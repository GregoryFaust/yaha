/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    QueryState.c    Includes code for Query State Structure.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include "Math.h"
#include "FragsClumps.inl"

/////
//     QueryState.
/////

// Note that the initialization and cleanup of QueryState objects is done is two phases.
// This is done to accomodate threads adjusting the size of their structures depending on the size of the queries they see.
// To make things simple, though, we will make all of the query read buffers big enough to handle the max possible query.
// This is currently 32KBases, so even six byte buffers this length is not a burden.
// This will keep us from having to change allocations and move the bases read so far in a query during query reading.
// Instead, we will be able to READ in any query up to the hard maximum, release the file lock, then change the other
//     structure sizes if they are two small.
// The structures that change a LOT with greater query length are those that are associated with fragments and DP.
// The former changes a lot because the max possible number of fragments is queryLength * maxHits.
// The latter change a lot because of the quadratic (or at least bandwidth * length) nature of the structures.
QueryState_t * makeQueryState(AlignmentArgs_t * AAs)
{
    QueryState_t * QS = (QueryState_t * )(malloc(sizeof(QueryState_t)));
    QS->AAs = AAs;
    QS->coverage = (BOOL *)(calloc(AAs->maxQueryLength, sizeof(BOOL)));
    QS->reverseBuf = (char *)(malloc(AAs->maxQueryLength));
    QS->reverseBufBase = QS->reverseBuf;
    QS->reverseBuf[AAs->maxQueryLength - 1] = '0';
    QS->reverseCodeBufBase = (char *)(malloc(AAs->maxQueryLength));
    QS->forwardBuf = (char *)(malloc(AAs->maxQueryLength));
    QS->qualBuf = (char *)(malloc(AAs->maxQueryLength));
    QS->queryID = (char *)(malloc(MAX_QUERY_ID_LEN));
    QS->forwardCodeBuf = (char *)(malloc(AAs->maxQueryLength));
    QS->clumps = NULL;
    //  Indicate we haven't read in a query yet.
    QS->queryLen = 0;
    QS->reversed = FALSE;
    // Initialize the Clump related structures.
    // This is done here because it doesn't allocate any memory,
    //     and never changed when the max size of a query changes for this thread.
    initClumpStructures(QS);
    // Similarly, the random number generator does not change with query length.
    QS->randGen = makeRandState();
    return QS;
}

void disposeQueryState(QueryState_t * QS)
{
    resetQueryState(QS);
    free(QS->reverseBufBase);
    free(QS->reverseCodeBufBase);
    free(QS->forwardCodeBuf);
    free(QS->forwardBuf);
    free(QS->qualBuf);
    free(QS->queryID);
    free(QS->coverage);
    // Since we initialize these in makeQueryState, we need to free them here.
    cleanUpClumpStructures(QS);
    if (QS->randGen != NULL) free(QS->randGen);
    free(QS);
}

// The next two are used when query length changes.
void initializeQueries (QueryState_t * QS)
{
    int maxQueryLength = QS->maxQueryLength;
    int maxRO_QO_PairCount = maxQueryLength * QS->AAs->maxHits;

    // Make these arrays here where we know how big it needs to be.
    QS->usedFrags = makeBitArray(maxRO_QO_PairCount);
    assert(QS->usedFrags);

    QS->fragArray = (Fragment_t *)malloc(maxRO_QO_PairCount*sizeof(Fragment_t));
    assert(QS->fragArray);

    QS->offsetCounts = (OffsetCount_t *)malloc(maxQueryLength*sizeof(OffsetCount_t));
    assert(QS->offsetCounts);

    QS->sortHeap = makeHeap(maxRO_QO_PairCount);

    // Initialize any thing the graph algorithms need.
    QS->tempFragArrayLen = 0;
    QS->tempFragArray = NULL;
    QS->fragGraphNodesLen = 0;
    QS->fragGraphNodes = NULL;

    // Do this last to increase chances that a realloc, if needed, will succeed without a copy.
    DPInit(QS);
}

void finalizeQueries(QueryState_t * QS)
{
    if (QS->usedFrags != NULL) disposeBitArray(QS->usedFrags);
    if (QS->fragArray != NULL) free(QS->fragArray);
    if (QS->offsetCounts != NULL) free(QS->offsetCounts);
    disposeHeap(QS->sortHeap);
    if (QS->tempFragArray != NULL) free(QS->tempFragArray);
    if (QS->fragGraphNodes != NULL) free(QS->fragGraphNodes);

    DPCleanUp(QS);
}

// Rest the coverage state for the query.
void setCoverage(QueryState_t * QS, int start, int len, BOOL value)
{
    // Reset the coverage information.
    memset(QS->coverage + start, value, len);
}

// Clean up the clump linked list and the coverage state for the query to get ready for the next one.
void resetQueryState(QueryState_t * QS)
{
    Clump_t * clump = QS->clumps;
    while (clump != NULL)
    {
        Clump_t * next = clump->next;
        disposeClump(clump, QS);
        clump = next;
    }
    QS->clumps = NULL;
    QS->clumpCount = 0;

    setCoverage(QS, 0, QS->queryLen, FALSE);
}

// Check for a contiguous run of covered locations on the query.
BOOL checkCoverage(QueryState_t * QS, int start, int end, int needed, BOOL value)
{
    int runLen = 0;
    int i;
    for (i = start; i<=end; i++)
    {
        if (QS->coverage[i] == value) runLen += 1;
        else runLen = 0;
        if (runLen >=needed) return TRUE;
    }
    return FALSE;
}

// The nest two add clumps into the linked list of clumps for the query.
void addClump(QueryState_t * QS, Clump_t * clump)
{
    setReversed(clump, QS->reversed);
    QS->clumps = pushClump(clump, QS->clumps);
    QS->clumpCount += 1;
}

void insertClumpAfter(QueryState_t * QS, Clump_t * prevClump, Clump_t * newClump)
{
    setReversed(newClump, QS->reversed);
    newClump->next = prevClump->next;
    prevClump->next = newClump;
    QS->clumpCount += 1;
}

// Generate the seed for the random number generator using the query codes for the query.
void generateRandomSeed(QueryState_t * QS)
{
    // two-bit compress the query codes, wrapping if necessary.
    int qoffset = 0;
    for (int i=0; i<5; i++)
    {
        UINT word = 0;
        for (int j=0; j<16; j++)
        {
            word = (word << 2) | (QS->forwardCodeBuf[qoffset] & 0x3);
            qoffset += 1;
            if (qoffset >= QS->queryLen) qoffset = 0;
        }
        QS->randGen->state[i] = word;
    }
}

void fprintQueryId(QueryState_t * QS, FILE * out)
{
    printChars(out, QS->queryID, 0, QS->queryIDLen);
    fprintf(out, " has length=%d, and strand=%c", QS->queryLen, (QS->reversed ? '-' : '+'));
}
