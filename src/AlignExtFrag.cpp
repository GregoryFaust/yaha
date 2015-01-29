/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    AlignExtFrag.cpp       Contains code to align and extend fragments.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include "Math.h"
#include "SW.inl"
#include "Math.inl"
#include "FragsClumps.inl"

/////
//     "Perfect" Fragment Extensions.
/////
//     These could be templated, but they are also called from C code.
/////

int extendFragmentForwardToStopPerfectly(Fragment_t *frag, char *refPtr, char *qPtr, int len)
{
    QOFF qOff = frag->endQueryOff + 1;
    ROFF rOff = fragEndRefOff(frag) + 1;
    int count = 0;
    while ((count < len) && codesEqual(qPtr[qOff+count], getFrom4Code(refPtr, rOff+count))) count += 1;
    if (count > 0) addLenToBack(frag, count);
    return count;
}

int extendFragmentBackwardToStopPerfectly(Fragment_t *frag, char *refPtr, char *qPtr, int len)
{
    QOFF qOff = frag->startQueryOff - 1;
    ROFF rOff = frag->startRefOff - 1;
    int count = 0;
    while ((count < len) && codesEqual(qPtr[qOff-count], getFrom4Code(refPtr, rOff-count))) count += 1;
    if (count > 0) addLenToFront(frag, count);
    return count;
}

/////
//     Clump Extensions.
/////
//     These are the only place the DP extension code is called.
//     We will always first do both forward and backward perfect extensions before DP extensions.
//     This allows us to skip some DP extensions if the extension length is not enough to overcome a mismatch.
//     Also, shorter DP extensions are faster DP extensions.
//     Consolidating clump extension here has also made it easier to be sure we don't run off the reference in either direction.
/////
//     Perhaps it is still possible to clean up the interface between this routine, and the ones in SW.cpp.
/////

// Extend the clump both forward and reverse directions.
// This will now do perfect extensions first, then DP extensions (possibly carefully).
template <BOOL goBack, BOOL goForw, BOOL goCarefully>
void extendClumpForwardReverseTemplated(Clump_t * clump, QueryState_t * QS)
{
    // Deconstruct the input args.
    SFragment_t * SFrag = SFragListPtr(clump)->head;
    Fragment_t * frag = getFragPtr(SFrag);
    EditOpList_t * list = EOLPtr(clump);
    AlignmentArgs_t * AAs = QS->AAs;
    char * refPtr = AAs->basePtr;
    char * queryPtr = getClumpQueryCodeBuf(QS, clump);
    int score = SFrag->score;

    int backLen = 0;
    int forwLen = 0;
    if (goBack)
    {
        backLen = MIN(frag->startQueryOff, frag->startRefOff);
        if (backLen > 0)
        {
            int newMatches = extendFragmentBackwardToStopPerfectly(frag, refPtr, queryPtr, backLen);
            if (newMatches > 0)
            {
                getFirstEditOp(list)->length += newMatches;
                score += newMatches * AAs->MScore;
                backLen -= newMatches;
            }
        }
    }
    if (goForw)
    {
        QOFF qlen = (QS->queryLen - 1) - frag->endQueryOff;
        ROFF rlen = AAs->maxROff - fragEndRefOff(frag);
        forwLen = MIN(qlen, rlen);
        if (forwLen > 0)
        {
            int newMatches = extendFragmentForwardToStopPerfectly(frag, refPtr, queryPtr, forwLen);
            if (newMatches > 0)
            {
                getLastEditOp(list)->length += newMatches;
                score += newMatches * AAs->MScore;
                forwLen -= newMatches;
            }
        }
    }

    SUINT addedQLen, addedRLen;
    if (goBack && backLen >= AAs->minExtLength)
    {
        int newScore;
        if(goCarefully)
            newScore = findAGSBackwardExtensionCarefully(QS, frag->startRefOff-1, queryPtr, frag->startQueryOff-1,
                                                         backLen, list, score, &addedQLen, &addedRLen);
        else
            newScore = findAGSBackwardExtension(QS, frag->startRefOff-1, queryPtr, frag->startQueryOff-1,
                                                backLen, list, &addedQLen, &addedRLen);
        if (newScore > 0)
        {
            score += newScore;
            addQLenToFront(frag, addedQLen);
            addRLenToFront(frag, addedRLen);
        }
    }
    if (goForw && forwLen >= AAs->minExtLength)
    {
        int newScore;
        if(goCarefully)
            newScore = findAGSForwardExtensionCarefully(QS, fragEndRefOff(frag)+1, queryPtr, frag->endQueryOff+1,
                                                        forwLen, list, score, &addedQLen, &addedRLen);
        else
            newScore = findAGSForwardExtension(QS, fragEndRefOff(frag)+1, queryPtr, frag->endQueryOff+1,
                                               forwLen, list, &addedQLen, &addedRLen);
        if (newScore > 0)
        {
            score += newScore;
            addQLenToBack(frag, addedQLen);
            addRLenToBack(frag, addedRLen);
        }
    }

    SFrag->score = score;
}

void extendClumpForwardReverse(Clump_t * clump, QueryState_t * QS)
{
    extendClumpForwardReverseTemplated<TRUE, TRUE, FALSE>(clump, QS);
}

void extendClumpForwardReverseCarefully(Clump_t * clump, QueryState_t * QS, BOOL goBack, BOOL goForw)
{
    if      (goBack && goForw) extendClumpForwardReverseTemplated<TRUE, TRUE, TRUE>(clump, QS);
    else if (goBack)           extendClumpForwardReverseTemplated<TRUE, FALSE, TRUE>(clump, QS);
    else                       extendClumpForwardReverseTemplated<FALSE, TRUE, TRUE>(clump, QS);
}


/////
//    The general routine for using DP to align a fragment.
//    This will use banded or unbanded DP depending on various sizes of things.
/////

SFragment_t * makeAndAlignSFragmentToFillGap(SFragment_t * SFrag1, SFragment_t * SFrag2, QueryState_t * QS, Clump_t * clump)
{
    Fragment_t * frag1 = getFragPtr(SFrag1);
    Fragment_t * frag2 = getFragPtr(SFrag2);

    // See if there is any space left, and if so, make an SFrag for it.
    SUINT qGap = calcQueryGap(frag1, frag2);
    SUINT rGap = calcRefGap(frag1, frag2);
    if (qGap == 0 && rGap == 0) return NULL;

    // Create SFrag and set offsets.
    SFragment_t * newSFrag = makeSFragment(QS);
    Fragment_t * newFrag = getFragPtr(newSFrag);
    newFrag->startQueryOff = frag1->endQueryOff + 1;
    newFrag->endQueryOff   = frag2->startQueryOff - 1;
    newFrag->startRefOff   = fragEndRefOff(frag1) + 1;
    fragSetEndRefOff(newFrag, frag2->startRefOff - 1);

    // Use a cascade for the case analysis, keeping the predicates as simple as possible.
    // We also want to avoid expensive DP whenever we can know the edit ops just by the gap lengths.
    // Due to SNPs, the 1,1 case is very common, and given we have perfectly extended, it has to be a mismatch.
    EditOpList_t * list = getEOLPtr(newSFrag);
    AlignmentArgs_t * AAs = QS->AAs;
    // We need to do these in the correct order to get the right edit ops.
    // Start with the lone indels
    if (qGap == 0)
    {
        addEditOpToFront(list, DeleteEditOp, rGap);
        newSFrag->score = calcGapCost(rGap, AAs);
    }
    else if (rGap == 0)
    {
        addEditOpToFront(list, InsertEditOp, qGap);
        newSFrag->score = calcGapCost(qGap, AAs);
    }
    // Now the very common single mismatch due to a snp or mis-sequencing.
    else if (rGap == 1 && qGap == 1)
    {
        addEditOpToFront(list, ReplaceEditOp, 1);
        newSFrag->score = -AAs->RCost;
    }
    // This is actually the most common case, but the test for it is complex.
    // Therefore, put it here at the end of the cascade.
    else //if (rGap > 0 && qGap > 0)
    {
        char * queryPtr = getClumpQueryCodeBuf(QS, clump);

        // Use banded DP when it requires fewer cells calculated than full DP.
        // This calculation is slightly oversimplified, and will slightly overuse full DP.
        // But full DP is faster than banded for similar number of cells calculated.
        int lenDiff = ABS(qGap - rGap);
        if (lenDiff + AAs->bandWidth * 2 + 1 < rGap)
        {
            newSFrag->score = findAGSAlignmentBanded(QS, newFrag->startRefOff, rGap,
                                                     queryPtr, newFrag->startQueryOff, qGap, list);
        }
        else
        {
            newSFrag->score = findAGSAlignment(QS, newFrag->startRefOff, rGap,
                                               queryPtr, newFrag->startQueryOff, qGap, list);
        }
    }

#ifdef DEBUG
    fprintf(stderr, "FRAG %d %d %d ", qGap, rGap, newSFrag->score);
    printEditOpList(list);
    fprintf(stderr, "\n");
#endif

    return newSFrag;
}
