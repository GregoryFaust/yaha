/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    BaseSeq.c       Code for Base Sequence and Base Sequences Structures.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include "Math.h"

/////
//     BaseSequence.
/////
BaseSequence_t * makeBaseSequence()
{
    BaseSequence_t * BS = (BaseSequence_t *)(malloc(sizeof(BaseSequence_t)));
    BS->name = NULL;
    BS->nameLen = 0;
    return BS;
}

void disposeBaseSequence(BaseSequence_t * BS)
{
    free(BS);
}

void printBaseSequence(BaseSequence_t * BS, FILE * out)
{
    fprintf(out, "Base Sequence \"");
    printChars(out, BS->name, 0, BS->nameLen);
    fprintf(out, "\" has %u bases starting at offset %u\n", BS->length, BS->startingOffset);
}

/////
//     BaseSequences.
/////
BaseSequences_t * makeBaseSequences(int maxCountArg)
{
    BaseSequences_t * BSs = (BaseSequences_t *)malloc(sizeof(BaseSequences_t));
    BSs->maxCount = maxCountArg;
    BSs->sequences = (BaseSequence_t **)malloc(maxCountArg*sizeof(BaseSequence_t *));
    BSs->curCount = 0;
    BSs->totalBaseCount = 0;
    BSs->basePtr = NULL;
    return BSs;
}

void disposeBaseSequences(BaseSequences_t * BSs)
{
    int i;
    for (i=0; i<BSs->curCount; i++) disposeBaseSequence(BSs->sequences[i]);
    free(BSs->sequences);
    free(BSs);
}

void addBaseSequence(BaseSequences_t * BSs, BaseSequence_t * BS)
{
    // If we are the current limit, expand the array.
    if (BSs->curCount == BSs->maxCount)
    {
        BSs->maxCount *= 2;
        // fprintf(stderr, "Reallocating to more base sequences.  Now can handle %d\n", BSs->maxCount);
        BSs->sequences = (BaseSequence_t **)realloc(BSs->sequences, BSs->maxCount*sizeof(BaseSequence_t *));
        if (BSs->sequences == NULL)  fatalError("Insufficeint memory available to allocate more BaseSequence Objects.");
    }
    BSs->sequences[BSs->curCount] = BS;
    BSs->curCount += 1;
    BSs->totalBaseCount += BS->length;
}

int findBaseSequenceNum(BaseSequences_t * BSs, ROFF offset)
{
    int i;
    for (i=0; i<BSs->curCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        if (offset >= BS->startingOffset && offset < BS->startingOffset + BS->length) return i;
    }
    return -1;
}

BaseSequence_t * findBaseSequence(BaseSequences_t * BSs, ROFF offset)
{
    int i;
    for (i=0; i<BSs->curCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        if (offset >= BS->startingOffset && offset < BS->startingOffset + BS->length) return BS;
    }
    return NULL;
}

void printBaseSequences(BaseSequences_t * BSs, char * gfilename, FILE * out)
{
    fprintf(out, "There are %d base sequences in %s.\n", BSs->curCount, gfilename);
    int i;
    for (i=0; i<BSs->curCount; i++)
    {
        printBaseSequence(BSs->sequences[i], out);
    }
}

// When reading the .nib2 file, the offsets are in terms of bytes.
// But for indexing and query processing, we need the offset in bases.
void normalizeBaseSequences(BaseSequences_t * BSs)
{
    int i;
    for (i=0; i<BSs->curCount; i++) BSs->sequences[i]->startingOffset *= 2;
}

ROFF baseSequencesMaxROff(BaseSequences_t * BSs)
{
    BaseSequence_t * BS = BSs->sequences[BSs->curCount - 1];
    return BS->startingOffset + BS->length;
}
