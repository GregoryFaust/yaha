/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    AlignArgs.c     Includes code for AlignmentArgs Structure.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include "Math.h"

/////
//     AlignmentArgs Structure.
/////

#define DEFAULT_ARG_VALUE -1

AlignmentArgs_t * makeAlignmentArgs ()
{
    AlignmentArgs_t * AAs = (AlignmentArgs_t *) malloc(sizeof(AlignmentArgs_t));

    // Set all the default values.
    AAs->gfileName = NULL;
    AAs->xfileName = NULL;
    AAs->qfileName = "stdin";
    AAs->ofileName = NULL;
    AAs->numThreads = 1;
    AAs->fastq = FALSE;

#ifdef QUERYSTATS
    AAs->qsfileName = NULL;
    AAs->queryStats = FALSE;
#endif

    // Any defaults that depend on other parameters must be set after the program args are processed.
    // So, mark them as DEFAULT_ARG_VALUE

    // Index parameters
    AAs->wordLen = 15;
    AAs->skipDist = 1;
    AAs->maxHits = DEFAULT_ARG_VALUE;

    // General Alignment Parameters
    AAs->maxGap = 50;
    AAs->maxIntron = DEFAULT_ARG_VALUE;
    AAs->minMatch = 25;
    AAs->minIdentity = 0.9;
    AAs->bandWidth = 5;
    AAs->maxDesert = 50;
    AAs->minRawScore = DEFAULT_ARG_VALUE;
    AAs->minNonOverlap = DEFAULT_ARG_VALUE;

    // These are the BWASW defaults.
    AAs->affineGapScoring = TRUE;
    AAs->GOCost = 5;
    AAs->GECost = 2;
    AAs->RCost  = 3;
    AAs->MScore = 1;
    AAs->XCutoff = 25;

    // These are for filtering of alignments against a "best" set of query spanning alignments.
    AAs->OQC = TRUE;
    AAs->OQCMinNonOverlap = DEFAULT_ARG_VALUE;
    AAs->BPCost = 5;
    AAs->maxBPLog = 5;
    // By default, output only primary alignments.
    // AAs->FBSMaxSimilar = 1;
    AAs->FBS = FALSE;
    AAs->FBS_PSLength = 0.90;
    AAs->FBS_PSScore = 0.90;

    // For now, just make the default the max.
    AAs->maxQueryLength = 32000;
    AAs->verbose = FALSE;

    AAs->outputBlast8 = FALSE;
    AAs->outputSAM = TRUE;
    AAs->hardClip  = TRUE;

    return AAs;
}

void disposeAlignmentArgs(AlignmentArgs_t * AAs)
{
    // The name of the index file and output file are always generated.
    if (AAs->gfileName != NULL) free(AAs->gfileName);
    if (AAs->ofileName != NULL) free(AAs->ofileName);
    free(AAs);
}

void printAlignmentArgs(AlignmentArgs_t * AAs, FILE * out)
{
    if (AAs->gfileName != NULL) fprintf(out, "Genome filename: %s\n", AAs->gfileName);
    if (AAs->xfileName != NULL) fprintf(out, "Index filename: %s\n", AAs->xfileName);
    if (AAs->qfileName != NULL) fprintf(out, "Query filename: %s\n", AAs->qfileName);
    if (AAs->ofileName != NULL) fprintf(out, "Output filename: %s\n", AAs->ofileName);
}

void postProcessAlignmentArgs(AlignmentArgs_t * AAs, BOOL query)
{
    // Make sure values not specified for AAs have reasonable values.
    if (AAs->maxIntron == DEFAULT_ARG_VALUE)
        AAs->maxIntron = AAs->maxGap;
    if (AAs->minRawScore   == DEFAULT_ARG_VALUE)
        AAs->minRawScore = AAs->minMatch;
    if (AAs->OQCMinNonOverlap == DEFAULT_ARG_VALUE)
        AAs->OQCMinNonOverlap = AAs->minMatch;
    if (AAs->OQCMinNonOverlap <= 0)
    {
        fprintf(stderr, "MNO parameter must be >=1.  MNO=1 will be used.\n");
        AAs->OQCMinNonOverlap = 1;
    }
    if (AAs->minNonOverlap == DEFAULT_ARG_VALUE)
    {
        AAs->minNonOverlap = AAs->OQCMinNonOverlap; // 1;
    }
    if (!AAs->affineGapScoring)
    {
        // Simulate edit distance by appropriate setting of the scoring paramenters.
        // This is not well tested.
        AAs->MScore = 1;
        AAs->RCost = AAs->GECost = 1;
        AAs->GOCost = 0;
    }
    // Calculate the min extension length requiring DP.
    // Since we always dp perfect extensions, the next bp is a mismatch.
    // DP is only needed if the length of the extension is long enough
    //      to have enough matches to make up for the cost of the mismatch.
    // The following is equivalent to ceiling(RCost/MScore) + 2;
    // We take the min of RCost and a single base gap cost.
    // But for any rational choice of the cost parameters, RCost <= AAs->GOCost + AAs->GECost.
    int len = 1;
    int score = 0;
    int target = MIN(AAs->RCost, AAs->GOCost + AAs->GECost);
    while (score <= target)
    {
        score += AAs->MScore;
        len += 1;
    }
    AAs->minExtLength = len;

    if (AAs->maxHits       == DEFAULT_ARG_VALUE)
    {
        if (query) AAs->maxHits = 650;
        else AAs->maxHits = SUINT_MAX_VALUE - 10;
    }
    else AAs->maxHits = MIN(AAs->maxHits, SUINT_MAX_VALUE - 10);
    // Only values between 1 and 9 make sense for this parameter.
    // For now, we will do the right thing with a warning.
    if (AAs->maxBPLog < 1)
    {
        fprintf(stderr, "MGDP parameter must be between 1 and 9 (inclusive). MGDP=1 will be used.\n");
        AAs->maxBPLog = 1;
    }
    if (AAs->maxBPLog > 9)
    {
        fprintf(stderr, "MGDP parameter must be between 1 and 9 (inclusive). MGDP=9 will be used.\n");
        AAs->maxBPLog = 9;
    }
}
