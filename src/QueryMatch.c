/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    QueryMatch.c    Query Helper Functions, especially top level routines for forming fragments, then clumps, then alignments.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include <time.h>
#include "Math.h"
#include "Math.inl"
#include "FragsClumps.inl"
// #include "Timing.inl"

// This next include brings in all the code for the binary heap used as a priority queue for hash seed hits.
// This is much faster than putting all the seed hits into an array and sorting them by RO to form fragments.
// The complexity for latter is nlog(n) in the number of seed hits.
// Using a priority queue reduces complexity to order nlog(QL).
// This is because only one hit (RO) per query offset (hash key) need be in the heap at once.
// This only works because the ROs for a given seed are sorted in increasing RO order in the index file.
#include "QueryHeap.inl"

////////////////////////////////////////////
// Query Helper Functions.
////////////////////////////////////////////

// TODO we need a strategy for these global timers given they won't work when using threads.
// setupTimers();
// defineTimer(alignClumpsTimer);
// defineTimer(scoreClumpsTimer);
// defineTimer(junkTimeTimer);

// Merge the seed matches into an implicit array sorted first by diagonal, then by starting query offset.
// The array is implicit because we condense the seeds on the same diagonal is overlapping or abutting query offsets
//     into fragments (extended seeds) as we go along.
// Unfortunately, diagonals that wrap around to very large numbers will block all other diags from entering the queue
//     until well after their rightful spot, thereby separating them from other fragments they should join.
// While this is rare, it is a case we must handle.
// To fix this, we stuff all of these aberrant ones into the heap early so the normal ones can get in.
// Thankfully, all diagonals that wrap will occur first in the list of ROFFS for a given query offset if at all.
// Therefore, we need not slow down the main sort loop with this problem (other than with a slightly larger heap size).
int findFragmentsSort (AlignmentArgs_t * AAs, QueryState_t * QS, int matchCount)
{
    ROFF * ROAPtr = AAs->ROAPtr;
    // Load up heap with initial values.
    for (int i=0;i<matchCount; i++)
    {
        OffsetCount_t * OC = &QS->offsetCounts[i];
        if (OC->count == 0) continue;
        ROFF roff = ROAPtr[OC->sOffset];
        // We only need stuff in the wrapping diags until we get the first normal one.
        while (roff < i)
        {
            heapAddAtEnd(QS->sortHeap, roff - i, i);
            OC->newCount += 1;
            roff = ROAPtr[OC->sOffset + OC->newCount];
        }
        heapAddAtEnd(QS->sortHeap, roff - i, i);
    }
    fixHeap(QS->sortHeap);

    // Now start main loop pulling a value off of the heap and adding a new one until all values used.
    // As the new RO-QO pairs are pulled off the heap, coelesce them into fragments.
    // This allows us to never build the RO-QO arrays!
    // We will start by half unwinding the first loop.
    int wordLen = AAs->wordLen;
    int fragCount = 0;
    Fragment_t * curFrag = QS->fragArray + fragCount;
    ROFF curDiag;
    QOFF curQO;
    decodeHeapTop(QS->sortHeap, &curDiag, &curQO);
    fragSetSOs(curFrag, curQO, curDiag + curQO);
    QOFF curEQO = curQO + wordLen;
    ROFF nextDiag;
    QOFF nextQO;
    while (!isHeapEmpty(QS->sortHeap))
    {
        decodeHeapTop(QS->sortHeap, &nextDiag, &nextQO);

        // Update the heap.
        OffsetCount_t * OC = &QS->offsetCounts[nextQO];
        OC->newCount += 1;
        // If we have another entry for this QO, then add it into the heap.
        // If not, just put the last thing on the heap at the head and let it percolate out to recreate heap invariant.
        fastHeapify(QS->sortHeap, (OC->newCount >= OC->count) ? removeHeapLeaf(QS->sortHeap)
                                                              : encodeHeapItem(ROAPtr[OC->sOffset + OC->newCount]-nextQO, nextQO));

        // Check if we need to start a new frag.
        if (nextDiag != curDiag || nextQO > curEQO)
        {
            // First finish the current frag.
            fragSetEOs(curFrag, curEQO-1);
            fragCount += 1;
            // Now start new frag.
            curFrag = QS->fragArray + fragCount;
            curDiag = nextDiag;
            curQO = nextQO;
            fragSetSOs(curFrag, curQO, curDiag + curQO);
            curEQO = curQO + wordLen;
        }
        else
        {
            // Continue with the current frag.
            curEQO = nextQO + wordLen;
        }
    }
    // Finish the last half unroll by completing the last frag.
    fragSetEOs(curFrag, curEQO-1);
    fragCount += 1;
    return fragCount;
}

/////
//     Start of helper routines for forming fragments into clumps (potential alignments).
/////

// The contract is that this will return exactly the frags that can be reached in hops
//     of less than maxGap.
// Apparently, this is no longer used.
int findAlignableFragsBack(Fragment_t * fragArray, int FNLB, int startFragNum, int maxGap)
{
    Fragment_t * curStartFrag = fragArray+startFragNum;
    ROFF curStartDiag = fragDiag(curStartFrag);
    for (int i=startFragNum; i>=FNLB; i--)
    {
        Fragment_t * curFrag = (fragArray + i);
        ROFF curDiag = fragDiag(curFrag);
        if (absDiagDiff(curStartDiag, curDiag) > maxGap) return i + 1;
        curStartDiag = curDiag;
    }
    return FNLB;
}

// The contract is that this will return exactly the frags that can be reached in hops
//     of less than maxGap.
int findAlignableFragsForw(Fragment_t * fragArray, int FNUB, int startFragNum, int maxGap)
{
    Fragment_t * curStartFrag = fragArray+startFragNum;
    ROFF curStartDiag = fragDiag(curStartFrag);
    for (int i=startFragNum; i<=FNUB; i++)
    {
        Fragment_t * curFrag = (fragArray + i);
        ROFF curDiag = fragDiag(curFrag);
        if (absDiagDiff(curStartDiag, curDiag) > maxGap) return i - 1;
        curStartDiag = curDiag;
    }
    return FNUB;
}

// This stategy is no longer used.
// checkStartEndCoverage is used instead.
static inline BOOL isNotOverLappingClump(Fragment_t * curFrag, int minLeft, Clump_t * clump)
{
    // We only keep frags that are not overlapping on the query.
    int SQO = clumpSQO(clump);
    return (curFrag->startQueryOff <= SQO - minLeft || curFrag->endQueryOff >= clumpEQO(clump) + minLeft);
}

// Helper for checkStartEndCoverage
static inline BOOL checkQueryRegion(BOOL * coverage, int start, int end)
{
    for (int i = start; i<=end; i++) if (coverage[i]) return FALSE;
    return TRUE;
}

// This is used to determine which frags to keep in a reference region that has already formed at least one clump.
static inline BOOL checkStartEndCoverage(QueryState_t * QS, Fragment_t * frag, int minLeft)
{
    // By subtracting 1 here, we eliminate a bunch of adds later.
    minLeft -= 1;
    int SQO = frag->startQueryOff;
    int EQO = frag->endQueryOff;
    if (EQO - SQO < minLeft) return FALSE;
#ifdef OLDWAY
    // This requires BOTH ends of the fragment to have MNO free bases in order to be kept.
    if (!checkQueryRegion(QS->coverage, SQO, SQO + minLeft)) return FALSE;
    if (!checkQueryRegion(QS->coverage, EQO-minLeft, EQO)) return FALSE;
    return TRUE;
#else
    // This requires at least ONE end of the fragment to have MNO free bases in order to be kept.
    // This produces better results in a few percentage of queries with aligments in complex regions.
    // It results in more overall alignments in ALL mode (as more remaining frags support more minMatch clumps).
    // It has very little impact on the Accuracy test resuls, and no impact on the SV results (as all the events are too large).
    if ((checkQueryRegion(QS->coverage, SQO, SQO + minLeft)) || (checkQueryRegion(QS->coverage, EQO-minLeft, EQO))) return TRUE;
    return FALSE;
#endif
}

// Mark fragments that are in query regions that are already covered by a clump (potential alignment) in the current ref region.
// Use the bit array to keep track of which fragments have already been used.
void eliminateFragments(AlignmentArgs_t * AAs, QueryState_t * QS, Fragment_t * fragArray,
                        BitArray_t * usedFrags, int startFragNum, int endFragNum, Clump_t * clump)
{
    // Make sure we have something to work with.
    if (isEmpty(clump)) return;

    // Decide which fragments to get rid of, based on the portion of the query and reference that has just been aligned.
    // "Expand" the area to exclude from further consideration by the minNonOverlap parameter.
    for (int i=startFragNum; i<=endFragNum; i++)
    {
        if (get(usedFrags, i)) continue;
        Fragment_t * curFrag = fragArray + i;
        if (!checkStartEndCoverage(QS, curFrag, AAs->minNonOverlap)) setTRUE(usedFrags, i);
    }
}

/////
//     Start of main workhorse routine for forming fragments into clumps (potential alignments).
/////
//     Use graph algorithm to find best collection of fragments to include into a clump
//          for any reference region with more than one fragment in it.
/////

void processFragmentsGapped(AlignmentArgs_t * AAs, QueryState_t * QS, int fragCount)
{
    // Apparently, some queries get no hits!??!
    if (fragCount == 0) return;

    // Go through all of the fragments looking for the largest ones.
    int nextFragNum = 0;
    while (nextFragNum < fragCount)
    {
        // Find the reference region.
        // We can't combine fragments if their diagonals differ by more than maxGap.
        // Therefore, that criteria defines a reference region.
        int startFragNum = nextFragNum;
        int endFragNum   = findAlignableFragsForw(QS->fragArray, fragCount-1, nextFragNum, AAs->maxGap);

        int numOfFrags = 1 + endFragNum - startFragNum;

#ifdef DEBUG
        fprintf(stderr, "There are %d frags near a diagonal for query ID: ", numOfFrags);
        fprintQueryId(QS, stderr);
        fprintf(stderr, "\n");
        for (int j=startFragNum; j<=endFragNum; j++)
        {
            fprintf(stderr, "%d: ", j);
            printDebugFragmentExtended(QS->fragArray+j, AAs);
        }
#endif

#ifdef QUERYSTATS
        if (AAs->queryStats)
        {
            QS->diagRegionCount += 1;
            QS->diagRegionTotal += numOfFrags;
#ifdef MOREQUERYSTATS
            if (numOfFrags >= 500)
            {
                ROFF curDiag = fragDiag(QS->fragArray+startFragNum);
                int  curDiagFragCount = 0;
                for (int i=startFragNum; i<=endFragNum; i++)
                {
                    Fragment_t * newFrag = QS->fragArray + i;
                    ROFF newDiag = fragDiag(newFrag);
                    if (newDiag == curDiag) curDiagFragCount += 1;
                    else
                    {
                        curDiag = newDiag;
                        curDiagFragCount = 1;
                    }
                }
            }
#endif
        }
#endif

        // The most common case is to just have one fragment to deal with.
        // And this case doesn't need any fancy machinery.
        // So, just take care of it here.
        if (numOfFrags == 1)
        {
            Fragment_t * frag = QS->fragArray + startFragNum;
            if (fragMatchCount(frag) >= AAs->minMatch)
            {
                Clump_t * clump = makeClump(QS);
                addFragment(clump, frag, QS);
                addClump(QS, clump);
            }
        }
        else
        {
            // We have more than one fragment in the region.
            // Use the graph algorithm to decide which fragments to use to form a clump.
            processFragmentRangeUsingGraph(AAs, QS, QS->fragArray, QS->usedFrags, startFragNum, endFragNum);
        }
        nextFragNum = endFragNum + 1;
    }

    // The contract is that usedFrags is zeroed at creation, and each query cleans up its own mess.
    // This way we don't have to clear the whole array, but just the used portion.
    clearStart(QS->usedFrags, fragCount);
}

// Perform alignment and scoring of clumps.
void postProcessClumps(QueryState_t * QS, AlignmentArgs_t * AAs)
{
    // Loop through the clumps aligning and scoring.
    Clump_t * clump = QS->clumps;
    QS->clumps = NULL;
    QS->clumpCount = 0;
    while (clump != NULL)
    {
        // addToTimer(junkTimeTimer);
        Clump_t * nextClump = clump->next;
        alignClump(clump, AAs, QS);
        // addToTimer(alignClumpsTimer);

        scoreClump(clump, AAs, QS);
        // addToTimer(scoreClumpsTimer);
        // To keep high water marks of especially EditOps down to a rational number....
        //    free resources for any clump that won't print anyway.
        if (isScored(clump))
        {
            QS->clumps = pushClump(clump, QS->clumps);
            QS->clumpCount += 1;
        }
        else disposeClump(clump, QS);
        clump = nextClump;
    }
}

int printClumps(QueryState_t * QS, AlignmentArgs_t * AAs)
{
    int printedClumps = 0;
    Clump_t * clump = QS->clumps;
    while (clump != NULL)
    {
        printClump(clump, AAs, QS);
        printedClumps += 1;
        clump = clump->next;
    }
    return printedClumps;
}
