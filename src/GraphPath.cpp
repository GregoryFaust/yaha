/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    GraphPath.c     File contains data structures and algorithms for finding a "best" path through a graph.
                             In paricular, the graph for picking best set of fragments for an alignment,
                             and the graph for picking best set of alignemnts for the optimal coverage set.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "Math.h"
#include "Math.inl"
#include "FragsClumps.inl"
#include "SW.inl"

////////////////////////////////
/////  Start of code for forming a best score path through a graph.
/////  For performance, we use separate but very similar code for two different cases.
////   We use a graph algorithm to find the best collection of fragments to include into clumps.
/////  And to find best collection of clumps that cover the query.
////   We also use a very simple sliding window approach to remove dups both before OQC and if OQC is not run.
////////////////////////////////

// The common aspects of the two graph algorithms/structures are as follows:
// There will be one graph node for each object (fragment or clump).
// All the "edges" will be implicit.
// One node will point to another node if it is possible
//     for the first node to come after the second node in the graph.
// We will act as if there is a dummy source node that will be pointed to by all possible
//     first nodes in the graph.  That is, to all the nodes
// But this node will also be implicit by just initializing each node score with their own score.
// The cost for an edge will be the sum of the cost of moving from the current node to the next node,
//     plus the advantage of adding the next node.
// Therefore, the cost can be negative.
// In addition, the node will know what object (fragment or clump) it corresponds to.
// It will know the cost of the best currently known path to this node.
// It will contain a pointer of the previous node on the optimal path.
// In general, we will cache all the information from the objects into the nodes that is needed
//    to run the graph algorithms so that we have the better cache locality of the nodes
//    vs. that of the objects they represent.

#define WorstScore (-(0x7fffff00))

///////////////////////////////////////////
// Start of the fragment specific code
///////////////////////////////////////////
// While forming clumps from fragments....
// We will assume that the fragments are sorted by diagonal, then by reference offset.
// We will assume that the fragments know their number of matching bases.
////////////////////////////////

struct fGraphNode
{
    // These are the items maintained by the graph algorithm.
    fGraphNode_t * bestPrev;       // The previous node in the best path to this node.
    Fragment_t   * frag;           // Pointer to the fragment represented by this node.  NULL if node is "dead".
    SINT           bestScore;      // The best score for this graph node.
    SINT           pathLength;     // The length of the best path to this node.
    QOFF           pathSQO;        // The SQO of the best path to this node.
    // The next items are just caches for information about the fragment to help cache locality.
    SINT          nodeScore;      // Cache of the starting score for the node cached from the fragment.
    ROFF          diag;           // The diagonal of a fragment.
    SINT          nodeLength;     // Cache of the total length for the node cached from the fragment.
    QOFF          SQO;            // Starting query offset of the fragment.
    QOFF          EQO;            // Ending query offset of the fragment.
};

// We have tried various criteria for the "best" fragment path.
// Emperically, the first one does better,
// presumably because it leaves more room for another alignment on what is left of the query.
// The second is kept around if we ever change our minds.
#ifdef TRUE
// For all fragment nodes with the best score, let's pick the one with the lowest end query off,
//     and if equal, the shorter query path length.
inline BOOL differentiateEqualFragNodesDuringBacktrack(fGraphNode_t * leftNode, fGraphNode_t * bestNode)
{
    // return (leftNode->EQO < bestNode->EQO);
    if (leftNode->EQO != bestNode->EQO) return (leftNode->EQO < bestNode->EQO);
    // Since the EQOs are equal, the one with the greater SQO has the shorter query path length.
    return (leftNode->pathSQO > bestNode->pathSQO);
}
#else
// For all fragment nodes with the best score, let's pick the one with the shortest query path length,
//     and if equal, the lower end query off.
inline BOOL differentiateEqualFragNodesDuringBacktrack(fGraphNode_t * leftNode, fGraphNode_t * bestNode)
{
    int lplen = ((int)leftNode->EQO - (int)leftNode->pathSQO);
    int bplen = ((int)bestNode->EQO - (int)bestNode->pathSQO);
    if (lplen != bplen) return (lplen < bplen);
    else return (leftNode->EQO < bestNode->EQO);
}
#endif

// Called for initialization of fragment graph nodes.
inline void initfGraphNode (fGraphNode_t * node, Fragment_t * frag, AlignmentArgs_t * AAs)
{
    // First do the generic node initialization.
    node->bestPrev = NULL;
    node->pathLength = 1;
    // Now the fragment specific initializations.
    node->frag = frag;
    node->diag = fragDiag(frag);
    node->nodeLength = fragMatchCount(frag);
    node->bestScore = node->nodeLength * AAs->MScore;
    node->SQO = frag->startQueryOff;
    node->EQO = frag->endQueryOff;
    node->pathSQO = node->SQO;
}

int compareFragsBySize(const void * arg1, const void * arg2)
{
    Fragment_t * frag1 = (Fragment_t *)arg1;
    Fragment_t * frag2 = (Fragment_t *)arg2;
    int f1Size = fragQueryLen(frag1);
    int f2Size = fragQueryLen(frag2);
    // If they are the same size, sort in ascending order of the starting query offset.
    if (f1Size == f2Size) return frag1->startQueryOff - frag2->startQueryOff;
    else return f2Size - f1Size;
}

void processBestFragmentPath(AlignmentArgs_t * AAs, QueryState_t * QS, fGraphNode_t * bestNode, int nodeCount, Clump_t * clump)
{
    for (fGraphNode_t * curNode = bestNode; curNode != NULL; curNode = curNode->bestPrev)
    {
        insertFragment(clump, curNode->frag, QS);
    }

    // Fragments might be resized during insertion.  So, need to get number of matching bases from clump.
    if (clump->matchedBases < AAs->minMatch)
        resetClump(clump, QS);
    else
        cleanUpClump(clump, QS);
}

int compareFragsByQueryOffsets(const void * arg1, const void * arg2)
{
    fGraphNode_t * node1 = (fGraphNode_t *)arg1;
    fGraphNode_t * node2 = (fGraphNode_t *)arg2;
    int SQO1 = node1->SQO;
    int SQO2 = node2->SQO;
    // Ascending in SQO
    int SQODiff = (SQO1 - SQO2);
    if (SQODiff != 0) return SQODiff;
    // Ascending in diag (SRO);
    return (diagDiff(node1->diag, node2->diag));
}

void buildBestClumpFromFragmentRange(AlignmentArgs_t * AAs, QueryState_t * QS, Fragment_t * fragArray,
                                     BitArray_t * usedFrags, int startFragNum, int endFragNum, Clump_t * clump)
{
    int numOfFrags = 1 + endFragNum - startFragNum;

    // See if the array needs to be larger.
    if (numOfFrags > QS->fragGraphNodesLen)
    {
        QS->fragGraphNodesLen = 2 * numOfFrags;
        QS->fragGraphNodes = (fGraphNode_t *)realloc(QS->fragGraphNodes, QS->fragGraphNodesLen*sizeof(fGraphNode_t));
    }

    // Initialize the nodes.
    // Because some fragments can already be excluded, we won't know the real count until we are done.
    int nodeCount = 0;
    for (int i=startFragNum; i<=endFragNum; i++)
    {
        if (get(usedFrags, i)) continue;
        fGraphNode_t * curNode = QS->fragGraphNodes + nodeCount;
        initfGraphNode(curNode, fragArray + i, AAs);
        nodeCount += 1;
    }
    // If all the fragments are used, we have no work to do.
    if (nodeCount == 0) return;

    // Sort fragments by increasing (plus strand normalized) query offsets.
    qsort(QS->fragGraphNodes, nodeCount, sizeof(fGraphNode_t), compareFragsByQueryOffsets);

    // Now start the graph algorithm.
    int bestScore = WorstScore;
    fGraphNode_t * bestNode = NULL;
    UINT maxGap = AAs->maxGap;
    int maxDesert = AAs->maxDesert;
    for (int i=0; i<nodeCount; i++)
    {
        fGraphNode_t * leftNode = QS->fragGraphNodes + i;
        int  lSQO = leftNode->SQO;
        int  lEQO = leftNode->EQO;
        ROFF lSRO = leftNode->diag + lSQO;
        // Equivalent to (lnode->diag + lnode->SQO) + (lnode->EQO - lnode->SQO);
        ROFF lERO = (leftNode->diag + leftNode->EQO);
        // We will go through the remainder of the list in reverse order so that we can stop once the SQO matches.
        for (int j=nodeCount-1; j>i; j--)
        {
            fGraphNode_t * rightNode = QS->fragGraphNodes + j;
            int rSQO = rightNode->SQO;
            if (rSQO == lSQO) break;
            // This simulates the absDiagDiff calculation.
            // But we can't break on this condition, but we can at least continue.
            // if (absDiagDiff(leftNode->diag, rightNode->diag) > maxGap) continue;
            UINT diagGap = absDiagDiff(leftNode->diag, rightNode->diag);
            if (diagGap > maxGap) continue;
            // Now simulate the canPrecede operation.
            // We've already tested the SQOs above, so just check the SROs.
            ROFF rSRO = rightNode->diag + rSQO;
            if (lSRO >= rSRO) continue;
            // Check if the desert is too big.
            int desert = MIN(calcGap(lEQO, rSQO), calcGap(lERO, rSRO));
            // if (desert > AAs->maxDesert) continue;
            if (desert > maxDesert) continue;

            // Check to make sure we have at least one new base being added.
            int maxOverlap = MAX(calcOverlap(lEQO, rSQO), calcOverlap(lERO, rSRO));
            int newbases = (rightNode->nodeLength - maxOverlap);
            if (newbases < 1) continue;

            // We've passed all the tests to decide if this would be a valid new path.
            // Now we need to estimate the score for the "edge".
            // This is the benefit of the (nonoverlapping) bases in the right node - the cost of a diagonal gap (indel).
            int newScore = leftNode->bestScore + (newbases * AAs->MScore) + calcGapCost(diagGap, AAs);

            // Now decide if we have a better connection or not.
            // If the newscore is worse, we don't.
            // If the newscore is better, we do.
            // If the scores are equal, and the right node alone has the same score as the two nodes combined,
            //     don't form the path (leave the right node to stand on its own).
            // Otherwise, use smallest diag diff, and if equal, the smallest desert,
            //     and if still equal, the shortest path query length (higher pathSQO).
            if      (rightNode->bestScore > newScore) continue;
            else if (rightNode->bestScore == newScore)
            {
                fGraphNode_t * prevBestNode = rightNode->bestPrev;
                if (prevBestNode == NULL) continue;
                int diagCompare = absDiagDiff(leftNode->diag, rightNode->diag) - absDiagDiff(prevBestNode->diag, rightNode->diag);
                if (diagCompare > 0) continue;
                else if (diagCompare == 0)
                {
                    int gapCompare = calcGap(leftNode->EQO, rightNode->SQO) - calcGap(prevBestNode->EQO, rightNode->SQO);
                    if (gapCompare > 0) continue;
                    else if (gapCompare == 0 && (leftNode->pathSQO <= prevBestNode->pathSQO)) continue;
                }
            }
            // We've passed the gauntlet, so we have a new best path for the right node.
            rightNode->bestScore = newScore;
            rightNode->bestPrev = leftNode;
            rightNode->pathLength = leftNode->pathLength + 1;
            rightNode->pathSQO = leftNode->pathSQO;
        }
        // Since this only cares about the left node, it may not matter if it is before or after the inner loop.
        if (leftNode->bestScore < bestScore) continue;
        if (leftNode->bestScore > bestScore || differentiateEqualFragNodesDuringBacktrack(leftNode, bestNode) > 0)
        {
            bestNode = leftNode;
            bestScore = leftNode->bestScore;
        }
    }

    processBestFragmentPath(AAs, QS, bestNode, nodeCount, clump);

}

void processFragmentRangeUsingGraph(AlignmentArgs_t * AAs, QueryState_t * QS, Fragment_t * fragArray,
                                    BitArray_t * usedFrags, int startFragNum, int endFragNum)
{
    // Reset query coverage to none.
    setCoverage(QS, 0, QS->queryLen, FALSE);
    while (TRUE)
    {
        Clump_t * clump = makeClump(QS);
        buildBestClumpFromFragmentRange(AAs, QS, fragArray, usedFrags, startFragNum, endFragNum, clump);
        if (isEmpty(clump))
        {
            disposeClump(clump, QS);
            return;
        }
        // Add the clump to the clump list, and mark its coverage on the query.
        setCoverage(QS, clumpSQO(clump), clumpQueryLen(clump), TRUE);
        // Eliminate any frags that dont cover an uncovered portion of query.
        eliminateFragments(AAs, QS, fragArray, usedFrags, startFragNum, endFragNum, clump);
        addClump(QS, clump);
    }
}

///////////////////////////////////////////////////
//     Start of clump specific code.
///////////////////////////////////////////////////

typedef struct cGraphNode cGraphNode_t;
struct cGraphNode
{
    // The graph algorithm fields
    cGraphNode_t * bestPrev;       // The previous node in the best path to this node.
    Clump_t      * clump;          // Pointer to the clump represented by this node.  NULL if node is "dead".
    SINT           bestScore;      // The best score for this graph node.
    SINT           pathLength;     // The length of the best path to this node.

    // The next items are just caches for information about the clump stored in the node to improve cache performance
    // The reference offsets are used for removing dups and to calculate break point penalty.
    // The query offsets are used for everything.
    ROFF           SRO;            // Starting reference offset.
    ROFF           ERO;            // Ending reference offset.
    QOFF           SQO;            // Starting query offset.  For clumps, this will always be relative to the + strand.
    QOFF           EQO;            // Ending query offset.    For clumps, this will always be relative to the + strand.
    // Now the cached values used during the graph algorithm
    SINT           nodeLength;     // Cache of the total length for the alignment cached from the clump.
    SINT           nodeScore;      // Cache of the starting score for the node cached from the clump.

    // This is used in the accurate overlap score calculations.
    SUINT          qLenInOQC;      // The number of bases from this clump in the best path the node is in.

    // These are needed for removing dups and to calculate break point penalty.
    BOOL          reversed;        // The strand that the clump is aligned to.
    UBYTE         seqNum;          // The sequence (chromosome) of the alignment.
};

inline void cleanUpGraph(cGraphNode_t * cGraphNodes, int nodeCount)
{
    free(cGraphNodes);
}

inline void setDead(cGraphNode_t * node)
{
    node->clump = NULL;
}

inline int isDead(cGraphNode_t * node)
{
    return (node->clump == NULL);
}

// Called for initialization of clump graph nodes.
inline void initcGraphNode (AlignmentArgs_t * AAs, QueryState_t * QS, cGraphNode_t * node, Clump_t * clump)
{
    // First do the generic node initialization.
    node->bestPrev = NULL;
    node->pathLength = 1;
    // Now the clump specific initializations.
    node->clump = clump;
    // We now keep track of the number of matched bases in a clump during fragment insertion.
    // TODO.  NO! We only get the clump score from fragment insertion for prefiltering, which we NO LONGER DO.
    //        Otherwise, this clumpScore comes from doing the alignment, where it has the real score.
    int clumpScore;
    clumpScore = clump->totScore;
    node->bestScore = node->nodeScore = clumpScore;
    node->nodeLength = clump->totLength;
    node->SQO = clumpPlusSQO(clump, QS);
    node->EQO = clumpPlusEQO(clump, QS);
    node->SRO = clumpSRO(clump);
    node->ERO = clumpERO(clump);
    node->reversed = isReversed(clump);
    node->qLenInOQC = clumpQueryLen(clump);
    node->seqNum = findBaseSequenceNum(AAs->BSs, node->SRO);
}

// We implement our own sort instead of using qsort so that we can use a random number generator as an argument.
// We can't use a global for this, as we wish to have the output of OQC stable for query file order and number of threads.
// So, each query/thread will have its own random number generator.

#ifdef TRUE
// This only works with 16 bit query offsets.
// So we will keep the old version around for when we go to 32 bit query offsets.
// The use of negatives allows for descendng sort.
// So, this is ascending on SQO, then descending on EQO, then descending on score.
// When otherwise equal, flip a coin to add stocastic (but query stable) result.
// The descending SQO sort allows easy subsumption test.
// Descending score tends to put the primary alignment first in a list of "equal" alignements.
inline UINT64 getCompareKey(cGraphNode_t * node)
{
    return (((((UINT64)node->SQO) << 16) + ((SUINT)-(SINT)node->EQO)) << 16) + ((SUINT)-node->nodeScore);
}

inline BOOL graphNodeLessThan(cGraphNode_t * node1, cGraphNode_t * node2, randState_t * randGen)
{
    UINT64 key1 = getCompareKey(node1);
    UINT64 key2 = getCompareKey(node2);
    if (key1 == key2) return (getRandBits(randGen) & 0x1);
    return (key1 < key2);
}

#else

// This would need to be converted to return a BOOL.
// Or leave this as is, and check in myQuickSortHelper for < 0 ala qsort compare routines.
// int graphNodeLessThan(cGraphNode_t * node1, cGraphNode_t * node2, randState randGen)
int CompareGraphNodeLessThan(cGraphNode_t * node1, cGraphNode_t * node2)
{
    // Sort in increasing query starting offsets.
    int SQO1 = node1->SQO;
    int SQO2 = node2->SQO;
    int SQODiff = (SQO1 - SQO2);
    if (SQODiff != 0) return SQODiff;
    // Sort in decreasing query ending offsets.
    // This allows an easy test for subsumption.
    int EQO1 = node1->EQO;
    int EQO2 = node2->EQO;
    int EQODiff = (EQO2 - EQO1);
    if (EQODiff != 0) return EQODiff;
    // Sort by decreasing node score.
    // This should tend to place the primary alignment first among a set of equal query offsets.
    // Unless they also have equal scores!
    int scoreDiff = (node2->nodeScore - node1->nodeScore);
    if (scoreDiff != 0) return scoreDiff;
    // Try to randomize the otherwise completely equal nodes.
    // This helps to get different "similar" nodes in the primary set for different queries.
    int rand = getRandBits(randGen) & 0x1;
    return rand ? 1 : -1;
}
#endif

inline void swap(void * one, void * two, void * temp, int len)
{
    memcpy(temp, one, len);
    memcpy(one, two, len);
    memcpy(two, temp, len);
}

void myQuickSortHelper(cGraphNode_t * array, int left, int right, cGraphNode_t * temp, randState_t * randGen)
{
    // Bottom out the recursion.
    if (left >= right) return;
    // Choose the pivot.
    // For now, just choose the middle point.
    // That will work well for sorted or reverse sorted input.
    // Consider using the median of left, right, and middle.
    // But the comparison is so costly, that this might not be worth the effort.
    int pivotIndex = (left + right)/2;
    // Perform the pivot.
    swap(array+pivotIndex, array+right, temp, sizeof(cGraphNode_t));
    cGraphNode_t * pivotNode = array + right;
    int storeIndex = left;
    for (int i=left; i<right; i++)
    {
        if (graphNodeLessThan(array+i, pivotNode, randGen))
        {
            swap(array+i, array+storeIndex, temp, sizeof(cGraphNode_t));
            storeIndex += 1;
        }
    }
    swap(array+storeIndex, array+right, temp, sizeof(cGraphNode_t));
    // Perform the recursion on the two sub-arrays.
    myQuickSortHelper(array, left, storeIndex-1, temp, randGen);
    myQuickSortHelper(array, storeIndex+1, right, temp, randGen);
}

inline void myQuickSort(cGraphNode_t * array, int length, randState_t * randGen)
{
    cGraphNode_t temp;
    myQuickSortHelper(array, 0, length-1, &temp, randGen);
}

inline BOOL nodesAreDups(cGraphNode_t * node1, cGraphNode_t * node2)
{
    // We'll check the SRO first, as they are least likely to be equal!
    return (node1->SRO == node2->SRO &&
            node1->ERO == node2->ERO &&
            node1->reversed == node2->reversed &&
            node1->SQO == node2->SQO &&
            node1->EQO == node2->EQO);
}

// We heuristically throw away alignments that are within the query offset boundaries of a longer alignment
//    and have less than 1/8 its score (which likely means it is about 1/8 its length given high identity thresholds.
// This tends to help performance a LOT when there is a small portion of a query that is very low complexity,
//     and has very little downside otherwise.
#define DELETE_SUBSUMED
#ifdef DELETE_SUBSUMED
inline BOOL nodeIsSubsumed(cGraphNode_t * node1, cGraphNode_t * node2, int node1ScoreThreshold)
{
    return (node1->EQO > node2->EQO && node2->nodeScore < node1ScoreThreshold);
}
#else
inline BOOL nodeIsSubsumed(cGraphNode_t * node1, cGraphNode_t * node2, int node1ScoreThreshold)
{
    return FALSE;
}
#endif

int deleteSubsumedDups (QueryState_t * QS, cGraphNode_t * graphNodes, int nodeCount)
{
    int count = 0;
    for (int i=0; i<nodeCount; i++)
    {
        cGraphNode_t * curNode = &graphNodes[i];
        if (isDead(curNode)) continue;
        // We will copy the good nodes to the lowest contiguous section of the array.
        // When many nodes are subsumed, this may significantly increase cache performance.
        // And it eliminates the check for dead nodes in the main graph algorithm.
        if (count != i) graphNodes[count] = graphNodes[i];
        count += 1;
        // Calculate this only once for the inner loop.
        int curScoreThreshold = curNode->nodeScore/8;
        for (int j=i+1; j<nodeCount; j++)
        {
            cGraphNode_t * nextNode = &graphNodes[j];
            if (isDead(nextNode)) continue;
            // Given the sort is decreasing in EQO, we need not go on if the next node has a higher EQO.
            if (nextNode->EQO > curNode->EQO) break;
            if (nodeIsSubsumed(curNode, nextNode, curScoreThreshold) || nodesAreDups(curNode, nextNode))
            {
                Clump_t * clump = nextNode->clump;
                disposeClump(clump, QS);
                setDead(nextNode);
            }
        }
    }
    return count;
}

/////
//     FilterBySimilarity
/////
//     The graph is assumed to already have had dups removed.
//     We start with three inlines for calculating mapping quality.
/////

struct primaryAttributes
{
    int alignedQueryLength;
    int numOutputSecondaries;
    SINT secondScore;
    SINT thirdScore;
};
typedef struct primaryAttributes PrimaryAttributes_t;

inline void initPAs(PrimaryAttributes_t * PAs, cGraphNode_t * node)
{
    // This is a cache for the aligned query length of the primary
    //     so we don't have to keep recalculating for each possible secondary.
    PAs->alignedQueryLength = (1 + node->EQO - node->SQO);
    PAs->secondScore = 0;
    PAs->thirdScore = 0;
    PAs->numOutputSecondaries = 0;
}

inline void  memoPAsFromOverlappingNode(PrimaryAttributes_t * PAs, cGraphNode_t * simNode)
{
    // Allow the secondary to contribute to the mapping quality calculation of the corresponding primary.
    if (simNode->nodeScore > PAs->secondScore)
    {
        PAs->thirdScore = PAs->secondScore;
        PAs->secondScore = simNode->nodeScore;
    }
    else if (simNode->nodeScore > PAs->thirdScore)
    {
        PAs->thirdScore = simNode->nodeScore;
    }
}

inline void calcMQfromPAs (Clump_t * clump, PrimaryAttributes_t * PAs)
{
    if (PAs->secondScore == 0) clump->mapQuality = 250;
    else
    {
        double scoreRatio = MAX(((double)clump->totScore - PAs->secondScore), 0.0)/((double)clump->totScore);
        scoreRatio = scoreRatio * (1.0 + MAX(((double)clump->totScore - PAs->thirdScore), 0.0)/clump->totScore)/2.0;
        clump->mapQuality = (UBYTE)((250.0 * scoreRatio) + 0.5);
    }
    clump->numSecondaries = PAs->numOutputSecondaries;
}

void filterBySimilarity (AlignmentArgs_t * AAs, QueryState_t * QS, cGraphNode_t * graphNodes,
                         int nodeCount, cGraphNode_t * bestNode)
{
    // We want to mark clumps in the primary query coverage set as "Primary" so we can write an optional SAM field saying that.
    // We will also create a small array of graph nodes for the primary alignments to improve cache performance,
    //    but more importantly give us an index to use in secondary arrays storing various things.
    Clump_t * newClumpList = NULL;
    int newClumpCount = 0;

    int primeCount = bestNode->pathLength;
    cGraphNode_t * primaries = (cGraphNode_t *)malloc(primeCount * sizeof(cGraphNode_t));
    PrimaryAttributes_t * PA_Array = (PrimaryAttributes_t *)malloc(primeCount * sizeof(PrimaryAttributes_t));
    // We want to have the indexes go from left to right on the query.
    // Yet, the best path starts at the right and goes left.
    // Therefore, we will add the primaries into the array in reverse order.
    int primeIndex = primeCount - 1;
    for (cGraphNode_t * pathNode = bestNode; pathNode != NULL; pathNode = pathNode->bestPrev)
    {
        // Instead of copying things one field at a time, it will be faster and more robust to use memcpy.
        memcpy(&primaries[primeIndex], pathNode, sizeof(cGraphNode_t));
        initPAs(&PA_Array[primeIndex], pathNode);
        // Mark the clump as primary and put it on the output list.
        // Mark the node as dead so it will be ignored below.
        Clump_t * clump = pathNode->clump;
        setPrimary(clump, TRUE);
        clump->matchedPrimary = primeIndex + 1;
        newClumpCount += 1;
        newClumpList = pushClump(clump, newClumpList);
        // Mark as dead so we don't re-add it to the clump list.
        setDead(pathNode);
        // Get ready for next iteration.
        primeIndex -= 1;
    }

    // Now let's filter out any clumps that don't match any of the clumps in the query coverage.
    // We will do this by starting a new list, then replacing the list in the QS when we are done.
    // It is easier/faster to filter by iterating through the graph nodes than the clump list.
    // This is because the nodes already have cached the appropriate start and end QOs relative to the plus strand.

    double targetOverlap = AAs->FBS_PSLength;
    for (int i=0; i<nodeCount; i++)
    {
        cGraphNode_t * curNode = &graphNodes[i];
        if (isDead(curNode)) continue;
        Clump_t * clump = curNode->clump;

        int curSQO = curNode->SQO;
        int curEQO = curNode->EQO;
        int curQLen = (1 + curEQO - curSQO);

        // Loop through the optimal path seeing which path clump the current clump has most overlap with.
        int maxOverlap = 0;
        int maxIndex = 0;
        for (int i=0; i<primeCount; i++)
        {
            cGraphNode_t * pathNode = &primaries[i];
            int pathSQO = pathNode->SQO;
            int pathEQO = pathNode->EQO;
            // This formula is from Aaron Quinlan as used in intersectBed.
            int overlap = 1 + MIN(curEQO, pathEQO) - MAX(curSQO, pathSQO);
            if (overlap > maxOverlap)
            {
                maxOverlap = overlap;
                maxIndex = i;
            }
        }
        // Make sure this alignment overlaps at least one of the primaries!
        if (maxOverlap > 0)
        {
            PrimaryAttributes_t * PAs = &PA_Array[maxIndex];
            // Now allow the secondary to contribute to the mapping quality calculation of the corresponding primary.
            memoPAsFromOverlappingNode(PAs, curNode);

            // Now see if this secondary passed the FBS tests.
            // We now do this even when not in FBS mode, so that we can output the correct secondary count.

            // First test the score criteria.
            cGraphNode_t * pathNode = &primaries[maxIndex];
            if (((double)curNode->nodeScore)/pathNode->nodeScore >= AAs->FBS_PSScore)
            {
                int pathSQO = pathNode->SQO;
                int pathEQO = pathNode->EQO;
                int overlap = 1 + MIN(curEQO, pathEQO) - MAX(curSQO, pathSQO);
                int pathQLen = PAs->alignedQueryLength;
                double overlapD = overlap;
                // Test the length criteria.
                if (overlapD/curQLen >= targetOverlap && overlapD/pathQLen >= targetOverlap)
                {
                    PAs->numOutputSecondaries += 1;
                    if (AAs->FBS)
                    {
                        // Mark the secondary with the index of the associated primary.
                        clump->matchedPrimary = maxIndex + 1;
                        newClumpCount += 1;
                        newClumpList = pushClump(clump, newClumpList);
                        continue;
                    }
                }
            }
        }
        // If we make it here, this secondary will not be output.
        // So delete the corresponding clump.
        disposeClump(clump, QS);
    }
    // Update the QS with the new list information.
    QS->clumps = newClumpList;
    QS->clumpCount = newClumpCount;
    QS->primaryCount = primeCount;

    // Now go back and fill in the mapping quality for the primary alignments.
    for (int i=0; i<primeCount; i++)
    {
        cGraphNode_t * pathNode = &primaries[i];
        Clump_t * clump = pathNode->clump;
        PrimaryAttributes_t * PAs = &PA_Array[i];
        calcMQfromPAs(clump, PAs);
    }

    // Free the temp structures.
    free(primaries);
    free(PA_Array);
}

////////////////////////////////
// Optimal Query Coverage
////////////////////////////////

///////////////////////////////
// Start with a lot of routines needed to calculate accurate overlap scores along the various "best paths".
// This is a complicated operation because of the fact that the nodes can be included in more than one "best path".
// So, we have to do a lot of recalculation of the overlapping segments as we examine different potential "best paths"
///////////////////////////////

// Calculate the score for a length along the alignment using the actual editops, not just length.
template<BOOL forward>
inline int calcScoreForLength(EditOpList_t * list, int length, AlignmentArgs_t * AAs)
{
    EditOp_t * item;
    int QLen = 0;
    int AGS = 0;
    if (forward) item = getFirstEditOp(list);
    else         item = getLastEditOp(list);
    while (isValidEONode(item, list) && QLen < length)
    {
        EditOpCode op = item->opcode;
        int len = item->length;
        // Delete has to be handled separately, as it adds no query length.
        if (op == DeleteEditOp) AGS -= (AAs->GOCost + AAs->GECost * len);
        else
        {
            // We chop off the last edit op to exactly match the desired length.
            if (QLen + len > length) len = length - QLen;
            QLen += len;
            if      (op == MatchEditOp  ) AGS += (AAs->MScore * len);
            else if (op == ReplaceEditOp) AGS -= (AAs->RCost * len);
            else if (op == InsertEditOp ) AGS -= (AAs->GOCost + AAs->GECost * len);
        }
        if (forward) item = getNextEditOp(item, list);
        else         item = getPrevEditOp(item, list);
    }
    return AGS;
}

int calcScoreForLengthForward(EditOpList_t * list, int length, AlignmentArgs_t * AAs)
{
    return calcScoreForLength<TRUE>(list, length, AAs);
}

int calcScoreForLengthBackward(EditOpList_t * list, int length, AlignmentArgs_t * AAs)
{
    return calcScoreForLength<FALSE>(list, length, AAs);
}

int calcAccurateOverlapScore(cGraphNode_t * leftNode, cGraphNode_t * rightNode, int overlap, AlignmentArgs_t * AAs, BOOL * rightBest)
{
    // We first calculate the accurate overlap score the for the rightNode.
    // This is the easy one, because we don't have to worry about multiple clumps being involved.
    // For reverse strand alignments, the edit op list is reversed relative to the query offsets.
    int rightOverlapScore;
    EditOpList_t * rightEOL = EOLPtr(rightNode->clump);
    if (rightNode->reversed)
        rightOverlapScore = calcScoreForLengthBackward(rightEOL, overlap, AAs);
    else
        rightOverlapScore = calcScoreForLengthForward(rightEOL, overlap, AAs);

    // The left side is harder, as we may have to navigate the entire best path to get enough
    // length of best scoring alignment to cover the needed overlap length.
    // We can't just assume that the leftNode itself will contribute the best segments
    //    for the best non-overlapping score.
    // However, because we only ever compare two alignments at a time to discover the best segments,
    //    and we do this entire calculation from a (normalized) query offset sorted list,
    //    best segments can only be at the (normalized) right most end of alignments, or out of reach
    //    of the right node overlap.
    int pathOverlapScore = 0;
    int remainingOverlap = overlap;
    cGraphNode_t * curNode = leftNode;
    while (TRUE)
    {
        EditOpList_t * curEOL = EOLPtr(curNode->clump);
        int curReverseQLen = MIN(remainingOverlap, curNode->qLenInOQC);
        remainingOverlap -= curReverseQLen;
        int curOverlapScore;
        if (curNode->reversed)
            curOverlapScore = calcScoreForLengthForward(curEOL, curReverseQLen, AAs);
        else
            curOverlapScore = calcScoreForLengthBackward(curEOL, curReverseQLen, AAs);
        pathOverlapScore += curOverlapScore;
        // Get ready for the next iteration.
        if (remainingOverlap <= 0) break;
        curNode = curNode->bestPrev;
    }

    // Now we have two overlap scores; one from the rightNode, and one from the left node's best path.
    // We will return the lower of the two scores.
    // But we also need to set the qLenInOQC values accordingly.
    // BUT ONLY IF WE ACTUALLY USE THIS SCORE AS THE NEW BEST BLOCK!!!!
    // So, return a second boolean value that says which way was best to be used later.
    int retval;
    if (pathOverlapScore > rightOverlapScore)
    {
        *rightBest = FALSE;
        retval = rightOverlapScore;
    }
    else
    {
        *rightBest = TRUE;
        retval = pathOverlapScore;
    }
    return retval;
}

void cacehQlenInOQCPathReverse(cGraphNode_t * leftNode, cGraphNode_t * rightNode, int overlap, BOOL rightBest)
{
    if (rightBest)
    {
        // This is the harder case, as we have to update the best path to reset their qLenInOQC values.
        rightNode->qLenInOQC = 1 + rightNode->EQO - rightNode->SQO;
        int remainingOverlap = overlap;
        cGraphNode_t * curNode = leftNode;
        while (TRUE)
        {
            int curReverseQLen = MIN(remainingOverlap, curNode->qLenInOQC);
            curNode->qLenInOQC -= curReverseQLen;
            remainingOverlap -= curReverseQLen;
            // Get ready for the next iteration.
            if (remainingOverlap <= 0) break;
            curNode = curNode->bestPrev;
            // if (curNode == NULL) fatalError("No nodes remaing on best path to cover overlap.");
        }
    }
    else
    {
        // The easy case.  Just remove the overlap from the full query length.
        rightNode->qLenInOQC = (1 + rightNode->EQO - rightNode->SQO) - overlap;
    }
}

// This is called for each node in the outer loop.
// By the time the node is the "left node" in comparisons,
//    its best prev path never changes again, but the qLens in it still can (see below).
// Therefore, we must set the path lengths for each outer loop iteration (left node).
// Unfortunately, we will have to reset the entire length of the path.
// This is because nodes towards the beginning of the path can be pointed to by more than one node.
// Each such path can have different qLens in the OQC because of the different nodes
//      towards the ends of the paths.
// As we only have backward path ptrs, we will do a recursive descent, and calculate the path on the way back up.
// This will closely mirror the full block process in sumCoverageScores.
// TODO!!  We need to think of a faster way to do this.  What data structure??
// Can we cache a dirty bit that controls if we can bypass all this?

cGraphNode_t * cacheQlenInOQCPath(cGraphNode_t * rightNode, AlignmentArgs_t * AAs)
{
    int qLen = 1 + rightNode->EQO - rightNode->SQO;
    // First check if this node is at the bottom of the recursion.
    // This will also be true of nodes that are not in a path.
    // If so, initialize its qLenInOQC, and return.
    if (rightNode->bestPrev == NULL)
    {
        rightNode->qLenInOQC = qLen;
        return rightNode;
    }
    // Now we need to get the next node.
    cGraphNode_t * leftNode = cacheQlenInOQCPath(rightNode->bestPrev, AAs);
    int overlap = calcOverlap(leftNode->EQO, rightNode->SQO);
    BOOL rightBest;
    if (overlap > 0)
    {
        calcAccurateOverlapScore(leftNode, rightNode, overlap, AAs, &rightBest);
        cacehQlenInOQCPathReverse(leftNode, rightNode, overlap, rightBest);
    }
    // If there is no overlap, we need to set rightNode cache.
    else
    {
        rightNode->qLenInOQC = qLen;
    }
    return rightNode;
}

// This is called when a right node is tentatively set as the new right node for a path.
// We only need to set its qLenInOQC field.
// The remainder of the path fixup will occur when the prevBest node becomes the left node..
// And therefore that part of the path will not change again.
inline void cacheQlenInRightNode(cGraphNode_t * rightNode, int overlap, BOOL rightBest)
{
    int qLen = 1 + rightNode->EQO - rightNode->SQO;
    if (rightBest) rightNode->qLenInOQC = qLen;
    else           rightNode->qLenInOQC = qLen - overlap;
}

///////////////////////////////
// End of calculate accurate overlap score code.
///////////////////////////////

///////////////////////////////
// The start of the main function that performs OQC algorithm.
///////////////////////////////

// Helper to distinguish between nodes in paths as right most node of best path
//     when they have the same path score.
// We will pick the one with the shortest path length in terms of nodes.
inline BOOL differentiateEqualClumpNodes(cGraphNode_t * leftNode, cGraphNode_t * bestNode)
{
    return (leftNode->pathLength < bestNode->pathLength);
}

// Main work horse routine.
void postFilterBySimilarity(AlignmentArgs_t * AAs, QueryState_t * QS)
{
    // The special cases are so common and can be handled so much faster,
    //     that they are worth handling separately.

    // Just in case some query has no clumps.
    int nodeCount = QS->clumpCount;
    if (nodeCount < 1) return;

    // Because it is so common, handle a single node as a special case.
    if (nodeCount == 1)
    {
        Clump_t * clump = QS->clumps;
        setPrimary(clump, TRUE);
        clump->mapQuality = 250;
        clump->numSecondaries = 0;
        clump->matchedPrimary = 1;
        QS->primaryCount = 1;
        return;
    }

#ifdef MORETIMING
    clock_t startTempClock, endTempClock;
    startTempClock = clock();
#endif

    // First build the graph of the clumps.
    // Start by initializing the nodes.
    // We need two iterators, one down the list of clumps, and the other down the graph array.
    // So, use an aux count variable for the latter.
    int curNodeCount = 0;
    cGraphNode_t * graphNodes = (cGraphNode_t *)malloc(nodeCount*sizeof(cGraphNode_t));
    for (Clump_t * clump = QS->clumps; clump != NULL; clump = clump->next)
    {
        cGraphNode_t * curNode = &graphNodes[curNodeCount];
        initcGraphNode(AAs, QS, curNode, clump);
        curNodeCount += 1;
    }

    // We will sort the graphNodes by plus-strand-normalized query coordinates.
    // This is essentially a topological sort of the nodes for our graph algorithm.
    // Therefore, we can skip actually forming the edges and doing the topo sort.
    // We will instead essentially form a virtual edge and compare scores as we go!
    // We need a RNG to randomize equal alignments.
    myQuickSort(graphNodes, curNodeCount, QS->randGen);

#ifdef MORETIMING
    endTempClock = clock();
    fprintf(stderr, "Time to build and sort graph with %d nodes is %lf seconds.\n",
            curNodeCount, ((double)(endTempClock-startTempClock))/CLOCKS_PER_SEC);
    startTempClock = endTempClock;
#endif

    // Filter out subsummed alignments.
    // This is much easier to do once the graph nodes are sorted.
    // As we might have fewer nodes/clumps after this process, reset the nodeCount as well.
    // We should consider not eliminating dups here.
    // That is because there are few dups and surprisingly, given the AccurateOverlapCalculation,
    //     seemingly duplicate alignments can lead to different primary sets due to differences in their EditOpLists.
    curNodeCount = deleteSubsumedDups(QS, graphNodes, curNodeCount);
    // fprintf(stderr, "Filtered out %d subsumed alignments out of %d leaving %d .\n", nodeCount-curNodeCount, nodeCount, curNodeCount);

#ifdef MORETIMING
    endTempClock = clock();
    fprintf(stderr, "Time to filter out %d subsumed alignments out of %d is %lf seconds.\n",
            nodeCount-curNodeCount, nodeCount, ((double)(endTempClock-startTempClock))/CLOCKS_PER_SEC);
    startTempClock = endTempClock;
#endif

    // We are going to relax edges, and find the "best" node all in one pass.
    int bestScore = WorstScore;
    cGraphNode_t * bestNode = NULL;
    int minNonOverlap = AAs->OQCMinNonOverlap;
    int BPCost = AAs->BPCost;
    int MBPL = AAs->maxBPLog;
    int startj = 1;
    for (int i=0; i<curNodeCount; i++)
    {
        cGraphNode_t * leftNode = graphNodes + i;
        cacheQlenInOQCPath(leftNode, AAs);
        int leftSQO = leftNode->SQO;
        int leftEQO = leftNode->EQO;

        // We will keep track of the first place there is enough left hand side non-overlap.
        // Therefore, next rightNode can start from there and go up.
        // This keeps a lot of similar clumps without enough left overlap from being compared.
        BOOL foundstartj = FALSE;
        for (int j=startj; j<curNodeCount; j++)
        {
            cGraphNode_t * rightNode = graphNodes + j;
            int rightSQO = rightNode->SQO;
            if ((rightSQO - leftSQO) >= minNonOverlap)
            {
                // Only change startj if this is the first one that passes the SQO test.
                if (!foundstartj)
                {
                    startj = j;
                    foundstartj = TRUE;
                }
                // Proceed to check the EQO test.
                int rightEQO = rightNode->EQO;
                if ((rightEQO - leftEQO) >= minNonOverlap)
                {
                    // How to order the tests to make things go as fast as possible in the common cases,
                    //     and only slow down in the complex cases?
                    // First calc the new score ignoring overlap.
                    // If this is not good enough, then subtracting overlap can only make it worse.
                    SINT newScore = leftNode->bestScore + rightNode->nodeScore;
                    if (rightNode->bestScore > newScore) continue;
                    // Now calculate the break point penalty.
                    int BPP;
                    if (leftNode->seqNum == rightNode->seqNum)
                    {
                        // Extra nonsense to avoid the silly unsigned int roll-over for negative distances.
                        ROFF distance;
                        if      (leftNode->SRO > rightNode->ERO) distance = leftNode->SRO - rightNode->ERO;
                        else if (rightNode->SRO > leftNode->ERO) distance = rightNode->SRO - leftNode->ERO;
                        else distance = 0;
                        if (distance <= 10) BPP = BPCost;
                        else
                        {
                            double log = log10(distance);
                            if (log > MBPL) log = (double)MBPL;
                            BPP = (int)(log * BPCost + 0.5);
                        }
                    }
                    else
                        BPP = MBPL * BPCost;
                    newScore -= BPP;
                    if (rightNode->bestScore > newScore) continue;
                    // Now calculate the overlap
                    // If it is 0, we have nothing more to do.
                    // If not, calculate the overlap score.
                    // And test again if the current best score is still better.
                    int overlap = calcOverlap(leftEQO, rightSQO);
                    BOOL rightBest;
                    if (overlap > 0)
                    {
                        newScore -= calcAccurateOverlapScore(leftNode, rightNode, overlap, AAs, &rightBest);
                        if (rightNode->bestScore > newScore) continue;
                    }
                    // When we get here, we already know that the rightNode->bestScore <= newScore.
                    if (rightNode->bestScore < newScore ||
                        (rightNode->bestPrev != NULL && differentiateEqualClumpNodes(leftNode, rightNode->bestPrev)))
                    {
                        if (overlap > 0) cacheQlenInRightNode(rightNode, overlap, rightBest);
                        rightNode->bestScore = newScore;
                        rightNode->bestPrev = leftNode;
                        rightNode->pathLength = leftNode->pathLength + 1;
                    }
                }
            }
            // If NONE pass the SQO test, we are done needing the inner loop.
            // However we need to keep going up on the outer loop to find the left node with the "best" score.
            if (!foundstartj)
            {
                startj = curNodeCount;
            }
        }
        // Since this only cares about the left node, it may not matter if it is before or after the inner loop.
        if (leftNode->bestScore < bestScore) continue;
        if (leftNode->bestScore > bestScore || (bestNode != NULL && differentiateEqualClumpNodes(leftNode, bestNode)))
        {
            bestNode = leftNode;
            bestScore = leftNode->bestScore;
        }
    }

#ifdef MORETIMING
    endTempClock = clock();
    fprintf(stderr, "Time to find best path through graph is %lf seconds.\n", ((double)(endTempClock-startTempClock))/CLOCKS_PER_SEC);
    startTempClock = endTempClock;
#endif

    filterBySimilarity (AAs, QS, graphNodes, curNodeCount, bestNode);

#ifdef MORETIMING
    endTempClock = clock();
    fprintf(stderr, "Time to filter alignments based on best path through graph is %lf seconds.\n",
            ((double)(endTempClock-startTempClock))/CLOCKS_PER_SEC);
    startTempClock = endTempClock;
#endif

    cleanUpGraph(graphNodes, curNodeCount);

#ifdef MORETIMING
    endTempClock = clock();
    fprintf(stderr, "Time to cleanup graph is %lf seconds.\n", ((double)(endTempClock-startTempClock))/CLOCKS_PER_SEC);
#endif
}

/////
//     Filter duplicate alignments in the case in which we are not running OQC
/////

// We will do this as a special case to make it go as fast as possible.
// It is worth having a small structure to make the sort go fast.
// We will visit each clump once to make the structure.
// And no clump accesses during the sort.
// Yet, clump accesses during the search for dups will be extremely rare.

typedef struct dupArrayElem DupArrayElem_t;
struct dupArrayElem
{
    Clump_t * clump;
    ROFF      SRO;
    int       score;
};

int compareDupArrayElems(const void * arg1, const void * arg2)
{
    DupArrayElem_t * elem1 = (DupArrayElem_t *)arg1;
    DupArrayElem_t * elem2 = (DupArrayElem_t *)arg2;
    // It's dangerous to subtract ROFFs, so do it this way.
    // Also, sort by decreasing score, so that the higher score of two dups will be kept.
    if      (elem1->SRO > elem2->SRO) return 1;
    else if (elem1->SRO < elem2->SRO) return -1;
    else return elem2->score - elem1->score;
}

inline BOOL clumpsAreDups(Clump_t * clump1, Clump_t * clump2)
{
    // This is intended to be cheapest to calculate.
    // It is a combination of starting with the least likely tests, to exit early
    // And cost to calculate
    return (clumpSRO(clump1) == clumpSRO(clump2) && clumpSQO(clump1) == clumpSQO(clump2) &&
            clumpEQO(clump1) == clumpEQO(clump2) && clumpERO(clump1) == clumpERO(clump2) &&
            isReversed(clump1) == isReversed(clump2));
}

void postFilterRemoveDups(AlignmentArgs_t * AAs, QueryState_t * QS)
{
    // First build the array of clump elements.
    int nodeCount = QS->clumpCount;
    if (nodeCount < 2) return;
    DupArrayElem_t * dupArray = (DupArrayElem_t *)malloc((nodeCount)*sizeof(DupArrayElem_t));
    int curNodeCount = 0;
    for (Clump_t * clump = QS->clumps; clump != NULL; clump = clump->next)
    {
        DupArrayElem_t * curElem = &dupArray[curNodeCount];
        curElem->clump = clump;
        curElem->SRO = clumpSRO(clump);
        curElem->score = clump->totScore;
        curNodeCount += 1;
    }

    qsort(dupArray, nodeCount, sizeof(DupArrayElem_t), compareDupArrayElems);

    // We will have to rebuild the clump list as we go.
    Clump_t * newClumpList = NULL;
    int newClumpCount = 0;
    for (int i=0; i<nodeCount; i++)
    {
        DupArrayElem_t * elem1 = &dupArray[i];
        Clump_t * clump1 = elem1->clump;
        if (clump1 == NULL) continue;
        for (int j=i+1; j<nodeCount; j++)
        {
            DupArrayElem_t * elem2 = &dupArray[j];
            if (elem1->SRO < elem2->SRO) break;
            Clump_t * clump2 = elem2->clump;
            if (clump2 == NULL) continue;
            if (clumpsAreDups(clump1, clump2))
            {
                disposeClump(clump2, QS);
                elem2->clump = NULL;
            }
        }
        newClumpCount += 1;
        newClumpList = pushClump(clump1, newClumpList);
    }

    // Set the filtered clump list back into the query state.
    QS->clumps = newClumpList;
    QS->clumpCount = newClumpCount;

    free(dupArray);
}
