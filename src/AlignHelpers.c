/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    AlignHelpers.c     Includes the code that add fragments to form clumps, align clumps,
                                   and split clumps if they are not proper local alignments.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include "Math.h"
#include "SW.inl"
#include "FragsClumps.c"

//////////////////////////////////
//  Alignment Operations.
//////////////////////////////////

//////////////////////////////////
// Form clumps from fragments.
//////////////////////////////////

// This uses several heuristics to try to build clumps that will lead to the best alignments.
// We first "believe" what the graph algorithm finds as the best set of fragments, and add them to the clump.
// While doing this, we are careful to trim fragments so that they do not overlap each other on either the
//     query or the reference.
// We also keep them "sorted" by SQO by inserting them at the front of the list each time.
// This works because the best path from the graph algorithm starts at the end node, and goes backwards.

// We then take a second pass over the clump, and eliminate most fragments that have been so chopped,
//     but only if the matches they contain can be found again during DP.
// That is, if they fall within the bandwidth for a SW operation.
// Doing this as a post pass has MUCH better results than a previous attempt to do this filtering
//     while adding the fragments.
// This is largely due to the fact that as a post pass we now know all of the diagonal differences
//     needed to make the decisions.

// Create SFragment for the frag and add to the beginning of the SFragment list.
void addFragment(Clump_t * clump, Fragment_t * frag, QueryState_t * QS)
{
    // Let's keep track of the number of matching bases for use later.
    clump->matchedBases += fragMatchCount(frag);
    // Now add the fragment to the clump.
    SFragment_t * SFrag = makeSFragment(QS);
    copyFragment(frag, getFragPtr(SFrag));
    addSFragmentToFront(SFragListPtr(clump), SFrag);
}

// The next two functions go hand in hand, first adding, then removing fragments using certain heuristics as described above.
// Since the best path of fragments is traversed from end to beginning, fragments will always go in the beginning of the list.
void insertFragment(Clump_t * clump, Fragment_t * frag1, QueryState_t * QS)
{
    // If the clump is empty, add the fragment unconditionally and return.
    if (isEmpty(clump))
    {
        addFragment(clump, frag1, QS);
        return;
    }

    SFragment_t * nextSFrag = getFirstSFrag(SFragListPtr(clump));
    Fragment_t * frag2 = getFragPtr(nextSFrag);
    int maxOverlap = calcMaxOverlap(frag1, frag2);
    if (maxOverlap > 0)
    {
        int frag1len = fragQueryLen(frag1);
        int frag2len = fragQueryLen(frag2);

        // By default, we will chop frag2 unless this is the second frag being added, then chop frag1.
        // This will tend to keep the first and last frags from being chopped.
        // But if one of the frags is shorter than the other, we will chop it instead.
        BOOL chopFrag1;
        if (frag1len != frag2len) chopFrag1 = (frag1len < frag2len);
        else                      chopFrag1 = (nextSFrag->next == NULL);

        if (chopFrag1)
            subLenFromBack(frag1, maxOverlap);
        else
            subLenFromFront(frag2, maxOverlap);
    }
    addFragment(clump, frag1, QS);
}

void cleanUpClump(Clump_t * clump, QueryState_t * QS)
{
    AlignmentArgs_t * AAs = QS->AAs;
    // Go out until the next frag >= wordLen or last frag.
    // Then blow away everything in between that is < WL and within the bandwidth.
    // This eliminates sub-optimal alignments caused by sub-optimal off diagonal fragments,
    //     while still ensuring that the matches in the fragment we are deleting will be picked up during SW
    //     if they contribute to the best score.
    // But we need to have one other precaution.
    // If the two frags >= WL are more than maxGap apart, we need to leave at least
    //    one frag in between them in order to maintain the constraint that the alignment
    //    will not have a single gap bigger than maxGap.
    // Since this happens rarely, we will leave in all fragments in this case (almost always 1 anyway).
    SFragment_t * SFrag1 = getFirstSFrag(SFragListPtr(clump));
    SFragment_t * SFrag2 = getNextSFrag(SFragListPtr(clump));
    SFragment_t * SFrag3 = getNextSFrag(SFragListPtr(clump));
    while (SFrag2 != NULL && SFrag3 != NULL)
    {
        Fragment_t * frag2 = getFragPtr(SFrag2);
        if (fragQueryLen(frag2) < AAs->wordLen)
        {
            // Find the next full frag, or the last frag.
            SFragment_t * anchorSFrag = SFrag3;
            while (fragQueryLen(getFragPtr(anchorSFrag)) < AAs->wordLen && anchorSFrag->next != NULL)
                anchorSFrag = anchorSFrag->next;

            ROFF f1Diag = fragDiag(getFragPtr(SFrag1));
            ROFF anchorDiag = fragDiag(getFragPtr(anchorSFrag));
            // Make sure the two anchors are within maxGap of each other.
            if (absDiagDiff(f1Diag, anchorDiag) <= AAs->maxGap)
            {
                SFragment_t * delSFrag = SFrag2;
                while (delSFrag != anchorSFrag)
                {
                    // Need to keep a pointer to the next one in case this one is deleted.
                    SFragment_t * delNextSFrag = delSFrag->next;
                    ROFF delDiag = fragDiag(getFragPtr(delSFrag));
                    // If the middle frag has a middle diagonal, or if it is above or below but within the bandWidth...
                    if (!((delDiag < f1Diag && delDiag < anchorDiag) || (delDiag > f1Diag && delDiag > anchorDiag)) ||
                        (MIN(absDiagDiff(f1Diag, delDiag), absDiagDiff(delDiag, anchorDiag)) <= AAs->bandWidth))
                    {
                        removeSFragment(SFragListPtr(clump), delSFrag);
                        disposeSFragment(delSFrag, QS);
                    }
                    delSFrag = delNextSFrag;
                }
            }
            // Fool the next iteration setup.
            SFrag1 = anchorSFrag;
            SFrag2 = anchorSFrag->next;
        }
        else
        {
            // Advance one frag for next iteration.
            SFrag1 = SFrag2;
            SFrag2 = SFrag3;
        }

        // Set up for next iteration.
        if (SFrag2 != NULL) SFrag3 = SFrag2->next;
    }

    // Handle the first and last fragments.
    // Unlike the above where being within the bandwidth ensures the matches will be found,
    //    here we have to worry about a desert causing XCutoff to miss the matches.
    // So, we are only going to remove frags is they are right next to the adjoining frag.
    // Also recall that extensions use twice the bandwidth as internal SW calls.
    // First the first frag.
    SFrag1 = getFirstSFrag(SFragListPtr(clump));
    Fragment_t * frag1 = getFragPtr(SFrag1);
    if (fragQueryLen(frag1) < AAs->wordLen)
    {
        Fragment_t * frag2 = getFragPtr(getNextSFrag(SFragListPtr(clump)));
        int qGap = calcGap(frag1->endQueryOff, frag2->startQueryOff);
        int rGap = calcGap(fragEndRefOff(frag1), frag2->startRefOff);
        if ((qGap == 0 && rGap <= 2 * AAs->bandWidth) ||
            (rGap == 0 && qGap <= 2 * AAs->bandWidth))
        {
            removeSFragment(SFragListPtr(clump), SFrag1);
            disposeSFragment(SFrag1, QS);
        }
    }

    // Now the last frag
    // We now have to worry about the 2 frag case and the first already being deleted.
    SFrag2 = getLastSFrag(SFragListPtr(clump));
    Fragment_t * frag2 = getFragPtr(SFrag2);
    if (fragQueryLen(frag2) < AAs->wordLen)
    {
        SFrag1 = getPrevSFrag(SFragListPtr(clump));
        if (SFrag1 == NULL) return;
        Fragment_t * frag1 = getFragPtr(SFrag1);
        int qGap = calcGap(frag1->endQueryOff, frag2->startQueryOff);
        int rGap = calcGap(fragEndRefOff(frag1), frag2->startRefOff);
        if ((qGap == 0 && rGap <= 2 * AAs->bandWidth) ||
            (rGap == 0 && qGap <= 2 * AAs->bandWidth))
        {
            removeSFragment(SFragListPtr(clump), SFrag2);
            disposeSFragment(SFrag2, QS);
        }
    }
}

//////////////////////////////////
// Complete the alignment of a clump
//////////////////////////////////

// A few forward signature definitions to allow the actual code to flow is logical order.
void collapseSFragments(Clump_t * clump, QueryState_t * QS);
int  splitClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS);

// This takes in a clump that hopefully only contains fragments we wish to keep in the alignment.
// Stitch them together to form the alignment, and calculate the score.
int alignClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS)
{
    if (isAligned(clump)) return 0;

#ifdef DEBUG
    fprintf(stderr, "I am in align clump.\n");
    fprintQueryId(QS, stderr);
    // Debug output of clump structure
    printDebugClump(clump, AAs, stderr);
#endif

    char * refPtr = AAs->basePtr;
    char * queryPtr = getClumpQueryCodeBuf(QS, clump);

    SFragment_t * SFrag1;
    SFragment_t * SFrag2;
    Fragment_t * frag2;
    Fragment_t * frag1;

    // Do perfect extensions of fragments towards each other.
    // We'll worry about the rest of the alignment later.
    SFrag1 = getFirstSFrag(SFragListPtr(clump));
    for (SFrag2 = getNextSFrag(SFragListPtr(clump)); SFrag2 != NULL; SFrag2 = getNextSFrag(SFragListPtr(clump)))
    {
        frag1 = getFragPtr(SFrag1);
        frag2 = getFragPtr(SFrag2);
        // Since we will now use DP to calculate the edit distance in the gap, only do perfect extensions.
        int gap = calcMinGap(frag1, frag2);
        gap -= extendFragmentBackwardToStopPerfectly(frag2, refPtr, queryPtr, gap);
        gap -= extendFragmentForwardToStopPerfectly(frag1, refPtr, queryPtr, gap);

        SFrag1 = SFrag2;
    }

    // Align each SFragment in turn.
    // Make an EditOpList and calculate score as well.
    for (SFrag1 = getFirstSFrag(SFragListPtr(clump)); SFrag1 != NULL; SFrag1 = getNextSFrag(SFragListPtr(clump)))
    {
        int qLen = fragQueryLen(getFragPtr(SFrag1));
        addEditOpToFront(getEOLPtr(SFrag1), MatchEditOp, qLen);
        SFrag1->score = AAs->MScore * qLen;
    }

    // Make SFragments for the remaining gaps, and align them.
    // The various alignment steps will provide an EditOpList for every SFragment.
    // Later we will combine these EditOpLists into one for total score calculation and output.
    for (SFrag1 = SFragListPtr(clump)->head; SFrag1 != NULL; SFrag1 = SFrag1->next)
    {
        if (SFrag1 == SFragListPtr(clump)->tail) break;
        SFrag2 = SFrag1->next;
        // Make an SFragment for any remaining gap between fragments.
        SFragment_t * newSFrag = makeAndAlignSFragmentToFillGap(SFrag1, SFrag2, QS, clump);

        if (newSFrag == NULL) continue;
        ASSERT(FALSE, "Adding new fragment for a gap.\n");
        addSFragmentAfter(SFragListPtr(clump), SFrag1, newSFrag);
    }

    collapseSFragments(clump, QS);
    extendClumpForwardReverse(clump, QS);
    setAligned(clump, TRUE);

#ifdef DEBUG
    printDebugClump(clump, AAs, stderr);
#endif

    return 1;
}

void collapseSFragments(Clump_t * clump, QueryState_t * QS)
{
    // Combine the SFragment EditOpLists into one.
    EditOpList_t * list = EOLPtr(clump);
    int tempTotalScore = 0;
    for (SFragment_t * SFrag1 = getFirstSFrag(SFragListPtr(clump)); SFrag1 != NULL; SFrag1 = getNextSFrag(SFragListPtr(clump)))
    {
        tempTotalScore += SFrag1->score;
        mergeEOLToBack(list, getEOLPtr(SFrag1));
    }
    // Combine the SFragments into one SFragment for the whole alignment!
    // Do this by extending the first SFragment to encompass the entire alignment length,
    //     then trim the tail of the SFragment list.
    // Note that from this point on, we will have no idea what the hit counts are.
    SFragment_t * SFrag0 = SFragListPtr(clump)->head;
    SFragment_t * SFragn = SFragListPtr(clump)->tail;
    Fragment_t * frag0 = getFragPtr(SFrag0);
    Fragment_t * fragn = getFragPtr(SFragn);
    frag0->endQueryOff = fragn->endQueryOff;
    fragSetEndRefOff(frag0, fragEndRefOff(fragn));
    SFrag0->score = tempTotalScore;
    // Reduce the SFragment list to just SFrag0.
    SFragmentList_t * SFList = SFragListPtr(clump);
    removeSFragment(SFList, SFrag0);
    emptySFragmentList(SFList, QS);
    addSFragmentToFront(SFList, SFrag0);
}

int scoreClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS)
{
    if (isScored(clump)) return 1;

#ifdef QUERYSTATS
    QS->alignCount += 1;
#endif

    // Now calculate the score and alignment length.
    // We can calc the score from the EditOpList very accurately.
    int AGS = 0;
    int maxAGS = 0;
    int matches = 0;
    int mismatches = 0;
    int inserts = 0;
    int deletes = 0;
    EditOpList_t * list = EOLPtr(clump);
    EditOp_t * lastItem = getLastEditOp(list);
    int alignedScore = SFragListPtr(clump)->head->score;

#ifdef NOTNOW
    // Believe it or not, this can happen, but it is rare.
    // The cases I have seen involve an indel > maxGap in length, but NOT for all such alignments!!
    // This is a result of finding two indels during DP which are each smaller than the maxGap,
    //   but then get merged in the backtrack code into one indel larger than maxGap.
    // In such cases, the original score will subtract two GOCs, but the merged/recalculated score only one.
    // See longer discussion in EverNote.
    int tempScore = EditOpList2AGS(list, AAs);
    if (alignedScore != tempScore)
    {
        fprintf(stderr, "Aligned score %d does not match recalculated score %d .\n", alignedScore, tempScore);
        // fprintQueryId(QS, stderr);
        // printDebugClump(clump, AAs, stderr);
    }
#endif

    forAllEditOpsInList(item, list)
    {
        EditOpCode op = item->opcode;
        int len = item->length;
        if      (op == MatchEditOp)   {matches    += len; AGS += (AAs->MScore * len);}
        else if (op == ReplaceEditOp) {mismatches += len; AGS -= (AAs->RCost  * len);}
        else if (op == InsertEditOp)  {inserts    += len; AGS -= (AAs->GOCost + AAs->GECost * len);}
        else if (op == DeleteEditOp)  {deletes    += len; AGS -= (AAs->GOCost + AAs->GECost * len);}
        if (AGS <= 0 || (AGS >= alignedScore && item != lastItem)) return splitClump(clump, AAs, QS);
        if (AGS > maxAGS) maxAGS = AGS;
    }
    if (matches >= AAs->minRawScore && maxAGS > AGS) return splitClump(clump, AAs, QS);
    if (matches < AAs->minRawScore) return 0;

    clump->matchedBases = matches;
    clump->mismatchedBases = mismatches;
    clump->gapBases = inserts + deletes;
    clump->totLength = matches + mismatches + inserts + deletes;
    clump->totScore = AGS;

    // Calculate the percentage and filter out unworthy alignments.
    double percent =  ((double)clump->matchedBases/clump->totLength);
    if (percent < AAs->minIdentity) return 0;

    // The clump has passed all tests.
    // Mark as successfully scored, and return;
    setScored(clump, TRUE);
    return 1;
}

//////////////////////////////////
// Handle the case in which we have found something that indicates we don't have a valid local alignment.
// We will break the alignment up into pieces each of which might be a valid local alignment.
// Not all portions will necessarily survive in isolation after split.
//////////////////////////////////

int splitClumpHelper(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS, int wSQO, int wEQO)
{
#ifdef DEBUG
    fprintf(stderr, "Entering splitClumpHelper, depth %d ", depth);
    fprintQueryId(QS, stderr);
    printDebugClump(clump, AAs, stderr);
    printEditOpList(EOLPtr(clump));
    fprintf(stderr, "\n");
#endif
    // Get the SFrag to store the offset and EOL as we find the surviving pieces.
    SFragment_t * curSFrag = SFragListPtr(clump)->head;
    Fragment_t * curFrag = getFragPtr(curSFrag);
    EditOpList_t * list = getEOLPtr(curSFrag);
    // At this point, the clump has the EOL, but the superFragment does not.
    // So, put the EOL back into the superFragment.
    mergeEOLToFront(list, EOLPtr(clump));

    // We will go through the EditOpList and find the maximum score.
    // Then we will go backwards and find the first zero.
    // Those endpoints will define the best alignment.
    // We will then have to process the remaining ends of the list.
    SUINT sQO = 0;
    SUINT eQO = 0;
    ROFF  sRO = 0;
    ROFF  eRO = 0;
    int matches = 0;
    int mismatches = 0;
    int inserts = 0;
    int deletes = 0;
    int AGS = 0;
    int totAGS = 0;
    int maxAGS = -10000;
    EditOp_t * maxItem = NULL;
    EditOp_t * minItem = NULL;
    forAllEditOpsInList(item, list)
    {
        EditOpCode op = item->opcode;
        int len = item->length;
        int newScore = 0;
        if      (op == MatchEditOp)   {matches    += len; newScore = (AAs->MScore * len);}
        else if (op == ReplaceEditOp) {mismatches += len; newScore = -(AAs->RCost  * len);}
        else if (op == InsertEditOp)  {inserts    += len; newScore = -(AAs->GOCost + AAs->GECost * len);}
        else if (op == DeleteEditOp)  {deletes    += len; newScore = -(AAs->GOCost + AAs->GECost * len);}
        totAGS += newScore;
        AGS += newScore;
        if (AGS < 0) AGS = 0;
        if (AGS > maxAGS)
        {
            maxAGS = AGS;
            maxItem = item;
            eQO = curFrag->startQueryOff + matches + mismatches + inserts - 1;
            eRO = curFrag->startRefOff + matches + mismatches + deletes - 1;
        }
    }

    // Now go backwards, looking for the first zero.
    // Note that this time we are subtracting the score from the max.
    AGS = maxAGS;
    matches = mismatches = inserts = deletes = 0;
    int maxMatch = 0;
    for (EditOp_t * item = maxItem; isValidEONode(item, list); item = getPrevEditOp(item, list))
    {
        EditOpCode op = item->opcode;
        int len = item->length;
        if      (op == MatchEditOp)   {matches    += len; AGS -= (AAs->MScore * len); if (len > maxMatch) maxMatch = len;}
        else if (op == ReplaceEditOp) {mismatches += len; AGS += (AAs->RCost  * len);}
        else if (op == InsertEditOp)  {inserts    += len; AGS += (AAs->GOCost + AAs->GECost * len);}
        else if (op == DeleteEditOp)  {deletes    += len; AGS += (AAs->GOCost + AAs->GECost * len);}
        // We have found the beginning of the "local" alignment.
        // Set the starting offsets of the frag we are building.
        if (AGS <= 0)
        {
            minItem = item;
            sQO = eQO - (matches + mismatches + inserts - 1);
            sRO = eRO - (matches + mismatches + deletes - 1);
            break;
        }
    }
    // Make sure we have at least one k-mer of the original seed size.
    if (maxMatch < AAs->wordLen) return 0;

    // Now we have the information to form the new head and/or tail SFragment(s).
    // But we must protect against the thread global EditOpList being reallocated between the two.
    // We can either first form both possible new clumps before recursing on either of them.
    // Or we can cheese out and just store the maxItem offset before the first recursive step, and restore below.
    EditOpPtr maxItemPtr = ptrFromNode(maxItem, list);

    int retval = 0;
    // First see if we need to form a head fragment.
    if (minItem != getFirstEditOp(list))
    {
        SFragment_t * newSFrag = makeSFragment(QS);
        EditOpList_t * newEOL = getEOLPtr(newSFrag);
        Clump_t * newClump = makeClump(QS);
        setReversed(newClump, isReversed(clump));
        addSFragmentToFront(SFragListPtr(newClump), newSFrag);
        mergeEOLToFront(newEOL, list);
        // NOTE!  This side effects the "list" to no longer include the head section.
        splitEditOpListBefore(newEOL, list, minItem);
        // Again make sure we have at least one k-mer of the original seed size.
        if (EditOpList2Maxmatch(newEOL, AAs->wordLen))
        {
            Fragment_t * newFrag = getFragPtr(newSFrag);
            newFrag->startQueryOff = curFrag->startQueryOff;
            newFrag->endQueryOff = sQO - 1;
            newFrag->startRefOff = curFrag->startRefOff;
            fragSetEndRefOff(newFrag, sRO - 1);
            retval += splitClumpHelper(newClump, AAs, QS, wSQO, wEQO);
        }
        // Check if this portion of the original made it through scoring.
        if (isScored(newClump))
        {
            // Add the clump into the clump list.
            setSplit(newClump, TRUE);
            setAligned(newClump, TRUE);
            addClump(QS, newClump);
            setReversed(newClump, isReversed(clump));
            // insertClumpAfter(QS, clump, newClump);
        }
        else
            disposeClump(newClump, QS);
    }

    // Restore the maxItem from the saved offset.
    maxItem = nodeFromPtr(maxItemPtr, list);
    // Now see if we need a tail fragment.
    if (maxItem != getLastEditOp(list))
    {
        SFragment_t * newSFrag = makeSFragment(QS);
        EditOpList_t * newEOL = getEOLPtr(newSFrag);
        Clump_t * newClump = makeClump(QS);
        setReversed(newClump, isReversed(clump));
        addSFragmentToFront(SFragListPtr(newClump), newSFrag);
        // NOTE!  This side effects the "list" to no longer include the tail section.
        splitEditOpListAfter(list, newEOL, maxItem);
        // Again make sure we have at least one k-mer of the original seed size.
        if (EditOpList2Maxmatch(newEOL, AAs->wordLen))
        {
            Fragment_t * newFrag = getFragPtr(newSFrag);
            newFrag->startQueryOff = eQO + 1;
            newFrag->endQueryOff = curFrag->endQueryOff;
            newFrag->startRefOff = eRO + 1;
            fragSetEndRefOff(newFrag, fragEndRefOff(curFrag));
            retval += splitClumpHelper(newClump, AAs, QS, wSQO, wEQO);
        }
        // Check if this portion of the original made it through scoring.
        if (isScored(newClump))
        {
            // Add the clump into the clump list.
            setSplit(newClump, TRUE);
            setAligned(newClump, TRUE);
            addClump(QS, newClump);
            setReversed(newClump, isReversed(clump));
            // insertClumpAfter(QS, clump, newClump);
        }
        else
            disposeClump(newClump, QS);
    }

    // Now we can do the actual good frag!
    curFrag->startQueryOff = sQO;
    curFrag->endQueryOff = eQO;
    curFrag->startRefOff = sRO;
    fragSetEndRefOff(curFrag, eRO);
    curSFrag->score = maxAGS;

    // We need to also put the list back into the clump;
    mergeEOLToFront(EOLPtr(clump), list);

    // We do careful extensions to avoid recreating the initial alignment and going into an infinite regress.
    // We start, as with initial alignment, with making perfect extensions in both directions.
    // This gives a higher starting score for the careful extensions.
    // Only extend if we are not at the beginning or end of the original alignment.
    // In such cases, the extensions have already been done in that direction.
    // Although perhaps with a lower score for careful extensions.
    BOOL goBack = (sQO != wSQO);
    BOOL goForw = (eQO != wEQO);
    extendClumpForwardReverseCarefully(clump, QS, goBack, goForw);

    // Mark the clump as split, and again score it to make sure it passes all filters.
    setSplit(clump, TRUE);
    retval += scoreClump(clump, AAs, QS);
    return retval;
}

// Just a wrapper function around the main workhorse above.
// Also allows us to pass in the original boundaries.
int splitClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS)
{
#ifdef DEBUG
    fprintf(stderr, "Splitting alignment to form one or more valid local alignments.");
    fprintQueryId(QS, stderr);
    printDebugClump(clump, AAs, stderr);
    printEditOpList(EOLPtr(clump));
#endif
    // Get the SQO and EQO of the entire clump plus its score to hand into the helper.
    SFragment_t * curSFrag = SFragListPtr(clump)->head;
    Fragment_t * curFrag = getFragPtr(curSFrag);
    int retval = splitClumpHelper(clump, AAs, QS, curFrag->startQueryOff, curFrag->endQueryOff);
#ifdef DEBUG
    fprintf(stderr, "Main clump to print is as follows.\n");
    printDebugClump(clump, AAs, stderr);
    fprintf(stderr, "\n");
#endif
    return retval;
}
