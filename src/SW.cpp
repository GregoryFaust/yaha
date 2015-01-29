/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    SW.c            Contains Dynamic Programming (Smith Waterman) algorithms for DNA only.

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
// #define TIMING
#include "Timing.inl"

// TODO: these timers can't be global in a threaded environment.
#ifdef TIMING
setupTimers();
defineTimer(junkTimer);
defineTimer(allocTimer);
defineTimer(DPINITTimer);
defineTimer(FullSWTimer);
defineTimer(ExtensionSWTimer);
defineTimer(FullBacktrackSWArraysTimer);
defineTimer(ExtensionBacktrackSWArraysTimer);
defineTimer(decompressTimer);
#endif

//////////////////////////////////
// EditOps and EditOpLists
/////////////////////////////////
// Just as with clumps and sFragments, we wish to use a slab allocator for editops.
// However, because we create so many more edit ops than any other structure,
//      and they are very small, we use a special size optimization for them.
// Instead of 8-byte forward and backward pointers, we use 4-byte for and back offsets.
// For this to work, we will need to keep ALL editops in a single array.
// Heuristically, each time the array is too small for the next neeeded editop,
//      we double the size of the array.
// (Alternatively, we could have each list be contained all in one array, and have the
//  array ptr in the list structure, but we choose not to addd that complication)
/////////////////////////////////

// If the EOLArray is full, we need to make a larger array, requiring reallocation.
extern "C" {void makeMoreEditOps(QueryState_t * QS);}

void makeMoreEditOps(QueryState_t * QS)
{
    addToTimer(junkTimer);

    // If this is the first time, alloc the array.
    // Otherwise, we will have to realloc.
    int startOffset;
    if (QS->EditOpArray == NULL)
    {
        QS->EditOpArray = (EditOp_t *)malloc(sizeof(EditOp_t) * QS->EditOpArrayMaxElements);
        // Never use the first element in the array so that we can use 0 as NULL ptr value.
        startOffset = 1;
    }
    else
    {
        // Check if we can make any more elements.
        if (QS->EditOpArrayMaxElements & 0x80000000) fatalError("Maximum number of EditOperation objects exceeded.");
        // Double the size of the array.
        startOffset = QS->EditOpArrayMaxElements;
        QS->EditOpArrayMaxElements *= 2;
        EditOp_t * EOLArrayt = (EditOp_t *)realloc(QS->EditOpArray, sizeof(EditOp_t) * QS->EditOpArrayMaxElements);
        // fprintf(stderr, "Reallocating EOLArray%s.  New size=%x\n", (QS->EditOpArray == EOLArrayt) ? "" : " to a new location", QS->EditOpArrayMaxElements);
        QS->EditOpArray = EOLArrayt;
    }
    // Check for success.
    if (QS->EditOpArray == NULL) fatalError("Insufficeint memory available to allocate (more) EditOperation objects.");
    // Now initialize the new area of the array.
    QS->EditOpFreeList = startOffset;
    for (int i=startOffset; i<QS->EditOpArrayMaxElements; i++)
    {
        QS->EditOpArray[i].next = i+1;
    }
    QS->EditOpArray[QS->EditOpArrayMaxElements-1].next = EditOpPtrNULL;
    addToTimer(allocTimer);
    fprintTimer(stderr, "Time to alloc EOLArray: ", allocTimer);
    resetTimer(allocTimer);
}

EditOpPtr allocEditOp (EditOpCode code, SUINT len, EditOpList_t * list)
{
    QueryState_t * QS = getEOLQSPtr(list);
    if (QS->EditOpFreeList == EditOpPtrNULL) makeMoreEditOps(QS);
    EditOpPtr ptr = QS->EditOpFreeList;
    EditOp_t * node = nodeFromPtr(ptr, list);
    QS->EditOpFreeList = node->next;
    node->next = EditOpPtrNULL;
    node->length = len;
    node->opcode = code;
    return ptr;
}

void freeEditOp(EditOpPtr ptr, EditOpList_t * list)
{
    QueryState_t * QS = getEOLQSPtr(list);
    EditOp_t * node = nodeFromPtr(ptr, list);
    node->next = QS->EditOpFreeList;
    QS->EditOpFreeList = ptr;
}

void emptyEditOpList(EditOpList_t * list)
{
    // If list is not empty, put onto free list.
    if (!(EOLisEmpty(list)))
    {
        QueryState_t * QS = getEOLQSPtr(list);
        EditOp_t * tailNode = nodeFromPtr(list->tail, list);
        tailNode->next = QS->EditOpFreeList;
        QS->EditOpFreeList = list->head;
    }
    list->head = list->tail = EditOpPtrNULL;
}

EditOpList_t * makeEditOpList(QueryState_t * QS)
{
    EditOpList_t * list = (EditOpList_t *)malloc(sizeof(EditOpList_t));
    initEditOpList(list, QS);
    return list;
}

void disposeEditOpList(EditOpList_t * list)
{
    emptyEditOpList(list);
    free(list);
}

void addEditOpToList(EditOpList_t * list, EditOpPtr prev, EditOpPtr next, EditOpPtr add)
{
    EditOp_t * newNode = nodeFromPtr(add, list);
    newNode->prev = prev;
    newNode->next = next;
    if (prev == EditOpPtrNULL) list->head = add;
    else nodeFromPtr(prev, list)->next = add;
    if (next == EditOpPtrNULL) list->tail = add;
    else nodeFromPtr(next, list)->prev = add;
}

void mergeEOLToFront (EditOpList_t * target, EditOpList_t * source)
{
    if (source->head == EditOpPtrNULL) return;
    if (target->head == EditOpPtrNULL)
    {
        target->head = source->head;
        target->tail = source->tail;
    }
    else
    {
        // Check for matching opcodes at the knit point before merging.
        EditOp_t * sourceTailNode = nodeFromPtr(source->tail, source);
        EditOp_t * targetHeadNode = nodeFromPtr(target->head, source);
        if (sourceTailNode->opcode == targetHeadNode->opcode)
        {
#ifdef NOTNOW
            // This code checks for the unusual case alluded to in AlignHelpers.c
            // The DP code can produce two ajoining indels of the same type.
            // Neither one is longer than maxGap, but when combined they are.
            // It is not clear how to get rid of this issue.
            // This code is left here for future investigations on the matter.
            if (sourceTailNode->opcode == InsertEditOp || sourceTailNode->opcode == DeleteEditOp)
            {
                int newLen = sourceTailNode->length + targetHeadNode->length;
                if (newLen > 50)
                {
                    fprintf(stderr, "EOLToFront: Mashing two Indels into an over long one.  %d + %d = %d \n", sourceTailNode->length, targetHeadNode->length, newLen);
                    fprintf(stderr, "Source ");
                    printEditOpList(source);
                    fprintf(stderr, "\nTarget ");
                    printEditOpList(target);
                    fprintf(stderr, "\n");
                }
            }
#endif
            // Do the merge
            sourceTailNode->length += targetHeadNode->length;
            EditOpPtr newHead = targetHeadNode->next;
            freeEditOp(target->head, target);
            // If we are deleting the only node in the target, the target tail must be set!
            if (newHead == EditOpPtrNULL)
            {
                target->tail = source->tail;
                targetHeadNode = NULL;
            }
            else targetHeadNode = nodeFromPtr(newHead, source);
            target->head = newHead;
        }
        // The normal case, or post merge.
        sourceTailNode->next = target->head;
        if (targetHeadNode != NULL) targetHeadNode->prev = source->tail;
        target->head = source->head;
    }
    source->head = source->tail = EditOpPtrNULL;
}

void mergeEOLToBack (EditOpList_t * target, EditOpList_t * source)
{
    if (source->head == EditOpPtrNULL) return;
    if (target->head == EditOpPtrNULL)
    {
        target->head = source->head;
        target->tail = source->tail;
    }
    else
    {
        // Check for matching opcodes at the knit point before merging.
        EditOp_t * targetTailNode = nodeFromPtr(target->tail, source);
        EditOp_t * sourceHeadNode = nodeFromPtr(source->head, source);
        if (targetTailNode->opcode == sourceHeadNode->opcode)
        {
#ifdef NOTNOW
            // This code checks for the unusual case eluded to in AlignHelpers.c
            // The DP code can produce two ajoining indels of the same type.
            // Neither one is longer than maxGap, but when combined they are.
            // It is not clear how to get rid of this issue.
            // This code is left here for future investigations on the matter.
            if (targetTailNode->opcode == InsertEditOp || targetTailNode->opcode == DeleteEditOp)
            {
                int newLen = targetTailNode->length + sourceHeadNode->length;
                if (newLen > 50)
                {
                    fprintf(stderr, "EOLToBack: Mashing two Indels into an over long one.  %d + %d = %d \n", targetTailNode->length, sourceHeadNode->length, newLen);
                    fprintf(stderr, "Target ");
                    printEditOpList(target);
                    fprintf(stderr, "\nSource ");
                    printEditOpList(source);
                    fprintf(stderr, "\n");
                }
            }
#endif
            // Do the merge
            targetTailNode->length += sourceHeadNode->length;
            EditOpPtr newHead = sourceHeadNode->next;
            freeEditOp(source->head, source);
            source->head = newHead;
            // If the source is now empty, we don't want to set the target tail to the free list!
            if (newHead == EditOpPtrNULL)
            {
                source->tail = target->tail;
                sourceHeadNode = NULL;
            }
            else sourceHeadNode = nodeFromPtr(newHead, source);
        }
        // The normal case, or post merge.
        targetTailNode->next = source->head;
        if (sourceHeadNode != NULL) sourceHeadNode->prev = target->tail;
        target->tail = source->tail;
    }
    source->head = source->tail = EditOpPtrNULL;
}

void splitEditOpListBefore (EditOpList_t * source2head, EditOpList_t * newTail, EditOp_t * node)
{
    EditOpPtr prev = node->prev;
    newTail->tail = source2head->tail;
    newTail->head = ptrFromNode(node, source2head);
    source2head->tail = prev;
    if (prev == EditOpPtrNULL) source2head->head = EditOpPtrNULL;
    else nodeFromPtr(prev, source2head)->next = EditOpPtrNULL;
    node->prev = EditOpPtrNULL;
}

void splitEditOpListAfter (EditOpList_t * source2head, EditOpList_t * newTail, EditOp_t * node)
{
    EditOpPtr next = node->next;
    newTail->tail = source2head->tail;
    source2head->tail = ptrFromNode(node, newTail);
    newTail->head = next;
    if (next == EditOpPtrNULL) newTail->tail = EditOpPtrNULL;
    else nodeFromPtr(next, newTail)->prev = EditOpPtrNULL;
    node->next = EditOpPtrNULL;
}

// This is a debugging aid.
void validateEditOpList(EditOpList_t * list)
{
    if (list->head == EditOpPtrNULL)
    {
        if (list->tail != EditOpPtrNULL) fprintf(stderr, "Error.  Tail ptr is non-null for null head ptr.\n");
        return;
    }
    if (list->tail == EditOpPtrNULL)
    {
        fprintf(stderr, "Error.  Tail ptr is null for non-null head ptr.\n");
        return;
    }
    if (nodeFromPtr(list->tail, list)->next != EditOpPtrNULL)
    {
        fprintf(stderr, "Error.  Tail node has non-null next ptr.\n");
        return;
    }
    EditOp_t * tail = nodeFromPtr(list->tail, list);
    BOOL foundtail = FALSE;
    forAllEditOpsInList(item, list)
    {
        if (item == tail) foundtail = TRUE;
    }
    if (foundtail == FALSE) fprintf(stderr, "Node pointed to by tail ptr not in list.\n");
}

// This is a debugging aid.
void printEditOpList (EditOpList_t * list)
{
    int count = 0;
    forAllEditOpsInList(item, list)
    {
        count += 1;
        fprintf(stderr, "%d%c", item->length, EditOpCode2char(item->opcode));
    }
}

/////
//     Start of SW code supporting Affine Gap Scoring
/////

// Here is the code to perform a modified version of SW needed for our algorithms.
// SW is defined to do "local" alignments, which are "unanchored" at either end.
// Meaning the start and end offsets for the best alignments on both the query and reference are unknown when we start.
// In use in this aligner, we will need to do "global" alignments in which both the start and end offsets are known.
// This is used to align between known runs of seed matches.
// In addition, we will use this code to do extensions.
// In which case, the start of query and reference are known, but the ending offsets for neither are known.

// The use of affine gap scoring means that inserts and deletes can "tunnel" past collections of matches and replaces
//     if the gap extension cost is less than a replacement cost.
// In particular, we need to store the so-called E and F arrays in additional to the usual V array.
// The good news is that we only need to store small parts of these arrays.
// We will keep track of the actual alignment with a "pointer" array,
//     which is stored as two parallel 2D arrays that store the edit operation and the lengths for indels.

// In addition, we will support banded SW as a very useful time/accuracy trade-off.
// The easiest way to imagine banded SW is to envision some space along the diagonal of the DP array(s)
//     that are calculated, while the rest of the values are not.
// However, while this saves time, it still requires the array(s) be full sized.
// This is impractical for very long queries.
// Therefore, instead, we will make the banded SW essentially use a vertical column that is the width of the sum of
//     the left and right bandwidths and the indel size (diagonal difference, if any), and height equal to the length of the DP run (along the query).
// It is possible that this will also improve cache performance, but this is yet to be seen.

// N.B.  To be correct, the integer lengths of scores must be at least as long as the length of queries.
// For short queries, we could therefore use 16 bit ints.
// However, emperically we have determined that full ints seem to run faster on X86 architecture.

// This is something less than -infinity so that any small - values added to it will not cause roll over to + value.
#define DPWorstScore (-(0x7fffff00))

// To save allocation and dealocation costs, we will allocate these once and keep throughout the run.
// Since the sizes (array layout) change with each banded run, we can't preassign values to the first row and column.
void DPInit(QueryState_t * QS)
{
    startTime();
    // Since we will use banded SW for any long alignments, we can VASTLY reduce the array sizes.
    // Need 4 times the bandWidth, as we double the bandWidth for extensions.
    // Add extra width of 2 for the overflow needed at begin and end of banded rows.
    AlignmentArgs_t * AAs = QS->AAs;
    int arrWidth = 4*AAs->bandWidth + AAs->maxGap + 3;
    int arrHeight = QS->maxQueryLength;
    QS->EOArray = (EditOpCode *) malloc (sizeof(EditOpCode) * arrWidth * arrHeight);
    QS->IDCount = (DPScore *) malloc (sizeof(DPScore) * arrWidth * arrHeight);

    // This is used as a temporary storage for the reference genome region of interest.
    // Having the temp helps with extensions because each base is examined multiple times,
    // Therefore, it is faster to avoid doing all of the calculations for the decompress each time.
    QS->refTemp = (char *)  malloc (arrHeight);

    QS->PVRowSave = (DPScore *) malloc (sizeof(DPScore) * arrWidth);
    QS->PFRowSave = (DPScore *) malloc (sizeof(DPScore) * arrWidth);
    QS->PIRowSave = (DPScore *) malloc (sizeof(DPScore) * arrWidth);

    // We use the unknown edit op to terminate the backtrace while finding the EditOpList for an alignment.
    QS->EOArray[0] = UnknownEditOp;
    QS->IDCount[0] = 0;

    // This is a hack for banded.
    // We will access these with -1 index for some rows of banded, and want DPWorstScore there.
    QS->PVRowSave[0] = DPWorstScore;
    QS->PVRowSave += 1;
    addToTimer(DPINITTimer);
}

void DPCleanUp(QueryState_t * QS)
{
    // This unwind of the hack for banded so free will work!
    QS->PVRowSave -= 1;

    // if (EOLArray != NULL) free(EOLArray);
    if (QS->EOArray != NULL) free(QS->EOArray);
    if (QS->IDCount != NULL) free(QS->IDCount);
    if (QS->PVRowSave != NULL) free(QS->PVRowSave);
    // if (VRowSave != NULL) free(VRowSave);
    if (QS->PFRowSave != NULL) free(QS->PFRowSave);
    if (QS->PIRowSave != NULL) free(QS->PIRowSave);
    if (QS->refTemp != NULL) free(QS->refTemp);

#ifdef TIMING
    // Print the timers here.
    endTime(junkTimer);
    fprintTimer(stderr, "Time to decompress Reference: ", decompressTimer);
    fprintTimer(stderr, "Time to perform Full SW: ", FullSWTimer);
    fprintTimer(stderr, "Time to perform SW for extensions: ", ExtensionSWTimer);
    fprintTimer(stderr, "Total time to perform SW: ", FullSWTimer+ExtensionSWTimer);
    fprintTimer(stderr, "Time to backtrack SW arrays for Full SW: ", FullBacktrackSWArraysTimer);
    fprintTimer(stderr, "Time to backtrack SW arrays for Extension SW: ", ExtensionBacktrackSWArraysTimer);
    fprintTimer(stderr, "Total time to backtrack SW arrays: ", FullBacktrackSWArraysTimer+ExtensionBacktrackSWArraysTimer);
    fprintTimer(stderr, "Total time in all above activities: ", FullBacktrackSWArraysTimer+ExtensionBacktrackSWArraysTimer+FullSWTimer+ExtensionSWTimer+decompressTimer);
#endif
}

// We need various versions of the AGS DP code to handle different aspects of an alignment.
// For aligning gaps between seeds, we will need to do a "global alignment" of the query
//     and reference regions that need to be aligned.
// For aligning extensions at the beginning and end of an alignment, we will need
//     to anchor just one end of the alignment, but let the termination point
//     be discovered by finding the max score and working back (ala SW).
// NOTE  We will keep track of the max value cell as we calculate the arrays to avoid
//     having to store all of the V array, and also to avoid going back through the
//     entire array to find the max value.
// At the moment, we have no case in which we want to have BOTH ends of the alignment
//     float (the so called "local aligment that SW does).
//     To add this, we would have to change both the array calculation and backtracking phases.
//     In particular, we would have to dis-allow negative numbers from the V array.
//     And also change backtracking to go as far as the first zero value, not all the way to (1,1).
//     This might make it mandatory to keep the entire V array during the forward phase.
// Also, to avoid the overhead of getting the right code from the reference genome over and over,
//     we will first copy the portion of interest into a char array.

// To maintain sanity in the code, we will use C++ templates for the various versions of the algorithm
//     in the same code base.
// We will define separate C function interfaces so we can just have this file be compiled via C++
//     while the remainder of the code stays in C.

// Decompress the portion of the compressed reference of interest to use one byte per code.
template <BOOL reverse>
void decompressRef(char * refPtr, char * decompressed, ROFF start, SUINT rLen)
{
    addToTimer(junkTimer);
    for (int i=0; i<rLen; i++)
    {
        if (reverse)
            decompressed[i] = getFrom4Code(refPtr, start-i);
        else
            decompressed[i] = getFrom4Code(refPtr, start+i);
    }
    addToTimer(decompressTimer);
}

// Some forward declarations of things later in the file.
template <BOOL banded, BOOL extension, BOOL reverse, BOOL XCutoff>
int findAffineGapScore(QueryState_t * QS, char * qStr, SUINT qLen, char * rStr, SUINT rLen, EditOpList_t * list, SUINT * addedQLen = NULL, SUINT * addedRLen = NULL);

int findAGSAlignment(QueryState_t * QS, ROFF rOff, SUINT rLen, char * qStr, SUINT qOff, SUINT qLen, EditOpList_t * list)
{
    decompressRef<FALSE>(QS->AAs->basePtr, QS->refTemp, rOff, rLen);
    // Since this is NOT banded, the bandwidth doesn't matter.
    int retval = findAffineGapScore<FALSE, FALSE, FALSE, FALSE>(QS, qStr+qOff, qLen, QS->refTemp, rLen, list);
    return retval;
}

int findAGSAlignmentBanded(QueryState_t * QS, ROFF rOff, SUINT rLen, char * qStr, SUINT qOff, SUINT qLen, EditOpList_t * list)
{
    decompressRef<FALSE>(QS->AAs->basePtr, QS->refTemp, rOff, rLen);
    int retval = findAffineGapScore<TRUE, FALSE, FALSE, FALSE>(QS, qStr+qOff, qLen, QS->refTemp, rLen, list);
    return retval;
}

// For extensions, the edit op list will be formed in the wrong order.
// Therefore, we will add to a temp list, and then append as appropriate.
template <BOOL reverse>
inline int findAGSExtension(QueryState_t * QS, ROFF rOff, char * qStr, SUINT qOff, SUINT qLenArg,
                            EditOpList_t * list, SUINT * addedQLen, SUINT * addedRLen)
{
    // Because of the below zero calcs, we need a full int for querylen.
    // TODO: Obviously once qLen comes in as an int, we can do away with this.
    int qLen = qLenArg;
#ifdef DEBUG
    if (reverse)
        fprintf(stderr, "Extending Backward from QO=%u, RO=%u, QL=%u\n", qOff, rOff, qLen);
    else
        fprintf(stderr, "Extending Forward from QO=%u, RO=%u, maxRO=%u, QL=%u\n", qOff, rOff, AAs->maxROff, qLen);
#endif
    *addedQLen = 0;
    *addedRLen = 0;
    if (qLen <= 0) return 0;
    AlignmentArgs_t * AAs = QS->AAs;
    int bandwidth = 2 * AAs->bandWidth;
    ROFF rLen = qLen + bandwidth;
    // Make sure we are not going to run off of the reference in either direction!!
    if (reverse && rLen > rOff)
    {
#ifndef COMPILE_USER_MODE
        fprintf(stderr, "In Backward Extension, trimming lengths to avoid running off beginning of reference.\n");
#endif
        rLen = rOff + 1;
        qLen = rLen - bandwidth;
        if (qLen <= 0) return 0;
    }
    if (!reverse && (rOff + rLen) > AAs->maxROff)
    {
#ifndef COMPILE_USER_MODE
        fprintf(stderr, "In Foreward Extension, trimming lengths to avoid running off end of reference.\n");
#endif
        rLen = AAs->maxROff - rOff;
        qLen = rLen - bandwidth;
        if (qLen <= 0) return 0;
    }

    // Perform the DP.
    decompressRef<reverse>(AAs->basePtr, QS->refTemp, rOff, rLen);
    EditOpList_t stackList;
    EditOpList_t * tempList = &stackList;
    int AGS = findAffineGapScore<TRUE, TRUE, reverse, TRUE>(QS, qStr+qOff, qLen, QS->refTemp, rLen, tempList, addedQLen, addedRLen);

    // If the extension has a negative or zero score, we are done.
    if (AGS <= 0) return 0;

    // Merge the editoplist to the list for the alignment.
    if (reverse)
        mergeEOLToFront(list, tempList);
    else
        mergeEOLToBack(list, tempList);
    return AGS;
}

// These are just convenience wrappers.
// They could be inlined, but what they are calling is so expensive anyway, why bother.
int findAGSBackwardExtension(QueryState_t * QS,  ROFF rOff, char * qStr, SUINT qOff, SUINT qLen,
                             EditOpList_t * list, SUINT * addedQLen, SUINT * addedRLen)
{
    return findAGSExtension<TRUE>(QS, rOff, qStr, qOff, qLen, list, addedQLen, addedRLen);
}

int findAGSForwardExtension(QueryState_t * QS, ROFF rOff, char * qStr, SUINT qOff, SUINT qLen,
                            EditOpList_t * list, SUINT * addedQLen, SUINT * addedRLen)
{
    return findAGSExtension<FALSE>(QS, rOff, qStr, qOff, qLen, list, addedQLen, addedRLen);
}

// #define TEMPDEBUG
// This alternate version will check that the extension does not make the entire alignment score go negative.
// This can't happen for backward alignments, as the backward alignment stops at a maximum.
// We could recalculate the score from the edit op list, but we are already doing such things too often.
int findAGSForwardExtensionCarefully(QueryState_t * QS, ROFF rOff, char * qStr, SUINT qOff, SUINT qLen,
                                     EditOpList_t * list, int score, SUINT * addedQLen, SUINT * addedRLen)
{
#ifdef TEMPDEBUG
    fprintf(stderr, "Entering extend forward carefully.  SRO=%u, Score is=%d:", rOff, score);
    printEditOpList(list);
    fprintf(stderr, "\n");
#endif

    // Find the extension
    EditOpList_t stackList;
    EditOpList_t * tempList = &stackList;
    initEditOpList(tempList, QS);
    int initAGS = findAGSExtension<FALSE>(QS, rOff, qStr, qOff, qLen, tempList, addedQLen, addedRLen);

#ifdef TEMPDEBUG
    fprintf(stderr, "qlen %d score %d ", qLen, initAGS);
    printEditOpList(tempList);
    fprintf(stderr, "\n");
#endif

    // If there is no positive extension, we are done.
    if (initAGS <= 0)
    {
        return 0;
    }

    // Now we need to go back through the extension making sure the score does not go below zero at any time.
    // Might as well calculate the lengths while we are at it.
    // The direction we need to loop through the editoplist is different depending on whether we are extending backward or forward.
    // The key is to start at the point farthest from the one that will attach to the rest of the query.

    // Going forward, we start at the end that is closest to the end of the query being extended.
    // Therefore, start with the existing score, and add as we go along.

    AlignmentArgs_t * AAs = QS->AAs;
    int QLen = 0;
    int RLen = 0;
    int AGS  = score;
    int maxAGS = score;
    EditOp_t * maxItem = NULL;
    int maxQLen = 0;
    int maxRLen = 0;
    forAllEditOpsInList(item, tempList)
    {
        EditOpCode op = item->opcode;
        int len = item->length;
        if      (op == MatchEditOp   )  { QLen += len; RLen += len; AGS += (AAs->MScore * len);               }
        else if (op == ReplaceEditOp )  { QLen += len; RLen += len; AGS -= (AAs->RCost * len);                }
        else if (op == InsertEditOp  )  { QLen += len;              AGS -= (AAs->GOCost + AAs->GECost * len); }
        else if (op == DeleteEditOp  )  {              RLen += len; AGS -= (AAs->GOCost + AAs->GECost * len); }
        if (AGS > maxAGS)
        {
            maxAGS = AGS;
            maxQLen = QLen;
            maxRLen = RLen;
            maxItem = item;
        }
        // If we hit zero, so we are done.
        // Decide whether any extension adds positive score.
        else if (AGS <= 0)
        {

#ifdef TEMPDEBUG
            fprintf(stderr, ":Before trim, added %u to QLen, and %u to RLen, and %d to score of %d.\n", *addedQLen, *addedRLen, initAGS, score);
            fprintf(stderr, "Trimming editoplist.  AGS=%d, maxAGS=%d, EOL: ", AGS, maxAGS);
            printEditOpList(tempList);
            fprintf(stderr, "\n");
#endif

            // The score went to zero, and the best was no better than what we stated with.
            // So, just keep what we had.
            if (maxAGS <= score)
            {
                emptyEditOpList(tempList);
                *addedQLen = 0;
                *addedRLen = 0;
                return 0;
            }
            // The score went to zero, but the best was better than what we stated with.
            // So, we want to truncate the extension at the best point.
            EditOpList_t unwantedTailStackList;
            EditOpList_t * unwantedTailList = &unwantedTailStackList;
            initEditOpList(unwantedTailList, QS);
            splitEditOpListAfter(tempList, unwantedTailList, maxItem);
            emptyEditOpList(unwantedTailList);
            *addedQLen = maxQLen;
            *addedRLen = maxRLen;
            initAGS = maxAGS - score;
#ifdef TEMPDEBUG
            fprintf(stderr, ":After trim, added %u to QLen, and %u to RLen, and %d to score.", *addedQLen, *addedRLen, initAGS);
            printEditOpList(tempList);
            fprintf(stderr, "\n");
#endif
            break;
        }
    }
#ifdef TEMPDEBUG
    fprintf(stderr, "sucessfully carefully extended. Adding: ");
    printEditOpList(tempList);
    fprintf(stderr, " to: ");
    printEditOpList(list);
#endif

    // Add the extension onto the alignment.
    // This also occurs with a "good" truncation of the original DP.
    mergeEOLToBack(list, tempList);

#ifdef TEMPDEBUG
    fprintf(stderr, " giving: ");
    printEditOpList(list);
    fprintf(stderr, "\n");
#endif

    // Return score.
    return initAGS;
}

int findAGSBackwardExtensionCarefully(QueryState_t * QS,  ROFF rOff, char * qStr, SUINT qOff, SUINT qLen,
                                      EditOpList_t * list, int score, SUINT * addedQLen, SUINT * addedRLen)
{
#ifdef TEMPDEBUG
    fprintf(stderr, "Entering extend backward carefully.  SRO=%u, qLen %d ", rOff, qLen);
    printEditOpList(list);
    fprintf(stderr, "\n");
#endif

    // Find the extension
    EditOpList_t stackList;
    EditOpList_t * tempList = &stackList;
    initEditOpList(tempList, QS);
    AlignmentArgs_t * AAs = QS->AAs;
    int initAGS = findAGSExtension<TRUE>(QS, rOff, qStr, qOff, qLen, tempList, addedQLen, addedRLen);

#ifdef TEMPDEBUG
    fprintf(stderr, "qlen %d score %d ", qLen, initAGS);
    printEditOpList(tempList);
    fprintf(stderr, "\n");
#endif

    // If there is no positive extension, we are done.
    if (initAGS <= 0)
    {
        return 0;
    }

    // Now we need to go back through the extension making sure the score does not go below zero at any time.
    // Might as well calculate the lengths while we are at it.
    // The direction we need to loop through the editoplist is different depending on whether we are extending backward or forward.
    // The key is to start at the point farthest from the one that will attach to the rest of the query.

    // Going backward, we start with the point of the extension farthest from the front of the query we are extending.
    // Therefore, we start with no score, and go all the way to the existing alignment, and see what the new and max score are.

    int QLen = 0;
    int RLen = 0;
    int AGS  = 0;
    int maxAGS = 0;
    EditOp_t * startItem = NULL;
    forAllEditOpsInList(item, tempList)
    {
        EditOpCode op = item->opcode;
        int len = item->length;
        if      (op == MatchEditOp   )  { QLen += len; RLen += len; AGS += (AAs->MScore * len);               }
        else if (op == ReplaceEditOp )  { QLen += len; RLen += len; AGS -= (AAs->RCost * len);                }
        else if (op == InsertEditOp  )  { QLen += len;              AGS -= (AAs->GOCost + AAs->GECost * len); }
        else if (op == DeleteEditOp  )  {              RLen += len; AGS -= (AAs->GOCost + AAs->GECost * len); }

        // If we hit zero, what is behind us is rubbish, so reset and look for a better place to start the extension.
        if (AGS <= 0)
        {
            AGS = 0;
            maxAGS = 0;
            QLen = 0;
            RLen = 0;
            startItem = item;
        }
        if (AGS > maxAGS) maxAGS = AGS;
    }

    // If the (remaining) extension has no value, stick with what we started with.
    // maxAGS - AGS is the amount that the extension drops in value between the max point and the point of attachment to the current alignment.
    // Therefore, if maxAGS - AGS > score, then the score would fall to zero if backtracking to the left from the current alignment.
    // If so, there might be a great alignment out to the left, but it is not part of a local extension of the current alignemnt.
    // So, let the caller find the bump to the left and leave the current alignment unextended.
    // Note that maxAGS >= AGS + score is just a better way to calculate maxAGS - AGS > score when AGS negative.
    if (AGS <= 0 || maxAGS >= AGS + score)
    {
        emptyEditOpList(tempList);
        *addedQLen = 0;
        *addedRLen = 0;
        return 0;
    }

#ifdef TEMPDEBUG
    fprintf(stderr, ":Before trim, added %u to QLen, and %u to RLen, and %d to score.\n", *addedQLen, *addedRLen, initAGS);
    fprintf(stderr, "Trimming editoplist.  AGS=%d, ", AGS);
    printEditOpList(tempList);
    fprintf(stderr, "\n");
#endif

    // If the startItem is not NULL, then we found a zero along the way.
    // So, trim the extension to the good part.
    if (startItem != NULL)
    {
        EditOpList_t wantedHeadStackList;
        EditOpList_t * wantedHeadList = &wantedHeadStackList;
        initEditOpList(wantedHeadList, QS);
        splitEditOpListAfter (tempList, wantedHeadList, startItem);

#ifdef TEMPDEBUG
        fprintf(stderr, ":After trim, added %u to QLen, and %u to RLen, and %d to score.", *addedQLen, *addedRLen, initAGS);
        printEditOpList(tempList);
        fprintf(stderr, "\n");
#endif

        emptyEditOpList(tempList);
        mergeEOLToFront(list, wantedHeadList);
    }
    // We are using the whole extension, so just add it on.
    else
    {
        mergeEOLToFront(list, tempList);
    }

#ifdef TEMPDEBUG
    fprintf(stderr, "sucessfully carefully extended backward. Resulting in: ");
    printEditOpList(list);
    fprintf(stderr, "\n");
#endif

    *addedQLen = QLen;
    *addedRLen = RLen;
    initAGS = AGS;
    return initAGS;
}
#undef TEMPDEBUG

////////////////////////
//     This is the main workhorse routine for SW calculations.
//     It is highly templated to handle the cases of global DP for portions of a query between seed hits in both full and banded versions.
//     It addition, it handles both bakcward and forward extensions in which one side of the SW is anchored, but the other side can float.
//     Such extensions are always banded, and heuristically, use XCutoff.
////////////////////////

template <BOOL banded, BOOL extension, BOOL reverse, BOOL XCutoff>
int findAffineGapScore(QueryState_t * QS, char * qStr, SUINT qLen, char * rStr, SUINT rLen, EditOpList_t * list, SUINT * addedQLen = NULL, SUINT * addedRLen = NULL)
{
    // As we will go row by row, then column by column, we do not need to store all of the V, E, and F arrays.
    // Instead we will store only the previous row and column for each.
    // However, we will have to store the entirety of the "walk back" arrays so that we can reclaim an editoplist.
    // As per standard literature usage, the following names are used for various quantities needed in score calculations.
    //  V    = the result score for a given cell.
    // PVRow = the previous row of V values.
    //  F    = the score for an insert at this cell.
    // PFRow = the previous row of F values.
    //  I    = insert length for current cell.
    // PIRow = previous row of insert length.
    //  E    = the score for a delete at this cell.
    // PECol = the previous column E value.
    //  D    = delete length for current cell.
    // PDCol = previous col of delete length.
    //  G    = the score for a match or replace at this cell.  This is really only for convenience, especially since we don't use a replacement matrix.
    // We also need the following to allow for proper reconstruction of the EditOpList from the score.
    // EOArray = an n*m array of "back pointers".  We use EditOpCodes instead, as they have the same information content.
    // IDCount = an n*m array of the count/length of inserts and deletes needed to handle insert/delete extension "tunneling".

    addToTimer(junkTimer);
    // Get local copies.
    AlignmentArgs_t * AAs = QS->AAs;
    int GOCost = AAs->GOCost;
    int GECost = AAs->GECost;
    int RCost  = AAs->RCost;
    int MScore = AAs->MScore;

    // Try putting array ptrs on the stack to avoid going indirect through the structure each time.
    DPScore * PVRow = QS->PVRowSave;
    DPScore * PFRow = QS->PFRowSave;
    DPScore * PIRow = QS->PIRowSave;
    DPScore * IDCount = QS->IDCount;
    EditOpCode * EOArray = QS->EOArray;
    DPScore   PVCol = DPWorstScore;

    // Let's calculate some parameters for the banded versions.
    // The left and right "bandwidth" can be different to account for the amount
    //     that the array differs from square, i.e. the net indel size.
    // This is because the bandwidth will need to be wider to
    //     to allow for all of the insert or delete to come in one chunk,
    int bandwidth, leftBandwidth, rightBandwidth;
    int arrWidth, maxi, maxj;
    int arrHeight = qLen + 1;
    if (banded)
    {
        // For extensions, we don't know which way to go.
        if (extension)
        {
            bandwidth = 2 * AAs->bandWidth;
            leftBandwidth = rightBandwidth = bandwidth;
            maxi = maxj = 0;
        }
        else
        {
            // We do know which direction an indel is, so be wider on that side.
            bandwidth = AAs->bandWidth;
            if (rLen > qLen)
            {
                rightBandwidth = bandwidth + (rLen - qLen);
                leftBandwidth = bandwidth;
            }
            else
            {
                leftBandwidth = bandwidth + (qLen - rLen);
                rightBandwidth = bandwidth;
            }
            // Tell the backtracker where to start.
            maxi = qLen;
            maxj = rightBandwidth;
        }
        arrWidth = leftBandwidth + rightBandwidth + 1;
    }
    else // Full DP
    {
        arrWidth = rLen + 1;
        // Tell the backtracker where to start.
        maxi = arrHeight - 1;
        maxj = arrWidth - 1;
    }

    // Array initializations.
    // NOTE:  We currently always use the first portion of the array (ie row widths are set per alignment).
    //        This demands we do the inits for each new size of alignment.
    int startInit, endInit;
    // Initialize the arrays with last ditch delete operations in the first row.
    if (banded)
    {
        startInit = leftBandwidth + 1;
        // As these change location for each banded call, we will need to initialize them here, not in DPInit.
        // As we don't need to set these for banded, we will set them here.
        EOArray[leftBandwidth] = UnknownEditOp;
        IDCount[leftBandwidth] = 0;
        // These are beyond the end of the calculated region, so won't be changed and can be safely set here.
        PFRow[arrWidth] = PVRow[arrWidth] = DPWorstScore;
        PIRow[arrWidth] = 0;
    }
    else startInit = 1;
    endInit = arrWidth;
    int deleteCount = 1;
    for (int j=startInit; j<endInit; j++)
    {
        // Cost of a delete from the first row.
        EOArray[j] = DeleteEditOp;
        IDCount[j] = deleteCount;
        PVRow[j] = -(GOCost + deleteCount * GECost);
        deleteCount += 1;
        // We can't do an insert from here.
        PFRow[j] = DPWorstScore;
        PIRow[j] = 0;
    }

    // Because banded will use column zero, these values cannot be set in the init routine.
    // These are the values that correspond to the place the backtrace will end.
    PFRow[startInit-1] = 0;
    PIRow[startInit-1] = 0;
    PVRow[startInit-1] = 0;
    // These are done in the init routine for non-banded, as unbanded always ends at location 0
    //     which is never used by banded.
    // EOArray[startInit-1] = UnknownEditOp;
    // IDCount[startInit-1] = 0;

    // Initialize the arrays with last ditch insert operations.
    // For full DP, this is the first column.
    // For banded, it is an anti-diagonal.
    if (banded) endInit = leftBandwidth;
    else endInit = arrHeight - 1;
    for (int i=1; i<=endInit; i++)
    {
        int loffset = (banded) ? leftBandwidth - i : 0;
        // Cost of an insert from the first column.
        (EOArray + i * arrWidth)[loffset] = InsertEditOp;
        (IDCount + i * arrWidth)[loffset] = i;
    }

    // If we are doing an extension, keep track of the maximum V value.
    // We will also keep track of the location of the maximum score cell from which to start the backtrack.
    DPScore maxScore;
    int Cutoff;
    if (extension)
    {
        maxScore = DPWorstScore;
        Cutoff = AAs->XCutoff;
    }

    // Set the start and end column locations for unbanded case.
    // For banded, these are set inside the loop.
    int startCol, endCol;
    if (! banded)
    {
        startCol = 1;
        endCol = arrWidth - 1;
    }
    ///////////////
    // Start of main loop.
    ///////////////
    // For each row................................
    int rowOffset = 0;
    DPScore V = 0;
    for (int i=1; i<arrHeight; i++)
    {
        // This is needed to treat the 1D arrays as if they were 2D.
        rowOffset += arrWidth;

        // We can't do a continued delete from the first column.
        DPScore PDCol = 0;
        DPScore PECol = DPWorstScore;

        // Now calculate the starting and ending column for the inner loop.
        // And set edge values of PVCol and PVRow accordingly.
        if (banded)
        {
            // Calculate the start and end columns for this iteration.
            // First the start column.
            // Make sure PVRow and PVCol has the proper insert cost.
            startCol = leftBandwidth + 1 - i;
            if (startCol <= 0)
            {
                startCol = 0;
                PVCol = DPWorstScore;
            }
            else PVCol = PVRow[startCol -1] = -(GOCost + i * GECost);
            // Now the end column.
            endCol = MIN(leftBandwidth + rLen - i, arrWidth - 1);
        }
        else
        {
            // Make sure PVCol has the correct insert cost for non-banded.
            PVCol =  -(GOCost + i * GECost);
        }

        // Keep track of the best score in the row for XCutoff
        int rowMaxScore;
        if (XCutoff)
        {
            rowMaxScore = DPWorstScore;
        }
        // The query offset does not change throughout the row, so calc once here.
        // We reversed the reference string when we decompressed it, but not the query.
        char qChar = (reverse) ? qStr[1 - i] : qStr[i - 1];
        ROFF rRowStartOff;
        if (banded) rRowStartOff = i - leftBandwidth - 1;

        ///////////////
        // Main inner loop.
        ///////////////
        // For each column...........
        for (int j=startCol; j<=endCol; j++)
        {
            // We need to set the relative locations of RM and I offsets.
            // This is what supports the columnar array shape (instead of diagonal) for banded case.
            int RMOff = (banded) ? j   : j-1;
            int IOff  = RMOff + 1;

            // Calculate the current E, F, G values, while maximizing V as we go.
            // Start with G, the score for the non-indel case.
            EditOpCode opcode;
            V = PVRow[RMOff];
            char rChar = (banded) ? rStr[rRowStartOff + j] : rStr[j - 1];
            if (codesEqual(qChar, rChar)) { V += MScore; opcode = MatchEditOp;   }
            else                          { V -= RCost;  opcode = ReplaceEditOp; }
            // if (extension) fprintf(stderr, "i %d j %d V %d rOff %d qOff %d rChar %d qChar %d \n", i, j, V, rRowStartOff + j, (reverse) ? 1-i : i-1, rChar, qChar);

            // This next line doesn't seem to really be needed, as we will never use the value for M or R.
            // (IDCount + rowOffset)[j] = 1;

            // Calculate E, the cost of a deletion.
            // D is the number of bases in the best deletion.
            // CE is cost of a continuing deletion.
            DPScore CE = PECol - GECost;
            // NE is cost of a newly started deletion.
            DPScore NE = PVCol - (GOCost + GECost);
            if (CE >= NE && (PDCol + 1) <= AAs->maxIntron) { PECol = CE; PDCol = PDCol + 1; }
            else                                           { PECol = NE; PDCol = 1;         }
            // In extensions, we want to greedily accept indels early in the array vs. later in the array.
            // This allows for more positive score after them to define a longer overall positive extension.
            if ((extension) ? (PECol >= V) : (PECol > V))
            {
                V = PECol;
                opcode = DeleteEditOp;
                (IDCount + rowOffset)[j] = PDCol;
            }

            // Calculate F, the cost of an insertion.
            // I is the number of bases in the best insertion.
            // CF is cost of a continuing insertion.
            DPScore F, I;
            DPScore CF = PFRow[IOff] - GECost;
            // NF is cost of a newly started insertion.
            DPScore NF = PVRow[IOff] - (GOCost + GECost);
            if (CF >= NF && (PIRow[IOff] + 1) <= AAs->maxGap) { F = CF; I = PIRow[IOff] + 1; }
            else                                              { F = NF; I = 1;               }
            // In extensions, we want to greedily accept indels early in the array vs. later in the array.
            // This allows for more positive score after them to define a longer overall positive extension.
            if (((extension) ? (F >= V) : (F > V)))
            {
                // fprintf(stderr, "Setting insertEditOp.  F=%d, V=%d, I=%d, PIRow[IOff]=%d\n", F, V, I, PIRow[IOff]);
                V = F;
                opcode = InsertEditOp;
                (IDCount + rowOffset)[j] = I;
            }
            // Keep track of these values for the next row.
            PFRow[j] = F;
            PIRow[j] = I;

            // We now have the new V value and the opcode.
            // Set the score in the backtrack array.
            (EOArray + rowOffset)[j] = opcode;
            // See if we have a new row max.
            if (XCutoff && V > rowMaxScore) rowMaxScore = V;
            // If we are extending, and the score is greater than the maxscore, keep track of this for start location of backtrack.
            if (extension)
            {
                if (V > maxScore)
                {
                    maxScore = V;
                    maxi = i;
                    maxj = j;
                }
            }
            // Get ready for next column.
            if (banded) PVRow[j] = V;
            else        PVRow[j-1] = PVCol;
            PVCol = V;
        }
#ifdef ALLDEBUG
        for (int k=0; k<arrWidth; k++) fprintf(stderr, "%4d ", VRow[k]);
        fprintf(stderr, "\n");
#endif
        // Get ready for next row.
        // First check for the XCutoff criteria to break out early.
        if (XCutoff && extension && rowMaxScore < (maxScore - Cutoff)) break;
        // For non banded, we need to put the last V value up for the next iteration.
        if (!banded) PVRow[endCol] = V;
    }
    // Done with the main DP loop.
    // Keep track of times and prepare for the backtrack loop.
    if (extension) {addToTimer(ExtensionSWTimer);}
    else           {addToTimer(FullSWTimer);}
    // The retval is the overall score for the DP run.
    int retval = (extension) ? maxScore : V;
    // No good extension found, so punt.
    if (extension && retval <= 0) return 0;

    // If we are extending, we know where the extension ends.
    // So, set the added query and reference extension lengths.
    if (extension)
    {
        initEditOpList(list, QS);
        *addedQLen = maxi;
        *addedRLen = maxi + (maxj - bandwidth);
    }

    ///////////////
    // Start of the backtrack loop to recover the edit operations for the best path.
    // This used to be a separate function, but was moved here to support threading,
    //     as too much state was needed between the DP loop and the backtrack loop.
    ///////////////

#ifdef ALLDEBUG
    fprintf(stderr, "Entering DP backtrack(%s %s %s). arrWidth=%d, arrHeight=%d, x=%d, y=%d.\n",
            (banded) ? "Banded" : "Full",
            (reverse) ? "Reverse" : "Forward",
            (extension) ? "Extension" : "Anchored",
            arrWidth, arrHeight, x, y);
    for (int y = 0; y < MIN(arrHeight, maxi); y++)
    {
        EditOpCode * EOPtr = EOArray + (y * arrWidth);
        DPScore * IDPtr = IDCount + (y * arrWidth);
        for(int x = 0; x < arrWidth; x++)
        {
            fprintf(stderr, "%d%c ", IDPtr[x], EditOpCode2char(EOPtr[x]));
        }
        fprintf(stderr, "\n");
    }
#endif

    // Set where to start the backtrack.
    int x = maxj;
    int y = maxi;

    // Start with a half unwind of the loop.
    // Get the next edit op and length for indel.
    EditOpCode * EOPtr = EOArray + (y * arrWidth);
    DPScore * IDPtr = IDCount + (y * arrWidth);

    // Define variables used in the loop.
    // Note that we need to know the type of the "previous" opcode to know when the opcode type has changed.
    EditOpCode EOCode;
    int EOLen;
    EditOpCode prevEOCode = EOPtr[x];
    int opLen = 0;
    // Now we are ready to start the main loop.
    // We iterate over both the edit op array and the indel length array in parallel.
    while ((EOCode = EOPtr[x]) != UnknownEditOp)
    {
        EOLen = IDPtr[x];
#ifdef ALLDEBUG
        fprintf(stderr, "Recovering EditOpList from DP arrays.  x=%d, y=%d, EOLen=%d, EOCode=%c, prevEOCode=%c\n",
                x, y, EOLen, EditOpCode2char(EOCode), EditOpCode2char(prevEOCode));
#endif
        // Both banded and not handle indels the same.
        // Deletes go back in the row.
        // Inserts for forward in the row, but backward up the columns.
        if (banded)
        {
            if      (EOCode == DeleteEditOp) { x -= EOLen;                                                                   }
            else if (EOCode == InsertEditOp) { x += EOLen; EOPtr -= (EOLen*arrWidth); IDPtr -= (EOLen*arrWidth);             }
            // Match or mismatch.
            // Go up one column.
            else                             {             EOPtr -= arrWidth;         IDPtr -= arrWidth;          EOLen = 1; }
        }
        else
        {
            if      (EOCode == DeleteEditOp) { x -= EOLen;                                                                   }
            else if (EOCode == InsertEditOp) {             EOPtr -= (EOLen*arrWidth); IDPtr -= (EOLen*arrWidth);             }
            // Match or mismatch.
            // Go back up the current diagonal.
            else                             { x -= 1;     EOPtr -= arrWidth;         IDPtr -= arrWidth;          EOLen = 1; }
        }

        // Now store the edit operation if at an inflection.
        if (prevEOCode != EOCode)
        {
            if (reverse) addEditOpToBack(list, prevEOCode, opLen);
            else         addEditOpToFront(list, prevEOCode, opLen);
            prevEOCode = EOCode;
            opLen = EOLen;
        }
        else opLen += EOLen;
    }

    // Partial unwind of end of loop.
    // Push on final operation.
    if (reverse) addEditOpToBack(list, prevEOCode, opLen);
    else         addEditOpToFront(list, prevEOCode, opLen);

#ifdef ALLDEBUG
    int score = EditOpList2AGS(list, AAs);
    fprintf(stderr, "%6d::", score);
    printEditOpList(list);
    fprintf(stderr, "\n");
#endif

    if (extension) {addToTimer(ExtensionBacktrackSWArraysTimer);}
    else           {addToTimer(FullBacktrackSWArraysTimer);}

    return retval;
}

///////////////////
// editoplist helper routines.
///////////////////

// Find out if there is a match edit op at least min in length.
BOOL EditOpList2Maxmatch(EditOpList_t * list, int min)
{
    forAllEditOpsInList(item, list)
    {
        if (item->opcode == MatchEditOp && item->length >= min) return TRUE;
    }
    return FALSE;
}

// Calculate the AGS for an edit op list.
int EditOpList2AGS(EditOpList_t * list, AlignmentArgs_t * AAs)
{
    int score = 0;
    forAllEditOpsInList(item, list)
    {
        EditOpCode op  = item->opcode;
        int        len = item->length;
        if      (op == MatchEditOp)                        score += (AAs->MScore * len);
        else if (op == ReplaceEditOp)                      score -= (AAs->RCost  * len);
        else if (op == InsertEditOp || op == DeleteEditOp) score -= (AAs->GOCost + AAs->GECost * len);
    }
    return score;
}
