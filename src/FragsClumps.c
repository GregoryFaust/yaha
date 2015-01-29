/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    FragsClumps.c   Includes code for fragment and clump manipulation.
                             As this file has so many inlines, we will #include it into AlignmentHelpers.c!!!
                             As such, it is not supposed to be separately compilable.
                             Separating the file into a .c and .inl might also work, but as one file,
                                  the compiler seems to put additional small functions inline that are not declared inline.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include "FragsClumps.inl"

/////
//     Fragments
/////

// Reset an existing fragment.
void resetFragment(Fragment_t *frag)
{
    frag->startQueryOff = SUINT_MAX_VALUE;
    frag->startRefOff = ROFF_MAX_VALUE;
    frag->endQueryOff = 0;
}

// Copy the data fields from one fragment to another.
// Sort of like a copy constructor, except we don't create a new fragment.
// So, caller is responsible for the memory management.
void copyFragment (Fragment_t * fragIn, Fragment_t * fragOut)
{
    fragOut->startQueryOff = fragIn->startQueryOff;
    fragOut->endQueryOff = fragIn->endQueryOff;
    fragOut->startRefOff = fragIn->startRefOff;
    fragOut->refLen = fragIn->refLen;
}

void printDebugFragment(Fragment_t * frag)
{
    fprintf(stderr, "R[%u-%u]=D:%u=Q[%d:%d-%d]\n", frag->startRefOff, fragEndRefOff(frag),
            frag->startRefOff - frag->startQueryOff, fragMatchCount(frag), frag->startQueryOff, frag->endQueryOff);
}

void printDebugFragmentExtended(Fragment_t * frag, AlignmentArgs_t * AAs)
{
    fprintf(stderr, "R[%u-%u]=D:%u=Q[%d:%d-%d]", frag->startRefOff, fragEndRefOff(frag),
            frag->startRefOff - frag->startQueryOff, fragMatchCount(frag), frag->startQueryOff, frag->endQueryOff);
    BaseSequence_t * BS = findBaseSequence(AAs->BSs, frag->startRefOff);
    fprintf(stderr, " starts at QO %d in ", frag->startQueryOff);
    printChars(stderr, BS->name, 0, BS->nameLen);
    fprintf(stderr, " at offset:%u\n", frag->startRefOff - BS->startingOffset);
}

/////
//     Slab Allocator common code
//     We will lazily allocate more slabs when needed, and only clean up at the end.
/////
//     Currently used for clump and SFrag linked lists.
/////

char * pushNewLBMallocBlock(int blockSize, LBMallocBlock_t **blockArrayPtr)
{
    char * newBlock = blockMalloc(blockSize, NULL);
    LBMallocBlock_t * newMallocBlock = (LBMallocBlock_t *)malloc(sizeof(LBMallocBlock_t));
    if (newMallocBlock == NULL) fatalError("Insufficeint memory available to allocate (more) objects.");
    newMallocBlock->size = blockSize;
    newMallocBlock->block = newBlock;
    newMallocBlock->next = *blockArrayPtr;
    *blockArrayPtr = newMallocBlock;
    return newBlock;
}

void freeLBMallocBlocks(LBMallocBlock_t * block)
{
    while (block != NULL)
    {
        LBMallocBlock_t * nextBlock = block->next;
        blockFree(block->block, block->size);
        free(block);
        block = nextBlock;
    }
}

/////
//     Now the clump related routines.
/////
//     A clump is a collection of non-overlapping SFragments.
//     All gaps in the clump should be between SFragments, or represented by SFragments,
//          depending where we are in the alignment process.
/////

void makeMoreClumps(QueryState_t * QS)
{
    int clumpCount = 4096;
    Clump_t * clumpArray = (Clump_t *)pushNewLBMallocBlock(sizeof(Clump_t) * clumpCount, &QS->clumpBlockList);
    for (int i=1; i<clumpCount; i++)
    {
        (clumpArray + (i - 1))->next = (clumpArray + i);
    }
    (clumpArray + (clumpCount - 1))->next = NULL;
    QS->freeClumpList = clumpArray;
}

Clump_t * makeClump(QueryState_t * QS)
{
    Clump_t * clump;
    if (QS->freeClumpList == NULL) makeMoreClumps(QS);
    clump = QS->freeClumpList;
    QS->freeClumpList = QS->freeClumpList->next;
    initEditOpList(EOLPtr(clump), QS);
    clump->SFragList.head = clump->SFragList.tail = NULL;
    resetClump(clump, QS);
    return clump;
}

void resetClump(Clump_t * clump, QueryState_t * QS)
{
    clump->totScore = 0;
    clump->totLength = 0;
    clump->matchedBases = 0;
    clump->status = 0;
    clump->mapQuality = 255;
    emptyEditOpList(EOLPtr(clump));
    emptySFragmentList(SFragListPtr(clump), QS);
    clump->next = NULL;
}

void disposeClump(Clump_t * clump, QueryState_t * QS)
{
    emptyEditOpList(EOLPtr(clump));
    emptySFragmentList(SFragListPtr(clump), QS);
    QS->freeClumpList = pushClump(clump, QS->freeClumpList);
}

/////
//     Super Fragments and SFragment Lists.
/////

void makeMoreSFrags(QueryState_t * QS)
{
    int SFragCount = 4096;
    SFragment_t * SFragArray = (SFragment_t *)pushNewLBMallocBlock(sizeof(SFragment_t) * SFragCount, &QS->SFragBlockList);
    for (int i=1; i<SFragCount; i++)
    {
        (SFragArray + (i - 1))->next = (SFragArray + i);
    }
    (SFragArray + (SFragCount - 1))->next = NULL;
    QS->freeSFragmentList = SFragArray;
}

SFragment_t * makeSFragment(QueryState_t * QS)
{
    SFragment_t * SFrag;
    if (QS->freeSFragmentList == NULL) makeMoreSFrags(QS);
    SFrag = QS->freeSFragmentList;
    QS->freeSFragmentList = QS->freeSFragmentList->next;
    resetFragment(getFragPtr(SFrag));
    initEditOpList(getEOLPtr(SFrag), QS);
    SFrag->score = 0;
    SFrag->prev = SFrag->next = NULL;
    return SFrag;
}

static inline void emptySFragment(SFragment_t * SFrag)
{
    emptyEditOpList(getEOLPtr(SFrag));
}

void disposeSFragment(SFragment_t * SFrag, QueryState_t * QS)
{
    emptyEditOpList(getEOLPtr(SFrag));
    SFrag->next = QS->freeSFragmentList;
    QS->freeSFragmentList = SFrag;
}

/////
//     SFragment Lists
/////

// The next ptr is to support iterators.
// We will make next point one ahead so that we can use this mechanism even for popping element off top.
SFragment_t * getFirstSFrag(SFragmentList_t * list)
{
    // Check if we have any nodes at all.
    if (list->head == NULL) return NULL;
    list->next = list->head->next;
    return list->head;
}

SFragment_t * getNextSFrag(SFragmentList_t * list)
{
    if (list->next == NULL) return NULL;
    SFragment_t * nextNode = list->next;
    list->next = nextNode->next;
    return nextNode;
}

SFragment_t * getLastSFrag(SFragmentList_t * list)
{
    // Check if we have any nodes at all.
    if (list->tail == NULL) return NULL;
    list->next = list->tail->prev;
    return list->tail;
}

SFragment_t * getPrevSFrag(SFragmentList_t * list)
{
    if (list->next == NULL) return NULL;
    SFragment_t * nextNode = list->next;
    list->next = nextNode->prev;
    return nextNode;
}

// As SFragmentLists are embedded in the clump structure, currently this is never called.
SFragmentList_t * makeSFragmentList()
{
    SFragmentList_t * list = (SFragmentList_t *)malloc(sizeof(SFragmentList_t));
    list->head = list->tail = NULL;
    return list;
}

void emptySFragmentList(SFragmentList_t * list, QueryState_t * QS)
{
    // We will first give each SFragment a chance to empty itself.
    // Then instead of freeing the SFragments, we will put them on the freeList.
    if (!SFragmentListIsEmpty(list))
    {
        for (SFragment_t * SFrag = getFirstSFrag(list); SFrag != NULL; SFrag = getNextSFrag(list))
        {
            emptySFragment(SFrag);
        }
        list->tail->next = QS->freeSFragmentList;
        QS->freeSFragmentList = list->head;
    }
    list->head = list->tail = NULL;
}

// As SFragmentLists are embedded in the clump structure, currently this is never called.
void disposeSFragmentList(SFragmentList_t * list, QueryState_t * QS)
{
    emptySFragmentList(list, QS);
    free(list);
}

// For functions that take nodes from the list as arguments,
//     we have to assume that the node actually is in the list.
// Otherwise, things would be too slow.
// So, we will check only in DEBUG mode.
BOOL SFragInList(SFragmentList_t * list, SFragment_t * SFrag)
{
    for (SFragment_t * currSFrag = getFirstSFrag(list); currSFrag != NULL; currSFrag = getNextSFrag(list))
    {
        if (currSFrag == SFrag) return TRUE;
    }
    return FALSE;
}

SFragment_t * findFollowingSFrag(SFragmentList_t * list, SUINT queryOff)
{
    for (SFragment_t * currSFrag = getFirstSFrag(list); currSFrag != NULL; currSFrag = getNextSFrag(list))
    {
        if (currSFrag->frag.startQueryOff > queryOff) return currSFrag;
    }
    return NULL;
}


// Helper functions to add and remove a node between two others.
// But it also looks for beginnning and end of list.
// Therefore, it can also be used for additions at the beginning and and of the list.
void addSFragmentToList(SFragmentList_t * list, SFragment_t * prev, SFragment_t * next, SFragment_t * newNode)
{
    newNode->prev = prev;
    newNode->next = next;
    if (prev == NULL) list->head = newNode;
    else prev->next = newNode;
    if (next == NULL) list->tail = newNode;
    else next->prev = newNode;
}

static inline void addSFragmentToFront(SFragmentList_t * list, SFragment_t * SFrag) __attribute__((always_inline));
static inline void addSFragmentToFront(SFragmentList_t * list, SFragment_t * SFrag)
{
    addSFragmentToList(list, NULL, list->head, SFrag);
}

static inline void addSFragmentToBack(SFragmentList_t * list, SFragment_t * SFrag) __attribute__((always_inline));
static inline void addSFragmentToBack(SFragmentList_t * list, SFragment_t * SFrag)
{
    addSFragmentToList(list, list->tail, NULL, SFrag);
}

static inline void addSFragmentBefore(SFragmentList_t * list, SFragment_t * next, SFragment_t * SFrag) __attribute__((always_inline));
static inline void addSFragmentBefore(SFragmentList_t * list, SFragment_t * next, SFragment_t * SFrag)
{
    ASSERT(SFragInList(list, next), "Attempt to add SFragment before node not in list.");
    addSFragmentToList(list, next->prev, next, SFrag);
}

static inline void addSFragmentAfter(SFragmentList_t * list, SFragment_t * prev, SFragment_t * SFrag) __attribute__((always_inline));
static inline void addSFragmentAfter(SFragmentList_t * list, SFragment_t * prev, SFragment_t * SFrag)
{
    ASSERT(SFragInList(list, prev), "Attempt to add SFragment after node not in list.");
    addSFragmentToList(list, prev, prev->next, SFrag);
}

void removeSFragment(SFragmentList_t * list, SFragment_t * curNode)
{
    ASSERT(SFragInList(list, curNode), "Attempt to remove SFragment not in list.");
    SFragment_t * next = curNode->next;
    SFragment_t * prev = curNode->prev;
    if (prev == NULL) list->head = next;
    else prev->next = next;
    if (next == NULL) list->tail = prev;
    else next->prev = prev;
    // We don't always want to do this.
    // So, when we do, we do it in the caller.
    // disposeSFragment(curNode);
}

void makeMoreEditOps(QueryState_t * QS);
void initClumpStructures(QueryState_t * QS)
{
    QS->clumpBlockList = NULL;
    QS->freeClumpList = NULL;
    makeMoreClumps(QS);
    QS->freeSFragmentList = NULL;
    QS->SFragBlockList = NULL;
    QS->EditOpArray = NULL;
    makeMoreSFrags(QS);
    QS->EditOpFreeList = EditOpPtrNULL;
    QS->EditOpArrayMaxElements = 0x100000;
    makeMoreEditOps(QS);
}

void cleanUpClumpStructures(QueryState_t * QS)
{
    // Clean up the Clump and SFragment free lists.
    // Since we are using a slab allocator, we need only clean up the blocks!
    freeLBMallocBlocks(QS->clumpBlockList);
    freeLBMallocBlocks(QS->SFragBlockList);
    free(QS->EditOpArray);
}

void printDebugClump(Clump_t * clump, AlignmentArgs_t * AAs, FILE * out)
{
    // First get info about the sequence the clump is in.
    SFragment_t * SFrag = SFragListPtr(clump)->head;
    if (SFrag == NULL)
    {
        fprintf(out, "Debug output of empty clump.\n");
        return;
    }
    Fragment_t * frag = getFragPtr(SFrag);
    BaseSequence_t * BS = findBaseSequence(AAs->BSs, frag->startRefOff);
    fprintf(out, "Clump [status:%02X] starts at QO[%c]:%d in ", clump->status,
            isReversed(clump) ? '-' : '+', frag->startQueryOff);
    printChars(out, BS->name, 0, BS->nameLen);
    fprintf(out, " at offset:%u\n", frag->startRefOff - BS->startingOffset);
    for (SFrag = getFirstSFrag(SFragListPtr(clump)); SFrag != NULL; SFrag = getNextSFrag(SFragListPtr(clump)))
    {
        frag = getFragPtr(SFrag);
        fprintf(out, "R[%u-%u]=D:%u=Q[%d:%d-%d] ", frag->startRefOff-BS->startingOffset,
                fragEndRefOff(frag)-BS->startingOffset, frag->startRefOff - frag->startQueryOff,
                fragMatchCount(frag), frag->startQueryOff, frag->endQueryOff);
        EditOpList_t * list = getEOLPtr(SFrag);
        if (!EOLisEmpty(list))
        {
            fprintf(stderr, "\n");
            printEditOpList(list);
            int AGS = EditOpList2AGS(list, AAs);
            fprintf(out, "[frag score=%d calc score = %d]", SFrag->score, AGS);
            validateEditOpList(list);
            fprintf(stderr, "\n");
        }
    }
    fprintf(out, "\n");
    EditOpList_t * list = EOLPtr(clump);
    if (list != NULL && list->head != EditOpPtrNULL)
    {
        printEditOpList(list);
        int AGS = EditOpList2AGS(list, AAs);
        fprintf(out, "[calc score = %d]", AGS);
        fprintf(stderr, "\n");
    }
}

#ifdef DEBUG
BOOL validateClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS, FILE * out, const char * outString)
{
    BOOL retval = TRUE;
    // First check to make sure the fragments have rational lengths and are in order
    for (SFragment_t * SFrag = getFirstSFrag(SFragListPtr(clump)); SFrag != NULL; SFrag = getNextSFrag(SFragListPtr(clump)))
    {
        Fragment_t * frag = getFragPtr(SFrag);
        if (frag->endQueryOff < frag->startQueryOff)
        {
            fprintf(out, "In validate Clump, frag has endQueryOff(%u) less than startQueryOff(%u).\n",
                    frag->endQueryOff, frag->startQueryOff);
            retval = FALSE;
        }
        if (fragEndRefOff(frag) < frag->startRefOff)
        {
            fprintf(out, "In validate Clump, frag has endRefOff(%u) less than startRefOff(%u) in ",
                    fragEndRefOff(frag), frag->startRefOff);
            BaseSequence_t * BS = findBaseSequence(AAs->BSs, frag->startRefOff);
            printChars(out, BS->name, 0, BS->nameLen);
            fprintf(out, " at offset=%u.\n", frag->startRefOff - BS->startingOffset);
            retval = FALSE;
        }
    }
    if (retval == FALSE)
    {
        fprintf(out, "Above validate Clump info about QueryID.");
        fprintQueryId(QS, out);
        fprintf(out, "Context=%s\n", outString);
    }
    return retval;
}
#else
static inline BOOL validateClump(Clump_t * clump, AlignmentArgs_t * AAs, FILE * out, const char * outString){return TRUE;}
#endif
