/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    QueryHeap.inl   Heap functions to sort the hash table results into contiguous fragments.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

/////
//     Binary Heap
/////

// We include these first two functions here even though they are not inline.
// This only works because this file is only included in one place: QueryMatch.c
Heap_t * makeHeap(int maxSizeArg)
{
    Heap_t * heap = (Heap_t *) malloc(sizeof(Heap_t));
    heap->contents = (UINT64 *) malloc(maxSizeArg*sizeof(UINT64));
    heap->size = 0;
    heap->maxSize = maxSizeArg;
    return heap;
}

void disposeHeap(Heap_t * heap)
{
    free(heap->contents);
    free(heap);
}

// For some reason these work faster as macros than as inlines.
// Common sense would indicate they should be identical!!
#ifdef NOTNOW
static inline UINT64 getItem(Heap_t * h, int index)
{
    return h->contents[index];
}
static inline void setItem(Heap_t * h, int index, UINT64 key)
{
    h->contents[index] = key;
}
#else
#define getItem(h, index) h->contents[index]
#define setItem(h, index, item) h->contents[index] = item
#endif

static inline int parent (int index)
{
    return index / 2;
}

static inline int left (int index)
{
    return (2 * index);
}

static inline int right (int index)
{
    return (2 * index) + 1;
}

static inline UINT64 encodeHeapItem(ROFF key, QOFF val)
{
    return ((((UINT64)key) << 32) + val);
}

static inline void decodeHeapItem(UINT64 temp, ROFF * key, QOFF * val)
{
    *val = temp & 0XFFFFFFFF;
    *key = temp >> 32;
}

static inline void decodeHeapTop(Heap_t * heap, ROFF * key, QOFF * val)
{
    decodeHeapItem(getItem(heap, 1), key, val);
}

static inline UINT64 removeHeapLeaf(Heap_t * heap)
{
    return getItem(heap, heap->size--);
}

static inline BOOL isHeapEmpty(Heap_t *  heap)
{
    return (heap->size == 0);
}

static inline void subHeapify(Heap_t * heap, int index)
{
    int curIndex = index;
    UINT64 bigItem = getItem(heap, index);
    int newIndex, minIndex;
    int hSize = heap->size;
    while (TRUE)
    {
        // Check if node to the left has smaller value than larger of right or current root.
        newIndex = left(curIndex);
        if (newIndex > hSize) break;  // Since left is smaller than right, we can exit if too large.
        if (getItem(heap, newIndex) < bigItem) minIndex = newIndex;
        else                                   minIndex = curIndex;
        // Check if node to the right has smaller value than current root.
        newIndex = right(curIndex);
        if (newIndex <= hSize && getItem(heap, newIndex) < getItem(heap, minIndex)) minIndex = newIndex;
        // Check if one of the left or right (or both) was smaller.
        // If so, swap with the smallest of them, and heapify what remains.
        if (minIndex == curIndex) break;
	// Swap items
        setItem(heap, curIndex, getItem(heap, minIndex));
        setItem(heap, minIndex, bigItem);
        curIndex = minIndex;
    }
 }

static inline void fastHeapify(Heap_t * heap, UINT64 bigItem)
{
    setItem(heap, 1, bigItem);
    subHeapify(heap, 1);
}

static inline void fixHeap(Heap_t * heap)
{
    for (int i=heap->size/2; i>0; i--)
    {
        subHeapify(heap, i);
    }
}

#define INVALID_ROFF 0
static inline void heapRemoveMinAndAdd(Heap_t * heap, ROFF key, QOFF val)
{
    UINT64 bigItem;
    if (key == INVALID_ROFF)
    {
        // The remove case
        bigItem = getItem(heap, heap->size);
        heap->size -= 1;
    }
    else
    {
        // The remove and add case
        bigItem = encodeHeapItem(key, val);
    }
    // Now we need to heapify to move the big item down in the tree.
    fastHeapify(heap, bigItem);
}

static inline BOOL heapPoll(Heap_t * heap, ROFF * key, QOFF * val)
{
    if (isHeapEmpty(heap)) return FALSE;

    // Get the return values.
    decodeHeapTop(heap, key, val);

    // Now we need to heapify to move the big item down in the tree.
    fastHeapify(heap, removeHeapLeaf(heap));
    return TRUE;
}

static inline void heapAddAtEnd(Heap_t * heap, ROFF key, QOFF val)
{
    UINT64 newItem = encodeHeapItem(key, val);
    int curIndex = (heap->size += 1);
    setItem(heap, curIndex, newItem);
}

static inline void heapAdd(Heap_t * heap, ROFF key, QOFF val)
{
    UINT64 newItem = encodeHeapItem(key, val);
    int curIndex = (heap->size += 1);
    setItem(heap, curIndex, newItem);
    int parentIndex;
    UINT64 parentItem;
    while (((parentIndex = parent(curIndex)) > 0) && ((parentItem = getItem(heap, parentIndex)) > newItem))
    {
        // Swap items.
        setItem(heap, parentIndex, newItem);
        setItem(heap, curIndex, parentItem);
        curIndex = parentIndex;
    }
}

static inline BOOL heapPeak(Heap_t * heap, ROFF * key, QOFF * val)
{
    if (isHeapEmpty(heap)) return FALSE;
    decodeHeapTop(heap, key, val);
    return TRUE;
}
