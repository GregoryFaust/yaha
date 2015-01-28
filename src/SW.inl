/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    SW.inl          Contains inline functions for EditOps used to represent alignments.
                             These are inlined for speed to support the unusual storage for EditOps.
                             See SW.cpp for more details.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

// We need this to allow this file to be included in both C and C++ code.

#ifdef __cplusplus
#define STATIC_DECL
#else
#define STATIC_DECL static
#endif

#ifdef __cplusplus
extern "C"
{
EditOpPtr allocEditOp (EditOpCode code, SUINT len, EditOpList_t * list);
void addEditOpToList(EditOpList_t * list, EditOpPtr prev, EditOpPtr next, EditOpPtr add);
}
#else
EditOpPtr allocEditOp (EditOpCode code, SUINT len, EditOpList_t * list);
void addEditOpToList(EditOpList_t * list, EditOpPtr prev, EditOpPtr next, EditOpPtr add);
#endif

STATIC_DECL inline QueryState_t * getEOLQSPtr(EditOpList_t * list)
{
    return list->QS;
}

STATIC_DECL inline EditOp_t * getEOLEOArrayPtr(EditOpList_t * list)
{
    return getEOLQSPtr(list)->EditOpArray;
}

STATIC_DECL inline char EditOpCode2char (EditOpCode code) { return code; }

STATIC_DECL inline EditOp_t * nodeFromPtr(EditOpPtr ptr, EditOpList_t * list)
{
    return (EditOp_t *)(&(getEOLEOArrayPtr(list)[ptr]));
}

STATIC_DECL inline EditOpPtr ptrFromNode(EditOp_t * node, EditOpList_t * list)
{
    return (EditOpPtr)(node - getEOLEOArrayPtr(list));
}

STATIC_DECL inline BOOL EOLisEmpty(EditOpList_t * list)
{
    return (list->head == EditOpPtrNULL);
}

STATIC_DECL inline void addEditOpToFront(EditOpList_t * list, EditOpCode code, SUINT len) __attribute__((always_inline));
STATIC_DECL inline void addEditOpToFront(EditOpList_t * list, EditOpCode code, SUINT len)
{
    EditOpPtr add = allocEditOp(code, len, list);
    addEditOpToList(list, EditOpPtrNULL, list->head, add);
}

STATIC_DECL inline void addEditOpToBack(EditOpList_t * list, EditOpCode code, SUINT len) __attribute__((always_inline));
STATIC_DECL inline void addEditOpToBack(EditOpList_t * list, EditOpCode code, SUINT len)
{
    EditOpPtr add = allocEditOp(code, len, list);
    addEditOpToList(list, list->tail, EditOpPtrNULL, add);
}

STATIC_DECL inline void initEditOpList(EditOpList_t * list, QueryState_t * QS)
{
    list->head = list->tail = EditOpPtrNULL;
    list->QS = QS;
}

#ifdef NOTNOW
STATIC_DECL inline EditOp_t * nodeFromPtrNull(EditOpPtr ptr)
{
    // return ptr ? (&((*(list->EOArrayPtr))[ptr])) : NULL;
    return ptr ? (EditOp_t *)(&(EOLArray[ptr])) : NULL;
}

STATIC_DECL inline int isValidEONode(EditOp_t * node)
{
    return (node != NULL);
}

STATIC_DECL inline EditOp_t * getFirstEditOp(EditOpList_t * list)
{
    return nodeFromPtrNull(list->head);
}

STATIC_DECL inline EditOp_t * getLastEditOp(EditOpList_t * list)
{
    return nodeFromPtrNull(list->tail);
}

STATIC_DECL inline EditOp_t * getNextEditOp(EditOp_t * item)
{
    return nodeFromPtrNull(item->next);
}

STATIC_DECL inline EditOp_t * getPrevEditOp(EditOp_t * item)
{
    return nodeFromPtrNull(item->prev);
}

#else

STATIC_DECL inline int isValidEONode(EditOp_t * node, EditOpList_t * list)
{
    return (node != getEOLEOArrayPtr(list));
}

STATIC_DECL inline EditOp_t * getFirstEditOp(EditOpList_t * list)
{
    return nodeFromPtr(list->head, list);
}

STATIC_DECL inline EditOp_t * getLastEditOp(EditOpList_t * list)
{
    return nodeFromPtr(list->tail, list);
}

STATIC_DECL inline EditOp_t * getNextEditOp(EditOp_t * item, EditOpList_t * list)
{
    return nodeFromPtr(item->next, list);
}

STATIC_DECL inline EditOp_t * getPrevEditOp(EditOp_t * item, EditOpList_t * list)
{
    return nodeFromPtr(item->prev, list);
}

// Define this as a macro.
#define forAllEditOpsInList(name, list) for (EditOp_t * name = getFirstEditOp(list); isValidEONode(name, list); name = getNextEditOp(name, list))

#endif
