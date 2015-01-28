/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    FragsClumps.inl Contains inline functions for various operations on clumps and especially frags.

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

/////
//     Inlines for Fragments.
/////

/////
//     Start with basic getter/setters for implicit fields such as endRefOff
/////

// Make up for the fact that these values are not stored in the fragment structure.
STATIC_DECL inline int fragQueryLen (Fragment_t * frag)
{
    return 1 + ((int)frag->endQueryOff - (int)frag->startQueryOff);
}

// While forming fragments, we won't bother to keep both the query length and ref length.
// We don't need to, as they will both be the same.
// Later, we allow them to differ for some fragments as we build the final alignment.
// This routine let's us set the ref length to the query length before they can diverge.
STATIC_DECL inline void setRefLen (Fragment_t * frag)
{
    frag->refLen = fragQueryLen(frag);
}

STATIC_DECL inline SUINT fragRefLen (Fragment_t * frag)
{
    return frag->refLen;
}

STATIC_DECL inline ROFF fragEndRefOff (Fragment_t * frag)
{
    return frag->startRefOff + frag->refLen - 1 ;
}

STATIC_DECL inline void fragSetEndRefOff (Fragment_t * frag, ROFF refOff)
{
    frag->refLen = 1 + refOff - frag->startRefOff;
}

STATIC_DECL inline void fragSetSOs(Fragment_t * frag, QOFF qOff, ROFF rOff)
{
    frag->startQueryOff = qOff;
    frag->startRefOff = rOff;
}

STATIC_DECL inline void fragSetEOs(Fragment_t * frag, QOFF qOff)
{
    frag->endQueryOff = qOff;
    setRefLen(frag);
}

STATIC_DECL inline void addQLenToFront (Fragment_t * frag, int len)
{
    frag->startQueryOff -= len;
}

STATIC_DECL inline void addRLenToFront (Fragment_t * frag, int len)
{
    frag->startRefOff -= len;
    frag->refLen += len;
}

STATIC_DECL inline void addLenToFront (Fragment_t * frag, int len)
{
    addQLenToFront (frag, len);
    addRLenToFront (frag, len);
}

STATIC_DECL inline void addQLenToBack (Fragment_t * frag, int len)
{
    frag->endQueryOff += len;
}

STATIC_DECL inline void addRLenToBack (Fragment_t * frag, int len)
{
    frag->refLen += len;
}

STATIC_DECL inline void addLenToBack (Fragment_t * frag, int len)
{
    addQLenToBack (frag, len);
    addRLenToBack (frag, len);
}

STATIC_DECL inline void subLenFromFront(Fragment_t * frag, int len)
{
    frag->startQueryOff += len;
    frag->startRefOff += len;
    frag->refLen -= len;
}

STATIC_DECL inline void subLenFromBack(Fragment_t * frag, int len)
{
    frag->endQueryOff -= len;
    frag->refLen -= len;
}

STATIC_DECL inline ROFF fragDiag(Fragment_t * frag)
{
    return frag->startRefOff - frag->startQueryOff;
}

STATIC_DECL inline int UINTDiff(UINT uint1, UINT uint2)
{
    // Machinations to avoid unsigned subtraction wrap.
    return (uint1 >= uint2) ? (int)(uint1 - uint2) : -((int)(uint2 - uint1));
}

STATIC_DECL inline UINT absUINTDiff(UINT uint1, UINT uint2)
{
    // Machinations to avoid unsigned subtraction wrap.
    return (uint1 > uint2) ? (uint1 - uint2) : (uint2 - uint1);
}

// These are just synonyms for the above.
// As diagonals are of type ROFF.
STATIC_DECL inline int diagDiff(ROFF diag1, ROFF diag2)
{
    return UINTDiff(diag1, diag2);
}

STATIC_DECL inline UINT absDiagDiff(ROFF diag1, ROFF diag2)
{
    return absUINTDiff(diag1, diag2);
}

/////
//     Various ways of calculating gaps and overlaps between two fragments.
/////

// These are macros instead of inlines because we need them for different data types.
// In C++ we would just have these be overloaded inlines.
// TODO  In future, if QOFF is ever implemented, these unsigned casts will no longer be needed.
#define calcGap(lowOffset, highOffset) ((highOffset > lowOffset) ? ((unsigned)(highOffset - lowOffset) - 1) : 0)
#define calcOverlap(lowOffset, highOffset) ((lowOffset >= highOffset) ? ((unsigned)(lowOffset - highOffset) + 1) : 0)

STATIC_DECL inline int calcMaxOverlap (Fragment_t * frag1, Fragment_t * frag2)
{
    return MAX(calcOverlap(frag1->endQueryOff, frag2->startQueryOff), calcOverlap(fragEndRefOff(frag1), frag2->startRefOff));
}

STATIC_DECL inline int calcQueryGap (Fragment_t * frag1, Fragment_t * frag2)
{
    return calcGap(frag1->endQueryOff, frag2->startQueryOff);
}

STATIC_DECL inline int calcRefGap (Fragment_t * frag1, Fragment_t * frag2)
{
    return calcGap(fragEndRefOff(frag1), frag2->startRefOff);
}

STATIC_DECL inline int calcMinGap (Fragment_t * frag1, Fragment_t * frag2)
{
    return MIN(calcQueryGap(frag1, frag2), calcRefGap(frag1, frag2));
}

STATIC_DECL inline int calcMaxGap (Fragment_t * frag1, Fragment_t * frag2)
{
    return MAX(calcQueryGap(frag1, frag2), calcRefGap(frag1, frag2));
}

/////
//     Inlines related to scoring.
/////

STATIC_DECL inline int calcGapCost(int len, AlignmentArgs_t * AAs)
{
    return (len > 0) ? -(AAs->GOCost + len * AAs->GECost) : 0;
}

// The match count for frags is now always the length.
STATIC_DECL inline int fragMatchCount(Fragment_t * frag)
{
    return frag->refLen;
}

/////
//     Inlines for Super Fragments.
/////

STATIC_DECL inline Fragment_t   * getFragPtr(SFragment_t * SFrag)
{
    return &(SFrag->frag);
}

STATIC_DECL inline EditOpList_t * getEOLPtr(SFragment_t * SFrag)
{
    return &(SFrag->EOList);
}

STATIC_DECL inline BOOL SFragmentListIsEmpty(SFragmentList_t * list)
{
    return (list->head == NULL);
}

/////
//     Inlines for Clumps.
/////

STATIC_DECL inline EditOpList_t * EOLPtr(Clump_t * clump)
{
    return &(clump->EOList);
}

STATIC_DECL inline SFragmentList_t * SFragListPtr(Clump_t * clump)
{
    return &(clump->SFragList);
}

// The status flags for clumps.
#define clumpReversed 0x01  // Indicates clump is aligned to inverse (3' to 5') strand.
#define clumpFormed   0x02  // Is this currently in use?
#define clumpAligned  0x04  // Indicates clump has been aligned.  This implies it will have a single SFrag with an edit op list.
#define clumpScored   0x08  // Indicates clump has successfully been scored.  This implies it has survived splitting and threshold tests.
#define clumpSplit    0x10  // Indicates clump was involved in a split during scoring.
#define clumpPrimary  0x20  // Indicates clump was in set of primary clumps spanning the query.

STATIC_DECL inline BOOL isReversed(Clump_t * clump)
{
    return ((clump->status & clumpReversed) != 0);
}

STATIC_DECL inline BOOL isFormed(Clump_t * clump)
{
    return ((clump->status & clumpFormed) != 0);
}

STATIC_DECL inline BOOL isAligned(Clump_t * clump)
{
    return ((clump->status & clumpAligned) != 0);
}

STATIC_DECL inline BOOL isScored(Clump_t * clump)
{
    return ((clump->status & clumpScored) != 0);
}

STATIC_DECL inline BOOL isSplit(Clump_t * clump)
{
    return ((clump->status & clumpSplit) != 0);
}

STATIC_DECL inline BOOL isPrimary(Clump_t * clump)
{
    return ((clump->status & clumpPrimary) != 0);
}

STATIC_DECL inline BOOL setReversed(Clump_t * clump, BOOL bit)
{
    if (bit)
        return clump->status |= clumpReversed;
    else
        return clump->status &= ~clumpReversed;
}

STATIC_DECL inline BOOL setFormed(Clump_t * clump, BOOL bit)
{
    if (bit)
        return clump->status |= clumpFormed;
    else
        return clump->status &= ~clumpFormed;
}

STATIC_DECL inline BOOL setAligned(Clump_t * clump, BOOL bit)
{
    if (bit)
        return clump->status |= clumpAligned;
    else
        return clump->status &= ~clumpAligned;
}

STATIC_DECL inline BOOL setScored(Clump_t * clump, BOOL bit)
{
    if (bit)
        return clump->status |= clumpScored;
    else
        return clump->status &= ~clumpScored;
}

STATIC_DECL inline BOOL setSplit(Clump_t * clump, BOOL bit)
{
    if (bit)
        return clump->status |= clumpSplit;
    else
        return clump->status &= ~clumpSplit;
}

STATIC_DECL inline BOOL setPrimary(Clump_t * clump, BOOL bit)
{
    if (bit)
        return clump->status |= clumpPrimary;
    else
        return clump->status &= ~clumpPrimary;
}

STATIC_DECL inline Fragment_t * clumpFirstFrag(Clump_t * clump)
{
    return getFragPtr(SFragListPtr(clump)->head);
}

STATIC_DECL inline Fragment_t * clumpLastFrag(Clump_t * clump)
{
    return getFragPtr(SFragListPtr(clump)->tail);
}

STATIC_DECL inline BOOL isEmpty(Clump_t * clump)
{
    return (SFragListPtr(clump)->head == NULL);
}

STATIC_DECL inline QOFF clumpSQO(Clump_t * clump)
{
    return (clumpFirstFrag(clump)->startQueryOff);
}

STATIC_DECL inline QOFF clumpEQO(Clump_t * clump)
{
    return (clumpLastFrag(clump)->endQueryOff);
}

STATIC_DECL inline ROFF clumpSRO(Clump_t * clump)
{
    return (clumpFirstFrag(clump)->startRefOff);
}

STATIC_DECL inline ROFF clumpERO(Clump_t * clump)
{
    return (fragEndRefOff(clumpLastFrag(clump)));
}

STATIC_DECL inline QOFF clumpPlusSQO(Clump_t * clump, QueryState_t * QS)
{
    if (isReversed(clump)) return (QS->queryLen - 1) - clumpEQO(clump);
    else return clumpSQO(clump);
}

STATIC_DECL inline QOFF clumpPlusEQO(Clump_t * clump, QueryState_t * QS)
{
    if (isReversed(clump)) return (QS->queryLen - 1) - clumpSQO(clump);
    else return clumpEQO(clump);
}

STATIC_DECL inline QOFF clumpQueryLen(Clump_t * clump)
{
    return 1 + clumpEQO(clump) - clumpSQO(clump);
}

STATIC_DECL inline QOFF clumpRefLen(Clump_t * clump)
{
    return 1 + clumpERO(clump) - clumpSRO(clump);
}

/////
//     Inline to push a clump onto a clump list.
/////
STATIC_DECL inline Clump_t * pushClump(Clump_t * clump, Clump_t * clumpList)
{
    clump->next = clumpList;
    return clump;
}

/////
//     Inlines for QueryState Accessors.
/////
//     Sneaking this in a file now inaptly named.
/////

STATIC_DECL inline char * getQueryBuf(QueryState_t * QS)
{
    return QS->queryBuf;
}

STATIC_DECL inline char * getQueryCodeBuf(QueryState_t * QS)
{
    return QS->queryCodeBuf;
}

STATIC_DECL inline char * getClumpQueryBuf(QueryState_t * QS, Clump_t * clump)
{
    return (isReversed(clump)) ? QS->reverseBuf : QS->forwardBuf;
}

STATIC_DECL inline char * getClumpQueryCodeBuf(QueryState_t * QS, Clump_t * clump)
{
    return (isReversed(clump)) ? QS->reverseCodeBuf : QS->forwardCodeBuf;
}
