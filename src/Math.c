/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    Math.c          Contains some of the common algorithmic code such as bitarrays.
                             Also contains the basic DNA compress and decompression helpers.
                             Although since many of these are inlines, also see Math.inl
                             Contains global read-only arrays for above.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include "Math.h"
#include "Math.inl"

// An empty function that can be called to do nothing,
//     but let compiler know it can't optimize away a variable is used as an argument to it!!!
void emptyFunction (int arg)
{
}

// Some basic integer math routines.
UINT64 iexp64(int base, int power)
{
    UINT64 retval = 1;
    for (int i=0; i<power; i++) retval *= base;
    return retval;
}

UINT iexp(int base, int power)
{
    UINT retval = 1;
    for (int i=0; i<power; i++) retval *= base;
    return retval;
}

int ilog(int base, UINT num)
{
    int retval=0;
    while (num > 1)
    {
        num = num/base;
        retval += 1;
    }
    return retval;
}

//////////////
//	Bit Array implementation.  Use unsigned bytes for the storage.
//	Bits start at the left and count to the right.
//	As usual for C, length is 1 based, but offset is 0 based.
//////////////
//      In an attempt to get things to go fast, I have replaced /8 with >>3
//      And likewise %8 with &0x07.  It is unclear if these
//      are faster than the divide operations or not, or whether the compiler
//      would do the right optimization since the 8 is always a constant.
//////////////

// Is it faster to access these prebuilt arrays?
// Or to just build the value we need each time from constants
// in the instruction stream via shifts?
UBYTE bits  [8] = {0x80,0x40,0x20,0x10,0x08,0x04,0x02,0x01};
UBYTE masks [8] = {0x7F,0xBF,0xDF,0xEF,0xF7,0xFB,0xFD,0xFE};

BitArray_t * makeBitArray(int numOfElements)
{
    BitArray_t *bArr = (BitArray_t *)(malloc(sizeof(BitArray_t)));
    bArr->numOfBits = numOfElements;
    bArr->numOfBytes = numOfElements >> 3;
    if ((bArr->numOfBits && 0x7) != 0) bArr->numOfBytes += 1;
    // Make sure this is all set to zero before we start.
    bArr->contents = (UBYTE *)(calloc(bArr->numOfBytes, 1));
    return bArr;
}

void disposeBitArray (BitArray_t * bArr)
{
    free(bArr->contents);
    free(bArr);
}

void setAll (BitArray_t *bArr, BOOL value)
{
    if (value) memset (bArr->contents, 0xff, bArr->numOfBytes);
    else       memset (bArr->contents, 0x00, bArr->numOfBytes);
}

#ifdef NOTNOW
// This is the only case that gets called, so optimize it.
// This is now defined in Math.inl
void setTRUE (BitArray_t *bArr, int offset)
{
    // Note transform to avoid division by 8.
    int byteOffset = offset >> 3;
    int bitOffset = offset & 0x07;
    bArr->contents[byteOffset] |= bits[bitOffset];
}
#endif

BOOL set (BitArray_t *bArr, int offset, BOOL value)
{
    ASSERT (offset >= 0 && offset < bArr->numOfBits, "Invalid offset in BitArray");

    // Note transform to avoid division by 8.
    int byteOffset = offset >> 3;
    int bitOffset = offset & 0x07;

    if (value) bArr->contents[byteOffset] |= bits[bitOffset];
    else       bArr->contents[byteOffset] &= masks[bitOffset];

    return (value);
}

BOOL get (BitArray_t *bArr, int offset)
{
    // Note transform to avoid division by 8.
    int byteOffset = offset >> 3;
    int bitOffset = offset & 0x07;

    BOOL retval = ((bArr->contents[byteOffset] & bits[bitOffset]) != 0);
    return (retval);
}

/////
//    Functions to encode and decode the base compressions.
//    Most of these are inlines, so look in Math.inl as well.
/////

// This covers all 128 ascii codes so that we need only do a lookup to get the 4-bit code.
// All nonsense codes now map to X instead of -1.
// TODO Consider moving X to be code 15.
char  fourBitCodes[]       = { 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
                               14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
                               14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
                               14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
//                              @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
                               14,  2,  5,  1,  6, 14, 14,  3,  7, 14, 14,  8, 14,  9,  4, 14,
//                              P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   _
                               14, 14, 10, 11,  0,  0, 12, 13, 14, 15, 14, 14, 14, 14, 14, 14,
//                              `   a   b   c   d   e   f   g   h   i   j   k   l   m   n   o
                               14,  2,  5,  1,  6, 14, 14,  3,  7, 14, 14,  8, 14,  9,  4, 14,
//                              p   q   r   s   t   u   v   w   x   y   z   {   |   }   ~ DEL
                               14, 14, 10, 11,  0,  0, 12, 13, 14, 15, 14, 14, 14, 14, 14, 14 };


char  fourBitChars[]       = {'T', 'C', 'A', 'G', 'N', 'B', 'D', 'H', 'K', 'M', 'R', 'S', 'V', 'W', 'X', 'Y'};
char  fourBitCompCodes[]   = { 2,   3,   0,   1,   4,  12,   7,   6,   9,   8,  15,  11,   5,  13,  14,  10 };
char  fourBitCompBases[]   = {'A', 'G', 'T', 'C', 'N', 'V', 'H', 'D', 'M', 'K', 'Y', 'S', 'B', 'W', 'X', 'R'};

// This is the BOOL equivalent way to map one 4-bit code to all possible 2 bit codes.
// For now we will only match codes that represent 1 or 2 bases, but not those that represent 3 or 4.
//                                 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
char  fourBitEquivalentCodes [] = {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1,  //  0 T
                                   0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1,  //  1 C
                                   0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0,  //  2 A
                                   0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0,  //  3 G
                                   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  //  4 N = All
                                   0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  //  5 B = G, T, C
                                   0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,  //  6 D = G, A, T
                                   0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,  //  7 H = A, C, T
                                   1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,  //  8 K = G, T
                                   0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  //  9 M = A, C
                                   0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,  // 10 R = G, A
                                   0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,  // 11 S = G, C
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  // 12 V = G, C, A
                                   1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,  // 13 W = A, T
                                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,  // 14 X = None
                                   1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}; // 15 Y = T, C

// Unpack a 4 bit value from a byte array (used for reading offset in .nib format)
char getFrom4Code(char* seq, ROFF offset)
{
    int byteOffset = offset >> 1;
    int low = offset & 0x1;
    UBYTE code = (UBYTE)seq[byteOffset];
    if (!low) code = (code >> 4) & 0xF;
    else code = code & 0xf;
    return code;
}

// This takes a hash code and recreates the code it came from.
// This is only used in debugging.
// seq had better be long enough to fit the code plus a 0 terminator!!
void decode (UINT code, char* seq, int len)
{
    int nextBase;
    for (int i=0; i<len; i++)
    {
        nextBase = code & 0x3;
        seq[len-(i+1)] = unmap2to8(nextBase);
        code = code >> 2;
    }
    // Null terminate the string.
    seq[len] = 0;
}

// Output the string of chars from the references
void output4Seq(FILE * stream, char* seq, ROFF offset, int len)
{
    for (int k=0; k<len; k++)
    {
        fprintf(stream, "%c", getCharFrom4Code(seq, offset+k));
    }
}

// Check the chars from the reference against the hash, and output if different.
// Output the string of chars from the references
void check4Seq(FILE * stream, char * hashSeq, char* seq, ROFF offset, int len)
{
    BOOL misMatch = FALSE;
    for (int k=0; k<len; k++)
    {
        if (hashSeq[k] != getCharFrom4Code(seq, offset+k)) misMatch = TRUE;
    }
    if (misMatch)
    {
        fprintf(stream, "Hash Sequence=\"%s\" ", hashSeq);
        fprintf(stream, "@R[%10d]=", offset);
        output4Seq(stream, seq, offset, len);
        fprintf(stream, "\n");
    }
}


/////
//    Random Number Generation
/////

// The following is due to George Marsaglia and has been adapted below.
/*
static unsigned long
x=123456789,y=362436069,z=521288629,w=88675123,v=886756453;
// replace defaults with five random seed values in calling program
unsigned long xorshift(void)
{
     unsigned long t;
     t=(x^(x>>7)); x=y; y=z; z=w; w=v;
     v=(v^(v<<6))^(t^(t<<13)); return (y+y+1)*v;
}
*/

randState_t * makeRandState()
{
    return (randState_t *) malloc (sizeof(randState_t));
}

// This will use the default init values from George Marsaglia
void initRandStateDefault(randState_t * rs)
{
    UINT init[5] = {123456789,362436069,521288629,88675123,886756453};
    initRandState(rs, init);
}

// This will use user specified init values.
void initRandState(randState_t * rs, UINT * instate)
{
    for (int i=0; i<5; i++) rs->state[i] = instate[i];
}

void disposeRandState(randState_t * rs)
{
    free(rs);
}

UINT getRandBits(randState_t * rs)
{
    UINT t;
    t = (rs->state[0]^(rs->state[0]>>7));
    rs->state[0] = rs->state[1];
    rs->state[1] = rs->state[2];
    rs->state[2] = rs->state[3];
    rs->state[3] = rs->state[4];
    rs->state[4] = (rs->state[4]^(rs->state[4]<<6))^(t^(t<<13));
    return (rs->state[1]+rs->state[1]+1)*rs->state[4];
}

// We need a function to return a float, and use it to return ints in a defined range.
// For sampling, we can use Floyd's algorithm, but to keep from using a hash table of values already selected,
//     we can use a bool array or bit array to keep track of the ones already hit.
double getRandDouble(randState_t * rs)
{
    return (double)(getRandBits(rs))/((double)(UINT_MAX_VALUE) + 1.0);
}

// This returns an unsigned int
UINT getRandUInt(randState_t * rs, UINT start, UINT end)
{
    return start + (UINT)(getRandDouble(rs)*(end-start));
}

// This is a modified version of Floyd's Algorithm.
// Instead of a HashSet we will use a BOOL array to keep track of which items have already been chosen.
// We will also "pick" one to keep or throw away
//     depending on whether or not we want to keep more or less than half of the input items.
void getRandSample(randState_t * rs, UINT * input, int inLen, UINT * output, int outLen)
{
    // Allocate a BOOL array to store whether or not an input item has already been marked.
    // Use memset to initialize to all FALSE.
    BOOL * marked = (BOOL *)malloc(inLen*sizeof(BOOL));
    memset(marked, FALSE, inLen);

    // Decide to mark the ones to keep or the ones not to keep.
    BOOL keepMarked = TRUE;
    int selectNum = outLen;
    if (outLen > inLen/2)
    {
        keepMarked = FALSE;
        selectNum = inLen-outLen;
    }
    // First select the items to mark.
    for (int i=inLen-selectNum; i<inLen; i++)
    {
        UINT pos = getRandUInt(rs, 0, i+1);
        if (marked[pos])
            marked[i] = TRUE;
        else
            marked[pos] = TRUE;
    }
    // Now put them into the output array.
    // Note that two passes are essential to keep the order the same as in the input,
    //     which is critical to our application.
    // We use two separate loops for keeping marked vs. unmarked, as we don't want to make the check for each item.
    int outOff = 0;
    for (int i=0; i<inLen; i++)
    {
        if (marked[i] == keepMarked)
        {
            output[outOff] = input[i];
            outOff += 1;
        }
    }
    if (outOff != outLen) fprintf(stderr, "Unmatched counts %d %d.\n", outOff, outLen);
    free(marked);
}
