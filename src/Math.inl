/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    Math.inl        Contains inline functions and macros for the genome compression and others.

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

// We need to declare the externs for the various mapping arrays.
// The actual arrays are in Math.c

extern char fourBitCodes[];
extern char fourBitChars[];
extern char fourBitCompCodes[];
extern char fourBitCompBases[];
extern char fourBitEquivalentCodes[];

// Here are all the mapping routine from chars to 4-bit and 2-bit and back again.

STATIC_DECL inline int map8to4(char base)
{
    return fourBitCodes[(int)base];
}

STATIC_DECL inline int map8to2(char base)
{
    int code = map8to4(base);
    ASSERT(code  >= 0 && code <= 3, "Code out of range in map8to2.");
    return code;
}

STATIC_DECL inline int complement2to2(int code)
{
    ASSERT(code  >= 0 && code <= 3, "Code out of range in complement.");
    return fourBitCompCodes[code];
}

STATIC_DECL inline int complement4to4(int code)
{
    ASSERT(code  >= 0 && code <= 15, "Code out of range in complement.");
    return fourBitCompCodes[code];
}

STATIC_DECL inline char complement2to8(int code)
{
    ASSERT(code  >= 0 && code <= 3, "Code out of range in complement.");
    return fourBitCompBases[code];
}

STATIC_DECL inline char complement4to8(int code)
{
    ASSERT(code  >= 0 && code <= 15, "Code out of range in complement.");
    return fourBitCompBases[code];
}

STATIC_DECL inline char complement8to8(char base)
{
    return complement4to8(map8to4(base));
}

STATIC_DECL inline char unmap2to8(int code)
{
    ASSERT(code >= 0 && code <= 3, "Code out of range in unmap.");
    return fourBitChars[code];
}

STATIC_DECL inline char unmap4to8(int code)
{
    ASSERT(code >= 0 && code <= 15, "Code out of range in unmap.");
    return fourBitChars[code];
}

STATIC_DECL inline char getCharFrom4Code (char* seq, ROFF offset)
{
    return unmap4to8(getFrom4Code(seq, offset));
}

// This one tests if two four codes are compatible for a "match" during alignment.
STATIC_DECL inline int codesEqual(int qCode, int rCode)
{
    // Currently unused hack to accept Ns in the query.
    // if (qCode == 4) return TRUE;
    // else return (qCode == rCode);

    // Currently unused hack to accept any non-2bit code from the query.
    // if (qCode > 3) return TRUE;
    // else return (qCode == rCode);

    // Currently unused method using a lookup to determine "equality".
    // Using a allowed equality array is arguably the "correct" thing to do.
    // But deciding what should be considered "equal" is far from clear.
    // The bit arithmetic is to turn two codes into a 1D array index.
    // The required equivalence array is in Match.c
    // return (int)(fourBitEquivalentCodes [(qCode<<4) | rCode]);

    // We now blindly accept the codes if they are equal.
    // This gets "called" a LOT in the inner loop of DP.
    // So, speed is of the essence.
    // Not matching non ACGT chars is not worth the runtime cost for the miniscule error rate induced.
    return (qCode == rCode);
}

// Bit Arrays
extern UBYTE  bits [8];
extern UBYTE  masks [8];

STATIC_DECL inline void clearStart (BitArray_t *bArr, int count)
{
    memset (bArr->contents, 0x00, (count/8) + 1);
}

STATIC_DECL inline void clearAll (BitArray_t *bArr)
{
    memset (bArr->contents, 0x00, bArr->numOfBytes);
}

STATIC_DECL inline BOOL clear (BitArray_t *bArr, int offset)
{
    return (set (bArr, offset, FALSE));
}

// For some crazy reason, the specific inlines for setting are SLOWER than the macro that calls a function to set.
// Even though the latter has to do the funcall AND decide at run time whether it is setting TRUE or FALSE.
// It is not even likely a code bloat issue, as only one setTRUE call appears in all the code.
// And there is NO call to setFALSE on a single byte (clearStart above is used for bulk resetting).
#ifdef NOTNOW
STATIC_DECL inline void setTRUE (BitArray_t *bArr, int offset)
{
    // Note transform to avoid division by 8.
    int byteOffset = offset >> 3;
    int bitOffset = offset & 0x07;
    // bArr->contents[byteOffset] |= bits[bitOffset];
    bArr->contents[byteOffset] |= (1 << bitOffset);
}

STATIC_DECL inline void setFALSE (BitArray_t *bArr, int offset)
{
    // Note transform to avoid division by 8.
    int byteOffset = offset >> 3;
    int bitOffset = offset & 0x07;
    bArr->contents[byteOffset] &= masks[bitOffset];
}

STATIC_DECL inline BOOL get (BitArray_t *bArr, int offset)
{
    // Note transform to avoid division by 8.
    int byteOffset = offset >> 3;
    int bitOffset = offset & 0x07;
    // Array indexing runs faster than creating the mask at runtime.
    // return ((bArr->contents[byteOffset] & bits[bitOffset]) != 0);
    return ((bArr->contents[byteOffset] & (1 << bitOffset)) != 0);
}
#else
#define setTRUE(BA,off) (set(BA, off, TRUE))
#define setFALSE(BA,off) (set(BA, off, FALSE))
#endif

