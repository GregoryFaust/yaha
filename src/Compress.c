/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    Compress.c     Compress .fna or .fasta file into a .nib2 file.
                          UnCompress .nib2 file into a .fasta file.


    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include "Math.h"
#include "Math.inl"

/////
//     .nib2 file format.
/////

//  The .nib2 file format is a variation on the .2bit file format.
//  But still keeps 4 bits per base to represent the non-ACGT and N codes
//     that are allowed in FASTA format files.
//
//  Unlike .2bit, we will put all the meta data at the start of the file
//     instead of interspersed with the sequence data.
//  The base sequences will be rounded up to the nearest 4 byte boundary,
//     with the code for 'X' used as filler as needed.
//  We will purposely make the early information 4 bytes at a time
//     to ease memory mapped IO.  In below, no length implies 4 bytes.
//
//  The file header will contain:
//  signature:            0x01020304
//  version:              The first two bytes are the major version, and lower two bytes the minor version.
//  byte offet of bases   The byte offset from start of file of the first base code.  This will be 32 bit aligned.
//  sequenceCount:        The number of sequences in the file.
//
//  For each sequence, we will then have:
//  sequence start        The byte offset from the start of the base sequence information.
//  sequence length       The count of the number of bases (not bytes) in the sequence.
//  seq name offset(16)   The offset of sequence name starting byte from start of file.
//  seq name length(16)   The length of the sequence name in bytes.
//
//  2/28/2013 Modified file structure (now "version 0.2");
//  Change seq name offset and seq name length to be full ints.
//  We don't need the extra size in the name length, but we do in the offset.
//  So, increase name length anyway to maintain seq block to be a multiple of 32 bits.
//  We will make it so that the code uses the file format version number to read in old as well as new files.
//
//  Number of Mask Blocks The number of mask blocks in the file.
//       For now we will continue to ignore masked blocks.  Half the mouse genome is masked!!
//       So, just put in zero for the number and go on.
//
//  We will then have all the sequence names with no spaces or nulls or line feeds.
//  We will then round to the nearest 32 bits padded with null chars as needed.
//  We will then start the base sequences as described above.

/////
//     Read in the Base Sequence information from a nib2 file.
/////

#define nib2HeaderSize 16
// Since this now depends on the file version, it can no longer be a compile time constant.
// #define BSBlockSize    12
#define MSBlockSize     0
#define MSHeaderSize    4
#define nib2marker     0x01020304

BaseSequences_t * loadBaseSequences(char * gfilePtr)
{

    // Get an integer pointer to the header information.
    UINT * gheadPtr = (UINT *) gfilePtr;
    // Check for expected header information.
    int versionNumber = gheadPtr[1];
    if (gheadPtr[0] != nib2marker || (versionNumber != 1 && versionNumber != 2))
    {
        fatalError("Input nib2 file bad header format.");
    }

    // Determine which file version we are dealing with.
    int BSBlockSize = 12;
    if (versionNumber == 2)
    {
        BSBlockSize = 16;
    }

    // Get the number of base sequences.
    int baseSeqCount = gheadPtr[3];
    // Make the return value.
    BaseSequences_t * BSs = makeBaseSequences(baseSeqCount);
    // Get a pointer to the start of the base array.
    BSs->basePtr = gfilePtr + gheadPtr[2];
    // Now get the start of the sequence name block.
    // TODO This calculation should take into account any masked blocks.
    char * nameStart = gfilePtr + nib2HeaderSize + (BSBlockSize * baseSeqCount) + MSHeaderSize +
        gheadPtr[nib2HeaderSize/4 + BSBlockSize/4 * baseSeqCount] * MSBlockSize;

    // Read in the Base Sequences data structure.
    gheadPtr += 4;
    UINT totalBaseCount = 0;
    for (int i=0; i<baseSeqCount; i++)
    {
        BaseSequence_t * BS = makeBaseSequence();
        BS->startingOffset = gheadPtr[0];
        BS->length = gheadPtr[1];
        totalBaseCount += BS->length;
        if (versionNumber == 1)
        {
            UINT nameInfo = gheadPtr[2];
            BS->name = nameStart + (uint16_t)(nameInfo >> 16);
            BS->nameLen = nameInfo & (uint16_t)0xFFFF;
            gheadPtr += 3;
        }
        else
        {
            // Version 2
            BS->name = nameStart + gheadPtr[2];
            BS->nameLen = gheadPtr[3];
            gheadPtr += 4;
        }
        addBaseSequence(BSs, BS);
    }
    BSs->totalBaseCount = totalBaseCount;

    return BSs;
}

/////
//     Output the Base Sequence information to a nib2 file.
/////

void outputBaseSequences(BaseSequences_t * BSs, FDES ofp)
{
    // Now we need to output the header for all the base sequence information.
    // Start by getting the total number of bytes in all of the sequence names
    //   rounded up to an even multiple of 4 bytes.
    int totNameLen = 0;
    int seqCount = BSs->curCount;
    for (int i=0; i<seqCount; i++)
    {
        totNameLen += BSs->sequences[i]->nameLen;
    }

    totNameLen = roundUpPow2(totNameLen, 2);
    int preambleSize = 20 + (16 * seqCount) + totNameLen;

    // Output the fixed file header.
    UINT value[4];
    char * valuePtr = (char *) &value;
    value[0] = nib2marker;
    value[1] = 2;
    value[2] = preambleSize;
    value[3] = seqCount;
    blockWrite(ofp, valuePtr, 16);

    // Output the sequence header information.
    int nameOffset = 0;
    for (int i=0; i<seqCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        value[0] = BS->startingOffset;
        value[1] = BS->length;
        value[2] = nameOffset;
        value[3] = BS->nameLen;
        blockWrite(ofp, valuePtr, 16);
        nameOffset += BS->nameLen;
    }

    // Output the mask information.
    // For now, just output zero as the count of mask blocks.
    value[0] = 0;
    blockWrite(ofp, valuePtr, 4);

    // Output the sequence name information
    for (int i=0; i<seqCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        blockWrite(ofp, BS->name, BS->nameLen);
    }
    // Pad out to multiple of 32 bytes with null chars.
    value[0] = 0;
    blockWrite(ofp, valuePtr, totNameLen - nameOffset);
}

/////
//     Compress Fasta file to .nib2 file
/////

// Define a little Helper function for finalizing a baseSequence.
// Used in compressFile.
ROFF finalizeBaseSequence(BaseSequence_t * BS, int baseCount, char * ofilePtr, ROFF ofileOffset)
{
    if (BS != NULL)
    {
        // We finally know the previous sequence length, so remember it.
        BS->length = baseCount;
        // Pad any extra base positions with 'X' codes.
        int Xcode = map8to4('X');
        if (baseCount & 0x1)
        {
            ofilePtr[ofileOffset] += Xcode;
            ofileOffset += 1;
        }
        // Let's round the output offset up to the nearest 32 bit boundary.
        int nextOfileOffset = roundUpPow2(ofileOffset, 2);
        for (int i=ofileOffset; i<nextOfileOffset; i++) ofilePtr[i] = (Xcode << 4) + Xcode;
        return nextOfileOffset;
    }
    return 0;
}

void compressFile(AlignmentArgs_t * AAs)
{
/////
//     Open the files.
/////

    char *gfilePtr;
    FSIZE gfileSize;
    FDES gfp = openForRead(AAs->gfileName, &gfileSize, &gfilePtr, TRUE);

    FDES ofp = openForWrite(AAs->ofileName);
    // Estimate the size we will need to store the compressed sequences.
    FSIZE ofileSize = gfileSize/2;
    char *ofilePtr = blockMalloc(ofileSize, AAs->ofileName);

    // Build the base sequences structure.
    // Start with 100, and it will automatically expand as needed.
    BaseSequences_t * BSs = makeBaseSequences(100);
    BaseSequence_t * BS = NULL;

/////
//     Main Read Loop.
/////

    // Place to put the current base code.
    int code;
    // Number of bases so far in the current sequence.
    int baseCount = 0;
    // Total number of non-ACGT bases in file.  Only used for stats at end.
    int offCodeCount = 0;
    // Offset in bytes in the input file.
    ROFF gfileOffset = 0;
    // Offset in bytes in the output file from the start of the compressed base sequences.
    ROFF ofileOffset = 0;
    while (gfileOffset < gfileSize)
    {
        // Throw away newlines.
        if (gfilePtr[gfileOffset] >= (char)0 && gfilePtr[gfileOffset] <= (char)31)
        {
            gfileOffset += 1;
            continue;
        }
        // Check for the start of a new sequence.
        if (gfilePtr[gfileOffset] == '>')
        {
            // If this is not the first sequence, then finalize data on previous sequence.
            ofileOffset = finalizeBaseSequence(BS, baseCount, ofilePtr, ofileOffset);
            baseCount = 0;
            // Make new base sequence.
            BS = makeBaseSequence();
            addBaseSequence(BSs, BS);
            // File in the name information.
            int newLine = findChar(gfilePtr, gfileOffset, '\n');
            BS->name = gfilePtr+gfileOffset+1;
            BS->nameLen = newLine - 1;
            gfileOffset += newLine + 1;
            // Check to see if there is a space in the sequence name, and if so chop there.
            for (int i=0; i<BS->nameLen; i++)
            {
                if (BS->name[i] == ' ')
                {
                    BS->nameLen = i;
                    break;
                }
            }
            // Store the starting offset for sequence.
            // We don't know the length yet, so update length at start of next seq or EOF.
            BS->startingOffset = ofileOffset;
            continue;
        }
        // We have a normal base code.
        code = map8to4(gfilePtr[gfileOffset]);
        gfileOffset += 1;
        if (code > 3)offCodeCount += 1;
        baseCount += 1;
        // Need different treatment for even and odd cases.
        if (baseCount & 0x1) ofilePtr[ofileOffset] = code << 4;
        else
        {
            ofilePtr[ofileOffset] += code;
            ofileOffset += 1;
        }
    }

    // We've hit the end of file.
    // Finalize the last BS, and pad to a round number of bytes.
    ofileOffset = finalizeBaseSequence(BS, baseCount, ofilePtr, ofileOffset);

    if (AAs->verbose) fprintf(stderr, "Compressed %d reference sequences from %s into %s.\n",
                              BSs->curCount, AAs->gfileName, AAs->ofileName);

    // Output the nib2 file header.
    outputBaseSequences(BSs, ofp);

    // Free the Base Sequence structures
    disposeBaseSequences(BSs);

/////
//     Output the array of bases and clean up
/////

    // Now output and free the array containing the base sequences.
    blockWrite(ofp, ofilePtr, ofileOffset);
    blockFree(ofilePtr, ofileSize);

    // Close the files.
    closeForRead(gfp, gfileSize, gfilePtr);
    closeForWrite(ofp);

}

/////
//     Uncompress a .nib2 file to a fasta file.
/////

#define FASTALineLength 50

void uncompressFile(AlignmentArgs_t * AAs)
{

/////
//    Open the input file and read in base sequence information.
/////

    char *gfilePtr;
    FSIZE gfileSize;
    FDES gfp = openForRead(AAs->gfileName, &gfileSize, &gfilePtr, TRUE);

    // Read in the .NIB2 header information.
    BaseSequences_t * BSs = loadBaseSequences(gfilePtr);
    if (AAs->verbose) fprintf(stderr, "Read in %d reference sequences from %s.\n", BSs->curCount, AAs->gfileName);

    // Calculate the size of the name sequences.
    int nameSize = 0;
    for (int i=0; i<BSs->curCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        // add 2 to account for the > and the newline.
        nameSize += BS->nameLen + 2;
    }

/////
//     Write output file
/////

    // We are ready to open the output file.
    // We will have to guess at the right size.
    // FASTA files seem to have a newline every 60 chars.
    // Add 1000 slop just because.....
    FSIZE ofileLen = (BSs->totalBaseCount/FASTALineLength)*(FASTALineLength+1) + nameSize + 1000;
    FDES ofp = openForWrite(AAs->ofileName);
    char * ofilePtr = blockMalloc(ofileLen, AAs->ofileName);

    ROFF ofileOffset = 0;
    for (int i=0; i<BSs->curCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        ofilePtr[ofileOffset++] = '>';
        for (int j=0; j<BS->nameLen; j++) ofilePtr[ofileOffset++] = BS->name[j];
        ofilePtr[ofileOffset++] = '\n';
        int newLineCount = 0;
        char * iPtr = BSs->basePtr + BS->startingOffset;
        for (int j=0; j<BS->length; j++)
        {
            if (newLineCount == FASTALineLength)
            {
                ofilePtr[ofileOffset++] = '\n';
                newLineCount = 0;
            }
            ofilePtr[ofileOffset++] = getCharFrom4Code(iPtr, j);
            newLineCount += 1;
        }
        if (newLineCount != 0)ofilePtr[ofileOffset++] = '\n';
    }

    blockWrite(ofp, ofilePtr, ofileOffset);

    // Close the files and other cleanup.
    disposeBaseSequences(BSs);
    closeForRead(gfp, gfileSize, gfilePtr);
    closeForWrite(ofp);

}
