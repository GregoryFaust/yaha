/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    Index.c         Build an index file for the input reference genome using the specified parameters.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include <assert.h>
#include "Math.h"

// Turn this on for a dump of counts per kmer, and counts per reference location
// #define FULLSTATS

// Predeclare some function signatures for function definitions appearing later in the file.
void outputCountStatistics(UINT HTSize, UINT * counts, AlignmentArgs_t * AAs);
void validateHT(UINT HTsize, UINT totalMatches, char *gbasePtr, ROFF *refPtr, UINT *sOffsPtr, int wordLen, int maxHits);

// Putting this inline here makes things much faster.
// The actual hash code is "returned" in the hash "out" parameter.
// The return value is zero on success, or 1 past the offending offset if fail.
inline int generateMatches4to2Fast(char * seq, ROFF offset, int length, UINT* hash)
{
    UINT retval = 0;
    for (int i=0; i<length; i++)
    {
        UBYTE code = (getFrom4Code(seq, offset+i));
        if (code > 3) return offset+i+1; // FALSE;
        retval = (retval << 2) + code;
    }
    *hash = retval;
    return 0; //TRUE;
}

// We will do this in two passes.
// The first pass just counts up the number of times each word appears in the reference.
// This information allows us to then build the rest of the data structure with the correct size information,
//   and to fill in that structure in the second pass.
void indexFile (AlignmentArgs_t * AAs)
{
    // First open the genome file.
    char *gfilePtr;
    FSIZE gfileSize;
    FDES gfp = openForRead(AAs->gfileName, &gfileSize, &gfilePtr, TRUE);

    // Read in the .nib2 file header.
    BaseSequences_t * BSs = loadBaseSequences(gfilePtr);
    char * gbasePtr = BSs->basePtr;
    normalizeBaseSequences(BSs);
    if (AAs->verbose) fprintf(stderr, "Read in %d reference sequences from %s.\n", BSs->curCount, AAs->gfileName);
    AAs->BSs = BSs;

#ifdef FULLSTATS
    int xfileNameLen = strlen(AAs->xfileName);
    char * kmerFileName = (char *)malloc(xfileNameLen+12);
    strcpy(kmerFileName, AAs->xfileName);
    strcpy(kmerFileName+xfileNameLen, ".kmerstats");
    char * refFileName = (char *)malloc(xfileNameLen+12);
    strcpy(refFileName, AAs->xfileName);
    strcpy(refFileName+xfileNameLen, ".refstats");
    fprintf(stderr, "Doing full stats production into %s and %s.\n", kmerFileName, refFileName);
    char * decodeBuf = (char *)malloc(AAs->wordLen+1);
#endif

////////////
//     First Pass
///////////
    // Cache these to make the loops go faster.
    int wordLen = AAs->wordLen;
    int skipDist = AAs->skipDist;

    // Since we can't yet allocate this in the output file, we will have to alloc on the heap to start.
    // TODO Perhaps we could allocate this and then do a realloc?  That has a chance to save a copy of this data later on.
    UINT   HTsize = iexp(4, wordLen);
    FSIZE  HTcountsSize = HTsize * sizeof(UINT);
    UINT   *counts = (UINT *)malloc(HTcountsSize);
    assert(counts);
    memset(counts, 0, HTcountsSize);

    // Calculate some masks for use later.
    UINT hashMask = ((UINT)0xFFFFFFFF) >> (32 - 2*wordLen);
    int  hashSaveLength = wordLen - skipDist;
    // fprintf(stderr, "HTsize is %u, and hashMask is %x.\n", HTsize, hashMask);

    // Iterate through each sequence.
    UINT hashCode = 0;
    UINT partialHashCode = 0;
    for (int i=0; i<BSs->curCount; i++)
    {
        // Iterate through each seed in the sequence.
        BaseSequence_t * BS = BSs->sequences[i];
        ROFF baseOff = BS->startingOffset;
        ROFF endingOffset = BS->startingOffset + BS->length - wordLen;

        int badOffset = generateMatches4to2Fast(gbasePtr, baseOff, wordLen, &hashCode);
        while (TRUE)
        {
            if (badOffset != 0)
            {
                // Get rid of all the bad ones in a row.
                while (getFrom4Code(gbasePtr, badOffset) > 3) badOffset += 1;
                // Renormalize to an offset that is a multiple of skipDist.
                baseOff = (((badOffset + (skipDist - 1))/skipDist) * skipDist);
                if (baseOff > endingOffset) break;
                // Try to load entire word.
                badOffset = generateMatches4to2Fast(gbasePtr, baseOff, wordLen, &hashCode);
                continue;
            }
            // The previous read was good, so add it to the counts.
            counts[hashCode] += 1;
            // and just get the next skipDist worth
            baseOff += skipDist;
            if (baseOff > endingOffset) break;
            badOffset = generateMatches4to2Fast(gbasePtr, baseOff + hashSaveLength, skipDist, &partialHashCode);
            // and tack it on to the rest from the previous get.
            hashCode = ((hashCode << (skipDist*2)) | partialHashCode) & hashMask;
        }
    }

#ifndef COMPILE_USER_MODE
    // The statistics only depend on the first pass information.
    if (AAs->verbose) outputCountStatistics(HTsize, counts, AAs);
#endif

#ifdef FULLSTATS
    FILE * kmerFile = openForPrint(kmerFileName);
    for (int i=0; i<HTsize; i++)
    {
        decode(i, decodeBuf, AAs->wordLen);
        fprintf(kmerFile, "%d\t%s\t%u\n", i, decodeBuf, counts[i]);
    }
    closeForPrint(kmerFile);
    free(kmerFileName);
    FILE * refFile = openForPrint(refFileName);
#endif

////////////
//     Second Pass
///////////
    // Find out how many matches we have, and allocate the new data structures.
    // Here is where we can cut off counts over the maximum.
    ROFF totalMatches = 0;
    for (int i=0; i<HTsize; i++)
    {
#ifndef RANDOM_SAMPLE_MAXHITS
        // If we are randomly sampling, we will need these high counts to remain for now.
        if (counts[i] > AAs->maxHits) counts[i] = 0;
#endif
        totalMatches += counts[i];
    }
    FSIZE xfileTotalSize = 4*4 + ((HTsize+1) * sizeof(ROFF)) + (totalMatches * sizeof(ROFF));

    // Now that we know the total size, we can open the file for write.
    FSIZE xfileAllocSize = xfileTotalSize;
    FDES xfp = openForWrite(AAs->xfileName);
    char *xfilePtr = blockMalloc(xfileAllocSize, AAs->xfileName);
    assert(xfilePtr);

    // Copy in the first few values in the file.
    // The word size, and the number of total references.
    UINT * ints = (UINT *)xfilePtr;
    // Start with a header that contains format version and creation parameters.
    ints[0] = CURRENT_INDEX_FILE_VERSION;
    ints[1] = wordLen;
    ints[2] = AAs->maxHits;
    ints[3] = totalMatches;

    // Account for the header
    FSIZE offset = 4*4;
    // Set up the pointers.
    UINT * sOffsPtr = (UINT *)(xfilePtr + offset);
    offset += (HTsize+1)*sizeof(UINT);
    ROFF * refPtr = (ROFF *)(xfilePtr + offset);

    // Go through the counts and set up the starting offsets.
    offset = 0;
    for (int i=0; i<HTsize; i++)
    {
        sOffsPtr[i] = offset;
        offset += counts[i];
    }
    // Add in the final "extra" offset.
    // This allows us to calculate the counts by subtracting offsets even for the last kmer.
    sOffsPtr[HTsize] = totalMatches;

    // Now we can re-use the counts array to store the number of counts seen so far in the second pass.
    memset(counts, 0, HTcountsSize);

    // Now we are ready to loop through the reference genome the second time.
    // Iterate through each sequence.
    for (int i=0; i<BSs->curCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        // Iterate through each seed in the sequence.
        ROFF baseOff = BS->startingOffset;
        ROFF endingOffset = BS->startingOffset + BS->length - wordLen;

        int badOffset = generateMatches4to2Fast(gbasePtr, baseOff, wordLen, &hashCode);
        while (TRUE)
        {
            if (badOffset != 0)
            {
                // Get rid of all the bad ones in a row.
                while (getFrom4Code(gbasePtr, badOffset) > 3) badOffset += 1;
                // Renormalize to an offset that is a multiple of skipDist.
                baseOff = (((badOffset + (skipDist - 1))/skipDist) * skipDist);
                if (baseOff > endingOffset) break;
                // Try to load entire word.
                badOffset = generateMatches4to2Fast(gbasePtr, baseOff, wordLen, &hashCode);
                continue;
            }
            // The previous read was good, so process it.
            // Max out before we overflow.
            UINT countLimit = sOffsPtr[hashCode+1] - sOffsPtr[hashCode];
            if (counts[hashCode] < countLimit)
            {
                refPtr[sOffsPtr[hashCode]+counts[hashCode]] = baseOff;
                counts[hashCode] += 1;
            }
#ifdef FULLSTATS
            // Output the kmer and count at this reference offset.
            // printChars(refFile, BS->name, 0, BS->nameLen);
            fprintf(refFile, "%u\t%u\n", baseOff, countLimit);
#endif
            // and just get the next skipDist worth
            baseOff += skipDist;
            if (baseOff > endingOffset) break;
            badOffset = generateMatches4to2Fast(gbasePtr, baseOff + hashSaveLength, skipDist, &partialHashCode);
            // and tack it on to the rest from the previous get.
            hashCode = ((hashCode << (skipDist*2)) | partialHashCode) & hashMask;
        }
    }

#ifdef FULLSTATS
    closeForPrint(refFile);
    free(refFileName);
    free(decodeBuf);
#endif

 #ifdef ALLDEBUG
    // Make sure the counts came out right.
    for (int i=0; i<HTsize; i++)
    {
        UINT countLimit = sOffsPtr[i+1] - sOffsPtr[i];
        if (counts[i] != countLimit) fprintf(stderr, "Mismatched before and after counts for %u.  Before=%d, After=%d.\n",
                                             i, countLimit, counts[i]);
    }
    // Make sure there aren't any 0 elements in the ROA (we should get at most one, and there will be none in hg18, for example)
    for (int i=0; i<totalMatches; i++)
    {
        if (refPtr[i] == 0) fprintf(stderr, "ROA entry is zero at location %u.\n", i);
    }
#endif
    // Now we can dispose of the counts.
    free (counts);

////////////
//     Third Pass -- Only if we are randomly sampling high hit k-mers.
///////////
    // If we are randomly sampling the hits to match maxHits, we need to go back through the k-mers and sample the hits.
#ifdef RANDOM_SAMPLE_MAXHITS
    // Initialize the random number generator.
    randState_t * randGen = makeRandState();
    initRandStateDefault(randGen);
    fprintf(stderr, "Randomly Sampling hits for %d-mers that occur more than %d times in the reference.\n", AAs->wordLen, AAs->maxHits);

    // Allocate a temporary array to hold the samples.
    ROFF * samples = (ROFF *)malloc (AAs->maxHits * sizeof(ROFF));
    ROFF newTotalMatches = 0;
    ROFF kmersOverLimit = 0;
    int maxHits = AAs->maxHits;
    for (int i=0; i<HTsize; i++)
    {
        UINT countLimit = sOffsPtr[i+1] - sOffsPtr[i];
        ROFF * srcPtr = &refPtr[sOffsPtr[i]];
        if (countLimit > maxHits)
        {
            kmersOverLimit += 1;
            // This will choose maxHits samples without replacement from the ROA section for this k-mer.
            // Also, the returned sample will retain the sort order of the original (i.e. ascending RO order).
            getRandSample(randGen, srcPtr, countLimit, samples, maxHits);
            srcPtr = samples;
            countLimit = maxHits;
        }
        // Now copy the (possibly reduced) hits to their (possibly) new location.
        ROFF * dstPtr = &refPtr[newTotalMatches];
        for (int j=0; j<countLimit; j++)
            dstPtr[j] = srcPtr[j];
        sOffsPtr[i] = newTotalMatches;
        newTotalMatches += countLimit;
    }
    // Free Structures
    disposeRandState(randGen);
    free(samples);

    // Now we need to resize the output file before writing.
    xfileTotalSize -= (sizeof(ROFF) * (totalMatches - newTotalMatches));
    sOffsPtr[HTsize] = newTotalMatches;
    ints[3] = newTotalMatches;
    // Some statistics
    fprintf(stderr, "%u %d-mers had more than %d hits.\n", kmersOverLimit, AAs->wordLen, maxHits);
    // fprintf(stderr, "Old ROA had %u entries.  New one has %u entries.\n", totalMatches, newTotalMatches);
    // fprintf(stderr, "Old xfileSize %zd, new xfileSize %zd\n", xfileAllocSize, xfileTotalSize);
    totalMatches = newTotalMatches;
#endif // RANDOM_SAMPLE_MAXHITS

#ifdef NOTNOW
    validateHT(HTsize, totalMatches, gbasePtr, refPtr, sOffsPtr, AAs->wordLen, AAs->maxHits);
#endif

    if (AAs->verbose) fprintf(stderr, "%u total %d-mer matches were indexed, resulting in a hash table size of %zd.\n",
                              totalMatches, AAs->wordLen, xfileTotalSize);

/////
//    Output the Index File.
/////

    // Now we are ready to clean up.
    // Close the genome file.
    closeForRead(gfp, gfileSize, gfilePtr);
    blockWrite(xfp, xfilePtr, xfileTotalSize);
    blockFree(xfilePtr, xfileAllocSize);
    closeForWrite(xfp);
    disposeBaseSequences(BSs);
}

void outputCountStatistics(UINT HTsize, UINT * counts, AlignmentArgs_t * AAs)
{
    int wordLen = AAs->wordLen;

    UINT totalMatches = 0;
    UINT zeroCount = 0;
    UINT overFlowCount = 0;
    UINT overFlowTotal = 0;
    int histBucketCount = 10000000;
    int histBucketSize = 1;
    UINT *buckets = (UINT *)malloc(histBucketCount*sizeof(UINT));
    memset(buckets, 0, histBucketCount*sizeof(UINT));
    for (int i=0; i<HTsize; i++)
    {
        totalMatches += counts[i];
        if (counts[i] == 0) zeroCount += 1;
        int bucket = counts[i]/histBucketSize;
        if (bucket >= histBucketCount)
        {
            overFlowCount += 1;
            overFlowTotal += counts[i];
        }
        else buckets[bucket] += 1;
    }

    fprintf(stderr, "Found %u total hits in %u bases.\n", totalMatches, AAs->BSs->totalBaseCount);
    fprintf(stderr, "Found %u %d-mers with zero hits.\n", zeroCount, AAs->wordLen);
    fprintf(stderr, "Found %u %d-mers with greater than %d hits, totaling %u.\n",
            overFlowCount, AAs->wordLen, histBucketCount * histBucketSize, overFlowTotal);

#define quartCount 8
    float percs[quartCount];
    int  kmerIndexes[quartCount];
    UINT kmerThresholds[quartCount];
    int  hitIndexes[quartCount];
    UINT hitThresholds[quartCount];
    percs[0] = 0.25;
    percs[1] = 0.50;
    percs[2] = 0.75;
    percs[3] = 0.90;
    percs[4] = 0.95;
    percs[5] = 0.99;
    percs[6] = 0.999;
    percs[7] = 0.9999;
    for (int i=0; i<quartCount; i++)
    {
        kmerIndexes[i] = -1;
        hitIndexes[i] = -1;
        kmerThresholds[i] = HTsize * percs[i] + 0.5;
        hitThresholds[i] = totalMatches * percs[i] + 0.5;
    }
    UINT newTotalkmers = 0;
    UINT newTotalMatches = 0;
    for (int i=0; i<histBucketCount; i++)
    {
        newTotalkmers += buckets[i];
        newTotalMatches += (buckets[i] * i);
        for (int j=0; j<quartCount; j++)
        {
            if (kmerIndexes[j] == -1 && newTotalkmers >= kmerThresholds[j]) kmerIndexes[j] = i;
            if (hitIndexes[j] == -1 && newTotalMatches >= hitThresholds[j]) hitIndexes[j] = i;
        }
        // fprintf(stderr, "%6d-%6d:%10u    %10u\n", i*histBucketSize, (i+1)*histBucketSize-1, buckets[i], newTotalMatches);
    }
    for (int i=0; i<quartCount; i++)
    {
        fprintf(stderr, "The %f percentile of %u %d-mers is %d, and of %u total matches is %d.\n",
                percs[i], HTsize, wordLen, kmerIndexes[i]*histBucketSize, totalMatches, hitIndexes[i]);
    }
    free(buckets);
}

void validateHT(UINT HTsize, UINT totalMatches, char *gbasePtr, ROFF *refPtr, UINT *sOffsPtr, int wordLen, int maxHits)
{
    // Output some values to see what is going on.
    char * decodeBuf = (char *)malloc(wordLen+1);
    assert(decodeBuf);
    ROFF newTotalMatches = 0;
    for (int i=0; i<HTsize; i++)
    {
        decode(i, decodeBuf, wordLen);
        // fprintf(stderr, "Hash code=%d, =>base seq\"%s\", matches %6d spots with starting offs=%u.\n", i, decodeBuf, counts[i], sOffsPtr[i]);
        int countLimit = sOffsPtr[i+1] - sOffsPtr[i];
        newTotalMatches += countLimit;
        if (countLimit < 0 || countLimit > maxHits)
            fprintf(stderr, "Found count (%d) not matching countlimit (%d) in hash %u => %s\n", countLimit, maxHits, i, decodeBuf);
        for (int j=0; j<countLimit; j++)
        {
            // fprintf(stderr, "@R[%10d]=", refPtr[sOffsPtr[i]+j]);
             check4Seq(stderr, decodeBuf, gbasePtr, refPtr[sOffsPtr[i]+j], wordLen);
            // fprintf(stderr, "\n");
        }
    }
    if (newTotalMatches != totalMatches)
        fprintf(stderr, "Unequal before (%u) and after (%u) totalMatches.\n", totalMatches, newTotalMatches);
    free(decodeBuf);
}
