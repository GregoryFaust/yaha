/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    Query.c         Process a Query file and align reads against the reference.
                             Do so in a mutli-threaded fashion with queries assigned to threads.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdlib.h>
#include <pthread.h>
#include <errno.h>
#include <sys/sysinfo.h>
#include "Math.h"
#include "Math.inl"
// Only output timing information in the experimental version.
#ifndef COMPILE_USER_MODE
#define TIMING
#endif
#include "Timing.inl"

////////////////////////////////////////////
// Perform query alignments.
////////////////////////////////////////////

// This is the one thing we can have global..
// The lock on the input file is used as the thread synchronization method for assigning queries in threads.
// Along with the similar lock on the output file, these are the only two thread syncronizations needed.

FILE * qFile;            // The query file.

//////////////
// Start with some small reading helper functions.
//////////////

// We used unlocked IO while processing a query, with lock and unlock surrounding the query read.
static inline int fgetchar(FILE * in)
{
    return getc_unlocked(in);
}

// Read from current location to the next occurence of specified char.
void readToChar (FILE * in, char fchar, BOOL preceedingNewline)
{
    char prevChar = 0;
    while (TRUE)
    {
        int curChar = fgetchar(in);
        if ((curChar == fchar && (!preceedingNewline || prevChar == '\n')) || curChar == EOF) return;
        prevChar = curChar;
    }
}

void openQueryFile(AlignmentArgs_t * AAs)
{
    // Open the query file.
    qFile = openForSeqRead(AAs->qfileName);

    // Look at the first character to decide if it is a fasta or fastq file.
    // This also positions the input where readNextQuery expects to be.
    flockfile(qFile);
    AAs->fastq = FALSE;
    if (fgetchar(qFile) == '@') AAs->fastq = TRUE;
    funlockfile(qFile);
}

void closeQueryFile(AlignmentArgs_t * AAs)
{
    closeForSeqRead(qFile);
}

void calcNewMaxQueryLength(QueryState_t * QS)
{
    // Leave a 20% buffer in case some later query in the file is larger than this one.
    int maxQueryLength = (int)(QS->queryLen * 1.2);
    // Make sure we don't go over the absolute limit for query length.
    if (maxQueryLength > QS->AAs->maxQueryLength)
        maxQueryLength = QS->AAs->maxQueryLength;
    // Set the new value back in the QS object.
    QS->maxQueryLength = maxQueryLength;
}

void reallocNewMaxQueryLength(QueryState_t * QS)
{
    // Deallocate the current structures.
    finalizeQueries(QS);
    // Calc the new size.
    calcNewMaxQueryLength(QS);
    // Reinitialize structures.
    initializeQueries(QS);
}

void readNextQuery(AlignmentArgs_t * AAs, QueryState_t * QS)
{
    // Lock the input file.
    flockfile(qFile);

    // Continue until we find a valid read or hit end of file.
    QS->queryLen = 0;
    while (TRUE)
    {
        // First read in the ">" or "@" line for the query ID.
        int charCount = 0;
        while (TRUE)
        {
            int curChar = fgetchar(qFile);
            // Check if we are done.
            if (curChar == '\n' || curChar == EOF) break;
            // Input the next char into the query ID buf (if it will fit),
            // Also convert spaces to underscores to aid further pipeline commands.
            if (charCount < MAX_QUERY_ID_LEN)
            {
                QS->queryID[charCount] = (curChar == ' ' ? '_' : curChar);
            }
            charCount += 1;
        }
        if (charCount > MAX_QUERY_ID_LEN)
        {
            fprintf(stderr, "Warning, Query Id length of %d exceeds maximum length %d.  Id will be truncated.\n",
                    charCount, MAX_QUERY_ID_LEN);
            QS->queryIDLen = MAX_QUERY_ID_LEN;
        }
        else
        {
            QS->queryIDLen = charCount;
        }

        // Now read in the query string.
        int reverseOff = AAs->maxQueryLength;
        char breakChar = AAs->fastq ? '+' : '>';
        int qbaseCount = 0;
        BOOL Fail = FALSE;
        while (TRUE)
        {
            int curChar = fgetchar(qFile);
            if (curChar == breakChar || curChar == EOF) break;
            if (curChar == '\n') continue;
            // Check to see that we are within bounds of maximum query length.
            if (qbaseCount >= AAs->maxQueryLength)
            {
                fprintf(stderr, "Warning.  Query sequence exceeds maximum length of %d.  Query will be skipped.\n",
                        AAs->maxQueryLength);
                // Throw away rest of query.
                readToChar(qFile, breakChar, FALSE);
                Fail = TRUE;
                break;
            }
            // process the next char of the query.
            // We fill both the forward and reverse query buffers as we go along.
            // In each case, store both the char for the base for output.
            // And also the 4-bit code for comparison to compressed reference.
            char forCode = map8to4(curChar);
            QS->forwardBuf[qbaseCount] = curChar;
            QS->forwardCodeBuf[qbaseCount] = forCode;
            reverseOff -= 1;
            char revCode = complement4to4(forCode);
            QS->reverseBufBase[reverseOff] = unmap4to8(revCode);
            QS->reverseCodeBufBase[reverseOff] = revCode;
            qbaseCount += 1;
        }
        // For a fastq file, we need to also find the quality string.
        int qualCount = 0;
        if (AAs->fastq)
        {
            // The '+' line can contain extra information.
            // First dump that.
            readToChar(qFile, '\n', FALSE);
            // Reading in the quality informtion is harder, as the @ needs to be preceded by a newline.
            // Given that '@' and newlines can appear in the quality string,
            // is even that good enough to ensure we are at the start of a new query?
            char prevChar = 0;
            while (TRUE)
            {
                int curChar = fgetchar(qFile);
                if ((curChar == '@' && prevChar == '\n') || curChar == EOF) break;
                prevChar = curChar;
                if (curChar == '\n') continue;
                if (qualCount >= AAs->maxQueryLength)
                {
                    fprintf(stderr, "Warning.  Quality score sequence exceeds maximum length of %d.  Query will be skipped.\n",
                            AAs->maxQueryLength);
                    // Throw away rest of quality string
                    readToChar(qFile, '@', TRUE);
                    Fail = TRUE;
                    break;
                }
                QS->qualBuf[qualCount] = curChar;
                qualCount += 1;
            }
            // Make sure the fastq file is well formed.
            if (qbaseCount != qualCount)
            {
                fprintf(stderr, "Warning.  Query sequence (%d) and quality score sequence (%d) have different lengths in fastq file."
                        "  Query will be skipped.\n", qbaseCount, qualCount);
                Fail = TRUE;
            }
        }
        if (qbaseCount > 0 && qbaseCount < AAs->wordLen)
        {
            fprintf(stderr, "Query length must be at least wordlen bases long. Query will be skipped.\n");
            Fail = TRUE;
        }
        // If anything has gone wrong, try again.
        if (Fail) continue;
        funlockfile(qFile);
        // Fill in the QS fields of interest.
        // Since reverse buf is filled in from back to front, now set the starting ptrs.
        QS->reverseCodeBuf = &QS->reverseCodeBufBase[reverseOff];
        QS->reverseBuf = &QS->reverseBufBase[reverseOff];
        QS->queryLen = qbaseCount;
        // Use the query sequence to generate a random seed for stocastically choosing amonst equal alignments.
        // By using the query sequence for the seed, the same one is chosen among equals regardless of file context of the query.
        generateRandomSeed(QS);
        return;
    }
    // This will never be executed, as EOF will result in qbaseCount of zero in while(TRUE) loop.
    // It is here just to avoid the compiler warning.
    return;
}

// We need 2-bit copression to calculate the hash key.
// Since we already have the query in unpacked 4 codes, this can be very fast.
// We only need to check for non-ACGT codes, and compress to 2 bits if OK.
static inline int generateMatches4UNPto2Fast(char * seq, ROFF offset, int length, UINT* hash)
{
    UINT retval = 0;
    for (int i=0; i<length; i++)
    {
        UBYTE code = seq[offset + i];
        if (code > 3) return offset + i + 1;
        retval = (retval << 2) + code;
    }
    *hash = retval;
    return 0;
}

/////////////////
//     Start of the main workhorse routine.
//     There is one call to this routine per thread.
//     We pass in the pre-allocated initialized structure for the thread to avoid malloc contention on startup.
/////////////////
//     A lot of the apparent complexity of this function is due to optionally keeping stats and timings.
//     Both of which have been tremendously useful to tune the algorithms/heuristics especially for speed.
/////////////////

void * processQueries(void * QSptr)
{
    setupTimers();
    defineTimer(setupSecs);
    defineTimer(findDiagonals);
    defineTimer(findClumps);
    defineTimer(filterClumps);
    defineTimer(processClumps);
    defineTimer(printClumpss);
    defineTimer(totalTimer);

    // Use skipDist of 1 for query, regardless of what was used to create the hash table.
    int skipDist = 1;
    // Cache the wordLen on the stack to speed loop access.
    QueryState_t * QS = (QueryState_t *) QSptr;
    AlignmentArgs_t * AAs = QS->AAs;
    int wordLen = AAs->wordLen;
    int maxHits = AAs->maxHits;
    ROFF *startingOffs = AAs->startingOffs;

#ifdef STATS
    // Let's keep some statistics about the queries we process.
    int queryCount = 0;
    int totalQueryLength = 0;
    uint64_t totalTotalCount = 0;
    int minQLength = 1000000000;
    int maxQLength = 0;
    int minTotalCount = 1000000000;
    int maxTotalCount = 0;

    uint64_t totalClumpsOut = 0;
    int minClumpsPerQuery = 1000000000;
    int maxClumpsPerQuery = -1;
    int totalNonAlignedQueries = 0;
#endif

/////
//     Now we are ready to go into the main query read loop.
/////

    // Calculate some masks for use later.
    UINT hashMask = ((UINT)0xFFFFFFFF) >> (32 - 2*wordLen);
    int  hashSaveLength = wordLen - skipDist;

    startTime();

    // To get started, make sure we have read a first query.
    // This is a partial unwind of the first time through the loop.
    // It allows the main thread to read the first query to get a query length estimate to use to create the other threads.
    if (QS->queryLen == 0) readNextQuery(AAs, QS);

    while (QS->queryLen > 0)
    {
///////
//     First Phase.  Read in the query, reverse it, and do the initial hashing.
///////

        // See if we need to reallocate the query length dependent structures.
        if (QS->queryLen > QS->maxQueryLength) reallocNewMaxQueryLength(QS);

        QOFF qbaseCount = QS->queryLen;

#ifdef STATS
        fprintQueryId(QS, stderr);
#endif
#ifdef QUERYSTATS
        struct timeval starttime, endtime;
        if (AAs->queryStats)
        {
            QS->alignCount = 0;
            QS->alignOutCount = 0;
            QS->DPCount = 0;
            QS->usec = 0;
            QS->seedMatches = 0;
            QS->diagRegionCount = 0;
            QS->diagRegionTotal = 0;
            gettimeofday(&starttime, NULL);
        }
#endif

/////
//     Now match the query against both the reference DNA and its reverse-complement strand.
//     Therefore the bases in the reversed query are also complemented.
/////

        // Do the forward query first.
        int offsetCount = ((qbaseCount - wordLen)/skipDist) + 1;
        for (int rev=0; rev<=1; rev++)
        {
            QS->reversed = (BOOL)rev;

            if (QS->reversed)
            {
                QS->queryCodeBuf = QS->reverseCodeBuf;
                QS->queryBuf = QS->reverseBuf;
            }
            else
            {
                QS->queryCodeBuf = QS->forwardCodeBuf;
                QS->queryBuf = QS->forwardBuf;
            }

/////
//     Find the hash seeds, filter overhits, and fill in structures from the hash table.
/////

#ifdef STATS
            UINT maxHitsOverCount = 0;
#endif
            // The loop complication is because we want to both encode each base once, and skip over non-ACGT codes.
            UINT totalCount = 0;
            int baseOff = 0;
            int endingOffset = qbaseCount - wordLen;
            UINT hashCode = 0;
            UINT partialHashCode = 0;
            int badOffset = generateMatches4UNPto2Fast(QS->queryCodeBuf, baseOff, wordLen, &hashCode);
            while (TRUE)
            {
                if (badOffset != 0)
                {
                    // Get rid of all the bad ones in a row.
                    while (badOffset <= endingOffset && QS->queryCodeBuf[badOffset] > 3) badOffset += 1;
                    // Make sure the ones we skipped over have their count set to zero.
                    // We need to handle the case where the above has gone too far.
                    badOffset = MIN(badOffset, endingOffset + 1);
                    // fprintf(stderr, "After throwing away bases, badOffset %d baseOff %d endingOffset %d \n", badOffset, baseOff, endingOffset);
                    for (int i=baseOff; i<badOffset; i++) QS->offsetCounts[i].count = 0;
                    // Move forward on the query.
                    baseOff = badOffset;
                    if (baseOff > endingOffset) break;
                    // Try to load entire word.
                    badOffset = generateMatches4UNPto2Fast(QS->queryCodeBuf, baseOff, wordLen, &hashCode);
                    continue;
                }
                // The previous read was good, so process it.
                // Max out before we overflow.
                UINT count = startingOffs[hashCode+1] - startingOffs[hashCode];
                if (count <= maxHits)
                {
                    totalCount += count;
                    QS->offsetCounts[baseOff].count = count;
                    QS->offsetCounts[baseOff].sOffset = startingOffs[hashCode];
                    QS->offsetCounts[baseOff].newCount = 0;
                }
                else
                {
#ifdef STATS
                    maxHitsOverCount += 1;
#endif
                    QS->offsetCounts[baseOff].count = 0;
                }
                // and just get the next skipDist worth
                baseOff += 1;
                if (baseOff > endingOffset) break;
                badOffset = generateMatches4UNPto2Fast(QS->queryCodeBuf, baseOff + hashSaveLength, 1, &partialHashCode);
                // and tack it on to the rest from the previous get.
                hashCode = ((hashCode << 2) | partialHashCode) & hashMask;
            }

#ifdef STATS
            fprintf(stderr, "K-mers above maxHits %u\n", maxHitsOverCount);
            totalTotalCount += totalCount;
            if (totalCount > 0 && totalCount < minTotalCount) minTotalCount = totalCount;
            if (totalCount > maxTotalCount) maxTotalCount = totalCount;
#endif
#ifdef QUERYSTATS
            if (AAs->queryStats)QS->seedMatches += totalCount;
#endif
            addToTimer(setupSecs);

            // If the totalCount is zero, we have no hits to process.
            if (totalCount == 0) continue;

/////
//     Second Phase.  Go through all hash matches, and look for diagonal fragments.
/////

            int fragCount = findFragmentsSort(AAs, QS, offsetCount);
            addToTimer(findDiagonals);

/////
//     Third Phase.  Put Fragments into potential alignments (clumps).
/////

            processFragmentsGapped(AAs, QS, fragCount);
            addToTimer(findClumps);

        }

        // Align and Score clumps, filtering those that don't meet thresholds.
        postProcessClumps(QS, AAs);
        addToTimer(processClumps);

        // Run OQC algorithm if requested.
        // If not, at least get rid of duplicates.
        if   (AAs->OQC) postFilterBySimilarity(AAs, QS);
        else            postFilterRemoveDups(AAs, QS);
        addToTimer(filterClumps);

        // The file lock and unlock here act as the thread synchromization mechanism.
        // This will ensure that all alignments for a given query are in contiguous lines in the output file.
        // Otherwise, many downstream tools such as samtools will be unhappy.
        flockfile(AAs->outFile);

        // This is to avoid compiler warning when STATS not on.
#ifdef STATS
        int tempClumpsPerQuery = printClumps(QS, AAs);
#else
        printClumps(QS, AAs);
#endif

        funlockfile(AAs->outFile);
        addToTimer(printClumpss);

#ifdef STATS
        queryCount += 1;
        totalQueryLength += qbaseCount;
        if (qbaseCount < minQLength) minQLength = qbaseCount;
        if (qbaseCount > maxQLength) maxQLength = qbaseCount;
        totalClumpsOut += tempClumpsPerQuery;
        if (tempClumpsPerQuery > maxClumpsPerQuery) maxClumpsPerQuery = tempClumpsPerQuery;
        if (tempClumpsPerQuery > 0 && tempClumpsPerQuery < minClumpsPerQuery) minClumpsPerQuery = tempClumpsPerQuery;
        if (tempClumpsPerQuery == 0) totalNonAlignedQueries += 1;
#endif

#ifdef QUERYSTATS
        if (AAs->queryStats)
        {
            gettimeofday(&endtime, NULL);
            QS->usec = ((endtime.tv_sec - starttime.tv_sec) * 1000000 +
                        (endtime.tv_usec - starttime.tv_usec));
            printChars(AAs->qsFile, QS->queryID, 0, QS->queryIDLen);
            double avgDiagRegionSize = (QS->diagRegionCount > 0) ? (((double)QS->diagRegionTotal/(double)QS->diagRegionCount)+0.0) : 0;
            fprintf(AAs->qsFile, "\t%u\t%u\t%f\t%u\t%u\t%u\t%u\n", QS->queryLen, QS->seedMatches, avgDiagRegionSize,
                    QS->alignCount, QS->DPCount, QS->alignOutCount, QS->usec);
        }
#endif

        // Clear query state to get ready for next query.
        resetQueryState(QS);
        // Read the next query to setup for next iteration.
        readNextQuery(AAs, QS);
    }

//////////////////////
//    Done with processing the query file.
//    Close files and clean up storage.
/////////////////////

#ifdef TIMING
    addToTimer(setupSecs);
    endTime(totalTimer);

    flockfile(stderr);
    // if (AAs->threadCount > 1) fprintf(stderr, "Thread %d statistics.\n", QS->threadNum);
    fprintTimerWithTotalPercent(stderr, "Query setup and hashing took: ", setupSecs, totalTimer);
    fprintTimerWithTotalPercent(stderr, "Sorting k-mers and forming fragments took: ", findDiagonals, totalTimer);
    fprintTimerWithTotalPercent(stderr, "Collecting fragments into clumps (using graph) took: ", findClumps, totalTimer);
    fprintTimerWithTotalPercent(stderr, "Aligning and scoring clumps took: ", processClumps, totalTimer);
    fprintTimerWithTotalPercent(stderr, "Filtering clumps took: ", filterClumps, totalTimer);
    fprintTimerWithTotalPercent(stderr, "Printing clumps took: ", printClumpss, totalTimer);
    funlockfile(stderr);
#endif

#ifdef STATS
    if (AAs->verbose)
    {
        fprintf(stderr, "%d queries processed.\n", queryCount);
        fprintf(stderr, "Query Lengths vary from %d to %d with average %d.\n",
                minQLength, maxQLength, (totalQueryLength/queryCount));
        fprintf(stderr, "Total Counts vary from %d to %d with average %"PRIu64".\n",
                minTotalCount, maxTotalCount, (totalTotalCount/(2 * queryCount)));
        fprintf(stderr, "There were %d queries with no Alignment.\n", totalNonAlignedQueries);
        if (totalClumpsOut <= 0) fprintf(stderr, "No Alignments found.\n");
        else
        {
            fprintf(stderr, "Total Alignments Output = %"PRIu64", average %4.2f per non-zero query.\n", totalClumpsOut,
                    ((double)totalClumpsOut/(queryCount-totalNonAlignedQueries)));
            fprintf(stderr, "Of those queries with an alignment, the min number of alignments was %d.\n", minClumpsPerQuery);
            fprintf(stderr, "The max number of alignments per query was %d.\n", maxClumpsPerQuery);
        }
    }
#endif

    // Deallocate thread local structures.
    finalizeQueries(QS);
    disposeQueryState(QS);
    return NULL;
}


/////////////////
//     Start of the routine called once to process the entire query file.
//     It then spawns threads as requested by user parameter.
/////////////////

void processQueryFile(AlignmentArgs_t * AAs)
{

/////
//     Open files, and initialize structures
//     The compressed genome and index structures are read only and shared amonst threads via AAs.
/////

    // Open the query file.
    openQueryFile(AAs);

    // Open the genome file.
    char *gfilePtr;
    FSIZE gfileSize;
    FDES gfp = openForRead(AAs->gfileName, &gfileSize, &gfilePtr, TRUE);

    // Read in the .nib2 file header.
    BaseSequences_t * BSs = loadBaseSequences(gfilePtr);
    normalizeBaseSequences(BSs);
    if (AAs->verbose) fprintf(stderr, "Read in %d reference sequences from %s.\n", BSs->curCount, AAs->gfileName);
    AAs->BSs = BSs;
    AAs->basePtr = BSs->basePtr;
    AAs->maxROff = baseSequencesMaxROff(BSs);

    // Open the output alignment file.
    // We have no idea how much output we will generate.
    // It will be all written sequentially.
    // So, no fancy memory mapped IO for this output.
    // Just do formatted ouput.
    AAs->outFile = openForPrint(AAs->ofileName);

#ifdef QUERYSTATS
    if (AAs->queryStats)
    {
        AAs->qsFile = openForPrint(AAs->qsfileName);
    }
#endif

    // Some of the file formats have headers.
    outputFileHeader(AAs);

    // Open the index file.
    // This is harder as it requires putting all pointers to various pieces back together.
    char *xfilePtr;
    FSIZE xfileSize;
    FDES xfp = openForRead(AAs->xfileName, &xfileSize, &xfilePtr, TRUE);

    // Start by getting the wordLen and totalCount from the file Header.
    UINT *xuintPtr = (UINT *)xfilePtr;
    int version = xuintPtr[0];
    if (version != CURRENT_INDEX_FILE_VERSION)
        fatalError("Index file version is out of date.\nPlease remake index file and try again.");
    AAs->wordLen = xuintPtr[1];
    int fileMaxHits = xuintPtr[2];
    if (fileMaxHits < AAs->maxHits)
    {
        fprintf(stderr, "WARNING: Index file made with maxHits of %d, while %d specified for this query run.\n"
                "Mimimum of two (%d) will be used.\n", fileMaxHits, AAs->maxHits, fileMaxHits);
        AAs->maxHits = fileMaxHits;
    }
    // Now we can recalculate the hash table size.
    UINT HTsize = iexp(4, AAs->wordLen);

    // Now we will walk through the various structures, assigning pointers of right type to right spots.
    FSIZE sizeofUINT = sizeof(UINT);
    FSIZE inputSize;
    // Space for four numbers in header.
    inputSize = 4*sizeofUINT;
    FSIZE startingOffBase = inputSize;
    // Space for reference starting Offsets array.
    inputSize += sizeofUINT * (HTsize + 1);
    FSIZE referenceBase = inputSize;

    // We've calculated all the offsets, now get the pointers.
    AAs->startingOffs = (ROFF *)(xfilePtr + startingOffBase);
    AAs->ROAPtr = (ROFF *)(xfilePtr + referenceBase);

#ifdef DEBUG
    fprintf(stderr, "Number of Hash Tables entries is %u, and ROA has %u entries.\n",
            HTsize, (UINT)((xfileSize-inputSize)/sizeof(ROFF)));
    fprintf(stderr, "Hash Table sizes are as follows: offs=%zd, refs=%zd.\n", sizeof(UINT)*(HTsize + 1), xfileSize-inputSize);
    fprintf(stderr, "The index pointers are: xfilePtr=%p, sOffsPtr=%p,  and refPtr=%p\n", xfilePtr, startingOffs, refs);
#endif

    // Read in first query to get an estimate of the size of the queries in the file.
    QueryState_t * QS = makeQueryState(AAs);
    readNextQuery(AAs, QS);
    calcNewMaxQueryLength(QS);
    QOFF maxQueryLength = QS->maxQueryLength;
    initializeQueries(QS);

    // Create extra threads.
    // We will use the main thread as well, so only need to spawn one less than requested.
    int numthreads = AAs->numThreads - 1;
    pthread_t * threads = NULL;
    if (numthreads > 0)
    {
        // Keep the thread count to something the machine can handle.
        int nprocs = get_nprocs();
        if (numthreads + 1 > nprocs)
        {
            fprintf(stderr, "Warning.  Requested number of threads (%d) is greater than number of processors."
                    "  %d threads will be used.\n", AAs->numThreads, nprocs);
            numthreads = nprocs - 1;
        }

        // Do the storage allocation for threads.
        // We allocate the thread local storage here as well to avoid malloc contention while spawning threads.
        threads = (pthread_t *)malloc(numthreads * sizeof(pthread_t));
        for (int i=0; i<numthreads; i++)
        {
            QueryState_t * QS = makeQueryState(AAs);
            // Set in an initial value for maxQueryLength based on first query read.
            QS->maxQueryLength = maxQueryLength;
            initializeQueries(QS);
            pthread_t * thread = threads+i;
            int err = pthread_create(thread, NULL, &processQueries, QS);

            // If we fail to create the thread for any reason, bail from trying to make more.
            if (err != 0)
            {
                fprintf(stderr, "Warning.  Thread creation resulted in error %d.  %d threads will be used.\n", err, i+1);
                numthreads = i;
                // Deallocate the unused thread local structures.
                finalizeQueries(QS);
                disposeQueryState(QS);
                break;
            }
        }
    }

    // Now also do work in this thread.
    // This won't return until the query file is exhausted.
    processQueries(QS);

    // Wait for all the rest of the threads (if any) to finish.
    for (int i=0; i<numthreads; i++)
    {
        pthread_join(threads[i], NULL);
    }
    if (threads != NULL) free(threads);

    // Free resources.
    // First close files.
    closeQueryFile(AAs);
    closeForRead(gfp, gfileSize, gfilePtr);
    closeForRead(xfp, xfileSize, xfilePtr);
    closeForPrint(AAs->outFile);

#ifdef QUERYSTATS
    if (AAs->queryStats)
    {
        closeForPrint(AAs->qsFile);
    }
#endif

    // Free allocated structures.
    disposeBaseSequences(BSs);
}
