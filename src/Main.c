/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    Main.c          Parse the command line arguments, and call the appropriate function to do the work.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "Math.h"
#define TIMING
#include "Timing.inl"

void printUsageString()
{
    char* useString = "Usage (Default parameter values shown in parenthesis):\n\n"
        "To create an index:\n"
        "yaha -g genomeFilename [-H maxHits (65525)] [-L wordLen (15)] [-S Skip-distance (1)]\n"
        "The genome file can be a FASTA file, or a nib2 file (created by a previous yaha index operation).\n\n"
        "To align queries:\n"
        "yaha -x yahaIndexFile [-q queryFile|(stdin)] [-o8|(-osh)|-oss outputFile|(stdout)][AdditionalOptions]\n"
        "The query file can be either a FASTA file or a FASTQ file.\n"
        "-o8     produces alignment output in modified Blast8 format.\n"
        "-osh    produces alignment output in SAM format with hard clipping.\n"
        "-oss    produces alignment output in SAM format with soft clipping.\n"
        "[-t numThreads (1)]\n\n"
        "Additional General Alignment Options:\n"
        "    [-BW BandWidth (5)] [-G maxGap (50)] [-H maxHits (650)] [-M minMatch (25)]\n"
        "    [-MD MaxDesert (50)] [-P minPercent-identity (0.9)] [-X Xdropoff (25)]\n\n"
        "[-AGS (Y)|N] controls use of Affine Gap Scoring.\n"
        "If -AGS is off, a simple edit distance calculation is done.  If on, the following are used:\n"
        "    [-GEC GapExtensionCost (2)] [-GOC GapOpenCost (5)] [-MS MatchScore (1)] [-RC ReplacementCost (3)]\n\n"
        "[-OQC (Y)|N] controls use of the Optimal Query Coverage Algorithm.\n"
        "If -OQC if off, all alignments meeting above criteria are output.\n"
        "If on, a set of \"primary\" alignments are found that optimally cover the length of the query, using the following options:\n"
        "    [-BP BreakpointPenalty (5)] [-MGDP MaxGenomicDistancePenalty (5)] [-MNO MinNonOverlap (minMatch)]\n"
        "The cost of inserting a breakpoint in the Optimal Coverage Set is BP*MIN(MGDP, Log10(genomic distance between reference loci)).\n\n"
        "[-FBS Y|(N)] controls inclusion of \"secondary\" alignments similar to a primary alignment found using OQC.\n"
        "If -FBS is on, the following are used.  A \"secondary\" alignemnt must satisfy BOTH criteria.\n"
        "    [-PRL PercentReciprocalLength (0.9)] [-PSS PercentSimilarScore (0.9)]\n\n"
#ifndef COMPILE_USER_MODE
        "Additional experimental parameters:\n"
        "To compress a FASTA file to a nib2 file without creating an index:\n"
        "yaha -g fastaGenomeFile -c\n"
        "To uncompress a nib2 file back into a FASTA file:\n"
        "yaha -g nib2GenomeFile -u\n"
        "For finer control of alignments:\n"
        "[-I maxIntron (maxGap)] allows separate control of max deletion size vs. maxGap for insertion size.\n"
#endif
        "For a more detailed help message, type yaha -xh.\n";
    fprintf(stderr, "%s", useString);
}

void printExtendedUsageString()
{
    char* useString = "Usage (Default parameter values shown in parenthesis):\n\n"
        "To create an index:\n"
        "yaha -g genomeFilename [-H maxHits (65525)] [-L wordLen (15)] [-S Skip-distance (1)]\n"
        "The genome file can be a FASTA file, or a nib2 file (created by a previous yaha index operation).\n\n"
        "[-H maxHits (655250] K-mers that appear > maxHits times in the genome will be randomly sampled down to max-hits locations.\n"
        "[-L wordLen (15)] the size of k-mers used for seed hits.\n"
        "[-S Skip-distance (1)] the distance to skip ahead in the genome before forming the next k-mer to index.\n"
        "  Skip-distance must be between 1 and wordLen (inclusive).\n\n"
        "To align queries:\n"
        "yaha -x yahaIndexFile [-q queryFile|(stdin)] [-o8|(-osh)|-oss outputFile|(stdout)][AdditionalOptions]\n"
        "The query file can be either a FASTA file or a FASTQ file.\n"
        "-o8     produces alignment output in modified Blast8 format.\n"
        "-osh    produces alignment output in SAM format with hard clipping.\n"
        "-oss    produces alignment output in SAM format with soft clipping.\n"
        "[-t numThreads (1)]\n\n"
        "Additional General Alignment Options:\n"
        "[-BW BandWidth (5)] the size of the bands used on either side of the diagonal during Smith-Waterman.\n"
        "[-G maxGap (50)] the max size of indels within alignments.  Larger indels will cause split alginments.\n"
        "[-H maxHits (650)] K-mers that appear > maxHits times in the index will be ignored.\n"
        "  To take advantage of k-mer sampling, use the same value of maxHits during index creation and alignment.\n"
        "[-M minMatch (25)] the minimum number of non-overlapping seed bases (not seed count) needed to initiate an alignment.\n"
        "[-MD MaxDesert (50)] the max stretch along the query with no index seed hits before an alignment will be broken in two.\n"
        "[-P minPercent-identity (0.9)] alignments with < this ratio of matches/alignment-length are thrown away.\n"
        "[-X Xdropoff (25)] alignment extension is stopped when an extension has score >= X-dropoff below the best extension.\n\n"
        "[-AGS (Y)|N] controls use of Affine Gap Scoring.\n"
        "If -AGS is off, a simple edit distance calculation is done.\n"
        "If on, the following parameters are used.  All are specified as positive integers.\n"
        "[-GEC GapExtensionCost (2)] cost for each base in an indel, including the first one.\n"
        "[-GOC GapOpenCost (5)] cost for starting a new indel.\n"
        "[-MS MatchScore (1)] positive score added for each matching base.\n"
        "[-RC ReplacementCost (3)] cost for each mis-matched base.\n\n"
        "[-OQC (Y)|N] controls use of the Optimal Query Coverage Algorithm.\n"
        "If -OQC if off, all alignments meeting above criteria are output.\n"
        "If on, a set of \"primary\" alignments are found that optimally cover the length of the query, using the following options:\n"
        "The cost of inserting a breakpoint in the Optimal Coverage Set is BP*MIN(MGDP, Log10(genomic distance between reference loci)).\n"
        "[-BP BreakpointPenalty (5)] Genomic distance penalty is multiplied by this factor.\n"
        "[-MGDP MaxGenomicDistancePenalty (5)] Log10 of the genomic distance is capped at this value.\n"
        "  Alignments on two different sequences (chromosomes) are always penalized by the max genomic distance penalty.\n"
        "[-MNO MinNonOverlap (minMatch)] Two alignments must each have >= this number of non-overlapping bases to be included together in the OCS.\n\n"
        "[-FBS Y|(N)] controls inclusion of \"secondary\" alignments similar to a primary alignment found using OQC.\n"
        "If -FBS is on, the following parameters are used.  A \"secondary\" alignemnt must satisfy BOTH criteria.\n"
        "[-PRL PercentReciprocalLength (0.9)] The overlap between the primary and secondary alignment must cover >= this %% of the length of BOTH alignemnts.\n"
        "[-PSS PercentSimilarScore (0.9)] The secondary alignment's score must be >= this percent of the primary alignment's score.\n"
#ifndef COMPILE_USER_MODE
        "\nAdditional experimental parameters:\n"
        "To compress a FASTA file to a nib2 file without creating an index:\n"
        "yaha -g fastaGenomeFile -c\n"
        "To uncompress a nib2 file back into a FASTA file:\n"
        "yaha -g nib2GenomeFile -u\n"
        "For finer control of alignments:\n"
        "[-I maxIntron (maxGap)] allows separate control of max deletion size vs. maxGap for insertion size.\n"
#endif
        "\n";
    fprintf(stderr, "%s", useString);
}

void printUsageStringAbort()
{
    printUsageString();
    exit(1);
}

// Helper functions for parsing parameter values by type.
inline BOOL streq(char * s1, const char * s2) __attribute__((always_inline));
inline BOOL streq(char * s1, const char * s2)
{
    return (strcmp(s1, s2) == 0);
}

BOOL parseBOOL(char * inString, char * key)
{
    if (strlen(inString) == 1)
    {
        if (strchr("YyTt", inString[0]) != NULL) return TRUE;
        if (strchr("NnFf", inString[0]) != NULL) return FALSE;
    }
    fprintf(stderr, "%s is not a valid value for parameter %s.\n", inString, key);
    fprintf(stderr, "Use one of 'YyTt' for Yes and 'NnFf' for No.\n\n.");
    printUsageString();
    exit(1);
}

int parseInt(char * inString, char * key)
{
    int retval = atoi(inString);
    if (retval < 0)
    {
        fprintf(stderr, "%s is not a valid value for parameter %s.\n", inString, key);
        fprintf(stderr, "Value must be a positive integer.\n\n");
        printUsageStringAbort();
    }
    return retval;
}

float parseFloat(char * inString, char * key)
{
    float retval = atof(inString);
    if (retval <= 0.0 || retval > 1.0)
    {
        fprintf(stderr, "%s is not a valid value for parameter %s.\n", inString, key);
        fprintf(stderr, "Value must be in the range 0<value<=1.0.\n\n");
        printUsageStringAbort();
    }
    return retval;
}

char * isFileNameOrStdin(char * inString)
{
    if (streq(inString, "-") || streq(inString, "-stdin") || streq(inString, "stdin"))
        return strdup("stdout");
    return strdup(inString);
}

char * isFileNameOrStdout(char * inString)
{
    if (streq(inString, "-stdout") || streq(inString, "stdout"))
        return strdup("stdout");
    return strdup(inString);
}

int main(int argc, char *argv[])
{
    // Output the version number.
    fprintf(stderr, "YAHA version 0.1.%d\n", BUILDNUM);

    // See if we have any arguments!
    if (argc <= 1)
    {
        printUsageString();
        return 0;
    }

    /*
     * Parse command line arguments.
     *
     * The gfile and the xfile will have the same name with different extensions.
     *
     */

    // Create structure to store arguments, and intialize default values.
    AlignmentArgs_t * AAs = makeAlignmentArgs();

    // Flags that represent which operation will be performed.
    // By default, we will make an index unless we are given a query file.
    BOOL query = FALSE;
    BOOL compress = FALSE;
    BOOL uncompress = FALSE;
    BOOL index = TRUE;

    for (int x = 1; x < argc; x++)
    {
        // First the general mode and I/O parameters.
        if (streq(argv[x], "-h") || streq(argv[x], "-?"))
        {
            printUsageString();
            return 0;
        }
        else if (streq(argv[x], "-xh"))
        {
            printExtendedUsageString();
            return 0;
        }
        else if (streq(argv[x], "-g"))
        {
            x++;
            AAs->gfileName = argv[x];
        }
        else if (streq(argv[x], "-q"))
        {
            x++;
            AAs->qfileName = isFileNameOrStdin(argv[x]);
            query = TRUE;
            index = FALSE;
        }
        else if (streq(argv[x], "-o8"))
        {
            AAs->outputBlast8 = TRUE;
            AAs->outputSAM = FALSE;
            x++;
            AAs->ofileName = isFileNameOrStdout(argv[x]);
        }
        // The H stands for "hard" clipping of query sequence in SAM file.
        // The S stands for "soft" clipping of qyery sequence in SAM file.
        else if (streq(argv[x], "-osh"))
        {
            AAs->outputBlast8 = FALSE;
            AAs->outputSAM = TRUE;
            AAs->hardClip = TRUE;
            x++;
            AAs->ofileName = isFileNameOrStdout(argv[x]);
        }
        else if (streq(argv[x], "-oss"))
        {
            AAs->outputBlast8 = FALSE;
            AAs->outputSAM = TRUE;
            AAs->hardClip = FALSE;
            x++;
            AAs->ofileName = isFileNameOrStdout(argv[x]);
        }
        else if (streq(argv[x], "-t"))
        {
            x++;
            AAs->numThreads = parseInt(argv[x], "-t");
        }
        else if (streq(argv[x], "-v"))
        {
             AAs->verbose = TRUE;
        }
        else if (streq(argv[x], "-x"))
        {
            x++;
            AAs->xfileName = argv[x];
            query = TRUE;
            index = FALSE;
        }
        // Now the parameters not available in USER_MODE.
#ifndef COMPILE_USER_MODE
        else if (streq(argv[x], "-c"))
        {
            compress = TRUE;
            index = FALSE;
        }
        else if (streq(argv[x], "-u"))
        {
            uncompress = TRUE;
            index = FALSE;
        }
#ifdef QUERYSTATS
        else if (streq(argv[x], "-qs"))
        {
            AAs->queryStats = TRUE;
            x++;
            AAs->qsfileName = argv[x];
        }
#endif
#endif

        // Now the parameters that control Algorithms.
        // Start with the ones that are used in both USER_MODE and experimental compiles.
        // Of those, start with the index related parameters.
        else if (streq(argv[x], "-H"))
        {
            x++;
            AAs->maxHits = parseInt(argv[x], "-H");
        }
        else if (streq(argv[x], "-L"))
        {
            x++;
            AAs->wordLen = parseInt(argv[x], "-L");
        }

        else if (streq(argv[x], "-S"))
        {
            x++;
            AAs->skipDist = parseInt(argv[x], "-S");
        }
        // Now the general alignment args.
        else if (streq(argv[x], "-BW"))
        {
            x++;
            AAs->bandWidth = parseInt(argv[x], "-BW");
        }
        else if (streq(argv[x], "-G"))
        {
            x++;
            AAs->maxGap = parseInt(argv[x], "-G");
        }
        else if (streq(argv[x], "-M"))
        {
            x++;
            AAs->minMatch = parseInt(argv[x], "-M");
        }
        else if (streq(argv[x], "-MD"))
        {
            x++;
            AAs->maxDesert = parseInt(argv[x], "-MD");
        }
        else if (streq(argv[x], "-P"))
        {
            x++;
            AAs->minIdentity = parseFloat(argv[x], "-P");
        }
        else if (streq(argv[x], "-X"))
        {
            x++;
            AAs->XCutoff = parseInt(argv[x], "-X");
        }
        // Now the Affine Gap Scoring args.
        else if (streq(argv[x], "-AGS"))
        {
            x++;
            AAs->affineGapScoring = parseBOOL(argv[x], "-AGS");
        }
        else if (streq(argv[x], "-GEC"))
        {
            x++;
            AAs->GECost = parseInt(argv[x], "-GEC");
        }
        else if (streq(argv[x], "-GOC"))
        {
            x++;
            AAs->GOCost = parseInt(argv[x], "-GOC");
        }
        else if (streq(argv[x], "-MS"))
        {
            x++;
            AAs->MScore = parseInt(argv[x], "-MS");
        }
        else if (streq(argv[x], "-RC"))
        {
            x++;
            AAs->RCost = parseInt(argv[x], "-RC");
        }
        // Now the OQC and FBS parameters.
        else if (streq(argv[x], "-OQC"))
        {
            x++;
            AAs->OQC = parseBOOL(argv[x], "-OQC");
        }
        else if (streq(argv[x], "-BP"))
        {
            x++;
            AAs->BPCost = parseInt(argv[x], "-BP");
        }
        else if (streq(argv[x], "-MGDP"))
        {
            x++;
            AAs->maxBPLog = parseInt(argv[x], "-MGDP");
        }
        else if (streq(argv[x], "-MNO"))
        {
            x++;
            AAs->OQCMinNonOverlap = parseInt(argv[x], "-MNO");
        }
        else if (streq(argv[x], "-FBS"))
        {
            x++;
            AAs->FBS = parseBOOL(argv[x], "-FBS");
        }
        else if (streq(argv[x], "-PRL"))
        {
            x++;
            AAs->FBS_PSLength = parseFloat(argv[x], "-PRL");
        }
        else if (streq(argv[x], "-PSS"))
        {
            x++;
            AAs->FBS_PSScore = parseFloat(argv[x], "-PSS");
        }
        // Now the ones only used in experimmental builds.
#ifndef COMPILE_USER_MODE
        else if (streq(argv[x], "-I"))
        {
            x++;
            AAs->maxIntron = parseInt(argv[x], "-I");
        }
#ifdef NOTNOW
        // This parameter is not currently supported.
        else if (streq(argv[x], "-MSA"))
        {
            x++;
            AAs->FBSMaxSimilar = parseInt(argv[x], "-MSA");
        }
#endif
        else if (streq(argv[x], "-R"))
        {
            x++;
            AAs->minRawScore = parseInt(argv[x], "-R");
        }
#endif
        else
        {
            fprintf(stderr, "%s is not a valid option.\n\n", argv[x]);
            printUsageStringAbort();
        }
    }

    //////
    // Check to see if we have a valid set of file name options.
    //////
    char * operString;
    if (index)
        operString = "index creation";
    else if (compress)
        operString = "genome compression";
    else if (uncompress)
        operString = "genome decompression";

    if ((index || compress || uncompress))
    {
        if (AAs->gfileName == NULL)
        {
            fprintf(stderr, "Genome file specification (-g) is required for %s.\n\n", operString);
            printUsageStringAbort();
        }
        if (AAs->ofileName != NULL)
        {
            fprintf(stderr, "Output file specification is not allowed during %s.\n\n", operString);
            printUsageStringAbort();
        }
    }
    if (query)
    {
        if (AAs->gfileName != NULL)
        {
            fprintf(stderr, "Genome file specification (-g) is not allowed for query alignment.\n");
            printUsageStringAbort();
        }
        if (AAs->xfileName == NULL)
        {
            fprintf(stderr, "Index file specification (-x) is required for query alignment.\n");
            printUsageStringAbort();
        }
        if (AAs->qfileName == NULL)
        {
            fprintf(stderr, "Query file specification (-q) is required for query alignment.\n");
            printUsageStringAbort();
        }
    }


    //////
    // Set up values for filenames.
    //////

    BOOL genomeIsFastaFile = FALSE;
    if (query)
    {
        // Here we need to create the .nib2 genome file name from the index file name.
        AAs->gfileName = (char *)malloc(strlen(AAs->xfileName)+10);
        strcpy(AAs->gfileName, AAs->xfileName);
        char * extLoc = strrchr(AAs->gfileName, '.');
        if (extLoc == NULL) fatalError("Specified index filename has improper or missing file extension.  Is it an index file?");
        strcpy(extLoc, ".nib2");
    }
    else
    {
        // We are given a genome file name.
        // Do some processing on it.

        // Get the extension of the gfile.
        char * gextActual = strrchr(AAs->gfileName, '.');

        // Check to see what type of genome file we have been given.
        if      (gextActual != NULL && (streq(gextActual, ".fna") || streq(gextActual, ".fa") || streq(gextActual, ".fasta"))) genomeIsFastaFile = TRUE;
        else if (gextActual != NULL && streq(gextActual, ".nib2")) genomeIsFastaFile = FALSE;
        else    fatalError("Expecting a \".fa\", \".fna\", \".fasta\", or \".nib2\" genome file.");

        // Also, empty the query file name.
        AAs->qfileName = NULL;
    }

    // Determine the output file name.
    if (compress && !genomeIsFastaFile)
    {
        fatalError("Expecting a \".fa\", \".fna\", or \".fasta\" genome file.");
    }
    else if (genomeIsFastaFile && (query || uncompress))
    {
        fatalError("Expecting a \".nib2\" genome file.");
    }
    else if (uncompress)
    {
        // Now generate the correct output file name.
        AAs->ofileName = (char *)malloc(strlen(AAs->gfileName)+10);
        strcpy(AAs->ofileName, AAs->gfileName);
        strcpy(strrchr(AAs->ofileName, '.'), ".fasta");
    }
    else if (compress)
    {
        AAs->ofileName = (char *)malloc(strlen(AAs->gfileName)+10);
        strcpy(AAs->ofileName, AAs->gfileName);
        strcpy(strrchr(AAs->ofileName, '.'), ".nib2");
    }
    else if (query && AAs->ofileName == NULL)
    {
        // Check if we are aligning queries, but have no specified output file or format.
        AAs->outputBlast8 = FALSE;
        AAs->outputSAM = TRUE;
        AAs->hardClip = TRUE;
        // Output to stdout  by default.
        AAs->ofileName = strdup("stdout");
    }

    // In case any of the default values of alignment arguments are based on other arguments.
    postProcessAlignmentArgs(AAs, query);

    if (index)
    {
        // Now generate the correct index file name.
        // Put the wordlen and skipDistance into the xfile extension,
        //     so we can have more than one index per genome file.
        char * oext = (char *)(malloc(20));
        sprintf(oext, ".X%02d_%02d_%05dS", AAs->wordLen, AAs->skipDist, AAs->maxHits);
        AAs->xfileName = (char *)malloc(strlen(AAs->gfileName)+20);
        strcpy(AAs->xfileName, AAs->gfileName);
        strcpy(strrchr(AAs->xfileName, '.'), oext);
        free(oext);
    }

    if (AAs->verbose) printAlignmentArgs(AAs, stderr);

    // Now that we have processed the input arguments, start the clock to time the operation.
    time_t startTime, endTime;
    startTime = time(NULL);

    if (compress)
    {
        compressFile(AAs);
        // The gfilename came in the parameter list, so should not be freed.
        AAs->gfileName = NULL;
    }
    else if (uncompress)
    {
        uncompressFile(AAs);
    }
    else if (query)
    {
        processQueryFile(AAs);
    }
    else if (index)
    {
        if (AAs->wordLen > 15)
        {
            fatalError("Word Length (-L) for index creation is currently restricted to < 16.");
        }
        if (AAs->skipDist < 1 || AAs->skipDist > AAs->wordLen)
        {
            fatalError("Skip Distance (-S) for index creation must be between 1 and WordLength (inclusive).");
        }

        // Creating an index.
        // Check if we should do a compression first.
        if (genomeIsFastaFile)
        {
            // We need to compress.
            // In this case, the gfilename is a fasta file, and ofileName will need a .nib2 extension.
            AAs->ofileName = (char *)malloc(strlen(AAs->gfileName)+10);
            strcpy(AAs->ofileName, AAs->gfileName);
            strcpy(strrchr(AAs->ofileName, '.'), ".nib2");
            if (fileNewerThan(AAs->gfileName, AAs->ofileName))
            {
                fprintf(stderr, "Compressing %s into %s.\n", AAs->gfileName, AAs->ofileName);
                compressFile(AAs);
                fprintf(stderr, "Finished compressing %s, now forming index.\n", AAs->ofileName);
                // The nib2 file we just created is now the genome input file.
            }
            else
            {
                fprintf(stderr, "%s already exists.  Creating the index file.\n", AAs->ofileName);
            }
            AAs->gfileName = AAs->ofileName;
            AAs->ofileName = NULL;
            if (AAs->verbose) printAlignmentArgs(AAs, stderr);
        }
        fprintf(stderr, "Creating index file %s.\n", AAs->xfileName);
        indexFile(AAs);
        fprintf(stderr, "Index %s created.\n", AAs->xfileName);
        free(AAs->xfileName);
        // We didn't create the genome file name, so we can't free it.
        if (!genomeIsFastaFile) AAs->gfileName = NULL;
    }
    else
    {
        fprintf(stderr, "Internal Error processing parameters.\\");
        printAlignmentArgs(AAs, stderr);
        abort();
    }

    // Output the timing information.
    struct rusage usagebuf;
    getrusage(RUSAGE_SELF, &usagebuf);
#ifndef COMPILE_USER_MODE
    fprintf(stderr, "Number of non-IO page faults: %zd, Number of IO page faults: %zd.\n", usagebuf.ru_minflt, usagebuf.ru_majflt);
#endif
    endTime = time(NULL);
    char * fname = AAs->gfileName;
    if (query) fname = AAs->qfileName;
    fprintf(stderr, "Operation on %s used ", fname);
    fprintTimeMicroSeconds(stderr, tv2usec(&usagebuf.ru_utime), 3);
    fprintf(stderr, " User, ");
    fprintTimeMicroSeconds(stderr, tv2usec(&usagebuf.ru_stime), 3);
    fprintf(stderr, " System, Total: ");
    fprintTimeMicroSeconds(stderr, (tv2usec(&usagebuf.ru_utime) + tv2usec(&usagebuf.ru_stime)), 3);
    fprintf(stderr, " in ");
    fprintTimeSeconds(stderr, (endTime-startTime), 0);
#ifdef COMPILE_USER_MODE
    fprintf(stderr, " wall time.\n");
#else
    fprintf(stderr, " wall time");
    if (query && AAs->numThreads > 1)
    {
        double totalTime = (tv2usec(&usagebuf.ru_utime) + tv2usec(&usagebuf.ru_stime))/1000000.0;
        double wallTime = (endTime-startTime);
        double speedUp = totalTime/wallTime;
        fprintf(stderr, " with thread speedup of %4.2f.\n", speedUp);
    }
    else fprintf(stderr, ".\n");
#endif

    // Clean up.
    disposeAlignmentArgs(AAs);

    return 0;
}
