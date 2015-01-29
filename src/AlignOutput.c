/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    AlignOutput.c    Includes code for output of alignments in various formats.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/


#include <stdlib.h>
#include "Math.h"
#include "Math.inl"
#include "FragsClumps.inl"
#include "SW.inl"

static inline int fputchar(FILE * out, char outchar)
{
    return putc_unlocked(outchar, out);
}

void outputFileHeader(AlignmentArgs_t * AAs)
{
    // Of the current output formats, only SAM has a header.
    if (!AAs->outputSAM) return;

    FILE * out = AAs->outFile;
    flockfile(out);
    // Start with the Header record.
    fprintf(out, "@HD\tVN:1.0\n");
    // Output the sequence dictionary.
    BaseSequences_t * BSs = AAs->BSs;
    for (int i=0; i<BSs->curCount; i++)
    {
        BaseSequence_t * BS = BSs->sequences[i];
        fprintf(out, "@SQ\tSN:");
        printChars(out, BS->name, 0, BS->nameLen);
        fprintf(out, "\tLN:%u\n", BS->length);
    }
    // Output the program name and command line.
    // The command line won't be what was actually put in, but the fully specified version of the parameters.
    fprintf(out, "@PG\tID:YAHA\tVN:0.1.%d\tCL:yaha", BUILDNUM);

    // I/O args
    // fprintf(out, " -g %s", AAs->gfileName);
    fprintf(out, " -q %s", AAs->qfileName);
    fprintf(out, " -x %s", AAs->xfileName);
    fprintf(out, " -os%c", AAs->hardClip ? 'h' : 's');
    fprintf(out, " %s", AAs->ofileName);

    // Thread arg.
    fprintf(out, " -t %d", AAs->numThreads);

    // General Alignment Args.
    fprintf(out, " -BW %d", AAs->bandWidth);
    // fprintf(out, " -DC %4.2f", AAs->DCost);
    fprintf(out, " -G %d", AAs->maxGap);
    fprintf(out, " -H %d", AAs->maxHits);
#ifndef COMPILE_USER_MODE
    fprintf(out, " -I %d", AAs->maxIntron);
#endif
    fprintf(out, " -M %d", AAs->minMatch);
    fprintf(out, " -MD %d", AAs->maxDesert);
    fprintf(out, " -P %4.2f", AAs->minIdentity);
#ifndef COMPILE_USER_MODE
    fprintf(out, " -R %d", AAs->minRawScore);
#endif
    fprintf(out, " -X %d",   AAs->XCutoff);

    // Affine Scoring Args.
    if (AAs->affineGapScoring)
    {
        fprintf(out, " -AGS Y");
        fprintf(out, " -GEC %d", AAs->GECost);
        fprintf(out, " -GOC %d", AAs->GOCost);
        fprintf(out, " -MS %d",  AAs->MScore);
        fprintf(out, " -RC %d",  AAs->RCost);
    }
    else fprintf(out, " -AGS N");

    // OQC and FBS Args.
    if (AAs->OQC)
    {
        fprintf(out, " -OQC Y");
        // fprintf(out, " -MSA %d", AAs->FBSMaxSimilar);
        fprintf(out, " -BP %d", AAs->BPCost);
        fprintf(out, " -MGDP %d", AAs->maxBPLog);
        fprintf(out, " -MNO %d", AAs->OQCMinNonOverlap);
        if (AAs->FBS)
        {
            fprintf(out, " -FBS Y");
            fprintf(out, " -PRL %4.2f", AAs->FBS_PSLength);
            fprintf(out, " -PSS %4.2f", AAs->FBS_PSScore);
        }
        else fprintf(out, " -FBS N");
    }
    else
    {
        fprintf(out, " -OQC N");
    }
    fputchar(out, '\n');
    funlockfile(out);
}

// This assumes the file is already locked.
// This allows us to keep all alignments from a query to be contiguous in the output file.
int printClump(Clump_t *clump, AlignmentArgs_t * AAs, QueryState_t * QS)
{
    double percent = .8;

    // Get the first and last fragments, and calculate preliminary reference offset range.
    SFragment_t * SFrag0 = SFragListPtr(clump)->head;
    SFragment_t * SFragn = SFragListPtr(clump)->tail;
    Fragment_t * frag0 = getFragPtr(SFrag0);
    Fragment_t * fragn = getFragPtr(SFragn);
    ROFF seqStartRefOff = frag0->startRefOff;
    ROFF seqEndRefOff = fragEndRefOff(fragn);

    // Now make sure we are within one base sequence
    // TODO Perhaps instead of bailing, we should split and check both ends separately.
    BaseSequence_t * BS = findBaseSequence(AAs->BSs, seqStartRefOff);
    if (!BS || (seqEndRefOff >= BS->startingOffset + BS->length))
    {
#ifndef COMPILE_USER_MODE
        fprintf(stderr, "query not within a base sequence.\n");
#endif
        return 0;
    }
    // Renormalize the start and end reference offsets to the base sequence.
    seqStartRefOff -= BS->startingOffset;
    seqEndRefOff -= BS->startingOffset;
    char * refPtr = AAs->basePtr;
    char * queryBuf = getClumpQueryBuf(QS, clump);

    // Now we can start the output.
#ifdef QUERYSTATS
    QS->alignOutCount += 1;
#endif

    FILE * out = AAs->outFile;

    // Output the SAM format.
    // This is the most complex of the output formats.
    // In particular, the CIGAR and MD strings together include all information in the "edit" sequence
    // to turn the query into the refernece.
    if (AAs->outputSAM)
    {
        printChars(out, QS->queryID, 0, QS->queryIDLen);
        fprintf(out, "\t%d\t", isReversed(clump) ? 0x10 : 0x00);
        printChars(out, BS->name, 0, BS->nameLen);
        fprintf(out, "\t%u\t%u\t",  seqStartRefOff + 1, clump->mapQuality);

        // We will just construct the CIGAR string here.
        EditOpList_t * list = EOLPtr(clump);

        // For now I have just translocated this code from the alignment function.
        // TODO I don't see why it is worth adding this to the editoplist at all!
        // If we get here, we have a real alignment that will not be split, nor thrown out.
        // So, it if finally safe to put in the clip regions.
        int clip = QS->queryLen - 1 - frag0->endQueryOff;
        if (clip > 0) addEditOpToBack(list, (AAs->hardClip ? HardClipEditOp : SoftClipEditOp), clip);
        clip = frag0->startQueryOff;
        if (clip > 0) addEditOpToFront(list, (AAs->hardClip ? HardClipEditOp : SoftClipEditOp), clip);

        int matches = 0;
        forAllEditOpsInList(item, list)
        {
            // CIGAR does not distinguish between match and replace!!
            if ((item->opcode == MatchEditOp || item->opcode == ReplaceEditOp))
            {
                matches += item->length;
                continue;
            }
            if (matches > 0)
            {
                fprintf(out, "%dM", matches);
                matches = 0;
            }
            fprintf(out, "%d%c", item->length, EditOpCode2char(item->opcode));
        }
        // Check if there are left over matches.
        if (matches > 0) fprintf(out, "%dM", matches);

        // Output MRNM, MPOS, and ISIZE
        fputs("\t*\t0\t0\t", out);
        // Output the query sequence.
        int qstart = 0;
        int qend = QS->queryLen - 1;
        if (AAs->hardClip)
        {
            qstart = frag0->startQueryOff;
            qend   = fragn->endQueryOff;
        }
        for (int i=qstart; i<=qend; i++) fputchar(out, queryBuf[i]);
        fputchar(out, '\t');
        // Output the qual string if we have the qual information.
        // For reversed strand, we have to reverse the string.
        if (AAs->fastq)
        {
            if (isReversed(clump))
                for (int i=qend; i>=qstart; i--) fputchar(out, QS->qualBuf[i]);
            else
                for (int i=qstart; i<=qend; i++) fputchar(out, QS->qualBuf[i]);
        }
        else fputchar(out, '*');
        fputchar(out, '\t');

        // Output TAG based information.
        // First the AS tag.
        fprintf(out, "AS:i:%d\t", clump->totScore);
        // The NM tag.
        fprintf(out, "NM:i:%d\t", clump->gapBases + clump->mismatchedBases);
        // TODO  If we know how many hits we have for the query, we can output H0, H1, H2, NH, IH, and HI.
        // The MD tag.
        fprintf(out, "MD:Z:");
        // Figure out the MD tag value.
        matches = 0;
        EditOpCode previous = UnknownEditOp;
        ROFF currRefOff = frag0->startRefOff;
        forAllEditOpsInList(item, list)
        {
            EditOpCode code = item->opcode;
            // This should never happen.
            if (code == previous)
            {
                fprintf(stderr, "Two identical codes in a row in EditOpList\n");
                printEditOpList(list);
                fprintf(stderr, "\n");
                fprintQueryId(QS, stderr);
                fprintf(stderr, "\n");
                abort();
            }
            // The ugly handling of matches is to allow for combine of matches across an Insert!!
            // This is what SAMTOOLS does for example!!
            if      (code == MatchEditOp)
            {
                matches += item->length;
                currRefOff += item->length;
            }
            else if (code == ReplaceEditOp)
            {
                if (matches > 0)
                {
                    fprintf(out, "%d", matches);
                    matches = 0;
                }
                // Agreed upon hack to distinguish these mismatches from any deleted reference bases before this.
                if (previous == DeleteEditOp) fputchar(out, '0');
                for (int i=0; i<item->length; i++) fputchar(out, getCharFrom4Code(refPtr, currRefOff + i));
                currRefOff += item->length;
            }
            else if (code == DeleteEditOp)
            {
                if (matches > 0)
                {
                    fprintf(out, "%d", matches);
                    matches = 0;
                }
                fputchar(out, '^');
                for (int i=0; i<item->length; i++) fputchar(out, getCharFrom4Code(refPtr, currRefOff + i));
                currRefOff += item->length;
            }
            previous = code;
        }
        if (matches > 0) fprintf(out, "%d", matches);
        // We will output the clump status information as a tag field.
        // In keeping with SAM format convention, any tag starting with X, Y, or Z is fair game.
        // We will use YS for "Yaha Status".
        // We will use a hex field so that it is easy to see the bits.
        fprintf(out, "\tYF:H:%02X", clump->status);
        // If we are using OQC, we will output the number of primary alignemnts.
        // And the position number of this alignment along the query (from left to right).
        if (AAs->OQC)
        {
            fprintf(out, "\tYI:i:%d", clump->matchedPrimary);
            fprintf(out, "\tYP:i:%d", QS->primaryCount);
            // If this is a primary alignment, also output the number of secondaries.
            if (isPrimary(clump)) fprintf(out, "\tYS:i:%d", clump->numSecondaries);
        }
        fputchar(out, '\n');
    }

    // Output the Blast 8 style one line per alignment format.
    // This produces one line of output per alignment.
    // The fields are tab delimited.
    // Our output is somewhat non-standard.
    // The first 10 fields are standard.
    // They are: QID, RSN, %ID, A-Len, Mismatches, GapBases, QSO, QEO, RO1, RO2
    // The start and end numbers are 1 based.   The range is inclusive.
    // The query start and end are relative to the inverted query if strand is "-".
    // If the reverse query is matched, then the reference start and end offsets are reversed!
    // We add four more fields.
    // First, '+' or '-' for which strand matched.
    // This is redundant with comparing RO1 and RO2, but convenient for post processing the file.
    // Second, the raw score.
    // Third, the query length.
    // Fourth, the %coverage of the alignment on the total query.
    // We don't output either normalized score or e-value as in normal Blast8 output.
    if (AAs->outputBlast8)
    {
        printChars(out, QS->queryID, 0, QS->queryIDLen);
        fputchar(out, '\t');
        printChars(out, BS->name, 0, BS->nameLen);
        fprintf(out, "\t%4.2f\t%d\t%d\t%d", percent*100, clump->totLength, clump->mismatchedBases, clump->gapBases);
        if (isReversed(clump)) fprintf(out, "\t%d\t%d\t%d\t%d\t%c", QS->queryLen-fragn->endQueryOff,
                                       QS->queryLen-frag0->startQueryOff, seqEndRefOff+1, seqStartRefOff+1, '-');
        else fprintf(out, "\t%d\t%d\t%d\t%d\t%c", frag0->startQueryOff+1, fragn->endQueryOff+1, seqStartRefOff+1,
                     seqEndRefOff+1, '+');
        fprintf(out, "\t%d\t%d\t%4.2f\n", clump->totScore, QS->queryLen, ((double)clump->matchedBases/QS->queryLen)*100);
    }

    return 1;
}
