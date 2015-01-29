/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    Math.h          All of the Header information for the entire project.
                             The name of the file is unfortunate for historic purposes.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

// We want to use the Large File System, especially for memmory mapped IO of the index.
#define _FILE_OFFSET_BITS 64

#include <stdint.h>
#include <stdio.h>
#include <string.h>

// Some online forums claim that memset as a macro can cause problems in pre 4.3 versions of redhat.
// Who knows if this is true.
// I continue to use it nonetheless.
// #undef memset

#include <assert.h>

// Control the usage of extra parameters for our internal investigations.
// Versus having a user friendly set of program options.
// However, we now set this option in the Makefile, so we don't have to change code to decide.
// #define COMPILE_USER_MODE

// Max length of a query id.
// Longer query names will be truncated.
#define MAX_QUERY_ID_LEN 200

// Standard Hack to tell C++ to leave my C function names alone!!!
#ifdef __cplusplus
extern "C" {
#endif

// Some basics for bools in case we go to C++
typedef uint8_t BOOL;
#define TRUE (1)
#define FALSE (0)

// Algorithmic switches

// This controls whether index creation will randomly sample the available hits to keep all k-mers within maxHits.
#define RANDOM_SAMPLE_MAXHITS

// Some debugging aids.
// This controls output statistics about the query file.
// TODO this should be moved to be controlled by the -v option at invocation.
// #define STATS
// This controls whether we output statistics about alignemnts per query.
// #define QUERYSTATS
// This controls some debug output and assert.
// #define DEBUG
// Use this if you REALLY want a lot of output!
// #define ALLDEBUG

#ifdef DEBUG
void outputDebugString(char * str);
#define ASSERT(exp,str) if (!(exp))outputDebugString(str);
#else
#define ASSERT(exp,str) ;
#endif

// Version constants
#define CURRENT_INDEX_FILE_VERSION (-1)

// Type used for file sizes and offsets.
// Needs to be 64 bits if file sizes are over 2GB.
typedef size_t   FSIZE;
typedef ssize_t SFSIZE;
// Use this type for filedescriptors.
typedef int FDES;
typedef int64_t SINT64;
#define SINT64_MAX_VALUE 0X7FFFFFFFFFFFFFFF
// Type for unsigned ints.
typedef uint64_t UINT64;
#define UINT64_MAX_VALUE 0XFFFFFFFFFFFFFFFF
// Needs to be bigger if using words over 16 bases long for hash.
typedef uint32_t UINT;
#define UINT_MAX_VALUE 0XFFFFFFFF
// Use this type for smaller unsigned ints, such as counts of words.
typedef uint16_t SUINT;
#define SUINT_MAX_VALUE 0XFFFF
// Use this type for smaller signed ints, such as alignment scores.
typedef int16_t SINT;
#define SINT_MAX_VALUE 0x7FFF
// Use this type for byte arrays.
typedef uint8_t UBYTE;
// Type used for offsets into the ROA and entries in the ROA.
typedef UINT ROFF;
#define ROFF_MAX_VALUE UINT_MAX_VALUE
// Type for Query Offsets.
typedef SUINT QOFF;
#define QOFF_MAX_VALUE SUINT_MAX_VALUE
typedef int DPScore;

//////////
//     Some generally useful functions/macros.
//////////
void   emptyFunction (int arg);
UINT   iexp(int base, int power);
UINT64 iexp64(int base, int power);
int    ilog(int base, UINT num);
// TODO: If the args are expressions, they may get evaled more than once?
#define roundUpPow2(num, pow) (((num + ((1 << pow)-1)) >> pow) << pow);
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a)   (((a) >= 0) ? (a) : -(a))

/////
//    Random Number Generation
//    Most of these are inlines, so look in Math.inl as well.
/////

struct randState
{
    UINT state[5];
};
typedef struct randState randState_t;

randState_t * makeRandState();
void disposeRandState(randState_t * rs);
void initRandStateDefault(randState_t * rs);
void initRandState(randState_t * rs,UINT * instate);
UINT getRandBits(randState_t * rs);
double getRandDouble(randState_t * rs);
UINT getRandUInt(randState_t * rs, UINT start, UINT end);
void getRandSample(randState_t * rs, UINT * input, int inLen, UINT * output, int outLen);

//////////
//     Bit Arrays are used as a memory efficient representation of sets.
//////////
struct bitArray
{
    int numOfBits;
    int numOfBytes;
    uint8_t * contents;
};
typedef struct bitArray BitArray_t;

BitArray_t * makeBitArray(int bitCount);
void disposeBitArray(BitArray_t *bitArr);
void setAll(BitArray_t *bitArr, BOOL value);
BOOL get(BitArray_t *bitArr, int offset);
BOOL set(BitArray_t *bitArr, int offset, BOOL value);
// void setTRUE(BitArray_t *bitArr, int offset);


//////////
//     Binary Heap is primary used for merge of ROA segments during seed match sort.
//////////
// The contents array and size (current size) is 1 based leaving the zeroeth element unused.
// This simplifies the parent child offset calculations.
/////////
struct heap
{
    UINT64 * contents;
    int      size;
    int      maxSize;
};
typedef struct heap Heap_t;

Heap_t * makeHeap(int maxSizeArg);
void disposeHeap(Heap_t * heap);

//////////
//     Signature for File Operations.
//////////
void fatalError(const char * str);
int  findChar(char * str, FSIZE offset, char search);
void printChars(FILE * out, char * str, int start, int len);
void fprintTimeSeconds (FILE * out, double seconds, int precision);
void fprintTimeMicroSeconds (FILE * out, UINT64 microSeconds, int precision);
void checkFSerr (SFSIZE returnCode);
BOOL fileNewerThan (char * file1Name, char * file2Name);
FDES openForRead(char * filename, FSIZE *fileSize, char **contentPtr, BOOL preload);
void closeForRead(FDES ifp, FSIZE fileSize, char *contentPtr);
FILE * openForSeqRead(char * filename);
void closeForSeqRead(FILE * inFile);
char * blockMalloc(FSIZE size, char * string);
void blockFree(char * ptr, FSIZE size);
void blockWrite(FDES ofp, char * ptr, FSIZE size);
FDES openForWrite(char * filename);
void closeForWrite(FDES ofp);
FILE * openForPrint(char * filename);
void closeForPrint(FILE * ofile);

//////////
//     Signature for Compression and DeCompression Operations.
//     These are now mostly inlines in Math.inl, so not much left here.
//////////

char getFrom4Code(char * seq, ROFF offset);

// These are used for debug output.
void decode (UINT hash, char* seq, int len);
void output4Seq(FILE * stream, char* seq, ROFF offset, int len);
void check4Seq(FILE * stream, char * hashSeq, char* seq, ROFF offset, int len);

//////////
//     Structure Definitions and related Functions for Alignment Data Structures.
//////////

//////////
//     Structure for a base sequence/chrom (typically from a file)
//////////
struct baseSequence
{
    char * name;
    SUINT  nameLen;
    ROFF   startingOffset;
    ROFF   length;
};
typedef struct baseSequence BaseSequence_t;

BaseSequence_t * makeBaseSequence();
void disposeBaseSequence(BaseSequence_t * BS);
void printBaseSequence(BaseSequence_t * BS, FILE * out);

// Structure for a collection of base sequences (typically from a file).
struct baseSequences
{
    BaseSequence_t ** sequences;
    int maxCount;
    int curCount;
    char * basePtr;
    UINT totalBaseCount;
};
typedef struct baseSequences BaseSequences_t;

BaseSequences_t * makeBaseSequences(int maxCountArg);
void addBaseSequence(BaseSequences_t * BSs, BaseSequence_t * BS);
BaseSequence_t * findBaseSequence(BaseSequences_t * BSs, ROFF offset);
int  findBaseSequenceNum(BaseSequences_t * BSs, ROFF offset);
ROFF baseSequencesMaxROff(BaseSequences_t * BSs);
void disposeBaseSequences(BaseSequences_t * BSs);
void printBaseSequences(BaseSequences_t * BSs, char * gfilename, FILE * out);
void outputBaseSequences(BaseSequences_t * BSs, FDES ofp);
BaseSequences_t * loadBaseSequences(char * gfilePtr);
void normalizeBaseSequences(BaseSequences_t * BSs);

//////////
//     Structure for the arguments to various alignment functions.
//     This will contain all of the parameters to either the index creation and/or the query process.
//////////
struct alignmentArgs
{
    // Filenames
    char * gfileName;        // The reference genome input file name.
    char * xfileName;        // The index file name.
    char * qfileName;        // The query file name.
    char * ofileName;        // The output file name.  (.nib2 if compress, otherwise, anything at all!)
    FILE * outFile;          // The output file pointer to open output file.

    // BaseSequences and a couple of cached values to keep from referencing the BS structure very often.
    BaseSequences_t * BSs;   // Base sequences structure read from the reference file.
    char *basePtr;           // Pointer to the base of the reference portion of the NIB2 file.

    // Controls the number of threads used.
    int   numThreads;

    // These are index related information, therefore read-only and sharable.
    ROFF  maxROff;           // Maximum reference offset that is within the reference.
    ROFF *startingOffs;      // Pointer to the SO portion of the index.
    ROFF *ROAPtr;            // Pointer to the ROA portion of the index.

    BOOL fastq;              // Is the input query file in fastq format?

    // Index related parameters.
    int   wordLen;           // The length of a k-mer to match.  Default is 8.
    int   skipDist;          // Distance between k-mer starting offsets for index.  Default is 1.  Query always uses 1.
    int   maxHits;           // Maximum hits allowed for a k-mer before information is thrown away.  Default is 1024.

    // General Alignment parameters.
    int   bandWidth;         // BandWidth to use for banded SW calls.  This is the band size on each side of the diagonal.
    int   maxIntron;         // Maximum allowed gap on the refernce.
    int   minMatch;          // Minimum number of k-mers needed on diagonal to trigger alignment.
    int   maxGap;            // Maximum allowed gap on the query, and in gapped extensions.
    int   maxDesert;         // Longest stretch between fragments on the diagonal without seed hits to include when putting clumps together.
    int   minNonOverlap;     // Minimum number of bases for a fragment that don't overlap a current alignment to reconsider a second alignment.
    float minIdentity;       // Percentage identity threshold for matching bases from start to end of alignment.  Default is .90
    int   minRawScore;       // Minimum raw score threshold.  Typical value is (wordLen+ ((minMatch-1) * skipdist)).

    // These are for Affine Gap Scoring.
    // All should be positive numbers.  The addition and subtraction is done properly in the algorithms.
    // AGS is on by default.
    BOOL  affineGapScoring;  // Boolean that controls if we are using affine gap scoring.
    UBYTE minExtLength;      // Number of bases an extension must be before DP is necessary.  minExtLength = Ceiling(RCost/MScore) + 1
    int   GOCost;            // The cost of opening a gap.
    int   GECost;            // The cost of extending a gap.
    int   RCost;             // The cost of a base replacement.  This takes the place of a replacement matrix which seems overkill for DNA.
    int   MScore;            // The score for a base match.
    int   XCutoff;           // X cutoff value for extensions, ala BLAST.

    // These are for Filtering of Alignment by Similarity to a set of query spanning alignments.
    // OQC is on by default.
    BOOL  OQC;               // Controls use of Optimal Query Coverage Algorithm.
    int   OQCMinNonOverlap;  // Number of bases that a clump must extend beyond a given alignment for a second alignment to be considered.
    BOOL  useAccOverlap;     // Controls use of more expensive, but more accurate overlap score calculations during OQC. (NO LONGER NEEDED.)
    int   BPCost;            // Cost for a breakpoint in a graph path.
    int   maxBPLog;          // Maximum allowed breakpoint distance as log base 10.  So, 7 is 10Mb max, while 1 means just the BPCost will be used.
    BOOL  FBS;               // Bool that controls output of secondary alignments.
    float FBS_PSLength;      // Percent of reciprocal overlap between two alignemnts to be considered similar.
    float FBS_PSScore;       // Percent of "best" score needed by secondary alignment to be considered similar.
    //int   FBSMaxSimilar;   // Maximum number of similar alignments to report.
                             // Currently, this is only tested with OQC on (After), and to see if it != 1.
                             // If one, only primaries are output.
                             // Otherwise, all secondaries are output.
                             // Eventually, this is supposed to put a limit on secondary alignments output with value > 1.

    int  maxQueryLength;     // The length of the maximum allowed query in the input file.  Default is 1300.
    BOOL verbose;            // Controls the reporting detail for run information.
    BOOL outputBlast8;       // Controls query output in blast8 file format.
    BOOL outputSAM;          // Controls query output in SAM file format.
    BOOL hardClip;           // Controls hard vs soft clipping of query read reporting in SAM output.

#ifdef QUERYSTATS
    char * qsfileName;       // The file name for the query stats output file.
    FILE * qsFile;           // The file pointer for same.
    BOOL   queryStats;       // A boolean to say that we are collecting query stats.
#endif
};
typedef struct alignmentArgs AlignmentArgs_t;

AlignmentArgs_t * makeAlignmentArgs ();
void disposeAlignmentArgs(AlignmentArgs_t * AAs);
void postProcessAlignmentArgs(AlignmentArgs_t * AAs, BOOL query);
void printAlignmentArgs(AlignmentArgs_t * AAs, FILE * out);

// Pre-declaration of the type of the queryState structure to be used in many other places before its definition below.
// Since the structure definitions are circular, we need to pre-declare something.
typedef struct queryState QueryState_t;

//////////
//     Dynamic Programming Interfaces.
//////////

// Data structure for storing lists of edit operations.

// Use these in place of an enum.
// The translation of these to the output char is then trivial!
typedef char  EditOpCode;
#define UnknownEditOp  'U'
#define HardClipEditOp 'H'
#define SoftClipEditOp 'S'
#define DeleteEditOp   'D'
#define InsertEditOp   'I'
#define MatchEditOp    'M'
#define ReplaceEditOp  'R'

// We will keep the offset of both the head and the tail.
// So that both add to front and back will fast.
// We have now gone to a linked list of sorts.
// But to avoid large ptrs (especially in 64-bit)
//     and lots of small allocs, we will still
//     use one big global pre-allocated array and use offsets for ptrs.
// Therefore, EditOpLists do not own their own storage.
typedef UINT  EditOpPtr;
#define EditOpPtrNULL ((EditOpPtr) 0)
struct editOp
{
    SUINT      length;
    EditOpCode opcode;
    EditOpPtr  next;
    EditOpPtr  prev;
};
typedef struct editOp EditOp_t;

struct editOpList
{
    QueryState_t * QS;      // Ptr to the QueryState that holds the allocation ptrs.
    EditOpPtr      head;    // Ptr to the head element.
    EditOpPtr      tail;    // Ptr to the tail element.
    // EditOpPtr  next;     // Ptr to next edit op in iterator.
};
typedef struct editOpList EditOpList_t;

EditOpList_t * makeEditOpList();
void           disposeEditOpList(EditOpList_t * list);
void           DPInit(QueryState_t * QS);
void           DPCleanUp(QueryState_t * QS);
void           emptyEditOpList(EditOpList_t * list);
void           mergeEOLToFront (EditOpList_t * target, EditOpList_t * source);
void           mergeEOLToBack (EditOpList_t * target, EditOpList_t * source);
void           splitEditOpListBefore (EditOpList_t * source2head, EditOpList_t * newTail, EditOp_t * node);
void           splitEditOpListAfter (EditOpList_t * source2head, EditOpList_t * newTail, EditOp_t * node);
void           printEditOpList (EditOpList_t * list);
void           validateEditOpList(EditOpList_t * list);

int            findAGSAlignment(QueryState_t * QS, ROFF rOff, QOFF rLen, char * qStr, QOFF qOff, QOFF qLen, EditOpList_t * list);
int            findAGSAlignmentBanded(QueryState_t * QS, ROFF rOff, QOFF rLen, char * qStr, QOFF qOff, QOFF qLen, EditOpList_t * list);
int            findAGSForwardExtension(QueryState_t * QS, ROFF rOff, char * qStr, QOFF qOff, QOFF qLen,
                                       EditOpList_t * list, QOFF * addedQLen, QOFF * adedRLen);
int            findAGSForwardExtensionCarefully(QueryState_t * QS, ROFF rOff, char * qStr, QOFF qOff, QOFF qLen,
                                                EditOpList_t * list, int score, QOFF * addedQLen, QOFF * addedRLen);
int            findAGSBackwardExtension(QueryState_t * QS, ROFF rOff, char * qStr, QOFF qOff, QOFF qLen,
                                        EditOpList_t * list, QOFF * addedQLen, QOFF * adedRLen);
int            findAGSBackwardExtensionCarefully(QueryState_t * QS, ROFF rOff, char * qStr, QOFF qOff, QOFF qLen,
                                                 EditOpList_t * list, int score, QOFF * addedQLen, QOFF * addedRLen);
BOOL           EditOpList2Maxmatch(EditOpList_t * list, int min);
int            EditOpList2AGS(EditOpList_t * list, AlignmentArgs_t * AAs);


//////////
//     Data Structures and Related Functions for Fragments, Clumps, and the Alignment Routines.
//////////


//////////
//     Offset Count
//////////
//     Used to keep track of the hits for a seed while we are forming fragments.
struct offsetCount
{
    ROFF  sOffset;
    SUINT count;
    SUINT newCount;
};
typedef struct offsetCount OffsetCount_t;

//////////
//     Fragment
//////////
// A fragment is a more or less perfect partial alignment of a query on the reference.
// It contains no insertions or deletions
// For now we won't store enough information to recreate the total alignment.
// Instead we will focus on length and score.
// It is debatable whether or not we should keep track of score as we go, or calculate in a post pass.

// 11/8/2010  Changing meaning of hitCount.
// It will now be the number of bases that match.
// We now always form fragments from contiguous runs of seed matches.
// Therefore, the number of matching bases is just the length of the fragment.
// TODO.  hitCount is now redundant with the fragment length.
//        We should consider eliminating it and refLen, and keep endRefOff directly.

struct fragment
{
    ROFF  startRefOff;
    SUINT startQueryOff;
    SUINT endQueryOff;
    SUINT hitCount;
    SUINT refLen;
};
typedef struct fragment Fragment_t;

void resetFragment(Fragment_t *frag);
void copyFragment(Fragment_t *fragIn, Fragment_t * fragOut);
int  extendFragmentForwardToStopPerfectly(Fragment_t *frag, char *refPtr, char *qPr, int len);
int  extendFragmentBackwardToStopPerfectly(Fragment_t *frag, char *refPtr, char *qPtr, int len);

//////////
//     Super Fragment.
//////////
//     Includes a fragment, but also has an editoplist for the fragment, and a score.
//     These are stored in a linked list for an alignment (and in the allocation pool.
//////////
typedef struct superFragment SFragment_t;
struct superFragment
{
    Fragment_t      frag;
    int             score;
    EditOpList_t    EOList;
    SFragment_t   * next;
    SFragment_t   * prev;
};

SFragment_t  * makeSFragment(QueryState_t * QS);
void           disposeSFragment(SFragment_t * SFrag, QueryState_t * QS);

struct sFragmentList
{
    SFragment_t * head;
    SFragment_t * tail;
    SFragment_t * next;
};
typedef struct sFragmentList SFragmentList_t;

// There are no calls to the first two of these, as the lists are all embedded in clumps.
SFragmentList_t * makeSFragmentList();
void              disposeSFragmentList(SFragmentList_t * list, QueryState_t * QS);
void              emptySFragmentList(SFragmentList_t * list, QueryState_t * QS);
void              removeSFragment(SFragmentList_t * list, SFragment_t * SFrag);
SFragment_t *     findFollowingSFrag(SFragmentList_t * list, SUINT queryOff);
SFragment_t *     getFirstSFrag(SFragmentList_t * list);
SFragment_t *     getNextSFrag(SFragmentList_t * list);
SFragment_t *     getLastSFrag(SFragmentList_t * list);
SFragment_t *     getPrevSFrag(SFragmentList_t * list);

//////////
//     Clump
////////////
//     A clump will be a collection of SFragments.
//     The collection will have different properties at differest stages of the alignment.
//     It will start as a bunch of seed strings potentially with space between them.
//     Later, the collection will represent an entire alignment with no space between them.
//     And finally, the list will be condensed into one SFragment for the entire alignment.
//     At no time will SFragments in the clump be allowed to overlap in either the reference or the query.
////////////
typedef struct clump Clump_t;
struct clump
{
    EditOpList_t      EOList;    // Pointer to the edit op array for the alignment.
    SFragmentList_t   SFragList; // Pointer to a linked list of super fragments.
    Clump_t         * next;      // Next clump in clump list.
    QOFF   totScore;              // Total Score for the clump.  Usually AGS score.
    // Why bother to cache this number?
    QOFF   totLength;             // Total alignment length.  Sum of matched, mismatched and gaps.
    QOFF   matchedBases;          // Count of matched bases in query and reference.  This is also the "raw" score.
    QOFF   mismatchedBases;       // Count of replacements.
    QOFF   gapBases;              // Sum of inserts and deletes.
    SUINT  numSecondaries;        // If this is a primary alignment, the number of matching secondary alignments.
    SUINT  matchedPrimary;        // The position number of the matching primary alignment along the query.
    UBYTE  status;                // Collection of flags that indicate status information about clump.
    UBYTE  mapQuality;            // Estimate of the mapping quality of the alignment.
};

Clump_t * makeClump(QueryState_t * QS);
void disposeClump(Clump_t * clump, QueryState_t * QS);
void resetClump(Clump_t * clump, QueryState_t * QS);
void addFragment(Clump_t * clump, Fragment_t * frag, QueryState_t * QS);
void insertFragment(Clump_t * clump, Fragment_t * frag, QueryState_t * QS);
void cleanUpClump(Clump_t * clump, QueryState_t * QS);
void extendClumpForwardReverse(Clump_t * clump, QueryState_t * QS);
void extendClumpForwardReverseCarefully(Clump_t * clump, QueryState_t * QS, BOOL goBack, BOOL goForw);
int  alignClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS);
int  scoreClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS);
void postProcessClumps(QueryState_t * QS, AlignmentArgs_t * AAs);
void outputFileHeader(AlignmentArgs_t * AAs);
int  printClumps(QueryState_t * QS, AlignmentArgs_t * AAs);
int  printClump(Clump_t * clump, AlignmentArgs_t * AAs, QueryState_t * QS);
SFragment_t * makeAndAlignSFragmentToFillGap(SFragment_t * SFrag1, SFragment_t * SFrag2, QueryState_t * QS, Clump_t * clump);
// Some debug output routines.
void printDebugClump(Clump_t * clump, AlignmentArgs_t * AAs, FILE * out);
void printDebugFragment(Fragment_t * frag);
void printDebugFragmentExtended(Fragment_t * frag, AlignmentArgs_t * AAs);

//////////
// Signatures for Query Helpers.
//////////
void initializeQueries (QueryState_t * QS);
void finalizeQueries (QueryState_t * QS);
int  findFragmentsSort (AlignmentArgs_t * AAs, QueryState_t * QS, int matchCount);
void processFragmentsGapped(AlignmentArgs_t * AAs, QueryState_t * QS, int fragCount);

void eliminateFragments(AlignmentArgs_t * AAs, QueryState_t * QS, Fragment_t * fragArray,
                        BitArray_t * usedFrags, int startFragNum, int endFragNum, Clump_t * clump);
void processFragmentRangeUsingGraph(AlignmentArgs_t * AAs, QueryState_t * QS, Fragment_t * fragArray,
                                    BitArray_t * usedFrags, int startFragNum, int endFragNum);
void postFilterBySimilarity(AlignmentArgs_t * AAs, QueryState_t * QS);
void postFilterRemoveDups(AlignmentArgs_t * AAs, QueryState_t * QS);

//////////
//     Slab Allocator decls used to hold thread local slab allocations in the querystate object.
//////////
typedef struct LBMallocBlock LBMallocBlock_t;
struct LBMallocBlock
{
    char *            block;     // Pointer to the payload block
    LBMallocBlock_t * next;      // Pointer to the next allocation block
    FSIZE             size;      // Size of the allocated block
};

char * pushNewLBMallocBlock(int blockSize, LBMallocBlock_t **blockArrayPtr);
void   freeLBMallocBlocks(LBMallocBlock_t * block);
void   cleanUpClumpStructures(QueryState_t * QS);
void   initClumpStructures(QueryState_t * QS);

/////
//     QueryState
/////
//     Store all of the thread local information for a query.
/////
// Forward declaration needed to define query state.
typedef struct fGraphNode fGraphNode_t;
struct queryState
{
    AlignmentArgs_t * AAs;        // Pointer to the alignment args (not owned by QS).

    char * queryID;               // Pointer to the start of the query id from the query file.
    char * queryBuf;              // Pointer to the start of the current query buffer, depending on reversed or not.
    char * forwardBuf;            // Pointer to the start of the forward buffer.  This will be allocated on heap.
    char * reverseBuf;            // Pointer to the start of the reverse buffer.  This will be allocated on heap.
    randState_t * randGen;        // Pointer to the state for the random number generator.
    SUINT  queryIDLen;            // The length of the query ID in chars.
    SUINT  queryLen;              // The length of the query in bases.
    char * qualBuf;               // Pointer to the start of the quality string for FASTQ query files.
                                  // We don't need a separate length, as it is the same len as the query.
                                  // We also won't bother to reverse this string for - strand, but just output backwards.
    SUINT  maxQueryLength;        // The maximum query length currently supported by the structures in this query state object.
                                  // This will always be <= the real maxQueryLength specified in the AAs object.
    BOOL   reversed;              // TRUE if this is the reverse complement of the original query in the input file.

    // We store both the char query string as well as its 4-bit encoding, one code per byte.
    char * forwardCodeBuf;        // Contains the forward query converted to 4-bit codes.  Allocated on heap.
    char * reverseCodeBuf;        // Contains the reverse query converted to 4-bit codes.  Allocated on heap.
    char * queryCodeBuf;          // Pointer to the start of the current code bufffer, depending on reverse or not.
    // The following two pointers allow us to build the reverse buffer while we are scanning the forward query.
    // In this case, we don't know the length of the query, so need to reset the main pointers to start of this query, not the buffer.
    char * reverseBufBase;        // Pointer to the base of the reverse buffer used for freeing.
    char * reverseCodeBufBase;    // Pointer to the base of the reverse code buffer used for freeing.

    // The clump list.
    // TODO: We should consider making this a real list with a tail pointer to make putting them back on the free list faster.
    Clump_t * clumps;             // A ptr to the head of a linked list of clumps.
    int    clumpCount;            // A count of the number of clumps in the linked list.
    int    primaryCount;          // A count of the number of primary alignments for this query.

    // This collections of structures are used to go from seed hits to clumps.
    BOOL * coverage;              // An array for storing which query positions are covered by the current clump set.
    OffsetCount_t * offsetCounts; // Array used during seed match merge. (strand).
    Fragment_t * fragArray;       // Array for those fragments.
    BitArray_t * usedFrags;       // Array parallel to the fragArray saying which frags have already been used.
    Heap_t * sortHeap;            // Binary heap used to merge ROA segments for sort of seed hits.
    fGraphNode_t * fragGraphNodes;// The array of frag graph nodes used to find best set of frags in reference region.
    Fragment_t * tempFragArray;   // Small frag array used to sort frags by decreasing size while adding to a clump.
    int fragGraphNodesLen;        // Size of frag graph array used to decide if we need to realloc to a larger size.
    int tempFragArrayLen;         // Size of small frag array used to decide if we need to realloc to a larger size.

    // Slab allocated objects.
    // We will lazily add slabs as they are needed.
    // Then only clean them up at the end of the program.

    // Clumps.
    LBMallocBlock_t * clumpBlockList;     // Ptr to head of linked list of allocated clump slabs.
    Clump_t *         freeClumpList;      // Ptr to head of linked list of free clump objects within the slabs.

    // SFragments.
    LBMallocBlock_t * SFragBlockList;     // Ptr to head of linked list of allocated SFragment slabs.
    SFragment_t *     freeSFragmentList;  // Ptr to head of linked list of free SFragment objects within the slabs.

    // EditOps
    EditOp_t * EditOpArray;               // Ptr to head of sole slab for editops.  This gets realloced if space runs out.
    EditOpPtr  EditOpFreeList;            // "Ptr" to the head of linked list of free EditOp objects with the slab.
    int        EditOpArrayMaxElements;    // Current allocated size of the EditOpArray;

    // DP/SW structures.
    char *     refTemp;                   // Ptr to decompressed refernece used during DP.
    EditOpCode * EOArray;                 // Ptr to the Edit Op array.
    DPScore    * IDCount;                 // Ptr to the count (length) for edit ops primarily used for I and D with AGS.
    DPScore    * PVRowSave;               // Ptr to the (previous) row of scores.
    DPScore    * PFRowSave;               // Ptr to the (previous) row of insert scores.
    DPScore    * PIRowSave;               // Ptr to the (previous) row of insert lengths.

#ifdef QUERYSTATS
    UINT   seedMatches;           // Keep track of the number of seed matches.
    UINT   alignCount;            // Keep the count of the alignments we tried.
    UINT   DPCount;               // Keep the count of the number of calls to DP.
    UINT   alignOutCount;         // Keep the count of the alignments we output.
    UINT   usec;                  // Keep track of the time spent in the alignment.
    // Next two used to calculate the avg size of the diagonal regions.
    UINT   diagRegionCount;       // Keep the count of the number of diagonal regions.
    UINT   diagRegionTotal;       // Keep the count of the total number of fragments in diagonal regions of interest.
#endif
};

QueryState_t * makeQueryState (AlignmentArgs_t * AAs);
void resetQueryState(QueryState_t * QS);
void addClump(QueryState_t * QS, Clump_t * clump);
void insertClumpAfter(QueryState_t * QS, Clump_t * prevClump, Clump_t * newClump);
void disposeQueryState (QueryState_t * QS);
void fprintQueryId(QueryState_t * QS, FILE * out);
void setCoverage(QueryState_t * QS, int start, int len, BOOL value);
BOOL checkCoverage(QueryState_t * QS, int start, int end, int needed, BOOL value);
void generateRandomSeed(QueryState_t * QS);

// End the cplusplus "C" extern declaration.
#ifdef __cplusplus
}
#endif


// Signatures for routines that perform the primary functions of the program.
void compressFile (AlignmentArgs_t * AAs);
void uncompressFile (AlignmentArgs_t * AAs);
void indexFile (AlignmentArgs_t * AAs);
void processQueryFile (AlignmentArgs_t * AAs);
