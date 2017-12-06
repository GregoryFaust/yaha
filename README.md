# *yaha*

**Written by:** Greg Faust (gf4ea@virginia.edu)  
[Ira Hall Lab, University of Virginia](http://faculty.virginia.edu/irahall/)

**Please cite:**  
[Faust G.G. and Hall I.M., "*YAHA*: fast and flexible long-read alignment with optimal breakpoint detection,"
*Bioinformatics* Oct. 2012; **28**(19): 2417-2424.](http://bioinformatics.oxfordjournals.org/content/28/19/2417)

There is an extensive [User Guide](https://www.dropbox.com/s/7j758vpbaskcq20/YAHA_User_Guide.0.1.83.pdf?dl=0)
that supplements the command information below.  
Click the preceeding link or download the file from this repository.

---

**Current version:** 0.1.83

Current support for Linux only.

## Summary
*yaha* is an open source, flexible, sensitive and accurate DNA aligner designed for single-end reads.
It supports three major modes of operation:

1. The default “Optimal Query Coverage” (**-OQC**) mode reports the best set of alignments that cover the length of each 
query. 
2. Using “Filter By Similarity” (**-FBS**), along with the best set of alignments,
*yaha* will also output alignments that are highly similar to an alignment in the best set. 
3. Finally, *yaha* can output all the alignments found for each query.  

The **-OQC** and **-FBS** modes are specifically tuned to form split read mappings that can be used to accurately 
identify structural variation events (deletions, duplications, insertions or inversions) 
between the subject query and the reference genome.

## Installation
*yaha* can be downloaded from the **_releases_** tab or manually downloaded via *git clone*.  For example:
~~~~~~~~~~~~~~~~~~
git clone git://github.com/GregoryFaust/yaha.git
cd yaha
make
cp bin/yaha /usr/local/bin/
~~~~~~~~~~~~~~~~~~

## Usage

**COMMON USAGE SCENARIOS:** 

To create an index for a reference genome:
```
yaha -g <genomeFilename> -H <maxHits> -L <seedLength> -S <Skipdistance>
```

To align sequencing data:
```
yaha -x <yahaIndexFile> -q <queryFile> -osh <outputFile> [Additional Options...]
```

---
**OPTIONS:**
Default values enclosed in square brackets []
```
Input/Output Options:
-g    FILE input genome file to use during index creation (FASTA or nib2)
-q    FILE input file of sequence reads to align (FASTA or FASTQ) [STDIN]
-osh  FILE output file for alignment output in SAM format with hard clipping(default) [STDOUT]
-oss  FILE output file for alignment output in SAM format with soft clipping [STDOUT]
-x    FILE reference index file to use during alignment
NOTE: At most one of -osh or -oss should be specified.

Index Creation Options:
-H    INT  maxHits: During index creation, seeds occuring more than maxHits times will be sampled [65565]
-L    INT  seedLength: Length of seed to use.  During alignment, seed length is taken from index file [15]
-S    INT  Skipdistance: Number of bases to skip ahead before forming next seed [1]

General Alignment Options:
-BW   INT  BandWidth: band size on each side of the diagonal of banded Smith Waterman [5]
-G    INT  maxGap: maximum indel size allowed with a single alignment [50]
-H    INT  maxHits: maximum times a seed is in the reference before it is ignored as too repetitive [650]
-M    INT  minMatch: minimum number of bases in seeds to start an alignment [25]
-MD   INT  MaxDesert: maximum number of contiguous bases without a seed before alignmment is split [50] 
-P    REAL minPercent-identity: minimum matching/alignment-length for a query to be included in output [0.9]
-X    INT  Xdropoff: maximum score dropoff before terminating alignment extensions [25]
-t    INT  numThreads: number of threads used to parallel process reads [1]

Affine Gap Scoring Options:
If -AGS is off, a simple edit distance calculation is done.
If on, the remaining options are used:
-AGS  BOOL (Y|N) controls use of Affine Gap Scoring [Y].  
-GEC  INT  GapExtensionCost: cost for extending a gap (indel) [2] 
-GOC  INT  GapOpenCost: cost for starting a new gap (indel) [5]
-MS   INT  MatchScore: score added for each matching base [1] 
-RC   INT  ReplacementCost: score subtracted for each mismatched base [3]

Optimal Query Coverage Options:
If -OQC if off, all alignments meeting above criteria are output.
If -OQC is on, a set of alignments are found that optimally cover the query, using the remaining options.
-OQC  BOOL (Y|N) controls use of the Optimal Query Coverage Algorithm.
-BP   INT BreakpointPenalty: penalty for inserting a breakpoint in split-read alignment [5]
-MGDP INT MaxGenomicDistancePenalty (5)] 
-MNO  INT MinNonOverlap: minimum number of unshared bases required in each split alignment [minMatch]
NOTE: The total cost of adding a breakpoint in a split-read mapping is:
  BP*MIN(MGDP, Log10(genomic distance between reference loci))

Filter By Similarity Options:
If -FBS is on, the remaining options are used.  An alignemnt must satisfy BOTH criteria to be "similar".
-FBS  BOOL (Y|N) controls output of alignments similar to best alignment found using OQC.
-PRL  REAL PercentReciprocalLength: minimum ratio of overlapping length between similar alignemnt [0.9] 
-PSS  REAL PercentSimilarScore: minimum ratio of scores between similar alignments [0.9]
```

See the User Guide for more details on all options and their usage.
