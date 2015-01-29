/* -*- Mode: C ; indent-tabs-mode: nil ; c-file-style: "stroustrup" ; column-number-mode: t -*-

    Project: YAHA, DNA alignment tool designed to find optimal split-read mappings on single-end queries.
    Author:  Greg Faust (gf4ea@virginia.edu)

    File:    FileHelpers.c     Helper routines for both mempory mapped and stream file IO, and outputting time.
                               Key elements of the file processing use GNU extensions, and not portable.

    License Information:

    Copyright 2009-2015 Gregory G. Faust

    Licensed under the MIT license (the "License");
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at http://opensource.org/licenses/MIT

*/

#include <time.h>
#include <stdlib.h>
#include <stdio_ext.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "Math.h"

// This is the function used in ASSERT.
#ifdef DEBUG
void outputDebugString (char * str)
{
    fprintf(stderr, str);
}
#endif

// Call this for fatal run time errors.
void fatalError(const char * errorStr)
{
    fprintf(stderr, "%s\n", errorStr);
    exit(1);
}

// Find the next occurrence of the char in the str.
// Assumes the char will exist, since it doesn't know when end of string will be reached.
// The contract is that offset + retval will point to the search char.
int findChar(char * str, FSIZE offset, char search)
{
    // Make sure we move at least one char ahead.
    int count = 1;
    while (TRUE)
    {
        if (str[offset+count] == search) return count;
        count += 1;
    }
}

static inline int fputchar(FILE * out, char outchar)
{
    return putc_unlocked(outchar, out);
}

void printChars(FILE * out, char * str, int start, int len)
{
    int i;
    for (i=start; i<(start+len); i++) fputchar(out, str[i]);
}

// A convenience function for outputing time is seconds in a more useful metric.
void fprintTimeSeconds (FILE * out, double seconds, int precision)
{
    double totalseconds = seconds;
    int hours = seconds/3600.;
    if (hours > 0)
    {
        seconds -= hours * 3600;
        fprintf(out, "%dH", hours);
    }
    int minutes = seconds/60.;
    if (minutes > 0)
    {
        seconds -= minutes * 60;
        fprintf(out, "%dM", minutes);
    }
    if (hours + minutes > 0)
    {
        fprintf(out, "%.0fS", seconds);
        fprintf(out, "(%.*fS)", precision, totalseconds);
    }
    else fprintf(out, "%.*fS", precision, totalseconds);
}

void fprintTimeMicroSeconds (FILE * out, UINT64 microSeconds, int precision)
{
    fprintTimeSeconds(out, ((double)microSeconds/1000000.0), precision);
}

// Check for a file system error.  If present, output message on stderr.
void checkFSerr (SFSIZE returnCode)
{
    if (returnCode == -1)
    {
        fprintf(stderr, "File system error: %d.\n", errno);
        exit(1);
    }
}

// Check for a file system error.  If present, output better message on stderr.
void checkFSerrWithFilename (SFSIZE returnCode, char * filename)
{
    if (returnCode == -1)
    {
        if (errno == ENOENT)
            fprintf(stderr, "File '%s' does not exist.\n", filename);
        else if (errno == ENOMEM)
            if (filename != NULL)
                fprintf(stderr, "Insufficient memory available to open '%s'.\n", filename);
            else
                fprintf(stderr, "Insufficient memory available to satisfy allocation request.\n");
        else
            fprintf(stderr, "File system error %d on file '%s'.\n", errno, filename);
        exit(1);
    }
}

// Check if one file is newer than another.
BOOL fileNewerThan (char * file1Name, char * file2Name)
{
    struct stat file1stat;
    struct stat file2stat;
    int retcode;
    retcode = stat(file2Name, &file2stat);
    // If file2 does not exist, then we will say file1 is newer.
    if (retcode == -1 && errno == ENOENT) return TRUE;
    checkFSerrWithFilename(retcode, file2Name);
    retcode = stat(file1Name, &file1stat);
    checkFSerrWithFilename(retcode, file1Name);
    // Both files exist, and we have the stats.
    // Now compare the modification times.
    struct timespec file1ts = file1stat.st_mtim;
    struct timespec file2ts = file2stat.st_mtim;
    if (file1ts.tv_sec < file2ts.tv_sec) return FALSE;
    if (file1ts.tv_sec > file2ts.tv_sec) return TRUE;
    if (file1ts.tv_nsec > file2ts.tv_nsec) return TRUE;
    return FALSE;
}


/////
//     INPUT FILES
/////

// Open the file for read.
// We will use memory mapped IO.
// And also check file stats to get the file size.
FDES openForRead(char * filename, FSIZE *fileSize, char **contentPtr, BOOL preload)
{
    // First open the file.
    FDES file;
    file = open(filename, O_RDONLY);
    checkFSerrWithFilename(file, filename);
    // Now get its size.
    struct stat * sbuf = (struct stat *)malloc(sizeof(struct stat));
    int err;
    err = fstat(file, sbuf);
    checkFSerrWithFilename(err, filename);
    FSIZE size = sbuf->st_size;
#ifdef DEBUG
    fprintf(stderr, "File %s is %d bytes long.\n", filename, size);
#endif
    // Now memory map the IO.
    char * filePtr = (char *)mmap(NULL, sbuf->st_size, PROT_READ, MAP_SHARED, file, 0);
    checkFSerrWithFilename((SFSIZE)filePtr, filename);

    if (preload)
    {
        // Try reading in each page to pre-load the page tables in sequential order.
        // This seems to vastly speed virtual page table creation or something
        //     relative to paging in things as they are used.
        int count = 0;
        for (char * ptr=filePtr; ptr<filePtr+size; ptr+=sbuf->st_blksize) count += ptr[0];
        // We call a silly function here to avoid the loop being optimized out for having no effect.
        emptyFunction(count);
    }
    // If we are here, we have successfully opened and mapped the file.
    // Fill in OUT parameters, and return.
    free(sbuf);
    *fileSize = size;
    *contentPtr = filePtr;
    return file;
}

// We will close the file, and unmap the memory mapped IO.
void closeForRead(FDES ifp, FSIZE fileSize, char *contentPtr)
{
    // First unmap the IO.
    int err = munmap(contentPtr, fileSize);
    checkFSerr(err);
    err = close(ifp);
    checkFSerr(err);
}

// Open the file for seq read, such as query files.
// We will not use memory mapped IO.
// We will use char IO.
FILE * openForSeqRead(char * filename)
{
    FILE * inFile;
    if (strcmp(filename, "stdin") == 0 || strcmp(filename, "-") == 0)
    {
        fprintf(stderr, "Reading queries from stdin.\n");
        inFile = stdin;
    }
    else
        inFile = fopen(filename, "r");
    if (inFile == NULL)
    {
        fprintf(stderr, "Failure to open input file: %s.  Error number:%d\n", filename, errno);
        exit(1);
    }

    // We will tell the OS that we will do our own locking and unlocking to improve input speed.
    __fsetlocking (inFile, FSETLOCKING_BYCALLER);

    return inFile;
}

void closeForSeqRead(FILE * inFile)
{
    fclose(inFile);
}

/////
//     OUTPUT FILES
/////
//     Output files are much more complicated in the flag settings and in the closing operation.
//     All of our output files are now either stream (for alignment output) or require a preamble.
//     Therefore, we will no longer mmap as part of opening the file, nor unmap or write as part of closing.
//     Instead we will provide some helper functions for mmapping and writing.
/////

// Allocate a big block of memory as anonymous memory mapped "file".
char * blockMalloc(FSIZE size, char * string)
{
    char * retval = (char *)mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANON, 0, 0);
    checkFSerrWithFilename((SFSIZE)retval, string);
    return retval;
}

// Free big blocks of memory.  Note the size is needed.
void blockFree(char * ptr, FSIZE size)
{
    int err = munmap(ptr, size);
    checkFSerrWithFilename(err, NULL);
}

// Write blocks of memory.
// As write may not always write out all the buffer in one call, repeat until success.
void blockWrite(FDES ofp, char * ptr, FSIZE size)
{
    FSIZE done = 0;
    while (done < size)
    {
        SFSIZE curr = write(ofp, ptr+done, size-done);
        checkFSerr(curr);
        done += curr;
    }
}

// Open for write will now only open the file and return the file descriptor.
FDES openForWrite(char * filename)
{
    // First we will open the file
    // Create it if it doesn't already exist.
    // Truncate it to zero length.
    // If it is created, set the permission bits.
    // User can RWX.  Everyone can read.
    FDES file;
    file = open(filename, O_RDWR | O_CREAT | O_TRUNC, S_IRWXU | S_IRGRP | S_IROTH);
    checkFSerrWithFilename(file, filename);
    return file;
}

// Just close the file and check for errors.
void closeForWrite(FDES ofp)
{
    // Now we can close the file.
    int err = close(ofp);
    checkFSerr(err);
}

// Open for stream IO.
FILE * openForPrint(char * filename)
{
    FILE * outFile;
    if (strcmp(filename, "stdout") == 0)
        outFile = stdout;
    else
        outFile = fopen(filename, "w");
    if (outFile == NULL)
    {
        fprintf(stderr, "Failure to open output file: %s.  Error number:%d\n", filename, errno);
        exit(1);
    }

    // We will tell the OS that we will do our own locking and unlocking to improve output speed.
    // By doing it here we presumably also speed up calls to functions such as fprintf and fputs,
    //     not just the calls to putc_unlocked.
    __fsetlocking (outFile, FSETLOCKING_BYCALLER);

    return outFile;
}

void closeForPrint(FILE * outFile)
{
    fclose(outFile);
}
