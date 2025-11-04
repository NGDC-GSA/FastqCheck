/*************************************************************************
    > File Name: file_read.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2022年07月26日 星期二 13时15分22秒
 ************************************************************************/


#ifndef __FILE_READ_H__
#define __FILE_READ_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <bzlib.h>

/* number of bytes readed from *.gz|*.bz2 file */
#define GZ_BUFF_SIZE (1048576<<2)
#define PATH_MAX 1024


#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif


/*! @typedef GzStream
 @abstract structure for gz/bz2 handle
 @field gz_fp        file handle of *.gz
 @field bz2_fp       file handle of *.bgz2
 @field out_fp       file handle of output file
 @field buf          used to store decompressed data, whose size is 'gz_stream_buff_size'
 @field begin, end   begin and end index in the buf
 @field is_write     is_write:1 -> (for file write)'w', otherwise -> 'r' (for file read)
 @ is_eof            is_eof:1 -> the end of the file
*/
typedef struct {
    gzFile gz_fp;
    BZFILE *bz2_fp;
    FILE *out_fp;
    int bzerror;
    char *buf;
    int begin;
    int end;
    int is_write;
    int is_eof;
} GzStream;


/*! @typedef FileObject
  @abstract the file object constructed by the files provided by user
  @field  n           the nubmer of files
  @field  n_max       the maximum memory allocated
  @field  file_name   the name of each file
  @field  gz_hd       file GzStream handle for each gz/bz2 file
 */
typedef struct {
    int n, n_max;
    char **file_name;
    GzStream **gz_hd;
} FileObject;


/*! @function: open a gz or bz2 file
  @param  file        input file name
  @param  mode        operation can only be 'r'(read) and 'w'(write)
  @return             GzStream object
 */
GzStream *gz_stream_open(char *file, char *mode);


/*! @function: read one line(until delimiter) from the compressed file
  @param  gz          GzStream object
  @param  delimiter   delimiter used to decide where to break, e.g. '\n'
  @param  ks_str      kstring_t type of string (must be NULL for the first time)
  @param  max_length  the maximum length allowed for one read sequence (default: 50MB)
  @return             operation status: 0->EOF; 1->OK; -1:truncate
 */
int gz_read_util(GzStream *gz, char delimiter, kstring_t *ks_str, int max_length);


/*! @function: destroy the GzStream object
  @param  gz          GzStream object
  @return             
 */
void gz_stream_destory(GzStream *gz);


/*! @function: read fastq file list
  @param  file_list   input file list, which contains the full path of each fastq file
  @return             pointer to FileObject
 */
FileObject *read_file_list(char *file_list);


/*! @function: get the basename of the given file path
  @param  file_path   file path
  @return             pointer to the start of the basename in the given file path
 */
char *get_path_basename(char *file_path);


#endif