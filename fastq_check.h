/*************************************************************************
    > File Name: fastq_check.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2022年07月26日 星期二 13时10分12秒
 ************************************************************************/


#ifndef __FASTQ_CHECK_H__
#define __FASTQ_CHECK_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bloom_filter.h"
#include "file_read.h"
#include "kqueue.h"


#define CACHE_SIZE 10000  /* number of cached read objects for each file */
#define BLOOM_ERROR 0.000000001  /* probability of false positive for bloomfilter */
#define READ_NAME_DUPLICATE_MAX 1000  /* maximum number of duplicate read name allowed */
#define MAX_INVALID_BASE_RATIO 0.01  /* maximum ratio of Non-ATCGN in the sequence */


#ifndef KSTRING_T
#define KSTRING_T kstring_t
/*! @typedef kstring_t
  @abstract the klib string object, which could be used to store any string buffer
  @field  l                the length of the string (without the null terminator)
  @field  m                the max memory allocated for the string buffer
  @field  s                the pointer to the string buffer
 */
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif


/*! @typedef read_t
  @abstract the read object
  @field  name           the name of the read, which will be truncated at the pair_marker (the marker will be detected and returned) when 'pair_check' is given,
                         eg. "@ST-E00126:HWM7:3173:1784 2:N:AAGAC" -> name: "@ST-E00126:HWM7:3173:1784" (marker: '2')
  @field  seq            the sequence of the read
  @field  comment        the comment of the read, usually is a character of '+'
  @field  qual           the quality of the read
 */
typedef struct read_t {
    kstring_t name;
    kstring_t seq;
    kstring_t comment;
    kstring_t qual;
} read_t;


/*! @typedef cache_t
  @abstract the fastq cache
  @field  n        the number of read objects
  @field  n_max    the max size of the cache
  @field  reads    the pointer to the reads object array
 */
typedef struct cache_t {
    size_t n, n_max;
    read_t *reads;
} cache_t;


/*! @typedef length_t
  @abstract prepare for the length distribution table
  @field  n_max     the max memory allocated for lengths array
  @field  lengths   the pointer to the lengths array
 */
typedef struct length_t {
    size_t n_max;
    uint64_t *lengths;
} length_t;


/*! @typedef result_t
  @abstract the result object, which contains all needed summary count
  @field  n_reads           the number of reads
  @field  n_bytes           the real bytes used by the fastq file
  @field  base_table        table of base distribution for all 256 characters ([ATCGNatcgn] are the valid bases)
  @field  qual_table        table of base quality distribution
  @field  length_obj        the length distribution table of the reads (init: 512 and could be resized)
 */
typedef struct result_t {
    uint64_t n_reads;
    uint64_t n_bytes;
    uint64_t base_table[256];
    uint64_t qual_table[256];
    length_t length_obj;
} result_t;


/*! @typedef message_t
  @abstract the message (error/warning/output result) of fastq checking
  @field  error_fp          pointer to the error message file
  @field  out_fp            pointer to the output file
  @field  n_duplicate       potential duplicate read name checked
  @field  name_check        check whether the read name is not in the same order (0 -> checked)
  @field  length_check      whether the read length has not been checked yet (0 -> checked)
  @field  pair_check        whether the read pair has not been checked yet (0 -> checked)
  @field  err_lock          the mutex used to protect the error file (only used in multithread mode)
 */
typedef struct message_t {
    FILE *err_fp;
    FILE *out_fp;
    int n_duplicate;
    int name_check;
    int length_check;
    int pair_check;
    pthread_mutex_t err_lock;
} message_t;


/*! @typedef thread_task_t
  @abstract the task object of fastq batch checking (each batch is processed by one thread)
  @field  status           the stage mask of the task (1: batch cached -> 11: content checked -> 111: duplicate checked)
  @field  n_file           the number of input files (1: single-end, 2: pair-end, 4: single-cell)
  @field  n_hash           the number of hash values calculated for the current batch (used by bloomfilter)
  @field  phred_value      the phred value of the fastq files (33 or 64)
  @field  max_length       the maximum length allowed for one read sequence
  @field  pair_check       check whether all reads have the same pair marker (1->check; 0->ignore)
  @field  total_reads      the total number of reads cached from the input files (accumulated over batches)
  @field  cache            the reads cache array of all files
  @field  result           the check result array of all files
  @field  hashes           the hash values array of the current batch (size: CACHE_SIZE x 2)
  @field  files_hd         the file handles of all input files
  @field  msg_obj          the message object (error file, output file and check flags)
  @field  status_lock      the mutex used to protect the status mask
  @field  batch_ready      the condition variable, signaled when all files of a batch are cached
 */
typedef struct thread_task_t {
    int32_t status;
    int32_t n_file;
    int32_t n_hash;
    int32_t phred_value;
    int32_t max_length;
    int32_t pair_check;
    uint64_t total_reads;
    cache_t *cache;
    result_t *result;
    uint64_t *hashes;
    GzStream **files_hd;
    message_t *msg_obj;
    pthread_mutex_t status_lock;
    pthread_cond_t batch_ready;
} thread_task_t;


/*! @typedef reader_job_t
  @abstract the reader job object, one job for one file of the batch
  @field  task             the pointer to the thread task object being processed
  @field  file_idx         the index of the file in the task
  @field  ret_value        the return value of the fastq reading (1: normal, 0: EOF, negative: error)
  @field  read_offset      the offset of the first read in the whole file
 */
typedef struct reader_job_t {
    thread_task_t *task;
    int32_t file_idx;
    int32_t ret_value;
    uint64_t read_offset;
} reader_job_t;


/*! @typedef thread_args_t
  @abstract the arguments object shared by the reader thread and the duplicate-check thread
  @field  task_index       the index of the current task (round-robin over the task array)
  @field  task_size        the total number of tasks (same as the thread number -t)
  @field  msg_obj          the message object (error file, output file and check flags)
  @field  task_obj         the array of the thread task objects
  @field  bloom_filter     the bloomfilter used to check duplicate read names
  @field  file_obj         the input file object (file names and file handles)
  @field  task_queue       the kqueue storing the tasks from the reader thread to the duplicate-check thread
 */
typedef struct thread_args_t {
    int32_t task_index;
    int32_t task_size;
    message_t *msg_obj;
    thread_task_t *task_obj;
    BaseBloomFilter *bloom_filter;
    FileObject *file_obj;
    kqueue_t *task_queue;
} thread_args_t;


/*! @typedef stat_file_t
  @abstract the statistical count for each file
  @field  length_avg      the average length of the file
  @field  length_stdev    the standard deviation of the file
  @field  n_bases         the number of bases of the file
 */
typedef struct stat_file_t {
    double length_avg;
    double length_stdev;
    uint64_t n_bases;
} stat_file_t;


/*! @typedef statistic_t
  @abstract the statistic count with the result_obj, which will be used in output .xml file
  @field  phred             the phred value of the fastq files
  @field  n_reads           the number of reads (all files have same reads number)
  @field  n_total_bases     the total number of bases for all files
  @field  n_total_bytes     the total number of bytes for all files
  @field  n_bases_num       the number of 'ACGTN' (in order)
  @field  qual_table        the quality distribution for all files
  @field  n_stat            the number of files (the size of the stat_obj array)
  @field  stat_obj          one stat object for one file
  @field  gc_content        the GC content of all the files
  @field  length_status     only 'fixed' and 'variable' is available for this variable
  @field  access_path       the access path of the file name
 */
typedef struct statistic_t {
    uint64_t phred;
    uint64_t n_reads;
    uint64_t n_total_bases;
    uint64_t n_total_bytes;
    uint64_t n_bases_num[8];
    uint64_t qual_table[256];
    uint64_t n_stat;
    stat_file_t *stat_obj;
    double gc_content;
    char length_status[16];
    char *access_path;
} statistic_t;


/* safely malloc memory for type */
#ifndef err_malloc
#define err_malloc(_p, _n, _type) do {                                  \
    _type *tem = (_type *)malloc((_n) * sizeof(_type));                 \
    if (!tem) {                                                         \
        fprintf(stderr, "\n[SysError:err_malloc]: failed to malloc memory!\n"); exit(-1); \
    } (_p) = tem;                                                       \
} while(0)
#endif


/* safely calloc memory for type */
#ifndef err_calloc
#define err_calloc(_p, _n, _type) do {                                  \
    _type *tem = (_type *)calloc((_n), sizeof(_type));                  \
    if (!tem) {                                                         \
        fprintf(stderr, "\n[SysError:err_calloc]: failed to calloc memory!\n"); exit(-1); \
    } (_p) = tem;                                                       \
} while(0)
#endif


/* safely realloc memory for type */
#ifndef err_realloc
#define err_realloc(_p, _n, _type) do {                                 \
    _type *tem = (_type *)realloc((_p), (_n) * sizeof(_type));          \
    if (!tem) {                                                         \
        fprintf(stderr, "\n[SysError:err_realloc]: failed to realloc memory!\n"); exit(-1); \
    } (_p) = tem;                                                       \
} while(0)
#endif


#endif
