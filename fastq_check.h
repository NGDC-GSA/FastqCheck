/*************************************************************************
    > File Name: fastq_check.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2022年07月26日 星期二 13时10分12秒
 ************************************************************************/


#ifndef __FASTQ_CHECK_H__
#define __FASTQ_CHECK_H__


#include <stdlib.h>
#include <string.h>


#define CACHE_SIZE 10000  /* number of cached read object for each file */
#define BLOOM_ERROR 0.000000001  /* probability of false positive for bloomfilter */
#define READ_NAME_DUPLICATE_MAX 1000  /* maximum number of duplicate read name allowed */
#define MAX_INVALID_BASE_RATIO 0.01  /* maximum ratio of Non-ATCGN in the sequence */


#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif


/*! @typedef read_t
  @abstract the read object
  @field  name           the name of read('@' will be repleaced with pair_marker if 'pair_check' is given). 
                         eg. @ST-E00126:HWM7:3173:1784 2:N:AAGAC -> 2ST-E00126:HWM7:3173:1784
  @field  seq            the sequence of the read
  @field  comment        the comment of the read, usually is a character of '+'
  @field  qual           the quality of the read
 */
typedef struct __read_t {
    kstring_t name;
    kstring_t seq;
    kstring_t comment;
    kstring_t qual;
} read_t;


/*! @typedef cache_t
  @abstract the fastq cache
  @field  n        the number of read object
  @field  n_max    the max size of the cache
  @field  reads    the pointer to the reads object array
 */
typedef struct __cache_t {
    size_t n, n_max;
    read_t *reads;
} cache_t;


/*! @typedef length_t
  @abstract prepare for the length distribution table
  @field  n_max     the max memory allocated for lengths array
  @field  lengths   the pointer to the lengths array
 */
typedef struct __length_t {
    size_t n_max;
    uint64_t *lengths;
} length_t;


/*! @typedef result_t
  @abstract the result object, which contains all needed summary count
  @field  phred             the phred value of the sequence
  @field  n_reads           the number of reads
  @field  n_bytes           the real bytes used by the fastq file
  @field  base_table        only focus on [ATCGNatcgn] and other characters also need to calculate
  @field  qual_table        table of base quality distribution
  @field  len_table         table of sequence length distribution
  @field  len_table_size    the size of the len_table (init: 512 and could be resized)
 */
typedef struct __result_t {
    uint64_t phred;
    uint64_t n_reads;
    uint64_t n_bytes;
    uint64_t base_table[256];
    uint64_t qual_table[256];
    length_t length_obj;
} result_t;


/*! @typedef stat_file_t
  @abstract the statistical count for each file
  @field  length_avg      the average length of the file
  @field  length_stdev    the standard deviation of the file
  @field  n_bases         the bases number of the file
 */
typedef struct __statistic_core_t {
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
  @field  n_bases_num       the number of 'ACGTN' (in oreder)
  @field  qual_table        the quality distribution for all files
  @field  n_stat            infact, it is the number of files
  @field  stat_obj          one stat object for one file
  @field  gc_content        the GC content of all the files
  @field  length_status     only 'fixed' and 'variable' is available for this variables
  @field  access_path       the access path of the file name
 */
typedef struct __statistic_t {
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


/*! @typedef message_t
  @abstract the message (error/warning/output result) of fastq checking
  @field  error_fp          pointer to the error message file
  @field  out_fp            pointer to the output file
  @field  cache_status      check whether the process of caching finised (!=1 means finished)
  @field  n_duplicate       potential duplicate read name checked
  @field  name_check        check whether the read name is not in same order (0 -> checked)
  @field  length_check      whether the read length has not been checked yet (0 -> checked)
  @field  pair_check        whether the read pair has not been checked yet (0 -> checked)
 */
typedef struct __message_t {
    FILE *err_fp;
    FILE *out_fp;
    int cache_status;
    int n_duplicate;
    int name_check;
    int length_check;
    int pair_check;
} message_t;


#endif
