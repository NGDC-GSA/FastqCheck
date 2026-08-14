/*************************************************************************
    > File Name: params_parse.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2022年07月27日 星期三 13时15分22秒
 ************************************************************************/


#ifndef __PARAMS_PARSE_H__
#define __PARAMS_PARSE_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>

#define __FASTQ_CHECK_VERSION__ "2.0.0"
#define __FASTQ_CHECK_CREATE_DATE__ "2022-07-25"
#define __FASTQ_CHECK_UPDATE_DATE__ "2026-07-30"


/*! @typedef arg_t
 @abstract structure for the command line args
 @field help            [0|1] 1: print the help information
 @field max_length      the maximum length allowed for one read sequence [default:50MB]
 @field pair_check      check whether all reads have the same pair marker (1->check; 0->ignore) [default:1]
 @field n_thread        number of thread used in fastq processing [default:1]
 @field compress_ratio  the compression ratio of the fastq file [default:10]
 @field max_reads       the estimated maximum number of reads (estimated by fastq_type, in million) [required]
 @field phred_value     the phred value of the fastq file (estimated by fastq_type, 33 or 64) [required]
 @field in_file         fastq file list (eg. samples.txt)
 @field out_file        output file of the given fastq files after checking (eg. output.xml)
 @field error_file      warning or error message about the given fastq files
*/
typedef struct arg_t {
    int32_t help;
    int32_t max_length;
    int32_t pair_check;
    int32_t n_thread;
    int32_t compress_ratio;
    uint64_t max_reads;
    int32_t phred_value;
    char *in_file;
    char *out_file;
    char *error_file;
} arg_t;


/*! @function
  * @abstract     command line parameters parsing.
  @param  argc    the number of args.
  @param  argv    the parameters.
  @return         args structure.
 */
arg_t *args_parse(int argc, char **argv);


#endif
