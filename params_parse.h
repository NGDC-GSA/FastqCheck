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
#include <getopt.h>

#define __FASTQ_CHECK_VERSION__ "1.1.21"
#define __FASTQ_CHECK_CREATE_DATE__ "2022-07-25"
#define __FASTQ_CHECK_UPDATE_DATE__ "2024-04-26"


/*! @typedef arg_t
 @abstract structure for the command line args
 @field help            [0|1] 1: print the help information
 @field max_length      the maximum length allowed for one read sequence [default:50MB]
 @field pair_check      check whether all reads have the same pair marker (1->check; 0->ignore) [default:1]
 @field n_thread        number of thread used in fastq processing [default:1]
 @field ratio           the compression ratio of the fastq file [default:10]
 @field in_file         fastq file list (eg. samples.txt)
 @field out_file        output file of the given fastq files after checking (eg. output.xml)
 @field error_file      warning or error message about the given fastq files
*/
typedef struct arg_t {
    int help;
    int max_length;
    int pair_check;
    int n_thread;
    int ratio;
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
