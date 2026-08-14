/*************************************************************************
    > File Name: file_type.h
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2022年11月03日 星期四 09时36分35秒
 ************************************************************************/

#ifndef __FILE_TYPE_H
#define __FILE_TYPE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include <bzlib.h>


/*! @typedef type_t
  @abstract the file type object, which stores the magic code of a file format
  @field  index          the offset of the magic code in the file header (0 for gz/bz2/zip/rar, 257 for tar)
  @field  len            the length of the magic code
  @field  code           the magic code of the file format
  @field  format         the name of the file format (eg. gz, bz2, zip, rar, tar)
 */
typedef struct type_t {
    uint32_t index;
    uint32_t len;
    uint8_t code[8];
    char format[8];
} type_t;


int file_type_check(char *file_list, FILE *err_fp);


#endif //__FILE_TYPE_H
