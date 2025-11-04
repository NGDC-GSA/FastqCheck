/*************************************************************************
    > File Name: file_type.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2022年11月03日 星期四 09时37分12秒
 ************************************************************************/

#include "file_type.h"


static const int N_TYPE = 5;

static const type_t TYPE_ARRAY[5] = {
        {0, 2, {0x1f, 0x8b}, "gz"},
        {0, 3, {0x42, 0x5a, 0x68}, "bz2"},
        {0, 2, {0x50, 0x4b}, "zip"},
        {0, 4, {0x52, 0x61, 0x72, 0x21}, "rar"},
        {257, 5, {0x75, 0x73, 0x74, 0x61, 0x72}, "tar"}
};


static int check_file_suffix(char *file_name, char *suffix)
{
    int idx1, idx2;

    idx1 = (int)strlen(file_name) - 1;
    idx2 = (int)strlen(suffix) - 1;

    if (idx1 < idx2) return 0;
    for (; idx1 >= 0; idx1--, idx2--) {
        if (file_name[idx1] != suffix[idx2]) return 0;
    }
    return 1;
}


static void line_break_trim(char *file_name)
{
    int l_str = (int)strlen(file_name) + 1;
    int idx;

    for (idx=0; idx < l_str; ++idx) {
        if (file_name[idx]=='\r' || file_name[idx]=='\n') break;
    }
    file_name[idx] = '\0';
}


static char *get_file_basename(char *file_path)
{
    char *start;

    start = strrchr(file_path, '/');
    if (start != NULL) return start+1;

    return file_path;
}


static const type_t *get_file_format(const unsigned char *cache, int n_byte)
{
    const type_t *type_obj;
    int flag;

    for (int i=0; i < N_TYPE; i++) {
        flag = 1;
        type_obj = &TYPE_ARRAY[i];

        if (n_byte < type_obj->index+type_obj->len)  /* it can not be 'tar' file */
            continue;

        for (int j=0; j < type_obj->len; j++) {
            if (cache[type_obj->index+j] != type_obj->code[j]) {
                flag = 0; break;
            }
        }
        if (flag == 1) return type_obj;
    }
    return NULL;
}


int file_type_check(char *file_list, FILE *err_fp)
{
    char buf[1024];
    int status = 0;

    FILE *file_fp = fopen(file_list, "r");
    if (file_fp == NULL) {
        char *err_fn = get_file_basename(file_list);
        fprintf(stderr, "[SysError:file_type_check:002] failed to open the file list of %s\n", err_fn);
        exit(-1);
    }

    while (fgets(buf, 1024, file_fp)) {
        if (buf[0] == '\r' || buf[0] == '\n') continue;  /* skip blank line */
        line_break_trim(buf);  /* trimming the '\r' and '\n' in the end of the filename */

        /* read 512 bytes from the file and check whether the suffix is accurate */
        unsigned char cache[512];
        char *err_fn = get_file_basename(buf);
        FILE *fp = fopen(buf, "rb");

        if (fp == NULL) {
            fprintf(stderr, "[SysError:file_type_check:007] failed to open file of %s\n", err_fn);
            exit(-1);
        }
        int n_byte = (int)fread(cache, 1, 512, fp);
        const type_t *type_obj = get_file_format(cache, n_byte);

        if (type_obj == NULL) {  /* not gz, bz2, zip, rar or tar file */
            status = -1;
            fprintf(err_fp, "[FileError:file_type_check:102] unsupported file format (.XXX), only file type of (.gz) and (.bz2) is available (Err: %s)!\n", err_fn);
            continue;
        }

        if (strcmp(type_obj->format, "gz") == 0) {  /* is .gz file */
            gzFile gz_fp = gzopen(buf, "r");
            if (gz_fp <= 0) {
                fprintf(stderr, "[SysError:file_type_check:005] failed to open gzip file of %s!\n", err_fn);
                exit(-1);
            }
            n_byte = gzread(gz_fp, cache, 512);
            if (n_byte < 0) {
                fprintf(err_fp, "[FileError:file_type_check:101] unexpected end of fastq file %s (truncated file) is detected!\n", err_fn); 
                exit(-1);
            }
            if (cache[0] != '@') {  /* the suffix of the filename has been changed */
                status = -1;
                type_obj = get_file_format(cache, n_byte);
                if (type_obj == NULL)  /* unsupported file format (eg. A.exe.gz -> A.gz) */
                    fprintf(err_fp, "[FileError:file_type_check:103] unsupported gzip file format of (.XXX.gz), please do not rename the original suffix of the filename (Err: %s)!\n", err_fn);
                else  /* known file format (eg. A.tar.gz->A.gz; A.rar.gz->A.gz; A.zip.gz->A.gz)*/
                    fprintf(err_fp, "[FileError:file_type_check:104] the file format of %s may be (.%s.gz), please do not rename the original suffix of the filename!\n", err_fn, type_obj->format);
            }
            gzclose(gz_fp); fclose(fp);
        }
        else if (strcmp(type_obj->format, "bz2") == 0) {  /* is .bz2 file */
            int bz_error=0;
            fseek(fp, 0, SEEK_SET);
            BZFILE *bz2_fp = BZ2_bzReadOpen(&bz_error, fp, 0, 0, NULL, 0);
            if (bz_error!=BZ_OK && bz_error!=BZ_STREAM_END) {
                fprintf(stderr, "[SysError:file_type_check:004] can not open bz2 file of %s!\n", err_fn);
                exit(-1);
            }
            n_byte = BZ2_bzRead(&bz_error, bz2_fp, cache, 512);
            if (bz_error!=BZ_OK && bz_error!=BZ_STREAM_END) {
                fprintf(err_fp, "[FileError:file_type_check:101] unexpected end of fastq file %s (truncated file) is detected!\n", err_fn); 
                exit(-1);
            }
            if (cache[0] != '@') {  /* the suffix of the filename has been changed */
                status = -1;
                type_obj = get_file_format(cache, n_byte);
                if (type_obj == NULL)  /* unsupported file format (eg. A.exe.gz -> A.gz) */
                    fprintf(err_fp, "[FileError:file_type_check:105] unsupported bzip2 file format of (.XXX.bz2), please do not rename the original suffix of the filename (Err: %s)!\n", err_fn);
                else  /* known file format (eg. A.tar.gz->A.gz; A.rar.gz->A.gz; A.zip.gz->A.gz)*/
                    fprintf(err_fp, "[FileError:file_type_check:106] the file format of %s may be (.%s.bz2), please do not rename the original suffix of the filename!\n", err_fn, type_obj->format);
            }
            BZ2_bzReadClose(&bz_error, bz2_fp);
        }
        else {  /* may be zip, rar, tar */
            status = -1;
            fprintf(err_fp, "[FileError:file_type_check:107] the file format of %s may be (.%s), please do not rename the original suffix of the filename!\n", err_fn, type_obj->format);
        }
    } fclose(file_fp);

    return status;
}