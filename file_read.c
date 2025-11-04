/*************************************************************************
    > File Name: file_read.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2022年07月26日 星期二 13时15分22秒
 ************************************************************************/

#include "file_read.h"


/* copy a string and allocate necessary memory */
static char *file_name_copy(const char *str)
{
    char *dest_str;
    int l_str = strlen(str) + 1, idx;

    dest_str = (char *)malloc(l_str * sizeof(char));
    if (!dest_str) {
        fprintf(stderr, "[SysError:file_name_copy:001] failed to malloc memory when copy file name %s!\n", str);
        exit(-1);
    }

    for (idx=0; idx < l_str; ++idx) {
        if (str[idx]=='\r' || str[idx]=='\n' || str[idx]=='\0') break;
        dest_str[idx] = str[idx];
    }
    dest_str[idx] = '\0';
    return dest_str;
}


FileObject *read_file_list(char *file_list)
{
    char buf[PATH_MAX];
    FILE *file_fp = fopen(file_list, "r");

    if (!file_fp) {
        char *err_fn = get_path_basename(file_list);
        fprintf(stderr, "[SysError:read_file_list:002] failed to open the file list of %s\n", err_fn);
        exit(-1);
    }

    FileObject *file_obj = calloc(1, sizeof(FileObject));
    char *fn;
    
    while (fgets(buf, PATH_MAX, file_fp)) {
        if (buf[0]=='\r' || buf[0]=='\n') continue;  /* skip blank line */

        if (file_obj->n == file_obj->n_max) {
            file_obj->n_max = file_obj->n_max ? file_obj->n_max<<1 : 4;
            file_obj->file_name = (char **)realloc(file_obj->file_name, sizeof(char*) * file_obj->n_max);
            file_obj->gz_hd = (GzStream **)realloc(file_obj->gz_hd, sizeof(GzStream*) * file_obj->n_max);
        }
        fn = file_name_copy(buf);
        file_obj->file_name[file_obj->n] = fn;
        file_obj->gz_hd[file_obj->n++] = gz_stream_open(fn, "r");
    } fclose(file_fp);

    return file_obj;
}


/* check whether the str is endwiths sub */
/* eg. get_str_ends('/home/xlzh/test.fq.gz', '.gz') -> true */
static int get_str_ends(char *str, char *sub)
{
    int i=strlen(str)-1;
    int j=strlen(sub)-1;

    if(i<j) return 0;
    for(; i>=0 && j>=0; i--,j--){
        if(str[i]!=sub[j]) return 0;
    }
    return 1;
}


/* get the basename of the given file path */
char *get_path_basename(char *file_path)
{
    #ifdef _WIN32
      #define _PDELIM_ 92 // '\'
    #else
      #define _PDELIM_ 47 // '/'
    #endif

    char *start;

    start = strrchr(file_path, _PDELIM_);
    if (start != NULL) return start + 1;

    return file_path;
}


GzStream *gz_stream_open(char *file, char *mode)
{
    if(strcmp(mode, "w") != 0 && strcmp(mode, "r") != 0) {
        fprintf(stderr, "[SysError:gz_stream_open:003] operate mode(%s) error, it should be \"w\" or \"r\".\n", mode);
        exit(0);
    }

    GzStream *gz = calloc(1, sizeof(GzStream));
    char *err_fn = get_path_basename(file);
    if(mode[0]=='r') {
        if(get_str_ends(file, ".bz2")){
            unsigned char buf[4];
            FILE *f=fopen(file, "rb");
            if(f==NULL) {
                fprintf(stderr, "[SysError:gz_stream_open:004] can not open bz2 file of (%s)!\n", err_fn);
                exit(0);
            }
            fread(buf, 1, 4, f);  /* read the magic number of the bzip2 file */
            fseek(f, 0, SEEK_SET);
            if (buf[0]==0x42 && buf[1]==0x5a && buf[2]==0x68 || buf[0]=='@') {
                gz->bz2_fp=BZ2_bzReadOpen(&(gz->bzerror), f, 0, 0, NULL, 0);
                if(gz->bzerror!=BZ_OK){
                    BZ2_bzReadClose(&(gz->bzerror), gz->bz2_fp);
                    fprintf(stderr, "[SysError:gz_stream_open:005] the file of (%s) does not seem to be a .bz2 file!\n", err_fn);
                    exit(0);
                }
            } else {
                fprintf(stderr, "[SysError:gz_stream_open:006] the file of (%s) is neither a real bzip2(.bz2) file nor a stander fastq file!\n", err_fn);
                exit(0);
            }
        }else{ /* open a .gz or a normal txt file for reading */
            unsigned char buf[4];
            FILE *f = fopen(file, "rb");
            if (f==NULL) {
                fprintf(stderr, "[SysError:gz_stream_open:007] failed to open file of (%s)!\n", err_fn);
                exit(0);
            }
            fread(buf, 1, 4, f);  /* read the magic number of the gzip file */
            if (buf[0]==0x1f && buf[1]==0x8b || buf[0]=='@') {
                gz->gz_fp=gzopen(file, "r");
                if(gz->gz_fp <= 0){
                    fprintf(stderr, "[SysError:gz_stream_open:008] failed to open gzip file of (%s)!\n", err_fn);
                    exit(0);
                }
            } else {
                fprintf(stderr, "[SysError:gz_stream_open:009] the file of (%s) is neither a real gzip(.gz) file nor a stander fastq file!\n", err_fn);
                exit(0);
            }
        }
        gz->is_write=0;
    }
    else {  /* open a gz or bz2 file for writing */
        if(get_str_ends(file, ".gz")) {
            gz->gz_fp=gzopen(file, "w");
            if(gz->gz_fp<=0){fprintf(stderr, "[SysError:gz_stream_open:010] failed to create file of (%s) in .gz format!\n", err_fn);exit(0);}
        }else if(get_str_ends(file, ".bz2")){
            FILE *f=fopen(file, "w");
            if(f==NULL){
                fprintf(stderr, "[SysError:gz_stream_open:011] failed to create file of (%s)!\n", err_fn);
                exit(0);
            }
            gz->bz2_fp=BZ2_bzWriteOpen(&(gz->bzerror), f, 9, 0, 30);
            if(gz->bzerror != BZ_OK){
                BZ2_bzWriteClose(&(gz->bzerror), gz->bz2_fp, 0, NULL, NULL);
                fprintf(stderr, "[SysError:gz_stream_open:012] failed to create file of (%s) in .bz2 format!\n", err_fn);
                exit(0);
            }
        }else{
            gz->out_fp=fopen(file, "w");
            if(gz->out_fp<=0){fprintf(stderr, "[SysError:gz_stream_open:013] failed to create a normal file of (%s)!\n", err_fn);
            exit(0);}
        }
        gz->is_write=1;
    }
    gz->buf=calloc(GZ_BUFF_SIZE, sizeof(char));
    gzbuffer(gz->gz_fp, GZ_BUFF_SIZE<<5);  /* increase the IO throughput (from 8KB to 128MB) */

    return gz;
}


int gz_read_util(GzStream *gz, char delimiter, kstring_t *ks_str, int max_length)
{
    int len = 0;
    char c;

    if(gz->is_eof && gz->begin>=gz->end) return 0;  /* end of the file */

    do {
        if (gz->begin >= gz->end){ /* gz->buf is full or the first time to read */
            gz->begin = 0;
            if (gz->gz_fp) {
                gz->end = gzread(gz->gz_fp, gz->buf, GZ_BUFF_SIZE);
                gzerror(gz->gz_fp, &(gz->bzerror));
                if (gz->bzerror < 0) return -1; /* truncated file detected */
            }
            else if (gz->bz2_fp) {
                gz->end = BZ2_bzRead(&(gz->bzerror), gz->bz2_fp, gz->buf, GZ_BUFF_SIZE);
                if (gz->bzerror!=BZ_OK && gz->bzerror!=BZ_STREAM_END) return -1; /* truncated file detected */
            }
            if (gz->end < GZ_BUFF_SIZE) gz->is_eof = 1;
        }
        while (gz->begin < gz->end) { /* copy the str to the user-provided ks_str */
            c = gz->buf[gz->begin++];
            if (len >= ks_str->m) {
                ks_str->m = ks_str->m ? ks_str->m<<1 : 512;
                if (ks_str->m > max_length) {
                    fprintf(stderr, "[SysError:gz_read_util:014] read length can not be longer than %d!\n", max_length);
                    fprintf(stderr, "[*] the main reason is that line breaks('\\n') can not be detected in the READ!\n");
                    return -2;
                }
                ks_str->s = (char *)realloc(ks_str->s, ks_str->m * sizeof(char));
                if (!ks_str->s) {
                    fprintf(stderr, "[SysError:gz_read_util:015] failed to reallocated memory!\n");
                    exit(-1);
                }
            }
            if (c == delimiter) {
                ks_str->s[len] = '\0'; ks_str->l = len;
                return 1;
            }
            ks_str->s[len++] = c;
        }
    } while (!gz->is_eof);

    return 0;  /* end of the file (EOF) */
}


void gz_stream_destory(GzStream *gz)
{
    if (gz->is_write){  /* write gz or gz2 file */
        if(gz->gz_fp){
            gzwrite(gz->gz_fp, gz->buf, gz->begin);
            gzclose(gz->gz_fp);
        }else if(gz->bz2_fp){
            BZ2_bzWrite(&(gz->bzerror), gz->bz2_fp, gz->buf, gz->begin);
            BZ2_bzWriteClose(&(gz->bzerror), gz->bz2_fp, 0, NULL, NULL);
        }else{
            fwrite(gz->buf, 1, gz->begin, gz->out_fp);
            fclose(gz->out_fp);
        }
    }
    else {  /* read gz or bz2 file */
        if(gz->bz2_fp) BZ2_bzReadClose(&(gz->bzerror), gz->bz2_fp);
        else gzclose(gz->gz_fp);
    }
    if(gz->buf) free(gz->buf);
    free(gz);
}
