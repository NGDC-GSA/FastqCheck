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
    const int l_str = (int)strlen(str) + 1;
    int idx;

    char *dest_str = (char *) malloc(l_str * sizeof(char));
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

    while (fgets(buf, PATH_MAX, file_fp)) {
        if (buf[0]=='\r' || buf[0]=='\n') continue;  /* skip blank line */

        if (file_obj->n == file_obj->n_max) {
            file_obj->n_max = file_obj->n_max ? file_obj->n_max<<1 : 4;
            file_obj->file_name = (char **)realloc(file_obj->file_name, sizeof(char*) * file_obj->n_max);
            file_obj->gz_hd = (GzStream **)realloc(file_obj->gz_hd, sizeof(GzStream*) * file_obj->n_max);
        }
        char *fn = file_name_copy(buf);
        file_obj->file_name[file_obj->n] = fn;
        file_obj->gz_hd[file_obj->n++] = gz_stream_open(fn, "r");
    } fclose(file_fp);

    return file_obj;
}


/* check whether the str is endswith sub */
/* eg. get_str_ends('/home/xlzh/test.fq.gz', '.gz') -> true */
static int get_str_ends(char *str, char *sub)
{
    int i = (int)strlen(str)-1;
    int j = (int)strlen(sub)-1;

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

    char *start = strrchr(file_path, _PDELIM_);
    if (start != NULL) return start + 1;

    return file_path;
}


/* ---------------- gz backend ---------------- */
static void *gz_handle_open(char *file, char *mode, char *err_fn)
{
    gzFile fp = gzopen(file, mode);
    if (fp == NULL) {
        if (mode[0] == 'r')
            fprintf(stderr, "[SysError:gz_stream_open:008] failed to open gzip file of (%s)!\n", err_fn);
        else
            fprintf(stderr, "[SysError:gz_stream_open:010] failed to create file of (%s) in .gz format!\n", err_fn);
        exit(0);
    }
    gzbuffer(fp, GZ_BUFF_SIZE<<5);
    return fp;
}

static int gz_read(void *fp, char *buf, int size)
{
    int zerr;
    int n = gzread((gzFile)fp, buf, size);
    gzerror((gzFile)fp, &zerr);

    /* gzread may return the last bytes with the error already flagged,
     * so the error state must be checked after every read */
    if (n < 0 || zerr < 0) return -1;  /* truncated file detected */
    return n;
}

static int gz_write(void *fp, const char *buf, int size)
{
    return gzwrite((gzFile)fp, buf, size);  /* the number of bytes written; <0: error */
}

static int gz_close(void *fp)
{
    return gzclose((gzFile)fp);  /* 0: normal; non-zero: error */
}

/* ---------------- bz2 backend ---------------- */
/* wrapper of BZFILE, recording the open mode so that a single
 * bz2_close callback can call the right close function */
typedef struct bz2_stream {
    BZFILE *fp;
    int is_write;  /* 1: opened for writing; 0: opened for reading */
} bz2_stream;

static void *bz2_handle_open(FILE *raw_fp, char *mode, char *err_fn)
{
    int bzerror;
    bz2_stream *bz = calloc(1, sizeof(bz2_stream));

    if (mode[0] == 'r') {
        bz->fp = BZ2_bzReadOpen(&bzerror, raw_fp, 0, 0, NULL, 0);
        if (bzerror != BZ_OK) {
            BZ2_bzReadClose(&bzerror, bz->fp);
            fprintf(stderr, "[SysError:gz_stream_open:005] the file of (%s) does not seem to be a .bz2 file!\n", err_fn);
            exit(0);
        }
    }
    else {
        bz->is_write = 1;
        bz->fp = BZ2_bzWriteOpen(&bzerror, raw_fp, 9, 0, 30);
        if (bzerror != BZ_OK) {
            BZ2_bzWriteClose(&bzerror, bz->fp, 0, NULL, NULL);
            fprintf(stderr, "[SysError:gz_stream_open:012] failed to create file of (%s) in .bz2 format!\n", err_fn);
            exit(0);
        }
    }
    return bz;
}

static int bz2_read(void *fp, char *buf, int size)
{
    bz2_stream *bz = (bz2_stream *)fp;
    int bzerror;
    int n = BZ2_bzRead(&bzerror, bz->fp, buf, size);
    if (bzerror != BZ_OK && bzerror != BZ_STREAM_END) return -1;  /* truncated file detected */
    return n;  /* n==0 at the end of the stream (EOF) */
}

static int bz2_write(void *fp, const char *buf, int size)
{
    bz2_stream *bz = (bz2_stream *)fp;
    int bzerror;
    BZ2_bzWrite(&bzerror, bz->fp, (void *)buf, size);
    return bzerror == BZ_OK ? size : -1;
}

static int bz2_close(void *fp)
{
    bz2_stream *bz = (bz2_stream *)fp;
    int bzerror;

    if (bz->is_write)
        BZ2_bzWriteClose(&bzerror, bz->fp, 0, NULL, NULL);
    else
        BZ2_bzReadClose(&bzerror, bz->fp);
    free(bz);
    return bzerror == BZ_OK ? 0 : -1;
}

/* ---------------- plain text backend ---------------- */
static void *plain_handle_open(char *file, char *mode, char *err_fn)
{
    FILE *fp = fopen(file, mode);
    if (fp == NULL) {
        if (mode[0] == 'r')
            fprintf(stderr, "[SysError:gz_stream_open:016] failed to open a normal file of (%s)!\n", err_fn);
        else
            fprintf(stderr, "[SysError:gz_stream_open:013] failed to create a normal file of (%s)!\n", err_fn);
        exit(0);
    }
    return fp;
}

/* plain text reading (currently plain files are read
 * transparently by the gz backend, so this callback is not bound)
 */
static int plain_read(void *fp, char *buf, int size)
{
    size_t n = fread(buf, 1, size, (FILE *)fp);
    if (n == 0 && ferror((FILE *)fp)) return -1;  /* truncated file detected */
    return (int)n;
}

static int plain_write(void *fp, const char *buf, int size)
{
    return (int)fwrite(buf, 1, size, (FILE *)fp);
}

static int plain_close(void *fp)
{
    return fclose((FILE *)fp);
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
                gz->stream.fp = bz2_handle_open(f, "r", err_fn);
                gz->stream.read = bz2_read;
                gz->stream.close = bz2_close;
            } else {
                fclose(f);
                fprintf(stderr, "[SysError:gz_stream_open:006] the file of (%s) is neither a real bzip2(.bz2) file nor a stander fastq file!\n", err_fn);
                exit(0);
            }
        }
        else{ /* open a .gz or a normal txt file for reading */
            unsigned char buf[4];
            FILE *f = fopen(file, "rb");
            if (f==NULL) {
                fprintf(stderr, "[SysError:gz_stream_open:007] failed to open file of (%s)!\n", err_fn);
                exit(0);
            }
            fread(buf, 1, 4, f);  /* read the magic number of the gzip file */
            fclose(f);  /* the FILE handle is used for magic number detection only */
            if (buf[0]==0x1f && buf[1]==0x8b || buf[0]=='@') {
                gz->stream.fp = gz_handle_open(file, "r", err_fn);
                gz->stream.read = gz_read;
                gz->stream.close = gz_close;
            } else {
                fprintf(stderr, "[SysError:gz_stream_open:009] the file of (%s) is neither a real gzip(.gz) file nor a stander fastq file!\n", err_fn);
                exit(0);
            }
        }
    }
    else {  /* open a gz or bz2 file for writing */
        if(get_str_ends(file, ".gz")) {
            gz->stream.fp = gz_handle_open(file, "w", err_fn);
            gz->stream.write = gz_write;
            gz->stream.close = gz_close;
        }
        else if(get_str_ends(file, ".bz2")){
            FILE *f=fopen(file, "w");
            if(f==NULL){
                fprintf(stderr, "[SysError:gz_stream_open:011] failed to create file of (%s)!\n", err_fn);
                exit(0);
            }
            gz->stream.fp = bz2_handle_open(f, "w", err_fn);
            gz->stream.write = bz2_write;
            gz->stream.close = bz2_close;
        }
        else{
            gz->stream.fp = plain_handle_open(file, "w", err_fn);
            gz->stream.write = plain_write;
            gz->stream.close = plain_close;
        }
    }
    gz->buf=calloc(GZ_BUFF_SIZE, sizeof(char));

    return gz;
}


int gz_read_util(GzStream *gz, char delimiter, kstring_t *ks_str, int max_length)
{
    size_t len = 0;

    if (gz->is_eof && gz->begin >= gz->end)  /* end of the file */
        return 0;

    /* the body runs at least once, so the remaining data
    * in the buffer is scanned even when is_eof is already set */
    do {
        if (gz->begin >= gz->end) {  /* gz->buf is full or the first time to read */
            gz->begin = 0;
            gz->end = gz->stream.read(gz->stream.fp, gz->buf, GZ_BUFF_SIZE);
            if (gz->end < 0) return -1; /* truncated file detected */
            if (gz->end < GZ_BUFF_SIZE) gz->is_eof = 1;
        }
        while (gz->begin < gz->end) {  /* search the delimiter in block via memchr */
            const int n_remain = gz->end - gz->begin;
            const char *delim = (char *)memchr(gz->buf + gz->begin, delimiter, n_remain);
            int n_bytes = delim ? (int)(delim - (gz->buf + gz->begin)) : n_remain;

            if (len + n_bytes >= (size_t)max_length) {  /* line too long */
                fprintf(stderr, "[SysError:gz_read_util:014] read length can not be longer than %d!\n", max_length);
                fprintf(stderr, "[*] the main reason is that line breaks('\\n') can not be detected in the READ!\n");
                return -2;
            }
            if (ks_str->m < len + n_bytes + 1) {  /* grow the user-provided ks_str */
                ks_str->m = ks_str->m ? len + n_bytes + 4 : 512;
                kroundup32(ks_str->m);
                ks_str->s = (char *)realloc(ks_str->s, ks_str->m * sizeof(char));
                if (!ks_str->s) {
                    fprintf(stderr, "[SysError:gz_read_util:015] failed to reallocated memory!\n");
                    exit(-1);
                }
            }

            /* copy the str to the user-provided ks_str */
            memcpy(ks_str->s + len, gz->buf + gz->begin, n_bytes);
            len += n_bytes;
            gz->begin += n_bytes;

            if (delim) {  /* delimiter found */
                gz->begin++;  /* skip the delimiter */
                ks_str->s[len] = '\0'; ks_str->l = len;
                return 1;
            }
        }
    } while (!gz->is_eof);

    /* end of the file (EOF) detected; the remaining bytes without a
     * trailing delimiter belong to a problematic file, so they are
     * discarded and reported as an incomplete read by the caller */
    return 0;
}


int gz_read_block(GzStream *gz)
{
    if(gz->is_eof && gz->begin>=gz->end)  /* end of the file */
        return 0;

    gz->end = gz->stream.read(gz->stream.fp, gz->buf, GZ_BUFF_SIZE);

    if (gz->end < 0)  /* truncated file detected */
        return -1;

    if (gz->end < GZ_BUFF_SIZE) {  /* end of the file (EOF) detected */
        gz->is_eof = 1;
        return 0;
    }

    return 1;  /* normal reading of the block */
}


void gz_stream_destroy(GzStream *gz)
{
    if (gz->stream.write)  /* write mode: flush the buffered data first */
        gz->stream.write(gz->stream.fp, gz->buf, gz->begin);

    gz->stream.close(gz->stream.fp);  /* close the underlying file (also flushes its internal buffer) */
    if(gz->buf) free(gz->buf);
    free(gz);
}
