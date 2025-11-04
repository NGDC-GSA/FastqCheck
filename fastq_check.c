/*************************************************************************
    > File Name: fastq_check.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com
    > Created Time: 2022年07月26日 星期二 13时09分32秒
 ************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <ctype.h>

#include "bloom_filter.h"
#include "file_read.h"
#include "params_parse.h"
#include "fastq_check.h"
#include "file_type.h"


/* Static structure of the program for the whole life cycle */
static BaseBloomFilter stBloomFilter = {0};

static long get_max_file_size(FileObject *file_obj, int *max_file_idx)
{
    long file_size = -1;
    FILE *fp;
    char *err_fn;

    for (int idx=0; idx < file_obj->n; idx++) {
        long f_size;
        fp = fopen(file_obj->file_name[idx], "r");
        if (fp == NULL) {
            err_fn = get_path_basename(file_obj->file_name[idx]);
            fprintf(stderr, "[SysError:get_file_size:019] faild to read the given file %s\n", err_fn);
            return -1;
        }
        fseek(fp, 0L, SEEK_END);
        f_size = ftell(fp);
        if (f_size > file_size) {
            file_size = f_size; *max_file_idx = idx;
        }
        fclose(fp);
    }
    return file_size;
}


static int bloomfilter_size_init(FileObject *file_obj, cache_t *fastq_cache, int compress_ratio)
{
    /* get the max_file_size */
    int max_file_idx;
    long file_size = get_max_file_size(file_obj, &max_file_idx);
    read_t *f_reads = fastq_cache[max_file_idx].reads;
    uint64_t f_reads_num = fastq_cache[max_file_idx].n;

    /* get the minimum and maximum read length*/
    uint64_t n_total_bytes = 0;
    uint32_t min_read_len=1<<30, max_read_len=0, min_read_byte=1<<30;
    for (uint64_t i=0; i < f_reads_num; i++) {
        read_t *read = &f_reads[i];
        uint32_t n_read_byte = read->name.l + read->seq.l + read->comment.l + read->qual.l + 4;
        
        n_total_bytes += n_read_byte;
        if (n_read_byte < min_read_byte) min_read_byte = n_read_byte;

        /* get the minimum and maximum length of the read */
        if (read->seq.l < min_read_len) min_read_len = read->seq.l;
        if (read->seq.l > max_read_len) max_read_len = read->seq.l;
    }

    /* evaluate the number of reads with best way */
    uint64_t max_item, avg_read_byte;

    avg_read_byte = max_read_len>min_read_len ? min_read_byte : (n_total_bytes/f_reads_num);
    max_item = file_size / avg_read_byte * compress_ratio;
    InitBloomFilter(&stBloomFilter, 0, max_item, BLOOM_ERROR);

    return 1;
}


static cache_t *fastq_cache_init(int n_file)
{
    cache_t *fastq_cache, *fc;

    fastq_cache = (cache_t *)calloc(n_file, sizeof(cache_t));
    for (int i=0; i < n_file; i++) {
        fc = &fastq_cache[i];
        fc->n_max = CACHE_SIZE;
        fc->reads = (read_t *)calloc(fc->n_max, sizeof(read_t));
    }
    return fastq_cache;
}


static result_t *fastq_result_init(int n_file)
{
    result_t *results, *rt;
    length_t *len_obj;

    results = (result_t *)calloc(n_file, sizeof(result_t));
    for (int i=0; i < n_file; i++) {
        rt = &results[i];
        rt->phred = 33; /* initiate to phred33 */
        len_obj = &(rt->length_obj);
        len_obj->n_max = 512;  /* it will be resized when seq_len > 512 */
        len_obj->lengths = (uint64_t *)calloc(len_obj->n_max, sizeof(uint64_t));
    }
    return results;
}


static message_t *check_message_init(char *error_file, char *out_file)
{
    char *err_fn;
    message_t *msg = (message_t *)calloc(1, sizeof(message_t));

    /* other value initialization */
    msg->cache_status = 1;
    msg->n_duplicate = 0;
    msg->name_check = msg->length_check = 1;
    msg->pair_check = 1;  /* it could be valid only if args->pair_check=1 */

    /* open the error file (*.err) */
    msg->err_fp = fopen(error_file, "w");
    if (msg->err_fp == NULL) {
        err_fn = get_path_basename(error_file);
        fprintf(stderr, "[SysError:check_message_init:020] failed to create error file of (%s)!\n", err_fn); exit(-1);
    }
    /* open the output file (*.xml) */
    msg->out_fp = fopen(out_file, "w");
    if (msg->out_fp == NULL) {
        err_fn = get_path_basename(out_file);
        fprintf(stderr, "[SysError:check_message_init:021] failed to create output file of (%s)!\n", err_fn); exit(-2);
    }
    return msg;
}


/* read four line from the input fastq file */
static int fastq_read_core(GzStream *gz, read_t *read, int max_length)
{
    int ret;  /* the normal return is: ret==1 */

    ret = gz_read_util(gz, '\n', &read->name, max_length);  /* get the read name*/
    switch (ret) {
        case -2: return -3;  /* failed to detect line breaks('\n') */
        case -1: return -1;  /* unexpected end of the file */
        case  0: return 0;  /* end of the file or empty file */
        case  1: break;  /* normal reading the file */
    }

    ret = gz_read_util(gz, '\n', &read->seq, max_length);  /* get the read sequence */
    switch (ret) {
        case -2: return -3;  /* failed to detect line breaks('\n') */
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only one line is readed */
        case  1: break;  /* normal reading the file */
    }

    ret = gz_read_util(gz, '\n', &read->comment, max_length); /* get the read comment, usually is '+' */
    switch (ret) {
        case -2: return -3;  /* failed to detect line breaks('\n') */
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only two lines are readed */
        case  1: break;  /* normal reading the file */
    }

    ret = gz_read_util(gz, '\n', &read->qual, max_length); /* get the read quality */
    switch (ret) {
        case -2: return -3;  /* failed to detect line breaks('\n') */
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only three lines are readed */
        case  1: break;  /* normal reading the file */
    }

    if (read->name.s[0]=='@' && read->comment.s[0]=='+')
        return ret;  /* means read status is OK */
    
    return -2; /* incomplete fastq reads */
}


static int fastq_cache_read(FileObject *file_obj, uint64_t n_total_reads, int max_length, cache_t *fastq_cache, message_t *msg)
{
    int ret, finish;
    GzStream *f_gz;
    cache_t *f_cache;
    char *err_fn;
    uint64_t line_num;

    for (int i=0; i < file_obj->n; i++) {
        finish = 0;
        f_gz = file_obj->gz_hd[i]; 
        f_cache = &fastq_cache[i]; f_cache->n = 0;

        for (int idx=0; idx < f_cache->n_max && finish != 1; idx++) {
            ret = fastq_read_core(f_gz, &f_cache->reads[idx], max_length);
            switch (ret) {
            case 1:
                f_cache->n++; break; /* normal reading of the fastq file */
            case 0:
                finish = 1; break;  /* end of the file (EOF) */
            case -1:  
                return -1;  /* unexpected end of fastq file */
            case -2:
                line_num = (n_total_reads+f_cache->n) * 4 + 1;
                fprintf(msg->err_fp, "[FormatError:fastq_cache_read:201] incomplete fastq read '%s' is detected!\n", f_cache->reads[idx].name.s);
                fprintf(msg->err_fp, "[*] File%d: the line number of the error read is %lld!\n", i+1, line_num);
                return -2;
            case -3: /* failed to detect line breaks('\n') */
                line_num = (n_total_reads+f_cache->n) * 4 + 1;
                fprintf(msg->err_fp, "[FormatError:fastq_cache_read:212] failed to detect line breaks('\\n') in the READ!\n");
                fprintf(msg->err_fp, "[*] File%d: the line number of the error read is %lld!\n", i+1, line_num);
                return -3;
            }
        }
    }
    /* check whether all caches have same number of reads */
    for (int i=1; i < file_obj->n; i++) {
        if(fastq_cache[i].n != fastq_cache[0].n) {
            fprintf(msg->err_fp, "[FormatError:fastq_cache_read:202] the number of reads cached is different!\n");
            fprintf(msg->err_fp, "[*] File1: %lld  <--> File%d: %lld\n", n_total_reads+fastq_cache[0].n, i+1, n_total_reads+fastq_cache[i].n);
            return -4;
        }
    }
    /* check whether all the files have been finished */
    if (fastq_cache[0].n < fastq_cache[0].n_max) return 0;
    return 1;  /* status is OK, could caching for next time */
}


static int windows_break_check(cache_t *fastq_cache, int n_file, message_t *msg)
{
    read_t *f_reads;
    kstring_t *r_seq;  /* sequence of the read */
    
    for (int idx=0; idx < n_file; idx++) {
        f_reads = fastq_cache[idx].reads;

        for (int i=0; i < fastq_cache[idx].n; i++) {
            r_seq = &(f_reads[i].seq);
            if (r_seq->s[r_seq->l-1] == '\r') goto error;
        }
    }
    return 1;

    error:
        fprintf(msg->err_fp, "[FormatError:windows_break_check:203] windows break ('\\r\\n') is detected in the fastq file!\n");
        return 0;
}


static int phred_check(cache_t *fastq_cache, result_t *results, int n_file, message_t *msg)
{
    uint64_t *f_qual_table;  /* quality table of the file */
    read_t *f_reads;
    kstring_t *r_qual;  /* quality of the read */

    for (int idx=0; idx < n_file; idx++) {
        f_qual_table = results[idx].qual_table;
        f_reads = fastq_cache[idx].reads;

        /* record the quality and its number */
        for (int i=0; i < fastq_cache[idx].n; i++) {
            r_qual = &(f_reads[i].qual);
            for (int j=0; j < r_qual->l; j++) f_qual_table[r_qual->s[j]]++;
        }
    }
    /* get the phred of each file */
    int phred33, phred64;
    
    for (int idx=0; idx < n_file; idx++) {
        phred33 = 0; phred64 = 0;
        f_qual_table = results[idx].qual_table;
        for (int i=33; i < 59; i++) phred33 += f_qual_table[i];  /* only occured in phred33 */
        for (int i=75; i < 104; i++) phred64 += f_qual_table[i];  /* only occured in phred64 */

        /* check phred and clean the temp value within  f_qual_table */
        if (phred64 > 0 && phred33 == 0) results[idx].phred = 64;
        memset(f_qual_table, 0, sizeof(uint64_t)*256);
    }
    /* check whether the phred of all files are same */
    for (int idx=1; idx < n_file; idx++) {
        /* different phred occured */
        if (results[idx].phred != results[0].phred) {
            fprintf(msg->err_fp, "[FormatError:phred_check:204] the phred value of the given files is different!\n");
            for (int idx=0; idx < n_file; idx++) {
                fprintf(msg->err_fp, "  [*] File %d: Phred%lld\n", idx, results[idx].phred);
            }
            return 0;
        }
    }
    return 1;
}


static int read_name_trim(kstring_t *read_name, int pair_check)
{
    char c, marker, edge;
    int idx;

    for (idx=0; idx < read_name->l; idx++) {
        c = read_name->s[idx];
        switch (c) {
        case ' ':  /* eg. "@ST-E00126:256:1784 1:N:AGGGACG" */
            marker = read_name->s[idx+1];
            goto finish;
        
        case '/':  /* eg. "@HWI-ST833:189:5:2200#0/2" */
            edge = read_name->s[idx+2];  /* could be ' ', '\t', '\n' */
            if (isspace(edge) || edge=='\0') {
                marker = read_name->s[idx+1];
                goto finish;
            }
            break;  /* Note: "@m54218_12717/43257/6645_7416" */

        default:
            break;
        }
    }
    return 0;  /* No need to handle the read name */

    finish:
        if (pair_check) read_name->s[0] = marker;
        read_name->s[idx] = '\0'; read_name->l = idx;  /* set the first character to pari marker */
        return 0;
}


static char get_first_read_marker(const kstring_t *read_name)
{
    char c, marker, edge;
    int idx;

    for (idx=0; idx < read_name->l; idx++) {
        c = read_name->s[idx];
        switch (c) {
        case ' ':  /* eg. "@ST-E00126:256:1784 1:N:AGGGACG" */
            marker = read_name->s[idx+1];
            return marker;
        
        case '/':  /* eg. "@HWI-ST833:189:5:2200#0/2" */
            edge = read_name->s[idx+2];  /* could be ' ', '\t', '\n' */
            if (isspace(edge) || edge=='\0') {
                marker = read_name->s[idx+1];
                return marker;
            }
            break;  /* Note: "@m54218_12717/43257/6645_7416" */

        default:
            break;
        }
    } return read_name->s[0];  /* set the '@' as the marker */
}


static int fastq_statistics_process(cache_t *fastq_cache, result_t *results, int n_file, int pair_check, message_t *msg)
{
    read_t *f_reads;
    uint64_t *f_base_table, *f_qual_table, line_num;
    length_t *f_len_obj;
    kstring_t *r_name, *r_seq, *r_comment, *r_qual;
    int n_reads = fastq_cache[0].n;  /* because the reads number of all files is the sample */
    
    for (int idx=0; idx < n_file; idx++) {
        results[idx].n_reads += n_reads;  /* count the number of the reads */

        /* process each read in the caches of the file */
        f_base_table = results[idx].base_table; f_qual_table = results[idx].qual_table;
        f_len_obj = &results[idx].length_obj;
        f_reads = fastq_cache[idx].reads;
        char pair_marker = get_first_read_marker(&(f_reads[0].name));
        char file_marker = '1' + idx;

        if (pair_check && pair_marker!=file_marker && msg->pair_check) {
            fprintf(msg->err_fp, "[FormatError:fastq_statistics_process:205] the pair marker of File%d should be %c instead of %c!\n", idx+1, file_marker, pair_marker);
            msg->pair_check = 0;
        }
        
        for (int i=0; i < n_reads; i++) {
            r_name = &(f_reads[i].name); r_seq = &(f_reads[i].seq);
            r_comment = &(f_reads[i].comment); r_qual = &(f_reads[i].qual);
            /* calculate the number of bytes of the read */
            results[idx].n_bytes += (r_name->l + r_seq->l + r_comment->l + r_qual->l + 4);

            /* check the format of the read (check start of '@' and '+') */
            if (r_name->s[0] != '@' || r_comment->s[0] != '+') {
                line_num = (results[idx].n_reads - n_reads + i) * 4 + 1;
                fprintf(msg->err_fp, "[FormatError:fastq_statistics_process:206] [F(%d):L(%lld)] the format of the read '%s' is wrong (read name not start with '@' or comment not start with '+')!\n", idx+1, line_num, r_name->s);
                fclose(msg->err_fp); exit(-1);
            }
            /* check whether the length of the seq and quality is the same */
            if (r_seq->l != r_qual->l && msg->length_check) {
                line_num = (results[idx].n_reads - n_reads + i) * 4 + 1;
                fprintf(msg->err_fp, "[FormatError:fastq_statistics_process:207] [F(%d):L(%lld)] the format of the read '%s' is wrong (length of the sequence and quality is different)!\n", idx+1, line_num, r_name->s);
                msg->length_check = 0;
            }
            /* calculate the quality and base distribution */
            for (int j=0; j < r_seq->l; j++) {
                f_base_table[r_seq->s[j]]++; f_qual_table[r_qual->s[j]]++;
            }
            /* only keep the read_name without the second part (like -> 1:N:0:AAGACG) */
            read_name_trim(r_name, pair_check);
            if (pair_check && r_name->s[0]!=pair_marker && msg->pair_check) {
                line_num = (results[idx].n_reads - n_reads + i) * 4 + 1;
                fprintf(msg->err_fp, "[FormatError:fastq_statistics_process:208] [F(%d):L(%lld)] the format of the read '%s' is wrong (its pair_marker '%c' is different to others '%c')!\n", idx+1, line_num, f_reads[i].name.s, r_name->s[0], pair_marker);
                msg->pair_check = 0;
            }
            /* calculate the length distribution of all the files */
            if (r_seq->l >= f_len_obj->n_max) {  /* resize the lengths array */
                while (r_seq->l >= f_len_obj->n_max) f_len_obj->n_max = f_len_obj->n_max<<2;
                f_len_obj->lengths = (uint64_t *)realloc(f_len_obj->lengths, f_len_obj->n_max * sizeof(uint64_t));
                if (!f_len_obj->lengths) {
                    fprintf(stderr, "[SysError:fastq_statistics_process:022] failed to allocated memory when resize lengths array!\n");
                    exit(-1);
                }
            }
            f_len_obj->lengths[r_seq->l]++;
        }
    }
    return 0;
}


static int read_name_check(cache_t *fastq_cache, uint64_t n_total_reads, int n_file, message_t *msg)
{
    kstring_t *r_name; /* get the read name of the first file */
    int ret, n_reads=fastq_cache[0].n;  /* because the reads number of all files is the sample */
    uint64_t line_num;

    for (int i=0; i < n_reads; i++) {
        r_name = &(fastq_cache[0].reads[i].name);

        for (int idx=1; idx < n_file; idx++) {
            ret = strcmp(r_name->s+1, fastq_cache[idx].reads[i].name.s+1);
            if (ret != 0 && msg->name_check) {
                line_num = (n_total_reads - n_reads + i) * 4 + 1;
                fprintf(msg->err_fp, "[FormatError:read_name_check:209] [L(%lld)] the read name is not in the same order between two fastq file!\n", line_num);
                fprintf(msg->err_fp, "  [*] File1: '@%s' <--> File%d: '@%s'\n", r_name->s+1, idx+1, fastq_cache[idx].reads[i].name.s+1);
                msg->name_check = 0;
            }
        }
        /* put the read_name into the global bloomfilter */
        ret = BloomFilter_Check(&stBloomFilter, r_name->s, r_name->l);
        switch (ret) {  /* ret -> 0:existed; 1:not existed; -N:failed */
        case 0:
            fprintf(stdout, "[Warning:read_name_check] the read name '@%s' may be a duplicate with a probability of 1e-9!\n", r_name->s+1);
            if (++msg->n_duplicate > READ_NAME_DUPLICATE_MAX) {
                fprintf(msg->err_fp, "[FormatError:read_name_check:210] the number of duplicate read name is larger than %d (details in STANDOUT or Screen)!\n", READ_NAME_DUPLICATE_MAX);
                fclose(msg->err_fp); exit(-1);
            }
            break;
        case 1:
            ret = BloomFilter_Add(&stBloomFilter, r_name->s, r_name->l);
            if (ret != 0 && msg->n_duplicate >= READ_NAME_DUPLICATE_MAX) {
                fprintf(stderr, "[SysError:read_name_check:023] the compression ratio of fastq file is greater than default parameter 10!\n");
                fprintf(stderr, "  [*] Please set the parameter -r (--ratio) greater than the default value 10, e.g. 15, 20, etc!\n");
                exit(-1);
            }
            break;
        default:  /* not possible to be here */
            fprintf(stderr, "[SysError:read_name_check:024] it is NOT POSSIBLE to be here, contact the ADMINISTRATOR to check why!\n");
            break;
        }
    }
    return 0;
}

static void strip_access_path(char *access_path)
{
    /* find the last path marker '/' */
    char *r_marker = strrchr(access_path, '/');
    if (r_marker == NULL) /* the access_path may be just a filename */
        return ;

    *r_marker = '\0'; /* trim the real filename */

    /* strip the prefix before CRA number */
    char *find = strstr(access_path, "/CRA");
    if (find != NULL) {
        size_t n = r_marker - find - 1;
        strncpy(access_path, find+1, n);
        access_path[n] = '\0';
    }
}


static statistic_t *calculate_result_summary(result_t *results, char *file_name, int n_file, message_t *msg)
{
    statistic_t *statis;

    { /* statis structure initate */
        statis = (statistic_t *)calloc(1, sizeof(statistic_t));
        statis->n_stat = n_file;
        statis->phred = results[0].phred;
        statis->stat_obj = (stat_file_t *)calloc(statis->n_stat, sizeof(stat_file_t));
        strcpy(statis->length_status, "fixed");
        statis->access_path = (char *)malloc(strlen(file_name)+1);
        strcpy(statis->access_path, file_name); 
    }
    for (int idx=0; idx < n_file; idx++) {
        stat_file_t *f_stat_obj = &(statis->stat_obj[idx]);
        result_t *f_result = &(results[idx]);
        length_t *f_len_obj = &(f_result->length_obj);
        uint64_t n_valid_bases = 0;

        /* calculate the number of reads and total bytes for all files */
        statis->n_reads = f_result->n_reads;
        statis->n_total_bytes += f_result->n_bytes;
        for (int i=0; i < 256; i++) f_stat_obj->n_bases += f_result->base_table[i];
        statis->n_total_bases += f_stat_obj->n_bases;

        /* calculate the "AaCcTtGgNn" number and quality table and check whether has bases other than 'ACTGN' */
        char target_base[16] = "AaCcGgTtNn";
        for (int i=0; i < 10; i++) {
            statis->n_bases_num[i/2] += f_result->base_table[target_base[i]];
            n_valid_bases += f_result->base_table[target_base[i]]; /* valid 'ACGTN' for this fastq file */
        }
        for (int i=0; i < 256; i++) statis->qual_table[i] += f_result->qual_table[i];
        if (n_valid_bases != f_stat_obj->n_bases) {
            float invalid_base_ratio = (float)(f_stat_obj->n_bases - n_valid_bases) / (float)f_stat_obj->n_bases;
            if (invalid_base_ratio > MAX_INVALID_BASE_RATIO)
                fprintf(msg->err_fp, "[FormatError:calculate_result_summary:211] the number of invalid bases (not 'ACGTN') in File%d exceeds %.3f!\n", idx+1, MAX_INVALID_BASE_RATIO);
        }

        /* calculate average length and standard deviation for each file */
        for (int len=0; len < f_len_obj->n_max; len++) f_stat_obj->length_avg += f_len_obj->lengths[len] * len;
        f_stat_obj->length_avg /= (double)statis->n_reads;

        double diff, sum=0.0;
        for (int len=0; len < f_len_obj->n_max; len++) {
            uint64_t len_count = f_len_obj->lengths[len];
            if (len_count == 0) continue;
            diff = (double)len - f_stat_obj->length_avg;
            sum += pow(diff, 2.0) * len_count;
        }
        f_stat_obj->length_stdev = pow(sum/(double)statis->n_reads, 0.5);
        if (f_stat_obj->length_stdev != 0.0) strcpy(statis->length_status, "variable");
    }
    /* calculate other filed or copy the content to statis object */
    statis->gc_content = (double)(statis->n_bases_num[1] + statis->n_bases_num[2]) / statis->n_total_bases;
    strip_access_path(statis->access_path);
    
    return statis;
}


static void write_result_summary(statistic_t *statistics, int n_file, message_t *msg)
{
    stat_file_t *f_stat_obj;

    /* basic statistics information */
    fprintf(msg->out_fp, "<Run accession=\"%s\" read_length=\"%s\" spot_count=\"%lld\" base_count=\"%lld\">\n", statistics->access_path, statistics->length_status, statistics->n_reads, statistics->n_total_bases);
    fprintf(msg->out_fp, "  <Size value=\"%lld\" units=\"bytes\"/>\n", statistics->n_total_bytes);
    fprintf(msg->out_fp, "  <Bases cs_native=\"false\" count=\"%lld\">\n", statistics->n_total_bases);
    fprintf(msg->out_fp, "    <Base value=\"A\" count=\"%lld\"/>\n", statistics->n_bases_num[0]);
    fprintf(msg->out_fp, "    <Base value=\"C\" count=\"%lld\"/>\n", statistics->n_bases_num[1]);
    fprintf(msg->out_fp, "    <Base value=\"G\" count=\"%lld\"/>\n", statistics->n_bases_num[2]);
    fprintf(msg->out_fp, "    <Base value=\"T\" count=\"%lld\"/>\n", statistics->n_bases_num[3]);
    fprintf(msg->out_fp, "    <Base value=\"N\" count=\"%lld\"/>\n", statistics->n_bases_num[4]);
    fprintf(msg->out_fp, "  </Bases>\n");
    fprintf(msg->out_fp, "  <GC-Content value=\"%.2f%%\"/>\n", statistics->gc_content*100.0);
    fprintf(msg->out_fp, "  <AlignInfo>\n");
    fprintf(msg->out_fp, "  </AlignInfo>\n");
    /* part of statistics */
    fprintf(msg->out_fp, "  <Statistics nreads=\"%lld\" nspots=\"%lld\">\n", statistics->n_stat, statistics->n_reads);
    for (int idx=0; idx < statistics->n_stat; idx++) {
        f_stat_obj = &(statistics->stat_obj[idx]);
        fprintf(msg->out_fp, "    <Read index=\"%d\" count=\"%lld\" bases=\"%lld\" average=\"%d\" stdev=\"%.2f\"/>\n", idx, statistics->n_reads, f_stat_obj->n_bases, (int)f_stat_obj->length_avg, f_stat_obj->length_stdev);
    }
    fprintf(msg->out_fp, "  </Statistics>\n");
    fprintf(msg->out_fp, "  <QualityCount>\n");
    /* quality distribution */
    for (int q=0; q < 256; q++) {
        if (statistics->qual_table[q] == 0) continue;
        fprintf(msg->out_fp, "<Quality value=\"%lld\" count=\"%lld\"/>\n", q - statistics->phred, statistics->qual_table[q]);
    }
    fprintf(msg->out_fp, "  </QualityCount>\n");
    /* part Dataase */
    fprintf(msg->out_fp, "  <Databases>\n");
    fprintf(msg->out_fp, "    <Database>\n");
    fprintf(msg->out_fp, "      <Table name=\"SEQUENCE\">\n");
    fprintf(msg->out_fp, "        <Statistics source=\"meta\">\n");
    fprintf(msg->out_fp, "          <Rows count=\"%lld\"/>\n", statistics->n_reads);
    fprintf(msg->out_fp, "          <Elements count=\"%lld\"/>\n", statistics->n_total_bases);
    fprintf(msg->out_fp, "        </Statistics>\n");
    fprintf(msg->out_fp, "      </Table>\n");
    fprintf(msg->out_fp, "    </Database>\n");
    fprintf(msg->out_fp, "  </Databases>\n");
    fprintf(msg->out_fp, "</Run>\n");
}


char *get_current_time(char *time_buf)
{
    time_t c_time;
    struct tm *tm_obj;
    int year, month, day, hour, minute, second;

    time(&c_time);
    tm_obj = gmtime(&c_time);

    year = tm_obj->tm_year + 1900;
    month = tm_obj->tm_mon + 1;
    day = tm_obj->tm_mday;
    hour = tm_obj->tm_hour + 8;
    minute = tm_obj->tm_min;
    second = tm_obj->tm_sec;

    sprintf(time_buf, "%d-%d-%d %d:%d:%d", year, month, day, hour, minute, second);
    return time_buf;
}


static int file_truncated_check(FileObject *file_obj, message_t *msg)
{
    GzStream *f_gz;
    char *err_fn;

    for (int i=0; i < file_obj->n; i++) {
        f_gz = file_obj->gz_hd[i];
        err_fn = get_path_basename(file_obj->file_name[i]);

        if (f_gz->gz_fp) {  /* is gziped file */
            while (1) {
                f_gz->end = gzread(f_gz->gz_fp, f_gz->buf, GZ_BUFF_SIZE);
                gzerror(f_gz->gz_fp, &(f_gz->bzerror));
                if (f_gz->bzerror < 0) {  /* truncated file detected */
                    fprintf(msg->err_fp, "[FileError:file_truncated_check:101] unexpected end of fastq file %s (truncated file) is detected!\n", err_fn); 
                    break;
                }
                if (f_gz->end < GZ_BUFF_SIZE) /* end of the file (EOF) */
                    break;
            }
        }
        else if (f_gz->bz2_fp) {  /* is bzip2 file */
            while (1) {
                f_gz->end = BZ2_bzRead(&(f_gz->bzerror), f_gz->bz2_fp, f_gz->buf, GZ_BUFF_SIZE);
                if (f_gz->bzerror == BZ_STREAM_END) /* end of the file (EOF) */
                    break;

                if (f_gz->bzerror != BZ_OK) { /* truncated file detected */
                    fprintf(msg->err_fp, "[FileError:file_truncated_check:101] unexpected end of fastq file %s (truncated file) is detected!\n", err_fn); 
                    break;
                }
            }
        }
    }
    return 0;
}


int main(int argc, char **argv)
{
    char time_buf[32];

    /* parse the parameters and initiate object of error message */
    arg_t *args = args_parse(argc, argv);
    message_t *msg = check_message_init(args->error_file, args->out_file);
    if (file_type_check(args->in_file, msg->err_fp) < 0) {
        fclose(msg->err_fp); exit(-1); /* something wrong with the file format */
    }

    /* base structure initiate */
    fprintf(stdout, "[%s] Reading the file list and initiating the basic data structure ...\n", get_current_time(time_buf));
    FileObject *file_obj = read_file_list(args->in_file);
    cache_t *fastq_cache = fastq_cache_init(file_obj->n);
    result_t *results = fastq_result_init(file_obj->n);
    int max_read_length = args->max_length;

    /* read the first cache and check the windows_break, phred, bloom_size */
    msg->cache_status = fastq_cache_read(file_obj, results[0].n_reads, max_read_length, fastq_cache, msg);
    if (msg->cache_status < 0) goto __final_process;

    windows_break_check(fastq_cache, file_obj->n, msg);
    phred_check(fastq_cache, results, file_obj->n, msg);
    bloomfilter_size_init(file_obj, fastq_cache, args->ratio);

    /* processing the fastq files for the rest reads */
    fprintf(stdout, "[%s] Reading and processing the fastq files ...\n", get_current_time(time_buf));
    do {
        /* process the fastq file */
        fastq_statistics_process(fastq_cache, results, file_obj->n, args->pair_check, msg);
        /* check the read name */
        read_name_check(fastq_cache, results[0].n_reads, file_obj->n, msg);
        /* continue read the fastq file except error occured */
        if (!msg->cache_status) break;
        msg->cache_status = fastq_cache_read(file_obj, results[0].n_reads, max_read_length, fastq_cache, msg);

        if (results[0].n_reads % 1000000 == 0)
            fprintf(stdout, "[*] Processing number of reads: %lld ...\n", results[0].n_reads);
            
    } while (msg->cache_status >= 0);
    fprintf(stdout, "[*] Processing total number of reads: %lld ...\n", results[0].n_reads);

    /* handle the error code from function of fastq_cache_read */
    __final_process:
    switch(msg->cache_status) {
    case 0: /* end of the file */
        /* calculate and write the results to out file when finished */
        fprintf(stdout, "[%s] Calculating and writing the result to the output file(.xml) ...\n", get_current_time(time_buf));
        statistic_t *statis = calculate_result_summary(results, file_obj->file_name[0], file_obj->n, msg);
        write_result_summary(statis, file_obj->n, msg);
        break;
    case -1: /* unexpected end of fastq file */
    case -2: /* incomplete fastq read detected */
    case -3: /* too long read or faild to detect line breaks (mainly) */
    case -4: /* number of reads cached is different */
        file_truncated_check(file_obj, msg);  /* ignore the file format, only checking whether the file is truncated */
    }

    fclose(msg->err_fp); fclose(msg->out_fp);
    fprintf(stdout, "[%s] Done!\n", get_current_time(time_buf));
}
