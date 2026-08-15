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
#include <math.h>
#include "bloom_filter.h"
#include "file_read.h"
#include "params_parse.h"
#include "fastq_check.h"
#include "file_type.h"
#include "threadpool.h"


static BaseBloomFilter *bloomfilter_memory_init(uint64_t max_reads, int compress_ratio)
{
    uint64_t est_max_reads = max_reads * 1000000;  /* the unit of max_reads is million */

    if (compress_ratio > 10)
        est_max_reads = est_max_reads / 10 * compress_ratio;

    BaseBloomFilter *bloom_filter = InitBloomFilter(est_max_reads, BLOOM_ERROR);
    if (bloom_filter == NULL) {
        fprintf(stderr, "[SysError:bloomfilter_size_init:022] failed to initialize bloom filter!\n");
        return NULL;
    }

    return bloom_filter;
}


static cache_t *fastq_cache_init(int n_file)
{
    cache_t *fastq_cache, *fc;

    fastq_cache = (cache_t *)calloc(n_file, sizeof(cache_t));
    if (fastq_cache == NULL) {
        fprintf(stderr, "[SysError:fastq_cache_init:020] failed to allocate fastq cache array!\n");
        return NULL;
    }

    for (int i=0; i < n_file; i++) {
        fc = &fastq_cache[i];
        fc->n_max = CACHE_SIZE;
        fc->reads = (read_t *)calloc(fc->n_max, sizeof(read_t));
        if (fc->reads == NULL) {
            fprintf(stderr, "[SysError:fastq_cache_init:021] failed to allocate reads cache for File%d!\n", i+1);
            return NULL;
        }
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
    msg->n_duplicate = 0;
    msg->name_check = msg->length_check = 1;
    msg->pair_check = 1;  /* it could be valid only if args->pair_check=1 */
    pthread_mutex_init(&msg->err_lock, NULL);

    /* open the error file (*.err) */
    msg->err_fp = fopen(error_file, "w");
    if (msg->err_fp == NULL) {
        err_fn = get_path_basename(error_file);
        fprintf(stderr,
                "[SysError:check_message_init:020] failed to create error file of (%s)!\n",
                err_fn);
        exit(-1);
    }
    /* open the output file (*.xml) */
    msg->out_fp = fopen(out_file, "w");
    if (msg->out_fp == NULL) {
        err_fn = get_path_basename(out_file);
        fprintf(stderr,
            "[SysError:check_message_init:021] failed to create output file of (%s)!\n",
            err_fn);
        exit(-2);
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
        default: break;  /* ret==1: normal reading the file */
    }

    ret = gz_read_util(gz, '\n', &read->seq, max_length);  /* get the read sequence */
    switch (ret) {
        case -2: return -3;  /* failed to detect line breaks('\n') */
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only one line is read */
        default: break;  /* ret==1: normal reading the file */
    }

    ret = gz_read_util(gz, '\n', &read->comment, max_length); /* get the read comment, usually is '+' */
    switch (ret) {
        case -2: return -3;  /* failed to detect line breaks('\n') */
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only two lines are readed */
        default: break;  /* ret==1: normal reading the file */
    }

    ret = gz_read_util(gz, '\n', &read->qual, max_length); /* get the read quality */
    switch (ret) {
        case -2: return -3;  /* failed to detect line breaks('\n') */
        case -1: return -1;  /* unexpected end of the file */
        case  0: return -2;  /* incomplete fastq reads, since only three lines are readed */
        default: break;  /* ret==1: normal reading the file */
    }

    if (read->name.s[0]=='@' && read->comment.s[0]=='+')
        return ret;  /* means read status is OK */
    
    return -2; /* incomplete fastq reads */
}


static char read_name_trim(kstring_t *read_name) {
    char c, edge;
    int idx;

    for (idx=0; idx < read_name->l; idx++) {
        c = read_name->s[idx];
        switch (c) {
        case ' ':  /* eg. "@ST-E00126:256:1784 1:N:AGGGACG" */
            read_name->s[idx] = '\0';
            read_name->l = idx;
            return read_name->s[idx+1];

        case '/':  /* eg. "@HWI-ST833:189:5:2200#0/2" */
            edge = read_name->s[idx+2];  /* could be ' ', '\t', '\n' */
            if (isspace(edge) || edge=='\0') {
                read_name->s[idx] = '\0';
                read_name->l = idx;
                return read_name->s[idx+1];
            }
            break;  /* Note: "@m54218_12717/43257/6645_7416" */

        default:
            break;
        }
    }
    return read_name->s[0];  /* No need to handle the read name */
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
    }
    return read_name->s[0];  /* set the '@' as the marker */
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


static statistic_t *calculate_result_summary(thread_args_t *args)
{
    statistic_t *statis;
    thread_task_t *tasks = args->task_obj;

    /* statis structure initiate */
    statis = calloc(1, sizeof(statistic_t));
    statis->n_stat = tasks->n_file;
    statis->phred = tasks->phred_value;
    statis->stat_obj = calloc(statis->n_stat, sizeof(stat_file_t));
    strcpy(statis->length_status, "fixed");
    statis->access_path = (char *)malloc(strlen(args->file_obj->file_name[0])+1);
    strcpy(statis->access_path, args->file_obj->file_name[0]);

    /* calculate the result for each file */
    for (int32_t file_idx=0; file_idx < tasks->n_file; file_idx++) {
        stat_file_t *f_stat = &statis->stat_obj[file_idx];
        uint64_t n_valid_base = 0;
        statis->n_reads = 0;  /* the reads count is same for all files */

        for (int32_t task_idx=0; task_idx < args->task_size; task_idx++) {
            thread_task_t *task = &tasks[task_idx];
            result_t *f_result = &task->result[file_idx];

            /* calculate the number of reads and total bytes for all files */
            for (int i=0; i < 256; i++) f_stat->n_bases += f_result->base_table[i];
            for (int i=0; i < 256; i++) statis->qual_table[i] += f_result->qual_table[i];
            statis->n_total_bytes += f_result->n_bytes;
            statis->n_reads += f_result->n_reads;

            /* calculate the "AaCcTtGgNn" number and quality table */
            for (int i=0; i < 10; i++) {
                const char base[16] = "AaCcGgTtNn";
                statis->n_bases_num[i/2] += f_result->base_table[base[i]];
                n_valid_base += f_result->base_table[base[i]];
            }

            /* calculate average length and standard deviation for each file */
            length_t *f_len = &f_result->length_obj;
            for (int len=0; len < f_len->n_max; len++) f_stat->length_avg += (double)f_len->lengths[len] * len;
        }
        if (n_valid_base != f_stat->n_bases) {
            fprintf(args->msg_obj->err_fp,
                    "[FormatError:calculate_result_summary_core:211] the number of invalid bases (not 'ACGTN') in File%d exceeds %.3f!\n",
                    file_idx + 1, MAX_INVALID_BASE_RATIO);
        }
        statis->n_total_bases += f_stat->n_bases;

        /* calculate the length distribution after all tasks were calculated */
        double square_sum = 0.0;
        f_stat->length_avg /= (double)statis->n_reads;

        for (int32_t task_idx=0; task_idx < args->task_size; task_idx++) {
            result_t *f_result = &tasks[task_idx].result[file_idx];
            length_t *f_len = &f_result->length_obj;

            for (int len=0; len < f_len->n_max; len++) {
                if (f_len->lengths[len] == 0) continue;
                double diff = (double)len - f_stat->length_avg;
                square_sum += pow(diff, 2.0) * (double)f_len->lengths[len];
            }
        }
        f_stat->length_stdev = pow(square_sum/(double)statis->n_reads, 0.5);
        if (f_stat->length_stdev != 0.0) strcpy(statis->length_status, "variable");
    }

    /* calculate other filed or copy the content to statis object */
    statis->gc_content = (double)(statis->n_bases_num[1] + statis->n_bases_num[2]) / (double)statis->n_total_bases;
    strip_access_path(statis->access_path);

    return statis;
}


static void write_result_summary(statistic_t *stat, message_t *msg)
{
    /* basic statistics information */
    fprintf(msg->out_fp, "<Run accession=\"%s\" read_length=\"%s\" spot_count=\"%lld\" base_count=\"%lld\">\n", stat->access_path, stat->length_status, stat->n_reads, stat->n_total_bases);
    fprintf(msg->out_fp, "  <Size value=\"%lld\" units=\"bytes\"/>\n", stat->n_total_bytes);
    fprintf(msg->out_fp, "  <Bases cs_native=\"false\" count=\"%lld\">\n", stat->n_total_bases);
    fprintf(msg->out_fp, "    <Base value=\"A\" count=\"%lld\"/>\n", stat->n_bases_num[0]);
    fprintf(msg->out_fp, "    <Base value=\"C\" count=\"%lld\"/>\n", stat->n_bases_num[1]);
    fprintf(msg->out_fp, "    <Base value=\"G\" count=\"%lld\"/>\n", stat->n_bases_num[2]);
    fprintf(msg->out_fp, "    <Base value=\"T\" count=\"%lld\"/>\n", stat->n_bases_num[3]);
    fprintf(msg->out_fp, "    <Base value=\"N\" count=\"%lld\"/>\n", stat->n_bases_num[4]);
    fprintf(msg->out_fp, "  </Bases>\n");
    fprintf(msg->out_fp, "  <GC-Content value=\"%.2f%%\"/>\n", stat->gc_content*100.0);
    fprintf(msg->out_fp, "  <AlignInfo>\n");
    fprintf(msg->out_fp, "  </AlignInfo>\n");

    /* part of statistics */
    fprintf(msg->out_fp, "  <Statistics nreads=\"%lld\" nspots=\"%lld\">\n", stat->n_stat, stat->n_reads);
    for (int idx=0; idx < stat->n_stat; idx++) {
        stat_file_t *f_stat = &(stat->stat_obj[idx]);
        fprintf(msg->out_fp,
            "    <Read index=\"%d\" count=\"%lld\" bases=\"%lld\" average=\"%d\" stdev=\"%.2f\"/>\n",
            idx, stat->n_reads, f_stat->n_bases, (int)f_stat->length_avg, f_stat->length_stdev);
    }
    fprintf(msg->out_fp, "  </Statistics>\n");
    fprintf(msg->out_fp, "  <QualityCount>\n");

    /* quality distribution */
    for (int q=0; q < 256; q++) {
        if (stat->qual_table[q] == 0) continue;
        fprintf(msg->out_fp, "<Quality value=\"%lld\" count=\"%lld\"/>\n", q - stat->phred, stat->qual_table[q]);
    }
    fprintf(msg->out_fp, "  </QualityCount>\n");

    /* part database */
    fprintf(msg->out_fp, "  <Databases>\n");
    fprintf(msg->out_fp, "    <Database>\n");
    fprintf(msg->out_fp, "      <Table name=\"SEQUENCE\">\n");
    fprintf(msg->out_fp, "        <Statistics source=\"meta\">\n");
    fprintf(msg->out_fp, "          <Rows count=\"%lld\"/>\n", stat->n_reads);
    fprintf(msg->out_fp, "          <Elements count=\"%lld\"/>\n", stat->n_total_bases);
    fprintf(msg->out_fp, "        </Statistics>\n");
    fprintf(msg->out_fp, "      </Table>\n");
    fprintf(msg->out_fp, "    </Database>\n");
    fprintf(msg->out_fp, "  </Databases>\n");
    fprintf(msg->out_fp, "</Run>\n");
}


static char *get_current_time(char *time_buf)
{
    const time_t now = time(NULL);
    struct tm tm_obj;

    if (localtime_r(&now, &tm_obj) == NULL) {
        return NULL;
    }

    if (strftime(time_buf, 32, "%Y-%m-%d %H:%M:%S", &tm_obj) == 0) {
        return NULL;
    }

    return time_buf;
}


static thread_args_t *thread_args_init(arg_t *args)
{
    char time_buf[32];
    thread_args_t *thr_args;

    err_calloc(thr_args, 1, thread_args_t);
    thr_args->task_index = 0;
    thr_args->task_size = args->n_thread;

    fprintf(stdout, "[%s] Initiating the base data struct ...\n", get_current_time(time_buf));
    thr_args->msg_obj = check_message_init(args->error_file, args->out_file);
    thr_args->bloom_filter = bloomfilter_memory_init(args->max_reads, args->compress_ratio);
    thr_args->file_obj = read_file_list(args->in_file);
    thr_args->task_queue = kqueue_init(thr_args->task_size);
    err_calloc(thr_args->task_obj, thr_args->task_size, thread_task_t);

    /* initiate the thread_tasks object */
    for (int i=0; i < thr_args->task_size; i++) {
        thread_task_t *thr_task = &thr_args->task_obj[i];
        thr_task->n_file = thr_args->file_obj->n;
        thr_task->phred_value = args->phred_value;
        thr_task->max_length = args->max_length;
        thr_task->pair_check = args->pair_check;

        /* struct initiate for all files */
        thr_task->cache = fastq_cache_init(thr_task->n_file);
        thr_task->result = fastq_result_init(thr_task->n_file);
        thr_task->files_hd = thr_args->file_obj->gz_hd;
        thr_task->msg_obj = thr_args->msg_obj;
        err_calloc(thr_task->hashes, CACHE_SIZE<<1, uint64_t);

        /* member of pthread initiate */
        pthread_mutex_init(&thr_task->status_lock, NULL);
        pthread_cond_init(&thr_task->batch_ready, NULL);
    }

    return thr_args;
}


static thread_task_t *thread_get_task(thread_args_t *thr_args)
{
    thread_task_t *thr_task;

    /* the kqueue_get_access is used to check whether kqueue is full */
    kqueue_get_access(thr_args->task_queue);  /* kqueue is not full */
    thr_task = &thr_args->task_obj[thr_args->task_index];
    thr_args->task_index = (thr_args->task_index + 1) % thr_args->task_size;

    return thr_task;
}


static void thread_fastq_reader_core(void *args)
{
    reader_job_t *job = (reader_job_t *)args;
    thread_task_t *task = job->task;
    cache_t *f_cache = &task->cache[job->file_idx];

    /* read the fastq and put them into cache object */
    GzStream *file_hd = task->files_hd[job->file_idx];
    const int32_t batch_mask = (1<<task->n_file) - 1;
    f_cache->n = 0;

    for (int32_t idx=0; idx < f_cache->n_max; idx++) {
        job->ret_value = fastq_read_core(file_hd, &f_cache->reads[idx], task->max_length);
        if (job->ret_value != 1)
            /*  0: end of the file (EOF)
             * -1: unexpected end of fastq file
             * -2: incomplete fastq read is detected
             * -3: failed to detect line breaks ('\n')
             */
            break;

        f_cache->n++;  /* normal reading of the fastq file */
    }

    /* set the file ready mask for the file batch (batch means pair-end files) */
    pthread_mutex_lock(&task->status_lock);
    task->status |= 1 << job->file_idx;
    if (task->status == batch_mask)  /* all files of the batch cached successfully */
        pthread_cond_broadcast(&task->batch_ready);
    pthread_mutex_unlock(&task->status_lock);
}


static int fastq_reader_error_process(reader_job_t *jobs, thread_task_t *task,  message_t *msg)
{
    int32_t ret_code = 1;  /* 1: normal caching of reads */
    int32_t n_file = task->n_file;

    for (int file_idx=0; file_idx < n_file; file_idx++) {
        reader_job_t *job = &jobs[file_idx];
        cache_t *f_cache = &task->cache[file_idx];

        if (job->ret_value == 1 && f_cache->n == task->cache->n)
            continue;

        pthread_mutex_lock(&msg->err_lock);
        ret_code = -1;  /* -1: abnormal caching of reads */

        switch (job->ret_value) {
            case 0:  /* EOF of the file */
                ret_code = f_cache->n > 0 ? 1 : 0;
                break;

            case 1:  /* cached different number of reads */
                fprintf(msg->err_fp,
                        "[FormatError:fastq_reader_error_handle:202] the number of reads cached is different!\n");
                fprintf(msg->err_fp,
                        "[*] File1: %lu <--> File%d: %lu\n",
                        task->cache->n, file_idx + 1, f_cache->n);
                break;

            case -2:
                fprintf(msg->err_fp,
                        "[FormatError:fastq_reader_error_handle:201] incomplete fastq read '%s' is detected!\n",
                        f_cache->reads[f_cache->n].name.s);
                fprintf(msg->err_fp,
                        "[*] File%d: the line number of the error read is %lld!\n",
                        file_idx + 1, (task->total_reads + f_cache->n) * 4 + 1);
                break;

            case -3:  /* failed to detect line breaks('\n') */
                fprintf(msg->err_fp,
                        "[FormatError:fastq_reader_error_handle:212] failed to detect line breaks('\\n') in the READ!\n");
                fprintf(msg->err_fp,
                        "[*] File%d: the line number of the error read is %lld!\n",
                        file_idx + 1, (task->total_reads + f_cache->n) * 4 + 1);
                break;

            default:  /* (ret_value == -1) unexpected end of the fastq file */
                break;
        }
        pthread_mutex_unlock(&msg->err_lock);
    }

    return ret_code;
}


static int32_t fastq_reader_cache_handle(threadpool_t *thr_pool, reader_job_t *jobs, thread_task_t *task)
{
    int32_t n_file = task->n_file;
    const int32_t batch_mask = (1<<n_file) - 1;
    task->status = task->n_hash = 0;  /* the task is ready to cache reads */

    /* caching reads for all files in parallel */
    for (int file_idx=0; file_idx < n_file; file_idx++) {
        reader_job_t *job = &jobs[file_idx];
        job->task = task;
        job->file_idx = file_idx;
        job->ret_value = 0;
        threadpool_add(thr_pool, thread_fastq_reader_core, job, 0);
    }

    /* submit fastq check task for read format checking */
    pthread_mutex_lock(&task->status_lock);
    while (task->status != batch_mask)
        pthread_cond_wait(&task->batch_ready, &task->status_lock);
    pthread_mutex_unlock(&task->status_lock);

    /* the total reads is the previous accumulated reads count */
    task->total_reads = jobs->read_offset;
    jobs->read_offset += task->cache->n;

    /* handle errors and EOF during fastq caching */
    return fastq_reader_error_process(jobs, task, task->msg_obj);
}


static void thread_truncate_check_core(void *args)
{
    reader_job_t *job = (reader_job_t *)args;
    thread_task_t *task = job->task;

    /* read the fastq in block */
    GzStream *file_hd = task->files_hd[job->file_idx];
    const int32_t batch_mask = (1<<task->n_file) - 1;

    while (1) {
        job->ret_value = gz_read_block(file_hd);
        if (job->ret_value > 0)  /* normal caching of the block */
            continue;

        /* truncated detected or EOF of the file */
        break;
    }

    /* set the file ready mask for the file batch (batch means pair-end files) */
    pthread_mutex_lock(&task->status_lock);
    task->status |= 1 << job->file_idx;
    if (task->status == batch_mask)  /* all files of the batch cached successfully */
        pthread_cond_broadcast(&task->batch_ready);
    pthread_mutex_unlock(&task->status_lock);
}


static void fastq_truncate_error_handle(threadpool_t *thr_pool, reader_job_t *jobs, thread_task_t *task)
{
    int32_t n_file = task->n_file;
    const int32_t batch_mask = (1<<n_file) - 1;
    task->status = 0;  /* the task is ready to check truncate */

    /* caching reads for all files in parallel */
    for (int file_idx=0; file_idx < n_file; file_idx++) {
        reader_job_t *job = &jobs[file_idx];
        job->task = task;
        job->file_idx = file_idx;
        job->ret_value = 0;
        threadpool_add(thr_pool, thread_truncate_check_core, job, 0);
    }

    /* submit fastq truncate checking task */
    pthread_mutex_lock(&task->status_lock);
    while (task->status != batch_mask)
        pthread_cond_wait(&task->batch_ready, &task->status_lock);
    pthread_mutex_unlock(&task->status_lock);

    /* handle errors from truncated checking */
    for (int file_idx=0; file_idx < n_file; file_idx++) {
        reader_job_t *job = &jobs[file_idx];
        if (job->ret_value == 0)  /* EOF of the file */
            continue;

        /* truncated file detected */
        pthread_mutex_lock(&task->msg_obj->err_lock);
        fprintf(task->msg_obj->err_fp,
                "[FileError:file_truncate_error_handle:] truncated fastq file (File%d) detected!\n",
                file_idx + 1);
        pthread_mutex_unlock(&task->msg_obj->err_lock);
    }
}


static void thread_fastq_check_core(void *args)
{
    thread_task_t *task = (thread_task_t*)args;
    message_t *msg = task->msg_obj;

    for (int file_idx=0; file_idx < task->n_file; file_idx++) {
        result_t *f_result = &task->result[file_idx];
        cache_t *f_cache = &task->cache[file_idx];

        f_result->n_reads += f_cache->n;
        char pair_marker = get_first_read_marker(&f_cache->reads->name);
        char file_marker = (char)('1' + file_idx);

        /* check the read pair marker and ensure they are the same as the file order */
        if (task->pair_check && pair_marker != file_marker && msg->pair_check) {
            pthread_mutex_lock(&msg->err_lock);
            fprintf(msg->err_fp,
                    "[FormatError:thread_fastq_check_core:205] "
                    "the pair marker of File%d should be %c instead of %c!\n",
                    file_idx + 1, file_marker, pair_marker);
            msg->pair_check = 0;
            pthread_mutex_unlock(&msg->err_lock);
        }

        for (int read_idx=0; read_idx < f_cache->n; read_idx++) {
            read_t *read = &f_cache->reads[read_idx];
            f_result->n_bytes += (read->name.l + read->seq.l + read->comment.l + read->qual.l + 4);

            /* check whether the pair_marker are same for all reads */
            char read_marker = read_name_trim(&read->name);
            if (task->pair_check && read_marker != pair_marker && msg->pair_check) {
                pthread_mutex_lock(&msg->err_lock);
                fprintf(msg->err_fp,
                        "[FormatError:thread_fastq_check_core:208] "
                        "[F(%d):L(%lld)] the format of the read '%s' is wrong (its pair_marker '%c' is different to others '%c')!\n",
                        file_idx + 1, (task->total_reads + read_idx) * 4 + 1, read->name.s, read_marker, pair_marker);
                msg->pair_check = 0;
                pthread_mutex_unlock(&msg->err_lock);
            }

            /* calculate the hash value for bloomfilter (only the first file) */
            if (file_idx == 0) {
                int32_t idx = task->n_hash << 1;
                task->hashes[idx] = MurmurHash2_x64(read->name.s, (int)read->name.l, 0);
                task->hashes[idx + 1] = MurmurHash2_x64(read->name.s, (int)read->name.l, MIX_UINT64(task->hashes[idx]));
                task->n_hash++;
            }
            /* check whether the read is the same order as the first file */
            else {
                read_t *f0_read = &task->cache->reads[read_idx];
                int32_t ret_value = strcmp(f0_read->name.s, read->name.s);

                if (ret_value != 0 && msg->name_check) {
                    pthread_mutex_lock(&msg->err_lock);
                    fprintf(msg->err_fp,
                            "[FormatError:thread_fastq_check_core:209] "
                            "[L(%lld)] the read name is not in the same order between two fastq file!\n",
                            (task->total_reads + read_idx) * 4 + 1);
                    fprintf(msg->err_fp,
                            "  [*] File1: '%s' <--> File%d: '%s'\n",
                            f0_read->name.s, file_idx + 1, read->name.s);
                    msg->name_check = 0;
                    pthread_mutex_unlock(&msg->err_lock);
                }
            }

            /* check whether the sequence length is the same as quality */
            if (read->seq.l != read->qual.l && msg->length_check) {
                pthread_mutex_lock(&msg->err_lock);
                fprintf(msg->err_fp,
                        "[FormatError:thread_fastq_check_core:207] "
                        "[F(%d):L(%lld)] the format of the read '%s' is wrong (length of the sequence and quality is different)!\n",
                        file_idx + 1, (task->total_reads + read_idx) * 4 + 1, read->name.s);
                msg->length_check = 0;
                pthread_mutex_unlock(&msg->err_lock);
            }

            /* calculate the quality and base distribution */
            for (int base_idx=0; base_idx < read->seq.l; base_idx++) {
                f_result->base_table[read->seq.s[base_idx]]++;
                f_result->qual_table[read->qual.s[base_idx]]++;
            }

            /* calculate the length distribution of all the files */
            if (read->seq.l >= f_result->length_obj.n_max) {
                size_t prev_size = f_result->length_obj.n_max;
                size_t new_size = read->seq.l + 1;

                kroundup32(new_size);
                f_result->length_obj.n_max = new_size;
                err_realloc(f_result->length_obj.lengths, new_size, uint64_t);
                memset(f_result->length_obj.lengths + prev_size, 0, new_size - prev_size);
            }
            f_result->length_obj.lengths[read->seq.l]++;
        }
    }

    /* the task status is a mask for different stages
     * single-end (1 file): 1 (1<<1 - 1) -> 11 (content checked) -> 111 (duplicate checked)
     * pair-end (2 files): 11 (1<<2 - 1) -> 111 (content checked) -> 1111 (duplicate checked)
     * single-cell (4 files): 1111 (1<<4 - 1) -> 11111 (content checked) -> 111111 (duplicate checked)
     */
    task->status = task->status << 1 | 0x1;
}


static void *thread_fastq_reader(void *args)
{
    thread_args_t *thr_args = (thread_args_t *)args;
    int32_t n_file = thr_args->file_obj->n;

    /* create the thread pool for file reading and processing */
    threadpool_t *reader_pool = threadpool_create(n_file, n_file, 0);
    threadpool_t *checker_pool = threadpool_create(thr_args->task_size, thr_args->task_size, 0);

    /* create the reader job for reads caching of each file */
    reader_job_t *reader_jobs;
    err_calloc(reader_jobs, n_file, reader_job_t);

    while (1) {
        /* get one task from kqueue when it is not full */
        thread_task_t *thr_task = thread_get_task(thr_args);
        int ret_value = fastq_reader_cache_handle(reader_pool, reader_jobs, thr_task);

        /* normal caching of the reads */
        if (ret_value > 0) {
            kqueue_push(thr_args->task_queue, thr_task);
            threadpool_add(checker_pool, thread_fastq_check_core, thr_task, 0);
            continue;
        }

        /* truncate checking when there have errors */
        if (ret_value < 0) {
            fastq_truncate_error_handle(reader_pool, reader_jobs, thr_task);
            kqueue_set_finish(thr_args->task_queue);
            break;
        }

        /* set finish signal when EOF (ret_value == 0) */
        kqueue_set_finish(thr_args->task_queue);
        break;
    }
    return NULL;
}


static void duplicate_error_process(thread_task_t *task, int32_t read_idx, int32_t ret_value)
{
    /* the item may existed in the bloomfilter */
    message_t *msg = task->msg_obj;

    if (ret_value == 1) {
        if (++msg->n_duplicate < READ_NAME_DUPLICATE_MAX)
            return;

        /* there may be duplicate in the files */
        pthread_mutex_lock(&msg->err_lock);
        cache_t *f_cache = task->cache;

        fprintf(msg->err_fp,
                "[FormatError:duplicate_error_process:210] "
                "[F(1):L(%lld)] the number of duplicate read name is larger than %d (details in STANDOUT or Screen)!\n",
                (task->total_reads + read_idx) * 4 + 1, READ_NAME_DUPLICATE_MAX);
        fprintf(msg->err_fp,
                "  [*] one of the duplicate read name: %s\n",
                f_cache->reads[read_idx].name.s);
        pthread_mutex_unlock(&msg->err_lock);
        return;
    }

    /* there is no space remain for more items */
    if (ret_value == -1) {
        pthread_mutex_lock(&msg->err_lock);
        fprintf(msg->err_fp,
                "[SysError:duplicate_error_process:023] "
                "the compression ratio of fastq file is greater than default parameter 10!\n");
        fprintf(msg->err_fp,
                "  [*] Please set the parameter -r (--ratio) greater than the default value 10, e.g. 15, 20, etc!\n");
        pthread_mutex_unlock(&msg->err_lock);
    }
}


static void *thread_duplicate_check(void *args)
{
    thread_args_t *thr_args = (thread_args_t *)args;
    thread_task_t *task;

    while (1) {
        /* To keep the order of FASTQ checking, process and pop out the task
         * from the head of the task queue each time
         * */
        if ((task = kqueue_get_front(thr_args->task_queue)) == NULL)
            break;

        const int32_t ready_mask = (1 << (task->n_file + 1)) - 1;
        if (task->status != ready_mask)
            continue;

        /* the front task is finished checking */
        size_t n_reads = task->cache->n;

        for (int read_idx=0; read_idx < n_reads; read_idx++) {
            int32_t hash_idx = read_idx << 1;
            int32_t ret_value = BloomFilterCheckAdd(thr_args->bloom_filter,task->hashes[hash_idx],task->hashes[hash_idx+1]);

            /* the item is insert into bloomfilter successfully */
            if (ret_value == 0)
                continue;

            duplicate_error_process(task, read_idx, ret_value);
        }

        /* show progress of the processing */
        if (task->total_reads % 1000000 == 0)
            fprintf(stdout, "[*] Processing number of reads: %lld ...\n", task->total_reads);

        /* pop out the finished task from task queue */
        kqueue_pop(thr_args->task_queue);
    }

    return NULL;
}


int main(int argc, char **argv)
{
    char time_buf[32];
    pthread_t reader_tid, checker_tid;

    /* parse the parameters and initiate object of error message */
    arg_t *args = args_parse(argc, argv);
    thread_args_t *thr_args = thread_args_init(args);

    /* create reader and writer threads */
    pthread_create(&reader_tid, NULL, thread_fastq_reader, (void*)thr_args);
    pthread_create(&checker_tid, NULL, thread_duplicate_check, (void*)thr_args);
    pthread_join(reader_tid, NULL);
    pthread_join(checker_tid, NULL);

    /* calculate and write the results to out file when finished */
    fprintf(stdout, "[%s] Calculating and writing the result to XML file ...\n", get_current_time(time_buf));
    statistic_t *statis = calculate_result_summary(thr_args);
    write_result_summary(statis, thr_args->msg_obj);

    /* close the file handles */
    fclose(thr_args->msg_obj->err_fp);
    fclose(thr_args->msg_obj->out_fp);
    fprintf(stdout, "[%s] Done!\n", get_current_time(time_buf));
}
