/*************************************************************************
    > File Name: params_parse.c
    > Author: xlzh
    > Mail: xiaolongzhang2015@163.com 
    > Created Time: 2022年07月27日 星期三 11时57分25秒
 ************************************************************************/

#include "params_parse.h"


static void show_usage(void)
{
    char *usage =
        "Usage: fastq_check -i <fastq_list.txt> -o <output_file.xml> -e <error_file.err> -x <max_reads> -P <phred_value>\n"
        "Note: you MUST RUN fastq_type for file type checking before fastq_check!\n\n"
        "Options:\n"
        "       -h|--help                 print help information\n\n"
        "[Required]\n"
        "       -i|--in_file        FILE  fastq file list, one sample path per line [.txt]\n"
        "       -o|--out_file       FILE  output xml file of the given fastq files after checking [.xml]\n"
        "       -e|--error_file     FILE  warning or error message about the given fastq files [.err]\n"
        "       -x|--max_reads      INT   the estimated maximum number of reads (estimated by fastq_type) [unit: million]\n"
        "       -P|--phred_value    INT   the phred value of the fastq file (estimated by fastq_type, 33 or 64) [required]\n\n"
        "[Optional]\n"
        "       -r|--ratio          INT   the estimated compression ratio of the fastq file [default:10]\n"
        "       -m|--max_length     INT   the maximum length (MB) allowed for one read sequence [default:50 (MB)]\n"
        "       -t|--thread         INT   number of thread used when processing the fastq file [default:1]\n"
        "       -p|--pair_check     INT   check whether all reads have the same pair marker (1->check; 0->ignore) [default:1]\n\n";

    fprintf(stderr, "Program: FastqCheck (v%s)\n", __FASTQ_CHECK_VERSION__);
    fprintf(stderr, "CreateDate: %s\n", __FASTQ_CHECK_CREATE_DATE__);
    fprintf(stderr, "UpdateDate: %s\n", __FASTQ_CHECK_UPDATE_DATE__);
    fprintf(stderr, "Author: XiaolongZhang (zhangxiaolong@big.ac.cn)\n\n");
    fprintf(stderr, "%s", usage);
    exit(-1);
}


/* copy a string and allocate necessary memory */
static char *x_strcopy(const char *str)
{
    char *new_str;

    size_t str_len = strlen(str) + 1;
    new_str = (char *)malloc(str_len * sizeof(char));
    if (!new_str) {
        fprintf(stderr, "[SysError:x_strcopy:016] failed to malloc memory for string %s!\n", str);
        exit(-1);
    }
    memcpy(new_str, str, str_len);

    return new_str;
}


static int get_file_number(char *file_list)
{
    char buf[1024];
    FILE *file_fp = fopen(file_list, "r");
    int file_number = 0;

    if (!file_fp) {
        fprintf(stderr, "[SysError:get_file_number:002] failed to open the file list when obtaining the number of files\n");
        exit(-1);
    }
    
    while (fgets(buf, 1024, file_fp)) {
        if (buf[0]=='\r' || buf[0]=='\n') continue;  /* skip blank line */
        ++file_number;
    }
    fclose(file_fp);

    return file_number;
}


static const struct option long_options[] =
{
    { "help", no_argument, NULL, 'h' },
    { "in_file", required_argument, NULL, 'i' },
    { "out_file", required_argument, NULL, 'o' },
    { "error_file", required_argument, NULL, 'e' },
    { "max_reads", required_argument, NULL, 'x' },
    { "phred_value", required_argument, NULL, 'P' },
    { "ratio", optional_argument, NULL, 'r' },
    { "max_length", optional_argument, NULL, 'm' },
    { "pair_check", optional_argument, NULL, 'p' },
    { "n_thread", optional_argument, NULL, 't' },
    { NULL, 0, NULL, 0 }
};


arg_t *args_parse(int argc, char **argv)
{
    int opt = 0;
    arg_t *args = calloc(1, sizeof(arg_t));

    /* set the default parameters */
    args->compress_ratio = 10;  // 10% of the original fastq file
    args->max_length = 50 * 1024 * 1024; // 50 MB
    args->pair_check = 1;
    args->n_thread = 1;

    if (argc < 11) show_usage();
    while ( (opt = getopt_long(argc, argv, "i:o:e:x:P:r:m:p:t:h", long_options, NULL)) != -1 )
    {
        switch (opt) {
            case 'h': args->help = 1; break;
            case 'i': args->in_file = x_strcopy(optarg); break;
            case 'o': args->out_file = x_strcopy(optarg); break;
            case 'e': args->error_file = x_strcopy(optarg); break;
            case 'x': args->max_reads = strtoull(optarg, NULL, 10); break;
            case 'P': args->phred_value = atoi(optarg); break;
            case 'r': args->compress_ratio = atoi(optarg); break;
            case 'm': args->max_length = atoi(optarg)*1024*1024; break;
            case 'p': args->pair_check = atoi(optarg); break;
            case 't': args->n_thread = atoi(optarg); break;
            case '?': fprintf(stderr, "[SysError:args_parse:017] failed to parse the parameter of %s!\n\n", optarg); args->help = 1;
        }
    }

    /* do not check the pair_marker in single-end fastq file */
    if (get_file_number(args->in_file) < 2) 
        args->pair_check = 0;

    /* the pari marker only could be 0 or 1 */
    if (args->pair_check != 0 && args->pair_check != 1) {
        fprintf(stderr, "[SysError:args_parse:018] the pair_check parameter must be 0 or 1!\n");
        args->help = 1;
    }

    /* the phred value only could be 33 or 64 */
    if (args->phred_value != 33 && args->phred_value != 64) {
        fprintf(stderr, "[SysError:args_parse:025] the phred_value parameter must be 33 or 64!\n");
        args->help = 1;
    }

    /* required parameters */
    if (args->help||!args->in_file||!args->out_file||!args->error_file||!args->max_reads||!args->phred_value)
        show_usage();

    return args;
}
