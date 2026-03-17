FastqCheck
========
Used for FASTQ format validation


__PROGRAM: fastq_check__<br>
__VERSION: 1.1.21__<br>
__PLATFORM: Linux__<br>
__COMPILER: gcc-4.8.5__<br>
__AUTHOR: xiaolong zhang__<br>
__EMAIL: xiaolongzhang2015@163.com__<br>
__DEPENDENCE__<br>
* __GNU make and gcc__<br>


Architecture
=========================
![fastq_check](images/fastq_check.png)


Building
=========================

```shell
cd FastqCheck
make
```

Usage
========================

```shell
Usage: fastq_check -i <fastq_list.txt> -o <output_file.xml> -e <error_file.err>
Options:
       -h|--help                 print help infomation

[Required]
       -i|--in_file        FILE  fastq file list, one sample path per line [.txt]
       -o|--out_file       FILE  output xml file of the given fastq files after checking [.xml]
       -e|--error_file     FILE  warning or error message about the given fastq files [.err]

[Optional]
       -r|--ratio          INT   the estimated compression ratio of the fastq file [default:10]
       -m|--max_length     INT   the maximum length (MB) allowed for one read sequence [default:50 (MB)]
       -t|--thread         INT   number of thread used when processing the fastq file [default:1]
       -p|--pair_check     INT   check whether all reads have the same pair marker (1->check; 0->ignore) [default:1]

```