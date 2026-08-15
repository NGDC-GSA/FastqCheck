#ifndef MICRO_BLOOMFILTER_H
#define MICRO_BLOOMFILTER_H

/**
 *  BloomFilter_x64实现: https://github.com/upbit/bloomfilter
 *
 *  仿照Cassandra中的BloomFilter实现，Hash选用MurmurHash2，通过双重散列公式生成散列函数
 *    Hash(key, i) = (H1(key) + i * H2(key)) % m
 *
 *  2012.12    完成初始版本
 *  2013.4.10  增加k/m的动态计算功能，参考：http://hur.st/bloomfilter
 *  2022.8.5   增加数据结构对大样本的支持(uint32_t -> uint64_t)
**/

#include <stdint.h>

#define BLOOMFILTER_VERSION "1.3"
#define MAGIC_CODE          (0x01464C42)

/**
 *  BloomFilter使用例子：
 *  static BaseBloomFilter *stBloomFilter = NULL;
 *
 *  初始化BloomFilter(最大100000元素，不超过0.00001的错误率)：
 *      stBloomFilter = InitBloomFilter(100000, 0.00001);
 *  重置BloomFilter：
 *      ResetBloomFilter(stBloomFilter);
 *  释放BloomFilter:
 *      FreeBloomFilter(stBloomFilter);
 *
 *  向BloomFilter中新增一个数值（0-正常，1-加入数值过多）：
 *      uint32_t dwValue;
 *      iRet = BloomFilterAdd(stBloomFilter, &dwValue, sizeof(uint32_t));
 *  检查数值是否在BloomFilter内（0-存在，1-不存在）：
 *      iRet = BloomFilterCheck(stBloomFilter, &dwValue, sizeof(uint32_t));
 *
 *  (1.1新增) 将生成好的BloomFilter写入文件:
 *      iRet = SaveBloomFilterToFile(&stBloomFilter, "dump.bin")
 *  (1.1新增) 从文件读取生成好的BloomFilter:
 *      iRet = LoadBloomFilterFromFile(&stBloomFilter, "dump.bin")
 *  (1.2新增) (xiaolongzhang: 20220805):
 *      (1) modify the code to allow more big memory set other than uint32_t
 *      (2) delete the function of "SaveBloomFilterToFile" and "LoadBloomFilterFromFile"
 *  (1.3新增) (xiaolongzhang: 20260518):
 *      (1) split the source code and header file to support subsequent parallel computing.
 *  (1.4新增) (xiaolongzhang: 20260816):
 *      (1) add function of BloomFilterCheckAdd for checking and adding item at the sametime
**/


#define BYTE_BITS           (8)


#pragma pack(1)

/*! @typedef BaseBloomFilter
  @abstract structure for BloomFilter object
  @field  dwMaxItems      maximum number of items expected in the filter
  @field  dProbFalse      expected false positive rate
  @field  dwFilterBits    total number of bits allocated for the filter
  @field  dwHashFuncs     number of hash functions used by the filter
  @field  dwSeed          seed used by MurmurHash2
  @field  dwCount         number of items added into the filter
  @field  dwFilterSize    size of filter memory in bytes
  @field  pstFilter       pointer to the bit array of BloomFilter
 */
typedef struct {
    uint64_t dwMaxItems;
    double dProbFalse;
    uint64_t dwFilterBits;
    uint64_t dwHashFuncs;

    uint64_t dwSeed;
    uint64_t dwCount;

    uint64_t dwFilterSize;
    unsigned char *pstFilter;
} BaseBloomFilter;

#pragma pack()


#define MIX_UINT64(v)       ((uint64_t)(((v) >> 32) ^ (v)))
#define SETBIT(filter, n)   ((filter)->pstFilter[(n) / BYTE_BITS] |= (1U << ((n) % BYTE_BITS)))
#define GETBIT(filter, n)   ((filter)->pstFilter[(n) / BYTE_BITS] & (1U << ((n) % BYTE_BITS)))


/*! @function: initialize BloomFilter by expected item count and false positive rate
  @param  dwMaxItems      maximum number of items expected in the filter
  @param  dProbFalse      expected false positive rate, should be between 0 and 1
  @return                 pointer to initialized BaseBloomFilter object; NULL->error
 */
BaseBloomFilter *InitBloomFilter(uint64_t dwMaxItems, double dProbFalse);


/*! @function: free the memory allocated for BloomFilter
  @param  pstBloomfilter  pointer to BaseBloomFilter object
  @return                 0->OK; negative value->error
 */
int FreeBloomFilter(BaseBloomFilter *pstBloomfilter);


/*! @function: reset BloomFilter and clear filter memory immediately
  @param  pstBloomfilter  pointer to BaseBloomFilter object
  @return                 0->OK; negative value->error
 */
int ResetBloomFilter(BaseBloomFilter *pstBloomfilter);


/*! @function: 64-bit MurmurHash2 implementation
  @param  key             pointer to input data
  @param  len             length of input data
  @param  seed            seed used for hash calculation
  @return                 64-bit hash value
 */
uint64_t MurmurHash2_x64(const void *key, int len, uint64_t seed);


/*! @function: add one item into BloomFilter
  @param  pstBloomfilter  pointer to BaseBloomFilter object
  @param  key             pointer to input data
  @param  len             length of input data
  @return                 0->OK; 1->item count exceeds limit; negative value->error
 */
int BloomFilterAdd(BaseBloomFilter *pstBloomfilter, const void *key, int len);


/*! @function: check whether one item may exist in BloomFilter
  @param  pstBloomfilter  pointer to BaseBloomFilter object
  @param  key             pointer to input data
  @param  len             length of input data
  @return                 0->possibly exists; 1->definitely not exists; negative value->error
 */
int BloomFilterCheck(const BaseBloomFilter *pstBloomfilter, const void *key, int len);


/*! @function: check whether one item may exist in BloomFilter and add it into the filter at the same time
  @param  pstBloomfilter  pointer to BaseBloomFilter object
  @param  hash1           the first hash value of the item (H1 in the double hashing formula)
  @param  hash2           the second hash value of the item (H2 in the double hashing formula)
  @return                 0->item added successfully; 1->item may already exist; negative->item count exceeds limit
 */
int BloomFilterCheckAdd(BaseBloomFilter *pstBloomfilter, uint64_t hash1, uint64_t hash2);


#endif
