#include "bloom_filter.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static void CalcBloomFilterParam(const uint64_t n, const double p, uint64_t *pm, uint64_t *pk)
{
    uint64_t m, k;

    m = (uint64_t)ceil(-1 * log(p) * n / 0.6185);
    m = (m - m % 64) + 64;

    k = (uint64_t)(0.6931 * m / n);
    k++;

    *pm = m;
    *pk = k;
}

BaseBloomFilter *InitBloomFilter(uint64_t dwMaxItems, double dProbFalse)
{
    BaseBloomFilter *pstBloomfilter;

    if (dwMaxItems == 0) {
        fprintf(stderr, "[BloomFilter:InitBloomFilter] invalid max item count 0\n");
        return NULL;
    }

    if ((dProbFalse <= 0) || (dProbFalse >= 1)) {
        fprintf(stderr, "[BloomFilter:InitBloomFilter] invalid false positive rate %.6f\n", dProbFalse);
        return NULL;
    }

    pstBloomfilter = (BaseBloomFilter *)calloc(1, sizeof(BaseBloomFilter));
    if (pstBloomfilter == NULL) {
        fprintf(stderr, "[BloomFilter:InitBloomFilter] failed to allocate bloom filter object\n");
        return NULL;
    }

    pstBloomfilter->dwMaxItems = dwMaxItems;
    pstBloomfilter->dProbFalse = dProbFalse;
    pstBloomfilter->dwSeed = 0;

    CalcBloomFilterParam(pstBloomfilter->dwMaxItems, pstBloomfilter->dProbFalse,
                         &pstBloomfilter->dwFilterBits, &pstBloomfilter->dwHashFuncs);

    pstBloomfilter->dwFilterSize = pstBloomfilter->dwFilterBits / BYTE_BITS;
    pstBloomfilter->pstFilter = (unsigned char *)malloc(pstBloomfilter->dwFilterSize);
    if (pstBloomfilter->pstFilter == NULL) {
        fprintf(stderr, "[BloomFilter:InitBloomFilter] failed to allocate filter buffer\n");
        free(pstBloomfilter);
        return NULL;
    }

    fprintf(stdout, "[*] BloomFilter initialization (n_item=%llu, false_positive=%.2e, n_bit=%llu, n_func=%llu, memory=%.2fMB)\n",
            pstBloomfilter->dwMaxItems, pstBloomfilter->dProbFalse, pstBloomfilter->dwFilterBits,
            pstBloomfilter->dwHashFuncs, (double)pstBloomfilter->dwFilterSize / 1024 / 1024);

    memset(pstBloomfilter->pstFilter, 0, pstBloomfilter->dwFilterSize);
    return pstBloomfilter;
}

int FreeBloomFilter(BaseBloomFilter *pstBloomfilter)
{
    if (pstBloomfilter == NULL)
        return -1;

    pstBloomfilter->dwCount = 0;

    free(pstBloomfilter->pstFilter);
    pstBloomfilter->pstFilter = NULL;
    free(pstBloomfilter);
    return 0;
}


int ResetBloomFilter(BaseBloomFilter *pstBloomfilter)
{
    if (pstBloomfilter == NULL)
        return -1;

    memset(pstBloomfilter->pstFilter, 0, pstBloomfilter->dwFilterSize);
    pstBloomfilter->dwCount = 0;
    return 0;
}

uint64_t MurmurHash2_x64(const void *key, const int len, const uint64_t seed)
{
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

    const uint64_t *data = (const uint64_t *)key;
    const uint64_t *end = data + (len / 8);

    while (data != end)
    {
        uint64_t k = *data++;

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const uint8_t *data2 = (const uint8_t *)data;

    switch (len & 7)
    {
    case 7: h ^= ((uint64_t)data2[6]) << 48;
    case 6: h ^= ((uint64_t)data2[5]) << 40;
    case 5: h ^= ((uint64_t)data2[4]) << 32;
    case 4: h ^= ((uint64_t)data2[3]) << 24;
    case 3: h ^= ((uint64_t)data2[2]) << 16;
    case 2: h ^= ((uint64_t)data2[1]) << 8;
    case 1: h ^= ((uint64_t)data2[0]);
            h *= m;
    }

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

int BloomFilterAdd(BaseBloomFilter *pstBloomfilter, const void *key, const int len)
{
    if ((pstBloomfilter == NULL) || (key == NULL) || (len <= 0))
        return -1;

    uint64_t hash1 = MurmurHash2_x64(key, len, pstBloomfilter->dwSeed);
    uint64_t hash2 = MurmurHash2_x64(key, len, MIX_UINT64(hash1));

    for (int i = 0; i < (int)pstBloomfilter->dwHashFuncs; i++)
    {
        const uint64_t pos = (hash1 + i * hash2) % pstBloomfilter->dwFilterBits;
        SETBIT(pstBloomfilter, pos);
    }

    pstBloomfilter->dwCount++;
    if (pstBloomfilter->dwCount <= pstBloomfilter->dwMaxItems)
        return 0;

    return 1;
}

int BloomFilterCheck(const BaseBloomFilter *pstBloomfilter, const void *key, const int len)
{
    if ((pstBloomfilter == NULL) || (key == NULL) || (len <= 0))
        return -1;

    uint64_t hash1 = MurmurHash2_x64(key, len, pstBloomfilter->dwSeed);
    uint64_t hash2 = MurmurHash2_x64(key, len, MIX_UINT64(hash1));

    for (int i = 0; i < (int)pstBloomfilter->dwHashFuncs; i++)
    {
        const uint64_t pos = (hash1 + i * hash2) % pstBloomfilter->dwFilterBits;
        if (GETBIT(pstBloomfilter, pos) == 0)
            return 1;
    }
    return 0;
}


int BloomFilterCheckAdd(BaseBloomFilter *pstBloomfilter, uint64_t hash1, uint64_t hash2)
{
    int is_existed = 1;

    for (int i = 0; i < (int)pstBloomfilter->dwHashFuncs; i++) {
        const uint64_t pos = (hash1 + i * hash2) % pstBloomfilter->dwFilterBits;
        if (GETBIT(pstBloomfilter, pos) != 0)
            continue;

        is_existed = 0;
        SETBIT(pstBloomfilter, pos);
    }

    if (is_existed)  /* the item may be existed in the bloomfilter */
        return 1;

    pstBloomfilter->dwCount++;
    if (pstBloomfilter->dwCount <= pstBloomfilter->dwMaxItems)
        return 0;

    return -1;  /* there is no space remain for more items */
}

