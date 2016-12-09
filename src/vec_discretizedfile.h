#ifndef VECTORDISCRETIZEDFILE_H
#define VECTORDISCRETIZEDFILE_H

#include <cstdint>
#include <immintrin.h>
#include "discretizedfile.h"


// Stored in VDO way

template <int VL>
class VectorDiscretizedFile {
public:
    VectorDiscretizedFile(DiscretizedFile * df);
    ~VectorDiscretizedFile();
    DiscretizedFileInfo info;
    void allocate();
    int32_t * data;
    int * decision;
    int32_t * getVDpack(int v, int d);
    int c1();
    int c0();
};

template <int VL>
VectorDiscretizedFile<VL>::VectorDiscretizedFile(DiscretizedFile * df) : info(df->info) {
    this->allocate();
    int cnt = 0;
    for (int v = 0; v < this->info.variableCount; ++v) {
        for (int dpack = 0; dpack < this->info.discretizations/VL; ++dpack) {
            for (int o = 0; o < this->info.objectCount; ++o) {
                for (int d = dpack * VL; d < dpack * VL + VL; d++) {
                   this->data[cnt] = df->getVD(v, d)[o];
                   cnt++;
                }
            }
        }
    }
    std::memcpy(this->decision, df->decision, sizeof(int) * this->info.objectCount);
}

template <int VL>
VectorDiscretizedFile<VL>::~VectorDiscretizedFile() {
    _mm_free(this->data);
    delete [] this->decision;
}

template <int VL>
void VectorDiscretizedFile<VL>::allocate() {
    this->data = (int32_t *) _mm_malloc(this->info.discretizations * this->info.objectCount * this->info.variableCount * sizeof(int32_t), VL * sizeof(int32_t));
//    this->data = new int32_t[this->info.discretizations * this->info.objectCount * this->info.variableCount];
    this->decision = new int[this->info.objectCount];
}

template <int VL>
int32_t * VectorDiscretizedFile<VL>::getVDpack(int v, int dpack) {
    std::size_t offset  = this->info.objectCount;
                offset *= (v * this->info.discretizations + (dpack * VL));
    return this->data + offset;
}

template <int VL>
int VectorDiscretizedFile<VL>::c1() {
    int c1 = 0;
    for (int i = 0; i < this->info.objectCount; ++i)
        if (this->decision[i] == 1)
            c1++;
    return c1;
}

template <int VL>
int VectorDiscretizedFile<VL>::c0() {
    int c0 = 0;
    for (int i = 0; i < this->info.objectCount; ++i)
        if (this->decision[i] == 0)
            c0++;
    return c0;
}
#endif
