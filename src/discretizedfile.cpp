#include <cstdint>
#include <cstdlib>
#include "discretizedfile.h"


DiscretizedFileInfo::DiscretizedFileInfo(int d, int o, int v) :
        discretizations(d),
        objectCount(o),
        variableCount(v) {}

DiscretizedFile::DiscretizedFile(DiscretizedFileInfo dfi) : info(dfi) {}
DiscretizedFile::~DiscretizedFile() {}

void DiscretizedFile::allocate() {
    this->data = new int32_t[this->info.discretizations * this->info.objectCount * this->info.variableCount];
    this->decision = new int[this->info.objectCount];
}

int32_t * DiscretizedFile::getVD(int v, int d) {
    std::size_t offset  = this->info.objectCount;
           offset *= (v * this->info.discretizations + d);
    return this->data + offset;
}

int DiscretizedFile::c1() {
    int c1 = 0;
    for (int i = 0; i < this->info.objectCount; ++i)
        if (this->decision[i] == 1)
            c1++;
    return c1;
}

int DiscretizedFile::c0() {
    int c0 = 0;
    for (int i = 0; i < this->info.objectCount; ++i)
        if (this->decision[i] == 0)
            c0++;
    return c0;
}
