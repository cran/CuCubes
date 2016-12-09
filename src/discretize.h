#ifndef DISCRETIZE_H
#define DISCRETIZE_H

#include <cstddef>
#include <cstdint>
#include "discretizedfile.h"
#include "datafile.h"

class DiscretizationInfo {
public:
    DiscretizationInfo(uint32_t seed, int disc, int div, float range);
    uint32_t seed;
    int disc;
    int div;
    float range;
};

void discretizeFile(DataFile *in,
                    DiscretizedFile *out,
                    DiscretizationInfo info);

#endif
