#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <random>
#include "discretize.h"

DiscretizationInfo::DiscretizationInfo(uint32_t seed, int disc, int div, float range) :
        seed(seed), disc(disc), div(div), range(range) {}

void discretize(uint32_t seed,
                uint32_t disc,
                uint32_t var,
                std::size_t div,
                std::size_t length,
                float *in_data,
                const std::vector<float>& sorted_in_data,
                int32_t *out_data,
                float range) {
    float* thr = new float[div];
    {
        float sum = 0.0f;
        {
            std::mt19937 seedGen0(seed);
            std::mt19937 seedGen1(seedGen0() ^ disc);
            std::mt19937 gen(seedGen1() ^ var);
            std::uniform_real_distribution<double> dis(1.0 - range, 1.0 + range);

            for (std::size_t d = 0; d < div; d++) {
                thr[d] = dis(gen);
                sum += thr[d];
            }

            sum += dis(gen);
        }

        std::size_t done = 0;

        for (std::size_t d = 0; d < div; d++) {
            done += std::lround(thr[d]/sum * length);
            if (done >= length) done = length-1;
            thr[d] = sorted_in_data[done];
        }
    }

    for (std::size_t i = 0; i < length; i++) {
        out_data[i] = 0;
        for (std::size_t d = 0; d < div; d++) {
            out_data[i] += in_data[i] > thr[d];
        }
    }

    delete[] thr;
}

void discretizeVar(DataFile *in,
                   DiscretizedFile *out,
                   int var,
                   DiscretizationInfo info) {
    float *in_data = in->getV(var);
    std::vector<float> sorted_in_data(in_data, in_data + in->info.objectCount);
    std::sort(sorted_in_data.begin(), sorted_in_data.end());
    for (int d = 0; d < info.disc; d ++) {
        discretize(info.seed, d, var, info.div, in->info.objectCount, in_data, sorted_in_data, out->getVD(var, d), info.range);
    }
}

void discretizeFile(DataFile *in,
                    DiscretizedFile *out,
                    DiscretizationInfo info) {
    memcpy(out->decision, in->decision, sizeof(int) * in->info.objectCount);
    for (int v = 0; v < in->info.variableCount; v++) {
        discretizeVar(in, out, v, info);
    }
}
