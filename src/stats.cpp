#include <cmath>
#include <cstring>

void reduceCounter(int div, float *in, int dim, float *out, int reduced) {
    div += 1;
    int rstride = std::pow(div, (reduced - 1));
    int size = std::pow(div, dim);
    int v = 0;
    std::memset(out, 0, sizeof(float) * std::pow(div, (dim - 1)));
    for (int c = 0; c < size; c += rstride * div) {
        for (int s = 0; s < rstride; s++, v++) {
            for (int d = 0; d < div; d++) {
                out[v] += in[c + s + (d * rstride)];
            }
        }
    }
}

float informationGain(int counters, float *c0, float *c1) {
    float ig = 0.0f;
    for (int i = 0; i < counters; i++) {
        float c = c0[i] + c1[i];
        if (c0[i] != 0.0f) ig += (c0[i]) * std::log2(c0[i]/c);
        if (c1[i] != 0.0f) ig += (c1[i]) * std::log2(c1[i]/c);
    }
    return ig;
}
