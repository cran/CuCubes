#ifndef STATS_VECTOR_H
#define STATS_VECTOR_H

template <typename T,
          T(*SET)(float),
          T(*MUL)(T a, T b),
          T(*ADD)(T a, T b),
          T (*FMA)(T a, T b, T c),
          T (*LOG)(T a)>
void vectorReduceCounter(int div, T *in, int dim, T *out, int reduced) {
    div += 1;
    int rstride = std::pow(div, (reduced - 1));
    int size = std::pow(div, dim);
    int v = 0;
    std::memset(out, 0, sizeof(T) * std::pow(div, (dim - 1)));
    for (int c = 0; c < size; c += rstride * div) {
        for (int s = 0; s < rstride; s++, v++) {
            for (int d = 0; d < div; d++) {
                out[v] = ADD(in[c + s + (d * rstride)], out[v]);
            }
        }
    }
}

template <typename T,
          T(*SET)(float),
          T(*MUL)(T a, T b),
          T(*ADD)(T a, T b),
          T (*FMA)(T a, T b, T c),
          T (*LOG)(T a)>
T vectorInformationGain(int counters, T *c0, T *c1) {
    T ig = SET(0.0f);
    for (int i = 0; i < counters; i++) {
        T c = ADD(c0[i], c1[i]);
        ig = FMA(c0[i], LOG(c0[i]/c), ig);
        ig = FMA(c1[i], LOG(c1[i]/c), ig);
    }
    return ig;
}

template <int VL,
          typename T>
T log_scalar(T in) {
    T out;
    float *fin = (float*)&in;
    float *fout = (float*)&out;
    for (int i = 0; i < VL; ++i) {
        fout[i] = std::log2(fin[i]);
    }
    return out;
}

#endif /* STATS_VECTOR_H */
