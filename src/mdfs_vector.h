#ifndef MDFS_VECTOR
#define MDFS_VECTOR

#include <algorithm>
#include <cstring>
#include <vector>
#include <numeric>
#include <cmath>

#include "mdfs_common.h"
#include "vec_stats.h"
#include "vec_discretizedfile.h"

#define CONTAINS(x, y) (std::find((x).begin(), (x).end(), (y)) != (x).end())

template<typename T,
	 T (*MUL) (T, T),
	 T (*ADD) (T, T)>
inline T my_fma(T a, T b, T c)
{
    return ADD(MUL(a, b), c);
}

template <int VL,
          typename T,
          T(*SET)(float),
          T(*MUL)(T a, T b),
          T(*ADD)(T a, T b),
          T(*SUB)(T a, T b),
          T (*FMA)(T a, T b, T c),
          T (*LOG)(T a),
          typename Td,
          Td(*SETd)(int a),
          Td(*MULd)(Td a, Td b),
          Td(*ADDd)(Td a, Td b)>
void vectorMdfs_scheme(AlgInfo ai,
                       VectorDiscretizedFile<VL> *in,
                       MDFSOutput &out,
                       int c0,
                       int c1)
{
    int cc = std::pow(ai.DIV + 1, ai.DIM);
    int cd = cc / (ai.DIV + 1);

    float sp0  = (float)c0 / (c1 + c0);
          sp0 *= ai.pseudo;
          sp0 /= cc;
    float sp1  = (float)c1 / (c1 + c0);
          sp1 *= ai.pseudo;
          sp1 /= cc;

    T p0 = SET(sp0);
    T p1 = SET(sp1);

    T* ig = (T*)_mm_malloc(sizeof(T) * ai.DISC/VL * ai.DIM, sizeof(T));
    float* dig = new float[ai.DIM];

    VarsTuple v(ai.DIM, in->info.variableCount);
    while (!v.done()) {
        std::list<int> current_interesting_vars;
        std::set_intersection(
            v.begin(), v.end(),
            ai.interesting_vars.begin(), ai.interesting_vars.end(),
            std::back_inserter(current_interesting_vars));
        if (!ai.interesting_vars.empty() && current_interesting_vars.empty()) {
            v.next();
            continue;
        }

        #pragma omp parallel for
        for (int d = 0; d < in->info.discretizations/VL; ++d) {
            T* counters = (T*)_mm_malloc(sizeof(T) * cc * 2, sizeof(T));
            T* reduced = (T*)_mm_malloc(sizeof(T) * cd * 2, sizeof(T));

            std::memset(counters, 0, sizeof(float) * VL * cc * 2);

            for (int o = 0; o < in->info.objectCount; ++o) {
                Td b = SETd(0);
                for (int vv = ai.DIM-1; vv >= 0; vv--) {
                    b = MULd(b, SETd(ai.DIV + 1));
                    b = ADDd(b, ((Td *)in->getVDpack(v.get(vv), d))[o]);
                }
                int dec = in->decision[o];
                int32_t * buckets = (int32_t*)&b;
                for (int bs = 0; bs < VL; bs++) {
                    ((float*)&(counters[dec * cc + buckets[bs]]))[bs] += 1.0f;
                }
            }

            for (int b = 0; b < cc; ++b) {
                counters[0 * cc + b] = ADD(counters[0 * cc + b], p0);
                counters[1 * cc + b] = ADD(counters[1 * cc + b], p1);
            }

            T ign = vectorInformationGain<T, SET, MUL, ADD, FMA, LOG>(cc, counters, counters+cc);

            for (int vv = 0; vv < ai.DIM; vv++) {
                vectorReduceCounter<T, SET, MUL, ADD, FMA, LOG>(ai.DIV, counters, ai.DIM, reduced, vv + 1);
                vectorReduceCounter<T, SET, MUL, ADD, FMA, LOG>(ai.DIV, counters+cc, ai.DIM, reduced+cd, vv + 1);
                T igg = vectorInformationGain<T, SET, MUL, ADD, FMA, LOG>(cd, reduced, reduced+cd);
                T igv = SUB(ign, igg);
                ig[vv * ai.DISC/VL + d] = igv;
            }

            _mm_free(counters);
            _mm_free(reduced);
        }

        switch (ai.rm) {
            case reduceMethod::RM_AVG:
                for (int vv = 0; vv < ai.DIM; vv++) {
                    T sum = SET(0.0f);
                    for (int i = 0; i < ai.DISC/VL; ++i) {
                        sum = ADD(sum, ig[vv * ai.DISC/VL + i]);
                    }
                    dig[vv] = 0.0f;
                    for (int i = 0; i < VL; i++) {
                        dig[vv] += ((float *)&sum)[i];
                    }
                    dig[vv] /= ai.DISC;
                }
                break;
            case reduceMethod::RM_MAX:
                for (int vv = 0; vv < ai.DIM; vv++) {
                    dig[vv] = 0.0f;
                    for (int i = 0; i < ai.DISC/VL; ++i) {
                        T tmp = ig[vv * ai.DISC/VL + i];
                        for (int ii = 0; ii < VL; ii++) {
                            dig[vv] = std::max(((float *)&tmp)[ii], dig[vv]);
                        }
                    }
                }
                break;
            default:
                break;
        }

        switch (out.type) {
            case MDFSOutputType::MaxIGs:
                for (int vv = 0; vv < ai.DIM; vv++) {
                    out.UpdateMaxIG(v.get(vv), dig[vv]);
                }
                break;
            case MDFSOutputType::MatchingTuples:
                for (int vv = 0; vv < ai.DIM; vv++) {
                    if (dig[vv] > ai.ig_thr && (current_interesting_vars.empty() || CONTAINS(current_interesting_vars, v.get(vv)))) {
                        out.AddTuple(v.get(vv), dig[vv], v);
                    }
                }
                break;
        }

        v.next();
    }

    _mm_free(ig);
    delete[] dig;
}

template <int VL,
          typename T,
          T(*SET)(float),
          T(*MUL)(T a, T b),
          T(*ADD)(T a, T b),
          T(*SUB)(T a, T b),
          T (*FMA)(T a, T b, T c),
          T (*LOG)(T a),
          typename Td,
          Td(*SETd)(int a),
          Td(*MULd)(Td a, Td b),
          Td(*ADDd)(Td a, Td b)>
void vectorMdfs(AlgInfo ai,
                VectorDiscretizedFile<VL> *in,
                MDFSOutput &out)
{
    vectorMdfs_scheme<VL, T, SET, MUL, ADD, SUB, FMA, LOG, Td, SETd, MULd, ADDd>(ai, in, out, in->c0(), in->c1());
}

#endif
