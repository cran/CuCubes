#include <algorithm>
#include <cstring>
#include <vector>
#include <numeric>
#include <cmath>

#include "mdfs_scalar.h"
#include "stats.h"

#define CONTAINS(x, y) (std::find((x).begin(), (x).end(), (y)) != (x).end())

void mdfs_scheme(AlgInfo ai,
                 DiscretizedFile *in,
                 MDFSOutput &out,
                 int c0,
                 int c1) {
    int cc = std::pow(ai.DIV + 1, ai.DIM);
    int cd = cc / (ai.DIV + 1);

    float p0  = (float)c0 / (c1 + c0);
          p0 *= ai.pseudo;
          p0 /= cc;
    float p1  = (float)c1 / (c1 + c0);
          p1 *= ai.pseudo;
          p1 /= cc;

    float* ig = new float[ai.DIM * ai.DISC];
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
        for (int d = 0; d < in->info.discretizations; ++d) {
            float* counters = new float[2 * cc];
            float* reduced = new float[2 * cd];

            std::memset(counters, 0, sizeof(float) * cc * 2);

            for (int o = 0; o < in->info.objectCount; ++o) {
                int b = 0;
                for (int vv = ai.DIM-1; vv >= 0; vv--) {
                    b *= ai.DIV + 1;
                    b += in->getVD(v.get(vv), d)[o];
                }

                int dec = in->decision[o];
                counters[dec * cc + b] += 1.0f;
            }

            for (int b = 0; b < cc; ++b) {
                counters[0 * cc + b] += p0;
                counters[1 * cc + b] += p1;
            }

            float ign = informationGain(cc, counters, counters+cc);

            for (int vv = 0; vv < ai.DIM; vv++) {
                reduceCounter(ai.DIV, counters, ai.DIM, reduced, vv + 1);
                reduceCounter(ai.DIV, counters+cc, ai.DIM, reduced+cd, vv + 1);
                float igg = informationGain(cd, reduced, reduced+cd);
                ig[vv * ai.DISC + d] = ign - igg;
            }

            delete[] counters;
            delete[] reduced;
        }

        switch (ai.rm) {
            case reduceMethod::RM_AVG:
                for (int vv = 0; vv < ai.DIM; vv++) {
                    dig[vv]  = std::accumulate(ig + vv * ai.DISC, ig + vv * ai.DISC + ai.DISC, 0.0f) / ai.DISC;
                }
                break;
            case reduceMethod::RM_MAX:
                for (int vv = 0; vv < ai.DIM; vv++)
                    dig[vv] = *std::max_element(ig + vv * ai.DISC, ig + vv * ai.DISC + ai.DISC);
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

    delete[] ig;
    delete[] dig;
}

void ScalarMDFS(AlgInfo ai,
                DiscretizedFile *in,
                MDFSOutput &out) {
    mdfs_scheme(ai, in, out, in->c0(), in->c1());
}
