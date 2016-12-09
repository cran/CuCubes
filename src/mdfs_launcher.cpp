#include "discretize.h"
#include "mdfs_common.h"
#include "mdfs_scalar.h"
#include "avxmdfs.h"
#include "avx2mdfs.h"

extern "C"
void CuCubes(MDFSAccelerationType *acceleration_type,
             MDFSOutputType *out_type,
                  int *n,                // obiekty
                  int *k,                // zmienne
                  int *dimension,        // wymiar [1..5]
                  int *divisions,        // podziały,
                  int *discretizations,  // liczba dyskretyzacji
                  int *seed,             // parametr do losowania dyskretyzacji
                  double *range,         // parametr dyskretyzacji,
                                         // {0-1, losowanie (1-range,1+range)}
                  double *pseudocount,   // suma wszystkich pseudozliczeń we wszystkich
                                         // kubełkach (vokselach)
                  int *reduce,           // 0 - max, 1 - avg
                  double *ig_thr,        // ig threshold value for the matching tuples output mode
                  int *interesting_vars, // interesting vars for the matching tuples output mode
                  int *interesting_vars_count,
                  double *data,          // długość n*k double, macierz - w formacie R, podajemy najpierw
                                         // wartości kolumny (czyli jednej zmiennej dla wszystkich obiektów)
                  int *decision,         // zmienna decyzyjna Boolowska - 0/1
                  double *IGmax)         // max IGs for each variable: array of length k
{
    int VAR = *k;
    int OBJ = *n;
    int DIM = *dimension;
    int DIV = *divisions;
    int DISC = *discretizations;
    int SEED = *seed;

    DataFile *df = new DataFile(DataFileInfo(OBJ, VAR), data, decision);

    DiscretizedFileInfo dfi(DISC, OBJ, VAR);
    DiscretizedFile *in = new DiscretizedFile(dfi);
    in->allocate();

    DiscretizationInfo di(SEED, DISC, DIV, (float)*range);
    discretizeFile(df, in, di);

    AlgInfo ai;
    ai.pseudo = (float) *pseudocount;
    ai.DIM = DIM;
    ai.DIV = DIV;
    ai.DISC = DISC;
    ai.rm = reduceMethod(*reduce);
    ai.ig_thr = *ig_thr;
    ai.interesting_vars = std::vector<int>(interesting_vars, interesting_vars + *interesting_vars_count);

    MDFSOutput out(*out_type, VAR);

    MDFSFunction mdfs = nullptr;
    switch (*acceleration_type) {
        case MDFSAccelerationType::Scalar:
            mdfs = ScalarMDFS;
            break;
        case MDFSAccelerationType::AVX:
            mdfs = AVXMdfs;
            break;
        case MDFSAccelerationType::AVX2:
            mdfs = AVX2Mdfs;
            break;
    }

    if (mdfs != nullptr) {
        mdfs(ai, in, out);
        switch (*out_type) {
            case MDFSOutputType::MaxIGs:
                out.CopyMaxIGsAsDouble(IGmax);
                break;
            case MDFSOutputType::MatchingTuples:
                out.Print();
                break;
        }
    }

    delete in;
    delete df;
}
