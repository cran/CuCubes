#ifndef MDFS_COMMON_H
#define MDFS_COMMON_H

#include <vector>
#include <list>

#include "discretizedfile.h"

enum class MDFSAccelerationType { Scalar, AVX, AVX2 };

enum class reduceMethod { RM_MAX, RM_AVG };

enum class MDFSOutputType { MaxIGs, MatchingTuples };

struct AlgInfo {
    int DIM;
    int DIV;
    int DISC;
    float pseudo;
    reduceMethod rm;
    float ig_thr;
    std::vector<int> interesting_vars;
};

class VarsTuple {
    private:
        const int dim;
        const int var_count;
        std::vector<int> v;
    public:
        VarsTuple(int dim, int var_count);
        void next();
        bool done();
        int get(int i);
        std::vector<int>::const_iterator begin() const;
        std::vector<int>::const_iterator end() const;
};

class MDFSTuple {
    int i;
    float ig;
    std::vector<int> v;
public:
    MDFSTuple(int i, float ig, std::vector<int>&& v);
    void Print();
};

class MDFSOutput {
    union {
        std::vector<float>* max_igs;
        std::list<MDFSTuple>* tuples;
    };
public:
    const MDFSOutputType type;
    MDFSOutput(MDFSOutputType type, int var_count);
    ~MDFSOutput();
    void Print();
    void UpdateMaxIG(int i, float v);
    void CopyMaxIGsAsDouble(double* copy);
    void AddTuple(int i, float ig, const VarsTuple &vt);
};

using MDFSFunction = void (*) (AlgInfo, DiscretizedFile*, MDFSOutput&);

#endif
