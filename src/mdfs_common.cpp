#include <R.h>

#include "mdfs_common.h"


VarsTuple::VarsTuple(int dim, int var_count) : dim(dim), var_count(var_count), v(dim+1) {
    v[0] = 0;
    for (int d = 1; d <= dim; d++)
        v[d] = d - 1;
}

void VarsTuple::next() {
    int d;
    for (d = dim; d >= 0; d--) {
        v[d] ++;
        if (v[d] < var_count-(dim-d))
            break;
    }
    for (d++; d <= dim; d++) {
        v[d] = v[d-1]+1;
    }
}

bool VarsTuple::done() {
    return (v[0] > 0);
}

int VarsTuple::get(int i) {
    return v[i + 1];
}

std::vector<int>::const_iterator VarsTuple::begin() const {
    return v.begin() + 1;
}

std::vector<int>::const_iterator VarsTuple::end() const {
    return v.end();
}


MDFSTuple::MDFSTuple(int i, float ig, std::vector<int>&& v): i(i), ig(ig), v(v) {
}

void MDFSTuple::Print() {
    Rprintf("%d:%f:%d", i, ig, v[0]);
    for (auto u = v.begin() + 1; u != v.end(); u++) {
        Rprintf(",%d", *u);
    }
    Rprintf("\n");
}


MDFSOutput::MDFSOutput(MDFSOutputType type, int var_count): type(type) {
    switch(type) {
        case MDFSOutputType::MaxIGs:
            max_igs = new std::vector<float>(var_count);
            break;
        case MDFSOutputType::MatchingTuples:
            tuples = new std::list<MDFSTuple>();
            break;
   }
}

MDFSOutput::~MDFSOutput() {
    switch(type) {
        case MDFSOutputType::MaxIGs:
            delete max_igs;
            break;
        case MDFSOutputType::MatchingTuples:
            delete tuples;
            break;
   }
}

void MDFSOutput::Print() {
    switch(type) {
        case MDFSOutputType::MaxIGs:
            Rprintf("%f", (*max_igs)[0]);
            for (auto v = max_igs->begin() + 1; v != max_igs->end(); v++) {
                Rprintf("\t%f", *v);
            }
            Rprintf("\n");
            break;
        case MDFSOutputType::MatchingTuples:
            for (auto v = tuples->begin(); v != tuples->end(); v++) {
                v->Print();
            }
            break;
   }
}

void MDFSOutput::UpdateMaxIG(int i, float v) {
    (*max_igs)[i] = std::max((*max_igs)[i], v);
}

void MDFSOutput::CopyMaxIGsAsDouble(double* copy) {
    std::copy(max_igs->begin(), max_igs->end(), copy);
}

void MDFSOutput::AddTuple(int i, float ig, const VarsTuple &vt) {
    tuples->emplace_back(i, ig, std::vector<int>(vt.begin(), vt.end()));
}
