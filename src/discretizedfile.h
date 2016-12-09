#ifndef DISCRETIZEDFILE_H
#define DISCRETIZEDFILE_H

#include <cstdint>

class DiscretizedFileInfo {
public:
    DiscretizedFileInfo(int d, int o, int v);
    int discretizations;
    int objectCount;
    int variableCount;
};

// Stored in VDO way

class DiscretizedFile {

public:
    DiscretizedFile(DiscretizedFileInfo dfi);
    ~DiscretizedFile();
    DiscretizedFileInfo info;
    void allocate();
    int32_t * data;
    int * decision;
    int32_t * getVD(int v, int d);
    int c1();
    int c0();
};

#endif
