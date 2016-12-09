#ifndef DATAFILE_H
#define DATAFILE_H

class DataFileInfo {

public:
    DataFileInfo(int o, int v);
    int objectCount;
    int variableCount;
};

//Stored in VO way

class DataFile {

public:
    DataFile(DataFileInfo dfi);
    DataFile(DataFileInfo dfi, double* data, int *decision);
    ~DataFile();
    DataFileInfo info;
    void allocate();
    float * data;
    int * decision;
    float *getV(int var);
};

#endif
