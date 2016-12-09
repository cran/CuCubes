
#include <cstdlib>
#include <cstring>
#include "datafile.h"

DataFileInfo::DataFileInfo(int o, int v) : objectCount(o), variableCount(v) {}


DataFile::DataFile(DataFileInfo dfi) : info(dfi) {}

DataFile::DataFile(DataFileInfo dfi, double *data, int *decision) : info(dfi) {
    this->allocate();
    for (int i = 0; i < dfi.objectCount * dfi.variableCount; i++)
        this->data[i] = (float)data[i];
    std::memcpy(this->decision, decision, sizeof(int) * dfi.objectCount);
}

void DataFile::allocate() {
    this->data = new float[this->info.objectCount * this->info.variableCount];
    this->decision = new int[this->info.objectCount];
}

float * DataFile::getV(int v) {
    std::size_t offset = v * this->info.objectCount;
    return this->data + offset;
}

DataFile::~DataFile() {}

