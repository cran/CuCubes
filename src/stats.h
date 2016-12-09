#ifndef STATS_H
#define STATS_H

void reduceCounter(int div, float *in, int dim, float *out, int stride);

float informationGain(int counters, float *c0, float *c1);

#endif
