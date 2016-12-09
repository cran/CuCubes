#include "avxmdfs.h"
#include "mdfs_vector.h"

void AVXMdfs(AlgInfo ai, DiscretizedFile *in, MDFSOutput &out)
{
    VectorDiscretizedFile<4> *vin = new VectorDiscretizedFile<4>(in);
    vectorMdfs<4,
               __m128,
               _mm_set1_ps,
               _mm_mul_ps,
               _mm_add_ps,
               _mm_sub_ps,
               my_fma<__m128, _mm_mul_ps, _mm_add_ps>,
               log_scalar<4, __m128>,//_mm128_log_ps,
               __m128i,
               _mm_set1_epi32,
               _mm_mullo_epi32,
               _mm_add_epi32>(ai, vin, out);
    delete vin;
}
