#include "avx2mdfs.h"
#include "mdfs_vector.h"

void AVX2Mdfs(AlgInfo ai, DiscretizedFile *in, MDFSOutput &out)
{
    VectorDiscretizedFile<8> *vin = new VectorDiscretizedFile<8>(in);
    vectorMdfs<8,
               __m256,
               _mm256_set1_ps,
               _mm256_mul_ps,
               _mm256_add_ps,
               _mm256_sub_ps,
               _mm256_fmadd_ps,
               log_scalar<8, __m256>,//_mm256_log_ps,
               __m256i,
               _mm256_set1_epi32,
               _mm256_mullo_epi32,
               _mm256_add_epi32>(ai, vin, out);

    delete vin;
}
