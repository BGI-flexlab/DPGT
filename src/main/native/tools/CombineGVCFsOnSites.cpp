#include <cstdint>
#include <vector>
#include <string>
#include "jni_md.h"
#include "tools/CombineGVCFsOnSites.hpp"
#include "tools/combine_gvcfs_on_sites.hpp"
#include "common/jni_utils.hpp"
#include "vcf/multi_vcf_reader.hpp"

JNIEXPORT void JNICALL Java_org_bgi_flexlab_dpgt_jointcalling_CombineGVCFsOnSites_Combine
    (JNIEnv *env, jobject java_this, jobjectArray vcfpaths, jobjectArray indexpaths, jint use_lix,
    jstring refpath, jstring outpath, jbyteArray bytes, jstring chrom, jlong start, jlong end)
{
    int vcf_count = env->GetArrayLength(vcfpaths);
    std::vector<std::string> vcfpaths_vec(vcf_count, "");
    JavaStringArrayToCppVec(env, vcfpaths, vcfpaths_vec);

    int index_count = env->GetArrayLength(indexpaths);
    std::vector<std::string> indexpaths_vec(index_count, "");
    JavaStringArrayToCppVec(env, indexpaths, indexpaths_vec);

    const char *chrom_cstr = env->GetStringUTFChars(chrom, NULL);
    const char *refpath_cstr = env->GetStringUTFChars(refpath, NULL);
    const char *outpath_cstr = env->GetStringUTFChars(outpath, NULL);

    MultiVcfReader reader{vcfpaths_vec, indexpaths_vec, true, use_lix == 0 ? false : true};
    reader.Querys(chrom_cstr, start, end);
    CachedRefSeq cached_ref_seq(refpath_cstr);
    VcfWriter vcf_writer(outpath_cstr);

    int byte_count = env->GetArrayLength(bytes);
    std::vector<uint8_t> byte_vec(byte_count, 0);
    int8_t *bytes_array = env->GetByteArrayElements(bytes, NULL);
    for (int i = 0; i < byte_count; ++i) {
        byte_vec[i] = (uint8_t)bytes_array[i];
    }
    env->ReleaseByteArrayElements(bytes, bytes_array, 0);
    VariantSiteSet vs_set(start, end);
    vs_set.fromUnit8Vector(byte_vec);

    CombineGVCFsOnSites combiner{&reader, &cached_ref_seq, &vs_set,
        chrom_cstr, &vcf_writer};
    
    combiner.run();

    env->ReleaseStringUTFChars(outpath, outpath_cstr);
    env->ReleaseStringUTFChars(refpath, refpath_cstr);
    env->ReleaseStringUTFChars(chrom, chrom_cstr);
}
