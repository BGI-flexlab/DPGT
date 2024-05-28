#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <cstdio>
#include "VariantSiteFinder.h"
#include "jni.h"
#include "jni_md.h"
#include "variant_site_finder.hpp"
#include "common/jni_utils.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

JNIEXPORT void JNICALL Java_org_bgi_flexlab_dpgt_jointcalling_VariantSiteFinder_FindVariantSite
    (JNIEnv * env, jobject java_this, jobjectArray vcfpaths, jobjectArray indexpaths,
    jint use_lix, jstring outpath, jstring chrom, jlong start, jlong end)
{
    int vcf_count = env->GetArrayLength(vcfpaths);
    std::vector<std::string> vcfpaths_vec(vcf_count, "");
    JavaStringArrayToCppVec(env, vcfpaths, vcfpaths_vec);

    int index_count = env->GetArrayLength(indexpaths);
    std::vector<std::string> indexpaths_vec(index_count, "");
    JavaStringArrayToCppVec(env, indexpaths, indexpaths_vec);

    const char *chrom_cstr = env->GetStringUTFChars(chrom, NULL);
    
    std::vector<uint8_t> site_bitset_as_bytes = FindVariantSite(
        vcfpaths_vec, indexpaths_vec, use_lix, chrom_cstr, start, end).toUint8Vector();
    
    // serialize to disk
    const char *outpath_cstr = env->GetStringUTFChars(outpath, NULL);
    FILE *outfp = fopen(outpath_cstr, "wb");
    if (outfp != NULL) {
        const int32_t size = site_bitset_as_bytes.size();
        fwrite(&size, sizeof(int32_t), 1, outfp);  // first 4-bytes is the size of bitset bytes
        fwrite(site_bitset_as_bytes.data(), 1, size, outfp);
        fclose(outfp);
    } else {
        std::string error_msg = "failed to open: " + std::string(outpath_cstr);
        perror(error_msg.c_str());
        std::exit(1);
    }
    
    env->ReleaseStringUTFChars(outpath, outpath_cstr);
    env->ReleaseStringUTFChars(chrom, chrom_cstr);
}
