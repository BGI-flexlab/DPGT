#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <cstdio>
#include "VariantSiteFinder.h"
#include "jni_md.h"
#include "variant_site_finder.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

JNIEXPORT void JNICALL Java_org_bgi_flexlab_dpgt_jointcalling_VariantSiteFinder_FindVariantSite
    (JNIEnv * env, jobject java_this, jobjectArray vcfpaths, jstring outpath, jstring chrom, jlong start, jlong end)
{
    int vcf_count = env->GetArrayLength(vcfpaths);
    std::vector<std::string> vcfpaths_vec(vcf_count, "");
    for (int i = 0; i < vcf_count; ++i) {
        jstring vcfpath = (jstring) (env->GetObjectArrayElement(vcfpaths, i));
        const char *vcfpath_cstr = env->GetStringUTFChars(vcfpath, NULL);
        vcfpaths_vec[i] = vcfpath_cstr;
        env->ReleaseStringUTFChars(vcfpath, vcfpath_cstr);
    }

    const char *chrom_cstr = env->GetStringUTFChars(chrom, NULL);
    
    std::vector<uint8_t> site_bitset_as_bytes = FindVariantSite(
        vcfpaths_vec, chrom_cstr, start, end).toUint8Vector();
    
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
