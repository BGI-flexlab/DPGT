#include <cstddef>
#include <cstdlib>
#include <vector>
#include "VariantSiteFinder.h"
#include "jni_md.h"
#include "variant_site_finder.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

JNIEXPORT jbyteArray JNICALL Java_org_bgi_flexlab_dpgt_jointcalling_VariantSiteFinder_FindVariantSite
    (JNIEnv * env, jobject java_this, jobjectArray vcfpaths, jstring chrom, jlong start, jlong end)
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
    
    // convert to java byte(signed char) type
    jbyte *bytes = new jbyte[site_bitset_as_bytes.size()];
    for (size_t i = 0; i < site_bitset_as_bytes.size(); ++i) {
        bytes[i] = site_bitset_as_bytes[i];
    }
    
    jbyteArray result = env->NewByteArray(site_bitset_as_bytes.size());
    env->SetByteArrayRegion(result, 0, site_bitset_as_bytes.size(),
        bytes);
    
    delete [] bytes;

    env->ReleaseStringUTFChars(chrom, chrom_cstr);

    return result;
}
