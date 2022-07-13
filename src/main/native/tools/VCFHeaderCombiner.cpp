#include <vector>
#include <string>
#include "jni_md.h"
#include "tools/VCFHeaderCombiner.hpp"
#include "vcf/multi_vcf_reader.hpp"
#include "vcf/vcf_writer.hpp"


JNIEXPORT void JNICALL Java_org_bgi_flexlab_dpgt_jointcalling_VCFHeaderCombiner_Combine
    (JNIEnv *env, jobject java_this, jobjectArray vcfpaths, jstring outpath)
{
    int vcf_count = env->GetArrayLength(vcfpaths);
    std::vector<std::string> vcfpaths_vec(vcf_count, "");
    for (int i = 0; i < vcf_count; ++i) {
        jstring vcfpath = (jstring) (env->GetObjectArrayElement(vcfpaths, i));
        const char *vcfpath_cstr = env->GetStringUTFChars(vcfpath, NULL);
        vcfpaths_vec[i] = vcfpath_cstr;
        env->ReleaseStringUTFChars(vcfpath, vcfpath_cstr);
    }

    const char *outpath_cstr = env->GetStringUTFChars(outpath, NULL);

    MultiVcfReader reader{vcfpaths_vec, false};
    bcf_hdr_t *header = reader.header();
    VcfWriter writer{outpath_cstr};
    writer.writeHeader(header);
    
    env->ReleaseStringUTFChars(outpath, outpath_cstr);
}
