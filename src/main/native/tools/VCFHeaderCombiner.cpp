#include <vector>
#include <string>
#include "jni_md.h"
#include "tools/VCFHeaderCombiner.hpp"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "vcf/vcf_writer.hpp"
#include "vcf/hts_vcf_utils.hpp"


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

    const int chunk_size = 100;
    std::vector<bcf_hdr_t *> mheaders; // merged header for each chunk
    std::vector<bcf_hdr_t *> headers;
    char fmode[5];
    strcpy(fmode, "r");
    int n = 1;
    for (int i = 0; i < (int)vcfpaths_vec.size(); ++i) {
        vcf_open_mode(fmode+1, vcfpaths_vec[i].c_str(), NULL);
        htsFile *fp = hts_open(vcfpaths_vec[i].c_str(), fmode);
        if ( ! fp ) {
            std::cerr << "[VCFHeaderCombiner_Combine] Error! Fail to open "
                << vcfpaths_vec[i] << std::endl;
            std::exit(1);
        }
        bcf_hdr_t *hdr = bcf_hdr_read(fp);
        headers.push_back(hdr);
        ++n;
        hts_close(fp);
        if (n >= chunk_size) {
            bcf_hdr_t *mhdr = bcf_hdr_merge_add_samples(headers);
            for (auto it: headers) {
                bcf_hdr_destroy(it);
            }
            headers.clear();
            mheaders.push_back(mhdr);
            n = 1;
        }
    }

    if (headers.size() > 0) {
        bcf_hdr_t *mhdr = bcf_hdr_merge_add_samples(headers);
        mheaders.push_back(mhdr);
    }
    bcf_hdr_t *result = bcf_hdr_merge_add_samples(mheaders);
    VcfWriter writer{outpath_cstr};
    writer.writeHeader(result);

    for (auto it: headers) {
        bcf_hdr_destroy(it);
    }
    headers.clear();

    for (auto it: mheaders) {
        bcf_hdr_destroy(it);
    }
    mheaders.clear();

    bcf_hdr_destroy(result);
    
    env->ReleaseStringUTFChars(outpath, outpath_cstr);
}
