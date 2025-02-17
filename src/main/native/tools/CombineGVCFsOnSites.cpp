/**
This file is part of DPGT.
Copyright (C) 2022 BGI.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// License End
#include <cstdint>
#include <vector>
#include <string>
#include "jni_md.h"
#include "tools/CombineGVCFsOnSites.hpp"
#include "tools/combine_gvcfs_on_sites.hpp"


JNIEXPORT void JNICALL Java_org_bgi_flexlab_dpgt_jointcalling_CombineGVCFsOnSites_Combine(
    JNIEnv *env, jobject java_this, jobjectArray vcfpaths, jstring refpath,
    jstring outpath, jbyteArray bytes, jstring chrom, jlong start, jlong end)
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
    const char *refpath_cstr = env->GetStringUTFChars(refpath, NULL);
    const char *outpath_cstr = env->GetStringUTFChars(outpath, NULL);

    MultiVcfReader reader{vcfpaths_vec, true};
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
