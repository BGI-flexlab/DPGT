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
#include <vector>
#include <string>
#include "jni_md.h"
#include "tools/vcf_header_combiner.hpp"
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

    combineVCFHeaders(vcfpaths_vec, outpath_cstr);
    
    env->ReleaseStringUTFChars(outpath, outpath_cstr);
}
