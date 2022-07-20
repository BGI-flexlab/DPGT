#include "tools/NativeResponse.hpp"
#include "spdlog/spdlog.h"

JNIEXPORT void JNICALL Java_org_bgi_flexlab_dpgt_utils_NativeResponse_apply
    (JNIEnv *, jobject)
{
    spdlog::get("cdpgt")->info("Response from dpgt c/c++ library.");
}
