#include "tools/NativeInitializer.hpp"
#include "common/check_jemalloc.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

JNIEXPORT void JNICALL Java_org_bgi_flexlab_dpgt_jointcalling_NativeInitializer_apply
    (JNIEnv *, jobject)
{
    spdlog::info("Initialize JNI environment.");
    // create a spdlog logger named 'cdpgt'
    auto cdpgt_logger = spdlog::stderr_color_mt("cdpgt");
    cdpgt_logger->set_level(spdlog::level::info);
    // check if jemalloc is in use
    spdlog::get("cdpgt")->info("Check if jemalloc is in use.");
    check_jemalloc();
}
