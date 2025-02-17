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
