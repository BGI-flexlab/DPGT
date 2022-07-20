#include "check_jemalloc.hpp"
#include "spdlog/spdlog.h"

bool check_jemalloc() {
    if (mallctl != nullptr) {
        const char *j;
        size_t s = sizeof(j);
        mallctl("version", &j,  &s, NULL, 0);
        spdlog::info("jemalloc was integrated, jemalloc version: {}\n", j);
        return true;
    } else {
        spdlog::warn("jemalloc not in use!");
        return false;
    }
}
