cmake_minimum_required(VERSION 3.10)
# project and version
project(dpgt LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 11)

set(JAVA_INCLUDE_PATH "/path/to/java8/include" CACHE STRING "java8 include path")

include_directories(vendor)
include_directories(vendor/htslib)
include_directories(vendor/eigen)
include_directories(vendor/cgranges/cpp)
include_directories(vendor/spdlog/include)
include_directories(vendor/ThreadPool)
include_directories(${JAVA_INCLUDE_PATH})
include_directories(${JAVA_INCLUDE_PATH}/linux)
include_directories(".")


# PROJECT_BINARY_DIR default is the build dir
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

if (NOT EXISTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    execute_process(COMMAND mkdir -p ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif (NOT EXISTS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

if (NOT EXISTS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND mkdir -p ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif (NOT EXISTS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(CMAKE_BUILD_TYPE Release CACHE STRING "build type is Release")

# CXX configure
if (CMAKE_BUILD_TYPE STREQUAL Debug)
    message(STATUS "[STATUS] build type is Debug")
    # -fprofile-arcs -ftest-coverage
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -pedantic -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall -Wextra -pedantic -fPIC")
elseif (CMAKE_BUILD_TYPE STREQUAL Release)
    message(STATUS "[STATUS] build type is Release")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC")
endif ()


set(MAKE_JOBS 4 CACHE STRING "number of jobs for run make in this cmake list")

set(USE_LIBDEFLATE ON CACHE STRING "if use libdeflate in htslib?")

# libdeflate library and include files in vendor
set(LIBDEFLATE ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate)
set(LIBDEFLATE_LIB ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate/libdeflate.so)

if (${USE_LIBDEFLATE})
    if (NOT EXISTS ${LIBDEFLATE})
        message(FATAL_ERROR "[FATAL_ERROR]: Can not find libdeflate directory in \
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor, Please check vendor directory \
    for existence of libdeflate.")
    endif (NOT EXISTS ${LIBDEFLATE})
endif (${USE_LIBDEFLATE})

if (${USE_LIBDEFLATE})
    if (EXISTS ${LIBDEFLATE_LIB})
        message(STATUS "[STATUS]: Find libdeflate library libdeflate.so")
    else (EXISTS ${LIBDEFLATE_LIB})
        message(STATUS "[STATUS]: libdeflate library libdeflate.so not exist in \
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate, build one ...")
        execute_process(COMMAND make -j ${MAKE_JOBS} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/libdeflate)
    endif (EXISTS ${LIBDEFLATE_LIB})
    execute_process(COMMAND cp ${LIBDEFLATE}/libdeflate.a ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND cp ${LIBDEFLATE}/libdeflate.so ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND cp ${LIBDEFLATE}/libdeflate.so.0 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif (${USE_LIBDEFLATE})


# htslib library and include files in vendor
set(HTSLIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
set(HTSLIB_LIB ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/libhts.so)

if (NOT EXISTS ${HTSLIB_DIR})
        message(FATAL_ERROR "[FATAL_ERROR]: Can not find htslib directory in \
${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib, Please check vendor directory \
for existence of htslib.")
endif (NOT EXISTS ${HTSLIB_DIR})


if (EXISTS ${HTSLIB_LIB})
    message(STATUS "[STATUS]: Find htslib shared library libhts.so")
else (EXISTS ${HTSLIB_LIB})
    message(STATUS "[STATUS]: htslib shared library libhts.so not exist in ${HTSLIB_DIR}, build one ...")
    execute_process(COMMAND autoheader WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    execute_process(COMMAND autoconf WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    if (${USE_LIBDEFLATE})
        execute_process(COMMAND bash -c "./configure CPPFLAGS=-I${LIBDEFLATE} LDFLAGS='-L${CMAKE_LIBRARY_OUTPUT_DIRECTORY} -Wl,-R${CMAKE_LIBRARY_OUTPUT_DIRECTORY}' --with-libdeflate" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    else (${USE_LIBDEFLATE})
        execute_process(COMMAND bash -c "./configure" WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    endif (${USE_LIBDEFLATE})
    execute_process(COMMAND make -j ${MAKE_JOBS} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
    execute_process(COMMAND cp ${HTSLIB_LIB} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    execute_process(COMMAND cp ${HTSLIB_LIB}.3 ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif (EXISTS ${HTSLIB_LIB})


find_package(Boost 1.74 REQUIRED COMPONENTS filesystem)
include_directories(${Boost_INCLUDE_DIRS})


add_subdirectory(vendor/spdlog)


add_executable(test_vcf_reader
    vcf/vcf_reader.cpp
    vcf/vcf_ibuffer.cpp
    test/test_vcf_reader.cpp)

target_link_libraries(test_vcf_reader
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    z m bz2 lzma pthread curl ssl crypto)


add_executable(grouper
    vcf/hts_vcf_utils.cpp
    vcf/vcf_reader.cpp
    vcf/vcf_ibuffer.cpp
    vcf/multi_vcf_reader.cpp
    vcf/vcf_constants.cpp
    vcf/variant_context.cpp
    vcf/allele.cpp
    test/test_multi_vcf_group_on_start.cpp)

target_link_libraries(grouper
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    z m bz2 lzma pthread curl ssl crypto)

add_executable(test_genotype_count_tables
    common/utils.cpp
    common/math_utils.cpp
    vcf/vcf_constants.cpp
    genotyper/genotype_allele_counts.cpp
    genotyper/genotype_likelihood_calculators.cpp
    test/test_genotype_count_tables.cpp)

target_link_libraries(test_genotype_count_tables
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    ${Boost_LIBRARIES})


add_executable(get_vcf_formats
    vcf/vcf_reader.cpp
    test/test_get_vcf_formats.cpp)

target_link_libraries(get_vcf_formats
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    z m bz2 lzma pthread curl ssl crypto)

add_executable(get_fmt_kvs
    vcf/vcf_reader.cpp
    test/test_get_fmt_kvs.cpp)

target_link_libraries(get_fmt_kvs
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    z m bz2 lzma pthread curl ssl crypto)

add_executable(test_vcf_attribute
    vcf/vcf_constants.cpp
    vcf/vcf_reader.cpp
    vcf/vcf_attribute.cpp
    test/test_vcf_attribute.cpp)

target_link_libraries(test_vcf_attribute
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    z m bz2 lzma pthread curl ssl crypto)


add_executable(combine-gvcfs
    common/cached_ref_seq.cpp
    common/histogram.cpp
    common/simple_interval.cpp
    common/fastx_seq.cpp
    common/fastx_suffix.cpp
    common/int_to_double_cache.cpp
    common/math_utils.cpp
    common/reference_context.cpp
    common/sequence_dictionary.cpp
    common/simple_interval.cpp
    common/utils.cpp
    genotyper/genotype_allele_counts.cpp
    genotyper/genotype_likelihood_calculator.cpp
    genotyper/genotype_likelihood_calculators.cpp
    genotyper/genotype_likelihoods.cpp
    genotyper/genotype.cpp
    genotyper/flat_genotype.cpp
    tools/combine_gvcfs.cpp
    tools/multi_variant_group_on_start.cpp
    tools/reference_confident_variant_merger.cpp
    vcf/vcf_id_table.cpp
    vcf/allele.cpp
    vcf/gatk_vcf_constants.cpp
    vcf/vcf_constants.cpp
    vcf/hts_vcf_utils.cpp
    vcf/multi_vcf_reader.cpp
    vcf/variant_builder.cpp
    vcf/variant_context_utils.cpp
    vcf/variant_context.cpp
    vcf/vcf_attribute.cpp
    vcf/vcf_shared_attribute.cpp
    vcf/vcf_ibuffer.cpp
    vcf/vcf_reader.cpp
    vcf/vcf_writer.cpp
    )


target_link_libraries(combine-gvcfs
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    ${Boost_LIBRARIES}
    z m bz2 lzma pthread curl ssl crypto profiler)


add_executable(test-genotype-allele-counts
    common/math_utils.cpp
    vcf/vcf_constants.cpp
    genotyper/genotype_allele_counts.cpp
    test/test_genotype_allele_counts.cpp)

target_link_libraries(test-genotype-allele-counts
    ${Boost_LIBRARIES})


add_executable(test-genotype-lk-calculator
    common/utils.cpp
    common/math_utils.cpp
    vcf/vcf_constants.cpp
    genotyper/genotype_allele_counts.cpp
    genotyper/genotype_likelihood_calculator.cpp
    genotyper/genotype_likelihood_calculators.cpp
    test/test_genotype_likelihood_calculator.cpp)

target_link_libraries(test-genotype-lk-calculator
    ${Boost_LIBRARIES})

add_library(cdpgt SHARED
    common/cached_ref_seq.cpp
    common/histogram.cpp
    common/simple_interval.cpp
    common/fastx_seq.cpp
    common/fastx_suffix.cpp
    common/int_to_double_cache.cpp
    common/math_utils.cpp
    common/reference_context.cpp
    common/sequence_dictionary.cpp
    common/simple_interval.cpp
    common/utils.cpp
    genotyper/genotype_allele_counts.cpp
    genotyper/genotype_likelihood_calculator.cpp
    genotyper/genotype_likelihood_calculators.cpp
    genotyper/genotype_likelihoods.cpp
    genotyper/genotype.cpp
    genotyper/flat_genotype.cpp
    tools/combine_gvcfs_on_sites.cpp
    tools/CombineGVCFsOnSites.cpp
    tools/reference_confident_variant_merger.cpp
    variant_site_finder/variant_site_finder.cpp
    vcf/vcf_id_table.cpp
    vcf/allele.cpp
    vcf/gatk_vcf_constants.cpp
    vcf/vcf_constants.cpp
    vcf/hts_vcf_utils.cpp
    vcf/multi_vcf_reader.cpp
    vcf/variant_builder.cpp
    vcf/variant_context_utils.cpp
    vcf/variant_context.cpp
    vcf/vcf_attribute.cpp
    vcf/vcf_shared_attribute.cpp
    vcf/vcf_ibuffer.cpp
    vcf/vcf_reader.cpp
    vcf/vcf_writer.cpp
    variant_site_finder/variant_site_finder.cpp
    variant_site_finder/VariantSiteFinder.cpp)

target_link_libraries(cdpgt
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    ${Boost_LIBRARIES}
    z m bz2 lzma pthread curl ssl crypto)


add_executable(test_combine_gvcfs_on_sites
    test/test_combine_gvcfs_on_sites.cpp
    common/cached_ref_seq.cpp
    common/histogram.cpp
    common/simple_interval.cpp
    common/fastx_seq.cpp
    common/fastx_suffix.cpp
    common/int_to_double_cache.cpp
    common/math_utils.cpp
    common/reference_context.cpp
    common/sequence_dictionary.cpp
    common/simple_interval.cpp
    common/utils.cpp
    genotyper/genotype_allele_counts.cpp
    genotyper/genotype_likelihood_calculator.cpp
    genotyper/genotype_likelihood_calculators.cpp
    genotyper/genotype_likelihoods.cpp
    genotyper/genotype.cpp
    genotyper/flat_genotype.cpp
    tools/combine_gvcfs_on_sites.cpp
    tools/reference_confident_variant_merger.cpp
    variant_site_finder/variant_site_finder.cpp
    vcf/vcf_id_table.cpp
    vcf/allele.cpp
    vcf/gatk_vcf_constants.cpp
    vcf/vcf_constants.cpp
    vcf/hts_vcf_utils.cpp
    vcf/multi_vcf_reader.cpp
    vcf/variant_builder.cpp
    vcf/variant_context_utils.cpp
    vcf/variant_context.cpp
    vcf/vcf_attribute.cpp
    vcf/vcf_shared_attribute.cpp
    vcf/vcf_ibuffer.cpp
    vcf/vcf_reader.cpp
    vcf/vcf_writer.cpp
    )

target_link_libraries(test_combine_gvcfs_on_sites
    ${HTSLIB_LIB}
    ${LIBDEFLATE_LIB}
    ${Boost_LIBRARIES}
    z m bz2 lzma pthread curl ssl crypto spdlog::spdlog)
