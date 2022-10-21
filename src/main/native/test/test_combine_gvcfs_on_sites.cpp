#include <boost/algorithm/string/join.hpp>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <future>
#include <string>
#include <vector>
#include "variant_site_finder/variant_site_finder.hpp"
#include "tools/combine_gvcfs_on_sites.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "ThreadPool.h"


static std::vector<std::vector<std::string>> chunkStringList
    (std::vector<std::string> string_list, int n)
{
    int p = string_list.size() / n;
    std::vector<std::vector<std::string>> result;
    if (p == 0) {
        spdlog::get("cdpgt")->error(
            "input string list length {} less than chunk number {}.", string_list.size(), n);
        std::exit(1);
    } else {
        for (int i = 0; i < n; ++i) {
            std::vector<std::string> r1;
            for (int j = 0; j < p; ++j) {
                r1.push_back(string_list[i*p+j]);
            }
            result.push_back(r1);
        }
        for (size_t i = p * n; i < string_list.size(); ++i) {
            result.back().push_back(string_list[i]);
        }
    }
    return result;
}

int main(int argc, char **argv) {
    auto cdpgt_logger = spdlog::stderr_color_mt("cdpgt");
    spdlog::set_level(spdlog::level::info);

    std::ifstream in(argv[1]);
    std::vector<std::string> vcfpaths;
    std::string line;
    while (std::getline(in, line)) {
        vcfpaths.push_back(line);
    }

    const char *chrom = "chr1";
    int64_t start = 2250900;
    int64_t end = 2260000;

    std::vector<std::vector<std::string>> vcfpaths_chunks = chunkStringList(vcfpaths, 23);
    ThreadPool *pool = new ThreadPool(10);
    std::vector<std::future<VariantSiteSet>> result;
    spdlog::get("cdpgt")->info("Finding variant sites...");
    for (auto &chunk: vcfpaths_chunks) {
        result.emplace_back(
                pool->enqueue(FindVariantSite, std::ref(chunk), chrom, start, end)
            );
    }

    VariantSiteSet vs_set{start, end};
    for (auto &&r: result) {
        vs_set.data() |= r.get().data();
    }

    delete pool;

    // spdlog::get("cdpgt")->info("Print variant sites...");
    for (int64_t i = 0; i < vs_set.size(); ++i) {
        if (vs_set.testIndex(i)) {
            std::cout << "VS\t" << vs_set.get(i) + 1 << std::endl;
        }
    }

    std::string prefix = argv[3];
    auto combine_func = [argv, chrom, start, end, &vs_set, &prefix](
        const std::vector<std::string> &vcfpaths1, int idx) -> void
    {
        if (idx == 11) {
            for (auto &s: vcfpaths1) {
                std::cout << s << std::endl;
            }
        }
        spdlog::get("cdpgt")->debug("chunk index: {}, vcfs: {}", idx,
            boost::join(vcfpaths1, " "));
        MultiVcfReader reader{vcfpaths1, true};

        reader.Querys(chrom, start, end);

        CachedRefSeq cached_ref_seq{argv[2]};

        VcfWriter vcf_writer{prefix+"."+std::to_string(idx)+".vcf.gz"};

        CombineGVCFsOnSites combiner{&reader, &cached_ref_seq, &vs_set,
            chrom, &vcf_writer};
        
        combiner.run();
    };

    ThreadPool *pool1 = new ThreadPool(1);

    std::vector<std::future<void>> result1;
    spdlog::get("cdpgt")->info("Combining gvcfs on variant sites...");
    for (int i = 0; i < (int)vcfpaths_chunks.size(); ++i) {
        result1.emplace_back(
                pool1->enqueue(combine_func, std::ref(vcfpaths_chunks[i]), i)
            );
    }

    for (auto &&r: result1) {
        r.get();
    }
    
    delete pool1;

    spdlog::get("cdpgt")->info("Done.");

    return 0;
}

