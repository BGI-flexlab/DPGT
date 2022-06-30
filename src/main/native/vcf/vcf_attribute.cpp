#include "vcf/vcf_attribute.hpp"



const std::map<int, std::string> VcfAttributeBase::TYPE_STRINGS = {
    {BCF_HT_FLAG, "FLAG"},
    {BCF_HT_INT, "INT"},
    {BCF_HT_REAL, "FLOAT"},
    {BCF_HT_STR, "STRING"},
    {BCF_HT_LONG, "LONG"}
};


const std::array<std::string, 5> VcfAttributeBase::SIZE_TYPE_STRINGS = {
    "FIXED", "VAR", "A", "G", "R"};
