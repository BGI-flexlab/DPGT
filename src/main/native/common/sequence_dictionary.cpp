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
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <limits.h>
#include <cassert>


#include "common/sequence_dictionary.hpp"
#include "common/fastx_seq.hpp"
#include "common/fastx_suffix.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>


std::string SequenceDictionaryRecord::Str() const
{
    std::ostringstream outstr;
    outstr << "@" << tag << "\t";
    std::vector<std::string> fields_vec;
    for (auto const &it: fields_keys)
    {
        fields_vec.push_back(it + ":" + fields.at(it));
    }
    std::string fields_str = boost::join(fields_vec, "\t");
    outstr << fields_str << std::endl; 
    return outstr.str();
}

const std::string SequenceDictionary::kVersion = "1.0";

void SequenceDictionary::InitSequenceRecords() {
    int i = 0;
    for (auto &it: records)
    {
        if (it.tag == "SQ")
        {
            sequence_records.push_back(&it);
            sn_indices[it.fields["SN"]] = i;
            ++i;
        }
    }
}

SequenceDictionary::SequenceDictionary(const sam_hdr_t *sam_header)
{
    records.push_back(MakeHeaderRecord());
    int64_t len;
    for (int32_t tid = 0; tid < sam_header->n_targets; ++tid)
    {
        SequenceDictionaryRecord record;
        record.tag = "SQ";
        record.fields["SN"] = sam_hdr_tid2name(sam_header, tid);
        len = sam_hdr_tid2len(sam_header, tid);
        record.fields["LN"] = std::to_string(len);
        record.seq_length = len;

        record.fields_keys = {"SN", "LN"};
        records.push_back(record);
    }
    InitSequenceRecords();
}

SequenceDictionary::SequenceDictionary(const SequenceDictionary &other)
{
    records = other.records;
    InitSequenceRecords();
}

SequenceDictionary::SequenceDictionary(SequenceDictionary &&other) noexcept
{
    records = std::move(other.records);
    sequence_records = std::move(other.sequence_records);
    sn_indices = std::move(other.sn_indices);
}

SequenceDictionary &SequenceDictionary::operator=(
    const SequenceDictionary &other)
{
    if (this != &other)
    {
        records = other.records;
        InitSequenceRecords();
    }
    return *this;
}

SequenceDictionary &SequenceDictionary::operator=(
    SequenceDictionary &&other) noexcept
{
    if (this != &other)
    {
        records = std::move(other.records);
        sequence_records = std::move(other.sequence_records);
        sn_indices = std::move(other.sn_indices);
    }
    return *this;
}


void SequenceDictionary::Clear()
{
    sn_indices.clear();
    sequence_records.clear();
    records.clear();
}


void SequenceDictionary::InitializeFromSequenceFile(
    const std::string &sequence_file)
{
    Clear();

    // add header
    records.push_back(MakeHeaderRecord());

    // file realpath
    char sequence_file_path[PATH_MAX];
    if (realpath(sequence_file.c_str(), sequence_file_path) == NULL)
    {
        char errMsg[4096];
        sprintf(errMsg, "[SequenceDictionary::InitializeFromSequenceFile] "
            "Error! Cannot get realpath for %s.", sequence_file.c_str());
        perror(errMsg);
        std::exit(1);
    }

    // md5 related variables
    hts_md5_context *md5;
    unsigned char digest[16];
    char hex[33];  // md5 hex string

    if (!(md5 = hts_md5_init()))
    {
        std::cerr << "[SequenceDictionary::InitializeFromSequenceFile] "
            << "Error! Can not initialize md5 context!" << std::endl;
        std::exit(1);
    }

    FastxReader reader{sequence_file};
    FastxSeq seq;
    while (reader.Read(seq) >= 0)
    {
        SequenceDictionaryRecord record;
        record.tag = "SQ";
        record.fields["SN"] = seq.name;  // sequence name
        record.fields["LN"] = std::to_string(seq.seq.size());  // sequence length
        record.seq_length = seq.seq.size();

        hts_md5_reset(md5);
        hts_md5_update(md5, seq.seq.data(), seq.seq.size());
        hts_md5_final(digest, md5);
        hts_md5_hex(hex, digest);

        record.fields["M5"] = hex;
        record.fields["UR"] = "file://" + std::string(sequence_file_path);

        record.fields_keys = {"SN", "LN", "M5", "UR"};
        records.push_back(record);

    }

    hts_md5_destroy(md5);

    InitSequenceRecords();
}


void SequenceDictionary::Load(const std::string &fasta_file)
{
    std::string dict_file = FastxSuffix::TrimSuffix(fasta_file) + ".dict";
    if (!boost::filesystem::exists(dict_file))
    {
        dict_file = fasta_file + ".dict";
    }

    if (!boost::filesystem::exists(dict_file))
    {
        std::cerr << "Can not find fasta dict file for reference: "
                  << fasta_file << std::endl;
        std::exit(1);
    }

    Clear();

    std::ifstream input;
    input.open(dict_file);

    if (!input.good())
    {
        std::cerr << "[SequenceDictionary::Load] Error! Can not read dict file"
                     ": " << dict_file << std::endl;
        std::exit(1);
    }

    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty()) continue;
        if (line[0] != '@') continue;
        SequenceDictionaryRecord record;
        std::vector<std::string> elements;
        boost::split(elements, line, boost::is_any_of("\t"));
        record.tag = std::string(elements[0].begin()+1, elements[0].end());
        for (int i = 1; i < (int)elements.size(); ++i)
        {
            std::string::size_type sep_n = elements[i].find(':');
            std::string key = elements[i].substr(0, sep_n);
            std::string value = elements[i].substr(sep_n+1);
            record.fields[key] = value;
            record.fields_keys.push_back(key);
        }
        if (record.tag == "SQ") {
            record.seq_length = atol(record.fields.at("LN").c_str());
        }
        records.push_back(record);
    }

    InitSequenceRecords();

    input.close();
}

void SequenceDictionary::Dump(const std::string &fasta_file) const
{
    std::string fasta_base = FastxSuffix::TrimFastaSuffix(fasta_file);
    std::string dict_file = fasta_base + ".dict";
    std::ofstream output;
    output.open(dict_file, std::ios_base::out);
    for (auto &record: records)
    {
        output << record.Str();
    }
    output.close();
}


SequenceDictionaryRecord SequenceDictionary::MakeHeaderRecord()
{
    SequenceDictionaryRecord header_record;
    header_record.tag = "HD";
    header_record.fields["VN"] = kVersion;
    header_record.fields["SO"] = "unsorted";
    header_record.fields_keys = {"VN", "SO"};
    return header_record;
}


int SequenceDictionary::SequenceIndexToName(int sequence_index,
    std::string &sequence_name) const
{
    if (records.empty())
    {
        return -1;      // can not convert index to name when this is empty
    }

    if (sequence_records.empty())
    {
        return -2;  // can not convert index to name when sequence records is empty
    }

    if (sequence_index >= (int)sequence_records.size() ||
        sequence_index < 0)
    {
        return -3;  // sequence index must >= 0 and < (int)sequence_records.size()
    }

    sequence_name = sequence_records[sequence_index]->fields["SN"];
    return 0;  // success
}


int SequenceDictionary::SequenceNameToIndex(const std::string &sequence_name,
    int &sequence_index) const
{
    if (records.empty())
    {
        return -1;      // can not convert index to name when this is empty
    }

    if (sn_indices.empty())
    {
        return -2;  // can not convert index to name when sequence records is empty
    }

    bool is_find = false;

    if (sn_indices.find(sequence_name) != sn_indices.end())
    {
        is_find = true;
        sequence_index = sn_indices.at(sequence_name);
    }

    if (is_find) return 0;

    return -3;  // sequence name not in this dictinary
}

int SequenceDictionary::SequenceIndexToLength(int sequence_index,
    int64_t &sequence_length) const
{
    if (records.empty())
    {
        return -1;      // can not convert index to name when this is empty
    }

    if (sequence_records.empty())
    {
        return -2;  // can not convert index to name when sequence records is empty
    }

    if (sequence_index >= (int)sequence_records.size() ||
        sequence_index < 0)
    {
        return -3;  // sequence index must >= 0 and < (int)sequence_records.size()
    }

    sequence_length = sequence_records[sequence_index]->seq_length;
    return 0;  // success
}


// #ifdef 0

// int main(int argc, char **argv)
// {
//     if ( strcmp(argv[1], "init") == 0 )
//     {
//         glat::SequenceDictionary seq_dict;
//         seq_dict.InitializeFromSequenceFile(argv[2]);
//         seq_dict.Dump(argv[3]);
//     } else if ( strcmp(argv[1], "load") == 0 )
//     {
//         glat::SequenceDictionary seq_dict;
//         seq_dict.Load(argv[2]);
//         seq_dict.Dump(argv[3]);
//     } else
//     {
//         std::cerr << "Error! Unknown subcommand." << std::endl;
//         std::exit(1);
//     }
    
//     return 0;
// }

// #endif // 0
