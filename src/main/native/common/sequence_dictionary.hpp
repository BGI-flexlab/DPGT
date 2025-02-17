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
#ifndef SEQUENCE_DICTIONARY_HPP
#define SEQUENCE_DICTIONARY_HPP

#include <cstdint>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <sstream>
#include "htslib/hts.h"
#include "htslib/sam.h"


struct SequenceDictionaryRecord
{
    std::string tag;   // type of the record. eg. "SQ"
    std::unordered_map<std::string, std::string> fields;
    std::vector<std::string> fields_keys;  // order of keys

    int64_t seq_length = -1;  // SEQ specific field

    std::string Str() const;
};


class SequenceDictionary
{
private:
    static const std::string kVersion;
    std::vector<SequenceDictionaryRecord> records;

    std::vector<SequenceDictionaryRecord *> sequence_records;  // @SQ

    // key: sequence name, value: sequence record idx in sequence_records
    std::unordered_map<std::string, int> sn_indices;

    SequenceDictionaryRecord MakeHeaderRecord();

    void InitSequenceRecords();

public:
    SequenceDictionary() = default;

    /**
     * Construct object from sam header.
     * 'M5' and 'UR' fields of the sequence record are absent.
     */
    SequenceDictionary(const sam_hdr_t *sam_header);

    SequenceDictionary(const SequenceDictionary &other);
    SequenceDictionary(SequenceDictionary &&other) noexcept;
    SequenceDictionary &operator=(const SequenceDictionary &other);
    SequenceDictionary &operator=(SequenceDictionary &&other) noexcept;

    void Clear();

    void InitializeFromSequenceFile(const std::string &sequence_file);
    void Load(const std::string &fasta_file);
    void Dump(const std::string &fasta_file) const;

    int32_t size() const {
        return static_cast<int32_t>(sequence_records.size());
    }

    int SequenceIndexToName(int sequence_index,
        std::string &sequence_name) const;
    int SequenceNameToIndex(const std::string &sequence_name,
        int &sequence_index) const;
    int SequenceIndexToLength(int sequence_index,
        int64_t &sequence_length) const;
};


#endif  // SEQUENCE_DICTIONARY_HPP
