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
