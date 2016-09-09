#pragma once

#include "PPData.h"
#include <fstream>
#include <vector>
#include <cassert>
#include <stdexcept>

class ProtData {
public:
    using Protein = PPData::Protein;

    ProtData(const char* filename, bool append_decoy) {
        // read data into memory
        std::basic_ifstream<char> file(filename, std::ios::binary);
        if (!file) { throw std::runtime_error("Fail to open fasta database file."); }
        file.unsetf(std::ios::skipws);
        file.seekg(0, std::ios::end);
        size_t size = file.tellg();
        file.seekg(0);
        std::vector<char> raw_data(size + 1);  // one more character for '\0'
        file.read(&raw_data.front(), static_cast<std::streamsize>(size));
        raw_data[size] = 0;

        // concatenate sequence in memory, by copying the data to a new place
        target_data_.resize(raw_data.size());
        auto state = ParseState::Name;
        auto index = 0;
        for (auto& c : raw_data) {
            switch (state) {
            case ParseState::Start:
                if (c == '>') { state = ParseState::Name; }
                break;
            case ParseState::Name:
                if (c == '\n') {
                    target_data_[index++] = '\0';
                    state = ParseState::Sequence;
                }
                else {
                    target_data_[index++] = c;
                }
                break;
            case ParseState::Sequence:
                if (c == '>') {
                    target_data_[index++] = '\0';
                    target_data_[index++] = '>';
                    state = ParseState::Name;
                }
                else if (c != ' ' && c != '\r' && c != '\n' && c != '\t') {
                    target_data_[index++] = c;
                }
                break;
            }
        }
        target_data_.resize(index);

        // build proteins
        state = ParseState::Start;  // reuse the same state
        const char* temp_name = nullptr;
        for (auto i = 0; i < target_data_.size(); ++i) {
            switch (state) {
            case ParseState::Start:
                if (target_data_[i] == '>') { state = ParseState::Name; }
                break;
            case ParseState::Name:
                temp_name = &target_data_[i];
                while (target_data_[i] != '\0') { ++i; }
                state = ParseState::Sequence;
                break;
            case ParseState::Sequence:
                const char* temp_sequence = &target_data_[i];
                while (target_data_[i] != '\0') { ++i; }
                size_t temp_length = &target_data_[i] - temp_sequence;
                // build protein
                proteins_.push_back(Protein(temp_name, temp_sequence, temp_length));
                state = ParseState::Start;
                break;
            }
        }

        // build decoy, if required
        if (append_decoy) {
            auto target_protein_num = proteins_.size();
            auto decoy_datamap_size = target_data_.size() + target_protein_num * 6;  // add DECOY_ prefix before protein name
            decoy_data_.resize(decoy_datamap_size);

            auto decoy_index = 0;
            for (auto i = 0; i < target_protein_num; ++i) {
                const char* decoy_name;
                const char* decoy_sequence;
                size_t decoy_seq_len;

                // build each decoy protein
                decoy_data_[decoy_index++] = '>';
                decoy_name = &decoy_data_[decoy_index];
                FillDecoyPrefix(decoy_index);

                auto& target_protein = proteins_[i];
                for (unsigned j = 0; target_protein.name[j] != '\0'; ++j) {
                    decoy_data_[decoy_index++] = target_protein.name[j];
                }
                decoy_data_[decoy_index++] = '\0';

                decoy_sequence = &decoy_data_[decoy_index];
                for (int j = target_protein.sequence_length - 1; target_protein.sequence[j] != '\0'; --j) {
                    decoy_data_[decoy_index++] = target_protein.sequence[j];
                }
                decoy_seq_len = &decoy_data_[decoy_index] - decoy_sequence;
                decoy_data_[decoy_index++] = '\0';

                // build decoy protein
                proteins_.push_back(Protein(decoy_name, decoy_sequence, decoy_seq_len));
            }
            assert(decoy_data_.size() == decoy_index);
        }
    }

    const std::vector<Protein>& GetProteins() const {
        return proteins_;
    }

private:
    std::vector<char> target_data_;
    std::vector<char> decoy_data_;
    std::vector<Protein> proteins_;

    enum class ParseState { Start, Name, Sequence };  // reuse twice

    void FillDecoyPrefix(int& index) {
        const char* const prefix = "DECOY_";
        for (auto i = 0; i < 6; ++i) {
            decoy_data_[index++] = prefix[i];
        }
    }
};