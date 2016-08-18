#pragma once

#include <fstream>
#include <vector>
#include <cassert>
#include <stdexcept>

class FastaData {
public:
    struct Protein {
        unsigned name_size;
        unsigned sequence_size;
        const char* name;
        const char* sequence;
    };  // arrange for memory alignment

    FastaData(const char* filename, bool append_decoy) : filename_(filename) {
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
        data_.resize(raw_data.size());
        auto state = ParseState::Name;
        auto index = 0;
        for (auto& c : raw_data) {
            switch (state) {
            case ParseState::Start:
                if (c == '>') { state = ParseState::Name; }
                break;
            case ParseState::Name:
                data_[index++] = c;
                if (c == '\n') { state = ParseState::Sequence; }
                break;
            case ParseState::Sequence:
                if (c == '>') {
                    data_[index++] = '\n';
                    data_[index++] = '>';
                    state = ParseState::Name;
                }
                else if (c != ' ' && c != '\r' && c != '\n' && c != '\t') {
                    data_[index++] = c;
                }
                break;
            }
        }
        data_.resize(index);

        // build proteins
        state = ParseState::Start;  // reuse the same state
        Protein temp;
        for (auto i = 0; i < data_.size(); ++i) {
            switch (state) {
            case ParseState::Start:
                if (data_[i] == '>') { state = ParseState::Name; }
                break;
            case ParseState::Name:
                temp.name = &data_[i];
                while (data_[i] != '\n') { ++i; }
                temp.name_size = &data_[i] - temp.name;
                state = ParseState::Sequence;
                break;
            case ParseState::Sequence:
                temp.sequence = &data_[i];
                while (data_[i] != '\n' && data_[i] != '\0') { ++i; }
                temp.sequence_size = &data_[i] - temp.sequence;
                // build protein
                proteins_.push_back(temp);
                state = ParseState::Start;
                break;
            }
        }

        // build decoy, if required
        if (append_decoy) {
            auto target_protein_num = proteins_.size();
            auto decoy_datamap_size = data_.size() + target_protein_num * 6;  // add DECOY_ prefix before each protein name
            decoy_data_.resize(decoy_datamap_size);

            proteins_.resize(target_protein_num * 2);
            auto decoy_index = 0;
            for (auto i = 0; i < target_protein_num; ++i) {
                Protein decoy_protein;

                // build each decoy protein
                decoy_data_[decoy_index++] = '>';
                decoy_protein.name = &decoy_data_[decoy_index];
                FillDecoyPrefix(decoy_index);

                auto& target_protein = proteins_[i];
                for (unsigned j = 0; j < target_protein.name_size; ++j) {
                    decoy_data_[decoy_index++] = target_protein.name[j];
                }
                decoy_protein.name_size = &decoy_data_[decoy_index] - decoy_protein.name;
                decoy_data_[decoy_index++] = '\n';

                decoy_protein.sequence = &decoy_data_[decoy_index];
                for (unsigned j = 0; j < target_protein.sequence_size; ++j) {
                    decoy_data_[decoy_index++] = target_protein.sequence[target_protein.sequence_size - 1 - j];
                }
                decoy_protein.sequence_size = &decoy_data_[decoy_index] - decoy_protein.sequence;
                decoy_data_[decoy_index++] = '\n';

                proteins_[target_protein_num + i] = decoy_protein;
            }
            assert(decoy_data_.size() == decoy_index);
            decoy_data_[decoy_index - 1] = '\0';  // the last position is already set '\n'
        }
    }

    std::vector<Protein> GetProteins() const {
        return proteins_;
    }

private:
    const char* filename_;
    std::vector<char> data_;
    std::vector<char> decoy_data_;
    std::vector<Protein> proteins_;

    enum class ParseState { Start, Name, Sequence };  // used in both ctor and load

    void FillDecoyPrefix(int& index) {
        const char* prefix = "DECOY_";
        for (auto i = 0; i < 6; ++i) {
            decoy_data_[index++] = prefix[i];
        }
    }
};
