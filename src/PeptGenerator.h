#pragma once

#include "FastaData.h"
#include <vector>
#include <unordered_map>

class PeptGenerator {
public:
    struct Peptide {
        char n_term;
        char c_term;
        double molecule_weight;
        std::string sequence;  // sequence after convert all L to I
        std::string protein_name;  // merely keep the first protein name

        // for reporting, require the original data memory
        const char* original_sequence;  // before convert L to I, point to the memory in fasta, because it is just for reporting
        unsigned original_sequence_size;

        bool operator==(const Peptide& another) const {
            return this->sequence == another.sequence;
        }
    };

    enum class EnzymeType { Trypsin };

    PeptGenerator(EnzymeType enzyme_type, unsigned max_miss_cleavage)
            : enzyme_type_(enzyme_type), max_miss_cleavage_(max_miss_cleavage) {}

    // customized mass table
    PeptGenerator(EnzymeType enzyme_type, unsigned max_miss_cleavage, const std::unordered_map<char, double> mass_table)
            : enzyme_type_(enzyme_type), max_miss_cleavage_(max_miss_cleavage), mass_table_(mass_table) {}

    std::vector<Peptide> Digest(FastaData::Protein protein) const {
        auto sequence = std::string(protein.sequence, protein.sequence_size);
        ConvertIL(sequence);  // in-place convert L to I
        auto protein_name = ExtractProteinName(protein);

        // find all cleavage sites
        std::vector<unsigned> cleavage_sites;
        switch (enzyme_type_) {  // to support more enzymes, simply add different cleavage rules here
        case EnzymeType::Trypsin:
            cleavage_sites = TrypsinCleavageSites(sequence);
            break;
        default:
            throw std::runtime_error("Enzyme type is not supported.");
        }

        // calculate segments' mass and adjust max miss cleavages
        std::vector<double> segments_mass;
        segments_mass = SegmentsMass(sequence, cleavage_sites);  // if segment equals to 0, then we ignore it
        auto max_miss_cleavage = cleavage_sites.size() - 1 < max_miss_cleavage_
                                 ? cleavage_sites.size() - 1
                                 : max_miss_cleavage_;

        // pre-allocate vector space
        auto vector_size = 0;
        for (auto loop = 0; loop <= max_miss_cleavage; ++loop) {
            vector_size += (cleavage_sites.size() - loop);
        }
        std::vector<Peptide> peptides(vector_size);

        // handle miss cleavage, put smaller loop inside to accelerate the computation
        auto peptide_index = 0;
        for (auto index = 0; index < cleavage_sites.size(); ++index) {
            for (unsigned miss_cleavage = 0; miss_cleavage <= max_miss_cleavage; ++miss_cleavage) {
                auto start = cleavage_sites[index];
                auto end = index + miss_cleavage + 1 < cleavage_sites.size()
                           ? cleavage_sites[index + miss_cleavage + 1]
                           : sequence.size();
                auto len = end - start;
                auto converted_sequence = sequence.substr(start, len);
                auto mass = water_;
                for (unsigned i = 0; i < miss_cleavage + 1; ++i) {
                    auto local_mass = segments_mass[index + i];
                    mass += local_mass;
                    if (local_mass == 0) {
                        mass = 0;  // if some mass equals to zero, ignore it
                        break;
                    }
                }
                if (mass == 0) {  // ignore it
                    if (end == sequence.size()) {
                        break;
                    }
                    continue;
                }
                // a filter of allowed mass range can be added here.

                Peptide p = {
                    start == 0 ? '-' : sequence[start - 1],  // n_term
                    end == sequence.size() ? '-' : sequence[end],  // c_term
                    mass,  // molecule weight
                    converted_sequence,  // sequence
                    protein_name,  // protein_name
                    protein.sequence + start,  // original sequence pointer
                    static_cast<unsigned int>(len)  // original sequence len
                };
                peptides[peptide_index++] = p;
                if (end == sequence.size()) {
                    break;  // break the small loop
                }
            }
        }
        peptides.resize(peptide_index);  // make it compact
        return peptides;
    }

private:
    const EnzymeType enzyme_type_;
    const unsigned max_miss_cleavage_;

    // internal mass table
    const double proton_ = 1.00727;
    const double water_ = 18.01528;
    const std::unordered_map<char, double> mass_table_ = {
        { 'G', 57.02147 },{ 'A', 71.03712 },{ 'S', 87.03203 },{ 'P', 97.05277 },
        { 'V', 99.06842 },{ 'T', 101.04768 },{ 'C', 103.00919 + 57.021464 /* Fixed Mod on C */ },
        { 'I', 113.08407 }, /*{ 'L', 113.08407 },*/
        { 'N', 114.04293 },{ 'D', 115.02695 },{ 'Q', 128.05858 },
        { 'K', 128.09497 },{ 'E', 129.04260 },{ 'M', 131.04049 },{ 'H', 137.05891 },
        { 'F', 147.06842 },{ 'R', 156.10112 },{ 'Y', 163.06333 },{ 'W', 186.07932 }
    };

    // helper functions
    static void ConvertIL(std::string& sequence) {  // convert all L to I, in-place modification
        for (auto& c : sequence) {
            if (c == 'L') { c = 'I'; }
            else if (c == 'J') { c = 'I'; }  // J means I or L
        }
    }

    // Trypsin KR
    static std::vector<unsigned> TrypsinCleavageSites(const std::string& sequence) {
        std::vector<unsigned> cleavage_sites;
        cleavage_sites.push_back(0);
        for (auto index = 1; index < sequence.size(); ++index) {
            if ((sequence[index - 1] == 'K' || sequence[index - 1] == 'R') && sequence[index] != 'P') {
                cleavage_sites.push_back(index);
            }
        }
        return cleavage_sites;
    }

    // return the mass value in each segment, so that we don't have to re-compute them
    std::vector<double> SegmentsMass(const std::string& sequence, std::vector<unsigned> cleavage_sites) const {
        std::vector<double> segments_mass;
        for (auto i = 0; i < cleavage_sites.size(); ++i) {
            auto start = cleavage_sites[i];
            auto end = i == cleavage_sites.size() - 1
                       ? sequence.size()
                       : cleavage_sites[i + 1];
            double segment = 0;
            try {
                for (auto j = start; j < end; ++j) {
                    segment += mass_table_.at(sequence[j]);  // will throw an exception if no such character in table
                }
            }
            catch (std::out_of_range e) {
                segment = 0;  // if there is an error, the segment equals to 0
            }
            segments_mass.push_back(segment);
        }
        return segments_mass;
    }

    static std::string ExtractProteinName(const FastaData::Protein& protein) {
        unsigned i = 0;
        for (; i < protein.name_size; ++i) {
            if (protein.name[i] == ' ') {
                break;
            }
        }
        return std::string(protein.name, i);
    }
};

// hashing support
namespace std {
    template<>
    struct hash<PeptGenerator::Peptide> {
        size_t operator()(const PeptGenerator::Peptide& p) const {
            return (hash<string>()(p.sequence));
        }
    };
}
