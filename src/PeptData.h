#pragma once

#include "PPData.h"
#include "ProtData.h"
#include "Hash.h"  // hash support for PPData::Peptide
#include <vector>
#include <numeric>
#include <unordered_set>
#include <algorithm>
#include <unordered_map>

class PeptData {
public:
    using Protein = PPData::Protein;
    using Peptide = PPData::Peptide;
    using EnzymeType = PPData::EnzymeType;

    // TODO: provide an API for assigning a different mass table (possibly with modifications)
    PeptData(const ProtData& proteins, EnzymeType enzyme_type, unsigned max_miss_cleavage, 
             double min_mass, double max_mass)
            : enzyme_type_(enzyme_type), max_miss_cleavage_(max_miss_cleavage), 
              min_mass_(min_mass), max_mass_(max_mass) {
        // calculate compact size
        auto sequences_len_sum = std::accumulate(proteins.begin(), proteins.end(), 0,
            [](const int& sum, const Protein& protein) {
                return sum + protein.sequence_length;
            }
        );
        auto compact_size = sequences_len_sum + proteins.size();  // append '\0' at the end of each sequence
        compact_sequences_.resize(compact_size);

        // build compact sequences and peptides
        std::unordered_set<Peptide> pool;
        int index = 0;
        for (auto& protein : proteins) {
            // fill compact sequence
            const char* compact_sequence = &compact_sequences_[index];
            for (unsigned i = 0; i < protein.sequence_length; ++i) {
                compact_sequences_[index++] = protein.sequence[i] == 'I' ? 'L' : protein.sequence[i];
            }  // prefer L
            compact_sequences_[index++] = '\0';

            // digest and build peptides into pool
            Digest(pool, protein, compact_sequence);
        }

        peptides_.insert(peptides_.end(), pool.begin(), pool.end());
        std::sort(peptides_.begin(), peptides_.end(),
                  [](const auto& one, const auto& another) { return one.mass < another.mass; });
    }

    size_t size() const { return peptides_.size(); }
    auto begin() const { return peptides_.cbegin(); }
    auto end() const { return peptides_.cend(); }

    auto lower_bound(double lower_mass) const {
        auto start = std::lower_bound(peptides_.begin(), peptides_.end(), lower_mass,
            [](const auto& peptide, const auto& val) {
                return peptide.mass < val;
            }
        );
        return start;
    }
    auto upper_bound(double upper_mass) const {
        auto end = std::upper_bound(peptides_.begin(), peptides_.end(), upper_mass,
            [](const auto& val, const auto& peptide) {
                return val < peptide.mass;
            }
        );
        return end;
    }
    
private:
    const EnzymeType enzyme_type_;
    const unsigned max_miss_cleavage_;
    const double min_mass_;
    const double max_mass_;

    std::vector<char> compact_sequences_;  // after convert IL
    std::vector<Peptide> peptides_;

    // internal mass table
    const double proton_ = 1.00727;
    const double water_ = 18.01528;

    const std::unordered_map<char, double> mass_table_ = {
        { 'G', 57.02147 },{ 'A', 71.03712 },{ 'S', 87.03203 },{ 'P', 97.05277 },
        { 'V', 99.06842 },{ 'T', 101.04768 },{ 'C', 103.00919 + 57.021464 /* Fixed Mod on C */ },
        /*{ 'I', 113.08407 },*/ { 'L', 113.08407 },
        { 'N', 114.04293 },{ 'D', 115.02695 },{ 'Q', 128.05858 },
        { 'K', 128.09497 },{ 'E', 129.04260 },{ 'M', 131.04049 },{ 'H', 137.05891 },
        { 'F', 147.06842 },{ 'R', 156.10112 },{ 'Y', 163.06333 },{ 'W', 186.07932 }
    };

    // builders
    void Digest(std::unordered_set<Peptide>& pool, 
                const Protein& protein, const char* compact_sequence) const {
        auto cleavage_sites = GenCleavageSites(compact_sequence, protein.sequence_length);
        auto segments_mass = SegmentsMass(compact_sequence, protein.sequence_length, cleavage_sites);  // if segment equals to 0, then we ignore it
        auto local_max_miss_cleavage = cleavage_sites.size() - 1 < max_miss_cleavage_
                                       ? cleavage_sites.size() - 1 : max_miss_cleavage_;

        // handle miss cleavage, put smaller loop inside to accelerate the computation
        for (unsigned index = 0; index < cleavage_sites.size(); ++index) {
            for (unsigned miss_cleavage = 0; miss_cleavage <= local_max_miss_cleavage; ++miss_cleavage) {
                auto start = cleavage_sites[index];
                auto end = index + miss_cleavage + 1 < cleavage_sites.size()
                           ? cleavage_sites[index + miss_cleavage + 1] 
                           : protein.sequence_length;  // next char of the end

                auto mass = water_;
                for (unsigned i = 0; i < miss_cleavage + 1; ++i) {
                    auto local_mass = segments_mass[index + i];
                    mass += local_mass;
                    if (local_mass == 0) {
                        mass = 0;  // if some mass equals to zero, ignore it
                        break;
                    }
                }
                if (mass == 0 /* contain intractable amino acid */
                        || mass < min_mass_ || max_mass_ < mass /* mass ourside range */) {
                    if (end == protein.sequence_length) { break; }
                    else { continue; }
                }

                pool.insert(Peptide(protein, compact_sequence, start, end, mass));
                if (end == protein.sequence_length) { break; }  // break the small loop
            }
        }
    }

    std::vector<unsigned> GenCleavageSites(const char* compact_sequence, size_t sequence_length) const {
        std::vector<unsigned> cleavage_sites;
        cleavage_sites.push_back(0);
        switch (enzyme_type_) {  // to support more enzymes, simply add different cleavage rules here
        case EnzymeType::Trypsin:  // Trypsin KR
            for (unsigned index = 1; index < sequence_length; ++index) {
                if ((compact_sequence[index - 1] == 'K' || compact_sequence[index - 1] == 'R')
                    && compact_sequence[index] != 'P') {
                    cleavage_sites.push_back(index);
                }
            }
            break;
        default:
            throw std::runtime_error("Enzyme type is not supported.");
        }
        return cleavage_sites;
    }

    // return the mass value in each segment, so that we don't have to re-compute them
    std::vector<double> SegmentsMass(const char* sequence, size_t sequence_length, 
                                     const std::vector<unsigned>& cleavage_sites) const {
        std::vector<double> segments_mass;
        for (unsigned i = 0; i < cleavage_sites.size(); ++i) {
            auto start = cleavage_sites[i];
            auto end = i == cleavage_sites.size() - 1
                       ? sequence_length
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
};
