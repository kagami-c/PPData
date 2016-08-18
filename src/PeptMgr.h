#pragma once

#include "FastaData.h"
#include "PeptGenerator.h"
#include <unordered_set>
#include <algorithm>

// the instance of this class is required to be valid until the program ends.
class PeptMgr {
public:
    PeptMgr(const char* filename, bool append_decoy, PeptGenerator::EnzymeType enzyme_type, unsigned max_miss_cleavage) 
            : PeptMgr(filename, append_decoy, enzyme_type, max_miss_cleavage, 700, 5000) {}  // default mass range [700 Da, 5000 Da]

    PeptMgr(const char* filename, 
            bool append_decoy, 
            PeptGenerator::EnzymeType enzyme_type, 
            unsigned max_miss_cleavage, 
            double min_mass, 
            double max_mass) : fasta_data_(filename, append_decoy), 
                               digester_(enzyme_type, max_miss_cleavage) {
        
        // use hash pool to remove duplicate
        std::unordered_set<PeptGenerator::Peptide> peptide_pool;
        auto proteins = fasta_data_.GetProteins();
        for (auto& p : proteins) {
            auto pts = digester_.Digest(p);
            for (auto& pt : pts) {
                if (min_mass <= pt.molecule_weight && pt.molecule_weight <= max_mass) {
                    peptide_pool.insert(std::move(pt));
                }
            }
        }

        // move into vector and sort it
        peptides_.resize(peptide_pool.size());
        auto index = 0;
        for (auto&& peptide : peptide_pool) {
            peptides_[index++] = peptide;
        }
        std::sort(peptides_.begin(), peptides_.end(), 
            [](const auto& one, const auto& another) {
                return one.molecule_weight < another.molecule_weight;
            }
        );
    }

    // return a sorted peptides vector
    const std::vector<PeptGenerator::Peptide>& Retrieve() const {
        return peptides_;
    }

    // return type is too long
    auto RetrieveMassRange(double min, double max) const {
        auto start = std::lower_bound(peptides_.begin(), peptides_.end(), min,
            [](const auto& peptide, const auto& val) {
                return peptide.molecule_weight < val;
            }
        );
        auto end = std::upper_bound(peptides_.begin(), peptides_.end(), max,
            [](const auto& val, const auto& peptide) {
                return val < peptide.molecule_weight;
            }
        );
        return make_pair(start, end);
    }

private:
    const FastaData fasta_data_;
    const PeptGenerator digester_;
    std::vector<PeptGenerator::Peptide> peptides_;
};