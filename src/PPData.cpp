// Copyright (C) 2016

#include "PPData.h"
#include "ProtData.h"
#include "PeptData.h"

// data interface ctor
PPData::Protein::Protein(const char* name, const char* sequence, size_t sequence_length)
        : name(name), sequence(sequence), sequence_length(sequence_length) {}

PPData::Peptide::Peptide(const Protein& protein, const char* compact_protein_sequence, 
                         size_t start_idx, size_t end_idx, double mass)
        : sequence(compact_protein_sequence + start_idx),
          sequence_length(end_idx - start_idx),
          n_term(start_idx == 0 ? '-' : protein.sequence[start_idx - 1]),
          c_term(end_idx == protein.sequence_length ? '-' : protein.sequence[end_idx]),
          mass(mass),
          protein(&protein),
          offset(start_idx) {}

// wrapper implementation
class PPData::Impl {
public:
    Impl(const char* filename, bool append_decoy, EnzymeType enzyme_type,
         unsigned max_miss_cleavage, double min_mass, double max_mass)
            : prot_data_(filename, append_decoy),
              pept_data_(prot_data_, enzyme_type, max_miss_cleavage, min_mass, max_mass) {}

    size_t size() const { return pept_data_.size(); }
    auto begin() const { return pept_data_.begin(); }
    auto end() const { return pept_data_.end(); }
    auto lower_bound(double lower_mass) const { return pept_data_.lower_bound(lower_mass); }
    auto upper_bound(double upper_mass) const { return pept_data_.upper_bound(upper_mass); }

private:
    ProtData prot_data_;
    PeptData pept_data_;
};

// container ctor
PPData::PPData(const char* filename, bool append_decoy, EnzymeType enzyme_type,
               unsigned max_miss_cleavage, double min_mass, double max_mass)
               : pImpl(std::make_unique<Impl>(filename, append_decoy, enzyme_type,
                                              max_miss_cleavage, min_mass, max_mass)) {}
PPData::~PPData() {}

// adapters
size_t PPData::size() const { return pImpl->size(); }
auto PPData::begin() const { return pImpl->begin(); }
auto PPData::end() const { return pImpl->end(); }
auto PPData::lower_bound(double mass) const { return pImpl->lower_bound(mass); }
auto PPData::upper_bound(double mass) const { return pImpl->upper_bound(mass); }

auto PPData::LoopWithin(double lower_mass, double upper_mass) const {
    return range<std::vector<Peptide>::const_iterator>(lower_bound(lower_mass), upper_bound(upper_mass));
}
