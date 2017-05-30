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
    const Peptide& operator[](const size_t index) const { return pept_data_[index]; }

private:
    ProtData prot_data_;
    PeptData pept_data_;
};

// container ctor
PPData::PPData(const char* filename, bool append_decoy, EnzymeType enzyme_type,
               unsigned max_miss_cleavage, double min_mass, double max_mass)
               : pImpl(std::make_unique<Impl>(filename, append_decoy, enzyme_type,
                                              max_miss_cleavage, min_mass, max_mass)) {}
PPData::PPData(const PPData& ppdata) : pImpl(new Impl(*ppdata.pImpl)) {}
PPData::~PPData() {}

// adapters
size_t PPData::size() const { return pImpl->size(); }
const PPData::Peptide& PPData::operator[](const size_t index) const { return pImpl->operator[](index); }
