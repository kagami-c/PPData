// Copyright (C) 2016

#include "PPData.h"
#include "PeptMgr.h"

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

// iterator implementation
class PPData::const_iterator::Impl {
public:
    Impl(const Peptide* ptr) : ptr_(ptr) {}
    void increment() { ++ptr_; }
    const Peptide* ptr() const { return ptr_; }

private:
    const Peptide* ptr_;
};

// const_iterator adapters
PPData::const_iterator::const_iterator(const Peptide& peptide) : pImpl(std::make_unique<Impl>(&peptide)) {}
PPData::const_iterator::const_iterator(const const_iterator& it) : pImpl(std::make_unique<Impl>(it.pImpl->ptr())) {}
PPData::const_iterator::~const_iterator() {}

PPData::const_iterator& PPData::const_iterator::operator++() { pImpl->increment(); return *this; }
PPData::const_iterator PPData::const_iterator::operator++(int) { const_iterator old(*pImpl->ptr()); ++(*this); return old; }

bool PPData::const_iterator::operator==(const_iterator other) const { return pImpl->ptr() == other.pImpl->ptr(); }
bool PPData::const_iterator::operator!=(const_iterator other) const { return pImpl->ptr() != other.pImpl->ptr(); }
const PPData::Peptide& PPData::const_iterator::operator*() const { return *pImpl->ptr(); }

// range object implementation
PPData::range::range(const_iterator begin, const_iterator end) : begin_(begin), end_(end) {}
PPData::const_iterator PPData::range::begin() const { return begin_; }
PPData::const_iterator PPData::range::end() const { return end_; }

// wrapper implementation
class PPData::Impl {
public:
    Impl(const char* filename, bool append_decoy, EnzymeType enzyme_type,
         unsigned max_miss_cleavage, double min_mass, double max_mass)
            : prot_data_(filename, append_decoy),
              pept_data_(prot_data_, enzyme_type, max_miss_cleavage, min_mass, max_mass) {}

    size_t size() const { return pept_data_.size(); }
    const_iterator begin() const { return const_iterator(*pept_data_.begin()); }
    const_iterator end() const {
        auto end_ptr = &*pept_data_.begin() + pept_data_.size();
        auto it = begin();
        it.pImpl->ptr = end_ptr;
        return const_iterator(it);
    }
    
    const_iterator lower_bound(double lower_mass) const {
        
    };
    const_iterator upper_bound(double upper_mass) const {
        
    };

private:
    ProtData prot_data_;
    PeptData pept_data_;
};

// ctors
PPData::PPData(const char* filename, bool append_decoy, EnzymeType enzyme_type,
               unsigned max_miss_cleavage, double min_mass, double max_mass)
               : pImpl(std::make_unique<Impl>(filename, append_decoy, enzyme_type,
                                              max_miss_cleavage, min_mass, max_mass)) {}
PPData::~PPData() {}

// adapters
size_t PPData::size() const { return pImpl->size(); }
PPData::const_iterator PPData::begin() const { return pImpl->begin(); }
PPData::const_iterator PPData::end() const { return pImpl->end(); }
PPData::const_iterator PPData::lower_bound(double mass) const { return pImpl->lower_bound(mass); }
PPData::const_iterator PPData::upper_bound(double mass) const { return pImpl->upper_bound(mass); }
PPData::range PPData::LoopWithin(double lower_mass, double upper_mass) const {
    return range(lower_bound(lower_mass), upper_bound(upper_mass));
}
