// Copyright (C) 2016

#include "PeptData.h"
#include "PeptMgr.h"

// data interface ctor
PeptData::Protein::Protein(const char* name, const char* sequence, size_t sequence_length)
        : name(name), sequence(sequence), sequence_length(sequence_length) {}

PeptData::Peptide::Peptide(const Protein& protein, size_t start_idx, size_t end_idx, double mass)
        : sequence(protein.sequence + start_idx),
          sequence_length(end_idx - start_idx + 1),
          n_term(start_idx == 0 ? '-' : protein.sequence[start_idx - 1]),
          c_term(end_idx == protein.sequence_length - 1 ? '-' : protein.sequence[end_idx + 1]),
          mass(mass),
          protein(&protein),
          offset(start_idx) {}

// iterator implementation
class PeptData::const_iterator::Impl {
public:
    Impl(const Peptide& peptide) : ptr(&peptide) {}
    void increment() { ++ptr; }
    const Peptide& element() const { return *ptr; }

private:
    const Peptide* ptr;
};

PeptData::const_iterator::const_iterator(const Peptide& peptide) : pImpl(std::make_unique<Impl>(peptide)) {}
PeptData::const_iterator::const_iterator(const const_iterator& it) : pImpl(std::make_unique<Impl>(*it)) {}
PeptData::const_iterator::~const_iterator() {}

PeptData::const_iterator& PeptData::const_iterator::operator++() { pImpl->increment(); return *this; }
PeptData::const_iterator PeptData::const_iterator::operator++(int) {
    const_iterator old(pImpl->element()); ++(*this); return old;
}

bool PeptData::const_iterator::operator==(const_iterator other) const {
    return &(pImpl->element()) == &(other.pImpl->element());
}

bool PeptData::const_iterator::operator!=(const_iterator other) const {
    return &(pImpl->element()) != &(other.pImpl->element());
}

const PeptData::Peptide& PeptData::const_iterator::operator*() const { return pImpl->element(); }

// range object implementation
PeptData::range::range(const_iterator begin, const_iterator end) : begin_(begin), end_(end) {}
PeptData::const_iterator PeptData::range::begin() const { return begin_; }
PeptData::const_iterator PeptData::range::end() const { return end_; }

// TODO: Rewrite the rest of the code
// container implementation
class PeptData::Impl {
public:
    Impl(const char* filename, bool append_decoy, EnzymeType enzyme_type,
         unsigned max_miss_cleavage, double min_mass, double max_mass)
            : database_name_(filename),
              // TODO: change to future version
              mgr_(filename, append_decoy, PeptGenerator::EnzymeType::Trypsin, max_miss_cleavage, min_mass, max_mass) {}

    // TODO: change the adapters
    const Peptide& lower_bound(double mass) const {
        auto peptides = mgr_.Retrieve();
        auto start = std::lower_bound(peptides.cbegin(), peptides.cend(), mass,
            [](const auto& peptide, const auto& val) {
                return peptide.mass < val;
            }
        );
        return *start;
    }

    const Peptide& upper_bound(double mass) const {
        auto peptides = mgr_.Retrieve();
        auto end = std::upper_bound(peptides.cbegin(), peptides.cend(), mass,
            [](const auto& peptide, const auto& val) {
                return val < peptide.mass;
            }
        );
        return *end;
    }

    const Peptide& begin() const { return *(mgr_.Retrieve().cbegin()); }
    const Peptide& end() const { return *(mgr_.Retrieve().cend()); }

    size_t size() const { return mgr_.Retrieve().size(); }

private:
    const char* const database_name_;

    // TODO: decompose mgr to two dissociated data structure
    PeptMgr mgr_;
};

PeptData::PeptData(const char* filename, bool append_decoy, EnzymeType enzyme_type,
                     unsigned max_miss_cleavage, double min_mass, double max_mass)
        : pImpl(std::make_unique<Impl>(filename, append_decoy, enzyme_type,
                                       max_miss_cleavage, min_mass, max_mass)) {}
PeptData::~PeptData() {}

// adapters
PeptData::const_iterator PeptData::lower_bound(double mass) const {
    return const_iterator(pImpl->lower_bound(mass));
}
PeptData::const_iterator PeptData::upper_bound(double mass) const {
    return const_iterator(pImpl->upper_bound(mass));
}
PeptData::range PeptData::LoopWithin(double lower_mass, double upper_mass) const {
    return range(lower_bound(lower_mass), upper_bound(upper_mass));
}

PeptData::const_iterator PeptData::begin() const { return const_iterator(pImpl->begin()); }
PeptData::const_iterator PeptData::end() const { return const_iterator(pImpl->end()); }

size_t PeptData::size() const { return pImpl->size(); }

