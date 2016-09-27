// Copyright (C) 2016

#pragma once

#include <memory>

class PPData {
public:
    struct Protein {
        const char* const name;
        const char* const sequence;
        const size_t sequence_length;

        Protein(const char* name, const char* sequence, size_t sequence_length);
    };

    struct Peptide {
        const char* sequence;
        size_t sequence_length;

        char n_term;
        char c_term;
        double mass;  // nonconst for adjusting modified amino acids

        const Protein* protein;
        size_t offset;  // offset in protein sequence

        Peptide(const Protein& protein, const char* compact_protein_sequence, 
                size_t start_idx, size_t end_idx, double mass);
    };

    enum class EnzymeType { Trypsin };

    // a simple range object for range based for loop
    template <typename iterator_type> 
    class range {  
    public:
        range(iterator_type begin, iterator_type end) : begin_(begin), end_(end) {}
        auto begin() const { return begin_; }
        auto end() const { return end_; }

    private:
        const iterator_type begin_;
        const iterator_type end_;
    };

    // ctors
    PPData(const char* filename, bool append_decoy, EnzymeType enzyme_type, 
           unsigned max_miss_cleavage, double min_mass, double max_mass);
    PPData(const char* filename) 
           : PPData(filename, false, EnzymeType::Trypsin, 0, 600.0, 5000.0) {}
    ~PPData();

    // access methods
    size_t size() const;
    const Peptide& operator[](const size_t index) const;

    auto begin() const;
    auto end() const;

    // subrange loop API
    auto lower_bound(double mass) const;  // including
    auto upper_bound(double mass) const;  // excluding
    auto LoopWithin(double lower_mass, double upper_mass) const;

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};
