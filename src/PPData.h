// Copyright (C) 2016

#pragma once

#include <iterator>
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
        const char* const sequence;
        const size_t sequence_length;

        const char n_term;
        const char c_term;
        const double mass;

        const Protein* const protein;
        const size_t offset;  // offset in protein sequence

        Peptide(const Protein& protein, const char* compact_protein_sequence, 
                size_t start_idx, size_t end_idx, double mass);
    };

    // iterator interface
    class const_iterator : public std::iterator<std::input_iterator_tag, Peptide> {
    public:
        explicit const_iterator(const Peptide& peptide);
        const_iterator(const const_iterator& it);
        ~const_iterator();
        const_iterator& operator++();
        const_iterator operator++(int);
        bool operator==(const_iterator other) const;
        bool operator!=(const_iterator other) const;
        const Peptide& operator*() const;

    private:
        class Impl;
        std::unique_ptr<Impl> pImpl;
    };

    // a simple range object for range based for loop
    class range {
    public:
        range(const_iterator begin, const_iterator end);
        const_iterator begin() const;
        const_iterator end() const;

    private:
        const const_iterator begin_;
        const const_iterator end_;
    };

    enum class EnzymeType { Trypsin };

    // ctors
    PPData(const char* filename,
              bool append_decoy,
              EnzymeType enzyme_type,
              unsigned max_miss_cleavage,
              double min_mass,
              double max_mass);
    PPData(const char* filename) : PPData(filename, false, EnzymeType::Trypsin, 0, 600.0, 5000.0) {}
    ~PPData();

    size_t size() const;
    const_iterator begin() const;
    const_iterator end() const;

    // subrange loop API
    const_iterator lower_bound(double mass) const;  // including
    const_iterator upper_bound(double mass) const;  // excluding
    range LoopWithin(double lower_mass, double upper_mass) const;

private:
    class Impl;
    std::unique_ptr<Impl> pImpl;
};
