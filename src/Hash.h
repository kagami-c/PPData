#pragma once

#include "PPData.h"
//#include <city.h>
#include <string>
#include <cstring>

// TODO: use cityhash
namespace std {
    template<> struct hash<PPData::Peptide> {
        size_t operator()(const PPData::Peptide& p) const {
            return (hash<string>()(string(p.sequence, p.sequence_length)));
//            return CityHash32(p.sequence, p.sequence_length);
        }
    };
}

inline bool operator==(const PPData::Peptide& one, const PPData::Peptide& another) {
    if (one.sequence_length != another.sequence_length) { return false; }
    return (0 == strncmp(one.sequence, another.sequence, one.sequence_length));
}
