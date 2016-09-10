#include "PeptMgr.h"
#include <gtest/gtest.h>

TEST(Unittest_PeptMgr, TrypsinDigest) {
    auto peptset1 = PeptMgr("uniprot-all.fasta", false, PeptGenerator::EnzymeType::Trypsin, 0);
    auto peptides = peptset1.Retrieve();
    EXPECT_EQ(554376, peptides.size());
    auto iters = peptset1.RetrieveMassRange(700, 5000);
    EXPECT_EQ(554376, std::distance(iters.first, iters.second));

    auto peptset2 = PeptMgr("uniprot-all.fasta", false, PeptGenerator::EnzymeType::Trypsin, 1);
    peptides = peptset2.Retrieve();
    EXPECT_EQ(1465360, peptides.size());
    iters = peptset2.RetrieveMassRange(700, 5000);
    EXPECT_EQ(1465360, std::distance(iters.first, iters.second));

    auto peptset3 = PeptMgr("uniprot-all.fasta", false, PeptGenerator::EnzymeType::Trypsin, 2);
    peptides = peptset3.Retrieve();
    EXPECT_EQ(2404590, peptides.size());
    iters = peptset3.RetrieveMassRange(700, 5000);
    EXPECT_EQ(2404590, std::distance(iters.first, iters.second));
}