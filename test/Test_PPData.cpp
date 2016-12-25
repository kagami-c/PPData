#include <PPData.h>
#include <ProtData.h>
#include <gtest/gtest.h>

TEST(Unittest_PPData, PPData_API) {
    PPData ppdata("uniprot-all.fasta", false, PPData::EnzymeType::Trypsin, 0, 700, 5000);
    EXPECT_EQ(554376, ppdata.size());
//    auto range = ppdata.LoopWithin(700, 5000);
//    EXPECT_EQ(554376, std::distance(range.begin(), range.end()));

    PPData ppdata2("uniprot-all.fasta", false, PPData::EnzymeType::Trypsin, 1, 700, 5000);
    EXPECT_EQ(1465360, ppdata2.size());
//    iters = peptset2.RetrieveMassRange(700, 5000);
//    EXPECT_EQ(1465360, std::distance(iters.first, iters.second));

    PPData ppdata3("uniprot-all.fasta", false, PPData::EnzymeType::Trypsin, 2, 700, 5000);
    EXPECT_EQ(2404590, ppdata3.size());
//    iters = peptset3.RetrieveMassRange(700, 5000);
//    EXPECT_EQ(2404590, std::distance(iters.first, iters.second));
}

TEST(Unittest_PPData, ProtData_API) {
    auto fasta = ProtData("uniprot-all.fasta", false);
    EXPECT_EQ(20198, fasta.size());

    auto decoy = ProtData("uniprot-all.fasta", true);
    EXPECT_EQ(20198 * 2, decoy.size());
}

TEST(Unittest_PPData, PeptNum) {
    PPData ppdata("C:\\Users\\Jiaan\\Desktop\\RandomDatabase\\random20000.fasta", true,
                  PPData::EnzymeType::Trypsin, 2, 1000, 5000);
    EXPECT_EQ(0, ppdata.size());
}