#include <gtest/gtest.h>
#include "iupac_codes.h"

using namespace dna_motif;

class IUPACCodesTest : public ::testing::Test {
protected:
    void SetUp() override {
        iupac_codes = &IUPACCodes::getInstance();
    }

    IUPACCodes* iupac_codes;
};

TEST_F(IUPACCodesTest, ValidIUPACCodes) {
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('A'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('N'));
    EXPECT_FALSE(iupac_codes->isValidIUPACCode('X'));
}

TEST_F(IUPACCodesTest, BitmasksCorrectness) {
    // Check Table 12 values
    // A: 0001 (1)
    EXPECT_EQ(iupac_codes->getBitmask('A'), 1);
    // T: 0010 (2)
    EXPECT_EQ(iupac_codes->getBitmask('T'), 2);
    // G: 0100 (4)
    EXPECT_EQ(iupac_codes->getBitmask('G'), 4);
    // C: 1000 (8)
    EXPECT_EQ(iupac_codes->getBitmask('C'), 8);
    
    // R (A/G): 0101 (5)
    EXPECT_EQ(iupac_codes->getBitmask('R'), 5);
    // N: 1111 (15)
    EXPECT_EQ(iupac_codes->getBitmask('N'), 15);
}

TEST_F(IUPACCodesTest, Hashing) {
    // Test Hash Generation
    // Pattern: AAAA TTTT
    // Hash: (1<<0)|(1<<4)|... 
    
    std::string s1 = "AAAAAAAA";
    uint32_t h1 = iupac_codes->hashSequence(s1);
    // 0x11111111
    EXPECT_EQ(h1, 0x11111111);
    
    std::string s2 = "TTTTTTTT";
    uint32_t h2 = iupac_codes->hashSequence(s2);
    // 0x22222222
    EXPECT_EQ(h2, 0x22222222);
}

TEST_F(IUPACCodesTest, MatchesBitwise) {
    // Motif: R (A/G)
    // Seq: A
    // Mask R: 0101, Mask A: 0001
    // (0001 & 0101) == 0001 -> Match
    EXPECT_TRUE(iupac_codes->matches('A', 'R'));
    
    // Seq: T (0010)
    // (0010 & 0101) == 0000 != 0010 -> No Match
    EXPECT_FALSE(iupac_codes->matches('T', 'R'));
}

TEST_F(IUPACCodesTest, MotifMatching) {
    std::string sequence = "ATGCATGC";
    std::string motif = "ATGC";

    EXPECT_TRUE(iupac_codes->matchesMotif(sequence, motif, 0));
    EXPECT_TRUE(iupac_codes->matchesMotif(sequence, motif, 4));
    
    // Ambiguous motif ATRC (R=A/G) should match ATGC
    std::string ambiguous = "ATRC";
    EXPECT_TRUE(iupac_codes->matchesMotif(sequence, ambiguous, 0));
}

TEST_F(IUPACCodesTest, FindMotifMatches) {
    std::string sequence = "ATGCATGCATGC";
    std::string motif = "ATGC";

    auto matches = iupac_codes->findMotifMatches(sequence, motif);

    EXPECT_EQ(matches.size(), 3);
    EXPECT_EQ(matches[0], 0);
    EXPECT_EQ(matches[1], 4);
    EXPECT_EQ(matches[2], 8);
}