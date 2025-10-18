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
    // Test standard nucleotides
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('A'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('T'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('G'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('C'));

    // Test ambiguous codes
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('R'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('Y'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('S'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('W'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('K'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('M'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('B'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('D'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('H'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('V'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('N'));

    // Test case insensitive
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('a'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('r'));
    EXPECT_TRUE(iupac_codes->isValidIUPACCode('n'));

    // Test invalid codes
    EXPECT_FALSE(iupac_codes->isValidIUPACCode('X'));
    EXPECT_FALSE(iupac_codes->isValidIUPACCode('Z'));
    EXPECT_FALSE(iupac_codes->isValidIUPACCode('1'));
    EXPECT_FALSE(iupac_codes->isValidIUPACCode('@'));
}

TEST_F(IUPACCodesTest, NucleotideMatching) {
    // Test standard nucleotide matching
    EXPECT_TRUE(iupac_codes->matches('A', 'A'));
    EXPECT_TRUE(iupac_codes->matches('T', 'T'));
    EXPECT_TRUE(iupac_codes->matches('G', 'G'));
    EXPECT_TRUE(iupac_codes->matches('C', 'C'));

    // Test ambiguous matching
    EXPECT_TRUE(iupac_codes->matches('A', 'R'));  // A matches R (A/G)
    EXPECT_TRUE(iupac_codes->matches('G', 'R'));  // G matches R (A/G)
    EXPECT_FALSE(iupac_codes->matches('T', 'R')); // T doesn't match R (A/G)
    EXPECT_FALSE(iupac_codes->matches('C', 'R')); // C doesn't match R (A/G)

    EXPECT_TRUE(iupac_codes->matches('T', 'Y'));  // T matches Y (T/C)
    EXPECT_TRUE(iupac_codes->matches('C', 'Y'));  // C matches Y (T/C)
    EXPECT_FALSE(iupac_codes->matches('A', 'Y')); // A doesn't match Y (T/C)
    EXPECT_FALSE(iupac_codes->matches('G', 'Y')); // G doesn't match Y (T/C)

    // Test N (any nucleotide)
    EXPECT_TRUE(iupac_codes->matches('A', 'N'));
    EXPECT_TRUE(iupac_codes->matches('T', 'N'));
    EXPECT_TRUE(iupac_codes->matches('G', 'N'));
    EXPECT_TRUE(iupac_codes->matches('C', 'N'));

    // Test case insensitive
    EXPECT_TRUE(iupac_codes->matches('a', 'A'));
    EXPECT_TRUE(iupac_codes->matches('A', 'a'));
    EXPECT_TRUE(iupac_codes->matches('g', 'R'));

    // Test invalid inputs
    EXPECT_FALSE(iupac_codes->matches('X', 'A'));
    EXPECT_FALSE(iupac_codes->matches('A', 'X'));
}

TEST_F(IUPACCodesTest, MotifMatching) {
    // Test simple motif matching
    std::string sequence = "ATGCATGC";
    std::string motif = "ATGC";

    EXPECT_TRUE(iupac_codes->matchesMotif(sequence, motif, 0));
    EXPECT_TRUE(iupac_codes->matchesMotif(sequence, motif, 4));
    EXPECT_FALSE(iupac_codes->matchesMotif(sequence, motif, 1));
    EXPECT_FALSE(iupac_codes->matchesMotif(sequence, motif, 5));

    // Test ambiguous motif matching
    std::string ambiguous_motif = "ATRC";  // R = A/G
    EXPECT_TRUE(iupac_codes->matchesMotif(sequence, ambiguous_motif, 0));  // ATGC matches ATRC
    EXPECT_TRUE(iupac_codes->matchesMotif(sequence, ambiguous_motif, 4));  // ATGC matches ATRC at pos 4 too

    // Test boundary conditions
    EXPECT_FALSE(iupac_codes->matchesMotif(sequence, motif, 5)); // Would go out of bounds
    EXPECT_FALSE(iupac_codes->matchesMotif("ATG", motif, 0));    // Sequence too short
}

TEST_F(IUPACCodesTest, FindMotifMatches) {
    std::string sequence = "ATGCATGCATGC";
    std::string motif = "ATGC";

    auto matches = iupac_codes->findMotifMatches(sequence, motif);

    EXPECT_EQ(matches.size(), 3);
    EXPECT_EQ(matches[0], 0);
    EXPECT_EQ(matches[1], 4);
    EXPECT_EQ(matches[2], 8);

    // Test with ambiguous motif
    std::string ambiguous_motif = "ATRC";  // R = A/G
    auto ambiguous_matches = iupac_codes->findMotifMatches(sequence, ambiguous_motif);

    EXPECT_EQ(ambiguous_matches.size(), 3);
    EXPECT_EQ(ambiguous_matches[0], 0);  // ATGC matches ATRC
    EXPECT_EQ(ambiguous_matches[1], 4);  // ATGC matches ATRC
    EXPECT_EQ(ambiguous_matches[2], 8);  // ATGC matches ATRC

    // Test with no matches
    std::string no_match_motif = "TTTT";
    auto no_matches = iupac_codes->findMotifMatches(sequence, no_match_motif);
    EXPECT_EQ(no_matches.size(), 0);
}

TEST_F(IUPACCodesTest, GetNucleotides) {
    // Test standard nucleotides
    auto a_nucs = iupac_codes->getNucleotides('A');
    EXPECT_EQ(a_nucs.size(), 1);
    EXPECT_TRUE(std::ranges::find(a_nucs, 'A') != a_nucs.end());

    // Test ambiguous nucleotides
    auto r_nucs = iupac_codes->getNucleotides('R');
    EXPECT_EQ(r_nucs.size(), 2);
    EXPECT_TRUE(std::ranges::find(r_nucs, 'A') != r_nucs.end());
    EXPECT_TRUE(std::ranges::find(r_nucs, 'G') != r_nucs.end());

    auto n_nucs = iupac_codes->getNucleotides('N');
    EXPECT_EQ(n_nucs.size(), 4);
    EXPECT_TRUE(std::ranges::find(n_nucs, 'A') != n_nucs.end());
    EXPECT_TRUE(std::ranges::find(n_nucs, 'T') != n_nucs.end());
    EXPECT_TRUE(std::ranges::find(n_nucs, 'G') != n_nucs.end());
    EXPECT_TRUE(std::ranges::find(n_nucs, 'C') != n_nucs.end());

    // Test invalid code
    auto invalid_nucs = iupac_codes->getNucleotides('X');
    EXPECT_EQ(invalid_nucs.size(), 0);
}
