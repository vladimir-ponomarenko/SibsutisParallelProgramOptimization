#include <gtest/gtest.h>
#include "motif_finder.h"
#include "iupac_codes.h"
#include <span>

using namespace dna_motif;

class MotifFinderTest : public ::testing::Test {
protected:
    void SetUp() override {
        iupac_codes = &IUPACCodes::getInstance();
        motif_finder = std::make_unique<MotifFinder>(*iupac_codes);

        sequences = {
            ChIPSequence("seq1", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"),
            ChIPSequence("seq2", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
            ChIPSequence("seq3", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
            ChIPSequence("seq4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
            ChIPSequence("seq5", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC")
        };

        motifs = {
            Motif("ATGCATGC", 10.5, 20.3, 30.1),
            Motif("TTTTTTTT", 15.2, 25.4, 35.6),
            Motif("GGGGGGGG", 12.8, 22.1, 32.9),
            Motif("ATRCATGC", 8.0, 18.0, 28.0)  // Ambiguous motif (R = A/G)
        };
    }

    IUPACCodes* iupac_codes;
    std::unique_ptr<MotifFinder> motif_finder;
    std::vector<ChIPSequence> sequences;
    std::vector<Motif> motifs;
};

TEST_F(MotifFinderTest, FindSingleMotif) {
    // Test finding exact match motif
    auto result = motif_finder->findSingleMotif(std::span<const ChIPSequence>(sequences), motifs[0]);  // ATGCATGC

    EXPECT_EQ(result.motif_pattern, "ATGCATGC");
    EXPECT_EQ(result.match_count, 2);  // seq1 and seq5 have matches
    EXPECT_DOUBLE_EQ(result.frequency, 0.4);  // 2 out of 5 sequences

    // Test finding TTTT motif
    auto result2 = motif_finder->findSingleMotif(std::span<const ChIPSequence>(sequences), motifs[1]);  // TTTTTTTT

    EXPECT_EQ(result2.motif_pattern, "TTTTTTTT");
    EXPECT_EQ(result2.match_count, 1);  // Only seq2 has matches
    EXPECT_DOUBLE_EQ(result2.frequency, 0.2);  // 1 out of 5 sequences
}

TEST_F(MotifFinderTest, FindMotifInSequence) {
    // Test exact match
    auto matches = motif_finder->findMotifInSequence(sequences[0], motifs[0], 0);
    EXPECT_EQ(matches.size(), 9);  // ATGCATGC appears 9 times in seq1 (40 chars / 8 chars = 5, but with overlaps)

    // Test no match
    auto no_matches = motif_finder->findMotifInSequence(sequences[1], motifs[0], 1);
    EXPECT_EQ(no_matches.size(), 0);

    // Test ambiguous motif
    auto ambiguous_matches = motif_finder->findMotifInSequence(sequences[0], motifs[3], 0);
    EXPECT_EQ(ambiguous_matches.size(), 9);  // ATRCATGC matches ATGCATGC (same as exact match)
}

TEST_F(MotifFinderTest, FindMotifs) {
    auto results = motif_finder->findMotifs(std::span<const ChIPSequence>(sequences), std::span<const Motif>(motifs));

    EXPECT_EQ(results.size(), motifs.size());

    // Check first result (ATGCATGC)
    EXPECT_EQ(results[0].motif_pattern, "ATGCATGC");
    EXPECT_EQ(results[0].match_count, 2);
    EXPECT_DOUBLE_EQ(results[0].frequency, 0.4);

    // Check second result (TTTTTTTT)
    EXPECT_EQ(results[1].motif_pattern, "TTTTTTTT");
    EXPECT_EQ(results[1].match_count, 1);
    EXPECT_DOUBLE_EQ(results[1].frequency, 0.2);

    // Check third result (GGGGGGGG)
    EXPECT_EQ(results[2].motif_pattern, "GGGGGGGG");
    EXPECT_EQ(results[2].match_count, 1);
    EXPECT_DOUBLE_EQ(results[2].frequency, 0.2);

    // Check fourth result (ATRCATGC - ambiguous)
    EXPECT_EQ(results[3].motif_pattern, "ATRCATGC");
    EXPECT_EQ(results[3].match_count, 2);
    EXPECT_DOUBLE_EQ(results[3].frequency, 0.4);
}

TEST_F(MotifFinderTest, CalculateFrequency) {
    EXPECT_DOUBLE_EQ(MotifFinder::calculateFrequency(0, 10), 0.0);
    EXPECT_DOUBLE_EQ(MotifFinder::calculateFrequency(5, 10), 0.5);
    EXPECT_DOUBLE_EQ(MotifFinder::calculateFrequency(10, 10), 1.0);
    EXPECT_DOUBLE_EQ(MotifFinder::calculateFrequency(3, 7), 3.0/7.0);

    // Edge case: division by zero
    EXPECT_DOUBLE_EQ(MotifFinder::calculateFrequency(5, 0), 0.0);
}

TEST_F(MotifFinderTest, PerformanceStats) {
    // Run some operations to generate performance stats
    motif_finder->findMotifs(std::span<const ChIPSequence>(sequences), std::span<const Motif>(motifs));

    auto stats = motif_finder->getPerformanceStats();

    EXPECT_GT(stats.size(), 0);
    EXPECT_TRUE(stats.find("find_motifs_total") != stats.end());
    EXPECT_TRUE(stats.find("process_single_motif") != stats.end());
}

TEST_F(MotifFinderTest, EmptyInput) {
    // Test with empty sequences
    std::vector<ChIPSequence> empty_sequences;
    auto result = motif_finder->findSingleMotif(std::span<const ChIPSequence>(empty_sequences), motifs[0]);

    EXPECT_EQ(result.match_count, 0);
    EXPECT_DOUBLE_EQ(result.frequency, 0.0);

    // Test with empty motifs
    std::vector<Motif> empty_motifs;
    auto results = motif_finder->findMotifs(std::span<const ChIPSequence>(sequences), std::span<const Motif>(empty_motifs));

    EXPECT_EQ(results.size(), 0);
}

TEST_F(MotifFinderTest, ShortSequences) {
    // Test with sequences shorter than motif
    std::vector<ChIPSequence> short_sequences = {
        ChIPSequence("short1", "ATG"),
        ChIPSequence("short2", "ATGC")
    };

    auto result = motif_finder->findSingleMotif(std::span<const ChIPSequence>(short_sequences), motifs[0]);  // ATGCATGC (8 chars)

    EXPECT_EQ(result.match_count, 0);
    EXPECT_DOUBLE_EQ(result.frequency, 0.0);
}

TEST_F(MotifFinderTest, MotifMatchDetails) {
    auto result = motif_finder->findSingleMotif(std::span<const ChIPSequence>(sequences), motifs[0]);

    // Check that we have the correct number of matches
    EXPECT_EQ(result.matches.size(), 2);  // Two sequences with matches

    // Check that matches are from correct sequences
    std::set<size_t> sequence_indices;
    for (const auto& match : result.matches) {
        sequence_indices.insert(match.sequence_index);
    }

    EXPECT_TRUE(sequence_indices.find(0) != sequence_indices.end());  // seq1
    EXPECT_TRUE(sequence_indices.find(4) != sequence_indices.end());  // seq5
}
