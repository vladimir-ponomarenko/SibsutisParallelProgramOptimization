#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include "dna_parser.h"

using namespace dna_motif;

class DNAParserTest : public ::testing::Test {
protected:
    void SetUp() override {
        parser = std::make_unique<DNAParser>();

        createTestChIPFile();
        createTestMotifsFile();
    }

    void TearDown() override {
        std::remove("test_chip.fst");
        std::remove("test_motifs.mot");
    }

    void createTestChIPFile() {
        std::ofstream file("test_chip.fst");
        file << ">seq1\tmetadata1\tmetadata2\n";
        file << "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
        file << ">seq2\tmetadata3\n";
        file << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";
        file << ">seq3\n";
        file << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n";
        file << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";
        file.close();
    }

    void createTestMotifsFile() {
        std::ofstream file("test_motifs.mot");
        file << "ATGCATGC\t10.5\t20.3\t30.1\n";
        file << "TTTTTTTT\t15.2\t25.4\t35.6\n";
        file << "# This is a comment\n";
        file << "GGGGGGGG\t12.8\t22.1\t32.9\n";
        file << "\n";
        file.close();
    }

    std::unique_ptr<DNAParser> parser;
};

TEST_F(DNAParserTest, ParseChIPSequences) {
    auto sequences_result = parser->parseChIPSequences("test_chip.fst");

    ASSERT_TRUE(sequences_result.has_value());
    const auto& sequences = sequences_result.value();

    EXPECT_EQ(sequences.size(), 3);

    // Check first sequence
    EXPECT_EQ(sequences[0].id, "seq1");
    EXPECT_EQ(sequences[0].sequence, "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC");
    EXPECT_EQ(sequences[0].metadata.size(), 2);
    EXPECT_EQ(sequences[0].metadata[0], "metadata1");
    EXPECT_EQ(sequences[0].metadata[1], "metadata2");

    // Check second sequence
    EXPECT_EQ(sequences[1].id, "seq2");
    EXPECT_EQ(sequences[1].sequence, "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
    EXPECT_EQ(sequences[1].metadata.size(), 1);
    EXPECT_EQ(sequences[1].metadata[0], "metadata3");

    // Check third sequence (with multi-line sequence)
    EXPECT_EQ(sequences[2].id, "seq3");
    EXPECT_EQ(sequences[2].sequence, "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
    EXPECT_EQ(sequences[2].metadata.size(), 0);
}

TEST_F(DNAParserTest, ParseMotifs) {
    auto motifs_result = parser->parseMotifs("test_motifs.mot");

    ASSERT_TRUE(motifs_result.has_value());
    const auto& motifs = motifs_result.value();

    EXPECT_EQ(motifs.size(), 3);

    // Check first motif
    EXPECT_EQ(motifs[0].pattern, "ATGCATGC");
    EXPECT_DOUBLE_EQ(motifs[0].score1, 10.5);
    EXPECT_DOUBLE_EQ(motifs[0].score2, 20.3);
    EXPECT_DOUBLE_EQ(motifs[0].score3, 30.1);

    // Check second motif
    EXPECT_EQ(motifs[1].pattern, "TTTTTTTT");
    EXPECT_DOUBLE_EQ(motifs[1].score1, 15.2);
    EXPECT_DOUBLE_EQ(motifs[1].score2, 25.4);
    EXPECT_DOUBLE_EQ(motifs[1].score3, 35.6);

    // Check third motif
    EXPECT_EQ(motifs[2].pattern, "GGGGGGGG");
    EXPECT_DOUBLE_EQ(motifs[2].score1, 12.8);
    EXPECT_DOUBLE_EQ(motifs[2].score2, 22.1);
    EXPECT_DOUBLE_EQ(motifs[2].score3, 32.9);
}

TEST_F(DNAParserTest, ValidateSequence) {
    // Valid sequences
    EXPECT_TRUE(parser->validateSequence("ATGC"));
    EXPECT_TRUE(parser->validateSequence("ATGCATGC"));
    EXPECT_TRUE(parser->validateSequence("AAAA"));
    EXPECT_TRUE(parser->validateSequence("TTTT"));
    EXPECT_TRUE(parser->validateSequence("GGGG"));
    EXPECT_TRUE(parser->validateSequence("CCCC"));

    // Case insensitive
    EXPECT_TRUE(parser->validateSequence("atgc"));
    EXPECT_TRUE(parser->validateSequence("AtGc"));

    // Invalid sequences
    EXPECT_FALSE(parser->validateSequence(""));
    EXPECT_FALSE(parser->validateSequence("ATGCX"));
    EXPECT_FALSE(parser->validateSequence("ATGC1"));
    EXPECT_FALSE(parser->validateSequence("ATG C"));
    EXPECT_FALSE(parser->validateSequence("ATGC@"));
}

TEST_F(DNAParserTest, ParseStatistics) {
    // Parse files to generate statistics
    auto sequences_result = parser->parseChIPSequences("test_chip.fst");
    auto motifs_result = parser->parseMotifs("test_motifs.mot");

    ASSERT_TRUE(sequences_result.has_value());
    ASSERT_TRUE(motifs_result.has_value());

    auto stats = parser->getStatistics();

    EXPECT_GT(stats["files_opened"], 0);
    EXPECT_GT(stats["files_closed"], 0);
    EXPECT_GT(stats["sequences_parsed"], 0);
    EXPECT_GT(stats["motifs_parsed"], 0);
}

TEST_F(DNAParserTest, FileNotFound) {
    auto sequences_result = parser->parseChIPSequences("nonexistent.fst");
    auto motifs_result = parser->parseMotifs("nonexistent.mot");

    EXPECT_FALSE(sequences_result.has_value());
    EXPECT_FALSE(motifs_result.has_value());
    EXPECT_EQ(sequences_result.error(), ParseError::FileNotFound);
    EXPECT_EQ(motifs_result.error(), ParseError::FileNotFound);
}

TEST_F(DNAParserTest, InvalidMotifFormat) {
    // Create file with invalid motif format
    std::ofstream file("invalid_motifs.mot");
    file << "ATGC\t10.5\n";  // Missing scores
    file << "TTTT\t15.2\t25.4\n";  // Missing one score
    file.close();

    auto motifs_result = parser->parseMotifs("invalid_motifs.mot");

    // Should parse valid motifs and skip invalid ones
    ASSERT_TRUE(motifs_result.has_value());
    const auto& motifs = motifs_result.value();
    EXPECT_EQ(motifs.size(), 0);  // All lines are invalid

    std::remove("invalid_motifs.mot");
}

TEST_F(DNAParserTest, EmptyFiles) {
    // Create empty files
    std::ofstream file1("empty_chip.fst");
    file1.close();

    std::ofstream file2("empty_motifs.mot");
    file2.close();

    auto sequences_result = parser->parseChIPSequences("empty_chip.fst");
    auto motifs_result = parser->parseMotifs("empty_motifs.mot");

    ASSERT_TRUE(sequences_result.has_value());
    ASSERT_TRUE(motifs_result.has_value());

    const auto& sequences = sequences_result.value();
    const auto& motifs = motifs_result.value();

    EXPECT_EQ(sequences.size(), 0);
    EXPECT_EQ(motifs.size(), 0);

    std::remove("empty_chip.fst");
    std::remove("empty_motifs.mot");
}
