#include <gtest/gtest.h>
#include <fstream>
#include "parallel_processor.h"

using namespace dna_motif;

class ParallelProcessorTest : public ::testing::Test {
protected:
    void SetUp() override {
        processor = std::make_unique<ParallelProcessor>();

        // Create test files
        createTestChIPFile();
        createTestMotifsFile();
    }

    void TearDown() override {
        // Clean up test files
        std::remove("test_chip_parallel.fst");
        std::remove("test_motifs_parallel.mot");
        std::remove("test_output.txt");

        // Finalize processor if initialized
        if (processor) {
            processor->finalize();
        }
    }

    void createTestChIPFile() {
        std::ofstream file("test_chip_parallel.fst");
        file << ">seq1\tmetadata1\n";
        file << "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
        file << ">seq2\tmetadata2\n";
        file << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";
        file << ">seq3\n";
        file << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n";
        file << ">seq4\n";
        file << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";
        file << ">seq5\n";
        file << "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n";
        file.close();
    }

    void createTestMotifsFile() {
        std::ofstream file("test_motifs_parallel.mot");
        file << "ATGCATGC\t10.5\t20.3\t30.1\n";
        file << "TTTTTTTT\t15.2\t25.4\t35.6\n";
        file << "GGGGGGGG\t12.8\t22.1\t32.9\n";
        file << "ATRCATGC\t8.0\t18.0\t28.0\n";  // Ambiguous motif
        file.close();
    }

    std::unique_ptr<ParallelProcessor> processor;
};

TEST_F(ParallelProcessorTest, Initialization) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));  // 2 OpenMP threads
    processor->finalize();
}

TEST_F(ParallelProcessorTest, ProcessMotifs) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));

    auto results = processor->processMotifs("test_chip_parallel.fst", "test_motifs_parallel.mot");

    // Should have results for all motifs
    EXPECT_EQ(results.size(), 4);

    // Check that results are valid
    for (const auto& result : results) {
        EXPECT_FALSE(result.motif_pattern.empty());
        EXPECT_GE(result.match_count, 0);
        EXPECT_GE(result.frequency, 0.0);
        EXPECT_LE(result.frequency, 1.0);
    }

    processor->finalize();
}

TEST_F(ParallelProcessorTest, SaveResults) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));

    auto results = processor->processMotifs("test_chip_parallel.fst", "test_motifs_parallel.mot");

    // Save results to file
    processor->saveResults(results, "test_output.txt");

    // Check that file was created and contains data
    std::ifstream file("test_output.txt");
    EXPECT_TRUE(file.is_open());

    std::string line;
    int line_count = 0;
    while (std::getline(file, line)) {
        line_count++;
    }

    // Should have header + 4 motif results
    EXPECT_EQ(line_count, 5);

    file.close();
    processor->finalize();
}

TEST_F(ParallelProcessorTest, PerformanceStats) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));

    processor->processMotifs("test_chip_parallel.fst", "test_motifs_parallel.mot");

    auto stats = processor->getPerformanceStats();

    // Should have some performance statistics
    EXPECT_GT(stats.size(), 0);

    // Check for expected performance metrics
    EXPECT_TRUE(stats.find("total_processing_time") != stats.end());
    EXPECT_TRUE(stats.find("file_loading_time") != stats.end());
    EXPECT_TRUE(stats.find("parallel_processing_time") != stats.end());

    processor->finalize();
}

TEST_F(ParallelProcessorTest, InvalidFiles) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));

    // Test with non-existent files
    EXPECT_THROW(processor->processMotifs("nonexistent.fst", "test_motifs_parallel.mot"), std::exception);
    EXPECT_THROW(processor->processMotifs("test_chip_parallel.fst", "nonexistent.mot"), std::exception);

    processor->finalize();
}

TEST_F(ParallelProcessorTest, EmptyFiles) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));

    // Create empty files
    std::ofstream file1("empty_chip.fst");
    file1.close();

    std::ofstream file2("empty_motifs.mot");
    file2.close();

    auto results = processor->processMotifs("empty_chip.fst", "empty_motifs.mot");

    // Should handle empty files gracefully
    EXPECT_EQ(results.size(), 0);

    std::remove("empty_chip.fst");
    std::remove("empty_motifs.mot");

    processor->finalize();
}

TEST_F(ParallelProcessorTest, DifferentThreadCounts) {
    // Test with different OpenMP thread counts
    std::vector<int> thread_counts = {1, 2, 4};

    for (int threads : thread_counts) {
        EXPECT_TRUE(processor->initialize(0, nullptr, threads));

        auto results = processor->processMotifs("test_chip_parallel.fst", "test_motifs_parallel.mot");

        // Results should be consistent regardless of thread count
        EXPECT_EQ(results.size(), 4);

        processor->finalize();
    }
}

TEST_F(ParallelProcessorTest, PrintResults) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));

    auto results = processor->processMotifs("test_chip_parallel.fst", "test_motifs_parallel.mot");

    // Test that printResults doesn't crash
    EXPECT_NO_THROW(processor->printResults(results));

    processor->finalize();
}

TEST_F(ParallelProcessorTest, MultipleRuns) {
    EXPECT_TRUE(processor->initialize(0, nullptr, 2));

    // Run multiple times to ensure no state issues
    for (int i = 0; i < 3; ++i) {
        auto results = processor->processMotifs("test_chip_parallel.fst", "test_motifs_parallel.mot");
        EXPECT_EQ(results.size(), 4);
    }

    processor->finalize();
}
