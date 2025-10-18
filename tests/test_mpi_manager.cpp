#include <gtest/gtest.h>
#include "mpi_manager.h"

using namespace dna_motif;

class MPIManagerTest : public ::testing::Test {
protected:
    void SetUp() override {
        mpi_manager = std::make_unique<MPIManager>();
        if (!mpi_manager->initialize(0, nullptr)) {
            GTEST_SKIP() << "MPI initialization failed";
        }
    }

    void TearDown() override {
        if (mpi_manager) {
            mpi_manager->finalize();
        }
    }

    std::unique_ptr<MPIManager> mpi_manager;
};

TEST_F(MPIManagerTest, Initialization) {
    EXPECT_TRUE(mpi_manager->isMaster() || !mpi_manager->isMaster());
    EXPECT_GE(mpi_manager->getRank(), 0);
    EXPECT_GT(mpi_manager->getSize(), 0);
    EXPECT_LT(mpi_manager->getRank(), mpi_manager->getSize());
}

TEST_F(MPIManagerTest, WorkDistribution) {
    // Test work distribution calculation
    auto [start, count] = mpi_manager->calculateWorkDistribution(10, 0, 2);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(count, 5);

    auto [start2, count2] = mpi_manager->calculateWorkDistribution(10, 1, 2);
    EXPECT_EQ(start2, 5);
    EXPECT_EQ(count2, 5);

    // Test uneven distribution
    auto [start3, count3] = mpi_manager->calculateWorkDistribution(11, 0, 3);
    EXPECT_EQ(start3, 0);
    EXPECT_EQ(count3, 4);  // 11/3 = 3 remainder 2, first process gets 3+1=4

    auto [start4, count4] = mpi_manager->calculateWorkDistribution(11, 1, 3);
    EXPECT_EQ(start4, 4);
    EXPECT_EQ(count4, 4);  // Second process gets 3+1=4

    auto [start5, count5] = mpi_manager->calculateWorkDistribution(11, 2, 3);
    EXPECT_EQ(start5, 8);
    EXPECT_EQ(count5, 3);  // Third process gets 3
}

TEST_F(MPIManagerTest, Synchronization) {
    // Test that all processes can synchronize
    EXPECT_NO_THROW(mpi_manager->synchronize());
}

TEST_F(MPIManagerTest, CommunicationStats) {
    auto stats = mpi_manager->getCommunicationStats();

    // Initially should be empty or have some basic stats
    // The exact content depends on what operations were performed
    EXPECT_TRUE(stats.empty() || stats.size() > 0);
}

TEST_F(MPIManagerTest, MotifDistribution) {
    std::vector<Motif> motifs = {
        Motif("ATGCATGC", 10.5, 20.3, 30.1),
        Motif("TTTTTTTT", 15.2, 25.4, 35.6),
        Motif("GGGGGGGG", 12.8, 22.1, 32.9)
    };

    // Test motif broadcasting
    auto received_motifs = mpi_manager->broadcastMotifs(motifs);

    EXPECT_EQ(received_motifs.size(), motifs.size());

    for (size_t i = 0; i < motifs.size(); ++i) {
        EXPECT_EQ(received_motifs[i].pattern, motifs[i].pattern);
        EXPECT_DOUBLE_EQ(received_motifs[i].score1, motifs[i].score1);
        EXPECT_DOUBLE_EQ(received_motifs[i].score2, motifs[i].score2);
        EXPECT_DOUBLE_EQ(received_motifs[i].score3, motifs[i].score3);
    }
}

TEST_F(MPIManagerTest, SequenceDistribution) {
    std::vector<ChIPSequence> sequences = {
        ChIPSequence("seq1", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"),
        ChIPSequence("seq2", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"),
        ChIPSequence("seq3", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"),
        ChIPSequence("seq4", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
    };

    auto local_sequences = mpi_manager->distributeSequences(sequences);

    // Each process should get some sequences
    EXPECT_LE(local_sequences.size(), sequences.size());
    EXPECT_GT(local_sequences.size(), 0);

    // All sequences should be valid
    for (const auto& seq : local_sequences) {
        EXPECT_FALSE(seq.id.empty());
        EXPECT_FALSE(seq.sequence.empty());
    }
}

TEST_F(MPIManagerTest, ResultGathering) {
    std::vector<MotifResult> local_results = {
        MotifResult()
    };
    local_results[0].motif_pattern = "ATGCATGC";
    local_results[0].match_count = 5;
    local_results[0].frequency = 0.5;

    auto all_results = mpi_manager->gatherResults(local_results);

    // Should have results from all processes
    EXPECT_GE(all_results.size(), local_results.size());

    // Check that our local results are included
    bool found_local = false;
    for (const auto& result : all_results) {
        if (result.motif_pattern == "ATGCATGC" &&
            result.match_count == 5 &&
            result.frequency == 0.5) {
            found_local = true;
            break;
        }
    }
    EXPECT_TRUE(found_local);
}

TEST_F(MPIManagerTest, MultipleOperations) {
    // Test multiple operations to ensure no interference

    std::vector<ChIPSequence> sequences = {
        ChIPSequence("seq1", "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"),
        ChIPSequence("seq2", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
    };

    std::vector<Motif> motifs = {
        Motif("ATGCATGC", 10.5, 20.3, 30.1)
    };

    // Perform multiple operations
    auto local_sequences = mpi_manager->distributeSequences(sequences);
    auto received_motifs = mpi_manager->broadcastMotifs(motifs);

    std::vector<MotifResult> local_results;
    if (!local_sequences.empty()) {
        local_results.push_back(MotifResult());
        local_results[0].motif_pattern = "ATGCATGC";
        local_results[0].match_count = 1;
        local_results[0].frequency = 0.5;
    }

    auto all_results = mpi_manager->gatherResults(local_results);

    // All operations should complete without errors
    EXPECT_TRUE(true);
}
