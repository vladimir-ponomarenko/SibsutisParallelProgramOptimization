#include <gtest/gtest.h>
#include "mpi_manager.h"

using namespace dna_motif;

class MPIManagerSimpleTest : public ::testing::Test {
protected:
    void SetUp() override {

    }

    void TearDown() override {

    }
};

TEST_F(MPIManagerSimpleTest, WorkDistributionCalculation) {
    // Test work distribution calculation without MPI initialization
    MPIManager manager;

    // Test even distribution
    auto [start, count] = manager.calculateWorkDistribution(10, 0, 2);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(count, 5);

    auto [start2, count2] = manager.calculateWorkDistribution(10, 1, 2);
    EXPECT_EQ(start2, 5);
    EXPECT_EQ(count2, 5);

    // Test uneven distribution
    auto [start3, count3] = manager.calculateWorkDistribution(11, 0, 3);
    EXPECT_EQ(start3, 0);
    EXPECT_EQ(count3, 4);  // 11/3 = 3 remainder 2, first process gets 3+1=4

    auto [start4, count4] = manager.calculateWorkDistribution(11, 1, 3);
    EXPECT_EQ(start4, 4);
    EXPECT_EQ(count4, 4);  // Second process gets 3+1=4

    auto [start5, count5] = manager.calculateWorkDistribution(11, 2, 3);
    EXPECT_EQ(start5, 8);
    EXPECT_EQ(count5, 3);  // Third process gets 3
}

TEST_F(MPIManagerSimpleTest, EdgeCases) {
    MPIManager manager;

    // Test with single process
    auto [start, count] = manager.calculateWorkDistribution(10, 0, 1);
    EXPECT_EQ(start, 0);
    EXPECT_EQ(count, 10);

    // Test with more processes than work
    auto [start2, count2] = manager.calculateWorkDistribution(3, 5, 10);
    EXPECT_EQ(start2, 3);  // Process 5 gets work starting at position 3
    EXPECT_EQ(count2, 0);  // But no work items

    // Test with zero work
    auto [start3, count3] = manager.calculateWorkDistribution(0, 0, 2);
    EXPECT_EQ(start3, 0);
    EXPECT_EQ(count3, 0);
}
