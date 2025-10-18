#include <gtest/gtest.h>
#include <cstring>
#include <iostream>
#include <string>


int test_main(int argc, char* argv[]) {
    std::string chip_seq_file;
    std::string motifs_file;
    std::string output_file;
    int num_threads = 0;
    bool verbose = false;
    (void)verbose;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            return 0;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                num_threads = std::stoi(argv[++i]);
                if (num_threads <= 0) {
                    return 1;
                }
            } else {
                return 1;
            }
        } else if (arg[0] != '-') {
            if (chip_seq_file.empty()) {
                chip_seq_file = arg;
            } else if (motifs_file.empty()) {
                motifs_file = arg;
            } else if (output_file.empty()) {
                output_file = arg;
            }
        } else {
            return 1;
        }
    }

    if (chip_seq_file.empty() || motifs_file.empty()) {
        return 1;
    }

    return 0;
}

class MainTest : public ::testing::Test {
protected:
    void SetUp() override {

    }

    void TearDown() override {

    }
};

TEST_F(MainTest, HelpOption) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("--help")};
    int argc = 2;

    EXPECT_EQ(test_main(argc, argv), 0);
}

TEST_F(MainTest, InvalidOptions) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("--invalid-option")};
    int argc = 2;

    EXPECT_NE(test_main(argc, argv), 0);
}

TEST_F(MainTest, MissingArguments) {
    char* argv[] = {const_cast<char*>("program")};
    int argc = 1;

    EXPECT_NE(test_main(argc, argv), 0);
}

TEST_F(MainTest, ThreadOption) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("--threads"), const_cast<char*>("4"), const_cast<char*>("test1.fst"), const_cast<char*>("test2.mot")};
    int argc = 5;

    EXPECT_EQ(test_main(argc, argv), 0);
}

TEST_F(MainTest, VerboseOption) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("--verbose"), const_cast<char*>("test1.fst"), const_cast<char*>("test2.mot")};
    int argc = 4;

    EXPECT_EQ(test_main(argc, argv), 0);
}

TEST_F(MainTest, ValidArguments) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("test1.fst"), const_cast<char*>("test2.mot"), const_cast<char*>("output.txt")};
    int argc = 4;

    EXPECT_EQ(test_main(argc, argv), 0);
}

TEST_F(MainTest, InvalidThreadCount) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("--threads"), const_cast<char*>("0"), const_cast<char*>("test1.fst"), const_cast<char*>("test2.mot")};
    int argc = 5;

    EXPECT_NE(test_main(argc, argv), 0);
}

TEST_F(MainTest, NegativeThreadCount) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("--threads"), const_cast<char*>("-1"), const_cast<char*>("test1.fst"), const_cast<char*>("test2.mot")};
    int argc = 5;

    EXPECT_NE(test_main(argc, argv), 0);
}

TEST_F(MainTest, ShortOptions) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("-t"), const_cast<char*>("2"), const_cast<char*>("-v"), const_cast<char*>("test1.fst"), const_cast<char*>("test2.mot")};
    int argc = 6;

    EXPECT_EQ(test_main(argc, argv), 0);
}

TEST_F(MainTest, MultipleOptions) {
    char* argv[] = {const_cast<char*>("program"), const_cast<char*>("--threads"), const_cast<char*>("4"), const_cast<char*>("--verbose"), const_cast<char*>("test1.fst"), const_cast<char*>("test2.mot"), const_cast<char*>("output.txt")};
    int argc = 7;

    EXPECT_EQ(test_main(argc, argv), 0);
}
