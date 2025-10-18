#include "parallel_processor.h"
#include <cstring>
#include <expected>
#include <filesystem>
#include <format>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <span>
#include <string>

using namespace dna_motif;

// Command line argument structure
struct CommandLineArgs {
  std::string chip_seq_file;
  std::string motifs_file;
  std::string output_file;
  int num_threads = 0;
  bool verbose = false;
  bool help = false;
};

// Error types for command line parsing
enum class ParseError {
  InvalidArgument,
  MissingRequired,
  InvalidValue,
  Unknown
};

using ParseResult = std::expected<CommandLineArgs, ParseError>;

void printUsage(std::string_view program_name) {
  std::cout << std::format(
      "Usage: {} [OPTIONS] <chip_seq_file> <motifs_file> [output_file]\n",
      program_name);
  std::cout << "\nOptions:\n";
  std::cout << "  -t, --threads <num>    Number of OpenMP threads per process "
               "(default: auto)\n";
  std::cout << "  -h, --help             Show this help message\n";
  std::cout << "  -v, --verbose          Enable verbose output\n";
  std::cout << "\nArguments:\n";
  std::cout << "  chip_seq_file          Path to ChIP-seq sequences file\n";
  std::cout << "  motifs_file            Path to motifs file\n";
  std::cout << "  output_file            Optional output file for results "
               "(default: stdout)\n";
  std::cout << "\nExample:\n";
  std::cout << std::format(
      "  mpirun -n 4 {} -t 8 sequences.fst motifs.mot results.txt\n",
      program_name);
}

// Parse command line arguments
ParseResult parseArguments(std::span<const char *> args) {
  CommandLineArgs result;
  for (size_t i = 1; i < args.size(); ++i) {
    std::string_view arg = args[i];

    if (arg == "-h" || arg == "--help") {
      result.help = true;
    } else if (arg == "-v" || arg == "--verbose") {
      result.verbose = true;
    } else if (arg == "-t" || arg == "--threads") {
      if (i + 1 < args.size()) {
        try {
          result.num_threads = std::stoi(std::string(args[++i]));
          if (result.num_threads <= 0) {
            return std::unexpected(ParseError::InvalidValue);
          }
        } catch (const std::exception &) {
          return std::unexpected(ParseError::InvalidValue);
        }
      } else {
        return std::unexpected(ParseError::InvalidArgument);
      }
    } else if (arg[0] != '-') {
      if (result.chip_seq_file.empty()) {
        result.chip_seq_file = arg;
      } else if (result.motifs_file.empty()) {
        result.motifs_file = arg;
      } else if (result.output_file.empty()) {
        result.output_file = arg;
      }
    } else {
      return std::unexpected(ParseError::Unknown);
    }
  }

  if (result.chip_seq_file.empty() || result.motifs_file.empty()) {
    return std::unexpected(ParseError::MissingRequired);
  }

  return result;
}

// Validate input files
bool validateInputFiles(const CommandLineArgs &args) {
  if (!std::filesystem::exists(args.chip_seq_file)) {
    std::cerr << std::format("Error: ChIP-seq file '{}' does not exist\n",
                             args.chip_seq_file);
    return false;
  }

  if (!std::filesystem::exists(args.motifs_file)) {
    std::cerr << std::format("Error: Motifs file '{}' does not exist\n",
                             args.motifs_file);
    return false;
  }

  return true;
}

// Print error message
void printError(ParseError error, std::string_view program_name) {
  switch (error) {
  case ParseError::InvalidArgument:
    std::cerr << "Error: Invalid argument\n";
    break;
  case ParseError::MissingRequired:
    std::cerr << "Error: Missing required arguments\n";
    break;
  case ParseError::InvalidValue:
    std::cerr << "Error: Invalid value for argument\n";
    break;
  case ParseError::Unknown:
    std::cerr << "Error: Unknown option\n";
    break;
  }
  printUsage(program_name);
}

int main(int argc, char *argv[]) {
  auto args_result = parseArguments(
      std::span<const char *>(const_cast<const char **>(argv), argc));

  if (!args_result) {
    printError(args_result.error(), argv[0]);
    return 1;
  }

  const auto &args = *args_result;

  if (args.help) {
    printUsage(argv[0]);
    return 0;
  }

  if (!validateInputFiles(args)) {
    return 1;
  }

  try {
    ParallelProcessor processor;
    if (!processor.initialize(argc, argv, args.num_threads)) {
      std::cerr << "Failed to initialize parallel processor\n";
      return 1;
    }

    auto results =
        processor.processMotifs(args.chip_seq_file, args.motifs_file);

    if (args.output_file.empty()) {
      processor.printResults(results);
    } else {
      processor.saveResults(results, args.output_file);
    }

    if (args.verbose) {
      auto stats = processor.getPerformanceStats();
      std::cout << "\n=== PERFORMANCE STATISTICS ===" << std::endl;
      for (const auto &[operation, time] : stats) {
        std::cout << std::format("{}: {:.4f} seconds\n", operation, time);
      }
    }

    processor.finalize();

  } catch (const std::exception &e) {
    std::cerr << std::format("Error: {}\n", e.what());
    return 1;
  }

  return 0;
}
