#pragma once

#include "common.h"
#include "concepts.h"
#include <expected>
#include <filesystem>
#include <format>
#include <ranges>

namespace dna_motif {

enum class ParseError {
  FileNotFound,
  InvalidFormat,
  InvalidSequence,
  InvalidMotif,
  IOError,
  Unknown
};

template <typename T> using ParseResult = std::expected<T, ParseError>;

/**
 * @brief Parser for ChIP-seq data files
 *
 * Parses input files containing DNA sequences in the specified format:
 * >id    metadata...
 * SEQUENCE_LINE_1
 * SEQUENCE_LINE_2
 * ...
 */
class DNAParser {
public:
  DNAParser() = default;
  ~DNAParser() = default;

  DNAParser(const DNAParser &) = default;
  DNAParser &operator=(const DNAParser &) = default;
  DNAParser(DNAParser &&) = default;
  DNAParser &operator=(DNAParser &&) = default;

  /**
   * @brief Parse ChIP-seq sequences from file
   * @param filename Path to input file
   * @return Expected vector of parsed ChIP sequences or error
   */
  [[nodiscard]] ParseResult<std::vector<ChIPSequence>>
  parseChIPSequences(std::string_view filename);

  /**
   * @brief Parse motifs from file
   * @param filename Path to motifs file
   * @return Expected vector of parsed motifs or error
   */
  [[nodiscard]] ParseResult<std::vector<Motif>>
  parseMotifs(std::string_view filename);

  /**
   * @brief Validate a DNA sequence
   * @param sequence Sequence to validate
   * @return true if sequence is valid
   */
  [[nodiscard]] bool validateSequence(std::string_view sequence) const noexcept;

  /**
   * @brief Get parsing statistics
   * @return Map with parsing statistics
   */
  [[nodiscard]] const std::unordered_map<std::string, size_t> &
  getStatistics() const noexcept {
    return stats_;
  }

  /**
   * @brief Reset parsing statistics
   */
  void resetStatistics() noexcept { stats_.clear(); }

  /**
   * @brief Check if file exists and is readable
   * @param filename Path to file
   * @return true if file exists and is readable
   */
  [[nodiscard]] bool isFileReadable(std::string_view filename) const noexcept;

  /**
   * @brief Get file size
   * @param filename Path to file
   * @return File size in bytes or 0 if error
   */
  [[nodiscard]] size_t getFileSize(std::string_view filename) const noexcept;

private:
  std::unordered_map<std::string, size_t> stats_;

  /**
   * @brief Parse a single ChIP sequence from lines
   * @param header_line Header line starting with '>'
   * @param sequence_lines Lines containing the DNA sequence
   * @return Parsed ChIP sequence
   */
  [[nodiscard]] ChIPSequence
  parseChIPSequence(std::string_view header_line,
                    std::span<const std::string> sequence_lines);

  /**
   * @brief Clean and concatenate sequence lines
   * @param sequence_lines Raw sequence lines
   * @return Cleaned concatenated sequence
   */
  [[nodiscard]] std::string
  cleanSequence(std::span<const std::string> sequence_lines) const;

  /**
   * @brief Parse motif line
   * @param line Line containing motif data
   * @return Parsed motif
   */
  [[nodiscard]] Motif parseMotifLine(std::string_view line);

  /**
   * @brief Update parsing statistics
   * @param key Statistic key
   * @param increment Value to add
   */
  void updateStats(std::string_view key, size_t increment = 1) noexcept;

  /**
   * @brief Read file content into string
   * @param filename Path to file
   * @return Expected file content or error
   */
  [[nodiscard]] ParseResult<std::string>
  readFile(std::string_view filename) const;

  /**
   * @brief Split text into lines
   * @param text Text to split
   * @return Vector of lines
   */
  [[nodiscard]] std::vector<std::string>
  splitLines(std::string_view text) const;

  /**
   * @brief Parse header line to extract ID and metadata
   * @param header_line Header line
   * @return Pair of (ID, metadata)
   */
  [[nodiscard]] std::pair<std::string, std::vector<std::string>>
  parseHeader(std::string_view header_line) const;

  /**
   * @brief Convert error to string message
   * @param error Error type
   * @return Error message
   */
  [[nodiscard]] static std::string errorToString(ParseError error) noexcept;
};

} // namespace dna_motif
