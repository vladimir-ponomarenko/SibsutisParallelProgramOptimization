#pragma once

#include "common.h"
#include "concepts.h"
#include <array>
#include <cctype>
#include <cstdint>
#include <ranges>
#include <span>
#include <unordered_set>

namespace dna_motif {

/**
 * @brief Manages IUPAC nucleotide codes for DNA motif matching
 *
 * IUPAC codes allow representation of ambiguous nucleotides:
 * A, T, G, C - standard nucleotides
 * Binary representation (Table 12 from task):
 * A: 0001, T: 0010, G: 0100, C: 1000
 * R = A/G (0101), Y = T/C (1010), etc.
 * N = 1111
 */
class IUPACCodes {
public:
  using nucleotide_set = std::array<char, 4>;
  using iupac_map_type = std::array<nucleotide_set, 256>;
  using bitmask_map_type = std::array<uint8_t, 256>;

  IUPACCodes();

  // Singleton
  static IUPACCodes &getInstance() {
    static IUPACCodes instance;
    return instance;
  }

  [[nodiscard]] bool isValidIUPACCode(char code) const noexcept {
    char upper_code = std::toupper(code);
    return valid_codes_[static_cast<unsigned char>(upper_code)];
  }

  [[nodiscard]] std::vector<char> getNucleotides(char code) const noexcept;

  [[nodiscard]] bool matches(char nucleotide, char iupac_code) const noexcept;

  /**
   * @brief Get 4-bit mask for a character (Table 12)
   * @param code Character (A, T, G, C, R, Y...)
   * @return 4-bit mask (e.g. A=1, T=2)
   */
  [[nodiscard]] uint8_t getBitmask(char code) const noexcept {
    return bitmask_map_[static_cast<unsigned char>(std::toupper(code))];
  }

  /**
   * @brief Create a 32-bit hash from an 8-char sequence/motif
   * Packs 8 chars * 4 bits into one uint32_t.
   * @param sequence String view of length 8
   * @return Packed hash
   */
  [[nodiscard]] uint32_t hashSequence(std::string_view sequence) const noexcept;

  // ----------------------------------------------------

  /**
   * @brief Check if a DNA sequence matches a motif pattern
   * @param sequence DNA sequence to check
   * @param motif Motif pattern with IUPAC codes
   * @param start_pos Starting position in sequence
   * @return true if sequence matches motif starting at start_pos
   */
  [[nodiscard]] bool matchesMotif(std::string_view sequence,

                                  std::string_view motif,

                                  size_t start_pos) const noexcept;

  /**
   * @brief Find all matches of a motif in a sequences
   * @param sequence DNA sequence to search
   * @param motif Motif pattern to find
   * @return Vector of starting positions where motif matches
   */
  [[nodiscard]] std::vector<size_t>
  findMotifMatches(std::string_view sequence, std::string_view motif) const;

  /**
   * @brief Get all valid IUPAC codes as a range
   * @return Range of all valid IUPAC code characters
   */
  [[nodiscard]] constexpr auto getAllCodes() const noexcept {
    return std::views::iota(0, 256) |
           std::views::filter([this](int i) { return valid_codes_[i]; }) |
           std::views::transform([](int i) { return static_cast<char>(i); });
  }

  /**
   * @brief Get count of valid IUPAC codes
   * @return Number of valid IUPAC codes
   */
  [[nodiscard]] constexpr size_t getCodeCount() const noexcept {
    return std::ranges::count(valid_codes_, true);
  }

  /**
   * @brief Check if a sequence contains only valid IUPAC codes
   * @param sequence Sequence to check
   * @return true if all characters are valid IUPAC code
   */
  [[nodiscard]] bool isValidSequence(std::string_view sequence) const noexcept {
    return std::ranges::all_of(sequence,
                               [this](char c) { return isValidIUPACCode(c); });
  }

  /**
   * @brief Get statistics about IUPAC code usage
   * @return Map with usage statistics
   */
  [[nodiscard]] std::unordered_map<char, size_t>
  getUsageStats(std::string_view sequence) const;

private:
  iupac_map_type iupac_map_;
  bitmask_map_type bitmask_map_;
  std::array<bool, 256> valid_codes_;
  constexpr void initializeIUPACMap() noexcept;
  constexpr void addMapping(char iupac_code,
                            std::initializer_list<char> nucleotides,
                            uint8_t mask) noexcept;
};

} // namespace dna_motif