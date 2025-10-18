#pragma once

#include "common.h"
#include "concepts.h"
#include <array>
#include <cctype>
#include <ranges>
#include <span>
#include <unordered_set>

namespace dna_motif {

/**
 * @brief Manages IUPAC nucleotide codes for DNA motif matching
 *
 * IUPAC codes allow representation of ambiguous nucleotides:
 * A, T, G, C - standard nucleotides
 * R = A/G, Y = T/C, S = G/C, W = A/T, K = G/T, M = A/C
 * B = C/G/T, D = A/G/T, H = A/C/T, V = A/C/G
 * N = A/T/G/C (any nucleotide)
 */
class IUPACCodes {
public:
  using nucleotide_set = std::array<char, 4>;
  using iupac_map_type = std::array<nucleotide_set, 256>;
  IUPACCodes();

  // Singleton
  static IUPACCodes &getInstance() {
    static IUPACCodes instance;
    return instance;
  }

  /**
   * @brief Check if a character is a valid IUPAC code
   * @param code Character to check
   * @return true if valid IUPAC code
   */
  [[nodiscard]] bool isValidIUPACCode(char code) const noexcept {
    char upper_code = std::toupper(code);
    return valid_codes_[static_cast<unsigned char>(upper_code)];
  }

  /**
   * @brief Get all possible nucleotides for an IUPAC code
   * @param code IUPAC code character
   * @return Vector of possible nucleotides
   */
  [[nodiscard]] std::vector<char> getNucleotides(char code) const noexcept {
    char upper_code = std::toupper(code);
    if (!isValidIUPACCode(upper_code)) {
      return {};
    }

    const auto &nucleotides =
        iupac_map_[static_cast<unsigned char>(upper_code)];
    std::vector<char> result;

    size_t count = 0;
    for (char nuc : nucleotides) {
      if (nuc != 0)
        count++;
    }

    result.reserve(count);
    for (char nuc : nucleotides) {
      if (nuc != 0)
        result.push_back(nuc);
    }

    return result;
  }

  /**
   * @brief Check if a nucleotide matches an IUPAC code
   * @param nucleotide Single nucleotide (A, T, G, C)
   * @param iupac_code IUPAC code character
   * @return true if nucleotide matches the IUPAC code
   */
  [[nodiscard]] bool matches(char nucleotide, char iupac_code) const noexcept {
    char upper_nucleotide = std::toupper(nucleotide);
    char upper_iupac_code = std::toupper(iupac_code);

    if (!isValidIUPACCode(upper_iupac_code))
      return false;

    const auto nucleotides = getNucleotides(upper_iupac_code);
    return std::ranges::find(nucleotides, upper_nucleotide) !=
           nucleotides.end();
  }

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
   * @return true if all characters are valid IUPAC codes
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
  std::array<bool, 256> valid_codes_;

  constexpr void initializeIUPACMap() noexcept;
  constexpr void addMapping(char iupac_code,
                            std::initializer_list<char> nucleotides) noexcept;
};

} // namespace dna_motif
