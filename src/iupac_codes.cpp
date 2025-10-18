#include "iupac_codes.h"
#include <algorithm>
#include <cctype>
#include <format>
#include <ranges>

namespace dna_motif {

IUPACCodes::IUPACCodes() { initializeIUPACMap(); }

constexpr void IUPACCodes::initializeIUPACMap() noexcept {
  iupac_map_.fill(nucleotide_set{});
  valid_codes_.fill(false);

  // Standard nucleotides
  addMapping('A', {'A'});
  addMapping('T', {'T'});
  addMapping('G', {'G'});
  addMapping('C', {'C'});

  // Two-way ambiguities
  addMapping('R', {'A', 'G'}); // puRine
  addMapping('Y', {'T', 'C'}); // pYrimidine
  addMapping('S', {'G', 'C'}); // Strong (3 H-bonds)
  addMapping('W', {'A', 'T'}); // Weak (2 H-bonds)
  addMapping('K', {'G', 'T'}); // Keto
  addMapping('M', {'A', 'C'}); // aMino

  // Three-way ambiguities
  addMapping('B', {'C', 'G', 'T'}); // not A
  addMapping('D', {'A', 'G', 'T'}); // not C
  addMapping('H', {'A', 'C', 'T'}); // not G
  addMapping('V', {'A', 'C', 'G'}); // not T

  // Four-way ambiguity
  addMapping('N', {'A', 'T', 'G', 'C'}); // aNy nucleotide
}

constexpr void
IUPACCodes::addMapping(char iupac_code,
                       std::initializer_list<char> nucleotides) noexcept {
  const auto index = static_cast<unsigned char>(iupac_code);
  valid_codes_[index] = true;

  auto &mapping = iupac_map_[index];
  size_t i = 0;
  for (const auto &nucleotide : nucleotides) {
    if (i < mapping.size()) {
      mapping[i++] = nucleotide;
    }
  }
}

bool IUPACCodes::matchesMotif(std::string_view sequence, std::string_view motif,
                              size_t start_pos) const noexcept {
  if (start_pos + motif.length() > sequence.length()) {
    return false;
  }

  return std::ranges::all_of(
      std::views::iota(0uz, motif.length()),
      [&](size_t i) { return matches(sequence[start_pos + i], motif[i]); });
}

std::vector<size_t> IUPACCodes::findMotifMatches(std::string_view sequence,
                                                 std::string_view motif) const {
  if (sequence.length() < motif.length()) {
    return {};
  }

  const size_t max_pos = sequence.length() - motif.length();

  std::vector<size_t> matches;
  auto match_positions =
      std::views::iota(0uz, max_pos + 1) | std::views::filter([&](size_t pos) {
        return matchesMotif(sequence, motif, pos);
      });

  std::ranges::copy(match_positions, std::back_inserter(matches));
  return matches;
}

std::unordered_map<char, size_t>
IUPACCodes::getUsageStats(std::string_view sequence) const {
  std::unordered_map<char, size_t> stats;

  for (const auto &code : sequence) {
    if (isValidIUPACCode(code)) {
      stats[code]++;
    }
  }

  return stats;
}

} // namespace dna_motif
