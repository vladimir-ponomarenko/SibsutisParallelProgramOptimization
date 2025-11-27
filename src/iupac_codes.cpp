#include "iupac_codes.h"
#include <algorithm>
#include <cctype>
#include <format>
#include <ranges>

namespace dna_motif {

IUPACCodes::IUPACCodes() { initializeIUPACMap(); }

constexpr void IUPACCodes::initializeIUPACMap() noexcept {
  iupac_map_.fill(nucleotide_set{});
  bitmask_map_.fill(0);
  valid_codes_.fill(false);

  // Bitmask values based on Table 12:
  // A=0001 (1), T=0010 (2), G=0100 (4), C=1000 (8)

  // Standard nucleotides
  addMapping('A', {'A'}, 1); // 0001
  addMapping('T', {'T'}, 2); // 0010
  addMapping('G', {'G'}, 4); // 0100
  addMapping('C', {'C'}, 8); // 1000

  // Two-way ambiguities
  addMapping('R', {'A', 'G'}, 1 | 4); // 0101
  addMapping('Y', {'T', 'C'}, 2 | 8); // 1010
  addMapping('M', {'A', 'C'}, 1 | 8); // 1001
  addMapping('K', {'G', 'T'}, 4 | 2); // 0110
  addMapping('W', {'A', 'T'}, 1 | 2); // 0011
  addMapping('S', {'G', 'C'}, 4 | 8); // 1100

  // Three-way ambiguities
  addMapping('B', {'C', 'G', 'T'}, 8 | 4 | 2); // 1110 (Not A)
  addMapping('D', {'A', 'G', 'T'}, 1 | 4 | 2); // 0111 (Not C)
  addMapping('H', {'A', 'C', 'T'}, 1 | 8 | 2); // 1011 (Not G)
  addMapping('V', {'A', 'C', 'G'}, 1 | 8 | 4); // 1101 (Not T)

  // Four-way ambiguity
  addMapping('N', {'A', 'T', 'G', 'C'}, 1 | 2 | 4 | 8); // 1111
}

constexpr void IUPACCodes::addMapping(char iupac_code,
                                      std::initializer_list<char> nucleotides,
                                      uint8_t mask) noexcept {
  const auto index = static_cast<unsigned char>(iupac_code);
  valid_codes_[index] = true;
  bitmask_map_[index] = mask;

  auto &mapping = iupac_map_[index];
  size_t i = 0;
  for (const auto &nucleotide : nucleotides) {
    if (i < mapping.size()) {
      mapping[i++] = nucleotide;
    }
  }
}

std::vector<char> IUPACCodes::getNucleotides(char code) const noexcept {
  char upper_code = std::toupper(code);
  if (!isValidIUPACCode(upper_code)) {
    return {};
  }

  const auto &nucleotides = iupac_map_[static_cast<unsigned char>(upper_code)];
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

bool IUPACCodes::matches(char nucleotide, char iupac_code) const noexcept {
  // Bitmasks for single char comparison too
  uint8_t nuc_mask = getBitmask(nucleotide);
  uint8_t iupac_mask = getBitmask(iupac_code);
  // (Seq & Motif) == Seq
  return (nuc_mask & iupac_mask) == nuc_mask;
}

uint32_t IUPACCodes::hashSequence(std::string_view sequence) const noexcept {
  uint32_t hash = 0;
  size_t len = std::min(sequence.length(), size_t(8));

  for (size_t i = 0; i < len; ++i) {
    uint32_t mask = static_cast<uint32_t>(getBitmask(sequence[i]));
    // Shift 4 bits per position (0, 4, 8, ... 28)
    hash |= (mask << (i * 4));
  }
  return hash;
}

bool IUPACCodes::matchesMotif(std::string_view sequence, std::string_view motif,
                              size_t start_pos) const noexcept {
  if (start_pos + motif.length() > sequence.length()) {
    return false;
  }

  uint32_t motif_hash = hashSequence(motif);

  std::string_view window = sequence.substr(start_pos, motif.length());
  uint32_t seq_hash = hashSequence(window);

  // (oligo_hash & motif_hash) == oligo_hash
  return (seq_hash & motif_hash) == seq_hash;
}

std::vector<size_t> IUPACCodes::findMotifMatches(std::string_view sequence,
                                                 std::string_view motif) const {
  if (sequence.length() < motif.length()) {
    return {};
  }

  // Pre-calculate Motif Hash
  const uint32_t motif_hash = hashSequence(motif);
  const size_t motif_len = motif.length();

  std::vector<size_t> matches;
  matches.reserve(sequence.length());
  const size_t max_pos = sequence.length() - motif_len;

  for (size_t i = 0; i <= max_pos; ++i) {
    uint32_t seq_hash = 0;
    bool possible = true;

    // Inline hash calculation
    for (size_t k = 0; k < motif_len; ++k) {
      uint8_t m = bitmask_map_[static_cast<unsigned char>(sequence[i + k])];
      if (m == 0) {
        possible = false;
        break;
      }
      seq_hash |= (static_cast<uint32_t>(m) << (k * 4));
    }

    if (possible && (seq_hash & motif_hash) == seq_hash) {
      matches.push_back(i);
    }
  }

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