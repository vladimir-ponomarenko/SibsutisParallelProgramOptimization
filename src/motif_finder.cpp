#include "motif_finder.h"
#include <algorithm>
#include <execution>
#include <format>
#include <iomanip>
#include <ranges>
#include <set>

namespace dna_motif {

MotifFinder::MotifFinder(const IUPACCodes &iupac_codes)
    : iupac_codes_(iupac_codes) {}

std::vector<MotifResult>
MotifFinder::findMotifs(std::span<const ChIPSequence> sequences,
                        std::span<const Motif> motifs) {
  Timer timer;
  std::vector<MotifResult> results;
  results.reserve(motifs.size());

  auto motif_results = motifs | std::views::transform([&](const Motif &motif) {
                         return processSingleMotif(sequences, motif);
                       });

  std::ranges::copy(motif_results, std::back_inserter(results));

  double total_time = timer.elapsed();
  updatePerformanceStats("find_motifs_total", total_time);

  return results;
}

MotifResult
MotifFinder::findSingleMotif(std::span<const ChIPSequence> sequences,
                             const Motif &motif) {
  Timer timer;
  MotifResult result(motif.pattern);

  auto sequence_indices = std::views::iota(0uz, sequences.size());

  for (const auto &[i, sequence] :
       std::views::zip(sequence_indices, sequences)) {
    auto matches = findMotifInSequence(sequence, motif, i);

    if (!matches.empty()) {
      result.match_count++;
      result.matches.push_back(std::move(matches[0]));
    }
  }

  result.calculateFrequency(sequences.size());

  double motif_time = timer.elapsed();
  updatePerformanceStats("find_single_motif", motif_time);

  return result;
}

std::vector<MotifMatch>
MotifFinder::findMotifInSequence(const ChIPSequence &sequence,
                                 const Motif &motif, size_t sequence_index) {
  if (sequence.sequence.length() < motif.pattern.length()) {
    return {};
  }

  auto match_positions =
      iupac_codes_.findMotifMatches(sequence.sequence, motif.pattern);

  std::vector<MotifMatch> matches;
  auto match_objects =
      match_positions | std::views::transform([&](size_t pos) {
        std::string matched_sequence =
            sequence.sequence.substr(pos, motif.pattern.length());
        return MotifMatch(sequence_index, pos, matched_sequence);
      });

  std::ranges::copy(match_objects, std::back_inserter(matches));
  return matches;
}

bool MotifFinder::matchesAtPosition(std::string_view sequence,
                                    std::string_view motif,
                                    size_t start_pos) const noexcept {
  return iupac_codes_.matchesMotif(sequence, motif, start_pos);
}

void MotifFinder::updatePerformanceStats(std::string_view operation,
                                         double time_seconds) noexcept {
  performance_stats_[std::string(operation)] = time_seconds;
}

MotifResult
MotifFinder::processSingleMotif(std::span<const ChIPSequence> sequences,
                                const Motif &motif) {
  Timer timer;
  MotifResult result(motif.pattern);

  auto sequence_indices = std::views::iota(0uz, sequences.size());

  for (const auto &[i, sequence] :
       std::views::zip(sequence_indices, sequences)) {
    auto matches = findMotifInSequence(sequence, motif, i);

    if (!matches.empty()) {
      result.match_count++;
      result.matches.push_back(std::move(matches[0]));
    }
  }

  result.calculateFrequency(sequences.size());

  double motif_time = timer.elapsed();
  updatePerformanceStats("process_single_motif", motif_time);

  return result;
}

Task<std::vector<MotifResult>>
MotifFinder::findMotifsAsync(std::span<const ChIPSequence> sequences,
                             std::span<const Motif> motifs) {
  co_return findMotifs(sequences, motifs);
}

std::vector<MotifResult>
MotifFinder::findMotifsParallel(std::span<const ChIPSequence> sequences,
                                std::span<const Motif> motifs) {
  Timer timer;
  std::vector<MotifResult> results;
  results.reserve(motifs.size());

  auto motif_results = motifs | std::views::transform([&](const Motif &motif) {
                         return processSingleMotif(sequences, motif);
                       });

  std::ranges::copy(motif_results, std::back_inserter(results));

  double total_time = timer.elapsed();
  updatePerformanceStats("find_motifs_parallel", total_time);

  return results;
}

} // namespace dna_motif
