#pragma once

#include "common.h"
#include "concepts.h"
#include "iupac_codes.h"
#include <coroutine>
#include <execution>
#include <future>
#include <ranges>

namespace dna_motif {

// Coroutine task for async motif finding
template <typename T> struct Task {
  struct promise_type {
    Task get_return_object() {
      return Task{std::coroutine_handle<promise_type>::from_promise(*this)};
    }
    std::suspend_never initial_suspend() { return {}; }
    std::suspend_never final_suspend() noexcept { return {}; }
    void return_value(T value) { result = std::move(value); }
    void unhandled_exception() { std::terminate(); }
    T result;
  };

  std::coroutine_handle<promise_type> handle;

  Task(std::coroutine_handle<promise_type> h) : handle(h) {}
  ~Task() {
    if (handle)
      handle.destroy();
  }

  Task(const Task &) = delete;
  Task &operator=(const Task &) = delete;

  Task(Task &&other) noexcept : handle(std::exchange(other.handle, {})) {}
  Task &operator=(Task &&other) noexcept {
    if (this != &other) {
      if (handle)
        handle.destroy();
      handle = std::exchange(other.handle, {});
    }
    return *this;
  }

  bool ready() const { return handle.done(); }
  T get() { return std::move(handle.promise().result); }
};

/**
 * @brief Core motif finding algorithm
 *
 * Implements motif matching using IUPAC codes
 * with support for parallel processing, coroutines, and ranges
 */
class MotifFinder {
public:
  explicit MotifFinder(const IUPACCodes &iupac_codes);
  ~MotifFinder() = default;

  MotifFinder(const MotifFinder &) = delete;
  MotifFinder &operator=(const MotifFinder &) = delete;
  MotifFinder(MotifFinder &&) = default;
  MotifFinder &operator=(MotifFinder &&) = default;

  /**
   * @brief Find all motif matches in a set of sequences
   * @param sequences Vector of ChIP sequences to search
   * @param motifs Vector of motifs to find
   * @return Vector of motif results with match counts and frequencies
   */
  [[nodiscard]] std::vector<MotifResult>
  findMotifs(std::span<const ChIPSequence> sequences,
             std::span<const Motif> motifs);

  /**
   * @brief Find matches for a single motif in all sequences
   * @param sequences Vector of ChIP sequences
   * @param motif Motif to find
   * @return Motif result with matches
   */
  [[nodiscard]] MotifResult
  findSingleMotif(std::span<const ChIPSequence> sequences, const Motif &motif);

  /**
   * @brief Find matches for a single motif in a single sequence
   * @param sequence ChIP sequence to search
   * @param motif Motif to find
   * @param sequence_index Index of sequence in the collection
   * @return Vector of motif matches
   */
  [[nodiscard]] std::vector<MotifMatch>
  findMotifInSequence(const ChIPSequence &sequence, const Motif &motif,
                      size_t sequence_index);

  /**
   * @brief Calculate frequency of motif matches
   * @param match_count Number of sequences with matches
   * @param total_sequences Total number of sequences
   * @return Frequency as a decimal (0.0 to 1.0)
   */
  [[nodiscard]] static constexpr double
  calculateFrequency(size_t match_count, size_t total_sequences) noexcept {
    return total_sequences > 0 ? static_cast<double>(match_count) /
                                     static_cast<double>(total_sequences)
                               : 0.0;
  }

  /**
   * @brief Get performance statistics
   * @return Map with performance metrics
   */
  [[nodiscard]] const std::unordered_map<std::string, double> &
  getPerformanceStats() const noexcept {
    return performance_stats_;
  }

  /**
   * @brief Reset performance statistics
   */
  void resetPerformanceStats() noexcept { performance_stats_.clear(); }

  /**
   * @brief Find motifs asynchronously
   * @param sequences Vector of ChIP sequences
   * @param motifs Vector of motifs to find
   * @return Task containing motif results
   */
  [[nodiscard]] Task<std::vector<MotifResult>>
  findMotifsAsync(std::span<const ChIPSequence> sequences,
                  std::span<const Motif> motifs);

  /**
   * @brief Find motifs in parallel
   * @param sequences Vector of ChIP sequences
   * @param motifs Vector of motifs to find
   * @return Vector of motif results
   */
  [[nodiscard]] std::vector<MotifResult>
  findMotifsParallel(std::span<const ChIPSequence> sequences,
                     std::span<const Motif> motifs);

private:
  const IUPACCodes &iupac_codes_;
  std::unordered_map<std::string, double> performance_stats_;

  /**
   * @brief Check if a sequence segment matches a motif
   * @param sequence DNA sequence
   * @param motif Motif pattern
   * @param start_pos Starting position in sequence
   * @return true if segment matches motif
   */
  [[nodiscard]] bool matchesAtPosition(std::string_view sequence,
                                       std::string_view motif,
                                       size_t start_pos) const noexcept;

  /**
   * @brief Update performance statistics
   * @param operation Operation name
   * @param time_seconds Time taken in seconds
   */
  void updatePerformanceStats(std::string_view operation,
                              double time_seconds) noexcept;

  /**
   * @brief Process a single motif with timing
   * @param sequences Sequences to search
   * @param motif Motif to find
   * @return Motif result
   */
  [[nodiscard]] MotifResult
  processSingleMotif(std::span<const ChIPSequence> sequences,
                     const Motif &motif);
};

} // namespace dna_motif
