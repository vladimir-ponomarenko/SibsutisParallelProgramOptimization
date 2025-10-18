#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <barrier>
#include <bit>
#include <cassert>
#include <chrono>
#include <concepts>
#include <condition_variable>
#include <coroutine>
#include <execution>
#include <expected>
#include <format>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <latch>
#include <memory>
#include <mpi.h>
#include <mutex>
#include <numbers>
#include <omp.h>
#include <optional>
#include <ranges>
#include <semaphore>
#include <shared_mutex>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <version>

namespace dna_motif {

template <typename T>
concept DNASequence = requires(T t) {
  { t.sequence } -> std::convertible_to<std::string>;
  { t.id } -> std::convertible_to<std::string>;
  requires std::ranges::range<decltype(t.metadata)>;
};

template <typename T>
concept MotifConcept = requires(T t) {
  { t.pattern } -> std::convertible_to<std::string>;
  { t.score1 } -> std::convertible_to<double>;
  { t.score2 } -> std::convertible_to<double>;
  { t.score3 } -> std::convertible_to<double>;
};

inline constexpr size_t CHIP_SEQ_LENGTH = 40;
inline constexpr size_t MOTIF_LENGTH = 8;
inline constexpr size_t IUPAC_CODE_SIZE = 15;
inline constexpr std::string_view VALID_DNA_NUCLEOTIDES = "ATGC";
inline constexpr std::string_view IUPAC_CODES = "ATGCWSRYMKBDHVN";

struct ChIPSequence {
  std::string id;
  std::string sequence;
  std::vector<std::string> metadata;

  ChIPSequence() = default;

  ChIPSequence(std::string_view seq_id, std::string_view seq)
      : id(seq_id), sequence(seq) {}

  ChIPSequence(ChIPSequence &&) = default;
  ChIPSequence &operator=(ChIPSequence &&) = default;

  ChIPSequence(const ChIPSequence &) = default;
  ChIPSequence &operator=(const ChIPSequence &) = default;

  ~ChIPSequence() = default;

  bool operator==(const ChIPSequence &other) const noexcept {
    return id == other.id && sequence == other.sequence &&
           metadata == other.metadata;
  }

  auto operator<=>(const ChIPSequence &other) const noexcept = default;

  std::span<const char> getSequenceSpan() const noexcept {
    return std::span<const char>(sequence.data(), sequence.size());
  }

  bool isValid() const noexcept {
    return !id.empty() && !sequence.empty() &&
           sequence.size() == CHIP_SEQ_LENGTH;
  }
};

struct Motif {
  std::string pattern;
  double score1;
  double score2;
  double score3;

  Motif() = default;

  Motif(std::string_view pat, double s1, double s2, double s3)
      : pattern(pat), score1(s1), score2(s2), score3(s3) {}

  Motif(Motif &&) = default;
  Motif &operator=(Motif &&) = default;

  Motif(const Motif &) = default;
  Motif &operator=(const Motif &) = default;

  ~Motif() = default;

  bool operator==(const Motif &other) const noexcept {
    return pattern == other.pattern &&
           std::abs(score1 - other.score1) <
               std::numeric_limits<double>::epsilon() &&
           std::abs(score2 - other.score2) <
               std::numeric_limits<double>::epsilon() &&
           std::abs(score3 - other.score3) <
               std::numeric_limits<double>::epsilon();
  }

  std::partial_ordering operator<=>(const Motif &other) const noexcept {
    if (auto cmp = pattern <=> other.pattern; cmp != 0)
      return cmp;
    if (auto cmp = score1 <=> other.score1; cmp != 0)
      return cmp;
    if (auto cmp = score2 <=> other.score2; cmp != 0)
      return cmp;
    return score3 <=> other.score3;
  }

  std::span<const char> getPatternSpan() const noexcept {
    return std::span<const char>(pattern.data(), pattern.size());
  }

  bool isValid() const noexcept {
    return !pattern.empty() && pattern.size() == MOTIF_LENGTH;
  }
};

struct MotifMatch {
  size_t sequence_index;
  size_t position;
  std::string matched_sequence;

  MotifMatch(size_t seq_idx, size_t pos, std::string_view matched)
      : sequence_index(seq_idx), position(pos), matched_sequence(matched) {}

  MotifMatch() = default;

  MotifMatch(MotifMatch &&) = default;
  MotifMatch &operator=(MotifMatch &&) = default;

  MotifMatch(const MotifMatch &) = default;
  MotifMatch &operator=(const MotifMatch &) = default;

  ~MotifMatch() = default;

  bool operator==(const MotifMatch &other) const noexcept {
    return sequence_index == other.sequence_index &&
           position == other.position &&
           matched_sequence == other.matched_sequence;
  }

  auto operator<=>(const MotifMatch &other) const noexcept {
    if (auto cmp = sequence_index <=> other.sequence_index; cmp != 0)
      return cmp;
    if (auto cmp = position <=> other.position; cmp != 0)
      return cmp;
    return matched_sequence <=> other.matched_sequence;
  }

  std::span<const char> getMatchedSequenceSpan() const noexcept {
    return std::span<const char>(matched_sequence.data(),
                                 matched_sequence.size());
  }
};

struct MotifResult {
  std::string motif_pattern;
  size_t match_count;
  double frequency;
  std::vector<MotifMatch> matches;

  MotifResult() : match_count(0), frequency(0.0) {}

  MotifResult(std::string_view pattern)
      : motif_pattern(pattern), match_count(0), frequency(0.0) {}

  MotifResult(MotifResult &&) = default;
  MotifResult &operator=(MotifResult &&) = default;

  MotifResult(const MotifResult &) = default;
  MotifResult &operator=(const MotifResult &) = default;

  ~MotifResult() = default;

  bool operator==(const MotifResult &other) const noexcept {
    return motif_pattern == other.motif_pattern &&
           match_count == other.match_count &&
           std::abs(frequency - other.frequency) <
               std::numeric_limits<double>::epsilon() &&
           matches == other.matches;
  }

  std::partial_ordering operator<=>(const MotifResult &other) const noexcept {
    if (auto cmp = motif_pattern <=> other.motif_pattern; cmp != 0)
      return cmp;
    if (auto cmp = match_count <=> other.match_count; cmp != 0)
      return cmp;
    if (auto cmp = frequency <=> other.frequency; cmp != 0)
      return cmp;
    return matches <=> other.matches;
  }

  std::span<const MotifMatch> getMatchesSpan() const noexcept {
    return std::span<const MotifMatch>(matches.data(), matches.size());
  }

  bool isValid() const noexcept {
    return !motif_pattern.empty() && motif_pattern.size() == MOTIF_LENGTH;
  }

  void calculateFrequency(size_t total_sequences) noexcept {
    if (total_sequences > 0) {
      frequency = static_cast<double>(match_count) /
                  static_cast<double>(total_sequences);
    } else {
      frequency = 0.0;
    }
  }
};

[[nodiscard]] std::string trim(std::string_view str);
[[nodiscard]] std::vector<std::string> split(std::string_view str,
                                             char delimiter);
[[nodiscard]] bool isValidDNASequence(std::string_view sequence);
[[nodiscard]] bool isValidIUPACCode(char code);
[[nodiscard]] std::string toUpperCase(std::string_view str);
[[nodiscard]] std::string toLowerCase(std::string_view str);
[[nodiscard]] bool startsWith(std::string_view str, std::string_view prefix);
[[nodiscard]] bool endsWith(std::string_view str, std::string_view suffix);
[[nodiscard]] std::string
replaceAll(std::string_view str, std::string_view from, std::string_view to);
[[nodiscard]] std::vector<std::string> splitLines(std::string_view str);
[[nodiscard]] std::string join(const std::vector<std::string> &strings,
                               std::string_view delimiter);
[[nodiscard]] std::string formatProgress(size_t current, size_t total,
                                         std::string_view operation);
void printProgress(size_t current, size_t total, std::string_view operation);

class Timer {
public:
  using clock_type = std::chrono::high_resolution_clock;
  using time_point = clock_type::time_point;
  using duration = clock_type::duration;

  Timer() : start_(clock_type::now()) {}

  void reset() noexcept { start_ = clock_type::now(); }

  [[nodiscard]] double elapsed() const noexcept {
    const auto end = clock_type::now();
    const auto duration_ms =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start_);
    return static_cast<double>(duration_ms.count()) / 1000.0;
  }

  [[nodiscard]] double elapsedMicroseconds() const noexcept {
    const auto end = clock_type::now();
    const auto duration_us =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start_);
    return static_cast<double>(duration_us.count());
  }

  [[nodiscard]] double elapsedNanoseconds() const noexcept {
    const auto end = clock_type::now();
    const auto duration_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_);
    return static_cast<double>(duration_ns.count());
  }

  [[nodiscard]] bool isRunning() const noexcept {
    return start_ != time_point{};
  }

private:
  time_point start_;
};

class ScopedTimer {
public:
  explicit ScopedTimer(std::string_view operation)
      : operation_(operation), timer_() {}

  ~ScopedTimer() {
    if (!operation_.empty()) {
      std::cout << std::format("{} completed in {:.3f} seconds\n", operation_,
                               timer_.elapsed());
    }
  }

  [[nodiscard]] double elapsed() const noexcept { return timer_.elapsed(); }

private:
  std::string operation_;
  Timer timer_;
};

class PerformanceCounter {
public:
  PerformanceCounter() = default;

  void increment() noexcept { count_.fetch_add(1, std::memory_order_relaxed); }

  void add(size_t value) noexcept {
    count_.fetch_add(value, std::memory_order_relaxed);
  }

  [[nodiscard]] size_t get() const noexcept {
    return count_.load(std::memory_order_relaxed);
  }

  void reset() noexcept { count_.store(0, std::memory_order_relaxed); }

private:
  std::atomic<size_t> count_{0};
};

} // namespace dna_motif
