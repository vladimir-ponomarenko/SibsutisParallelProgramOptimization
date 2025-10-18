#pragma once

#include <concepts>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

namespace dna_motif::concepts {

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

template <typename T>
concept MotifMatch = requires(T t) {
  { t.sequence_index } -> std::convertible_to<size_t>;
  { t.position } -> std::convertible_to<size_t>;
  { t.matched_sequence } -> std::convertible_to<std::string>;
};

template <typename T>
concept MotifResult = requires(T t) {
  { t.motif_pattern } -> std::convertible_to<std::string>;
  { t.match_count } -> std::convertible_to<size_t>;
  { t.frequency } -> std::convertible_to<double>;
  requires std::ranges::range<decltype(t.matches)>;
};

template <typename T>
concept Parser = requires(T t, const std::string &filename) {
  { t.parseChIPSequences(filename) } -> std::ranges::range;
  { t.parseMotifs(filename) } -> std::ranges::range;
  { t.validateSequence(filename) } -> std::convertible_to<bool>;
};

template <typename T>
concept MotifFinderConcept = requires(
    T t, const std::vector<ChIPSequence> &sequences, const Motif &motif) {
  { t.findSingleMotif(sequences, motif) } -> MotifResult;
  {
    t.findMotifInSequence(sequences[0], motif, size_t{})
  } -> std::ranges::range;
};

template <typename T>
concept MPIManagerConcept =
    requires(T t, const std::vector<ChIPSequence> &data) {
      { t.initialize(0, nullptr) } -> std::convertible_to<bool>;
      { t.getRank() } -> std::convertible_to<int>;
      { t.getSize() } -> std::convertible_to<int>;
      { t.isMaster() } -> std::convertible_to<bool>;
      { t.distributeSequences(data) } -> std::ranges::range;
      { t.broadcastMotifs(data) } -> std::ranges::range;
      { t.gatherResults(data) } -> std::ranges::range;
    };

template <typename T>
concept ParallelProcessor =
    requires(T t, const std::string &file1, const std::string &file2) {
      { t.initialize(0, nullptr, 0) } -> std::convertible_to<bool>;
      { t.processMotifs(file1, file2) } -> std::ranges::range;
      { t.finalize() } -> std::same_as<void>;
    };

template <typename T>
concept IUPACCodes = requires(T t, char code, char nucleotide) {
  { t.isValidIUPACCode(code) } -> std::convertible_to<bool>;
  { t.getNucleotides(code) } -> std::ranges::range;
  { t.matches(nucleotide, code) } -> std::convertible_to<bool>;
  {
    t.matchesMotif(std::string{}, std::string{}, size_t{})
  } -> std::convertible_to<bool>;
};

template <typename T>
concept ParallelProcessable = std::ranges::range<T> && requires(T t) {
  { std::ranges::size(t) } -> std::convertible_to<size_t>;
  { std::ranges::begin(t) } -> std::input_iterator;
  { std::ranges::end(t) } -> std::sentinel_for<decltype(std::ranges::begin(t))>;
};

template <typename T>
concept PerformanceStats = requires(T t, const std::string &key, double value) {
  { t[key] = value } -> std::same_as<double &>;
  { t.at(key) } -> std::convertible_to<double>;
  { t.find(key) } -> std::convertible_to<bool>;
};

template <typename T>
concept FileOperation = requires(T t, const std::string &filename) {
  { t.open(filename) } -> std::convertible_to<bool>;
  { t.is_open() } -> std::convertible_to<bool>;
  { t.close() } -> std::same_as<void>;
};

template <typename T>
concept Timer = requires(T t) {
  { t.reset() } -> std::same_as<void>;
  { t.elapsed() } -> std::convertible_to<double>;
};

} // namespace dna_motif::concepts
