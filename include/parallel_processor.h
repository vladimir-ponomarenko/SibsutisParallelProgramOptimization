#pragma once

#include "common.h"
#include "motif_finder.h"
#include "mpi_manager.h"

namespace dna_motif {

/**
 * @brief Main parallel processor coordinating MPI and OpenMP
 *
 * Orchestrates the entire motif finding process using both
 * MPI for inter-process communication and OpenMP for intra-process parallelism
 */
class ParallelProcessor {
public:
  ParallelProcessor();
  ~ParallelProcessor();

  /**
   * @brief Initialize the parallel processor
   * @param argc Command line arguments count
   * @param argv Command line arguments
   * @param num_threads Number of OpenMP threads per process
   * @return true if initialization successful
   */
  bool initialize(int argc, char *argv[], int num_threads = 0);

  /**
   * @brief Process motif finding with given input files
   * @param chip_seq_file Path to ChIP-seq sequences file
   * @param motifs_file Path to motifs file
   * @return Vector of motif results
   */
  std::vector<MotifResult> processMotifs(const std::string &chip_seq_file,
                                         const std::string &motifs_file);

  /**
   * @brief Print results to console
   * @param results Motif results to print
   */
  void printResults(const std::vector<MotifResult> &results) const;

  /**
   * @brief Save results to file
   * @param results Motif results to save
   * @param output_file Output file path
   */
  void saveResults(const std::vector<MotifResult> &results,
                   const std::string &output_file) const;

  /**
   * @brief Get performance statistics
   * @return Map with performance metrics
   */
  std::unordered_map<std::string, double> getPerformanceStats() const;

  /**
   * @brief Finalize the parallel processor
   */
  void finalize();

private:
  std::unique_ptr<MPIManager> mpi_manager_;
  std::unique_ptr<MotifFinder> motif_finder_;
  std::unique_ptr<IUPACCodes> iupac_codes_;
  std::unordered_map<std::string, double> performance_stats_;
  bool initialized_;

  /**
   * @brief Load and parse input files
   * @param chip_seq_file ChIP-seq file path
   * @param motifs_file Motifs file path
   * @return Pair of (sequences, motifs)
   */
  std::pair<std::vector<ChIPSequence>, std::vector<Motif>>
  loadInputFiles(const std::string &chip_seq_file,
                 const std::string &motifs_file);

  /**
   * @brief Process motifs in parallel using OpenMP
   * @param sequences Sequences to process
   * @param motifs Motifs to find
   * @return Local results from current process
   */
  std::vector<MotifResult>
  processMotifsParallel(const std::vector<ChIPSequence> &sequences,
                        const std::vector<Motif> &motifs);

  /**
   * @brief Update performance statistics
   * @param operation Operation name
   * @param time_seconds Time taken in seconds
   */
  void updatePerformanceStats(const std::string &operation,
                              double time_seconds);
};

} // namespace dna_motif
