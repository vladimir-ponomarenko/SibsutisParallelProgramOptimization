#pragma once

#include "common.h"

namespace dna_motif {

/**
 * @brief MPI manager for distributed motif finding
 *
 * Handles MPI communication and data distribution
 * for parallel motif finding across multiple processes
 */
class MPIManager {
public:
  MPIManager();
  ~MPIManager();

  /**
   * @brief Initialize MPI environment
   * @param argc Command line arguments count
   * @param argv Command line arguments
   * @return true if initialization successful
   */
  bool initialize(int argc, char *argv[]);

  /**
   * @brief Finalize MPI environment
   */
  void finalize();

  /**
   * @brief Get current process rank
   * @return Process rank (0 to size-1)
   */
  int getRank() const { return rank_; }

  /**
   * @brief Get total number of processes
   * @return Total number of MPI processes
   */
  int getSize() const { return size_; }

  /**
   * @brief Check if current process is master
   * @return true if rank is 0
   */
  bool isMaster() const { return rank_ == 0; }

  /**
   * @brief Distribute sequences among processes
   * @param all_sequences All sequences to distribute
   * @return Sequences assigned to current process
   */
  std::vector<ChIPSequence>
  distributeSequences(const std::vector<ChIPSequence> &all_sequences);

  /**
   * @brief Broadcast motifs to all processes
   * @param motifs Motifs to broadcast
   * @return Motifs received by current process
   */
  std::vector<Motif> broadcastMotifs(const std::vector<Motif> &motifs);

  /**
   * @brief Gather motif results from all processes
   * @param local_results Results from current process
   * @return Combined results from all processes
   */
  std::vector<MotifResult>
  gatherResults(const std::vector<MotifResult> &local_results);

  /**
   * @brief Synchronize all processes
   */
  void synchronize();

  /**
   * @brief Get communication statistics
   * @return Map with communication metrics
   */
  std::unordered_map<std::string, double> getCommunicationStats() const;

  /**
   * @brief Calculate work distribution for sequences
   * @param total_sequences Total number of sequences
   * @param process_rank Current process rank
   * @param total_processes Total number of processes
   * @return Pair of (start_index, count) for this process
   */
  std::pair<size_t, size_t>
  calculateWorkDistribution(size_t total_sequences, int process_rank,
                            int total_processes) const;

private:
  int rank_;
  int size_;
  bool initialized_;
  std::unordered_map<std::string, double> comm_stats_;

  /**
   * @brief Update communication statistics
   * @param operation Communication operation name
   * @param bytes Number of bytes transferred
   * @param time_seconds Time taken in seconds
   */
  void updateCommStats(const std::string &operation, size_t bytes,
                       double time_seconds);
};

} // namespace dna_motif
