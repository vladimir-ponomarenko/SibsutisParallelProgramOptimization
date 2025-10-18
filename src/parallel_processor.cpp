#include "parallel_processor.h"
#include "dna_parser.h"
#include <fstream>
#include <iomanip>
#include <set>

namespace dna_motif {

ParallelProcessor::ParallelProcessor() : initialized_(false) {}

ParallelProcessor::~ParallelProcessor() { finalize(); }

bool ParallelProcessor::initialize(int argc, char *argv[], int num_threads) {
  // Initialize MPI
  mpi_manager_ = std::make_unique<MPIManager>();
  if (!mpi_manager_->initialize(argc, argv)) {
    std::cerr << "Failed to initialize MPI" << std::endl;
    return false;
  }

  if (num_threads > 0) {
    omp_set_num_threads(num_threads);
  }

  // Initialize IUPAC codes
  iupac_codes_ = std::make_unique<IUPACCodes>();

  // Initialize motif finder
  motif_finder_ = std::make_unique<MotifFinder>(*iupac_codes_);

  initialized_ = true;

  if (mpi_manager_->isMaster()) {
    std::cout << "ParallelProcessor initialized with "
              << mpi_manager_->getSize() << " MPI processes and "
              << omp_get_max_threads() << " OpenMP threads per process"
              << std::endl;
  }

  return true;
}

std::vector<MotifResult>
ParallelProcessor::processMotifs(const std::string &chip_seq_file,
                                 const std::string &motifs_file) {
  if (!initialized_) {
    throw std::runtime_error("ParallelProcessor not initialized");
  }

  Timer total_timer;

  auto [sequences, motifs] = loadInputFiles(chip_seq_file, motifs_file);

  if (mpi_manager_->isMaster()) {
    std::cout << "Loaded " << sequences.size() << " sequences and "
              << motifs.size() << " motifs" << std::endl;
  }

  // Distribute work among MPI processes
  std::vector<ChIPSequence> local_sequences =
      mpi_manager_->distributeSequences(sequences);
  std::vector<Motif> local_motifs = mpi_manager_->broadcastMotifs(motifs);

  if (mpi_manager_->isMaster()) {
    std::cout << "Work distributed. Processing motifs..." << std::endl;
  }

  // Process motifs in parallel using OpenMP
  std::vector<MotifResult> local_results =
      processMotifsParallel(local_sequences, local_motifs);

  // Gather results from all processes
  std::vector<MotifResult> all_results =
      mpi_manager_->gatherResults(local_results);

  double total_time = total_timer.elapsed();
  updatePerformanceStats("total_processing_time", total_time);

  if (mpi_manager_->isMaster()) {
    std::cout << "Processing completed in " << std::fixed
              << std::setprecision(2) << total_time << " seconds" << std::endl;
  }

  return all_results;
}

void ParallelProcessor::printResults(
    const std::vector<MotifResult> &results) const {
  if (!mpi_manager_->isMaster()) {
    return;
  }

  std::cout << "\n=== MOTIF FINDING RESULTS ===" << std::endl;
  std::cout << std::setw(20) << "Motif Pattern" << std::setw(15)
            << "Match Count" << std::setw(15) << "Frequency" << std::endl;
  std::cout << std::string(50, '-') << std::endl;

  for (const auto &result : results) {
    std::cout << std::setw(20) << result.motif_pattern << std::setw(15)
              << result.match_count << std::setw(15) << std::fixed
              << std::setprecision(4) << result.frequency << std::endl;
  }

  std::cout << std::endl;
}

void ParallelProcessor::saveResults(const std::vector<MotifResult> &results,
                                    const std::string &output_file) const {
  if (!mpi_manager_->isMaster()) {
    return;
  }

  std::ofstream file(output_file);
  if (!file.is_open()) {
    std::cerr << "Cannot open output file: " << output_file << std::endl;
    return;
  }

  file << "Motif_Pattern\tMatch_Count\tFrequency\n";

  for (const auto &result : results) {
    file << result.motif_pattern << "\t" << result.match_count << "\t"
         << std::fixed << std::setprecision(6) << result.frequency << "\n";
  }

  file.close();

  std::cout << "Results saved to: " << output_file << std::endl;
}

std::unordered_map<std::string, double>
ParallelProcessor::getPerformanceStats() const {
  return performance_stats_;
}

void ParallelProcessor::finalize() {
  if (initialized_ && mpi_manager_) {
    mpi_manager_->finalize();
    initialized_ = false;
  }
}

std::pair<std::vector<ChIPSequence>, std::vector<Motif>>
ParallelProcessor::loadInputFiles(const std::string &chip_seq_file,
                                  const std::string &motifs_file) {
  Timer timer;

  DNAParser parser;
  std::vector<ChIPSequence> sequences;
  std::vector<Motif> motifs;

  try {
    auto sequences_result = parser.parseChIPSequences(chip_seq_file);
    auto motifs_result = parser.parseMotifs(motifs_file);

    if (!sequences_result) {
      throw std::runtime_error(
          "Failed to parse ChIP sequences: " +
          std::to_string(static_cast<int>(sequences_result.error())));
    }

    if (!motifs_result) {
      throw std::runtime_error(
          "Failed to parse motifs: " +
          std::to_string(static_cast<int>(motifs_result.error())));
    }

    sequences = std::move(sequences_result.value());
    motifs = std::move(motifs_result.value());

    auto stats = parser.getStatistics();
    if (mpi_manager_->isMaster()) {
      std::cout << "Parsing statistics:" << std::endl;
      for (const auto &stat : stats) {
        std::cout << "  " << stat.first << ": " << stat.second << std::endl;
      }
    }

  } catch (const std::exception &e) {
    std::cerr << "Error loading input files: " << e.what() << std::endl;
    throw;
  }

  double load_time = timer.elapsed();
  updatePerformanceStats("file_loading_time", load_time);

  return {sequences, motifs};
}

std::vector<MotifResult> ParallelProcessor::processMotifsParallel(
    const std::vector<ChIPSequence> &sequences,
    const std::vector<Motif> &motifs) {
  Timer timer;
  std::vector<MotifResult> results;
  results.reserve(motifs.size());

// Process motifs in parallel using OpenMP
#pragma omp parallel
  {
    std::vector<MotifResult> local_results;
    local_results.reserve(motifs.size());

#pragma omp for schedule(dynamic)
    for (size_t i = 0; i < motifs.size(); ++i) {
      MotifResult result = motif_finder_->findSingleMotif(sequences, motifs[i]);
      local_results.push_back(result);
    }

// Merge results from all threads
#pragma omp critical
    {
      results.insert(results.end(), local_results.begin(), local_results.end());
    }
  }

  double parallel_time = timer.elapsed();
  updatePerformanceStats("parallel_processing_time", parallel_time);

  return results;
}

void ParallelProcessor::updatePerformanceStats(const std::string &operation,
                                               double time_seconds) {
  performance_stats_[operation] = time_seconds;
}

} // namespace dna_motif
