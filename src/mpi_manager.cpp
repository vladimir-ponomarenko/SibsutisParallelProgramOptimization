#include "mpi_manager.h"
#include <stdexcept>

namespace dna_motif {

MPIManager::MPIManager() : rank_(0), size_(1), initialized_(false) {}

MPIManager::~MPIManager() {
  if (initialized_) {
    finalize();
  }
}

bool MPIManager::initialize(int argc, char *argv[]) {
  int mpi_initialized;
  MPI_Initialized(&mpi_initialized);

  if (!mpi_initialized) {
    int result = MPI_Init(&argc, &argv);
    if (result != MPI_SUCCESS) {
      return false;
    }
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
  initialized_ = true;

  return true;
}

void MPIManager::finalize() {
  if (initialized_) {
    MPI_Finalize();
    initialized_ = false;
  }
}

std::vector<ChIPSequence> MPIManager::distributeSequences(
    const std::vector<ChIPSequence> &all_sequences) {
  Timer timer;
  std::vector<ChIPSequence> local_sequences;

  if (isMaster()) {
    // Master process distributes sequences
    size_t total_sequences = all_sequences.size();

    // Send total count to all processes
    for (int dest = 1; dest < size_; ++dest) {
      MPI_Send(&total_sequences, 1, MPI_UNSIGNED_LONG, dest, 0, MPI_COMM_WORLD);
    }

    // Distribute sequences
    for (int dest = 1; dest < size_; ++dest) {
      auto work_dist = calculateWorkDistribution(total_sequences, dest, size_);
      size_t start_idx = work_dist.first;
      size_t count = work_dist.second;

      // Send count
      MPI_Send(&count, 1, MPI_UNSIGNED_LONG, dest, 1, MPI_COMM_WORLD);

      // Send sequences
      for (size_t i = 0; i < count; ++i) {
        const auto &seq = all_sequences[start_idx + i];

        // Send ID length and ID
        int id_len = seq.id.length();
        MPI_Send(&id_len, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
        MPI_Send(seq.id.c_str(), id_len, MPI_CHAR, dest, 3, MPI_COMM_WORLD);

        // Send sequence length and sequence
        int seq_len = seq.sequence.length();
        MPI_Send(&seq_len, 1, MPI_INT, dest, 4, MPI_COMM_WORLD);
        MPI_Send(seq.sequence.c_str(), seq_len, MPI_CHAR, dest, 5,
                 MPI_COMM_WORLD);

        // Send metadata count and metadata
        int meta_count = seq.metadata.size();
        MPI_Send(&meta_count, 1, MPI_INT, dest, 6, MPI_COMM_WORLD);

        for (const auto &meta : seq.metadata) {
          int meta_len = meta.length();
          MPI_Send(&meta_len, 1, MPI_INT, dest, 7, MPI_COMM_WORLD);
          MPI_Send(meta.c_str(), meta_len, MPI_CHAR, dest, 8, MPI_COMM_WORLD);
        }
      }
    }

    // Master keeps its own portion
    auto work_dist = calculateWorkDistribution(total_sequences, 0, size_);
    size_t start_idx = work_dist.first;
    size_t count = work_dist.second;

    local_sequences.assign(all_sequences.begin() + start_idx,
                           all_sequences.begin() + start_idx + count);
  } else {
    // Worker processes receive sequences
    size_t total_sequences;
    MPI_Recv(&total_sequences, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    size_t count;
    MPI_Recv(&count, 1, MPI_UNSIGNED_LONG, 0, 1, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    local_sequences.reserve(count);

    for (size_t i = 0; i < count; ++i) {
      ChIPSequence seq;

      // Receive ID
      int id_len;
      MPI_Recv(&id_len, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      seq.id.resize(id_len);
      MPI_Recv(&seq.id[0], id_len, MPI_CHAR, 0, 3, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      // Receive sequence
      int seq_len;
      MPI_Recv(&seq_len, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      seq.sequence.resize(seq_len);
      MPI_Recv(&seq.sequence[0], seq_len, MPI_CHAR, 0, 5, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      // Receive metadata
      int meta_count;
      MPI_Recv(&meta_count, 1, MPI_INT, 0, 6, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      seq.metadata.reserve(meta_count);

      for (int j = 0; j < meta_count; ++j) {
        int meta_len;
        MPI_Recv(&meta_len, 1, MPI_INT, 0, 7, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        std::string meta(meta_len, '\0');
        MPI_Recv(&meta[0], meta_len, MPI_CHAR, 0, 8, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        seq.metadata.push_back(meta);
      }

      local_sequences.push_back(seq);
    }
  }

  double comm_time = timer.elapsed();
  updateCommStats("distribute_sequences",
                  all_sequences.size() * sizeof(ChIPSequence), comm_time);

  return local_sequences;
}

std::vector<Motif>
MPIManager::broadcastMotifs(const std::vector<Motif> &motifs) {
  Timer timer;
  std::vector<Motif> received_motifs;

  if (isMaster()) {
    // Master broadcasts motifs
    int motif_count = motifs.size();
    MPI_Bcast(&motif_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    for (const auto &motif : motifs) {
      // Broadcast pattern
      int pattern_len = motif.pattern.length();
      MPI_Bcast(&pattern_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(const_cast<char *>(motif.pattern.c_str()), pattern_len,
                MPI_CHAR, 0, MPI_COMM_WORLD);

      // Broadcast scores
      double scores[3] = {motif.score1, motif.score2, motif.score3};
      MPI_Bcast(scores, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    received_motifs = motifs;
  } else {
    // Workers receive motifs
    int motif_count;
    MPI_Bcast(&motif_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    received_motifs.reserve(motif_count);

    for (int i = 0; i < motif_count; ++i) {
      Motif motif;

      // Receive pattern
      int pattern_len;
      MPI_Bcast(&pattern_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
      motif.pattern.resize(pattern_len);
      MPI_Bcast(&motif.pattern[0], pattern_len, MPI_CHAR, 0, MPI_COMM_WORLD);

      // Receive scores
      double scores[3];
      MPI_Bcast(scores, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      motif.score1 = scores[0];
      motif.score2 = scores[1];
      motif.score3 = scores[2];

      received_motifs.push_back(motif);
    }
  }

  double comm_time = timer.elapsed();
  updateCommStats("broadcast_motifs", motifs.size() * sizeof(Motif), comm_time);

  return received_motifs;
}

std::vector<MotifResult>
MPIManager::gatherResults(const std::vector<MotifResult> &local_results) {
  Timer timer;
  std::vector<MotifResult> all_results;

  if (isMaster()) {
    // Master gathers results
    all_results = local_results;

    for (int src = 1; src < size_; ++src) {
      // Receive result count
      int result_count;
      MPI_Recv(&result_count, 1, MPI_INT, src, 0, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      // Receive results
      for (int i = 0; i < result_count; ++i) {
        MotifResult result;

        // Receive pattern
        int pattern_len;
        MPI_Recv(&pattern_len, 1, MPI_INT, src, 1, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        result.motif_pattern.resize(pattern_len);
        MPI_Recv(&result.motif_pattern[0], pattern_len, MPI_CHAR, src, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Receive counts
        MPI_Recv(&result.match_count, 1, MPI_UNSIGNED_LONG, src, 3,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&result.frequency, 1, MPI_DOUBLE, src, 4, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);

        all_results.push_back(result);
      }
    }
  } else {
    // Workers send results
    int result_count = local_results.size();
    MPI_Send(&result_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

    for (const auto &result : local_results) {
      // Send pattern
      int pattern_len = result.motif_pattern.length();
      MPI_Send(&pattern_len, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
      MPI_Send(result.motif_pattern.c_str(), pattern_len, MPI_CHAR, 0, 2,
               MPI_COMM_WORLD);

      // Send counts
      MPI_Send(&result.match_count, 1, MPI_UNSIGNED_LONG, 0, 3, MPI_COMM_WORLD);
      MPI_Send(&result.frequency, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
    }
  }

  double comm_time = timer.elapsed();
  updateCommStats("gather_results", local_results.size() * sizeof(MotifResult),
                  comm_time);

  return all_results;
}

void MPIManager::synchronize() { MPI_Barrier(MPI_COMM_WORLD); }

std::unordered_map<std::string, double>
MPIManager::getCommunicationStats() const {
  return comm_stats_;
}

std::pair<size_t, size_t>
MPIManager::calculateWorkDistribution(size_t total_sequences, int process_rank,
                                      int total_processes) const {
  size_t base_work = total_sequences / total_processes;
  size_t extra_work = total_sequences % total_processes;

  size_t start_idx = process_rank * base_work +
                     std::min(static_cast<size_t>(process_rank), extra_work);
  size_t count =
      base_work + (process_rank < static_cast<int>(extra_work) ? 1 : 0);

  return {start_idx, count};
}

void MPIManager::updateCommStats(const std::string &operation, size_t bytes,
                                 double time_seconds) {
  comm_stats_[operation + "_bytes"] += bytes;
  comm_stats_[operation + "_time"] += time_seconds;
}

} // namespace dna_motif
