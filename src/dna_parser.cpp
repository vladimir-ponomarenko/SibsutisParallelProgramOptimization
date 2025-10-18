#include "dna_parser.h"
#include <algorithm>
#include <filesystem>
#include <format>
#include <iomanip>
#include <ranges>
#include <set>
#include <stdexcept>

namespace dna_motif {

ParseResult<std::vector<ChIPSequence>>
DNAParser::parseChIPSequences(std::string_view filename) {
  if (!isFileReadable(filename)) {
    return std::unexpected(ParseError::FileNotFound);
  }

  // Read file content
  auto file_content = readFile(filename);
  if (!file_content) {
    return std::unexpected(file_content.error());
  }

  updateStats("files_opened");

  auto lines = splitLines(*file_content);

  std::vector<ChIPSequence> sequences;
  std::string current_header;
  std::vector<std::string> current_sequence_lines;

  for (const auto &line : lines) {
    const auto trimmed_line = trim(line);

    if (trimmed_line.empty()) {
      continue;
    }

    if (trimmed_line[0] == '>') {
      if (!current_header.empty() && !current_sequence_lines.empty()) {
        try {
          ChIPSequence seq =
              parseChIPSequence(current_header, current_sequence_lines);
          if (validateSequence(seq.sequence)) {
            sequences.push_back(std::move(seq));
            updateStats("sequences_parsed");
          } else {
            updateStats("sequences_invalid");
          }
        } catch (const std::exception &e) {
          updateStats("sequences_parse_errors");
          std::cerr << std::format("Warning: Failed to parse sequence: {}\n",
                                   e.what());
        }
      }

      current_header = trimmed_line;
      current_sequence_lines.clear();
    } else {
      current_sequence_lines.push_back(trimmed_line);
    }
  }

  if (!current_header.empty() && !current_sequence_lines.empty()) {
    try {
      ChIPSequence seq =
          parseChIPSequence(current_header, current_sequence_lines);
      if (validateSequence(seq.sequence)) {
        sequences.push_back(std::move(seq));
        updateStats("sequences_parsed");
      } else {
        updateStats("sequences_invalid");
      }
    } catch (const std::exception &e) {
      updateStats("sequences_parse_errors");
      std::cerr << std::format("Warning: Failed to parse sequence: {}\n",
                               e.what());
    }
  }

  updateStats("files_closed");

  return sequences;
}

ParseResult<std::vector<Motif>>
DNAParser::parseMotifs(std::string_view filename) {
  if (!isFileReadable(filename)) {
    return std::unexpected(ParseError::FileNotFound);
  }

  auto file_content = readFile(filename);
  if (!file_content) {
    return std::unexpected(file_content.error());
  }

  updateStats("files_opened");

  auto lines = splitLines(*file_content);

  std::vector<Motif> motifs;

  for (const auto &line : lines) {
    const auto trimmed_line = trim(line);

    if (trimmed_line.empty() || trimmed_line[0] == '#') {
      continue;
    }

    try {
      Motif motif = parseMotifLine(trimmed_line);
      motifs.push_back(std::move(motif));
      updateStats("motifs_parsed");
    } catch (const std::exception &e) {
      updateStats("motifs_parse_errors");
      std::cerr << std::format("Warning: Failed to parse motif line: {} - {}\n",
                               trimmed_line, e.what());
    }
  }

  updateStats("files_closed");

  return motifs;
}

bool DNAParser::validateSequence(std::string_view sequence) const noexcept {
  if (sequence.empty()) {
    return false;
  }

  return std::ranges::all_of(sequence, [](char c) {
    const char upper_c = std::toupper(c);
    return upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C';
  });
}

bool DNAParser::isFileReadable(std::string_view filename) const noexcept {
  try {
    return std::filesystem::exists(filename) &&
           std::filesystem::is_regular_file(filename);
  } catch (...) {
    return false;
  }
}

size_t DNAParser::getFileSize(std::string_view filename) const noexcept {
  try {
    return std::filesystem::file_size(filename);
  } catch (...) {
    return 0;
  }
}

ParseResult<std::string> DNAParser::readFile(std::string_view filename) const {
  try {
    std::ifstream file{std::string(filename)};
    if (!file.is_open()) {
      return std::unexpected(ParseError::FileNotFound);
    }

    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());

    if (file.bad()) {
      return std::unexpected(ParseError::IOError);
    }

    return content;
  } catch (const std::exception &) {
    return std::unexpected(ParseError::IOError);
  }
}

std::vector<std::string> DNAParser::splitLines(std::string_view text) const {
  std::vector<std::string> lines;

  auto line_views =
      text | std::views::split('\n') | std::views::transform([](auto &&range) {
        return std::string(range.begin(), range.end());
      });

  std::ranges::copy(line_views, std::back_inserter(lines));
  return lines;
}

std::pair<std::string, std::vector<std::string>>
DNAParser::parseHeader(std::string_view header_line) const {
  auto header_parts = split(header_line, '\t');

  if (header_parts.empty()) {
    throw std::runtime_error("Invalid header format");
  }

  // Remove '>' from the beginning
  std::string id = header_parts[0].substr(1);

  std::vector<std::string> metadata;
  std::ranges::copy(header_parts | std::views::drop(1),
                    std::back_inserter(metadata));

  return {id, metadata};
}

ChIPSequence
DNAParser::parseChIPSequence(std::string_view header_line,
                             std::span<const std::string> sequence_lines) {
  auto [id, metadata] = parseHeader(header_line);

  std::string sequence = cleanSequence(sequence_lines);

  ChIPSequence chip_seq(id, sequence);
  chip_seq.metadata = std::move(metadata);

  return chip_seq;
}

std::string
DNAParser::cleanSequence(std::span<const std::string> sequence_lines) const {
  std::string sequence;
  sequence.reserve(sequence_lines.size() * 40);

  for (const auto &line : sequence_lines) {
    const auto clean_line = trim(line);
    auto cleaned = clean_line |
                   std::views::filter([](char c) { return !std::isspace(c); });
    std::ranges::copy(cleaned, std::back_inserter(sequence));
  }

  return sequence;
}

Motif DNAParser::parseMotifLine(std::string_view line) {
  auto parts = split(line, '\t');

  if (parts.size() < 4) {
    throw std::runtime_error(
        std::format("Invalid motif line format: {}", line));
  }

  std::string pattern = trim(parts[0]);
  double score1 = std::stod(parts[1]);
  double score2 = std::stod(parts[2]);
  double score3 = std::stod(parts[3]);

  return Motif(pattern, score1, score2, score3);
}

void DNAParser::updateStats(std::string_view key, size_t increment) noexcept {
  stats_[std::string(key)] += increment;
}

std::string DNAParser::errorToString(ParseError error) noexcept {
  switch (error) {
  case ParseError::FileNotFound:
    return "File not found";
  case ParseError::InvalidFormat:
    return "Invalid file format";
  case ParseError::InvalidSequence:
    return "Invalid DNA sequence";
  case ParseError::InvalidMotif:
    return "Invalid motif format";
  case ParseError::IOError:
    return "I/O error";
  case ParseError::Unknown:
  default:
    return "Unknown error";
  }
}

std::string trim(std::string_view str) {
  const auto first = str.find_first_not_of(" \t\n\r\f\v");
  if (first == std::string_view::npos) {
    return "";
  }
  const auto last = str.find_last_not_of(" \t\n\r\f\v");
  return std::string(str.substr(first, last - first + 1));
}

std::vector<std::string> split(std::string_view str, char delimiter) {
  std::vector<std::string> tokens;

  auto token_views = str | std::views::split(delimiter) |
                     std::views::transform([](auto &&range) {
                       return std::string(range.begin(), range.end());
                     });

  std::ranges::copy(token_views, std::back_inserter(tokens));
  return tokens;
}

bool isValidDNASequence(std::string_view sequence) {
  if (sequence.empty()) {
    return false;
  }

  return std::ranges::all_of(sequence, [](char c) {
    const char upper_c = std::toupper(c);
    return upper_c == 'A' || upper_c == 'T' || upper_c == 'G' || upper_c == 'C';
  });
}

bool isValidIUPACCode(char code) {
  const char upper_code = std::toupper(code);
  return VALID_DNA_NUCLEOTIDES.find(upper_code) != std::string_view::npos ||
         IUPAC_CODES.find(upper_code) != std::string_view::npos;
}

std::string toUpperCase(std::string_view str) {
  std::string result;
  result.reserve(str.size());

  std::ranges::transform(str, std::back_inserter(result),
                         [](char c) { return std::toupper(c); });

  return result;
}

std::string toLowerCase(std::string_view str) {
  std::string result;
  result.reserve(str.size());

  std::ranges::transform(str, std::back_inserter(result),
                         [](char c) { return std::tolower(c); });

  return result;
}

bool startsWith(std::string_view str, std::string_view prefix) {
  return str.size() >= prefix.size() && str.substr(0, prefix.size()) == prefix;
}

bool endsWith(std::string_view str, std::string_view suffix) {
  return str.size() >= suffix.size() &&
         str.substr(str.size() - suffix.size()) == suffix;
}

std::string replaceAll(std::string_view str, std::string_view from,
                       std::string_view to) {
  std::string result(str);
  size_t pos = 0;

  while ((pos = result.find(from, pos)) != std::string::npos) {
    result.replace(pos, from.length(), to);
    pos += to.length();
  }

  return result;
}

std::vector<std::string> splitLines(std::string_view str) {
  std::vector<std::string> lines;

  auto line_views =
      str | std::views::split('\n') | std::views::transform([](auto &&range) {
        return std::string(range.begin(), range.end());
      });

  std::ranges::copy(line_views, std::back_inserter(lines));
  return lines;
}

std::string join(const std::vector<std::string> &strings,
                 std::string_view delimiter) {
  if (strings.empty()) {
    return "";
  }

  std::string result = strings[0];
  for (size_t i = 1; i < strings.size(); ++i) {
    result += delimiter;
    result += strings[i];
  }

  return result;
}

std::string formatProgress(size_t current, size_t total,
                           std::string_view operation) {
  if (total == 0)
    return "";

  const double percentage = (static_cast<double>(current) / total) * 100.0;
  return std::format("{}: {}/{} ({:.1f}%)", operation, current, total,
                     percentage);
}

void printProgress(size_t current, size_t total, std::string_view operation) {
  if (total == 0)
    return;

  const auto progress = formatProgress(current, total, operation);
  std::cout << "\r" << progress << std::flush;

  if (current == total) {
    std::cout << std::endl;
  }
}

} // namespace dna_motif
