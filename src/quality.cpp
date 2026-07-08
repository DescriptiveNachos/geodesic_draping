#include "geodesic_draping/quality.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace geodesic_draping {
namespace {

struct BboxXY {
  double minX = 0.0;
  double maxX = 0.0;
  double minY = 0.0;
  double maxY = 0.0;
};

double triangleArea(const Vec3& p0, const Vec3& p1, const Vec3& p2) {
  return 0.5 * (p1 - p0).cross(p2 - p0).norm();
}

double percentile(std::vector<double> values, double fraction) {
  if (values.empty()) {
    return 0.0;
  }
  std::sort(values.begin(), values.end());
  const double clamped = std::clamp(fraction, 0.0, 1.0);
  const size_t index = static_cast<size_t>(std::round(clamped * static_cast<double>(values.size() - 1)));
  return values[index];
}

BboxXY computeBboxXY(const SurfaceMeshData& mesh) {
  if (mesh.vertices.empty()) {
    throw std::runtime_error("quality analysis requires at least one vertex");
  }

  BboxXY bbox;
  bbox.minX = bbox.maxX = mesh.vertices.front().x();
  bbox.minY = bbox.maxY = mesh.vertices.front().y();
  for (const Vec3& vertex : mesh.vertices) {
    bbox.minX = std::min(bbox.minX, vertex.x());
    bbox.maxX = std::max(bbox.maxX, vertex.x());
    bbox.minY = std::min(bbox.minY, vertex.y());
    bbox.maxY = std::max(bbox.maxY, vertex.y());
  }
  return bbox;
}

void raiseLevel(SolveQualityReport& report, SolveQualityLevel level, const std::string& warning) {
  if (static_cast<int>(level) > static_cast<int>(report.level)) {
    report.level = level;
  }
  report.warnings.push_back(warning);
}

std::string formatWarning(const std::string& label, double value, double threshold) {
  std::ostringstream out;
  out << label << " " << value << " exceeds threshold " << threshold;
  return out.str();
}

std::string trim(const std::string& text) {
  const size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

bool setThresholdValue(SolveQualityThresholds& thresholds, const std::string& key, double value) {
  if (key == "edgeLengthP99ToMedianWarning") {
    thresholds.edgeLengthP99ToMedianWarning = value;
  } else if (key == "edgeLengthP99ToMedianPoor") {
    thresholds.edgeLengthP99ToMedianPoor = value;
  } else if (key == "triangleAspectRatioWarning") {
    thresholds.triangleAspectRatioWarning = value;
  } else if (key == "triangleAspectRatioPoor") {
    thresholds.triangleAspectRatioPoor = value;
  } else if (key == "faceAreaRelativeToBboxWarning") {
    thresholds.faceAreaRelativeToBboxWarning = value;
  } else if (key == "faceAreaRelativeToBboxPoor") {
    thresholds.faceAreaRelativeToBboxPoor = value;
  } else if (key == "generatorMinLengthFraction") {
    thresholds.generatorMinLengthFraction = value;
  } else if (key == "generatorMaxLengthMultiplier") {
    thresholds.generatorMaxLengthMultiplier = value;
  } else if (key == "shearMaxMinusP99Warning") {
    thresholds.shearMaxMinusP99Warning = value;
  } else if (key == "shearMaxMinusP99Poor") {
    thresholds.shearMaxMinusP99Poor = value;
  } else if (key == "localShearJumpP99Warning") {
    thresholds.localShearJumpP99Warning = value;
  } else if (key == "localShearJumpP99Poor") {
    thresholds.localShearJumpP99Poor = value;
  } else if (key == "localShearJumpMaxWarning") {
    thresholds.localShearJumpMaxWarning = value;
  } else if (key == "localShearJumpMaxPoor") {
    thresholds.localShearJumpMaxPoor = value;
  } else {
    return false;
  }
  return true;
}

MeshQualityStats analyzeMeshQuality(const SurfaceMeshData& mesh) {
  if (mesh.faces.empty()) {
    throw std::runtime_error("quality analysis requires at least one face");
  }

  const BboxXY bbox = computeBboxXY(mesh);
  MeshQualityStats stats;
  stats.bboxWidth = bbox.maxX - bbox.minX;
  stats.bboxHeight = bbox.maxY - bbox.minY;
  const double bboxArea = std::max(stats.bboxWidth * stats.bboxHeight, std::numeric_limits<double>::epsilon());

  std::vector<double> edgeLengths;
  edgeLengths.reserve(mesh.faces.size() * 3);
  std::vector<double> aspectRatios;
  aspectRatios.reserve(mesh.faces.size());
  stats.minFaceAreaRelativeToBbox = std::numeric_limits<double>::infinity();

  for (const Face& face : mesh.faces) {
    const Vec3& p0 = mesh.vertices[face[0]];
    const Vec3& p1 = mesh.vertices[face[1]];
    const Vec3& p2 = mesh.vertices[face[2]];
    const std::array<double, 3> lengths{
        (p1 - p0).norm(),
        (p2 - p1).norm(),
        (p0 - p2).norm(),
    };
    edgeLengths.insert(edgeLengths.end(), lengths.begin(), lengths.end());

    const double area = triangleArea(p0, p1, p2);
    const double relativeArea = area / bboxArea;
    stats.minFaceAreaRelativeToBbox = std::min(stats.minFaceAreaRelativeToBbox, relativeArea);
    if (relativeArea <= 0.0) {
      ++stats.nearDegenerateFaceCount;
    }

    const double longest = *std::max_element(lengths.begin(), lengths.end());
    const double aspect = area > 0.0 ? (0.5 * longest * longest) / area
                                     : std::numeric_limits<double>::infinity();
    aspectRatios.push_back(aspect);
  }

  stats.edgeLengthMedian = percentile(edgeLengths, 0.50);
  stats.edgeLengthP99 = percentile(edgeLengths, 0.99);
  stats.edgeLengthP99ToMedian =
      stats.edgeLengthMedian > 0.0 ? stats.edgeLengthP99 / stats.edgeLengthMedian
                                   : std::numeric_limits<double>::infinity();
  stats.triangleAspectRatioMax = percentile(aspectRatios, 1.0);
  stats.triangleAspectRatioP99 = percentile(aspectRatios, 0.99);
  return stats;
}

GeneratorQualityStats analyzeGeneratorQuality(const SurfaceMeshData& mesh,
                                              const DrapeResult& result,
                                              const SolveQualityThresholds& thresholds) {
  if (!result.origin.extrinsicPoint) {
    throw std::runtime_error("generator quality analysis requires an extrinsic result origin");
  }
  const BboxXY bbox = computeBboxXY(mesh);
  const Vec3& seed = *result.origin.extrinsicPoint;
  const std::array<double, 4> wallDistances{
      seed.x() - bbox.minX,
      bbox.maxX - seed.x(),
      seed.y() - bbox.minY,
      bbox.maxY - seed.y(),
  };

  GeneratorQualityStats stats;
  stats.nearestWallDistance = *std::min_element(wallDistances.begin(), wallDistances.end());
  stats.farthestWallDistance = *std::max_element(wallDistances.begin(), wallDistances.end());
  stats.plausibleMinLength = thresholds.generatorMinLengthFraction * stats.nearestWallDistance;
  stats.plausibleMaxLength = thresholds.generatorMaxLengthMultiplier * stats.farthestWallDistance;

  stats.lengths.reserve(4);
  for (const TraceFamily& family : result.traces) {
    for (const DrapeTrace* trace : {&family.positive, &family.negative}) {
      stats.lengths.push_back(trace->length);
      if (!trace->hitBoundary) {
        ++stats.missedBoundaryCount;
      }
      if (trace->length < stats.plausibleMinLength) {
        ++stats.shortGeneratorCount;
      }
      if (trace->length > stats.plausibleMaxLength) {
        ++stats.longGeneratorCount;
      }
    }
  }
  return stats;
}

std::vector<std::pair<size_t, size_t>> scalarAdjacency(const SurfaceMeshData& mesh, DrapeSolveMode mode) {
  std::vector<std::pair<size_t, size_t>> pairs;
  if (mode == DrapeSolveMode::Complete) {
    std::set<std::pair<size_t, size_t>> edges;
    for (const Face& face : mesh.faces) {
      for (size_t i = 0; i < 3; ++i) {
        const size_t a = face[i];
        const size_t b = face[(i + 1) % 3];
        edges.insert(std::minmax(a, b));
      }
    }
    pairs.assign(edges.begin(), edges.end());
    return pairs;
  }

  std::map<std::pair<size_t, size_t>, size_t> edgeOwner;
  for (size_t faceIndex = 0; faceIndex < mesh.faces.size(); ++faceIndex) {
    const Face& face = mesh.faces[faceIndex];
    for (size_t i = 0; i < 3; ++i) {
      const auto edge = std::minmax(face[i], face[(i + 1) % 3]);
      const auto inserted = edgeOwner.emplace(edge, faceIndex);
      if (!inserted.second) {
        pairs.emplace_back(inserted.first->second, faceIndex);
      }
    }
  }
  return pairs;
}

ShearQualityStats analyzePrimaryShearQuality(const SurfaceMeshData& mesh, const DrapeResult& result) {
  const std::vector<double>* shear = nullptr;
  if (result.mode == DrapeSolveMode::Complete) {
    if (!result.vertexShearAnglesDegrees) {
      throw std::runtime_error("complete result is missing primary vertex shear");
    }
    shear = &*result.vertexShearAnglesDegrees;
  } else {
    if (!result.faceShearAnglesDegrees) {
      throw std::runtime_error("fast/hybrid result is missing primary face shear");
    }
    shear = &*result.faceShearAnglesDegrees;
  }

  ShearQualityStats stats;
  stats.sampleCount = shear->size();
  std::vector<double> finiteValues;
  finiteValues.reserve(shear->size());
  for (double value : *shear) {
    if (std::isfinite(value)) {
      finiteValues.push_back(value);
    } else {
      ++stats.nonFiniteCount;
    }
  }

  stats.p95 = percentile(finiteValues, 0.95);
  stats.p99 = percentile(finiteValues, 0.99);
  stats.max = percentile(finiteValues, 1.0);
  stats.maxMinusP99 = stats.max - stats.p99;

  std::vector<double> jumps;
  for (const auto& pair : scalarAdjacency(mesh, result.mode)) {
    const double a = (*shear)[pair.first];
    const double b = (*shear)[pair.second];
    if (std::isfinite(a) && std::isfinite(b)) {
      jumps.push_back(std::abs(a - b));
    }
  }
  stats.localJumpP95 = percentile(jumps, 0.95);
  stats.localJumpP99 = percentile(jumps, 0.99);
  stats.localJumpMax = percentile(jumps, 1.0);
  return stats;
}

void evaluateMesh(SolveQualityReport& report, const SolveQualityThresholds& thresholds) {
  const MeshQualityStats& mesh = report.mesh;
  if (mesh.edgeLengthP99ToMedian > thresholds.edgeLengthP99ToMedianPoor) {
    raiseLevel(report, SolveQualityLevel::Poor,
               formatWarning("edge length p99/median", mesh.edgeLengthP99ToMedian, thresholds.edgeLengthP99ToMedianPoor));
  } else if (mesh.edgeLengthP99ToMedian > thresholds.edgeLengthP99ToMedianWarning) {
    raiseLevel(report, SolveQualityLevel::Warning,
               formatWarning("edge length p99/median", mesh.edgeLengthP99ToMedian, thresholds.edgeLengthP99ToMedianWarning));
  }

  if (mesh.triangleAspectRatioMax > thresholds.triangleAspectRatioPoor) {
    raiseLevel(report, SolveQualityLevel::Poor,
               formatWarning("max triangle aspect ratio", mesh.triangleAspectRatioMax, thresholds.triangleAspectRatioPoor));
  } else if (mesh.triangleAspectRatioMax > thresholds.triangleAspectRatioWarning) {
    raiseLevel(report, SolveQualityLevel::Warning,
               formatWarning("max triangle aspect ratio", mesh.triangleAspectRatioMax, thresholds.triangleAspectRatioWarning));
  }

  if (mesh.minFaceAreaRelativeToBbox < thresholds.faceAreaRelativeToBboxPoor) {
    raiseLevel(report, SolveQualityLevel::Poor,
               "minimum face area relative to XY bounding box is extremely small");
  } else if (mesh.minFaceAreaRelativeToBbox < thresholds.faceAreaRelativeToBboxWarning) {
    raiseLevel(report, SolveQualityLevel::Warning,
               "minimum face area relative to XY bounding box is very small");
  }
}

void evaluateGenerators(SolveQualityReport& report) {
  const GeneratorQualityStats& generators = report.generators;
  if (generators.missedBoundaryCount > 0) {
    raiseLevel(report, SolveQualityLevel::Warning,
               "one or more generators did not hit the mesh boundary");
  }
  if (generators.shortGeneratorCount > 0) {
    raiseLevel(report, SolveQualityLevel::Warning,
               "one or more generators are short relative to seed-to-wall distances");
  }
  if (generators.longGeneratorCount > 0) {
    raiseLevel(report, SolveQualityLevel::Warning,
               "one or more generators are long relative to seed-to-wall distances");
  }
}

void evaluateShear(SolveQualityReport& report, const SolveQualityThresholds& thresholds) {
  const ShearQualityStats& shear = report.primaryShear;
  if (shear.nonFiniteCount > 0) {
    raiseLevel(report, SolveQualityLevel::Poor,
               "primary shear contains non-finite values");
  }

  if (shear.maxMinusP99 > thresholds.shearMaxMinusP99Poor) {
    raiseLevel(report, SolveQualityLevel::Poor,
               formatWarning("primary shear max-p99", shear.maxMinusP99, thresholds.shearMaxMinusP99Poor));
  } else if (shear.maxMinusP99 > thresholds.shearMaxMinusP99Warning) {
    raiseLevel(report, SolveQualityLevel::Warning,
               formatWarning("primary shear max-p99", shear.maxMinusP99, thresholds.shearMaxMinusP99Warning));
  }

  if (shear.localJumpP99 > thresholds.localShearJumpP99Poor) {
    raiseLevel(report, SolveQualityLevel::Poor,
               formatWarning("primary shear local-jump p99", shear.localJumpP99, thresholds.localShearJumpP99Poor));
  } else if (shear.localJumpP99 > thresholds.localShearJumpP99Warning) {
    raiseLevel(report, SolveQualityLevel::Warning,
               formatWarning("primary shear local-jump p99", shear.localJumpP99, thresholds.localShearJumpP99Warning));
  }

  if (shear.localJumpMax > thresholds.localShearJumpMaxPoor) {
    raiseLevel(report, SolveQualityLevel::Poor,
               formatWarning("primary shear local-jump max", shear.localJumpMax, thresholds.localShearJumpMaxPoor));
  } else if (shear.localJumpMax > thresholds.localShearJumpMaxWarning) {
    raiseLevel(report, SolveQualityLevel::Warning,
               formatWarning("primary shear local-jump max", shear.localJumpMax, thresholds.localShearJumpMaxWarning));
  }
}

} // namespace

const char* toString(SolveQualityLevel level) {
  switch (level) {
    case SolveQualityLevel::Good:
      return "Good";
    case SolveQualityLevel::Warning:
      return "Warning";
    case SolveQualityLevel::Poor:
      return "Poor";
  }
  return "Unknown";
}

SolveQualityThresholds loadSolveQualityThresholds(const std::filesystem::path& path,
                                                  SolveQualityThresholds thresholds) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("could not open quality threshold config: " + path.string());
  }

  std::string line;
  size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    const size_t comment = line.find('#');
    if (comment != std::string::npos) {
      line = line.substr(0, comment);
    }
    line = trim(line);
    if (line.empty()) {
      continue;
    }

    const size_t equals = line.find('=');
    if (equals == std::string::npos) {
      throw std::runtime_error("quality threshold config line " + std::to_string(lineNumber) +
                               " is missing '='");
    }

    const std::string key = trim(line.substr(0, equals));
    const std::string valueText = trim(line.substr(equals + 1));
    if (key.empty() || valueText.empty()) {
      throw std::runtime_error("quality threshold config line " + std::to_string(lineNumber) +
                               " has an empty key or value");
    }

    size_t parsedLength = 0;
    const double value = std::stod(valueText, &parsedLength);
    if (parsedLength != valueText.size()) {
      throw std::runtime_error("quality threshold config line " + std::to_string(lineNumber) +
                               " has a non-numeric value");
    }
    if (!setThresholdValue(thresholds, key, value)) {
      throw std::runtime_error("quality threshold config line " + std::to_string(lineNumber) +
                               " uses unknown key '" + key + "'");
    }
  }
  return thresholds;
}

SolveQualityReport analyzeSolveQuality(const SurfaceMeshData& mesh,
                                       const DrapeResult& result,
                                       const SolveQualityThresholds& thresholds) {
  SolveQualityReport report;
  report.mesh = analyzeMeshQuality(mesh);
  report.generators = analyzeGeneratorQuality(mesh, result, thresholds);
  report.primaryShear = analyzePrimaryShearQuality(mesh, result);
  evaluateMesh(report, thresholds);
  evaluateGenerators(report);
  evaluateShear(report, thresholds);
  return report;
}

} // namespace geodesic_draping
