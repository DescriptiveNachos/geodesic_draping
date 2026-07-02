#pragma once

#include "geodesic_draping/geodrape.h"

#include <cstddef>
#include <string>
#include <vector>

namespace geodesic_draping {

enum class SolveQualityLevel {
  Good,
  Warning,
  Poor,
};

struct SolveQualityThresholds {
  double edgeLengthP99ToMedianWarning = 8.0;
  double edgeLengthP99ToMedianPoor = 15.0;
  double triangleAspectRatioWarning = 50.0;
  double triangleAspectRatioPoor = 150.0;
  double faceAreaRelativeToBboxWarning = 1e-12;
  double faceAreaRelativeToBboxPoor = 1e-14;
  double generatorMinLengthFraction = 0.20;
  double generatorMaxLengthMultiplier = 2.50;
  double shearMaxMinusP99Warning = 25.0;
  double shearMaxMinusP99Poor = 45.0;
  double localShearJumpP99Warning = 25.0;
  double localShearJumpP99Poor = 45.0;
  double localShearJumpMaxWarning = 45.0;
  double localShearJumpMaxPoor = 70.0;
};

struct MeshQualityStats {
  double bboxWidth = 0.0;
  double bboxHeight = 0.0;
  double edgeLengthMedian = 0.0;
  double edgeLengthP99 = 0.0;
  double edgeLengthP99ToMedian = 0.0;
  double triangleAspectRatioMax = 0.0;
  double triangleAspectRatioP99 = 0.0;
  double minFaceAreaRelativeToBbox = 0.0;
  size_t nearDegenerateFaceCount = 0;
};

struct GeneratorQualityStats {
  double nearestWallDistance = 0.0;
  double farthestWallDistance = 0.0;
  double plausibleMinLength = 0.0;
  double plausibleMaxLength = 0.0;
  std::vector<double> lengths;
  size_t shortGeneratorCount = 0;
  size_t longGeneratorCount = 0;
  size_t missedBoundaryCount = 0;
};

struct ShearQualityStats {
  size_t sampleCount = 0;
  size_t nonFiniteCount = 0;
  double p95 = 0.0;
  double p99 = 0.0;
  double max = 0.0;
  double maxMinusP99 = 0.0;
  double localJumpP95 = 0.0;
  double localJumpP99 = 0.0;
  double localJumpMax = 0.0;
};

struct SolveQualityReport {
  SolveQualityLevel level = SolveQualityLevel::Good;
  std::vector<std::string> warnings;
  MeshQualityStats mesh;
  GeneratorQualityStats generators;
  ShearQualityStats primaryShear;
};

const char* toString(SolveQualityLevel level);

SolveQualityReport analyzeSolveQuality(const SurfaceMeshData& mesh,
                                       const DrapeResult& result,
                                       const SolveQualityThresholds& thresholds = {});

} // namespace geodesic_draping
