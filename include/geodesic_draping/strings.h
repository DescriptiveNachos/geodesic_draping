#pragma once

#include "geodesic_draping/geodrape.h"

#include <stdexcept>
#include <string>

namespace geodesic_draping {

inline std::string drapeSolveModeName(DrapeSolveMode mode) {
  switch (mode) {
  case DrapeSolveMode::Fast:
    return "fast";
  case DrapeSolveMode::Hybrid:
    return "hybrid";
  case DrapeSolveMode::Complete:
    return "complete";
  }
  return "unknown";
}

inline DrapeSolveMode parseDrapeSolveMode(const std::string& value) {
  if (value == "fast") return DrapeSolveMode::Fast;
  if (value == "hybrid") return DrapeSolveMode::Hybrid;
  if (value == "complete") return DrapeSolveMode::Complete;
  throw std::invalid_argument("mode must be one of: 'fast', 'hybrid', 'complete'");
}

inline std::string retrievalDomainName(RetrievalDomain domain) {
  switch (domain) {
  case RetrievalDomain::Intrinsic:
    return "intrinsic";
  case RetrievalDomain::Extrinsic:
    return "extrinsic";
  case RetrievalDomain::Subdivision:
    return "subdivision";
  }
  return "unknown";
}

inline RetrievalDomain parseRetrievalDomain(const std::string& value) {
  if (value == "intrinsic") return RetrievalDomain::Intrinsic;
  if (value == "extrinsic") return RetrievalDomain::Extrinsic;
  if (value == "subdivision") return RetrievalDomain::Subdivision;
  throw std::invalid_argument("retrieval domain must be one of: 'intrinsic', 'extrinsic', 'subdivision'");
}

inline IntrinsicTriangulationBackend parseIntrinsicBackend(const std::string& value) {
  if (value == "signpost") return IntrinsicTriangulationBackend::Signpost;
  if (value == "integer" || value == "integer-coordinates") {
    return IntrinsicTriangulationBackend::IntegerCoordinates;
  }
  throw std::invalid_argument("intrinsic_backend must be one of: 'signpost', 'integer'");
}

inline RefinementMode parseRefinementMode(const std::string& value) {
  if (value == "none") return RefinementMode::None;
  if (value == "flip") return RefinementMode::DelaunayFlip;
  if (value == "refine") return RefinementMode::DelaunayRefine;
  throw std::invalid_argument("refinement must be one of: 'none', 'flip', 'refine'");
}

} // namespace geodesic_draping
