#pragma once

#include "geodesic_draping/mesh.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace geodesic_draping::fixture_io {

inline std::string readText(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open " + path.string());
  }
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

inline std::string arrayTextForKey(const std::string& text, const std::string& key) {
  const std::string quotedKey = "\"" + key + "\"";
  const size_t keyPos = text.find(quotedKey);
  if (keyPos == std::string::npos) {
    throw std::runtime_error("missing JSON key " + key);
  }
  const size_t arrayStart = text.find('[', keyPos + quotedKey.size());
  if (arrayStart == std::string::npos) {
    throw std::runtime_error("missing array for JSON key " + key);
  }

  int depth = 0;
  for (size_t i = arrayStart; i < text.size(); ++i) {
    if (text[i] == '[') {
      ++depth;
    } else if (text[i] == ']') {
      --depth;
      if (depth == 0) {
        return text.substr(arrayStart, i - arrayStart + 1);
      }
    }
  }
  throw std::runtime_error("unterminated array for JSON key " + key);
}

inline std::vector<double> numbersInArrayText(const std::string& array) {
  std::vector<double> values;
  const char* cursor = array.c_str();
  while (*cursor != '\0') {
    char* end = nullptr;
    const double value = std::strtod(cursor, &end);
    if (end != cursor) {
      values.push_back(value);
      cursor = end;
    } else {
      ++cursor;
    }
  }
  return values;
}

inline std::vector<double> numbersInArrayForKey(const std::string& text, const std::string& key) {
  return numbersInArrayText(arrayTextForKey(text, key));
}

inline double numberForKey(const std::string& text, const std::string& key) {
  const std::string quotedKey = "\"" + key + "\"";
  const size_t keyPos = text.find(quotedKey);
  if (keyPos == std::string::npos) {
    throw std::runtime_error("missing JSON key " + key);
  }
  const size_t colon = text.find(':', keyPos + quotedKey.size());
  if (colon == std::string::npos) {
    throw std::runtime_error("missing JSON value for key " + key);
  }
  const char* cursor = text.c_str() + colon + 1;
  char* end = nullptr;
  const double value = std::strtod(cursor, &end);
  if (end == cursor) {
    throw std::runtime_error("expected numeric JSON value for key " + key);
  }
  return value;
}

inline SurfaceMeshData loadMesh(const std::filesystem::path& fixtureDir) {
  const std::string text = readText(fixtureDir / "mesh.json");
  const std::vector<double> vertexValues = numbersInArrayForKey(text, "vertices");
  const std::vector<double> faceValues = numbersInArrayForKey(text, "faces");

  if (vertexValues.size() % 3 != 0 || faceValues.size() % 3 != 0) {
    throw std::runtime_error("mesh arrays are not triples in " + fixtureDir.string());
  }

  SurfaceMeshData mesh;
  mesh.vertices.reserve(vertexValues.size() / 3);
  for (size_t i = 0; i < vertexValues.size(); i += 3) {
    mesh.vertices.emplace_back(vertexValues[i], vertexValues[i + 1], vertexValues[i + 2]);
  }

  mesh.faces.reserve(faceValues.size() / 3);
  for (size_t i = 0; i < faceValues.size(); i += 3) {
    mesh.faces.push_back(Face{static_cast<size_t>(faceValues[i]),
                              static_cast<size_t>(faceValues[i + 1]),
                              static_cast<size_t>(faceValues[i + 2])});
  }
  return mesh;
}

inline Vec2 loadSeedXY(const std::filesystem::path& fixtureDir) {
  const std::vector<double> values = numbersInArrayForKey(readText(fixtureDir / "inputs.json"), "seed_xy");
  if (values.size() != 2) {
    throw std::runtime_error("seed_xy must have two entries");
  }
  return Vec2(values[0], values[1]);
}

inline double loadAngleDegrees(const std::filesystem::path& fixtureDir) {
  return numberForKey(readText(fixtureDir / "inputs.json"), "angle_degrees");
}

} // namespace geodesic_draping::fixture_io
