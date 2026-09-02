#include "geodesic_draping/geodrape.h"
#include "geodesic_draping/strings.h"

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
namespace gd = geodesic_draping;
namespace gcs = geometrycentral::surface;

namespace {

[[noreturn]] void throwNotImplemented(const std::string& message) {
  PyErr_SetString(PyExc_NotImplementedError, message.c_str());
  throw py::error_already_set();
}

gd::RetrievalDomain parseRetrieval(const std::string& value) {
  if (value == "intrinsic") {
    throwNotImplemented("retrieval='intrinsic' is not exposed in Python v1");
  }
  if (value == "extrinsic" || value == "subdivision") return gd::parseRetrievalDomain(value);
  throw py::value_error("retrieval must be one of: 'extrinsic', 'subdivision'");
}

geometrycentral::LevelSetConstraint parseLevelSetConstraint(const std::string& value) {
  if (value == "none") return geometrycentral::LevelSetConstraint::None;
  throw py::value_error("level_set_constraint must be 'none' in Python v1");
}

gd::SurfaceMeshData meshDataFromArrays(py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
                                       py::array_t<long long, py::array::c_style | py::array::forcecast> faces) {
  py::buffer_info vertexInfo = vertices.request();
  py::buffer_info faceInfo = faces.request();

  if (vertexInfo.ndim != 2 || vertexInfo.shape[1] != 3) {
    throw py::value_error("vertices must have shape (V, 3)");
  }
  if (faceInfo.ndim != 2 || faceInfo.shape[1] != 3) {
    throw py::value_error("faces must have shape (F, 3)");
  }

  const auto vertexView = vertices.unchecked<2>();
  const auto faceView = faces.unchecked<2>();
  gd::SurfaceMeshData mesh;
  mesh.vertices.reserve(static_cast<size_t>(vertexInfo.shape[0]));
  for (py::ssize_t i = 0; i < vertexInfo.shape[0]; ++i) {
    mesh.vertices.emplace_back(vertexView(i, 0), vertexView(i, 1), vertexView(i, 2));
  }

  mesh.faces.reserve(static_cast<size_t>(faceInfo.shape[0]));
  for (py::ssize_t i = 0; i < faceInfo.shape[0]; ++i) {
    std::array<size_t, 3> face{};
    for (py::ssize_t j = 0; j < 3; ++j) {
      const long long index = faceView(i, j);
      if (index < 0 || static_cast<size_t>(index) >= mesh.vertices.size()) {
        throw py::value_error("faces contain an out-of-range vertex index");
      }
      face[static_cast<size_t>(j)] = static_cast<size_t>(index);
    }
    mesh.faces.push_back(face);
  }

  return mesh;
}

gd::Vec2 vec2FromArray(py::array_t<double, py::array::c_style | py::array::forcecast> values,
                       const std::string& name) {
  py::buffer_info info = values.request();
  if (info.ndim != 1 || info.shape[0] != 2) {
    throw py::value_error(name + " must have shape (2,)");
  }
  const auto view = values.unchecked<1>();
  return gd::Vec2(view(0), view(1));
}

gd::SignedHeatSolveOptions makeHeatOptions(bool preserveSourceNormals,
                                           const std::string& levelSetConstraint,
                                           double softLevelSetWeight,
                                           double diffusionTimeCoefficient) {
  gd::SignedHeatSolveOptions options;
  options.preserveSourceNormals = preserveSourceNormals;
  options.levelSetConstraint = parseLevelSetConstraint(levelSetConstraint);
  options.softLevelSetWeight = softLevelSetWeight;
  options.diffusionTimeCoefficient = diffusionTimeCoefficient;
  return options;
}

gd::IntrinsicConstructionOptions makeIntrinsicOptions(const std::string& backend) {
  gd::IntrinsicConstructionOptions options;
  options.backend = gd::parseIntrinsicBackend(backend);
  return options;
}

gd::RefinementOptions makeRefinementOptions(const std::string& refinement,
                                            std::optional<double> angleThreshold,
                                            std::optional<double> circumradiusThreshold,
                                            std::optional<size_t> maxInsertions) {
  gd::RefinementOptions options;
  options.mode = gd::parseRefinementMode(refinement);
  options.angleThreshold = angleThreshold;
  options.circumradiusThreshold = circumradiusThreshold;
  options.maxInsertions = maxInsertions;
  return options;
}

gd::DrapeSolveOptions makeSolveOptions(const std::string& mode,
                                       double fiberAngle,
                                       std::optional<double> traceLength,
                                       std::optional<size_t> maxTraceIterations) {
  gd::DrapeSolveOptions options;
  options.mode = gd::parseDrapeSolveMode(mode);
  options.fiberAngle = fiberAngle;
  options.advanced.trace.traceLength = traceLength;
  options.advanced.trace.maxIterations = maxTraceIterations;
  return options;
}

gd::RetrievalOptions makeRetrievalOptions(const std::string& retrieval, bool sampleVertexShear) {
  gd::RetrievalOptions options;
  options.domain = parseRetrieval(retrieval);
  options.sampleVertexShear = sampleVertexShear;
  return options;
}

py::array_t<double> vector3Array(const std::vector<geometrycentral::Vector3>& values) {
  py::array_t<double> array({static_cast<py::ssize_t>(values.size()), py::ssize_t{3}});
  auto out = array.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(values.size()); ++i) {
    out(i, 0) = values[static_cast<size_t>(i)].x;
    out(i, 1) = values[static_cast<size_t>(i)].y;
    out(i, 2) = values[static_cast<size_t>(i)].z;
  }
  return array;
}

py::array_t<double> verticesArray(const gd::DrapeResult& result) {
  if (!result.mesh) {
    throw std::runtime_error("result does not contain a mesh");
  }
  auto& mesh = const_cast<gcs::SurfaceMesh&>(*result.mesh);
  py::array_t<double> array({static_cast<py::ssize_t>(mesh.nVertices()), py::ssize_t{3}});
  auto out = array.mutable_unchecked<2>();

  if (result.vertexPositions) {
    for (gcs::Vertex vertex : mesh.vertices()) {
      const geometrycentral::Vector3 value = (*result.vertexPositions)[vertex];
      const py::ssize_t i = static_cast<py::ssize_t>(vertex.getIndex());
      out(i, 0) = value.x;
      out(i, 1) = value.y;
      out(i, 2) = value.z;
    }
    return array;
  }

  if (result.extrinsicGeometry) {
    auto& geometry = const_cast<gcs::VertexPositionGeometry&>(*result.extrinsicGeometry);
    for (gcs::Vertex vertex : mesh.vertices()) {
      const geometrycentral::Vector3 value = geometry.vertexPositions[vertex];
      const py::ssize_t i = static_cast<py::ssize_t>(vertex.getIndex());
      out(i, 0) = value.x;
      out(i, 1) = value.y;
      out(i, 2) = value.z;
    }
    return array;
  }

  throw std::runtime_error("result has no Python-exportable vertex positions");
}

py::array_t<long long> facesArray(const gd::DrapeResult& result) {
  if (!result.mesh) {
    throw std::runtime_error("result does not contain a mesh");
  }
  auto& mesh = const_cast<gcs::SurfaceMesh&>(*result.mesh);
  py::array_t<long long> array({static_cast<py::ssize_t>(mesh.nFaces()), py::ssize_t{3}});
  auto out = array.mutable_unchecked<2>();
  for (gcs::Face face : mesh.faces()) {
    const py::ssize_t i = static_cast<py::ssize_t>(face.getIndex());
    py::ssize_t j = 0;
    for (gcs::Vertex vertex : face.adjacentVertices()) {
      if (j >= 3) {
        throw std::runtime_error("Python export expects triangular faces");
      }
      out(i, j++) = static_cast<long long>(vertex.getIndex());
    }
    if (j != 3) {
      throw std::runtime_error("Python export expects triangular faces");
    }
  }
  return array;
}

py::array_t<double> faceScalarArray(const gcs::FaceData<double>& field) {
  auto& mesh = *field.getMesh();
  py::array_t<double> array(static_cast<py::ssize_t>(mesh.nFaces()));
  auto out = array.mutable_unchecked<1>();
  for (gcs::Face face : mesh.faces()) {
    out(static_cast<py::ssize_t>(face.getIndex())) = field[face];
  }
  return array;
}

py::array_t<double> vertexScalarArray(const gcs::VertexData<double>& field) {
  auto& mesh = *field.getMesh();
  py::array_t<double> array(static_cast<py::ssize_t>(mesh.nVertices()));
  auto out = array.mutable_unchecked<1>();
  for (gcs::Vertex vertex : mesh.vertices()) {
    out(static_cast<py::ssize_t>(vertex.getIndex())) = field[vertex];
  }
  return array;
}

py::array_t<double> distancesArray(const std::array<gcs::VertexData<double>, 2>& fields) {
  auto& mesh = *fields[0].getMesh();
  py::array_t<double> array({py::ssize_t{2}, static_cast<py::ssize_t>(mesh.nVertices())});
  auto out = array.mutable_unchecked<2>();
  for (size_t family = 0; family < fields.size(); ++family) {
    for (gcs::Vertex vertex : mesh.vertices()) {
      out(static_cast<py::ssize_t>(family), static_cast<py::ssize_t>(vertex.getIndex())) =
          fields[family][vertex];
    }
  }
  return array;
}

py::array_t<double> directionFieldsArray(
    const std::array<gcs::FaceData<geometrycentral::Vector3>, 2>& fields) {
  auto& mesh = *fields[0].getMesh();
  py::array_t<double> array(
      {py::ssize_t{2}, static_cast<py::ssize_t>(mesh.nFaces()), py::ssize_t{3}});
  auto out = array.mutable_unchecked<3>();
  for (size_t family = 0; family < fields.size(); ++family) {
    for (gcs::Face face : mesh.faces()) {
      const geometrycentral::Vector3 value = fields[family][face];
      const py::ssize_t i = static_cast<py::ssize_t>(face.getIndex());
      out(static_cast<py::ssize_t>(family), i, 0) = value.x;
      out(static_cast<py::ssize_t>(family), i, 1) = value.y;
      out(static_cast<py::ssize_t>(family), i, 2) = value.z;
    }
  }
  return array;
}

py::list generatorArrays(const gd::DrapeResult& result) {
  py::list families;
  if (!result.extrinsicGenerators) {
    families.append(py::list());
    families.append(py::list());
    return families;
  }

  for (size_t family = 0; family < 2; ++family) {
    py::list pair;
    pair.append(vector3Array((*result.extrinsicGenerators)[2 * family].points));
    pair.append(vector3Array((*result.extrinsicGenerators)[2 * family + 1].points));
    families.append(pair);
  }
  return families;
}

py::dict resultDict(const gd::DrapeResult& result) {
  py::dict dict;
  dict["vertices"] = verticesArray(result);
  dict["faces"] = facesArray(result);
  dict["domain"] = gd::retrievalDomainName(result.domain);
  dict["mode"] = gd::drapeSolveModeName(result.mode);
  dict["generators"] = generatorArrays(result);
  dict["direction_fields"] = result.directionFields
                                  ? py::object(directionFieldsArray(*result.directionFields))
                                  : py::object(py::none());
  dict["distances"] = result.distances
                          ? py::object(distancesArray(*result.distances))
                          : py::object(py::none());
  dict["face_shear"] = result.faceShear
                           ? py::object(faceScalarArray(*result.faceShear))
                           : py::object(py::none());
  dict["vertex_shear"] = result.vertexShear
                             ? py::object(vertexScalarArray(*result.vertexShear))
                             : py::object(py::none());
  return dict;
}

class PyGeoDrapeSolver {
public:
  PyGeoDrapeSolver(py::array_t<double, py::array::c_style | py::array::forcecast> vertices,
                   py::array_t<long long, py::array::c_style | py::array::forcecast> faces,
                   const std::string& intrinsicBackend,
                   const std::string& refinement,
                   std::optional<double> angleThreshold,
                   std::optional<double> circumradiusThreshold,
                   std::optional<size_t> maxInsertions,
                   bool preserveSourceNormals,
                   const std::string& levelSetConstraint,
                   double softLevelSetWeight,
                   double diffusionTimeCoefficient)
      : solver_(meshDataFromArrays(vertices, faces),
                makeHeatOptions(preserveSourceNormals,
                                levelSetConstraint,
                                softLevelSetWeight,
                                diffusionTimeCoefficient),
                makeIntrinsicOptions(intrinsicBackend),
                makeRefinementOptions(refinement,
                                      angleThreshold,
                                      circumradiusThreshold,
                                      maxInsertions)) {}

  py::dict solve(py::array_t<double, py::array::c_style | py::array::forcecast> seedXY,
                 double fabricAngle,
                 const std::string& mode,
                 double fiberAngle,
                 const std::string& retrieval,
                 bool sampleVertexShear,
                 std::optional<double> traceLength,
                 std::optional<size_t> maxTraceIterations) {
    const gd::DrapeSolveOptions solveOptions =
        makeSolveOptions(mode, fiberAngle, traceLength, maxTraceIterations);
    const gd::RetrievalOptions retrievalOptions =
        makeRetrievalOptions(retrieval, sampleVertexShear);
    const gd::Vec2 seed = vec2FromArray(seedXY, "seed_xy");
    gd::DrapeResult result;
    {
      py::gil_scoped_release release;
      result = solver_.solve(seed, fabricAngle, solveOptions, retrievalOptions);
    }
    return resultDict(result);
  }

  py::dict retrieve(const std::string& retrieval, bool sampleVertexShear) {
    const gd::RetrievalOptions retrievalOptions =
        makeRetrievalOptions(retrieval, sampleVertexShear);
    gd::DrapeResult result;
    {
      py::gil_scoped_release release;
      result = solver_.retrieve(retrievalOptions);
    }
    return resultDict(result);
  }

private:
  gd::GeoDrapeSolver solver_;
};

} // namespace

PYBIND11_MODULE(_core, m) {
  py::class_<PyGeoDrapeSolver>(m, "GeoDrapeSolver")
      .def(py::init<py::array_t<double, py::array::c_style | py::array::forcecast>,
                    py::array_t<long long, py::array::c_style | py::array::forcecast>,
                    const std::string&,
                    const std::string&,
                    std::optional<double>,
                    std::optional<double>,
                    std::optional<size_t>,
                    bool,
                    const std::string&,
                    double,
                    double>(),
           py::arg("vertices"),
           py::arg("faces"),
           py::kw_only(),
           py::arg("intrinsic_backend") = "signpost",
           py::arg("refinement") = "none",
           py::arg("angle_threshold") = py::none(),
           py::arg("circumradius_threshold") = py::none(),
           py::arg("max_insertions") = py::none(),
           py::arg("preserve_source_normals") = false,
           py::arg("level_set_constraint") = "none",
           py::arg("soft_level_set_weight") = -1.0,
           py::arg("diffusion_time_coefficient") = 1.0)
      .def("solve",
           &PyGeoDrapeSolver::solve,
           py::arg("seed_xy"),
           py::arg("fabric_angle"),
           py::kw_only(),
           py::arg("mode") = "complete",
           py::arg("fiber_angle") = 90.0,
           py::arg("retrieval") = "extrinsic",
           py::arg("sample_vertex_shear") = false,
           py::arg("trace_length") = py::none(),
           py::arg("max_trace_iterations") = py::none())
      .def("retrieve",
           &PyGeoDrapeSolver::retrieve,
           py::kw_only(),
           py::arg("retrieval") = "extrinsic",
           py::arg("sample_vertex_shear") = false);
}
