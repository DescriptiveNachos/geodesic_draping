# geodesic_draping

Compute geodesic draping approximations on triangle meshes as described in https://papers.ssrn.com/sol3/papers.cfm?abstract_id=6553829. The draping of a fabric is computed as the global crossing of two equidistant fiber families constrained by generators created as straightest geodesics from a drape origin and initial fiber directions.  

While the paper used a python implementation with a customized potpourri3D fork, the present code implements the solver directly in C++ and makes it available through python bindings. This allows the solver core to use intrinsic triangulations which can improve robustness through intrinsic edge flips and refinement.

The implementation supports different `intrinsic backends`, `retrieval domains` and `solve modes`. 
The solver core computes the draping on intrinsic geomtry. There are multiple valid choices to represent this domain. We currently provide a choice between geometry-centrals `signpost` and `integer coordinates`. The former is faster, the later provides added robustness, which is especially valuable when using `intrinsic refinement`. 
As the internal solver core works on intrinsic geometry there are multiple possible return domains. The direct return of the intrinsic domain is only available in C++. The bindings offer to return the result on the extrinsic input domain or the common subdivision of intrinsic and extrinsic domain. The former is the original input mesh, the latter also allows the return of face shear and direction fields.\
For solve modes both C++ and bindings allow `Fast`, `Hybrid` and `Complete` solves. `Fast` computes face shear from signed-heat direction fields. Does not
return distance fields. `Hybrid` computes face shear from signed-heat direction fields and also
integrates vertex distance fields. `Complete` integrates vertex distance fields and computes shear from the
complete solve path.\

The main public header is:

```cpp
#include "geodesic_draping/geodrape.h"
```

You can find example use case of this repo in https://github.com/DescriptiveNachos/geodesic_draping_visualizer which provides an interactive UI for visualizing the draping approximation on a mesh.

# Python

## Python Install

The Python package is built with `scikit-build-core`, the simplest install is via one of the available wheels Download a matching wheel from the GitHub release if one is available. Check compatible tags with
```powershell
python -m pip debug --verbose
```
and install with:

```powershell
python -m pip install geodesic_draping-0.1.0-cp312-cp312-win_amd64.whl
```

To install from source:

```powershell
git clone https://github.com/DescriptiveNachos/geodesic_draping.git
cd geodesic_draping
git submodule update --init --recursive
python -m pip install .
```

For editable development installs:

```powershell
python -m pip install -e .[test]
pytest tests/python
```

The Python extension is built from the same CMake target as the C++ library.
Polyscope support in the C++ debug viewer is disabled for Python wheels by
default.

## Python Examples

Small Python examples live in `examples/`, which show the basic usage and can be used as a quick install verification. To execute the examples just navigate to the repo root and run:

```powershell
python examples/basic_solve.py
python examples/persistent_solver.py
```

For a Polyscope debug views use:

```powershell
python -m pip install polyscope
python tools/plot_drape_result.py
python examples/subdivision_debug_plot.py
```

## Python API

solve_drape() provides a convenience one-shot solver:
```python
import numpy as np
import geodesic_draping as gd

vertices = np.asarray(..., dtype=float)  # shape (V, 3)
faces = np.asarray(..., dtype=np.int64)  # shape (F, 3)

result = gd.solve_drape(
    vertices,
    faces,
    seed_xy=np.array([0.0, 0.0]),
    fabric_angle=20.0,
    mode="complete",
    sample_vertex_shear=True,
)

print(result.distances.shape)
print(result.vertex_shear.shape)
```
`DrapeResult` is a dataclass-style object with Python-owned NumPy arrays:

```python
result.vertices         # (V, 3)
result.faces            # (F, 3)
result.domain           # "extrinsic" or "subdivision"
result.mode             # "fast", "hybrid", or "complete"
result.generators       # [[family0_plus, family0_minus], [family1_plus, family1_minus]]
result.direction_fields # None or (2, F, 3)
result.distances        # None or (2, V)
result.face_shear       # None or (F,)
result.vertex_shear     # None or (V,)
```

Use `GeoDrapeSolver` for repeated solves on the same mesh. This keeps the C++
solver and its factorizations alive and exelerates repeated solves drastically. 

```python
solver = gd.GeoDrapeSolver(
    vertices,
    faces,
    intrinsic_backend="integer",  # "signpost" or "integer"
    refinement="none",            # "none", "flip", or "refine"
)

fast = solver.solve(
    seed_xy=np.array([0.0, 0.0]),
    fabric_angle=20.0,
    mode="fast",
)

subdivision = solver.retrieve(retrieval="subdivision")
print(subdivision.face_shear.shape)
```

## Build A Wheel

```powershell
python -m pip install scikit-build-core pybind11 numpy
python -m pip wheel . --no-build-isolation -w dist
```

The project wheel should contain only the Python package, the compiled `_core`
extension, and package metadata.

# C++

## C++ Build

Source builds from the Git repository require the Geometry Central submodule.
Clone normally, then initialize Geometry Central:

```powershell
git clone https://github.com/DescriptiveNachos/geodesic_draping.git
cd geodesic_draping
git submodule update --init --recursive extern/geometry-central
```

If you want to use the debug plotting tests you will also need polyscope

```powershell
git submodule update --init --recursive extern/polyscope
```

The repo provides both basic tests and some tools for benchmarking and debug plots, disable anything you do not need, but make sure to enable polyscope if you want to use debug plotting tools. To get all provided tests and tools choose:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 `
  -DGEODESIC_DRAPING_BUILD_TESTS=ON `
  -DGEODESIC_DRAPING_BUILD_TOOLS=ON 
  -DGEODESIC_DRAPING_ENABLE_POLYSCOPE=ON

cmake --build build --config Release
```
If you enabled tests you may run the following to confirm a successful build:
```powershell
ctest --test-dir build --build-config Release --output-on-failure
```
If you also enabled tools and polyscope you can also run:
```powershell
.\build\Release\debug_drape_result.exe demo_part
```
and test the following useful options:

```powershell
--mode fast|hybrid|complete
--domain extrinsic|subdivision
--backend signpost|integer
--refinement none|flip|refine
--sample-vertex-shear
```

The active fixture is `test_data/fixtures/demo_part`.

## Mesh Input

```cpp
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Face = std::array<size_t, 3>;

struct SurfaceMeshData {
  std::vector<Vec3> vertices;
  std::vector<Face> faces;
};
```

`vertices` are 3D positions. `faces` are triangular vertex-index triples.

## One-Shot Solve

Use `solveDrape()` when you only need a single solve:

```cpp
geodesic_draping::SurfaceMeshData meshData = ...;
geodesic_draping::Vec2 seedXY(0.0, 0.0);
double fabricAngle = 20.0; // degrees

geodesic_draping::DrapeResult result =
    geodesic_draping::solveDrape(meshData, seedXY, fabricAngle);
```

The seed is supplied in input XY coordinates and projected to the mesh. Public
facing angles are degrees.

## Persistent Solver

Use `GeoDrapeSolver` for repeated solves on the same mesh. This keeps solver
state and lazy factorizations alive across calls.

```cpp
geodesic_draping::GeoDrapeSolver solver(meshData);

geodesic_draping::DrapeSolveOptions solveOptions;
solveOptions.mode = geodesic_draping::DrapeSolveMode::Complete;

geodesic_draping::RetrievalOptions retrievalOptions;
retrievalOptions.domain = geodesic_draping::RetrievalDomain::Extrinsic;

geodesic_draping::DrapeResult result =
    solver.solve(seedXY, fabricAngle, solveOptions, retrievalOptions);
```

After a solve, `solver.retrieve(retrievalOptions)` can retrieve the last core
result in another domain.

## Modes & Options
The solve behavior is configured with several option structs:\
`DrapeSolveMode`:
```cpp
enum class DrapeSolveMode {
  Fast,
  Hybrid,
  Complete,
};
```
`DrapeSolveOptions`:
```cpp
struct DrapeSolveOptions {
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  double fiberAngle = 90.0;
  AdvancedTraceOptions trace;
};
```
`fabricAngle` is the global fabric orientaion and passed to `solve()` / `solveDrape()`. `fiberAngle` is the
second fabric-family angle relative to the first, the default 90 degrees corresponds to an orthotropic fabric.

```cpp
enum class RetrievalDomain {
  Intrinsic,
  Extrinsic,
  Subdivision,
};
```
`RetrievalOptions`:
```cpp
struct RetrievalOptions {
  RetrievalDomain domain = RetrievalDomain::Extrinsic;
  bool sampleVertexShear = false;
};
```
`sampleVertexShear` optionally samples primary face shear to vertices.

## Intrinsic Backend And Refinement

```cpp
enum class IntrinsicTriangulationBackend {
  Signpost,
  IntegerCoordinates,
};

enum class RefinementMode {
  None,
  DelaunayFlip,
  DelaunayRefine,
};

struct RefinementOptions {
  RefinementMode mode = RefinementMode::None;
  std::optional<double> angleThreshold;
  std::optional<double> circumradiusThreshold;
  std::optional<size_t> maxInsertions;
};
```

Constructor with explicit intrinsic backend and refinement:

```cpp
geodesic_draping::RefinementOptions refinementOptions;
refinementOptions.mode = geodesic_draping::RefinementMode::DelaunayFlip;

geodesic_draping::GeoDrapeSolver solver(
    meshData,
    geodesic_draping::SignedHeatSolveOptions{},
    geodesic_draping::IntrinsicTriangulationBackend::IntegerCoordinates,
    refinementOptions);
```

## Intrinsic Input

If the seed and fabric direction already live on the active intrinsic mesh:

```cpp
geometrycentral::surface::SurfacePoint seed = ...;
geometrycentral::surface::BarycentricVector fabricDirection = ...;

geodesic_draping::DrapeResult result =
    solver.solveFromIntrinsic(seed, fabricDirection, solveOptions, retrievalOptions);
```

## Result

`DrapeResult` is a single result type with optional fields. A field is populated
only when the selected solve mode and retrieval domain can provide it.

Always present when successful:

```cpp
RetrievalDomain domain;
DrapeSolveMode mode;
const geometrycentral::surface::SurfaceMesh* mesh;
```

Geometry/domain handles:

```cpp
const geometrycentral::surface::IntrinsicGeometryInterface* intrinsicGeometry;
const geometrycentral::surface::VertexPositionGeometry* extrinsicGeometry;
std::optional<geometrycentral::surface::VertexData<geometrycentral::Vector3>>
    vertexPositions;
```

Seed, directions, and generator traces:

```cpp
std::optional<geometrycentral::surface::SurfacePoint> intrinsicSeed;
std::optional<std::array<geometrycentral::surface::BarycentricVector, 4>>
    intrinsicDirections;
std::optional<std::array<IntrinsicGeneratorTrace, 4>> intrinsicGenerators;

std::optional<geometrycentral::Vector3> extrinsicSeed;
std::optional<std::array<geometrycentral::Vector3, 4>> extrinsicDirections;
std::optional<std::array<ExtrinsicGeneratorTrace, 4>> extrinsicGenerators;
```

Fields:

```cpp
std::optional<std::array<geometrycentral::surface::FaceData<geometrycentral::Vector3>, 2>>
    directionFields;
std::optional<std::array<geometrycentral::surface::VertexData<double>, 2>>
    distances;
std::optional<geometrycentral::surface::FaceData<double>> faceShear;
std::optional<geometrycentral::surface::VertexData<double>> vertexShear;
```

For `solveDrape()`, `DrapeResult::storageOwner` keeps solver-owned geometry
alive. Keep the result object alive while using pointers or Geometry Central
data views stored in it.

## Benchmark

The repo also includes a very basic benchmark tool, which can be called like:

```powershell
.\build\Release\benchmark_drape.exe demo_part --warmup-runs 1 --measured-runs 10
```

This reports cold one-shot `solveDrape()` timings and warm repeated
`GeoDrapeSolver::solve()` timings for fast, hybrid, and complete modes.

To include subdivision retrieval cost:

```powershell
.\build\Release\benchmark_drape.exe demo_part --domain subdivision --warmup-runs 1 --measured-runs 10
```

## Notes

- The solver depends on Geometry Central.
- Polyscope is only needed for the debug viewer.
- The C++ API is Geometry-Central-native internally; the Python API converts at
  the package boundary to and from NumPy arrays.

## License

This project is licensed under the GNU General Public License v3.0 or later.
See `LICENSE`.

Geometry Central and Polyscope are MIT-licensed third-party dependencies.
