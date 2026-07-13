# geodesic_draping

C++ implementation of geodesic draping on triangle meshes. The solver traces
two generator families, computes signed-heat direction fields, optionally
integrates distances, and returns shear data on the requested mesh domain.

The main public header is:

```cpp
#include "geodesic_draping/geodrape.h"
```

## Build

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 `
  -DGEODESIC_DRAPING_BUILD_TESTS=ON `
  -DGEODESIC_DRAPING_ENABLE_POLYSCOPE=ON

cmake --build build --config Release
ctest --test-dir build --build-config Release --output-on-failure
```

Polyscope is optional. Disable it if you only need the library and benchmark:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 `
  -DGEODESIC_DRAPING_ENABLE_POLYSCOPE=OFF
```

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

## Solve Modes

```cpp
enum class DrapeSolveMode {
  Fast,
  Hybrid,
  Complete,
};
```

- `Fast`: computes face shear from signed-heat direction fields. Does not
  return distance fields.
- `Hybrid`: computes face shear from signed-heat direction fields and also
  integrates vertex distance fields.
- `Complete`: integrates vertex distance fields and computes shear from the
  complete solve path.

`DrapeSolveOptions`:

```cpp
struct DrapeSolveOptions {
  DrapeSolveMode mode = DrapeSolveMode::Complete;
  double fiberAngle = 90.0;
  AdvancedSolveOptions advanced;
};
```

`fabricAngle` is passed to `solve()` / `solveDrape()`. `fiberAngle` is the
second fabric-family angle relative to the first.

## Retrieval Domains

```cpp
enum class RetrievalDomain {
  Intrinsic,
  Extrinsic,
  Subdivision,
};
```

- `Intrinsic`: returns data on the active intrinsic triangulation.
- `Extrinsic`: returns data transferred to the original input mesh when
  available.
- `Subdivision`: returns data on the common subdivision of input and intrinsic
  meshes. This is the best domain for debugging face fields after flips or
  refinement.

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

struct IntrinsicConstructionOptions {
  IntrinsicTriangulationBackend backend =
      IntrinsicTriangulationBackend::Signpost;
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

Constructor with explicit intrinsic options:

```cpp
geodesic_draping::IntrinsicConstructionOptions intrinsicOptions;
intrinsicOptions.backend =
    geodesic_draping::IntrinsicTriangulationBackend::IntegerCoordinates;

geodesic_draping::RefinementOptions refinementOptions;
refinementOptions.mode = geodesic_draping::RefinementMode::DelaunayFlip;

geodesic_draping::GeoDrapeSolver solver(
    meshData,
    geodesic_draping::SignedHeatSolveOptions{},
    intrinsicOptions,
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

## Debug Viewer

When Polyscope is enabled:

```powershell
.\build\Release\debug_drape_result.exe demo_part --mode complete --domain subdivision
```

Useful options:

```powershell
--mode fast|hybrid|complete
--domain extrinsic|subdivision
--backend signpost|integer
--refinement none|flip|refine
--sample-vertex-shear
```

The active fixture is `test_data/fixtures/demo_part`.

## Benchmark

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
- The current C++ API is native-first. Python bindings are not included yet.
