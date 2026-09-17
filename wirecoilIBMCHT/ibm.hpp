#ifndef IBM_HPP
#define IBM_HPP

// Forward declarations for everything defined (out of call order) in ibm.udf.
// Include this AFTER the geometry struct definitions (Point_3, Vertex, Vector,
// FacetNormal, Point, Triangle, IntegrationPoint, LagrangePoint, EmbreeMesh, ...)
// and BEFORE the first call site (IBM_Setup / the parse_stl/classify.../
// identifyCuttingTrianglesWithIntersections block), so every function below is
// visible to the compiler no matter where its body actually sits in the file.

#include <string>
#include <vector>
#include <array>
#include <utility>
#include <occa.hpp>

// --- geometry types (defined earlier in ibm.udf; declared here so the
//     prototypes below are self-contained if this header is ever moved) ---
struct Point_3;
struct Vertex;
struct Vector;
struct FacetNormal;
struct Point;
struct Triangle;
struct IntegrationPoint;
struct LagrangePoint;
struct EmbreeMesh;

// --- misc scalar/vector helpers ---
bool try_parse_float(const std::string &token, dfloat &value);

Vector vectorSubtract(const Vertex &a, const Vertex &b);
Vector crossProduct(const Vector &a, const Vector &b);
Vector normalize(const Vector &v);
dfloat Dot(const Vector &vector1, const Vector &vector2);
dfloat magnitude(const Vector &v);

dfloat computeTriangleArea(const Point_3 &v1, const Point_3 &v2, const Point_3 &v3);
dfloat computeSTLSurfaceArea(const std::vector<Vertex> &vertices);

// --- STL I/O ---
bool parse_stl(const std::string &filename,
                std::vector<Vertex> &vertices,
                std::vector<FacetNormal> &normals);
std::vector<Vertex> gatherSTLVertices();

// --- point / triangle / segment geometry ---
dfloat signed_distance_to_plane(const Point &point, const Vertex &v0, const Vertex &v1,
                                 const Vertex &v2, FacetNormal &normal);
Point projectPointOnPlane(const Point &point, const Vertex &v0, const Vertex &v1,
                           const Vertex &v2, FacetNormal &Normal);
bool pointInTriangle(const Point &p, const Vertex &v0, const Vertex &v1, const Vertex &v2);
dfloat pointToLineSegmentDistance(const Point &p, const Vertex &a, const Vertex &b);
dfloat signedDistanceToTriangle(const Point &p, const Vertex &v0, const Vertex &v1,
                                 const Vertex &v2, FacetNormal &Normal);
dfloat signedDistanceToTriangles(const std::vector<dfloat> &sdfVector);
dfloat intersectionSDF(const std::vector<dfloat> &sdfVector);

bool doLineSegmentsIntersect(const Point &p1, const Point &p2, const Point &q1, const Point &q2);
bool findLineSegmentIntersection(const Point &p1, const Point &p2, const Point &q1,
                                  const Point &q2, Point &intersectionPoint);
Point lineSegmentPlaneIntersection(const Point &p1, const Point &p2, const Vertex &v0,
                                    const FacetNormal &normal);

Point convertVertexToPoint(const Vertex &vertex);
std::vector<std::pair<Point, Point>> getElementSegments(const std::vector<Point> &quadPoints,
                                                          size_t elementIdx, int p);
bool arePointsSimilar(const Point &a, const Point &b, dfloat threshold = 1e-6);
dfloat distance(const Point &p1, const Point &p2);

// --- classification (Embree path) ---
EmbreeMesh buildEmbreeMeshFromSTL(const std::string &stlFilename);
std::vector<double> classifyPointsWithEMBREE(const std::vector<Point> &testPoints,
                                              const std::string &stlFilename);

// --- classification (VTK path) ---
std::vector<double> classifyPointsWithVTK(const std::vector<Point> &testPoints,
                                           const std::string &stlFilename);

// --- integration / Lagrange point generation ---
std::vector<IntegrationPoint> computeIntegrationPointsForTriangle(
    const Point_3 &v1, const Point_3 &v2, const Point_3 &v3);
std::vector<LagrangePoint> computeLagrangePointsForTriangle(
    const Point_3 &v1, const Point_3 &v2, const Point_3 &v3, int elementNumber);

// --- surface mesh output ---
void writeSurfaceMeshToVTK(const std::string &filename,
                            const std::vector<Point_3> &points,
                            const std::vector<std::array<std::size_t, 3>> &faces);

// --- cut-cell identification ---
void identifyCuttingTrianglesWithIntersections(
    const std::vector<Point> &quadPoints,
    int pointsPerElement,
    const std::vector<Vertex> &vertices,
    const std::vector<FacetNormal> &normals,
    const std::vector<dfloat> &sdfValues);

// --- IBM driver entry points (called from kgj.udf) ---
void IBM_Setup();
void IBM_PenaltyNoSlipBCCustomSource(double time);
void IBM_FillProp(occa::memory o_prop, double propSolid);

#endif // IBM_HPP

