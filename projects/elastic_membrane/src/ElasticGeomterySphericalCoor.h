/**
 * @file ElasticGeomterySphericalCoor.h
 * @brief Spherical-coordinate specialization of the elastic shell geometry.
 *
 * Extends ElasticGeometry to work in spherical (theta, phi) coordinates
 * rather than Cartesian edge-based coordinates. Metrics and curvatures are
 * computed from coordinate differences on the sphere, which avoids the
 * dihedral-angle approach used by the base class and gives smoother
 * gradients for near-spherical shells.
 */
#pragma once
#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/surface_mesh.h"

#include <Eigen/SparseCore>

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "colormap.h"
#include <map>
#include <numeric>

#include "ElasticGeometry.h"
#include <fstream>
#include <limits>

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace geometrycentral {
namespace surface {

/**
 * @class ElasticGeometrySphericalCoor
 * @brief Elastic shell geometry using spherical coordinates for metric/curvature computation.
 *
 * Overrides the metric and curvature callbacks from ElasticGeometry so that
 * the first and second fundamental forms are evaluated from (theta, phi)
 * coordinate differences on the unit sphere, rather than from edge lengths
 * and dihedral angles. The gradient is also overridden to use central finite
 * differences in Cartesian space with local re-evaluation of the spherical
 * quantities around each perturbed vertex.
 */
class ElasticGeometrySphericalCoor : public ElasticGeometry {

  public:

    /** @brief Minimal constructor — delegates to the fully specified one with zero material constants. */
    ElasticGeometrySphericalCoor(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_);

    /**
     * @brief Fully specified constructor.
     *
     * Initializes the base ElasticGeometry, registers spherical-coordinate
     * dependency quantities (centroids, vertex coordinates, face centroid
     * coordinates), allocates all data arrays, and triggers the initial
     * computation of reference and actual metrics/curvatures.
     */
    ElasticGeometrySphericalCoor(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                                 const double& THICKNESS_, const double& YOUNGs_, const double& POISSONs_,
                                 const double& PRESSURE_);

    virtual ~ElasticGeometrySphericalCoor() {};

    // ========== Overridden Callbacks ==========

    /** @brief Computes reference metric from spherical coordinate differences and reference edge lengths. */
    void computeReferenceMetric() override;
    /** @brief Computes actual metric from spherical coordinate differences and current edge lengths. */
    void computeActualMetric() override;
    /** @brief Computes reference curvature from face normals and centroid coordinate differences. */
    void computeReferenceCurvature() override;
    /** @brief Computes actual curvature from face normals and centroid coordinate differences. */
    void computeActualCurvature() override;
    /** @brief Computes the elastic gradient via central finite differences with local re-evaluation. */
    void computeGradient() override;

    // ========== Spherical Coordinate Data ==========

    /** @brief Per-vertex spherical coordinates (theta, phi). */
    VertexData<Vector2> vertexCoordinates;
    /** @brief Registers a dependency on vertexCoordinates. */
    void requireVertexCoordinates();
    /** @brief Releases a dependency on vertexCoordinates. */
    void unrequireVertexCoordinates();

    /** @brief Per-face centroid spherical coordinates (theta, phi). */
    FaceData<Vector2> faceCentroidCoordinates;
    /** @brief Registers a dependency on faceCentroidCoordinates. */
    void requireFaceCentroidCoordinates();
    /** @brief Releases a dependency on faceCentroidCoordinates. */
    void unrequireFaceCentroidCoordinates();

    /** @brief Per-face centroid 3D position. */
    FaceData<Vector3> faceCentroidPosition;
    /** @brief Registers a dependency on faceCentroidPosition. */
    void requireFaceCentroidPosition();
    /** @brief Releases a dependency on faceCentroidPosition. */
    void unrequireFaceCentroidPosition();

    /** @brief Forces re-computation of the elastic Cauchy tensor by clearing and re-requiring it. */
    void updateElasticCauchyTensor();
    /** @brief Overwrites face centroid coordinates from external data and recomputes metrics/curvatures. */
    void updateFaceCentroidCoordinates(const FaceData<Vector2> faceCoordinates);

  protected:

    // ========== Dependency Handles ==========

    /** @brief Dependency handle for faceCentroidPosition. */
    DependentQuantityD<FaceData<Vector3>> faceCentroidPositionQ;
    /** @brief Computes face centroid positions by averaging vertex positions. */
    virtual void computeCentroids();

    /** @brief Dependency handle for vertexCoordinates. */
    DependentQuantityD<VertexData<Vector2>> vertexCoordinatesQ;
    /** @brief Converts vertex 3D positions to spherical (theta, phi) coordinates. */
    virtual void computeVertexCoordinates();

    /** @brief Dependency handle for faceCentroidCoordinates. */
    DependentQuantityD<FaceData<Vector2>> faceCentroidCoordinatesQ;
    /** @brief Converts face centroid 3D positions to spherical (theta, phi) coordinates. */
    virtual void computeFaceCentroidCoordinates();

  private:

    // ========== Per-Face Metric/Curvature Helpers ==========

    /** @brief Computes the actual metric tensor for a single face from current edge lengths. */
    void calculateFaceActualMetric(Face f);
    /** @brief Computes the reference metric tensor for a single face from reference edge lengths. */
    void calculateFaceReferenceMetric(Face f);
    /** @brief Computes the actual curvature tensor for a single face from normals and centroids. */
    void calcualteFaceActualCurvature(Face f);
    /** @brief Computes the reference curvature tensor for a single face from normals and centroids. */
    void calcualteFaceReferenceCurvature(Face f);

    // ========== Local Update Helpers ==========

    /** @brief Re-evaluates energy for all faces adjacent to vertex v after a perturbation. */
    void localEnergyChange(const Vertex v);
    /** @brief Recomputes normals and centroid positions for faces adjacent to vertex v. */
    void calculateAdjacentNormalAndCentroids(const Vertex v);
    /** @brief Recomputes the actual metric for faces adjacent to vertex v. */
    void calculateAdjacentMetric(const Vertex v);
    /** @brief Recomputes the actual curvature for faces adjacent to vertex v. */
    void calculateAdjacentCurvature(const Vertex v);

    // ========== Core Computation ==========

    /** @brief Computes the curvature tensor for face f from centroid coordinate differences and face normals. */
    Eigen::Vector3f calculateCurvature(Face f);
    /** @brief Computes the metric tensor for face f from squared edge lengths and coordinate differences. */
    Eigen::Vector3f calculateMetric(Face f, float lengths[3]);
    /** @brief Computes the shortest-path coordinate difference between two spherical points, handling periodicity. */
    Vector2 calculateCoordinateDiff(Vector2 p1, Vector2 p2);

    /** @brief True after the reference metric has been initialized (prevents re-computation). */
    bool isReferenceMetricInitializedF = false;
    /** @brief True after the reference curvature has been initialized (prevents re-computation). */
    bool isReferenceCurvatureInitializedF = false;
};

} // namespace surface
} // namespace geometrycentral