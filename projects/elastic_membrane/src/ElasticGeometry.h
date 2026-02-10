#pragma once

/**
 * @file ElasticGeometry.h
 * @brief Defines the ElasticGeometry class, extending VertexPositionGeometry with
 *        elastic membrane phenomenology (reference/actual metrics, curvatures,
 *        elastic energy, pressure, and gradient-based solving).
 */

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

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace geometrycentral {
namespace surface {

/**
 * @brief Elastic membrane geometry built on top of VertexPositionGeometry.
 *
 * Extends geometry-central's VertexPositionGeometry with elastic shell quantities:
 * reference and actual metrics/curvatures, material parameters (thickness, Young's
 * modulus, Poisson's ratio), elastic energy computation (stretching + bending),
 * pressure loading, and energy-gradient computation for optimization.
 *
 * Quantities follow the geometry-central dependency system: call requireXxx()
 * before accessing, and unrequireXxx() when no longer needed.
 */
class ElasticGeometry : public VertexPositionGeometry {

  public:

    /** @brief Construct with mesh only (no positions, all quantities zeroed). */
    ElasticGeometry(SurfaceMesh& mesh_);

    /** @brief Construct from vertex positions; reference state equals the initial configuration. */
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_);

    /** @brief Construct with uniform material parameters; reference state equals the initial configuration. */
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_, const double& THICKNESS_,
                    const double& YOUNGs_, const double& POISSONs_, const double& PRESSURE_);

    /** @brief Construct with uniform material parameters and explicitly specified reference lengths and angles. */
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                    const EdgeData<double>& L_bar_, const EdgeData<double>& B_bar_, const double& THICKNESS_,
                    const double& YOUNGs_, const double& POISSONs_, const double& PRESSURE_);

    /** @brief Fully specified constructor with per-face thickness, elastic tensor, and pressure. */
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                    const EdgeData<double>& L_bar_, const EdgeData<double>& B_bar_, const FaceData<double>& THICKNESS_,
                    const FaceData<Eigen::Matrix3f>& ElasticTensor_, const double PRESSURE_);

    virtual ~ElasticGeometry() {};

    // ========== Geometric Quantities (Reference State) ==========

    /** @brief Reference edge lengths (bar{l}_e), defining the stress-free metric. */
    EdgeData<double> referenceLengths;
    void requireReferenceLegths();
    void unrequireReferenceLegths();

    /** @brief Reference dihedral angles along edges, defining the stress-free curvature. */
    EdgeData<double> referenceEdgeDihedralAngles;
    void requireReferenceEdgeDihedralAngles();
    void unrequireReferenceEdgeDihedralAngles();

    /** @brief Reference metric per face, vectorized as (l1^2, l2^2, (l1^2+l2^2-l3^2)/2). */
    FaceData<Eigen::Vector3f> referenceMetric;
    void requireReferenceMetric();
    void unrequireReferenceMetric();

    /** @brief Reference curvature tensor per face, vectorized as (b11, b22, b12). */
    FaceData<Eigen::Vector3f> referenceCurvature;
    void requireReferenceCurvature();
    void unrequireReferenceCurvature();

    // ========== Geometric Quantities (Actual/Current State) ==========

    /** @brief Actual (current) metric per face, same vectorization as referenceMetric. */
    FaceData<Eigen::Vector3f> actualMetric;
    void requireActualMetric();
    void unrequireActualMetric();

    /** @brief Actual (current) curvature tensor per face, vectorized as (b11, b22, b12). */
    FaceData<Eigen::Vector3f> actualCurvature;
    void requireActualCurvature();
    void unrequireActualCurvature();

    /** @brief Actual shape operator per face, stored as (S11, S22, S12, S21). */
    FaceData<Eigen::Vector4f> actualShape;

    /** @brief Indices of the two base halfedges used to define the local frame per face. */
    FaceData<Eigen::Vector2i> baseEdges;

    // ========== Material Properties ==========

    /** @brief 3x3 elastic Cauchy tensor per face, for computing A*(g - bar{g}). */
    FaceData<Eigen::Matrix3f> elasticCauchyTensor;
    void requireElasticCauchyTensor();
    void unrequireElasticCauchyTensor();

    /** @brief Shell thickness per face. */
    FaceData<double> thickness;
    void requireThickness();
    void unrequireThickness();

    /** @brief Young's modulus per face. */
    FaceData<double> youngsModulus;
    void requireYoungsModulus();
    void unrequireYoungsModulus();

    /** @brief Poisson's ratio per face. */
    FaceData<double> poissonsRatio;
    void requirePoissonsRatio();
    void unrequirePoissonsRatio();

    /** @brief Uniform pressure applied normal to the surface. */
    double pressure;
    void requirePressure();
    void unrequirePressure();

    /** @brief Global coordinate scaling factor. */
    double coordinate_scale = 1;

    // ========== Energy Quantities ==========

    /** @brief Total elastic energy per face (stretching + bending). */
    FaceData<double> elasticEnergy;
    void requireElasticEnergy();
    void unrequireElasticEnergy();

    /** @brief Total energy per face (elastic + pressure work). */
    FaceData<double> totalEnergy;
    void requireTotalEnergy();
    void unrequireTotalEnergy();

    /** @brief Signed volume contribution per face (for pressure energy). */
    FaceData<double> faceVolume;
    void requireFaceVolume();
    void unrequireFaceVolume();

    /** @brief Stretching (in-plane) energy density per face. */
    FaceData<double> stretchingEnergy;
    void requireStretchingEnergy();
    void unrequireStretchingEnergy();

    /** @brief Bending (out-of-plane) energy density per face. */
    FaceData<double> bendingEnergy;
    void requireBendingEnergy();
    void unrequireBendingEnergy();

    // ========== Constraints ==========

    /** @brief Per-vertex region label for specifying regions with different properties. */
    VertexData<int> regions;
    void requireRegions();
    void unrequireRegions();

    /** @brief Per-vertex flag indicating whether the vertex position is fixed. */
    VertexData<bool> fixedVertexes;
    void requireFixedVertexes();
    void umrequireFixedVertexes();

    /** @brief Per-edge flag indicating whether the dihedral angle is fixed. */
    EdgeData<bool> fixedAngles;
    void requireFIxedAngles();
    void unrequireFIxedAngles();

    // ========== Solver Interface ==========

    /** @brief Per-vertex elastic energy gradient (force direction for minimization). */
    VertexData<Vector3> elasticGradient;

    /** @brief Computes the elastic energy gradient with respect to vertex positions. */
    void virtual computeGradient();

    /** @brief Sets reference dihedral angles from the current geometry. */
    void setReferenceAngles();

    /** @brief Computes a single descent step. */
    void computeStep();

    /** @brief Runs the full energy minimization solver. */
    void solve();

    /** @brief Evolves the geometry forward in time. */
    void evolve();

    /** @brief Returns the in-plane stress tensor per face. */
    FaceData<Eigen::Vector3f> getStress();

    /** @brief Returns the bending moment tensor per face. */
    FaceData<Eigen::Vector3f> getMoment();

    // ========== Curvature Queries ==========

    /** @brief Returns the mean curvature from the reference metric/curvature at face f. */
    double getReferenceMeanCurvautre(Face f);

    /** @brief Returns the Gaussian curvature from the reference metric/curvature at face f. */
    double getReferenceGaussianCurvautre(Face f);

    /** @brief Returns the mean curvature from the actual metric/curvature at face f. */
    double getActualMeanCurvautre(Face f);

    /** @brief Returns the Gaussian curvature from the actual metric/curvature at face f. */
    double getActualGaussianCurvautre(Face f);

    /** @brief Returns the mean curvature using actual curvature with reference metric at face f. */
    double getSemiActualMeanCurvautre(Face f);

    /** @brief Returns the Gaussian curvature using actual curvature with reference metric at face f. */
    double getSemiActualGaussianCurvautre(Face f);

    /** @brief Computes the actual shape operator for all faces from metric and curvature. */
    void getActualShape();

    /** @brief Computes the local 2D frame basis (two edge vectors) for all faces. */
    void getFaceBasis();

  protected:

    // ========== Dependent Quantity Handles (geometry-central dependency system) ==========

    /** @brief Dependency handle for referenceLengths. */
    DependentQuantityD<EdgeData<double>> referenceLengthsQ;
    virtual void computeReferenceLengths();

    /** @brief Dependency handle for referenceEdgeDihedralAngles. */
    DependentQuantityD<EdgeData<double>> referenceEdgeDihedralAnglesQ;
    virtual void computeReferenceEdgeDihedralAngles();

    /** @brief Dependency handle for referenceMetric. */
    DependentQuantityD<FaceData<Eigen::Vector3f>> referenceMetricQ;
    virtual void computeReferenceMetric();

    /** @brief Dependency handle for actualMetric. */
    DependentQuantityD<FaceData<Eigen::Vector3f>> actualMetricQ;
    virtual void computeActualMetric();

    /** @brief Dependency handle for referenceCurvature. */
    DependentQuantityD<FaceData<Eigen::Vector3f>> referenceCurvatureQ;
    virtual void computeReferenceCurvature();

    /** @brief Dependency handle for actualCurvature. */
    DependentQuantityD<FaceData<Eigen::Vector3f>> actualCurvatureQ;
    virtual void computeActualCurvature();

    /** @brief Dependency handle for elasticCauchyTensor. */
    DependentQuantityD<FaceData<Eigen::Matrix3f>> elasticCauchyTensorQ;
    virtual void computeElasticCauchyTensor();

    /** @brief Dependency handle for thickness. */
    DependentQuantityD<FaceData<double>> thicknessQ;
    virtual void computeThickness();

    /** @brief Dependency handle for youngsModulus. */
    DependentQuantityD<FaceData<double>> youngsModulusQ;
    virtual void computeYoungsModulus();

    /** @brief Dependency handle for poissonsRatio. */
    DependentQuantityD<FaceData<double>> poissonsRatioQ;
    virtual void computePoissonsRatio();

    /** @brief Dependency handle for elasticEnergy. */
    DependentQuantityD<FaceData<double>> elasticEnergyQ;
    virtual void computeElasticEnergy();

    /** @brief Dependency handle for totalEnergy. */
    DependentQuantityD<FaceData<double>> totalEnergyQ;
    virtual void computeTotalEnergy();

    /** @brief Dependency handle for faceVolume. */
    DependentQuantityD<FaceData<double>> faceVolumeQ;
    virtual void computeFaceVolume();

    /** @brief Dependency handle for stretchingEnergy. */
    DependentQuantityD<FaceData<double>> stretchingEnergyQ;
    virtual void computeStretchingEnergy();

    /** @brief Dependency handle for bendingEnergy. */
    DependentQuantityD<FaceData<double>> bendingEnergyQ;
    virtual void computeBendingEnergy();

    /** @brief Dependency handle for pressure. */
    DependentQuantityD<double> pressureQ;
    virtual void computePressure();

    /** @brief Dependency handle for regions. */
    DependentQuantityD<VertexData<int>> regionsQ;
    virtual void computeRegions();

    /** @brief Dependency handle for fixedVertexes. */
    DependentQuantityD<VertexData<bool>> fixedVertexesQ;
    virtual void computeFixedVertexs();

    /** @brief Dependency handle for fixedAngles. */
    DependentQuantityD<EdgeData<bool>> fixedAnglesQ;
    virtual void computeFixedAngles();

    // ========== Initialization Flags ==========

    bool isActualMetricInitializedF = false;
    bool isActualCurvatureInitializedF = false;
    bool isElasticTensorInitializedF = false;

    // ========== Per-Face Energy Computation Helpers ==========

    /** @brief Computes the stretching energy contribution for a single face. */
    void calculate_stretching_energy(const Face& f);

    /** @brief Computes the bending energy contribution for a single face. */
    void calculate_bending_energy(const Face& f);

    /** @brief Computes the total elastic energy for a single face. */
    void calculateFaceEnergy(const Face& f);

    /** @brief Computes the signed volume contribution for a single face. */
    void calculateFaceVolume(const Face& f);

    /** @brief Computes the total energy (elastic + pressure) for a single face. */
    void calculateFaceTotalEnergy(const Face& f);

    /** @brief Recomputes edge lengths for all edges adjacent to vertex v. */
    void calculate_adjacent_edges_lenght(const Vertex& v);

    /** @brief Recomputes face areas for all faces adjacent to vertex v. */
    void ElasticGeometry::calculate_adjucent_faces_area(const Vertex& v);

    /** @brief Recomputes face volumes for all faces adjacent to vertex v. */
    void calculate_adjacent_faces_volume(const Vertex& v);

    /** @brief Recomputes elastic energy for all faces adjacent to vertex v. */
    void calculate_adjacent_faces_energy(const Vertex& v);

    /** @brief Recomputes total energy for all faces adjacent to vertex v. */
    void calculate_adjacent_faces_total_energy(const Vertex& v);

  private:

    /** @brief Returns the curvature vector for face f (reference if _ref_or_act=0, actual if 1). */
    Vector3 get_curvature(Face& _f, const int& _ref_or_act);

    /** @brief Recomputes all local energy quantities around vertex v after a position change. */
    void updateLocalEnergy(const Vertex& v);

    /** @brief Recomputes the actual metric for all faces adjacent to vertex v. */
    void calculate_adjacent_faces_metric(const Vertex& v);

    /** @brief Recomputes the actual curvature for all faces adjacent to vertex v. */
    void calculate_adjacent_faces_curvature(const Vertex& v);

    /** @brief Recomputes dihedral angles for all edges adjacent to vertex v. */
    void calculate_adjacent_edges_dihedral_angles(const Vertex& v);

    /** @brief Computes the actual metric tensor for a single face. */
    void calculate_metric(const Face& f);

    /** @brief Computes the actual curvature tensor for a single face. */
    void calculate_curvature(const Face& f);

    /** @brief Computes the reference metric tensor for a single face from reference lengths. */
    void calculate_reference_metric(const Face& f);

    /** @brief Returns the local orthonormal frame basis vectors for face f. */
    std::vector<Vector3> ElasticGeometry::getFrameBasis(Face& f);

    /** @brief Computes the reference curvature for a single face from reference dihedral angles. */
    void calculate_reference_curvature(const Face& f);

    /** @brief Returns the mean curvature tr(a^{-1} b) / 2 from metric a and curvature b. */
    double getMean(Eigen::Vector3f a, Eigen::Vector3f b);

    /** @brief Returns the Gaussian curvature det(a^{-1} b) from metric a and curvature b. */
    double getDet(Eigen::Vector3f a, Eigen::Vector3f b);
};

} // namespace surface
} // namespace geometrycentral
