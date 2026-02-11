/**
 * @file ElasticGeometry.cpp
 * @brief Implementation of the ElasticGeometry class — elastic membrane
 *        geometry with reference/actual metrics, curvatures, energies,
 *        and gradient computation.
 */

#include "ElasticGeometry.h"
#include <fstream>
#include <limits>


using namespace geometrycentral;
using namespace geometrycentral::surface;


namespace geometrycentral {
namespace surface {


// clang-format off


/**
 * @brief Principal constructor: fully specified elastic geometry.
 *
 * Initializes all dependency-tracked quantities, binds compute callbacks,
 * stores the supplied reference lengths, dihedral angles, thickness, elastic
 * tensor, and pressure. Computes reference and actual metrics/curvatures
 * when vertex positions are non-zero.
 *
 * @param mesh_                  Surface mesh.
 * @param inputVertexPositions_  Initial vertex positions.
 * @param L_bar_                 Reference edge lengths.
 * @param B_bar_                 Reference dihedral angles.
 * @param THICKNESS_             Per-face shell thickness.
 * @param ElasticTensor_         Per-face 3x3 elastic Cauchy tensor.
 * @param PRESSURE_              Uniform pressure.
 */
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_,const VertexData<Vector3>& inputVertexPositions_, const EdgeData<double>& L_bar_,
                    const EdgeData<double>& B_bar_,const FaceData<double>& THICKNESS_,
                    const FaceData<Eigen::Matrix3f>& ElasticTensor_, const double PRESSURE_) : 
    VertexPositionGeometry(mesh_),

    referenceLengthsQ                   (&referenceLengths,                 std::bind(&ElasticGeometry::computeReferenceLengths, this),                         quantities),
    referenceEdgeDihedralAnglesQ        (&referenceEdgeDihedralAngles,      std::bind(&ElasticGeometry::computeReferenceEdgeDihedralAngles,this),               quantities),
    referenceMetricQ                    (&referenceMetric,                  std::bind(&ElasticGeometry::computeReferenceMetric, this),                          quantities),
    actualMetricQ                       (&actualMetric,                     std::bind(&ElasticGeometry::computeActualMetric, this),                             quantities),
    referenceCurvatureQ                 (&referenceCurvature,               std::bind(&ElasticGeometry::computeReferenceCurvature, this),                       quantities),
    actualCurvatureQ                    (&actualCurvature,                  std::bind(&ElasticGeometry::computeActualCurvature, this),                          quantities),
    elasticCauchyTensorQ                (&elasticCauchyTensor,              std::bind(&ElasticGeometry::computeElasticCauchyTensor, this),                      quantities),
    thicknessQ                          (&thickness,                        std::bind(&ElasticGeometry::computeThickness, this),                                quantities),
    youngsModulusQ                      (&youngsModulus,                    std::bind(&ElasticGeometry::computeYoungsModulus, this),                            quantities),
    poissonsRatioQ                      (&poissonsRatio,                    std::bind(&ElasticGeometry::computePoissonsRatio, this),                            quantities),
    elasticEnergyQ                      (&elasticEnergy,                    std::bind(&ElasticGeometry::computeElasticEnergy, this),                            quantities),
    totalEnergyQ                        (&totalEnergy,                      std::bind(&ElasticGeometry::computeTotalEnergy, this),                              quantities),
    faceVolumeQ                         (&faceVolume,                       std::bind(&ElasticGeometry::computeFaceVolume, this),                               quantities),
    stretchingEnergyQ                   (&stretchingEnergy,                 std::bind(&ElasticGeometry::computeStretchingEnergy, this),                         quantities),
    bendingEnergyQ                      (&bendingEnergy,                    std::bind(&ElasticGeometry::computeBendingEnergy, this),                            quantities),
    pressureQ                           (&pressure,                         std::bind(&ElasticGeometry::computePressure,this),                                  quantities),
    regionsQ                            (&regions,                          std::bind(&ElasticGeometry::computeRegions, this),                                  quantities),
    fixedVertexesQ                      (&fixedVertexes,                    std::bind(&ElasticGeometry::computeFixedVertexs,this),                              quantities),
    fixedAnglesQ                        (&fixedAngles,                      std::bind(&ElasticGeometry::computeFixedAngles,this),                               quantities)    
    {
         vertexPositions = inputVertexPositions_;
         requireVertexPositions();
         vertexPositionsQ.clearable = false;

         referenceLengths=L_bar_;
         requireReferenceLegths();
         referenceLengthsQ.clearable = false;
        
         referenceEdgeDihedralAngles=B_bar_;
         requireReferenceEdgeDihedralAngles();
         referenceEdgeDihedralAnglesQ.clearable=false;

         thickness=THICKNESS_;
         requireThickness();
         thicknessQ.clearable=false;

         youngsModulus= FaceData<double>(mesh_, 0.);
         requireYoungsModulus();

         poissonsRatio= FaceData<double>(mesh_, 0.);
         requirePoissonsRatio();

         coordinate_scale = 1;

         elasticCauchyTensor=ElasticTensor_;
         requireElasticCauchyTensor();

         pressure=PRESSURE_;
         requirePressure();
         pressureQ.clearable=false;

         referenceMetric = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(1., 2., 3.));
         requireReferenceMetric();

         referenceCurvature = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(1., 2., 3.));
         requireReferenceCurvature();

         actualMetric = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(0., 0., 0.));
         actualCurvature = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(0., 0., 0.));

         double xsum=0;
         double ysum=0;
         double zsum=0;
         for (Vertex v : this->mesh.vertices()) {
           xsum += inputVertexPositions[v].x;
           ysum += inputVertexPositions[v].y;
           zsum += inputVertexPositions[v].z;
         }
         if(xsum+ysum+zsum!=0){
            requireActualMetric();
            requireActualCurvature();
         }
          
         
         stretchingEnergy = FaceData<double>(mesh_, 0.);
         bendingEnergy = FaceData<double>(mesh_, 0.);
         elasticEnergy = FaceData<double>(mesh_, 0.);
         requireStretchingEnergy();
         requireBendingEnergy();
         requireElasticEnergy();

      }

// clang-format on


/** @brief Mesh-only constructor: all quantities zeroed, no positions. */
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_)
    : ElasticGeometry::ElasticGeometry(mesh_, VertexData<Vector3>(mesh_, Vector3{0., 0., 0.}),
                                       EdgeData<double>(mesh_, 0.), EdgeData<double>(mesh_, 0.),
                                       FaceData<double>(mesh_, 0.), FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()),
                                       0.) {}

/** @brief Positions-only constructor: reference state set from current geometry, no material params. */
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_)
    : ElasticGeometry::ElasticGeometry(mesh_, inputVertexPositions_, EdgeData<double>(mesh_, 0.),
                                       EdgeData<double>(mesh_, 0.), FaceData<double>(mesh_, 0.),
                                       FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()), 0.) {
    this->requireReferenceLegths();
    this->requireReferenceEdgeDihedralAngles();
}


/**
 * @brief Uniform isotropic constructor: builds elastic tensor from Young's modulus
 *        and Poisson's ratio, reference state from current geometry.
 */
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                                 const double& THICKNESS_, const double& YOUNGs_, const double& POISSONs_,
                                 const double& PRESSURE_)
    : ElasticGeometry::ElasticGeometry(mesh_, inputVertexPositions_, EdgeData<double>(mesh_, 0.),
                                       EdgeData<double>(mesh_, 0.), FaceData<double>(mesh_, THICKNESS_),
                                       FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()), PRESSURE_) {
    youngsModulus = FaceData<double>(mesh_, YOUNGs_);
    requireYoungsModulus();
    youngsModulusQ.clearable = false;

    poissonsRatio = FaceData<double>(mesh_, POISSONs_);
    requirePoissonsRatio();
    poissonsRatioQ.clearable = false;

    unrequireElasticCauchyTensor();
    elasticCauchyTensorQ.clearIfNotRequired();
    elasticCauchyTensor = FaceData<Eigen::Matrix3f>(this->mesh, Eigen::Matrix3f());
    requireElasticCauchyTensor();
    elasticCauchyTensorQ.clearable = false;

    requireActualMetric();
    requireActualCurvature();

    unrequireStretchingEnergy();
    unrequireBendingEnergy();
    unrequireElasticEnergy();
    stretchingEnergyQ.clearIfNotRequired();
    bendingEnergyQ.clearIfNotRequired();
    elasticEnergyQ.clearIfNotRequired();
    stretchingEnergy = FaceData<double>(mesh_, 0.);
    bendingEnergy = FaceData<double>(mesh_, 0.);
    elasticEnergy = FaceData<double>(mesh_, 0.);
    requireStretchingEnergy();
    requireBendingEnergy();
    requireElasticEnergy();
    requireEdgeDihedralAngles();
    edgeDihedralAnglesQ.clearable = false;
}


// ========== Require/Unrequire Pairs ==========
// Each pair delegates to the DependentQuantityD handle, triggering lazy
// computation on first require and reference-counting for cleanup.

void ElasticGeometry::requireReferenceLegths() {
    referenceLengthsQ.require();
}
void ElasticGeometry::unrequireReferenceLegths() {
    referenceLengthsQ.unrequire();
}

void ElasticGeometry::requireReferenceEdgeDihedralAngles() {
    referenceEdgeDihedralAnglesQ.require();
}
void ElasticGeometry::unrequireReferenceEdgeDihedralAngles() {
    referenceEdgeDihedralAnglesQ.unrequire();
}

void ElasticGeometry::requireReferenceMetric() {
    referenceMetricQ.require();
}
void ElasticGeometry::unrequireReferenceMetric() {
    referenceMetricQ.unrequire();
}

void ElasticGeometry::requireActualMetric() {
    actualMetricQ.require();
}
void ElasticGeometry::unrequireActualMetric() {
    actualMetricQ.unrequire();
}

void ElasticGeometry::requireReferenceCurvature() {
    referenceCurvatureQ.require();
}
void ElasticGeometry::unrequireReferenceCurvature() {
    referenceCurvatureQ.unrequire();
}


void ElasticGeometry::requireActualCurvature() {
    actualCurvatureQ.require();
}
void ElasticGeometry::unrequireActualCurvature() {
    actualCurvatureQ.unrequire();
}


void ElasticGeometry::requireElasticCauchyTensor() {
    elasticCauchyTensorQ.require();
}
void ElasticGeometry::unrequireElasticCauchyTensor() {
    elasticCauchyTensorQ.unrequire();
}


void ElasticGeometry::requireThickness() {
    thicknessQ.require();
}
void ElasticGeometry::unrequireThickness() {
    thicknessQ.unrequire();
}


void ElasticGeometry::requireYoungsModulus() {
    youngsModulusQ.require();
}
void ElasticGeometry::unrequireYoungsModulus() {
    youngsModulusQ.unrequire();
}


void ElasticGeometry::requirePoissonsRatio() {
    poissonsRatioQ.require();
}
void ElasticGeometry::unrequirePoissonsRatio() {
    poissonsRatioQ.unrequire();
}


void ElasticGeometry::requireElasticEnergy() {
    elasticEnergyQ.require();
}
void ElasticGeometry::unrequireElasticEnergy() {
    elasticEnergyQ.unrequire();
}


void ElasticGeometry::requireTotalEnergy() {
    totalEnergyQ.require();
}
void ElasticGeometry::unrequireTotalEnergy() {
    totalEnergyQ.unrequire();
}


void ElasticGeometry::requireFaceVolume() {
    faceVolumeQ.require();
}
void ElasticGeometry::unrequireFaceVolume() {
    faceVolumeQ.unrequire();
}

void ElasticGeometry::requireStretchingEnergy() {
    stretchingEnergyQ.require();
}
void ElasticGeometry::unrequireStretchingEnergy() {
    stretchingEnergyQ.unrequire();
}

void ElasticGeometry::requireBendingEnergy() {
    bendingEnergyQ.require();
}
void ElasticGeometry::unrequireBendingEnergy() {
    bendingEnergyQ.unrequire();
}


void ElasticGeometry::requirePressure() {
    pressureQ.require();
}
void ElasticGeometry::unrequirePressure() {
    pressureQ.unrequire();
}


void ElasticGeometry::requireRegions() {
    regionsQ.require();
}
void ElasticGeometry::unrequireRegions() {
    regionsQ.unrequire();
}


void ElasticGeometry::requireFixedVertexes() {
    fixedVertexesQ.require();
}
void ElasticGeometry::umrequireFixedVertexes() {
    fixedVertexesQ.unrequire();
}


void ElasticGeometry::requireFIxedAngles() {
    fixedAnglesQ.require();
}
void ElasticGeometry::unrequireFIxedAngles() {
    fixedAnglesQ.unrequire();
}




// ========== Helper Functions ==========

/**
 * @brief Checks whether any element in a MeshData container is exactly zero.
 * @tparam data_type MeshData type.
 * @param data The mesh data to check.
 * @return true if at least one element is zero.
 */
template <typename data_type>
static bool is_illegal(data_type& data) {
    bool any_zeros = false;
    for (int index = 0; index < data.size(); index++) {
        any_zeros += data[index] == 0.0;
    }
    return any_zeros;
}



/**
 * @brief Computes the dot product of two Vector3 values.
 * @param vec1 First vector.
 * @param vec2 Second vector.
 * @return The dot product vec1 . vec2.
 */
double project3(const Vector3& vec1, const Vector3& vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

/**
 * @brief Builds a local 2D frame from the first two edge directions of a face.
 * @param f The face whose frame basis is computed.
 * @return A vector of two normalized edge direction vectors.
 */
std::vector<Vector3> ElasticGeometry::getFrameBasis(Face& f) {
    int baseIndex = 0;
    Vector3 base_vectors[2];
    Vector3 _temp_edge_vec;
    for (Edge e : f.adjacentEdges()) { // As basis we use the same basis
        if (baseIndex < 2) {
            _temp_edge_vec = this->vertexPositions[e.firstVertex()] - this->vertexPositions[e.secondVertex()];
            _temp_edge_vec = _temp_edge_vec.normalize();
            base_vectors[baseIndex] = _temp_edge_vec;
        }
        baseIndex++;
    }
    

    std::vector<Vector3> frame = {base_vectors[0], base_vectors[1]};
    return frame;
}



/**
 * @brief Computes the curvature tensor components (b11, b22, b12) for a face.
 *
 * Uses a fixed reference coordinate system on the triangle (vertices at
 * {0,0}, {1,0}, {0,1}) to compute edge normals, centroid-to-midedge vectors,
 * and curvature from dihedral angles. The curvature is measured relative to
 * the reference metric even for the actual curvature case.
 *
 * @param f           Face to compute curvature for.
 * @param ref_or_act  0 = reference curvature (from reference angles),
 *                    1 = actual curvature (from current dihedral angles).
 * @return Vector3 of curvature components (b11, b22, b12).
 */
Vector3 ElasticGeometry::get_curvature(Face& f, const int& ref_or_act) {
    EdgeData<double> angles;
    Eigen::Vector3f metric;
    if (ref_or_act == 0) {
        angles = this->referenceEdgeDihedralAngles;
        metric = this->referenceMetric[f];
    } else {
        this->edgeDihedralAnglesQ.ensureHave();
        angles = this->edgeDihedralAngles;
        metric = this->referenceMetric[f];
    }

    double totLength = 0.;
    double curve_comp1 = 0.; // in the 1 direction
    double curve_comp2 = 0.; // in the 2 direction
    double curve_comp3 = 0.;
    double curvatureMagnitude;
    double tmpEdgeLength;
    double centerToMidLength;
    Vector2 edgesCoordinates[3] = {{1. * coordinate_scale, 0.},
                                   {-1. * coordinate_scale, 1. * coordinate_scale},
                                   {0., -1. * coordinate_scale}}; // edges in coordinates
    Vector2 centroidMidEdgeVec[3] = {{.1666667 * coordinate_scale, -.3333333 * coordinate_scale},
                                     {.1666667 * coordinate_scale, .1666667 * coordinate_scale},
                                     {-.3333333 * coordinate_scale, .1666667 * coordinate_scale}};

    Vector2 edgeNormal;
    double factor =1/std::sqrt(metric[0] * metric[1] - metric[2] * metric[2]);
    int edgeindex = 0;
    Vector3 angs = {0., 0., 0.};
    double CML2[3] = {0., 0., 0.};
    Vector3 edgeNormals[3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
    
    for (Edge e : f.adjacentEdges()) {
        totLength += referenceLengths[e];
        tmpEdgeLength = referenceLengths[e];
        edgeNormals[edgeindex] = {
            factor / tmpEdgeLength *
                          (metric[2] * edgesCoordinates[edgeindex][0] + metric[1] * edgesCoordinates[edgeindex][1]),
                      -factor / tmpEdgeLength *
                          (metric[0] * edgesCoordinates[edgeindex][0] + metric[2] * edgesCoordinates[edgeindex][1])};
        centerToMidLength =
            centroidMidEdgeVec[edgeindex][0] * edgeNormal[0] * metric[0] +
            centroidMidEdgeVec[edgeindex][1] * edgeNormal[1] * metric[1] +
            (centroidMidEdgeVec[edgeindex][0] * edgeNormal[1] + centroidMidEdgeVec[edgeindex][1] * edgeNormal[0]) *
                metric[2];
        CML2[edgeindex] = std::sqrt(metric[0] * centroidMidEdgeVec[edgeindex][0] * centroidMidEdgeVec[edgeindex][0] +
                          metric[1] * centroidMidEdgeVec[edgeindex][1] * centroidMidEdgeVec[edgeindex][1] +
                          2 * metric[2] * centroidMidEdgeVec[edgeindex][0] * centroidMidEdgeVec[edgeindex][1]);
        curvatureMagnitude = 0.5 * angles[e];
        angs[edgeindex] = curvatureMagnitude;
        curve_comp1 += tmpEdgeLength * curvatureMagnitude / centerToMidLength * 
                       (metric[0] * edgeNormal[0] + metric[2] * edgeNormal[1]) *
                       (metric[0] * edgeNormal[0] + metric[2] * edgeNormal[1]);
        curve_comp2 += tmpEdgeLength * curvatureMagnitude / centerToMidLength *
                       (metric[2] * edgeNormal[0] + metric[1] * edgeNormal[1]) *
                       (metric[2] * edgeNormal[0] + metric[1] * edgeNormal[1]);
        curve_comp3 += tmpEdgeLength * curvatureMagnitude / centerToMidLength *
                       (metric[0] * edgeNormal[0] + metric[2] * edgeNormal[1]) *
                       (metric[2] * edgeNormal[0] + metric[1] * edgeNormal[1]);
        edgeindex++;
    }

   Vector3 res{curve_comp1 / totLength, curve_comp2 / totLength, curve_comp3 / totLength};
   res = {(-4. * angs[0] * CML2[0] + 8. * angs[1] * CML2[1] + 8. * angs[2] * CML2[2]) / std::pow(coordinate_scale, 2.),
          (8. * angs[0] * CML2[0] + 8. * angs[2] * CML2[2] - 4. * angs[1] * CML2[1]) / std::pow(coordinate_scale, 2.),
          (-2. * angs[0] * CML2[0] + 10. * angs[2] * CML2[2] - 2. * angs[1] * CML2[1]) /
              std::pow(coordinate_scale, 2.)};

    return res;
}

// ========== Compute Callbacks ==========

/**
 * @brief Resets the reference curvature from current reference angles and
 *        recomputes the elastic energy.
 */
void ElasticGeometry::setReferenceAngles() {
    referenceMetricQ.ensureHave();
    elasticEnergyQ.unrequire();
    referenceCurvatureQ.unrequire();
    elasticEnergyQ.clearIfNotRequired();
    for (Edge e : mesh.edges()) {
        referenceEdgeDihedralAngles[e] *= 1.;
    }
    refreshQuantities();
    computeElasticEnergy();    
    elasticEnergyQ.require();
   

}


/**
 * @brief Computes reference edge lengths from current geometry if they are all zero.
 *
 * Only runs once at initialization when no reference lengths were supplied;
 * copies the actual edge lengths as the stress-free reference.
 */
void ElasticGeometry::computeReferenceLengths() {
    if (referenceLengths.toVector().isZero()) {
        this->requireEdgeLengths();
        for (Edge e : this->mesh.edges()) {
            referenceLengths[e] = this->edgeLength(e);
        }
    }
}


/**
 * @brief Computes reference dihedral angles from current geometry if they are all zero.
 *
 * Same lazy-initialization logic as computeReferenceLengths.
 */
void ElasticGeometry::computeReferenceEdgeDihedralAngles() {
    if (referenceEdgeDihedralAngles.toVector().isZero()) {
        this->requireEdgeDihedralAngles();
        for (Edge e : this->mesh.edges()) {
            referenceEdgeDihedralAngles[e] = this->edgeDihedralAngle(e);
        }
    }
}


/** @brief Computes the reference metric for every face from reference edge lengths. */
void ElasticGeometry::computeReferenceMetric() {
    for (Face f : this->mesh.faces()) {
        calculate_reference_metric(f);
    }
}

/** @brief Computes the actual metric for every face from current edge lengths. */
void ElasticGeometry::computeActualMetric() {
    for (Face f : this->mesh.faces()) {
        calculate_metric(f);
    }
}

/** @brief Computes the reference curvature tensor for every face from reference dihedral angles. */
void ElasticGeometry::computeReferenceCurvature() {
    this->faceNormalsQ.ensureHave();

    for (Face f : this->mesh.faces()) {
        Vector3 _curve_comp = get_curvature(f, 0);
        this->referenceCurvature[f][0] = _curve_comp[0];
        this->referenceCurvature[f][1] = _curve_comp[1];
        this->referenceCurvature[f][2] = _curve_comp[2];
    }
} 

/** @brief Computes the actual curvature tensor for every face from current dihedral angles. */
void ElasticGeometry::computeActualCurvature() {
    this->faceNormalsQ.ensureHave();
    for (Face f : this->mesh.faces()) {
        Vector3 _curve_comp = get_curvature(f, 1);
        this->actualCurvature[f][0] = _curve_comp[0];
        this->actualCurvature[f][1] = _curve_comp[1];
        this->actualCurvature[f][2] = _curve_comp[2];
    }
}

/**
 * @brief Builds the 3x3 elastic Cauchy tensor per face from Young's modulus,
 *        Poisson's ratio, and the inverse reference metric.
 *
 * Only runs once (guarded by isElasticTensorInitializedF) and only if
 * Young's modulus is non-zero. The tensor encodes the isotropic linear
 * elastic response A such that stress ~ A * (g - bar{g}).
 */
void ElasticGeometry::computeElasticCauchyTensor() {
    if (!isElasticTensorInitializedF && !youngsModulus.toVector().isZero()) {
        double _invmet[3];
        double _det;
        double _coef;
        for (Face f : this->mesh.faces()) {
            _det = referenceMetric[f][0] * referenceMetric[f][1] - referenceMetric[f][2] * referenceMetric[f][2];
            _invmet[0] =  referenceMetric[f][1] / _det;
            _invmet[1] =  referenceMetric[f][0] / _det;
            _invmet[2] = -referenceMetric[f][2] / _det;
            _coef = youngsModulus[f] / (1 - poissonsRatio[f] * poissonsRatio[f]) / 8.;
            elasticCauchyTensor[f](0, 0) = _coef * _invmet[0] * _invmet[0];
            elasticCauchyTensor[f](0, 1) =
                _coef * (_invmet[2] * _invmet[2] * (1 - poissonsRatio[f]) + _invmet[0] * _invmet[1] * poissonsRatio[f]);
            elasticCauchyTensor[f](1, 0) = elasticCauchyTensor[f](0, 1);
            elasticCauchyTensor[f](0, 2) = _coef * 2 * _invmet[0] * _invmet[2];
            elasticCauchyTensor[f](2, 0) = .5* elasticCauchyTensor[f](0, 2);
            elasticCauchyTensor[f](1, 1) = _coef * _invmet[1] * _invmet[1];
            elasticCauchyTensor[f](1, 2) = _coef * 2* _invmet[1] * _invmet[2];
            elasticCauchyTensor[f](2, 1) = 0.5 * elasticCauchyTensor[f](1, 2);
            elasticCauchyTensor[f](2, 2) = _coef * (_invmet[2] * _invmet[2] * (1 + poissonsRatio[f]) +
                                                    _invmet[0] * _invmet[1] * (1 - poissonsRatio[f]));
        }
        isElasticTensorInitializedF = true;
    }
}

/** @brief No-op: thickness is set externally and does not need recomputation. */
void ElasticGeometry::computeThickness() {}

/** @brief No-op: Young's modulus is set externally. */
void ElasticGeometry::computeYoungsModulus() {}

/** @brief No-op: Poisson's ratio is set externally. */
void ElasticGeometry::computePoissonsRatio() {}

/**
 * @brief Computes total elastic energy (stretching + bending) for every face.
 *
 * Requires bending and stretching energies, then delegates per-face
 * computation to calculateFaceEnergy.
 */
void  ElasticGeometry::computeElasticEnergy() {
    this->requireBendingEnergy();
    this->requireStretchingEnergy();
    if (elasticEnergy.size() == 0) elasticEnergy = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces())
    {
        calculateFaceEnergy(f);
    }
}


/**
 * @brief Computes total energy (elastic - pressure * volume) for every face.
 */
void ElasticGeometry::computeTotalEnergy() {
    this->requireElasticEnergy();
    this->requireFaceVolume();
    if (totalEnergy.size() == 0) totalEnergy = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces()) {
        calculateFaceTotalEnergy(f);
    }    
}


/**
 * @brief Computes the signed volume contribution for every face (assumes centered object).
 */
void ElasticGeometry::computeFaceVolume() {
    this->requireFaceNormals();
    if (faceVolume.size() == 0) faceVolume = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces()) {
        calculateFaceVolume(f);                
       
    }
}

/** @brief Computes stretching (in-plane) energy for every face. */
void  ElasticGeometry::computeStretchingEnergy() {
    this->requireFaceAreas();
    if (stretchingEnergy.size() == 0) stretchingEnergy = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces()) {
        calculate_stretching_energy(f);
    }
}

/**
 * @brief Computes bending (out-of-plane) energy for every face.
 *
 * Ensures all prerequisite quantities (metrics, Cauchy tensor, curvatures)
 * are available, forces an actual curvature recomputation, then delegates
 * per-face computation to calculate_bending_energy.
 */
void ElasticGeometry::computeBendingEnergy() {
    actualMetricQ.ensureHave();
    referenceMetricQ.ensureHave();
    elasticCauchyTensorQ.ensureHave();
    referenceCurvatureQ.ensureHave();
    actualCurvatureQ.ensureHave();
    computeActualCurvature();

    if (bendingEnergy.size() == 0) bendingEnergy = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces()) {
        calculate_bending_energy(f);
    }
}



/**
 * @brief Computes the total energy gradient via central finite differences.
 *
 * For each vertex and each spatial direction (x, y, z), perturbs the vertex
 * position by +/- epsilon, recomputes local energy, and estimates the
 * partial derivative as -(E+ - E-) / (2*epsilon). The result is stored
 * in elasticGradient (pointing downhill).
 */
void ElasticGeometry::computeGradient() {
    vertexDualAreasQ.ensureHave();
    vertexNormalsQ.ensureHave();
    totalEnergyQ.ensureHave();
    elasticEnergyQ.ensureHave();
    elasticGradient = VertexData<Vector3>(this->mesh, Vector3{0., 0., 0.});
    double _epsilon = 1.e-6;    
    for (Vertex v : this->mesh.vertices()) {
        for (int _direction = 0; _direction < 3; _direction++) {
            double _ePlus=0.;
            double _eMinus=0.;

            vertexPositions[v][_direction] += _epsilon;
            updateLocalEnergy(v);
            for (Face f : v.adjacentFaces()) {
                _ePlus += this->totalEnergy[f];
            }

            vertexPositions[v][_direction] -= 2.*_epsilon;
            updateLocalEnergy(v);
            for (Face f : v.adjacentFaces()) {
                _eMinus += this->totalEnergy[f];
            }

            vertexPositions[v][_direction] += _epsilon;
            updateLocalEnergy(v);

            elasticGradient[v][_direction] += -(_ePlus - _eMinus) / 2. / _epsilon;
        }
    }
}

/**
 * @brief Recomputes all derived quantities in the 1-ring around vertex v.
 *
 * Called after a vertex position perturbation during gradient computation.
 * Updates edge lengths, face areas, dihedral angles, face normals,
 * metrics, curvatures, energies, volumes, and total energies locally.
 *
 * @param v The vertex whose neighborhood needs updating.
 */
void ElasticGeometry::updateLocalEnergy(const Vertex& v) {
    faceVolumeQ.ensureHave();
    edgeLengthsQ.ensureHave();
    edgeDihedralAnglesQ.ensureHave();
    actualMetricQ.ensureHave();
    actualCurvatureQ.ensureHave();
    bendingEnergyQ.ensureHave();
    stretchingEnergyQ.ensureHave();
    elasticEnergyQ.ensureHave();
    totalEnergyQ.ensureHave();
    calculate_adjacent_edges_lenght(v);
    calculate_adjucent_faces_area(v);
    calculate_adjacent_edges_dihedral_angles(v);
    calculate_adjacent_faces_metric(v);
    calculate_adjacent_faces_curvature(v);
    calculate_adjacent_faces_energy(v);
    calculate_adjacent_faces_volume(v);
    calculate_adjacent_faces_total_energy(v);    
}

/** @brief Recomputes total energy for all faces adjacent to vertex v. */
void ElasticGeometry::calculate_adjacent_faces_total_energy(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculateFaceTotalEnergy(f);
    }
}

/** @brief Recomputes face volume for all faces adjacent to vertex v. */
void ElasticGeometry::calculate_adjacent_faces_volume(const Vertex& v) {
    for (Face f :v.adjacentFaces()) {
        calculateFaceVolume(f);
    }
}

/**
 * @brief Computes the signed volume of the tetrahedron formed by face f and the origin.
 *
 * V_f = (vertex . normal) * area / 3. Used for pressure work computation.
 *
 * @param f The face.
 */
void ElasticGeometry::calculateFaceVolume(const Face& f) {
    faceAreasQ.ensureHave();
    faceNormalsQ.ensureHave();
    this->faceVolume[f] = dot(this->vertexPositions[f.halfedge().vertex()], this->faceNormals[f]) * this->faceAreas[f]/3.;
}

/** @brief Recomputes edge lengths for all edges incident on vertex v. */
void ElasticGeometry::calculate_adjacent_edges_lenght(const Vertex& v) {
    Vector3 _edgeVec;
    for (Edge e : v.adjacentEdges()) {
        _edgeVec = this->vertexPositions[v] - this->vertexPositions[e.otherVertex(v)];
        this->edgeLengths[e] = _edgeVec.norm();
    }
}


/** @brief Recomputes face areas for all faces adjacent to vertex v via cross product. */
void ElasticGeometry::calculate_adjucent_faces_area(const Vertex& v) {
    Vector3 _edgeVec;
    for (Face f : v.adjacentFaces()) {
        Vector3 normalSum = Vector3::zero();
        for (Halfedge heF : f.adjacentHalfedges()) {
            Halfedge he = heF;
            Vector3 pA = vertexPositions[he.vertex()];
            he = he.next();
            Vector3 pB = vertexPositions[he.vertex()];
            he = he.next();
            Vector3 pC = vertexPositions[he.vertex()];
            
            normalSum += cross(pB - pA, pC - pA); 
            if (he.next() == heF) break;
        }
        faceAreas[f] = 0.5 * normalSum.norm();
    }
}

/**
 * @brief Recomputes face normals and dihedral angles around vertex v.
 *
 * First pass: recomputes face normals for all adjacent faces from the
 * cross product of edge vectors. Second pass: recomputes dihedral angles
 * for all edges of adjacent faces using the updated normals.
 *
 * @param v The vertex whose neighborhood normals/angles need updating.
 */
void ElasticGeometry::calculate_adjacent_edges_dihedral_angles(const Vertex& v) {
    vertexPositionsQ.ensureHave();
    faceNormalsQ.ensureHave();

    for (Face f : v.adjacentFaces()) {
        Vector3 normalSum = Vector3::zero();
        for (Halfedge heF : f.adjacentHalfedges()) {
            Halfedge he = heF;
            Vector3 pA = vertexPositions[he.vertex()];
            he = he.next();
            Vector3 pB = vertexPositions[he.vertex()];
            he = he.next();
            Vector3 pC = vertexPositions[he.vertex()];

            normalSum += cross(pB - pA, pC - pA);
            if (he.next() == heF) break;
        }

        Vector3 normal = unit(normalSum);
        faceNormals[f] = normal;
    }

    for (Face f : v.adjacentFaces()) {
        for (Edge e : f.adjacentEdges()) {
            if (e.isBoundary()) continue;

            if (!e.isManifold()) {
                continue;
            }

            Vector3 N1 = faceNormals[e.halfedge().face()];
            Vector3 N2 = faceNormals[e.halfedge().sibling().face()];
            Vector3 pTail = vertexPositions[e.halfedge().vertex()];
            Vector3 pTip = vertexPositions[e.halfedge().next().vertex()];
            Vector3 edgeDir = unit(pTip - pTail);

            edgeDihedralAngles[e] = atan2(dot(edgeDir, cross(N1, N2)), dot(N1, N2));
        }
    }

}

/** @brief Recomputes the actual metric for all faces adjacent to vertex v. */
void ElasticGeometry::calculate_adjacent_faces_metric(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculate_metric(f);
    }
}

/**
 * @brief Computes the actual metric tensor for a single face from current edge lengths.
 *
 * Metric is stored as (l1^2, l3^2, (l1^2+l3^2-l2^2)/2), all divided by
 * coordinate_scale^2. Edge ordering follows f.adjacentEdges() iteration.
 *
 * @param f The face.
 */
void ElasticGeometry::calculate_metric(const Face& f) {
   Eigen::Vector3f _faceEdgesLengths(3);
   int ind = 0;
   for (Edge e : f.adjacentEdges()) {
       _faceEdgesLengths(ind) = this->edgeLength(e);
        ind += 1;
   }
   this->actualMetric[f][0] = std::pow(_faceEdgesLengths(0), 2.)/ std::pow(coordinate_scale, 2.);
   this->actualMetric[f][1] = std::pow(_faceEdgesLengths(2), 2.)/ std::pow(coordinate_scale, 2.);
   this->actualMetric[f][2] = 0.5 * (std::pow(_faceEdgesLengths(0), 2.) + std::pow(_faceEdgesLengths(2), 2.) -
                                     std::pow(_faceEdgesLengths(1), 2.)) /  std::pow(coordinate_scale, 2.);
}


/**
 * @brief Computes the reference metric tensor for a single face from reference edge lengths.
 *
 * Same vectorization as calculate_metric but uses referenceLengths instead
 * of current edge lengths.
 *
 * @param f The face.
 */
void ElasticGeometry::calculate_reference_metric(const Face& f) {
    Eigen::Vector3f _faceEdgesLengths(3);
    int ind = 0;
    for (Edge e : f.adjacentEdges()) {
        _faceEdgesLengths(ind) = this->referenceLengths[e];
        ind += 1;
    }
    this->referenceMetric[f][0] = std::pow(_faceEdgesLengths(0), 2.) / std::pow(coordinate_scale, 2.);
    this->referenceMetric[f][1] = std::pow(_faceEdgesLengths(2), 2.) / std::pow(coordinate_scale, 2.);
    this->referenceMetric[f][2] = 0.5 * (std::pow(_faceEdgesLengths(0), 2.) + std::pow(_faceEdgesLengths(2), 2.) - std::pow(_faceEdgesLengths(1), 2.)) /
        std::pow(coordinate_scale, 2.);
    
}


/** @brief Recomputes the actual curvature for all faces adjacent to vertex v. */
void ElasticGeometry::calculate_adjacent_faces_curvature(const Vertex& v) {
    Vector3 _curve_comp;
    for (Face f : v.adjacentFaces()) {
        _curve_comp = get_curvature(f, 1);
        this->actualCurvature[f][0] = _curve_comp[0];
        this->actualCurvature[f][1] = _curve_comp[1];
        this->actualCurvature[f][2] = _curve_comp[2];
    }
}




/** @brief Recomputes elastic energy for all faces adjacent to vertex v. */
void ElasticGeometry::calculate_adjacent_faces_energy(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculateFaceEnergy(f);
    }
}

/**
 * @brief Computes the total elastic energy for a single face.
 *
 * E_f = area * (thickness * E_stretch + thickness^3 / 3 * E_bend).
 *
 * @param f The face.
 */
void ElasticGeometry::calculateFaceEnergy(const Face& f) {
    calculate_stretching_energy(f);
    calculate_bending_energy(f);
    elasticEnergy[f] =  1.0* thickness[f] * stretchingEnergy[f];
    elasticEnergy[f] += 1.0 * thickness[f] * thickness[f] * thickness[f] * bendingEnergy[f] / 3.0;
    elasticEnergy[f] *= faceAreas[f];    
}
 
/**
 * @brief Computes total energy (elastic - pressure * volume) for a single face.
 * @param f The face.
 */
void ElasticGeometry::calculateFaceTotalEnergy(const Face& f) {
    requireElasticEnergy();
    requireFaceVolume();
    totalEnergy[f] =1.*elasticEnergy[f] - 1. * this->pressure * faceVolume[f];
}

/**
 * @brief Computes the 2D stretching energy density for a single face.
 *
 * E_stretch = (g - bar{g})^T A (g - bar{g}), where g and bar{g} are
 * the actual and reference metric vectors and A is the elastic Cauchy tensor.
 *
 * @param f The face.
 */
void ElasticGeometry::calculate_stretching_energy(const Face& f) {
    actualMetricQ.ensureHave();
    referenceMetricQ.ensureHave();
    elasticCauchyTensorQ.ensureHave();
    Eigen::Vector3f _metricDiff = actualMetric[f] - referenceMetric[f];

    this->stretchingEnergy[f] = (elasticCauchyTensor[f](0, 0) * _metricDiff[0] * _metricDiff[0] +
                          elasticCauchyTensor[f](1, 1) * _metricDiff[1] * _metricDiff[1] +
                          2. * elasticCauchyTensor[f](2, 2) * _metricDiff[2] * _metricDiff[2] +
                          2. * elasticCauchyTensor[f](1, 0) * _metricDiff[0] * _metricDiff[1] +
                          4. * elasticCauchyTensor[f](2, 0) * _metricDiff[0] * _metricDiff[2] +
                          4. * elasticCauchyTensor[f](2, 1) * _metricDiff[1] * _metricDiff[2]);
}


/**
 * @brief Computes the bending energy density for a single face.
 *
 * E_bend = (b - bar{b})^T A (b - bar{b}), where b and bar{b} are
 * the actual and reference curvature vectors and A is the elastic Cauchy tensor.
 *
 * @param f The face.
 */
void ElasticGeometry::calculate_bending_energy(const Face& f) {
    actualMetricQ.ensureHave();
    referenceMetricQ.ensureHave();
    elasticCauchyTensorQ.ensureHave();
    referenceCurvatureQ.ensureHave();
    actualCurvatureQ.ensureHave();

    Eigen::Vector3f _curvDiff = actualCurvature[f] - referenceCurvature[f];

    this->bendingEnergy[f] = elasticCauchyTensor[f](0, 0) * _curvDiff[0] * _curvDiff[0] +
                          elasticCauchyTensor[f](1, 1) * _curvDiff[1] * _curvDiff[1] +
                          2.0 * elasticCauchyTensor[f](2, 2) * _curvDiff[2] * _curvDiff[2] +
                          2.0 * elasticCauchyTensor[f](1, 0) * _curvDiff[0] * _curvDiff[1] +
                          4.0 * elasticCauchyTensor[f](2, 0) * _curvDiff[0] * _curvDiff[2] +
                          4.0 * elasticCauchyTensor[f](2, 1) * _curvDiff[1] * _curvDiff[2];
}

/** @brief No-op: pressure is set externally. */
void ElasticGeometry::computePressure() {}

/** @brief No-op placeholder: region computation not yet implemented. */
void ElasticGeometry::computeRegions() {}

/** @brief No-op placeholder: fixed vertex computation not yet implemented. */
void ElasticGeometry::computeFixedVertexs() {}

/** @brief No-op placeholder: fixed angle computation not yet implemented. */
void ElasticGeometry::computeFixedAngles() {}


// ========== Curvature Query Functions ==========

/** @brief Returns the mean curvature at face f computed from the reference metric and curvature. */
double ElasticGeometry::getReferenceMeanCurvautre(Face f) {
    referenceCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getMean(referenceMetric[f], referenceCurvature[f]);
}
/** @brief Returns the Gaussian curvature at face f from reference metric and curvature. */
double ElasticGeometry::getReferenceGaussianCurvautre(Face f) {
    referenceCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getDet(referenceMetric[f], referenceCurvature[f]);
}
/** @brief Returns the mean curvature at face f computed from the actual metric and curvature. */
double ElasticGeometry::getActualMeanCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    actualMetricQ.ensureHave();
    return getMean(actualMetric[f], actualCurvature[f]);
}
/** @brief Returns the Gaussian curvature at face f from the actual metric and curvature. */
double ElasticGeometry::getActualGaussianCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    actualMetricQ.ensureHave();
    return getDet(actualMetric[f], actualCurvature[f]);
}
/**
 * @brief Returns the "semi-actual" mean curvature at face f.
 *
 * Uses the reference metric with the actual curvature tensor, giving
 * a hybrid measure useful for comparing deformed curvature against
 * the undeformed metric.
 */
double ElasticGeometry::getSemiActualMeanCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getMean(referenceMetric[f], actualCurvature[f]);
}
/**
 * @brief Returns the "semi-actual" Gaussian curvature at face f.
 *
 * Uses the reference metric with the actual curvature tensor, giving
 * a hybrid determinant measure for comparing deformed curvature against
 * the undeformed metric.
 */
double ElasticGeometry::getSemiActualGaussianCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getDet(referenceMetric[f], actualCurvature[f]);
}


/**
 * @brief Computes the shape operator for every face from the actual state.
 *
 * The shape operator S = a^{-1} b is computed per-face as a 2x2 matrix
 * stored in a Vector4f (row-major: S00, S11, S01, S10), where a is the
 * actual metric tensor and b is the actual curvature tensor. The inverse
 * is computed analytically using the 2x2 metric determinant.
 */
void ElasticGeometry::getActualShape() {
    actualShape = FaceData<Eigen::Vector4f>(mesh, Eigen::Vector4f(0, 0, 0, 0));
    actualCurvatureQ.ensureHave();
    actualMetricQ.ensureHave();
    Eigen::Vector3f a;
    Eigen::Vector3f b;
    for (Face f : this->mesh.faces()) {
        a = actualMetric[f];
        b = actualCurvature[f];
        actualShape[f] = Eigen::Vector4f((a[1] * b[0] - a[2] * b[2]) / (a[0] * a[1] - a[2] * a[2]),
                                         (a[0] * b[1] - a[2] * b[2]) / (a[0] * a[1] - a[2] * a[2]),
                                         -(a[2] * b[1] - a[1] * b[2]) / (a[0] * a[1] - a[2] * a[2]),
                                         -(a[2] * b[0] - a[0] * b[2]) / (a[0] * a[1] - a[2] * a[2]));
    }
}

/**
 * @brief Stores the first two halfedge indices of each face as its local basis edges.
 *
 * For each face, records the indices of the first two adjacent halfedges
 * into baseEdges[f] as an (ind1, ind2) pair. These two edges define the
 * local tangent frame used for metric and curvature computations.
 */
void ElasticGeometry::getFaceBasis() {
    baseEdges = FaceData<Eigen::Vector2i>(mesh, Eigen::Vector2i(-1,-1));    
    int ind1;
    int ind2;
    for (Face f : this->mesh.faces()) {
        ind1 = -1;
        ind2 = -1;
        for (Halfedge HE : f.adjacentHalfedges()) {
            if (ind1 == -1) {
                ind1 = HE.getIndex();
            }
            else if (ind2 == -1) {
                ind2 = HE.getIndex();
            }            
        }
        baseEdges[f] = Eigen::Vector2i(ind1, ind2);
    }
}

/**
 * @brief Computes the mean curvature H = tr(a^{-1} b) / 2 from metric a and curvature b.
 *
 * @param a Metric tensor stored as (a11, a22, a12).
 * @param b Curvature tensor stored as (b11, b22, b12).
 * @return The mean curvature scalar.
 */
double ElasticGeometry::getMean(Eigen::Vector3f a, Eigen::Vector3f b) {
    return (a[1] * b[0] - 2 * a[2] * b[2] + a[0] * b[1]) / (a[0] * a[1] - a[2] * a[2]) / 2.0;
}
/**
 * @brief Computes the Gaussian curvature K = det(b) / det(a) from metric a and curvature b.
 *
 * @param a Metric tensor stored as (a11, a22, a12).
 * @param b Curvature tensor stored as (b11, b22, b12).
 * @return The Gaussian curvature scalar.
 */
double ElasticGeometry::getDet(Eigen::Vector3f a, Eigen::Vector3f b) {
    return (b[0] * b[1] - b[2] * b[2]) / (a[0] * a[1] - a[2] * a[2]);
}

} // namespace surface
} // namespace geometrycentral