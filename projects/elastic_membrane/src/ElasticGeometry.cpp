#include "ElasticGeometry.h"
#include <fstream>
#include <limits>


using namespace geometrycentral;
using namespace geometrycentral::surface;


namespace geometrycentral {
namespace surface {


// clang-format off


/// <summary>
/// Main, principal constructor.
/// </summary>
/// <param name="mesh_"> the  surfave mesh object relatedto this geometry</param>
/// <param name="inputVertexPositions_"> as he name implies the vertex position, should be a VertexData calss  </param>
/// <param name="L_bar_"> EdgdeData - reference legnths </param>
/// <param name="B_bar_"> EdgdeData - reference curvatures </param>
/// <param name="THICKNESS_"> FaceData - thickness of each face</param>
/// <param name="ElasticTensor_"> FaceData - Elastic tensor, in form of a 3X3</param>
/// <param name="PRESSURE_"> pressure</param>
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


ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_)
    : ElasticGeometry::ElasticGeometry(mesh_, VertexData<Vector3>(mesh_, Vector3{0., 0., 0.}),
                                       EdgeData<double>(mesh_, 0.), EdgeData<double>(mesh_, 0.),
                                       FaceData<double>(mesh_, 0.), FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()),
                                       0.) {}

ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_)
    : ElasticGeometry::ElasticGeometry(mesh_, inputVertexPositions_, EdgeData<double>(mesh_, 0.),
                                       EdgeData<double>(mesh_, 0.), FaceData<double>(mesh_, 0.),
                                       FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()), 0.) {
    this->requireReferenceLegths();
    this->requireReferenceEdgeDihedralAngles();
}


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




template <typename data_type>
static bool is_illegal(data_type& data) {
    bool any_zeros = false;
    for (int index = 0; index < data.size(); index++) {
        any_zeros += data[index] == 0.0;
    }
    return any_zeros;
}



double project3(const Vector3& vec1, const Vector3& vec2) {
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

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


void ElasticGeometry::computeReferenceLengths() {
    if (referenceLengths.toVector().isZero()) {
        this->requireEdgeLengths();
        for (Edge e : this->mesh.edges()) {
            referenceLengths[e] = this->edgeLength(e);
        }
    }
}


void ElasticGeometry::computeReferenceEdgeDihedralAngles() {
    if (referenceEdgeDihedralAngles.toVector().isZero()) {
        this->requireEdgeDihedralAngles();
        for (Edge e : this->mesh.edges()) {
            referenceEdgeDihedralAngles[e] = this->edgeDihedralAngle(e);
        }
    }
}


void ElasticGeometry::computeReferenceMetric() {
    for (Face f : this->mesh.faces()) {
        calculate_reference_metric(f);
    }
}

void ElasticGeometry::computeActualMetric() {
    for (Face f : this->mesh.faces()) {
        calculate_metric(f);
    }
}

void ElasticGeometry::computeReferenceCurvature() {
    this->faceNormalsQ.ensureHave();

    for (Face f : this->mesh.faces()) {
        Vector3 _curve_comp = get_curvature(f, 0);
        this->referenceCurvature[f][0] = _curve_comp[0];
        this->referenceCurvature[f][1] = _curve_comp[1];
        this->referenceCurvature[f][2] = _curve_comp[2];
    }
} 

void ElasticGeometry::computeActualCurvature() {
    this->faceNormalsQ.ensureHave();
    for (Face f : this->mesh.faces()) {
        Vector3 _curve_comp = get_curvature(f, 1);
        this->actualCurvature[f][0] = _curve_comp[0];
        this->actualCurvature[f][1] = _curve_comp[1];
        this->actualCurvature[f][2] = _curve_comp[2];
    }
}

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

void ElasticGeometry::computeThickness() {}

void ElasticGeometry::computeYoungsModulus() {}

void ElasticGeometry::computePoissonsRatio() {}

void  ElasticGeometry::computeElasticEnergy() {
    this->requireBendingEnergy();
    this->requireStretchingEnergy();
    if (elasticEnergy.size() == 0) elasticEnergy = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces())
    {
        calculateFaceEnergy(f);
    }
}


void ElasticGeometry::computeTotalEnergy() {
    this->requireElasticEnergy();
    this->requireFaceVolume();
    if (totalEnergy.size() == 0) totalEnergy = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces()) {
        calculateFaceTotalEnergy(f);
    }    
}


void ElasticGeometry::computeFaceVolume() { //assuming centered objects    
    this->requireFaceNormals();
    if (faceVolume.size() == 0) faceVolume = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces()) {
        calculateFaceVolume(f);                
       
    }
}

void  ElasticGeometry::computeStretchingEnergy() {
    this->requireFaceAreas();
    if (stretchingEnergy.size() == 0) stretchingEnergy = FaceData<double>(this->mesh, 0.);
    for (Face f : this->mesh.faces()) {
        calculate_stretching_energy(f);
    }
}

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

void ElasticGeometry::calculate_adjacent_faces_total_energy(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculateFaceTotalEnergy(f);
    }
}

void ElasticGeometry::calculate_adjacent_faces_volume(const Vertex& v) {
    for (Face f :v.adjacentFaces()) {
        calculateFaceVolume(f);
    }
}

void ElasticGeometry::calculateFaceVolume(const Face& f) {
    faceAreasQ.ensureHave();
    faceNormalsQ.ensureHave();
    this->faceVolume[f] = dot(this->vertexPositions[f.halfedge().vertex()], this->faceNormals[f]) * this->faceAreas[f]/3.;
}

void ElasticGeometry::calculate_adjacent_edges_lenght(const Vertex& v) {
    Vector3 _edgeVec;
    for (Edge e : v.adjacentEdges()) {
        _edgeVec = this->vertexPositions[v] - this->vertexPositions[e.otherVertex(v)];
        this->edgeLengths[e] = _edgeVec.norm();
    }
}


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

void ElasticGeometry::calculate_adjacent_faces_metric(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculate_metric(f);
    }
}

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


void ElasticGeometry::calculate_adjacent_faces_curvature(const Vertex& v) {
    Vector3 _curve_comp;
    for (Face f : v.adjacentFaces()) {
        _curve_comp = get_curvature(f, 1);
        this->actualCurvature[f][0] = _curve_comp[0];
        this->actualCurvature[f][1] = _curve_comp[1];
        this->actualCurvature[f][2] = _curve_comp[2];
    }
}




void ElasticGeometry::calculate_adjacent_faces_energy(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculateFaceEnergy(f);
    }
}

void ElasticGeometry::calculateFaceEnergy(const Face& f) {
    calculate_stretching_energy(f);
    calculate_bending_energy(f);
    elasticEnergy[f] =  1.0* thickness[f] * stretchingEnergy[f];
    elasticEnergy[f] += 1.0 * thickness[f] * thickness[f] * thickness[f] * bendingEnergy[f] / 3.0;
    elasticEnergy[f] *= faceAreas[f];    
}
 
void ElasticGeometry::calculateFaceTotalEnergy(const Face& f) {
    requireElasticEnergy();
    requireFaceVolume();
    totalEnergy[f] =1.*elasticEnergy[f] - 1. * this->pressure * faceVolume[f];
}

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

void ElasticGeometry::computePressure() {}

void ElasticGeometry::computeRegions() {}

void ElasticGeometry::computeFixedVertexs() {}

void ElasticGeometry::computeFixedAngles() {}


double ElasticGeometry::getReferenceMeanCurvautre(Face f) {
    referenceCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getMean(referenceMetric[f], referenceCurvature[f]);
}
double ElasticGeometry::getReferenceGaussianCurvautre(Face f) {
    referenceCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getDet(referenceMetric[f], referenceCurvature[f]);
}
double ElasticGeometry::getActualMeanCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    actualMetricQ.ensureHave();
    return getMean(actualMetric[f], actualCurvature[f]);
}
double ElasticGeometry::getActualGaussianCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    actualMetricQ.ensureHave();
    return getDet(actualMetric[f], actualCurvature[f]);
}
double ElasticGeometry::getSemiActualMeanCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getMean(referenceMetric[f], actualCurvature[f]);
}
double ElasticGeometry::getSemiActualGaussianCurvautre(Face f) {
    actualCurvatureQ.ensureHave();
    referenceMetricQ.ensureHave();
    return getDet(referenceMetric[f], actualCurvature[f]);
}


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

double ElasticGeometry::getMean(Eigen::Vector3f a, Eigen::Vector3f b) {
    return (a[1] * b[0] - 2 * a[2] * b[2] + a[0] * b[1]) / (a[0] * a[1] - a[2] * a[2]) / 2.0;
}
double ElasticGeometry::getDet(Eigen::Vector3f a, Eigen::Vector3f b) {
    return (b[0] * b[1] - b[2] * b[2]) / (a[0] * a[1] - a[2] * a[2]);
}

} // namespace surface
} // namespace geometrycentral