#include "ElasticGeometry.h" // my class to be implemented
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
    referenceEdgeDihedralAngleQ         (&referenceEdgeDihedralAngle,       std::bind(&ElasticGeometry::computeReferenceEdgeDihedralAngle,this),                quantities),
    referenceMetricQ                    (&referenceMetric,                  std::bind(&ElasticGeometry::computeReferenceMetric, this),                          quantities),
    actualMetricQ                       (&actualMetric,                     std::bind(&ElasticGeometry::computeActualMetric, this),                             quantities),
    referenceCurvatureQ                 (&referenceCurvature,               std::bind(&ElasticGeometry::computeReferenceCurvature, this),                       quantities),
    actualCurvatureQ                    (&actualCurvature,                  std::bind(&ElasticGeometry::computeActualCurvature, this),                          quantities),
    elasticCauchyTensorQ                (&elasticCauchyTensor,              std::bind(&ElasticGeometry::computeElasticCauchyTensor, this),                      quantities),
    thicknessQ                          (&thickness,                        std::bind(&ElasticGeometry::computeThickness, this),                                quantities),
    youngsModulusQ                      (&youngsModulus,                    std::bind(&ElasticGeometry::computeYoungsModulus, this),                            quantities),
    poissonsRatioQ                      (&poissonsRatio,                    std::bind(&ElasticGeometry::computePoissonsRatio, this),                            quantities),
    elasticEnergyQ                      (&elasticEnergy,                    std::bind(&ElasticGeometry::computeElasticEnergy, this),                            quantities),
    pressureQ                           (&pressure,                         std::bind(&ElasticGeometry::computePressure,this),                                  quantities),
    regionsQ                            (&regions,                          std::bind(&ElasticGeometry::computeRegions, this),                                  quantities),
    fixedVertexesQ                      (&fixedVertexes,                    std::bind(&ElasticGeometry::computeFixedVertexs,this),                              quantities),
    fixedAnglesQ                        (&fixedAngles,                      std::bind(&ElasticGeometry::computeFixedAngles,this),                               quantities)    
    {

         //std::cout << "Elastic Geometry!";
         vertexPositions = inputVertexPositions_;// VertexData<Vector3>(mesh_, Vector3{0., 0., 0.});
         // The input vertex positions share storage with vertexPositions, incremented the required counter and make sure they never get cleared
         requireVertexPositions();
         vertexPositionsQ.clearable = false;

         referenceLengths=L_bar_;
         requireReferenceLegths();
         referenceLengthsQ.clearable = false;
        
         referenceEdgeDihedralAngle=B_bar_;
         requireEdgeDihedralAngles();
         referenceEdgeDihedralAngleQ.clearable=false;

         thickness=THICKNESS_;
         requireThickness();
         thicknessQ.clearable=false;

         elasticCauchyTensor=ElasticTensor_;
         requireElasticCauchyTensor();
         elasticCauchyTensorQ.clearable=false;

         pressure=PRESSURE_;
         requirePressure();
         pressureQ.clearable=false;


         // Also compute essential stuff? (consider not)//
         requireReferenceMetric();
         requireReferenceCurvature();
         requireActualMetric();
         requireActualCurvature();

         requireElasticCauchyTensor();
         
         //Dont require energy of not requested specifically.
      }

// clang-format on


// Simplest of all contructors - only mesh. no position, no nothin (everything is 0)
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_)
    : ElasticGeometry::ElasticGeometry(mesh_, VertexData<Vector3>(mesh_, Vector3{0., 0., 0.}),
                                       EdgeData<double>(mesh_, 0), EdgeData<double>(mesh_, 0),
                                       FaceData<double>(mesh_, 0), FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()),
                                       0) {
}

// Not so basic after all - mesh with coordinates. Creates a basic compatible without energy or thickness
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_)
    : ElasticGeometry::ElasticGeometry(mesh_, inputVertexPositions_,
                                       EdgeData<double>(mesh_, 0), EdgeData<double>(mesh_, 0),
                                       FaceData<double>(mesh_, 0), FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()),
                                       0) {
    this->requireReferenceLegths();
    this->requireReferenceEdgeDihedralAngle();
 }


// Notbasic - mesh with coordinates, thickness and Youngs modulus - creating a compatible isotropic geometry, no pressure though need to implement (at least partially) the constructor)
 ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                                  const double& THICKNESS_, const double& YOUNGs_, const double& POISSONs_)
     : ElasticGeometry::ElasticGeometry(mesh_, inputVertexPositions_, EdgeData<double>(mesh_, 0),
                                        EdgeData<double>(mesh_, 0), FaceData<double>(mesh_, THICKNESS_),
                                        FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()), 0),
       youngsModulus()
       {
     this->requireReferenceLegths();
     this->requireReferenceEdgeDihedralAngle();
     this->requireReferenceCurvature();
     this->requireReferenceMetric();



     



     
 }




 
void ElasticGeometry::requireReferenceLegths() {
    referenceLengthsQ.require();
}
void ElasticGeometry::unrequireReferenceLegths() {
    referenceLengthsQ.unrequire();
}

void ElasticGeometry::requireReferenceEdgeDihedralAngle() {
    referenceEdgeDihedralAngleQ.require();
}
void ElasticGeometry::unrequireReferenceEdgeDihedralAngle() {
    referenceEdgeDihedralAngleQ.unrequire();
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


void ElasticGeometry::requirePoissonsratio() {
    poissonsRatioQ.require();
}
void ElasticGeometry::unrequirePoissonsratio() {
    poissonsRatioQ.unrequire();
}


void ElasticGeometry::requireElasticEnergy() {
    elasticEnergyQ.require();
}
void ElasticGeometry::unrequireElasticEnergy() {
    elasticEnergyQ.unrequire();
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


/// HELPER functions

// Checks if the a general MeshData has legal values or not.  Since it is a general MeshData we use a template "data_type", to catch'em all.
template <typename data_type>
static bool is_illegal(data_type& data) {
    bool any_zeros = false;
    for (int index = 0; index < data.size(); index++) {
        any_zeros += data[index] == 0;
    }
    return any_zeros;
}

// "COMPUTE" functions. implement!



void ElasticGeometry::computeReferenceLengths() {
    //     If reference lengths are not given (i.e. they are zero) conpute them once at initialization only, set them to
    //     be actual values.  Otherwise, we are not yet changing reference values, so there is no need to compute
    //     anything.

     if (is_illegal(referenceLengths)){//any_zeros) { // Indicating illegal data.
         // Calculate all reference lengths, not just those that are zero.  ####### CONSIDER CHANGING THIS #######
         this->requireEdgeLengths();
         for (Edge e : this->mesh.edges()) {
             referenceLengths[e] = this->edgeLengths(e);
         }        
     }
}


void ElasticGeometry::computeReferenceEdgeDihedralAngle() {
    // Same logic as for reference lengths ubove.
    if (is_illegal(referenceEdgeDihedralAngle)) {
        this->requireEdgeDihedralAngles();
        for (Edge e : this->mesh.edges()) {
            referenceEdgeDihedralAngle[e] = this->edgeDihedralAngle(e);
        }
    }
} 




void ElasticGeometry::computeReferenceMetric() { //CONSIDER delegating the calculation inside to an external, more general and morr readable function.
    referenceMetric = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(1., 2., 3.));
    Eigen::Vector3f _faceEdgesLengths(3);
    for (Face f : this->mesh.faces()) {
        int ind = 0;
        for (Edge e : f.adjacentEdges()) {
            _faceEdgesLengths(ind) = this->referenceLengths[e];
            ind += 1;
        }
        referenceMetric[f][0] = std::pow(_faceEdgesLengths(0), 2);
        referenceMetric[f][1] = std::pow(_faceEdgesLengths(1), 2);
        referenceMetric[f][2] = 0.5 * (std::pow(_faceEdgesLengths(0), 2) + std::pow(_faceEdgesLengths(1), 2) -
                                       std::pow(_faceEdgesLengths(2), 2));
    }
}

void ElasticGeometry::computeActualMetric() {
    Eigen::Vector3f _faceEdgesLengths(3);
    for (Face f : this->mesh.faces()) {
        int ind = 0;
        for (Edge e : f.adjacentEdges()) {
            _faceEdgesLengths(ind) = this->edgeLength(e);
            ind += 1;
        }
        actualMetric[f][0] = std::pow(_faceEdgesLengths(0), 2);
        actualMetric[f][1] = std::pow(_faceEdgesLengths(1), 2);
        actualMetric[f][2] = 0.5 * (std::pow(_faceEdgesLengths(0), 2) + std::pow(_faceEdgesLengths(1), 2) -
                                       std::pow(_faceEdgesLengths(2), 2));
    }
}

void ElasticGeometry::computeReferenceCurvature() {}

void ElasticGeometry::computeActualCurvature() {}

void ElasticGeometry::computeElasticCauchyTensor() { ///
    bool _isIlliegal =
        this->elasticCauchyTensor.toVector().isZero(); // simple test currently if all entries are zero. ## WILL NEED TO CHANGE THIS FOR A VALIDITY TEST! ##
    if (_isIlliegal) {
        float _Atensor[6];
        float _invmet[3];
        float _det;
        for (Face f : this->mesh.faces()) {
            _det = referenceMetric[f][0] * referenceMetric[f][1] - referenceMetric[f][2] * referenceMetric[f][2];
            _invmet[0] = referenceMetric[f][1]/_det;
            _invmet[1] = referenceMetric[f][0] / _det;
            _invmet[2] = -referenceMetric[f][3] / _det;
            elasticCauchyTensor[f](0, 0) = youngsModulus * _invmet[0] * _invmet[0];
            elasticCauchyTensor[f](0, 1) = _invmet[2] * _invmet[2] * (1- poissonsRatio);
            elasticCauchyTensor[f](1, 0) = ;
            elasticCauchyTensor[f](0, 2) = ;
            elasticCauchyTensor[f](2, 0) = ;
            elasticCauchyTensor[f](1, 1) = ;
            elasticCauchyTensor[f](1, 2) = ;
            elasticCauchyTensor[f](2, 1) = ;
            elasticCauchyTensor[f](2, 2) = ;

        }
    }    
}

void ElasticGeometry::computeThickness() {}

void ElasticGeometry::computeYoungsModulus() {}

void ElasticGeometry::computePoissonsRatio() {}

void ElasticGeometry::computeElasticEnergy() {}

void ElasticGeometry::computePressure() {}


void ElasticGeometry::computeRegions() {}

void ElasticGeometry::computeFixedVertexs() {}

void ElasticGeometry::computeFixedAngles() {}




} // namespace surface
} // namespace geometrycentral