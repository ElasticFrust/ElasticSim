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
    referenceEdgeDihedralAnglesQ        (&referenceEdgeDihedralAngles,      std::bind(&ElasticGeometry::computeReferenceEdgeDihedralAngles,this),                quantities),
    referenceMetricQ                    (&referenceMetric,                  std::bind(&ElasticGeometry::computeReferenceMetric, this),                          quantities),
    actualMetricQ                       (&actualMetric,                     std::bind(&ElasticGeometry::computeActualMetric, this),                             quantities),
    referenceCurvatureQ                 (&referenceCurvature,               std::bind(&ElasticGeometry::computeReferenceCurvature, this),                       quantities),
    actualCurvatureQ                    (&actualCurvature,                  std::bind(&ElasticGeometry::computeActualCurvature, this),                          quantities),
    elasticCauchyTensorQ                (&elasticCauchyTensor,              std::bind(&ElasticGeometry::computeElasticCauchyTensor, this),                      quantities),
    thicknessQ                          (&thickness,                        std::bind(&ElasticGeometry::computeThickness, this),                                quantities),
    youngsModulusQ                      (&youngsModulus,                    std::bind(&ElasticGeometry::computeYoungsModulus, this),                            quantities),
    poissonsRatioQ                      (&poissonsRatio,                    std::bind(&ElasticGeometry::computePoissonsRatio, this),                            quantities),
    elasticEnergyQ                      (&elasticEnergy,                    std::bind(&ElasticGeometry::computeElasticEnergy, this),                            quantities),
    stretchingEnergyQ                   (&stretchingEnergy,                 std::bind(&ElasticGeometry::computeStretchingEnergy, this),                         quantities),
    bendingEnergyQ                      (&bendingEnergy,                    std::bind(&ElasticGeometry::computeBendingEnergy, this),                            quantities),
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
        
         referenceEdgeDihedralAngles=B_bar_;
         requireEdgeDihedralAngles();
         referenceEdgeDihedralAnglesQ.clearable=false;

         thickness=THICKNESS_;
         requireThickness();
         thicknessQ.clearable=false;

         youngsModulus= FaceData<double>(mesh_, 0);
         requireYoungsModulus();

         poissonsRatio= FaceData<double>(mesh_, 0);
         requirePoissonsRatio();

         elasticCauchyTensor=ElasticTensor_;
         requireElasticCauchyTensor();
         //elasticCauchyTensorQ.clearable=false;

         pressure=PRESSURE_;
         requirePressure();
         pressureQ.clearable=false;

        




         // Also compute essential stuff? (consider not)//
         referenceMetric = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(1., 2., 3.));
         requireReferenceMetric();

         referenceCurvature = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(1., 2., 3.));
         requireReferenceCurvature();

         actualMetric = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(0., 0., 0.));
         actualCurvature = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(0., 0., 0.));
         //elasticCauchyTensor=FaceData<Eigen::Matrix3f>(this->mesh, Eigen::Matrix3f());





         float xsum=0;
         float ysum=0;
         float zsum=0;
         for (Vertex v : this->mesh.vertices()) {
           xsum += inputVertexPositions[v].x;
           ysum += inputVertexPositions[v].y;
           zsum += inputVertexPositions[v].z;
         }
         if(xsum+ysum+zsum!=0){
            requireActualMetric();
            requireActualCurvature();
            //requireElasticCauchyTensor();
         }
          
         
         stretchingEnergy = FaceData<double>(mesh_, 0);
         bendingEnergy = FaceData<double>(mesh_, 0);
         elasticEnergy = FaceData<double>(mesh_, 0);
         requireStretchingEnergy();
         requireBendingEnergy();
         requireElasticEnergy();

      }

// clang-format on


// Simplest of all contructors - only mesh. no position, no nothin (everything is 0)
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_)
    : ElasticGeometry::ElasticGeometry(mesh_, VertexData<Vector3>(mesh_, Vector3{0., 0., 0.}),
                                       EdgeData<double>(mesh_, 0), EdgeData<double>(mesh_, 0),
                                       FaceData<double>(mesh_, 0), FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()),
                                       0) {}

// Not so basic after all - mesh with coordinates. Creates a basic compatible without energy or thickness
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_)
    : ElasticGeometry::ElasticGeometry(mesh_, inputVertexPositions_, EdgeData<double>(mesh_, 0),
                                       EdgeData<double>(mesh_, 0), FaceData<double>(mesh_, 0),
                                       FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()), 0) {
    this->requireReferenceLegths();
    this->requireReferenceEdgeDihedralAngles();
}


// Notbasic - mesh with coordinates, thickness and Youngs modulus - creating a compatible isotropic geometry, no
// pressure though need to implement (at least partially) the constructor)
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                                 const double& THICKNESS_, const double& YOUNGs_, const double& POISSONs_,
                                 const double& PRESSURE_)
    : ElasticGeometry::ElasticGeometry(mesh_, inputVertexPositions_, EdgeData<double>(mesh_, 0),
                                       EdgeData<double>(mesh_, 0), FaceData<double>(mesh_, THICKNESS_),
                                       FaceData<Eigen::Matrix3f>(mesh_, Eigen::Matrix3f()), PRESSURE_) {
    // THE ABOVE CALL toto the general contructor creates a compatible elasticmembrane withtout any rigidity. Following,
    // we implement an elastic tensor

    youngsModulus = FaceData<double>(mesh_, YOUNGs_);
    requireYoungsModulus();
    youngsModulusQ.clearable = false;

    poissonsRatio = FaceData<double>(mesh_, POISSONs_);
    requirePoissonsRatio();
    poissonsRatioQ.clearable = false;

    // After creating the relevant poisson ratio and young modulus values. Creat the defult elastic tensor-
    // std::cout << "\nCalling Cauchy calc...\n";
    unrequireElasticCauchyTensor();
    elasticCauchyTensorQ.clearIfNotRequired();
    elasticCauchyTensor = FaceData<Eigen::Matrix3f>(this->mesh, Eigen::Matrix3f());
    requireElasticCauchyTensor();
    elasticCauchyTensorQ.clearable = false;



    //Genereate the energy
    unrequireStretchingEnergy();
    unrequireBendingEnergy();
    unrequireElasticEnergy();
    stretchingEnergyQ.clearIfNotRequired();
    bendingEnergyQ.clearIfNotRequired();
    elasticEnergyQ.clearIfNotRequired();
    stretchingEnergy = FaceData<double>(mesh_, 0);
    bendingEnergy = FaceData<double>(mesh_, 0);
    elasticEnergy = FaceData<double>(mesh_, 0);
    requireStretchingEnergy();
    requireBendingEnergy();
    requireElasticEnergy();
    stretchingEnergyQ.clearable = false;
    bendingEnergyQ.clearable = false;
    elasticEnergyQ.clearable = false;

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


/// HELPER functions

// Checks if the a general MeshData has legal values or not.  Since it is a general MeshData we use a template
// "data_type", to catch'em all.
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

Vector3  ElasticGeometry::get_curvature(Face& _f, const int& _ref_or_act) {
    EdgeData<double> angles;
    if (_ref_or_act == 0)
    {
       angles = this->referenceEdgeDihedralAngles;
    }
    else
    {
       angles = this->edgeDihedralAngles;
    }

    float _totLength = 0;
    float _curve_comp1 = 0; // in the 1 direction
    float _curve_comp2 = 0; // in the 2 direction
    float _curve_comp3 = 0; // in the 1-2 "direction"
    Vector3 _temp_edge_vec;
    Vector3 _dual_edge_vec;
    Vector3 _basis1 = this->faceTangentBasis[_f][0].normalize(); // Do we need to normalize? Be sure.
    Vector3 _basis2 = this->faceTangentBasis[_f][1].normalize();
    double _proj1;
    double _proj2;
    for (Edge e : _f.adjacentEdges()) {
        _totLength += referenceLengths[e];
        double _curvature_magnitude = (angles[e] - std::_Pi_val) / this->edgeLength(e) / this->edgeCotanWeight(e);
        _temp_edge_vec = this->vertexPositions[e.firstVertex()] - this->vertexPositions[e.secondVertex()];
        _temp_edge_vec = _temp_edge_vec.normalize();
        _dual_edge_vec = _temp_edge_vec.rotateAround(this->faceNormals[_f], std::_Pi_val / 2);
        _proj1 = project3(_basis1, _dual_edge_vec);
        _proj2 = project3(_basis2, _dual_edge_vec);
        _curve_comp1 += _curvature_magnitude * _proj1 * _proj1;
        _curve_comp2 += _curvature_magnitude * _proj2 * _proj2;
        _curve_comp3 += _curvature_magnitude * _proj1 * _proj2;
    }
    Vector3 res{_curve_comp1, _curve_comp2, _curve_comp3};
    return res;    
}

// "COMPUTE" functions. implement!


void ElasticGeometry::computeReferenceLengths() {
    //     If reference lengths are not given (i.e. they are zero) conpute them once at initialization only, set them to
    //     be actual values.  Otherwise, we are not yet changing reference values, so there is no need to compute
    //     anything.

    if (referenceLengths.toVector().isZero()) { // any_zeros) { // Indicating illegal data.
        // Calculate all reference lengths, not just those that are zero.  ####### CONSIDER CHANGING THIS #######
        this->requireEdgeLengths();
        for (Edge e : this->mesh.edges()) {
            referenceLengths[e] = this->edgeLength(e);
        }
    }
}


void ElasticGeometry::computeReferenceEdgeDihedralAngles() {
    // Same logic as for reference lengths ubove.
    if (referenceEdgeDihedralAngles.toVector().isZero()) {
        this->requireEdgeDihedralAngles();
        for (Edge e : this->mesh.edges()) {
            referenceEdgeDihedralAngles[e] = this->edgeDihedralAngle(e);
        }
    }
}


void ElasticGeometry::computeReferenceMetric() { // CONSIDER delegating the calculation inside to an external, more
                                                 // general and morr readable function.
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
        calculate_metric(f);
        /*int ind = 0;
        for (Edge e : f.adjacentEdges()) {
            _faceEdgesLengths(ind) = this->edgeLength(e);
            ind += 1;
        }
        actualMetric[f][0] = std::pow(_faceEdgesLengths(0), 2);
        actualMetric[f][1] = std::pow(_faceEdgesLengths(1), 2);
        actualMetric[f][2] = 0.5 * (std::pow(_faceEdgesLengths(0), 2) + std::pow(_faceEdgesLengths(1), 2) -
                                    std::pow(_faceEdgesLengths(2), 2));*/
    }
}

void ElasticGeometry::computeReferenceCurvature() { // For a dihedral angle \theta, the cruvature along that direction
                                                    // is \theta/(length dual). and lengof dual = lengh* (edge cotan
                                                    // weight). Direction is easy perp. to edge. Curvature is then sum o
                                                    // direct_curve * eXe single edge (e) true curvatrue  =
                                                    // Dihedralangle/length of dual (l*).  ==> k= D/l* a single edge
                                                    // curvature tensor  =   k  e X e total curvature tensor of face S =
                                                    // \sum_e l/L_tot  k e X e  (wieghted)
                                                    //               L_tot = \sum_e l
                                                    // second fundamantal  form = a*S
    this->requireEdgeCotanWeights();
    this->requireFaceNormals();
    this->requireFaceTangentBasis();

    for (Face f : this->mesh.faces()) {
        Vector3 _curve_comp = get_curvature(f, 0);
        this->referenceCurvature[f][0] = _curve_comp[0];
        this->referenceCurvature[f][1] = _curve_comp[1];
        this->referenceCurvature[f][2] = _curve_comp[2];
    }
} 

void ElasticGeometry::computeActualCurvature() {
    this->requireEdgeCotanWeights();
    this->requireFaceNormals();
    this->requireFaceTangentBasis();

    for (Face f : this->mesh.faces()) {
        Vector3 _curve_comp = get_curvature(f, 1);
        this->referenceCurvature[f][0] = _curve_comp[0];
        this->referenceCurvature[f][1] = _curve_comp[1];
        this->referenceCurvature[f][2] = _curve_comp[2];
    }
}

void ElasticGeometry::computeElasticCauchyTensor() { ///
    if (!isElasticTensorInitializedF && !youngsModulus.toVector().isZero()) {
        // std::cout << "\n Executing Cauchy calc...\n";
        // float _Atensor[6];
        float _invmet[3];
        float _det;
        float _coef;
        for (Face f : this->mesh.faces()) {
            _det = referenceMetric[f][0] * referenceMetric[f][1] - referenceMetric[f][2] * referenceMetric[f][2];
            _invmet[0] =  referenceMetric[f][1] / _det;
            _invmet[1] =  referenceMetric[f][0] / _det;
            _invmet[2] = -referenceMetric[f][2] / _det;
            _coef = youngsModulus[f] / (1 - poissonsRatio[f] * poissonsRatio[f]) / 8;
            elasticCauchyTensor[f](0, 0) = _coef * _invmet[0] * _invmet[0];
            elasticCauchyTensor[f](0, 1) =
                _coef * (_invmet[2] * _invmet[2] * (1 - poissonsRatio[f]) + _invmet[0] * _invmet[1] * poissonsRatio[f]);
            elasticCauchyTensor[f](1, 0) = elasticCauchyTensor[f](0, 1);
            elasticCauchyTensor[f](0, 2) = _coef * 2 * _invmet[0] * _invmet[2];
            elasticCauchyTensor[f](2, 0) = .5 * elasticCauchyTensor[f](0, 2);
            elasticCauchyTensor[f](1, 1) = _coef * _invmet[1] * _invmet[1];
            elasticCauchyTensor[f](1, 2) = _coef * 2* _invmet[1] * _invmet[2];
            elasticCauchyTensor[f](2, 1) = 0.5 * elasticCauchyTensor[f](1, 2);
            elasticCauchyTensor[f](2, 2) = _coef * (_invmet[2] * _invmet[2] * (1 + poissonsRatio[f]) +
                                                    _invmet[0] * _invmet[1] * (1 - poissonsRatio[f]));
        }
        isElasticTensorInitializedF = true;
    }
}

void ElasticGeometry::computeThickness() {
    // Nothing to compute. Currently thickness is constant
}

void ElasticGeometry::computeYoungsModulus() {
    // Nothing to compute. Currently  is constant
}

void ElasticGeometry::computePoissonsRatio() {
    // Nothing to compute. Currently  is constant
}

void inline ElasticGeometry::computeElasticEnergy() {
    this->requireBendingEnergy();
    this->requireStretchingEnergy();
    if (elasticEnergy.size() == 0) elasticEnergy = FaceData<double>(this->mesh, 0);
    for (Face f : this->mesh.faces())
    {
        elasticEnergy[f] = thickness[f] * stretchingEnergy[f];
        //+1 / 3 * thickness[f] * thickness[f] * thickness[f] * bendingEnergy[f];
    }
}

void inline ElasticGeometry::computeStretchingEnergy() {
    this->requireFaceAreas();
    if (stretchingEnergy.size() == 0) stretchingEnergy = FaceData<double>(this->mesh, 0);
    for (Face f : this->mesh.faces()) {
        calculate_stretching_energy(f); // energy content  
        if (stretchingEnergy[f] < 0) {
            std::cout << "ERROR! Negative stretching energy!  at face: " << f.getIndex() << ".\n";
            std::cout << "\n  Reference lengths:  {";
            int edgecount = 0;
            for (Edge e : f.adjacentEdges()) {
                edgecount++;
                std::cout << referenceLengths[e];
                if (edgecount == 3)
                    std::cout << "}\n";
                else
                    std::cout << ",";
            }
            
            std::cout << "\n  Actual lengths:  {";
            edgecount = 0;
            for (Edge e : f.adjacentEdges()) {
                edgecount++;
                std::cout << edgeLengths[e];
                if (edgecount == 3)
                    std::cout << "}\n";
                else
                    std::cout << ",";
            }

            std::cout << "\n  Reference Metric: \n";
            std::cout << referenceMetric[f][0] << ", \t";
            std::cout << referenceMetric[f][1] << ", \t";
            std::cout << referenceMetric[f][2] << "\n";

            std::cout << "\n  Actual Metric: \n";
            std::cout << actualMetric[f][0] << ", \t";
            std::cout << actualMetric[f][1] << ", \t";
            std::cout << actualMetric[f][2] << "\n";


            std::cout << "\n  Cauchy Tensor:\n";
            std::cout <<elasticCauchyTensor[f](0, 0) << ", \t";
            std::cout <<elasticCauchyTensor[f](0, 1) << ", \t";
            std::cout <<elasticCauchyTensor[f](0, 2);
            std::cout << "\n";
            std::cout <<elasticCauchyTensor[f](1, 0) << ", \t";
            std::cout <<elasticCauchyTensor[f](1, 1) << ", \t";
            std::cout <<elasticCauchyTensor[f](1, 2);
            std::cout << "\n";
            std::cout << elasticCauchyTensor[f](2, 0) << ", \t";
            std::cout << elasticCauchyTensor[f](2, 1) << ", \t";
            std::cout << elasticCauchyTensor[f](2, 2) << "\n";

        }
    }
}

void ElasticGeometry::computeBendingEnergy() { // Currently not implemented. Here as a placeholder   benedicte wants to see something!
    if (bendingEnergy.size() == 0) bendingEnergy = FaceData<double>(this->mesh, 0);
}


void ElasticGeometry::computeGradient() {
    elasticGradient = VertexData<Vector3>(this->mesh, Vector3{0, 0, 0});
    double _epsilon = 1e-6;    
    for (Vertex v : this->mesh.vertices()) {
        for (int _direction = 0; _direction < 3; _direction++) {
            double _ePlus=0;
            double _eMinus=0;
            double _eOrig = 0;

           /* double _eInit = 0;
            double _eFinit = 0;
            Vector3 _posInit;
            Vector3 _posFinit;
            Eigen::Vector3f _metricInit;
            Eigen::Vector3f _metricFinit;*/  //Debuggers



           /*_posInit = vertexPositions[v];           
            for (Face f : v.adjacentFaces()) {
                _eInit += this->elasticEnergy[f];
                _metricInit = actualMetric[f];
            }*/
            
            vertexPositions[v][_direction] += _epsilon;
            //_posFinit = vertexPositions[v];
            updateLocalEnergy(v);
            for (Face f : v.adjacentFaces())
            {
                _ePlus += this->elasticEnergy[f];
                //_metricFinit = actualMetric[f];
            }

            vertexPositions[v][_direction] -= 2*_epsilon;
            //_posFinit = vertexPositions[v];
            updateLocalEnergy(v);
            for (Face f : v.adjacentFaces()) {
                _eMinus += this->elasticEnergy[f];
                //_metricFinit = actualMetric[f];
            }
              
            vertexPositions[v][_direction] += _epsilon;
            //_posFinit = vertexPositions[v];
            updateLocalEnergy(v);

             for (Face f : v.adjacentFaces()) {
                //_eFinit += this->elasticEnergy[f];
                //_metricFinit = actualMetric[f];
            }

            elasticGradient[v][_direction] +=  -(_ePlus - _eMinus) / 2 / _epsilon;
        }        
    }
}

// TO DO: Break down the COMPUTEGRADIENT routine to local calculation. To this end we need to update local edge lengths
// (easy), and then calculate local metric - > for this creat a new, small routine given a face (also implement in get
// actual and referenc metric fucntions)
// 
// Psudo Code:
// void update_local_energy (vertex v) {
// calculate_adjacent_edges_lenght(v);
// calculate_adjacent_faces_metric(v);
// calculate_adjacent_faces_curvature(v); 
// calculate_adjacent_faces_energy(v); 
// }
// 
// void calculate_adjacent_edges_lenght(vertex v)  {
// for (adjacent edge) calculatelength(edge);
// }
// 
// void calculate_adjacent_faces_metric(vertex v){
// for (adjacent face) calculatemetric(face,v);
// }
// 
// calculatemetric(face,v) - calculates the metric given a change in v, we need v to  compare if it is l_1 l_2 or l_3;
// similar expression for the curvature.
// finally, we need to calculate the local energy;
// 
// calculate_adjacent_faces_energy(v);  is just a simple run over the adjacent faces and calculating (no need to know vertex position with relation to face)
// 

void ElasticGeometry::updateLocalEnergy(const Vertex& v) {   
    calculate_adjacent_edges_lenght(v);
    calculate_adjacent_faces_metric(v);
    //calculate_adjacent_faces_curvature(v); // to implement
    calculate_adjacent_faces_energy(v); 
}

void ElasticGeometry::calculate_adjacent_edges_lenght(const Vertex& v) {
    Vector3 _edgeVec;
    for (Edge e : v.adjacentEdges()) {
        _edgeVec = this->vertexPositions[v] - this->vertexPositions[e.otherVertex(v)];
        this->edgeLengths[e] = _edgeVec.norm();
       /* std::cout << "\n  edge diff: " << this->edgeLengths[e] - _edgeVec.norm();*/
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
   this->actualMetric[f][0] = std::pow(_faceEdgesLengths(0), 2);
   this->actualMetric[f][1] = std::pow(_faceEdgesLengths(1), 2);
   this->actualMetric[f][2] = 0.5 * (std::pow(_faceEdgesLengths(0), 2) + std::pow(_faceEdgesLengths(1), 2) -
                                    std::pow(_faceEdgesLengths(2), 2));
    
}

void ElasticGeometry::calculate_adjacent_faces_curvature(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculate_curvature(f);
    }
}

void ElasticGeometry::calculate_curvature(const Face& f) { // To implement    
}


void ElasticGeometry::calculate_adjacent_faces_energy(const Vertex& v) {
    for (Face f : v.adjacentFaces()) {
        calculate_stretching_energy(f);
        calculate_bending_energy(f);
        elasticEnergy[f] =
            thickness[f] * stretchingEnergy[f] + 1 / 3 * thickness[f] * thickness[f] * thickness[f] * bendingEnergy[f];
    }
}

void ElasticGeometry::calculate_stretching_energy(const Face& f) {
    Eigen::Vector3f _metricDiff = actualMetric[f] - referenceMetric[f];
    stretchingEnergy[f] = elasticCauchyTensor[f](0, 0) * _metricDiff[0] * _metricDiff[0] +
                          elasticCauchyTensor[f](1, 1) * _metricDiff[1] * _metricDiff[1] +
                          2 * elasticCauchyTensor[f](2, 2) * _metricDiff[2] * _metricDiff[2] +
                          2 * elasticCauchyTensor[f](1, 0) * _metricDiff[0] * _metricDiff[1] +
                          4 * elasticCauchyTensor[f](2, 0) * _metricDiff[0] * _metricDiff[2] +
                          4 * elasticCauchyTensor[f](2, 1) * _metricDiff[1] * _metricDiff[2]; // energy content
    stretchingEnergy[f] *= faceAreas[f]; // 2D energy (not including thicness)
}


void ElasticGeometry::calculate_bending_energy(const Face& f) {
    bendingEnergy[f] *= 0; // 2D energy (not including thicness)
}

void ElasticGeometry::computePressure() {
    // Nothing to compute. Currently  is constant
}


void ElasticGeometry::computeRegions() {} // NOT YET IMPLEMENTED


void ElasticGeometry::computeFixedVertexs() {} // NOT YET IMPLEMENTED


void ElasticGeometry::computeFixedAngles() {} // NOT YET IMPLEMENTED






} // namespace surface
} // namespace geometrycentral