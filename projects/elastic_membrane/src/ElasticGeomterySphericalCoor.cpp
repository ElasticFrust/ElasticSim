
//#include "ElasticGeometry.cpp" // my class to be implemented
#include "ElasticGeomterySphericalCoor.h"
#include <fstream>
#include <limits>



using namespace geometrycentral;
using namespace geometrycentral::surface;


namespace geometrycentral {
namespace surface {
   
    ElasticGeometrySphericalCoor::ElasticGeometrySphericalCoor(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_)
                                : ElasticGeometrySphericalCoor::ElasticGeometrySphericalCoor(mesh_, inputVertexPositions_, 0, 0, 0, 0) {};
 
    //clang-format off
    ElasticGeometrySphericalCoor::ElasticGeometrySphericalCoor(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                                                           const double& THICKNESS_, const double& YOUNGs_, const double& POISSONs_, const double& PRESSURE_): 
        ElasticGeometry(mesh_, inputVertexPositions_, THICKNESS_, YOUNGs_, POISSONs_, PRESSURE_),
        faceCentroidPositionQ               (&faceCentroidPosition,         std::bind(&ElasticGeometrySphericalCoor::computeCentroids, this),                   quantities),
        vertexCoordinatesQ                  (&vertexCoordinates,            std::bind(&ElasticGeometrySphericalCoor::computeVertexCoordinates, this),           quantities),
        faceCentroidCoordinatesQ            (&faceCentroidCoordinates,      std::bind(&ElasticGeometrySphericalCoor::computeFaceCentroidCoordinates,this),      quantities)
    {

    faceCentroidPosition = FaceData<Vector3>(this->mesh, Vector3({0, 0, 0}));
    faceCentroidCoordinates = FaceData<Vector2>(this->mesh, Vector2({0, 0}));
    vertexCoordinates = VertexData<Vector2>(this->mesh, Vector2({0, 0}));
    requireFaceCentroidPosition();
    requireFaceCentroidCoordinates();
    requireVertexCoordinates();
    faceCentroidPositionQ.clearable = false;
    faceCentroidCoordinatesQ.clearable = false;
    vertexCoordinatesQ.clearable = false;
    isElasticTensorInitializedF = false;

     unrequireActualCurvature();
     unrequireActualMetric();
     unrequireReferenceCurvature();
     unrequireReferenceMetric();
     unrequireElasticCauchyTensor();
     actualMetricQ.clearable = true;
     referenceMetricQ.clearable = true;
     actualCurvatureQ.clearable = true;
     referenceCurvatureQ.clearable = true;
     elasticCauchyTensorQ.clearable = true;

     actualMetricQ.clearIfNotRequired();
     actualCurvatureQ.clearIfNotRequired();
     referenceCurvatureQ.clearIfNotRequired();
     referenceMetricQ.clearIfNotRequired();
     elasticCauchyTensorQ.clearIfNotRequired();

     
     referenceMetric = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(1., 2., 3.));
     requireReferenceMetric();

     referenceCurvature = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(1., 2., 3.));
     requireReferenceCurvature();

     actualMetric = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(0., 0., 0.));
     actualCurvature = FaceData<Eigen::Vector3f>(this->mesh, Eigen::Vector3f(0., 0., 0.));
     requireActualMetric();     
     requireActualCurvature();
     
     elasticCauchyTensor = FaceData<Eigen::Matrix3f>(this->mesh, Eigen::Matrix3f());
     requireElasticCauchyTensor();
     elasticCauchyTensorQ.clearable = false;

    }
    //clang-format on


    void ElasticGeometrySphericalCoor::requireVertexCoordinates() {
        vertexCoordinatesQ.require();
    };

    void ElasticGeometrySphericalCoor::unrequireVertexCoordinates() {
        vertexCoordinatesQ.unrequire();
    };


    void ElasticGeometrySphericalCoor::requireFaceCentroidCoordinates() {
        faceCentroidCoordinatesQ.require();
    };
    void ElasticGeometrySphericalCoor::unrequireFaceCentroidCoordinates() {
        faceCentroidCoordinatesQ.unrequire();
    };


    void ElasticGeometrySphericalCoor::requireFaceCentroidPosition() {
        faceCentroidPositionQ.require();
    };
    void ElasticGeometrySphericalCoor::unrequireFaceCentroidPosition() {
        faceCentroidPositionQ.unrequire();
    };

    /// <summary>
    ///
    /// </summary>
    void ElasticGeometrySphericalCoor::computeReferenceMetric(){
        computeCentroids();
        if (!isReferenceMetricInitializedF) {
            referenceLengthsQ.ensureHave();
            for (Face f : mesh.faces()) {
                calculateFaceReferenceMetric(f);
                if (referenceMetric[f][0] < 0) {
                    double ls[3] = {0, 0, 0};
                    int iter = 0;
                    for (Edge e : f.adjacentEdges()) {
                        ls[iter] = edgeLengths[e];
                        iter++;
                    }
                    std::cout << "\nNegative Metric!\nFace: " << f.getIndex()
                              << "\nNeighboring Edges Lenghts: " << ls[0] << ", " << ls[1] << ", " << ls[2]
                              << "\nCentroid Coordinates: (" << faceCentroidCoordinates[f][0] << ", "
                              << faceCentroidCoordinates[f][1] << ")\n";
                }
            }            
        }
        isReferenceMetricInitializedF = true;
    };


    void ElasticGeometrySphericalCoor::computeActualMetric()  {
        computeCentroids();
        edgeLengthsQ.ensureHave();
        for (Face f : mesh.faces()) {
            calculateFaceActualMetric(f);
            calculateFaceReferenceMetric(f);
            if (actualMetric[f][0] < 0 || actualMetric[f][1] < 0) {
                double ls[3] = {0, 0, 0};
                int iter = 0;
                for (Edge e : f.adjacentEdges()) {
                    ls[iter] = edgeLengths[e];
                    iter++;
                }
                std::cout << "\nNegative Metric!\nFace: " << f.getIndex() << "\nNeighboring Edges Lenghts: " << ls[0]
                          << ", " << ls[1] << ", " << ls[2] << "\nCentroid Coordinates: ("
                          << faceCentroidCoordinates[f][0] << ", " << faceCentroidCoordinates[f][1] << ")\n";
            }
        }
    };

    void ElasticGeometrySphericalCoor::computeReferenceCurvature() {
        if (!isReferenceCurvatureInitializedF) {
            faceNormalsQ.ensureHave();
            for (Face f : mesh.faces()) {
                calcualteFaceReferenceCurvature(f);
            }
            isReferenceCurvatureInitializedF = true;
        }
    };

    void ElasticGeometrySphericalCoor::computeActualCurvature() {
        faceNormalsQ.ensureHave();
        for (Face f : mesh.faces()) {
            calcualteFaceActualCurvature(f);
        }
    };

    //void ElasticGeometrySphericalCoor::computeElasticEnergy()  {};

    void ElasticGeometrySphericalCoor::computeGradient() {
        vertexDualAreasQ.ensureHave();
        vertexNormalsQ.ensureHave();
        totalEnergyQ.ensureHave();
        elasticEnergyQ.ensureHave();
        computeCentroids();
        elasticGradient = VertexData<Vector3>(this->mesh, Vector3{0., 0., 0.});
        double _epsilon = 1.e-6;
        for (Vertex v : this->mesh.vertices()) {
            for (int _direction = 0; _direction < 3; _direction++) {
                double _ePlus = 0.;
                double _eMinus = 0.;
                double _eOrig = 0.;


                vertexPositions[v][_direction] += _epsilon;
                localEnergyChange(v);
                for (Face f : v.adjacentFaces()) {
                    _ePlus += this->totalEnergy[f];
                }

                vertexPositions[v][_direction] -= 2. * _epsilon;
                localEnergyChange(v);
                for (Face f : v.adjacentFaces()) {
                    _eMinus += this->totalEnergy[f];
                }

                vertexPositions[v][_direction] += _epsilon;
                
                localEnergyChange(v);

                elasticGradient[v][_direction] += -(_ePlus - _eMinus) / 2. / _epsilon;
                
            }
        }

     refreshQuantities();
    };

    void ElasticGeometrySphericalCoor::updateElasticCauchyTensor() {
        isElasticTensorInitializedF = false;
        elasticCauchyTensorQ.clearable = true;
        elasticCauchyTensorQ.clearIfNotRequired();
        elasticCauchyTensor = FaceData<Eigen::Matrix3f>(this->mesh, Eigen::Matrix3f());
        requireElasticCauchyTensor();
        elasticCauchyTensorQ.clearable = false;
    }

    /// <summary>
    /// same as ElasticGeometry::updateLocalEnergy(const Vertex& v)  up to difference in geomtery
    /// </summary>
    /// <param name="v"></param>
    void ElasticGeometrySphericalCoor::localEnergyChange(const Vertex v) {
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
        calculateAdjacentNormalAndCentroids(v);
        //calculate_adjacent_edges_dihedral_angles(v); // THIS ALSO UPDATES NORMAL!
        calculateAdjacentMetric(v);
        calculateAdjacentCurvature(v);
        calculate_adjacent_faces_energy(v);
        calculate_adjacent_faces_volume(v);
        calculate_adjacent_faces_total_energy(v);
    };


    void ElasticGeometrySphericalCoor::calculateAdjacentNormalAndCentroids(const Vertex v) {
        for (Face f : v.adjacentFaces()) {
            // For general polygons, take the sum of the cross products at each corner
            Vector3 normalSum = Vector3::zero();
            Vector3 centroid = Vector3::zero();
            for (Halfedge heF : f.adjacentHalfedges()) {

                // Gather vertex positions for next three vertices
                Halfedge he = heF;
                Vector3 pA = vertexPositions[he.vertex()];
                he = he.next();
                Vector3 pB = vertexPositions[he.vertex()];
                he = he.next();
                Vector3 pC = vertexPositions[he.vertex()];

                normalSum += cross(pB - pA, pC - pA);
                centroid += (pA + pB + pC) / 3.; /// sketchy (wrong) line if not triangle 

                // In the special case of a triangle, there is no need to to repeat at all three corners; the result
                // will be the same
                if (he.next() == heF) break;
            }

            Vector3 normal = unit(normalSum);
            faceNormals[f] = normal;
            faceCentroidPosition[f] = Vector3(centroid);
        }
    };

    void ElasticGeometrySphericalCoor::calculateAdjacentMetric(const Vertex v) {
        for (Face f : v.adjacentFaces()) {
            calculateFaceActualMetric(f);
        }
    };

    void ElasticGeometrySphericalCoor::calculateAdjacentCurvature(const Vertex v) {
        for (Face f : v.adjacentFaces()) {
            calcualteFaceActualCurvature(f);
        }
    };

    void ElasticGeometrySphericalCoor::calculateFaceActualMetric(Face f) {
        edgeLengthsQ.ensureHave();
        float lengths2[3] = {0, 0, 0};
        int iter = 0;
        for (Edge e : f.adjacentEdges()) {
            lengths2[iter] = (float) std::pow(edgeLengths[e], 2);
            iter += 1;
        }
        actualMetric[f] = calculateMetric(f, lengths2);
    }


    void ElasticGeometrySphericalCoor::calculateFaceReferenceMetric(Face f) {
        if (!isReferenceMetricInitializedF) {
            referenceLengthsQ.ensureHave();
            float lengths2[3] = {0, 0, 0};
            int iter = 0;
            for (Edge e : f.adjacentEdges()) {
                lengths2[iter] = (float)std::pow(referenceLengths[e], 2);
                iter += 1;
            }
            referenceMetric[f] = calculateMetric(f, lengths2);
        }
    };

    void ElasticGeometrySphericalCoor::calcualteFaceActualCurvature(Face f) {
        faceNormalsQ.ensureHave();
        actualCurvature[f] = calculateCurvature(f);
    };
    void ElasticGeometrySphericalCoor::calcualteFaceReferenceCurvature(Face f) {
        faceNormalsQ.ensureHave();
        referenceCurvature[f] = calculateCurvature(f);
    };



    Eigen::Vector3f ElasticGeometrySphericalCoor::calculateMetric(Face f, float lengs2[3]) {   
        float lengths2[3] = {lengs2[0], lengs2[1], lengs2[2]};
        Vector2 coor_diffs[3] = {{0, 0}, {0, 0}, {0, 0}};
        int iter = 0;
        int singularV = -1;
        Vertex v1, v2;
        // NEW way - using distances from centroid:
        float tolerance = 5e-2;

        if (false && f.getIndex() == 3871) {
            std::cout << "\nFace 3871 Centroid coordinates: (" << faceCentroidCoordinates[f][0] << ", "
                      << faceCentroidCoordinates[f][1] << ")\n";
            std::cout << "\nFace 3871 Centroid position: (" << faceCentroidPosition[f][0] << ", "
                      << faceCentroidPosition[f][1] << ", " << faceCentroidPosition[f][2] << ")\n";
        }
        
        for (Vertex v : f.adjacentVertices()) {                
            coor_diffs[iter] = calculateCoordinateDiff(faceCentroidCoordinates[f], vertexCoordinates[v]);
            if (vertexCoordinates[v][0] < tolerance)
                if (abs(vertexCoordinates[v][1] - faceCentroidCoordinates[f][1])>PI/3) /// NASTY PATCH need a better criteria (basically if 0,0 is IN this  face)
                    coor_diffs[iter] = Vector2({-vertexCoordinates[v][0] - faceCentroidCoordinates[f][0], 0});
                else
                    coor_diffs[iter] = Vector2({vertexCoordinates[v][0] - faceCentroidCoordinates[f][0], 0});
            if (PI - vertexCoordinates[v][0] < tolerance)
                if (abs(vertexCoordinates[v][1] - faceCentroidCoordinates[f][1]) > PI / 3)
                    coor_diffs[iter] = Vector2({2 * PI - vertexCoordinates[v][0] - faceCentroidCoordinates[f][0], 0});
                else 
                    coor_diffs[iter] = Vector2({ vertexCoordinates[v][0] - faceCentroidCoordinates[f][0], 0});    
            lengths2[iter] =
                (float)dot(faceCentroidPosition[f] - vertexPositions[v], faceCentroidPosition[f] - vertexPositions[v]);
            iter += 1;
            if (false && f.getIndex() == 3871) {
                std::cout << "\nvertex " << iter << "coordinates: (" << vertexCoordinates[v][0] << ", "
                          << vertexCoordinates[v][1] << ")\n";
                std::cout << "Vertex position: (" << vertexPositions[v][0] << ", " << vertexPositions[v][1] << ", "
                          << vertexPositions[v][2] << ")\n";
            }

        }

        // OLD way  using triangle edges. 
        //for (Vertex v : f.adjacentVertices()) {
        //    if (false && (vertexCoordinates[v][0] < 1e-3 || PI - vertexCoordinates[v][0] < 1e-3)) {                
        //        singularV = v.getIndex();
        //    }
        //}
        //for (Edge e : f.adjacentEdges()) {
        //    v1 = e.firstVertex();
        //    v2 = e.secondVertex();
        //    if (v1.getIndex() == singularV) {
        //        coor_diffs[iter] = calculateCoordinateDiff(faceCentroidCoordinates[f], vertexCoordinates[v2]);
        //        lengths2[iter] = (float)dot(faceCentroidPosition[f] - vertexPositions[v2],
        //                                    faceCentroidPosition[f] - vertexPositions[v2]);
        //    } else if (v2.getIndex() == singularV) {
        //        coor_diffs[iter] = calculateCoordinateDiff(vertexCoordinates[v1], faceCentroidCoordinates[f]);
        //        lengths2[iter] = (float)dot(- faceCentroidPosition[f] + vertexPositions[v1],
        //                                    - faceCentroidPosition[f] + vertexPositions[v1]);
        //    } else { //none is singular or close to it
        //        coor_diffs[iter] = calculateCoordinateDiff(vertexCoordinates[v1], vertexCoordinates[v2]);
        //    }
        //    iter += 1;
        //    if (false &&  f.getIndex() == 32) {
        //        if (v1.getIndex() == singularV) {
        //            std::cout << "\nvertex " << iter << " is pole, centroid coordinates: ("
        //                      << faceCentroidCoordinates[f][0] << ", " << faceCentroidCoordinates[f][1] << ")\n";
        //            std::cout << "\nvertex " << iter << " is pole, centroid position: (" << faceCentroidPosition[f].x
        //                      << ", " << faceCentroidPosition[f].y << ", " << faceCentroidPosition[f].z
        //                      << ")\n";
        //            std::cout << "\nvertex " << iter << " is pole, original coordinates: ("
        //                      << vertexCoordinates[e.firstVertex()][0] << ", " << vertexCoordinates[e.firstVertex()][1]
        //                      << ")\n";
        //            std::cout << "\nvertex " << iter << " is pole, original position: ("
        //                      << vertexPositions[e.firstVertex()][0] << ", " << vertexPositions[e.firstVertex()][1]
        //                      << ", " << vertexPositions[e.firstVertex()][2] << ")\n";

        //        } else {
        //            std::cout << "\nvertex " << iter << " coordinates: (" << vertexCoordinates[e.firstVertex()][0]
        //                      << ", " << vertexCoordinates[e.firstVertex()][1] << ")\n";
        //            std::cout << "\nvertex " << iter << "  position: (" << vertexPositions[e.firstVertex()][0] << ", "
        //                      << vertexPositions[e.firstVertex()][1] << ", " << vertexPositions[e.firstVertex()][2]
        //                      << ")\n";
        //        }
        //    }
        //}
        double ddet = (coor_diffs[0][1] * coor_diffs[1][0] - coor_diffs[0][0] * coor_diffs[1][1]) *
                      (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) *
                      (coor_diffs[2][1] * coor_diffs[1][0] - coor_diffs[2][0] * coor_diffs[1][1]);

       double res[3] = {
            (lengths2[2] * coor_diffs[0][1] * coor_diffs[1][1] *
                        (-coor_diffs[0][1] * coor_diffs[1][0] + coor_diffs[0][0] * coor_diffs[1][1]) +
                    lengths2[1] * coor_diffs[0][1] * coor_diffs[2][1] *
                        (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) +
                    lengths2[0] * coor_diffs[1][1] * coor_diffs[2][1] *
                        (-coor_diffs[1][1] * coor_diffs[2][0] + coor_diffs[1][0] * coor_diffs[2][1]))/ddet,

               (lengths2[2] * coor_diffs[0][0] * coor_diffs[1][0] *
                        (-coor_diffs[0][1] * coor_diffs[1][0] + coor_diffs[0][0] * coor_diffs[1][1]) +
                    lengths2[1] * coor_diffs[0][0] * coor_diffs[2][0] *
                        (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) +
                    lengths2[0] * coor_diffs[1][0] * coor_diffs[2][0] *
                        (-coor_diffs[1][1] * coor_diffs[2][0] + coor_diffs[1][0] * coor_diffs[2][1]))/ddet,

               (lengths2[2] * (coor_diffs[0][1] * coor_diffs[1][0] + coor_diffs[0][0] * coor_diffs[1][1]) *
                        (coor_diffs[0][1] * coor_diffs[1][0] - coor_diffs[0][0] * coor_diffs[1][1]) +
                    -lengths2[1] * (coor_diffs[0][1] * coor_diffs[2][0] + coor_diffs[0][0] * coor_diffs[2][1]) *
                        (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) +
                    lengths2[0] * (coor_diffs[1][1] * coor_diffs[2][0] + coor_diffs[1][0] * coor_diffs[2][1]) *
                        (coor_diffs[1][1] * coor_diffs[2][0] - coor_diffs[1][0] * coor_diffs[2][1]))/ddet};


       if (false && f.getIndex() == 3871) {
           
           std::cout << "\nFace 3871 edge 0 length^2: " << lengths2[0] << " coordinate difference: (" << coor_diffs[0][0]<< ", " << coor_diffs[0][1] << ")\n";
           std::cout << "\nFace 3871 edge 1 length^2: " << lengths2[1] << " coordinate difference: (" << coor_diffs[1][0]<< ", " << coor_diffs[1][1] << ")\n";
           std::cout << "\nFace 3871 edge 2 length^2: " << lengths2[2] << " coordinate difference: (" << coor_diffs[2][0]<< ", " << coor_diffs[2][1] << ")\n";
           std::cout << "";
       }


        return Eigen::Vector3f((float) res[0], (float) res[1], (float) res[2]);
    };

    Eigen::Vector3f ElasticGeometrySphericalCoor::calculateCurvature(Face f) {
        double angles[3] = {0, 0, 0};
        Vector2 coor_diffs[3] = {{0, 0}, {0, 0}, {0, 0}};
        Vector3 dr[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        int iter = 0;
        int singularV = -1;
        Vertex v1, v2;
        for (Vertex v : f.adjacentVertices()) {
            if (vertexCoordinates[v][0] < 5e-2 || PI - vertexCoordinates[v][0] < 5e-2) {
                singularV = v.getIndex();
            }
        }
        for (Face af : f.adjacentFaces()) {
            coor_diffs[iter] = calculateCoordinateDiff(faceCentroidCoordinates[f],faceCentroidCoordinates[af]);
            if ( true && singularV > -1 ) {
                coor_diffs[iter] += Vector2({1e-3, 1e-3});
            }

            
            dr[iter] = faceCentroidPosition[af] - faceCentroidPosition[f] ;
            angles[iter] = -dot(faceNormals[f], dr[iter]); // a matter of convention -is 
            iter += 1;
        }
        double ddet = (coor_diffs[0][1] * coor_diffs[1][0] - coor_diffs[0][0] * coor_diffs[1][1]) *
                      (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) *
                      (coor_diffs[2][1] * coor_diffs[1][0] - coor_diffs[2][0] * coor_diffs[1][1]);

         double res[3] = {
               (angles[2] * coor_diffs[0][1] * coor_diffs[1][1] *
                        (-coor_diffs[0][1] * coor_diffs[1][0] + coor_diffs[0][0] * coor_diffs[1][1]) +
                    angles[1] * coor_diffs[0][1] * coor_diffs[2][1] *
                        (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) +
                    angles[0] * coor_diffs[1][1] * coor_diffs[2][1] *
                        (-coor_diffs[1][1] * coor_diffs[2][0] + coor_diffs[1][0] * coor_diffs[2][1]))/ddet,

               (angles[2] * coor_diffs[0][0] * coor_diffs[1][0] *
                        (-coor_diffs[0][1] * coor_diffs[1][0] + coor_diffs[0][0] * coor_diffs[1][1]) +
                    angles[1] * coor_diffs[0][0] * coor_diffs[2][0] *
                        (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) +
                    angles[0] * coor_diffs[1][0] * coor_diffs[2][0] *
                        (-coor_diffs[1][1] * coor_diffs[2][0] + coor_diffs[1][0] * coor_diffs[2][1]))/ddet,

               (angles[2] * (coor_diffs[0][1] * coor_diffs[1][0] + coor_diffs[0][0] * coor_diffs[1][1]) *
                        (coor_diffs[0][1] * coor_diffs[1][0] - coor_diffs[0][0] * coor_diffs[1][1]) +
                    -angles[1] * (coor_diffs[0][1] * coor_diffs[2][0] + coor_diffs[0][0] * coor_diffs[2][1]) *
                        (coor_diffs[0][1] * coor_diffs[2][0] - coor_diffs[0][0] * coor_diffs[2][1]) +
                    angles[0] * (coor_diffs[1][1] * coor_diffs[2][0] + coor_diffs[1][0] * coor_diffs[2][1]) *
                        (coor_diffs[1][1] * coor_diffs[2][0] - coor_diffs[1][0] * coor_diffs[2][1]))/ddet};



         if (false && (f.getIndex() == 3871)) {
             std::cout << "\n Face No: " << f.getIndex() << "\nCoordinate difference 0: (" << coor_diffs[0][0] << ", "
                       << coor_diffs[0][1] << ")\nCoordinate difference 1: (" << coor_diffs[1][0] << ", "
                       << coor_diffs[1][1] << ")\nCoordinate difference 2: (" << coor_diffs[2][0] << ", "
                       << coor_diffs[2][1] << ")\n dr 0: (" << dr[0][0] << ", " << dr[0][1] << ", " << dr[0][2]
                       << ")\n dr 1: (" << dr[1][0] << ", " << dr[1][1] << ", " << dr[1][2] << ")\n dr 2: (" << dr[2][0]
                       << ", " << dr[2][1] << ", " << dr[2][2] << ")\nAngles: (" << angles[0] << ", " << angles[1]
                       << ", " << angles[2] << ")\n";
             std::cout << "\n Curvature: (" << res[0] << ", " << res[1] << ", " << res[2] << ")\n\n";
         }

        return Eigen::Vector3f((float)res[0], (float)res[1], (float)res[2]);
    };



    /// <summary>
    /// compute face centroid porisition, shoulde be done whenver
    /// </summary>
    void ElasticGeometrySphericalCoor::computeCentroids() {
        Vector3 tpos;
        for (Face f : mesh.faces()) {
            tpos = {0, 0, 0};
            for (Vertex v : f.adjacentVertices()) {
                tpos += vertexPositions[v] / 3.;
            }
            faceCentroidPosition[f].x = tpos.x;
            faceCentroidPosition[f].y = tpos.y;
            faceCentroidPosition[f].z = tpos.z;
            if (false && f.getIndex() == 32) {
                std::cout << "\nCentroid Pos: " << tpos[0] << ", " << tpos[1] << ", " << tpos[2] << "\n";
            }
        }
    }

    /// <summary>
    /// Coordinates are stereoscopic projections on the unit sphere. {1,0} is the \theta (polar) direction {0,1} is the
    /// \phi (azimuthal) direction. Calculated once at the begining, or for an additional face/vertex after refininmernt
    /// (not implemented)
    /// </summary>
    void ElasticGeometrySphericalCoor::computeVertexCoordinates() { // add test to see if already have coordinates
        if (isReferenceMetricInitializedF) return;
        vertexPositionsQ.ensureHave();
        Vector3 tempPos = {0, 0, 0};
        double rho = 0;

        // vertex coordinates first
        for (Vertex v : mesh.vertices()) {
            tempPos = vertexPositions[v];
            rho = std::sqrt(tempPos.x * tempPos.x + tempPos.z * tempPos.z);
            vertexCoordinates[v] = Vector2({atan2(rho, tempPos.y), atan2(tempPos.z, tempPos.x)});
        }
    };


    void ElasticGeometrySphericalCoor::computeFaceCentroidCoordinates() { // add test to see if already have coordinates
        requireFaceCentroidPosition();
        if (isReferenceMetricInitializedF) return;
        faceCentroidPositionQ.ensureHave();
        Vector3 tempPos = {0, 0, 0};
        double rho = 0;

        // the face centroid coordinates
        for (Face f : mesh.faces()) {
            tempPos = faceCentroidPosition[f];
            rho = std::sqrt(tempPos.x * tempPos.x + tempPos.z * tempPos.z);
            faceCentroidCoordinates[f] = Vector2({atan2(rho, tempPos.y), atan2(tempPos.z, tempPos.x)});
        }
    };


    void ElasticGeometrySphericalCoor::updateFaceCentroidCoordinates(const FaceData<Vector2> faceCoordinates) {
        faceCentroidCoordinatesQ.ensureHave();
        for (Face f : mesh.faces()) {            
            faceCentroidCoordinates[f] = Vector2(faceCoordinates[f]);
        }        
        computeActualMetric();
        computeActualCurvature();
    }


    Vector2 ElasticGeometrySphericalCoor::calculateCoordinateDiff(Vector2 p1, Vector2 p2) { 
        Vector2 res = p2 - p1;
        if (true && (p1[0] < 1e-3 || p2[0] < 1e-3 || PI - p1[0] < 1e-3 || PI - p2[0]<1e-3)) {
            if (abs(res[1]) > abs(res[1] - 1 * PI)) res[1] = res[1] - 1 * PI;
            else if (abs(res[1]) > abs(res[1] + 1 * PI))  res[1] = res[1] + 1* PI;
            return res;
        }
        else if (abs(res[1]) > abs(res[1] - 2 * PI)) res[1] = res[1] - 2 * PI;
        else if (abs(res[1]) > abs(res[1] + 2 * PI)) res[1] = res[1] + 2 * PI;
        return res + 0*Vector2({1e-3, 1e-3});

    };

    
} // namespace surface
} // namespace geometrycentral