/**
 * @file ElasticGeomterySphericalCoor.cpp
 * @brief Implementation of the spherical-coordinate elastic shell geometry.
 *
 * Provides metric and curvature computation using (theta, phi) coordinate
 * differences on the unit sphere, overriding the base class edge-based
 * approach. Includes central finite-difference gradient evaluation with
 * local re-evaluation around each perturbed vertex.
 */
#include "ElasticGeomterySphericalCoor.h"
#include <fstream>
#include <limits>



using namespace geometrycentral;
using namespace geometrycentral::surface;


namespace geometrycentral {
namespace surface {
   
    /** @brief Minimal constructor — delegates to the fully specified constructor with zero material constants. */
    ElasticGeometrySphericalCoor::ElasticGeometrySphericalCoor(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_)
                                : ElasticGeometrySphericalCoor::ElasticGeometrySphericalCoor(mesh_, inputVertexPositions_, 0, 0, 0, 0) {};
 
    /**
     * @brief Fully specified constructor.
     *
     * Initializes the base ElasticGeometry, registers three spherical-coordinate
     * dependency quantities (centroids, vertex coordinates, face centroid
     * coordinates), allocates all data arrays, unrequires/clears inherited
     * quantities so they can be re-evaluated in spherical form, then triggers
     * initial computation of reference and actual metrics, curvatures, and the
     * elastic Cauchy tensor.
     */
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


    // ========== Require / Unrequire ==========

    /** @brief Registers a dependency on vertexCoordinates. */
    void ElasticGeometrySphericalCoor::requireVertexCoordinates() {
        vertexCoordinatesQ.require();
    };
    /** @brief Releases a dependency on vertexCoordinates. */
    void ElasticGeometrySphericalCoor::unrequireVertexCoordinates() {
        vertexCoordinatesQ.unrequire();
    };

    /** @brief Registers a dependency on faceCentroidCoordinates. */
    void ElasticGeometrySphericalCoor::requireFaceCentroidCoordinates() {
        faceCentroidCoordinatesQ.require();
    };
    /** @brief Releases a dependency on faceCentroidCoordinates. */
    void ElasticGeometrySphericalCoor::unrequireFaceCentroidCoordinates() {
        faceCentroidCoordinatesQ.unrequire();
    };

    /** @brief Registers a dependency on faceCentroidPosition. */
    void ElasticGeometrySphericalCoor::requireFaceCentroidPosition() {
        faceCentroidPositionQ.require();
    };
    /** @brief Releases a dependency on faceCentroidPosition. */
    void ElasticGeometrySphericalCoor::unrequireFaceCentroidPosition() {
        faceCentroidPositionQ.unrequire();
    };

    // ========== Overridden Metric / Curvature Callbacks ==========

    /**
     * @brief Computes the reference metric tensor for all faces.
     *
     * On first call, iterates over all faces, computes the metric from
     * reference edge lengths and spherical coordinate differences via
     * calculateFaceReferenceMetric(). Warns on negative metric components.
     * Subsequent calls are no-ops (guarded by isReferenceMetricInitializedF).
     */
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


    /**
     * @brief Computes the actual metric tensor for all faces.
     *
     * Updates centroids, then computes both actual and reference metrics for
     * every face. Warns on negative metric components (degenerate triangles).
     */
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

    /**
     * @brief Computes the reference curvature tensor for all faces (one-time).
     *
     * On first call, evaluates curvature from face normals and centroid
     * coordinate differences via calcualteFaceReferenceCurvature().
     * Guarded by isReferenceCurvatureInitializedF.
     */
    void ElasticGeometrySphericalCoor::computeReferenceCurvature() {
        if (!isReferenceCurvatureInitializedF) {
            faceNormalsQ.ensureHave();
            for (Face f : mesh.faces()) {
                calcualteFaceReferenceCurvature(f);
            }
            isReferenceCurvatureInitializedF = true;
        }
    };

    /** @brief Computes the actual curvature tensor for all faces from current normals and centroids. */
    void ElasticGeometrySphericalCoor::computeActualCurvature() {
        faceNormalsQ.ensureHave();
        for (Face f : mesh.faces()) {
            calcualteFaceActualCurvature(f);
        }
    };

    /**
     * @brief Computes the elastic gradient via central finite differences.
     *
     * For each vertex and each Cartesian direction (x, y, z), perturbs the
     * vertex position by +/- epsilon, locally re-evaluates energy for
     * adjacent faces via localEnergyChange(), and approximates the gradient
     * component as -(E+ - E-) / (2 * epsilon). Calls refreshQuantities()
     * at the end to restore globally consistent state.
     */
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

    /**
     * @brief Forces re-computation of the elastic Cauchy tensor.
     *
     * Resets the initialization flag, temporarily makes the quantity clearable,
     * clears it, re-allocates, re-requires it, and locks it again.
     */
    void ElasticGeometrySphericalCoor::updateElasticCauchyTensor() {
        isElasticTensorInitializedF = false;
        elasticCauchyTensorQ.clearable = true;
        elasticCauchyTensorQ.clearIfNotRequired();
        elasticCauchyTensor = FaceData<Eigen::Matrix3f>(this->mesh, Eigen::Matrix3f());
        requireElasticCauchyTensor();
        elasticCauchyTensorQ.clearable = false;
    }

    /**
     * @brief Re-evaluates all energy-related quantities for faces adjacent to vertex v.
     *
     * Ensures all dependencies are available, then locally recomputes edge
     * lengths, face areas, normals, centroids, metrics, curvatures, energies,
     * volumes, and total energies for the 1-ring of v. Used during finite-
     * difference gradient evaluation.
     */
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
        calculateAdjacentMetric(v);
        calculateAdjacentCurvature(v);
        calculate_adjacent_faces_energy(v);
        calculate_adjacent_faces_volume(v);
        calculate_adjacent_faces_total_energy(v);
    };


    /**
     * @brief Recomputes face normals and centroid positions for all faces adjacent to vertex v.
     *
     * For each adjacent face, accumulates cross products of edge vectors to
     * compute the face normal, and averages vertex positions for the centroid.
     * Short-circuits after the first halfedge iteration for triangles.
     */
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

    /** @brief Recomputes the actual metric for all faces adjacent to vertex v. */
    void ElasticGeometrySphericalCoor::calculateAdjacentMetric(const Vertex v) {
        for (Face f : v.adjacentFaces()) {
            calculateFaceActualMetric(f);
        }
    };

    /** @brief Recomputes the actual curvature for all faces adjacent to vertex v. */
    void ElasticGeometrySphericalCoor::calculateAdjacentCurvature(const Vertex v) {
        for (Face f : v.adjacentFaces()) {
            calcualteFaceActualCurvature(f);
        }
    };

    /**
     * @brief Computes the actual metric tensor for a single face.
     *
     * Squares the current edge lengths and passes them to calculateMetric()
     * which solves for the metric using spherical coordinate differences.
     */
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


    /**
     * @brief Computes the reference metric tensor for a single face (one-time).
     *
     * Squares the reference edge lengths and passes them to calculateMetric().
     * Guarded by isReferenceMetricInitializedF.
     */
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

    /** @brief Computes the actual curvature tensor for a single face from normals and centroids. */
    void ElasticGeometrySphericalCoor::calcualteFaceActualCurvature(Face f) {
        faceNormalsQ.ensureHave();
        actualCurvature[f] = calculateCurvature(f);
    };
    /** @brief Computes the reference curvature tensor for a single face from normals and centroids. */
    void ElasticGeometrySphericalCoor::calcualteFaceReferenceCurvature(Face f) {
        faceNormalsQ.ensureHave();
        referenceCurvature[f] = calculateCurvature(f);
    };



    // ========== Core Metric / Curvature Computation ==========

    /**
     * @brief Computes the metric tensor for face f from squared edge lengths and spherical coordinates.
     *
     * Uses centroid-to-vertex coordinate differences in spherical (theta, phi) space.
     * Applies pole-proximity patches when vertices are near theta=0 or theta=pi
     * to handle coordinate singularities. Replaces the input squared edge lengths
     * with squared Euclidean centroid-to-vertex distances and solves for the three
     * independent metric components (a11, a22, a12) via a Cramer-like formula.
     *
     * @param f     The face.
     * @param lengs2 Squared edge lengths (overwritten internally with centroid distances).
     * @return Metric tensor as (a11, a22, a12).
     */
    Eigen::Vector3f ElasticGeometrySphericalCoor::calculateMetric(Face f, float lengs2[3]) {   
        float lengths2[3] = {lengs2[0], lengs2[1], lengs2[2]};
        Vector2 coor_diffs[3] = {{0, 0}, {0, 0}, {0, 0}};
        int iter = 0;
        float tolerance = 5e-2;

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
        }

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


        return Eigen::Vector3f((float) res[0], (float) res[1], (float) res[2]);
    };

    /**
     * @brief Computes the curvature tensor for face f from centroid coordinate differences and normals.
     *
     * Uses coordinate differences between the face centroid and its adjacent face
     * centroids. Projects the displacement vectors onto the face normal to obtain
     * the "angles" (normal height differences). Near-pole faces receive a small
     * perturbation to avoid singular coordinate differences. Solves for the three
     * independent curvature components (b11, b22, b12) via a Cramer-like formula.
     *
     * @param f The face.
     * @return Curvature tensor as (b11, b22, b12).
     */
    Eigen::Vector3f ElasticGeometrySphericalCoor::calculateCurvature(Face f) {
        double angles[3] = {0, 0, 0};
        Vector2 coor_diffs[3] = {{0, 0}, {0, 0}, {0, 0}};
        Vector3 dr[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        int iter = 0;
        int singularV = -1;
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



        return Eigen::Vector3f((float)res[0], (float)res[1], (float)res[2]);
    };



    // ========== Coordinate Computation ==========

    /**
     * @brief Computes face centroid 3D positions by averaging the three vertex positions.
     *
     * Iterates over all faces and sets faceCentroidPosition[f] to the mean
     * of the adjacent vertex positions (assumes triangular faces).
     */
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
        }
    }

    /**
     * @brief Converts vertex 3D positions to spherical (theta, phi) coordinates.
     *
     * theta = atan2(rho, y) where rho = sqrt(x^2 + z^2), phi = atan2(z, x).
     * Only runs once — returns immediately if the reference metric is already initialized.
     */
    void ElasticGeometrySphericalCoor::computeVertexCoordinates() {
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


    /**
     * @brief Converts face centroid 3D positions to spherical (theta, phi) coordinates.
     *
     * Same formula as computeVertexCoordinates() but applied to face centroids.
     * Only runs once — returns immediately if the reference metric is already initialized.
     */
    void ElasticGeometrySphericalCoor::computeFaceCentroidCoordinates() {
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


    /**
     * @brief Overwrites face centroid coordinates from external data and recomputes derived quantities.
     *
     * Copies the provided faceCoordinates into faceCentroidCoordinates, then
     * triggers recomputation of the actual metric and actual curvature.
     */
    void ElasticGeometrySphericalCoor::updateFaceCentroidCoordinates(const FaceData<Vector2> faceCoordinates) {
        faceCentroidCoordinatesQ.ensureHave();
        for (Face f : mesh.faces()) {            
            faceCentroidCoordinates[f] = Vector2(faceCoordinates[f]);
        }        
        computeActualMetric();
        computeActualCurvature();
    }


    /**
     * @brief Computes the shortest-path coordinate difference between two spherical points.
     *
     * Handles phi-periodicity: near the poles (theta near 0 or pi), wraps by +/- PI;
     * otherwise wraps by +/- 2*PI to find the shortest azimuthal path.
     *
     * @param p1 First point (theta, phi).
     * @param p2 Second point (theta, phi).
     * @return Coordinate difference (delta_theta, delta_phi).
     */
    Vector2 ElasticGeometrySphericalCoor::calculateCoordinateDiff(Vector2 p1, Vector2 p2) { 
        Vector2 res = p2 - p1;
        if (true && (p1[0] < 1e-3 || p2[0] < 1e-3 || PI - p1[0] < 1e-3 || PI - p2[0]<1e-3)) {
            if (abs(res[1]) > abs(res[1] - 1 * PI)) res[1] = res[1] - 1 * PI;
            else if (abs(res[1]) > abs(res[1] + 1 * PI))  res[1] = res[1] + 1* PI;
            return res;
        }
        else if (abs(res[1]) > abs(res[1] - 2 * PI)) res[1] = res[1] - 2 * PI;
        else if (abs(res[1]) > abs(res[1] + 2 * PI)) res[1] = res[1] + 2 * PI;
        return res;

    };

    
} // namespace surface
} // namespace geometrycentral