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

#include "ElasticGeometry.h" // my class to be implemented
//#include "ElasticGeometry.cpp" // my class to be implemented
#include <fstream>
#include <limits>

using namespace geometrycentral;
using namespace geometrycentral::surface;


namespace geometrycentral {
namespace surface {

class ElasticGeometrySphericalCoor : public ElasticGeometry {


  public:

    

    ElasticGeometrySphericalCoor(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_);

    ElasticGeometrySphericalCoor(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions_,
                                 const double& THICKNESS_, const double& YOUNGs_, const double& POISSONs_,
                                 const double& PRESSURE_);

    virtual ~ElasticGeometrySphericalCoor() {};


    void computeReferenceMetric() override;
    void computeActualMetric() override;
    void computeReferenceCurvature() override;
    void computeActualCurvature() override;
   // void computeElasticEnergy() override; //no need we are changing the references, but caluclation should be the same
    void computeGradient() override; //we do need  to override gradient - the base class calculation is horrible 
    
   
    VertexData<Vector2> vertexCoordinates; //    = VertexData<Vector2>(this->mesh, Vector2({0, 0}));
    void requireVertexCoordinates();
    void unrequireVertexCoordinates();
    FaceData<Vector2> faceCentroidCoordinates; //  = FaceData<Vector2>(this->mesh, Vector2({0, 0}));
    void requireFaceCentroidCoordinates();
    void unrequireFaceCentroidCoordinates();
    FaceData<Vector3> faceCentroidPosition; //  = FaceData<Vector3>(this->mesh, Vector3({0, 0, 0}));
    void requireFaceCentroidPosition();
    void unrequireFaceCentroidPosition();

    void updateElasticCauchyTensor();
    void updateFaceCentroidCoordinates(const FaceData<Vector2> faceCoordinates);


  protected:
    // == Quantities
    DependentQuantityD<FaceData<Vector3>> faceCentroidPositionQ;
    virtual void computeCentroids();
    DependentQuantityD<VertexData<Vector2>> vertexCoordinatesQ;
    virtual void computeVertexCoordinates();
    DependentQuantityD<FaceData<Vector2>> faceCentroidCoordinatesQ;
    virtual void computeFaceCentroidCoordinates();



  private:
    void calculateFaceActualMetric(Face f);
    void calculateFaceReferenceMetric(Face f);
    void calcualteFaceActualCurvature(Face f);
    void calcualteFaceReferenceCurvature(Face f);
        
    void localEnergyChange(const Vertex v);

    void calculateAdjacentNormalAndCentroids(const Vertex v);
    void calculateAdjacentMetric(const Vertex v);
    void calculateAdjacentCurvature(const Vertex v);
    

    Eigen::Vector3f calculateCurvature(Face f);
    Eigen::Vector3f calculateMetric(Face f, float lengths[3]);
    Vector2 calculateCoordinateDiff(Vector2 p1, Vector2 p2);
    
    bool isReferenceMetricInitializedF = false;
    bool isReferenceCurvatureInitializedF = false;




};
} // namespace surface
} // namespace geometrycentral