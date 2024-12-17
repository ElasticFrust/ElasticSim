#include "ElasticGeometry.h" // my class to be implemented
#include <fstream>
#include <limits>


using namespace geometrycentral;
using namespace geometrycentral::surface;


namespace geometrycentral {
namespace surface {
ElasticGeometry::ElasticGeometry(SurfaceMesh& mesh_)
    : VertexPositionGeometry(mesh_) {

          std::cout << "Elastic Geometry!";
          // vertexPositions = VertexData<Vector3>(mesh_, Vector3{0., 0., 0.});

          // The input vertex positions share storage with vertexPositions, incremented the required counter and make
          // sure
          //// they never get cleared
          // requireVertexPositions();
          // vertexPositionsQ.clearable = false;
      };
} // namespace surface
} // namespace geometrycentral