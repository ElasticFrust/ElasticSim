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




using namespace geometrycentral;
using namespace geometrycentral::surface;


 
// ELASTICGEOMTRY is VERTEXPOSITIONGEOMTRY with implemented phenomelogy capable of associating an elastic energy to a memebrane. It includes the following properties:
//
//	Quantities (to be refereshed and declared correctly), and properties, those marked with {} are "refernces" - they do no require updating, and are given as a propertiy :
//		* {REFERENCE METRIC} [edges] - Reference metric (\bar g) is implemented as reference lengths of edges.
//		* {REFERNCE CURVATURE} [edges?] - Prefered angles (\bar g^{-1} \bar{b}) between triangles along an edge. Measured by the reference metric.
//		* ACTUAL CURVATURE [edges? same as reference in any cas] - In principle shoud be actual curvature from the POV of the reference metric, but might prove difficult to calculate. Especially since the angles currespond to actual curvature.
//																Possible. Note that S ~\Delta \theta / \ell -> \bar{a}^{-1} b= \Delta \theta / \bar \ell ??. 
//																Another option - use S -\bar{S}  as this will chagne only quantiatively. (And in anycase this is negligible typically).     $$$ DO THE DETAILED CALC $$$
//		* {PRESSURE} scalar constant in volume - Q: Easy to implement in EOMs, hard to implement in energy, since -> E_p = -\int P dV = -\int P \sqrt{g} d^3 x, but what is d^3 x?  
//												 A: For a convex shape, we can assign an arbitrary "center", to which we know how to calculate the volum of the tetrahedron, created by the triangle "A" and the "center" C.
//													Volume is then V=1/3 A h, where:   A- area of triangle-base, h- the height of C from A.  Even though the center is arbitrary we know V_t=1/3 A_t h_t, and we also know that 
//													V_tot = sum_t V_t = 1/3 sum A_t h_t.  We also know that the total volume is independent of choice of "C" -->? I want to say that it can be written as Atot*H	
//		* {THICKNESS} [faces]
//		* {YOUNGS MODULUS} [faces]
//		* {POISSON'S RATIO}[faces]
//		* Elastic ENERGY [faces] - while the curvatures and lengths might be given on edges,energy is a ssociated with faces: E= \int (g-g\barg)^2 dV ~ \sum_ A (g-g0)^2.  Energy is typically defined in coordinated basis. Translating this to cordinatless is by  noting that for a trangle, it is enoug to consider two of its sides.
//                                  then the metic d_{\mu} f \cdot d_{\nu} f is approximated bu \Delta_{\mu} f \cdot \Delta_{\nu}f.  Furthermor we can work in the natural frame along two of the triangles edges {A,B}, so that e_{A} (intrinsic) is {1,0} and e_{B}= {0,1}, then:
//                                  g_{\mu \nu} --> g_{A B}= \Delta_{A} f \cdot \Delta _{B} f  = {{ l_{A}^2 , \frac{1}{2}(l_{A}^2+l_{B}^2-l_{C}^2)} ,  {\frac{1}{2}(l_{A}^2+l_{B}^2-l_{C}^2) , l_{B}^2}   }, Where l_{n}  is the length (reference or actual) of the n-th edge.
//                                  Note that \frac{1}{2}(l_{A}^2+l_{B}^2-l_{C}^2) = l_{A} l_{B} \cos(\thaeta_{AB}), whre \thaeta_{AB} is the internal angle between A and B.
//		*//									
//	
//
// 
//
//

		
namespace geometrycentral {
namespace surface {

class ElasticGeometry : public VertexPositionGeometry {


  public:

    // REMEBER TO ADD ABILITY TO FIX POSITION/ANGLES WHATERVER. What data structures needed? add those to contrcutor!


    // Construct empty? no geometry, everywhere zero no reference metric
    ElasticGeometry(SurfaceMesh& mesh_); // Basically inherits that of vertex position geometry.

    // Construct from positions (Reference lengths values chosen to be the current values. Set thickness and Y modulus, and poisson ratio to defaults [what defaults?, prob zeros])
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions);

    //Simple Uniform elastics, reference is current)
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions,  const double& THICKNESS,
                    const double& YOUNGs, const double& POISSONs);

    // Simple Uniform elastics, reference is specified)
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions, const EdgeData<double>& L_bar,
                    const EdgeData<double>& B_bar, const double& THICKNESS,
                    const double& YOUNGs, const double& POISSONs);

    // Regions with different elasticit, reference is specified)
    /*ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions, const EdgeData<double>& L_bar,
                    const EdgeData<double>& B_bar, const double& THICKNESS, const double& YOUNGs,
                    const double& POISSONs);*/

    //Constructor with additional fields maybe? Prob. A subclass - Living Geometry? Adapting Geometry? Process Geometry?

    // Detailed constructor (called by all the above basically)
    ElasticGeometry(SurfaceMesh& mesh_, const VertexData<Vector3>& inputVertexPositions, const EdgeData<double>& L_bar,
                    const EdgeData<double>& B_bar, const FaceData<double>& THICKNESS,
                    const FaceData<Eigen::Matrix3f>& ElasticTensor);

    // Destructor
    virtual ~ElasticGeometry() {};


    //Members
    EdgeData<double> Ref_lengths;
    EdgeData<double> Ref_dihedral_angle;
    FaceData<Eigen::Vector3f> Ref_Metric; // Reference Matrix, vectorized.  R[1]= l_1^2 , R[2]=l_2^2, R[3]= (l_1^2+l_2^2-l_3^2)/2;
    FaceData<Eigen::Vector3f> Act_Metric; // Actual Matrix, vectorized (updated based on actual lenths, same scheme.
    FaceData<Eigen::Vector3f> Ref_Curvature;
    FaceData<Eigen::Vector3f> Act_curvature;
    FaceData<Eigen::Matrix3f> Cauchy; //The elastic tensor, reporesneted as a 3X3 matrix.  ***We need to define multiplication rules.***
    FaceData<double> thickness;
    FaceData<double> YoungsModulus;
    FaceData<double> PoissonsRatio;
    FaceData<double> Energy; // \Delta g  A \Delta g  \sqrt{G}
    // MISSINg DATA STRUCTS for constraints! ///
    VertexData<int> Regions; // If we want to specify simple regions with different properties
    VertexData<bool> Fixed_Ver; // If we want to fix specific vertexs
    EdgeData<bool> Fixed_angles;  // If we want to fix specific dihedral  angles.
    //fixing lengths? faces? what else?



  protected:
  private:  
      
};
} // namespace surface
} // namespace geometrycentral