//Simulation MAIN_FILE

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "ElasticGeometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "colormap.h"
#include <map>
#include <numeric>

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<ElasticGeometry> EG;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

//Read vertexPosition, and creates a (weighted) confifguration scalar(s).
std::tuple<std::unique_ptr<Eigen::VectorXd>, std::unique_ptr<Eigen::VectorXd>, std::unique_ptr<Eigen::VectorXd>>  
                  to_configuration(const std::unique_ptr<VertexPositionGeometry>& GM, const std::unique_ptr<ManifoldSurfaceMesh>& SM) {


   
   std::unique_ptr<Eigen::VectorXd> x = std::make_unique<Eigen::VectorXd>(SM->nVertices());
   std::unique_ptr<Eigen::VectorXd> y = std::make_unique<Eigen::VectorXd>(SM->nVertices());
   std::unique_ptr<Eigen::VectorXd> z = std::make_unique<Eigen::VectorXd>(SM->nVertices());
   VertexData<Vector3> VD = GM->vertexPositions;
   Eigen::SparseMatrix<double>& Mass = GM->vertexGalerkinMassMatrix;

    //double temp;
    //double temp2;
    size_t iV=0;
    for (Vertex v : SM->vertices()) {
        //temp = VD[v].x;
        //temp2 = x[iV];
        (*x)[iV] = VD[v].x;
        //temp2 = x[iV];
        (*y)[iV] = VD[v].y;
        (*z)[iV] = VD[v].z;
        iV++;
    }
    
    *x = Mass * (*x);
    *y = Mass * (*y);
    *z = Mass * (*z);
    
    return std::make_tuple(std::move(x),std::move(y), std::move(z));
}

//Set new configurations into vertexPositions
void updategeometry(const std::unique_ptr<VertexPositionGeometry>&geometry,
    const std::unique_ptr<ManifoldSurfaceMesh>&mesh, const std::unique_ptr<Eigen::VectorXd>&newx,
    const std::unique_ptr<Eigen::VectorXd>&newy, const std::unique_ptr<Eigen::VectorXd>&newz)
{

    VertexData<Vector3> VP = geometry->vertexPositions;
    size_t iV = 0;
    for (Vertex v : mesh->vertices()) {
        // temp = VD[v].x;
        // temp2 = x[iV];
        VP[v].x = (*newx)[iV];
        VP[v].y = (*newy)[iV];
        VP[v].z = (*newz)[iV];
        iV++;
    }
    VertexData<Vector3> dif;
    Vector3 difvec;
    for (Vertex v : mesh->vertices()) {
        difvec = VP[v];
        // dif[v] = difvec;
    }
    geometry->vertexPositions = VP;

    geometry->refreshQuantities();
}

//Advance Lapalacian solution by one step
void one_step_solve(const double stepsize, const std::unique_ptr<VertexPositionGeometry>& GM, const std::unique_ptr<ManifoldSurfaceMesh>& SM,int i)
{

    // This  IF is for debug
    if (i<=-3) //Debugging.
    {
        // Creat configurations {x,y,z} before advancement of lapalacian, and {newx,newy,newz} after.
        std::cout << "creating vectors\n";
        std::unique_ptr<Eigen::VectorXd> mass_x;
        std::unique_ptr<Eigen::VectorXd> mass_y;
        std::unique_ptr<Eigen::VectorXd> mass_z;
        std::unique_ptr<Eigen::VectorXd> newx = std::make_unique<Eigen::VectorXd>(mesh->nVertices());
        std::unique_ptr<Eigen::VectorXd> newy = std::make_unique<Eigen::VectorXd>(mesh->nVertices());
        std::unique_ptr<Eigen::VectorXd> newz = std::make_unique<Eigen::VectorXd>(mesh->nVertices());
        VertexData<Vector3> VP = geometry->vertexPositions;


        // This should be reiteratated
        // geometry->unrequireCotanLaplacian();
        std::cout << "getting laplacian\n";
        geometry->requireCotanLaplacian();
        // geometry->unrequireVertexGalerkinMassMatrix();
        std::cout << "getting mass\n";
        geometry->requireVertexGalerkinMassMatrix();
        std::cout << "getting laplacian values\n";
        Eigen::SparseMatrix<double>& Lap =
            geometry->cotanLaplacian; // as a reference they should be automatically updated
        std::cout << "getting mass values\n";
        Eigen::SparseMatrix<double>& Mass = geometry->vertexGalerkinMassMatrix;
        std::cout << "getting ID values\n";
        Eigen::SparseMatrix<double> Id = Eigen::SparseMatrix<double>::SparseMatrix(Lap);
        Id.setIdentity();


        double h = stepsize;

        std::cout << "getting operator values\n";
        Eigen::SparseMatrix<double> opMat = Mass + h * Lap;
        // Eigen::SparseMatrix<double> disc_lap_op = Id - h * Lap;
        std::cout << "getting vertex positions values\n";
        geometry->requireVertexPositions();
        std::cout << "getting configuration\n";
        std::tie(mass_x, mass_y, mass_z) = to_configuration(geometry, mesh);

        // solving;
        std::cout << "creating solver\n";
        //Eigen::SparseLUEigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, 1, Eigen::COLAMDOrdering<int>> solver;
        std::cout << "solver analyze\n";
        solver.analyzePattern(opMat);
        std::cout << "solver factorization\n";
        solver.factorize(opMat);
        std::cout << "solve_x\n";
        *newx = solver.solve(*mass_x);
        std::cout << "solve_y\n";
        *newy = solver.solve(*mass_y);
        std::cout << "solve_z\n";
        *newz = solver.solve(*mass_z);

        std::cout << "updating geometry\n";
        updategeometry(geometry, mesh, newx, newy, newz);

    }
    
    // Creat configurations {x,y,z} before advancement of lapalacian, and {newx,newy,newz} after.
    std::unique_ptr<Eigen::VectorXd> mass_x;
    std::unique_ptr<Eigen::VectorXd> mass_y;
    std::unique_ptr<Eigen::VectorXd> mass_z;
    std::unique_ptr<Eigen::VectorXd> newx = std::make_unique<Eigen::VectorXd>(mesh->nVertices());
    std::unique_ptr<Eigen::VectorXd> newy = std::make_unique<Eigen::VectorXd>(mesh->nVertices());
    std::unique_ptr<Eigen::VectorXd> newz = std::make_unique<Eigen::VectorXd>(mesh->nVertices());
    VertexData<Vector3> VP = geometry->vertexPositions;


    //Solve one step - first by reasing all the data

    
    //geometry->unrequireCotanLaplacian();
    geometry->requireCotanLaplacian();
    //geometry->unrequireVertexGalerkinMassMatrix();
    geometry->requireVertexGalerkinMassMatrix();
    Eigen::SparseMatrix<double>& Lap = geometry->cotanLaplacian; // as a reference they should be automatically updated
    Eigen::SparseMatrix<double>& Mass = geometry->vertexGalerkinMassMatrix;
    Eigen::SparseMatrix<double> Id = Eigen::SparseMatrix<double>::SparseMatrix(Lap);
    Id.setIdentity();
     

    double h = stepsize;
    Eigen::SparseMatrix<double> opMat = Mass + h*Lap;
    // Eigen::SparseMatrix<double> disc_lap_op = Id - h * Lap;

    geometry->requireVertexPositions();
    std::tie(mass_x, mass_y, mass_z) = to_configuration(geometry, mesh);

    // solving;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(opMat);
    solver.factorize(opMat);
    *newx = solver.solve(*mass_x);
    *newy = solver.solve(*mass_y);
    *newz = solver.solve(*mass_z);

    geometry->refreshQuantities();

    updategeometry(geometry, mesh, newx, newy, newz);

}
 
// The whole sovler  - calls "ONE_STEP_SOLVE"  as much as needed
int solver(const double stepsize, const int iter_num, const std::unique_ptr<ManifoldSurfaceMesh>& mesh,
           const std::unique_ptr<VertexPositionGeometry>& geometry, polyscope::SurfaceMesh* psMesh) {
    std::cout << "Calculation beginning\n";
    std::cout << "0/" << iter_num << " steps passed ";

    for (int i = 0; i < iter_num; i++) {
        one_step_solve(stepsize, geometry, mesh,i);
        std::cout <<"\r" << i + 1 << "/"<< iter_num << " steps passed";

       psMesh->updateVertexPositions(geometry->vertexPositions);
        geometry->requireFaceAreas();
        geometry->requireEdgeLengths();
        geometry->requireVertexNormals();
        geometry->requireVertexMeanCurvatures();
        geometry->requireVertexNormals();

        psMesh->addFaceScalarQuantity("Areas", geometry->faceAreas);
        psMesh->addEdgeScalarQuantity("Lengths", geometry->edgeLengths);
        psMesh->addVertexVectorQuantity("Normals", geometry->vertexNormals);
        psMesh->addVertexScalarQuantity("Mean Curvature", geometry->vertexMeanCurvatures);

       // polyscope::frameTick();

         polyscope::view::lookAt(glm::vec3{0, 2, 5}, glm::vec3{0., 0., 0.});

        //Adjust some screenshot default settings if you'd like
        std::string filename = "D:/code output/geometry/screenshots_raw/bunny_screenshot_00" + std::to_string(i) + ".png";
        //std::cout << filename;
        polyscope::screenshot(filename, true);

    }
    std::cout << "\r" << iter_num << "/" << iter_num << " steps passed\n";
    std::cout << "Calculation Complete";
    return 0;
}




void mySubroutine() {
    std::string folder = "D:/code output/geometry/screenshots_raw/";
    double gradAmp = 0;
    double height = 0;
    for (Vertex v : mesh->vertices()) gradAmp += EG->elasticGradient[v].norm2();
    for (Vertex v : mesh->vertices()) height = std::max(height, EG->vertexPositions[v].y);
    double stepsize = .02;
    int max_steps = 10000;
    int count = 0;
    int ss_count = 0;
    double last_ener = EG->elasticEnergy.toVector().sum();
    EG->requireVertexDualAreas();
    EG->requireVertexNormals();
    VertexData<Vector3> forces = VertexData<Vector3>(EG->mesh, Vector3{0, 0, 0});
    while (gradAmp > -1e-3 &&  count <= max_steps) {
        count++;
        for (Vertex v : mesh->vertices()) {
            forces[v] = stepsize * (EG->elasticGradient[v] + EG->vertexNormals[v] * EG->vertexDualAreas[v] * EG->pressure);
        }
        EG->vertexPositions += forces;
        EG->refreshQuantities();
        EG->computeGradient();
        gradAmp = 0;
        height = 0;
        for (Vertex v : mesh->vertices()) gradAmp += EG->elasticGradient[v].norm2();
        for (Vertex v : mesh->vertices()) height = std::max(height, EG->vertexPositions[v].y);
        if (count % 10 == 0) {
            polyscope::refresh();
            std::cout << "Count: " << count << "\t" << "Total Gradient: " << gradAmp << "\t"
                      << "Average Gradient: " << gradAmp / mesh->nVertices() << "\t"
                      << "dEnergy: " << EG->elasticEnergy.toVector().sum() - last_ener << "\n";
            std::string file= folder + "bunny_grad_" + std::to_string(ss_count) + ".png";
            psMesh->updateVertexPositions(EG->vertexPositions);
            psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy);
            //psMesh->addVertexVectorQuantity("Gradient", EG->elasticGradient);
            if (count % 10 == 0) {
                polyscope::screenshot(file, true);
                ss_count += 1;
            }
            last_ener = EG->elasticEnergy.toVector().sum();
        }
    }
}
   
       


void myCallback() {

    // Since options::openImGuiWindowForUserCallback == true by default,
    // we can immediately start using ImGui commands to build a UI

    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
                               // instead of full width. Must have
                               // matching PopItemWidth() below.

    if (ImGui::Button("run subroutine")) {
        // executes when button is pressed
        mySubroutine();
    }
    /*ImGui::SameLine();
    if (ImGui::Button("hi")) {
        polyscope::warning("hi");
    }*/

    ImGui::PopItemWidth();
}


int main(int argc, char** argv) {

        // Configure the argument parser
    args::ArgumentParser parser("15-458 HW2");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help&) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // If a mesh name was not given, use default mesh.
    std::string filepath = "C:/Users/dgrossma/Documents/GitHub/DG-DDG/input/sphere.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }

    // Load mesh
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filepath);
    //std::unique_ptr<ElasticGeometry> EG1(new ElasticGeometry(*mesh));
   // std::unique_ptr<VertexPositionGeometry> EG2(new VertexPositionGeometry(*mesh));
    //std::unique_ptr<ElasticGeometry> EG3(new ElasticGeometry(
        //*mesh, geometry->inputVertexPositions, EdgeData<double>(*mesh, 0), EdgeData<double>(*mesh, 0),
        //FaceData<double>(*mesh, 0), FaceData<Eigen::Matrix3f>(*mesh, Eigen::Matrix3f()), 0));
    std::unique_ptr<ElasticGeometry> BG(new ElasticGeometry(*mesh, geometry->inputVertexPositions, 1, 1, .5, 1));
    EG = std::move(BG);
    EG->requireElasticEnergy();

    
    
    std::cout << "\n  Reference Metric [0]: \n";
    std::cout << EG->referenceMetric[0][0] << ", \t";
    std::cout << EG->referenceMetric[0][1] << ", \t";
    std::cout << EG->referenceMetric[0][2] << "\n";

    std::cout << "\n  Actual Metric [0]: \n";
    std::cout << EG->actualMetric[0][0] << ", \t";
    std::cout << EG->actualMetric[0][1] << ", \t";
    std::cout << EG->actualMetric[0][2] << "\n";
    
    
    std::cout << "\n  Cauchy Tensor [0]:\n";
    std::cout << EG->elasticCauchyTensor[0](0, 0) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](0, 1) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](0, 2) ;
    std::cout << "\n";
    std::cout << EG->elasticCauchyTensor[0](1, 0) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](1, 1) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](1, 2) ;
    std::cout << "\n";
    std::cout << EG->elasticCauchyTensor[0](2, 0) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](2, 1) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](2, 2) << "\n";

    std::cout << "\n  Energy [0]:\n";
    std::cout << EG->elasticEnergy[0] << "\n";



    std::cout << "\n  Inflating\n";
   
    VertexData<Vector3> VP = EG->vertexPositions;
    for (Vertex v : mesh->vertices()) {
        VP[v].y *= 1.1;
        VP[v].x *= 1.1;
        VP[v].z *= 1.1;
    }
    EG->vertexPositions = VP;
    EG->refreshQuantities();

     std::cout << "\n  Reference Metric [0]: \n";
    std::cout << EG->referenceMetric[0][0] << ", \t";
    std::cout << EG->referenceMetric[0][1] << ", \t";
    std::cout << EG->referenceMetric[0][2] << "\n";

    std::cout << "\n  Actual Metric [0]: \n";
    std::cout << EG->actualMetric[0][0] << ", \t";
    std::cout << EG->actualMetric[0][1] << ", \t";
    std::cout << EG->actualMetric[0][2] << "\n";


    std::cout << "\n  Cauchy Tensor [0]:\n";
    std::cout << EG->elasticCauchyTensor[0](0, 0) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](0, 1) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](0, 2);
    std::cout << "\n";
    std::cout << EG->elasticCauchyTensor[0](1, 0) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](1, 1) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](1, 2);
    std::cout << "\n";
    std::cout << EG->elasticCauchyTensor[0](2, 0) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](2, 1) << ", \t";
    std::cout << EG->elasticCauchyTensor[0](2, 2) << "\n";

    std::cout << "\n  Energy [0]:\n";
    std::cout << EG->elasticEnergy[0] << "\n";


    EG->computeGradient();

   

    // Initialize polyscope
    polyscope::init();
    polyscope::options::alwaysRedraw = true;

    // Set the callback function
    //polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), EG->vertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy);
    psMesh->addVertexVectorQuantity("Gradient", EG->elasticGradient);
    


    //polyscope::state::userCallback = myCallback;  
    
    polyscope::show();

   /* VP = EG->vertexPositions;
    for (Vertex v : mesh->vertices()) {
        VP[v].y *= 1/1.01;
        VP[v].x *= 1 / 1.01;
        VP[v].z *= 1 / 1.01;
    }
    EG->vertexPositions = VP;
    EG->refreshQuantities();
    psMesh->updateVertexPositions(EG->vertexPositions);
    psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy);

    polyscope::show();*/


    mySubroutine();


    return EXIT_SUCCESS;
}