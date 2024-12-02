//Simulation MAIN_FILE

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

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;

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

double my_sum() {
    double SUM = 0.0;
    for (Edge e : mesh->edges()) {
        SUM += geometry->edgeLength(e);
    }
    return SUM / mesh->nEdges();
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
    std::string filepath = "C:/Users/dgrossma/Documents/GitHub/DG-DDG/input/bunny.obj";
    if (inputFilename) {
        filepath = args::get(inputFilename);
    }

    // Load mesh
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(filepath);

    // Initialize polyscope
    polyscope::init();
    polyscope::options::alwaysRedraw = false;

    // Set the callback function
    //polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    

    geometry->requireFaceAreas();
    geometry->requireEdgeLengths();
    geometry->requireVertexNormals();
    geometry->requireVertexMeanCurvatures();
    geometry->requireVertexNormals();

    psMesh->addFaceScalarQuantity("Areas", geometry->faceAreas);
    psMesh->addEdgeScalarQuantity("Lengths", geometry->edgeLengths);
    psMesh->addVertexVectorQuantity("Normals", geometry->vertexNormals);
    psMesh->addVertexScalarQuantity("Mean Curvature", geometry->vertexMeanCurvatures);


    polyscope::view::lookAt(glm::vec3{0, 2, 5}, glm::vec3{0., 0., 0.});
    polyscope::show();


    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->inputVertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));


    geometry->requireFaceAreas();
    geometry->requireEdgeLengths();
    geometry->requireVertexNormals();
    geometry->requireVertexMeanCurvatures();
    geometry->requireVertexNormals();

    psMesh->addFaceScalarQuantity("Areas", geometry->faceAreas);
    psMesh->addEdgeScalarQuantity("Lengths", geometry->edgeLengths);
    psMesh->addVertexVectorQuantity("Normals", geometry->vertexNormals);
    psMesh->addVertexScalarQuantity("Mean Curvature", geometry->vertexMeanCurvatures);
   
    solver(.0001, 1000, mesh, geometry,psMesh);





   psMesh = polyscope::registerSurfaceMesh("changed mesh", geometry->vertexPositions, mesh->getFaceVertexList(),
                                            polyscopePermutations(*mesh));

   geometry->requireFaceAreas();
   geometry->requireEdgeLengths();
   geometry->requireVertexNormals();
   geometry->requireVertexMeanCurvatures();
   geometry->requireVertexNormals();

   psMesh->addFaceScalarQuantity("Areas", geometry->faceAreas);
   psMesh->addEdgeScalarQuantity("Lengths", geometry->edgeLengths);
   psMesh->addVertexVectorQuantity("Normals", geometry->vertexNormals);
   psMesh->addVertexScalarQuantity("Mean Curvature", geometry->vertexMeanCurvatures);
    
   polyscope::show();

    
    

    // Add visualization options.
    //flipZ();
    //psMesh->setSmoothShade(true);
    //psMesh->setSurfaceColor({1.0, 0.45, 0.0});
    //computeNormals();
    //normalVectors = psMesh->addVertexVectorQuantity("Normals", EW_vectors);
    //normalVectors->setEnabled(false);
    //computeShaded();
    //vertexColors = psMesh->addVertexColorQuantity("Plot", shaded_colors);
    //vertexColors->setEnabled(true);

    //// Initialize quantities.
    //TOTAL_ANGLE_DEFECT = 0; //  geometry->totalAngleDefect();
    //EULER_CHARACTERISTIC = 0; //  geometry->eulerCharacteristic();

    // Give control to the polyscope gui
    //polyscope::show();

    return EXIT_SUCCESS;
}