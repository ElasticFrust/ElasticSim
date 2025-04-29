//Simulation MAIN_FILE

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "ElasticGeometry.h"
#include "geometrycentral/surface/rich_surface_mesh_data.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "colormap.h"
#include <map>
#include <numeric>

#include <iostream>
#include <fstream> 

#include <Windows.h>


using namespace geometrycentral;
using namespace geometrycentral::surface;




// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<ElasticGeometry> EG;
std::unique_ptr<RichSurfaceMeshData> richData;

std::string workingFolder;
int snapshotEvery;
int printinEvery;

double Tinit, Ttarget;



// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

//Read vertexPosition, and creates a (weighted) confifguration scalar(s).    // GOOD code, no use,
/* std::tuple<std::unique_ptr<Eigen::VectorXd>, std::unique_ptr<Eigen::VectorXd>, std::unique_ptr<Eigen::VectorXd>>  
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
*/

//Not important can be deleted
/*
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
*/    // Not important can delete       

double logfunc(double x) {
    return std::log10(x+1e-10);
}


template <typename datatype, typename meshtype>
double mean_func(MeshData<meshtype,datatype>& arg) {
    return arg.toVector().mean();
}


double mean_abs(VertexData<Vector3>& vec) {
    VertexData<double> absvec = VertexData<double>(*vec.getMesh(), 0);
    for (Vertex v : vec.getMesh()->vertices()) absvec[v] = vec[v].norm();
    return mean_func(absvec);
}


int saveLog(std::vector<std::string> heads, std::vector<std::string> data, std::string folder, std::string file,
            bool append = false) {
    std::ofstream logfile;
    if ((data.size() != heads.size()) && !data.empty()) {
        std::cout << "Header and data not the same size, Log may be corrupted";
        return 1;
    }
    if (append) {
        logfile.open(folder + file, std::ios_base::app);
        for (std::string datum : data) {
            logfile << datum << "\t";
        }
        logfile << "\n";
    }
    else { // First time = write only header
        logfile.open(folder + file);
        for (std::string header : heads) {
            logfile << header << "\t";
        }
        logfile << "\n";
    }
    return 0;
}


void printline(std::vector<std::string>& headers, std::vector<std::string>& data) {
    for (int index = 0; index < headers.size(); index++) {
        std::cout << headers[index] << ": " << data[index] << "\t";
    }
    std::cout << "\n";
}


int writeRichData(RichSurfaceMeshData& RD, ElasticGeometry& geo, std::string file) {
    // commented out qquanitties don't have the correct conversion.  Since we were in any case about the chagne some of
    // therir formats and since they are currenyl *not* basic quantities, we live withoutit for the meanwhile
    // 
    //RD.addMeshConnectivity(); // Calculate once, if that changes - recreat richdata
    RD.addGeometry(geo);

    RD.addFaceProperty("Thickness", geo.thickness);
    RD.addFaceProperty("Pressure", FaceData<double>(geo.mesh, geo.pressure));
    RD.addFaceProperty("Youngs_Modulus", geo.youngsModulus);
    RD.addFaceProperty("Poissons_Ratio", geo.poissonsRatio);
    // RD.addFaceProperty("Elastic Tensor", geo.elasticCauchyTensor);

    RD.addEdgeProperty("Reference_Lengths", geo.referenceLengths);
    RD.addEdgeProperty("Reference_Dihedral_Angles", geo.referenceEdgeDihedralAngles);
    // RD.addFaceProperty("Reference Metric", geo.referenceMetric);
    // RD.addFaceProperty("Reference Curvature", geo.referenceCurvature);

    RD.addEdgeProperty("Actual_Lengths", geo.edgeLengths);
    RD.addEdgeProperty("Actual_Dihedral_Angles", geo.edgeDihedralAngles);
    // RD.addFaceProperty("Actual Metric", geo.actualMetric);
    // RD.addFaceProperty("Actual Curvature", geo.actualCurvature);

    //RD.addVertexProperty("Elastic Gradient", geo.elasticGradient); //seemingly should work but breaks code 

    RD.addFaceProperty("Bending_Energy_density", geo.bendingEnergy);
    RD.addFaceProperty("Stretching_Energy_density", geo.stretchingEnergy);
    RD.addFaceProperty("Elastic_Energy", geo.elasticEnergy);


    // RD.addFaceProperty("Reference Area",geo.referenceAreas);


    RD.write(file);

    return 0;
}

int readRichaData(RichSurfaceMeshData& RD, ElasticGeometry& geo) {
    geo.thickness = RD.getFaceProperty<double>("Thickness");
    geo.pressure = RD.getFaceProperty<double>("Pressure")[0];
    geo.youngsModulus = RD.getFaceProperty<double>("Youngs_Modulus");
    geo.poissonsRatio = RD.getFaceProperty<double>("Poissons_Ratio");
    // geo.elasticCauchyTensor = RD.getFaceProperty<Eigen::Matrix3f>("Elastic Tensor");

    geo.referenceLengths = RD.getEdgeProperty<double>("Reference_Lengths");
    geo.referenceEdgeDihedralAngles = RD.getEdgeProperty<double>("Reference_Dihedral_Angles");
    // geo.referenceMetric = RD.getFaceProperty<Eigen::Vector3f>("Reference Metric");
    // geo.referenceCurvature = RD.getFaceProperty<Eigen::Vector3f>("Reference Curvature");

    geo.edgeLengths = RD.getEdgeProperty<double>("Actual_Lengths");
    geo.edgeDihedralAngles = RD.getEdgeProperty<double>("Actual_Dihedral_Angles");
    // geo.actualMetric = RD.getFaceProperty<Eigen::Vector3f>("Actual Metric");
    // geo.actualCurvature = RD.getFaceProperty<Eigen::Vector3f>("Actual Curvature");

    //geo.elasticGradient = RD.getVertexProperty<Vector3>("Elastic Gradient");

    geo.bendingEnergy = RD.getFaceProperty<double>("Bending_Energy_density");
    geo.stretchingEnergy = RD.getFaceProperty<double>("Stretching_Energy_density");
    geo.elasticEnergy = RD.getFaceProperty<double>("Elastic_Energy");


    // geo.referenceAreas = D.getFaceProperty<double>("Reference Area");

    geo.refreshQuantities();
    geo.computeGradient();

    return 0;
}


bool stabilitycheck(VertexData<Vector3> current, VertexData<Vector3> old) {
    bool stable = true;
    for (Vertex v : mesh->vertices()) {
        Vector3 diff = current[v] - old[v];
        double currentnorm = current[v].norm();
        double diffnorm = diff.norm();
        double oldnorm = old[v].norm();
        if (diffnorm / (oldnorm + 1e-6) > 1) 
            stable = false;
    }
    return stable;
}

void mySubroutine() {
    std::string screen_folder = workingFolder;
    std::string Log_folder = workingFolder;
    std::string log_file = "log.log";
    std::string datafile_w;
       
    double gradAmp = 0;
    double height = 0;
    for (Vertex v : mesh->vertices()) gradAmp += EG->elasticGradient[v].norm2();
    for (Vertex v : mesh->vertices()) height = std::max(height, EG->vertexPositions[v].y);
    double stepsize;
    bool stepsize_updated_flag =false;
    int max_steps = 1000000;
    int count = 0;
    int reg_count = 0;
    int ss_count = 0;
    double last_ener = EG->elasticEnergy.toVector().sum();
    EG->requireVertexDualAreas();
    EG->requireVertexNormals();
    EG->requireVertexMeanCurvatures();
    EG->requireFaceAreas();
    double d_ener = 1;
    double rel_ener = 1;
    double tot_ener = EG->elasticEnergy.toVector().sum();
    std::time_t* stopwatch =  0;
    VertexData<Vector3> forces_last = VertexData<Vector3>(EG->mesh, Vector3{10,10,10});
    VertexData<Vector3> forces = VertexData<Vector3>(EG->mesh, Vector3{10, 10, 10}); 
    std::vector<std::string> headers = {"Iteration",
                                        "Epoch",
                                        "Gradient",
                                        "Total Energy",
                                        "dEnergy",
                                        "rel dEnergy",
                                        "Average Mean Curvature",
                                        "pressure",
                                        "step size",
                                        "Number of triangles"};
    std::vector<std::string> data;    
    saveLog(headers, data, Log_folder, log_file);
    int press_reg_steps = 5000;
    /*for (int reg_step = 1; reg_step <= press_reg_steps; reg_step++) {
        std::cout << "\n \n  Current reg_step: " << reg_step << "/" << press_reg_steps << "\n \n";*/
        stepsize = .01;
        reg_count = 0;
        double pres_func = 0;
        int thickness_reg_steps =20000;
        double thick_func;   

        do {
            count++;
            reg_count++;
            if (stepsize < 0.01 && count % 3000 == 1 && count > 2900) stepsize *=2;
            pres_func = std::min(count * 1.0, press_reg_steps * 1.0) / press_reg_steps/1.0 * EG->pressure;
            thick_func =
                std ::min(count * 1.0, thickness_reg_steps * 1.0) * (Ttarget - Tinit) / thickness_reg_steps + Tinit;
           // std::cout << pres_func << "\n";
            stepsize_updated_flag = false;
            // std::cout << count << "\n";
            
            for (Vertex v : mesh->vertices()) {
                forces[v] =
                    stepsize * (EG->elasticGradient[v] + pres_func* EG->vertexNormals[v] * EG->vertexDualAreas[v]  );
            }
            if (reg_count > 1) { /// CHAGNE
                if (stabilitycheck(forces, forces_last)) {
                    forces_last = forces;
                } else {
                    for (Vertex v : mesh->vertices()) {
                        forces[v] = -forces_last[v];
                    }
                    stepsize *= 0.5;
                    std::cout << "\n \n  Unstable, adapting step size... Recalculating last step with new stepsize: "
                              << stepsize << "\n \n";
                    stepsize_updated_flag = true;
                }
            } else {
                forces_last = forces;
            }
            EG->vertexPositions += forces;
            for (Face f : mesh->faces()) EG->thickness[f] = thick_func;
            EG->refreshQuantities();
            EG->computeGradient();
            gradAmp = 0;
            height = 0;
            for (Vertex v : mesh->vertices()) gradAmp += EG->elasticGradient[v].norm2();
            for (Vertex v : mesh->vertices()) height = std::max(height, EG->vertexPositions[v].y);
            last_ener = tot_ener;
            tot_ener = EG->elasticEnergy.toVector().sum();
            if (!stepsize_updated_flag) {
                d_ener = (tot_ener - last_ener)/stepsize;
                rel_ener = d_ener / last_ener;
            }
            
            if ((count - 1) % printinEvery == 0) {
                data.assign({std::to_string(count), 
                             std::to_string(std::time(stopwatch)),
                             std::to_string(mean_abs(forces)),
                             std::to_string(tot_ener),
                             std::to_string(d_ener), 
                             std::to_string(rel_ener),
                             std::to_string(mean_func(EG->vertexMeanCurvatures / EG->vertexDualAreas)),
                             std::to_string(pres_func),
                             std::to_string(stepsize),
                             std::to_string(mesh->nFaces())});
                printline(headers, data);
                saveLog(headers, data, Log_folder, log_file, true);
            }

            if ((count - 1) % snapshotEvery == 0) {
                polyscope::view::resetCameraToHomeView();
                polyscope::refresh();
                std::string file1 =
                    screen_folder + "energy_" + std::to_string(count) + "order_" + std::to_string(ss_count) + ".png";
                std::string file2 =
                    screen_folder + "curvature_" + std::to_string(count) + "order_" + std::to_string(ss_count) + ".png";

                datafile_w =
                    screen_folder + "RichData_" + std::to_string(count) + "order_" + std::to_string(ss_count) + ".ply";

                psMesh->updateVertexPositions(EG->vertexPositions);
                psMesh->setAllQuantitiesEnabled(false);
                auto energy_int = psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy);
                /*auto energy_int_log =
                    psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy.toVector().unaryExpr(&logfunc));*/
                auto energy_quantity = psMesh->addFaceScalarQuantity(
                    "Elastic Energy Content", (EG->elasticEnergy / EG->faceAreas)); //.toVector().unaryExpr(&logfunc));
                energy_quantity->setMapRange(
                    std::make_pair(0, 0.005)); //(-5, -3)); // EG->elasticEnergy.toVector().maxCoeff()));
                energy_quantity->setColorMap("coolwarm");
                energy_quantity->setEnabled(true);
                energy_quantity->draw();
                // auto energy_quantity_log = psMesh->addFaceScalarQuantity(
                //"Elastic Energy Content (log)", (EG->elasticEnergy / EG->faceAreas).toVector().unaryExpr(&logfunc));
                // energy_quantity_log->setMapRange(std::make_pair(-5, -1)); //
                // EG->elasticEnergy.toVector().maxCoeff())); energy_quantity_log->setColorMap("coolwarm");
                // energy_quantity_log->setEnabled(true);
                // energy_quantity_log->draw();

                polyscope::screenshot(file1, true);

                auto stretch_ener = psMesh->addFaceScalarQuantity("Stretching Energy density", EG->stretchingEnergy);
                auto BEND_ener = psMesh->addFaceScalarQuantity("bending Energy density", EG->bendingEnergy);
                auto gradient_quantity = psMesh->addVertexVectorQuantity("Gradient", EG->elasticGradient);
                gradient_quantity->draw();
                auto normals_quantity = psMesh->addVertexVectorQuantity("Normal (pressur)", EG->vertexNormals);
                normals_quantity->draw();
                auto angles_quant_vis = psMesh->addEdgeScalarQuantity(
                    "dihedral angles difference", EG->edgeDihedralAngles - EG->referenceEdgeDihedralAngles);
                auto ref_angles_vis =
                    psMesh->addEdgeScalarQuantity("reference dihedral angles", EG->referenceEdgeDihedralAngles);
                auto angles_vis = psMesh->addEdgeScalarQuantity("actual dihedral angles", EG->edgeDihedralAngles);
                // auto BEND_ener2 = psMesh->addFaceScalarQuantity("bending Energy", EG->bendingEnergy);

                // if(d_ener>0)  polyscope::show();
                auto curv_quantity =
                    psMesh->addVertexScalarQuantity("MeanCurvature", EG->vertexMeanCurvatures / EG->vertexDualAreas);
                curv_quantity->setMapRange(std::make_pair(0, .5)); // EG->elasticEnergy.toVector().maxCoeff()));
                curv_quantity->setColorMap("viridis");
                curv_quantity->setEnabled(true);
                curv_quantity->draw();
                // psMesh->addVertexVectorQuantity("Gradient", EG->elasticGradient);


               


                /* for (Face f : mesh->faces()) {
                    if (EG->stretchingEnergy[f] < 0) {
                        std::cout << "\nError in iteration: " << count << ".\n";
                        std::cout << "ERROR! Negative stretching energy!  at face: " << f.getIndex()
                                  << ", Energy: " << EG->stretchingEnergy[f] << ".\n ";
                        std::cout << "\n  Reference lengths:  {";
                        int edgecount = 0;
                        for (Edge e : f.adjacentEdges()) {
                            edgecount++;
                            std::cout << EG->referenceLengths[e];
                            if (edgecount == 3)
                                std::cout << "}\n";
                            else
                                std::cout << ",";
                        }

                        std::cout << "\n  Actual lengths:  {";
                        edgecount = 0;
                        for (Edge e : f.adjacentEdges()) {
                            edgecount++;
                            std::cout << EG->edgeLengths[e];
                            if (edgecount == 3)
                                std::cout << "}\n";
                            else
                                std::cout << ",";
                        }

                        std::cout << "\n  Reference Metric: \n";
                        std::cout << EG->referenceMetric[f][0] << ", \t";
                        std::cout << EG->referenceMetric[f][1] << ", \t";
                        std::cout << EG->referenceMetric[f][2] << "\n";

                        std::cout << "\n  Actual Metric: \n";
                        std::cout << EG->actualMetric[f][0] << ", \t";
                        std::cout << EG->actualMetric[f][1] << ", \t";
                        std::cout << EG->actualMetric[f][2] << "\n";


                        std::cout << "\n  Cauchy Tensor:\n";
                        std::cout << EG->elasticCauchyTensor[f](0, 0) << ", \t";
                        std::cout << EG->elasticCauchyTensor[f](0, 1) << ", \t";
                        std::cout << EG->elasticCauchyTensor[f](0, 2);
                        std::cout << "\n";
                        std::cout << EG->elasticCauchyTensor[f](1, 0) << ", \t";
                        std::cout << EG->elasticCauchyTensor[f](1, 1) << ", \t";
                        std::cout << EG->elasticCauchyTensor[f](1, 2);
                        std::cout << "\n";
                        std::cout << EG->elasticCauchyTensor[f](2, 0) << ", \t";
                        std::cout << EG->elasticCauchyTensor[f](2, 1) << ", \t";
                        std::cout << EG->elasticCauchyTensor[f](2, 2) << "\n";
                        // polyscope::show();
                    }
                }*/
                // if (ss_count == 7000) polyscope::show();


                psMesh->addVertexScalarQuantity("MeanCurvature", EG->vertexMeanCurvatures / EG->vertexDualAreas);

                writeRichData(*richData, *EG, datafile_w);
                polyscope::screenshot(file2, true);
                ss_count += 1;
                // polyscope::show();
            }
        } while ((std::abs(rel_ener) > 1e-4 && count <= max_steps && stepsize > 1e-9) ||
                 (count <= press_reg_steps || count <= thickness_reg_steps));
    //}
    std::cout << "\n  \n \t \t SIM COMPLETE! \n \n";    
    if (count > max_steps) std::cout << "max interation exceeded \n";
    if (stepsize <= 1e-6) std::cout << "unstable \n";
    if (std::abs(d_ener) <= 1e-6) std::cout << "minimum found! \n";
    std::cout << "Total iterations: " << count << "\n ";
    std::cout << "Final log: \n";
    printline(headers, data);

    datafile_w = screen_folder + "RichData_Final.ply";   
    writeRichData(*richData, *EG, datafile_w);
    polyscope::show();
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



//int main(int argc, const char** argv) {
//    args::ArgumentParser p("parser");
//    args::CompletionFlag c(p, {"complete"});
//    args::ValueFlag<std::string> f(p, "name", "description", {'f', "foo"}, "abc");
//    args::ValueFlag<std::string> b(p, "name", "description", {'b', "bar"}, "abc");
//    args::MapFlag<std::string, int> m(p, "mappos", "mappos", {'m', "map"}, {{"1", 1}, {"2", 2}});
//    args::Positional<std::string> pos(p, "name", "desc");
//
//    try {
//        p.ParseCLI(argc, argv);
//    } catch (args::Completion& e) {
//        std::cout << e.what();
//    }
//    std::cout << p;
//    std::cout << pos.Get() << "\n";
//    if (f.Get()== "")
//        std::cout << "blabla";
//    else
//        std::cout << f.Get();
//
//    return 0;
//}

int main(int argc, char** argv) {

        // Configure the argument parser
    args::ArgumentParser parser("15-458 HW2");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
    args::ValueFlag<double> thickness(parser, "thickness", "thickness value (uniform)", {'t', "thickness"}, 1);
    args::ValueFlag<double> pressure(parser, "pressure", "pressure value", {'p', 'P', "pressure"}, 10);
    args::ValueFlag<double> Youngs(parser, "Young's Modulus", "Young's Modulus", {'Y', 'E', "youngs"}, 1);
    args::ValueFlag<double> Poissons(parser, "Poisson's ratio", "Poisson's ratio", {'u', "poisson", "poissons"}, .5);  
    args::ValueFlag<std::string> SaveFolderName(parser, "folder name", "Folder to save snapshots and final result", {'f', "folder"}, "D:/code_output/geometry/screenshots_raw/");
    args::ValueFlag<int> printingCoutner(parser, "print counter","print meta data after every [print counter] iterations", {'l', "log", "print"},100);
    args::ValueFlag<int> snapshotCoutner(parser, "snapshot counter", "snapshot after every [snapshot counter] iterations", {'s', "snap", "snapshot"}, 100);
    args::ValueFlag<double> otherVal(parser, "other", "other", {'o', "other"}, 1);
    /*args::MapFlag<std::string, bool> readfromfileflag(parser, "input", "an input file .ply",
                                                      {'r', "read", "in", "input", "readfile", "inputfile"},
                                 {{"true", true}, {"false", false}});*/

    
   
   


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
    std::string filepath = "C:/Users/dgrossma/Documents/GitHub/ElasticSim/input/short_cilinder.obj";//sphere.obj"; //
    if (inputFilename) {
        filepath = args::get(inputFilename);       
    }


     std::cout << "\n Acquired paramerers: \n"
              << "mesh file: " << inputFilename.Get() << "\n"
              << "working folder: " << SaveFolderName.Get() << "\n"
              << "thickness: " << thickness.Get() << "\n"
              << "pressure: " << pressure.Get() << "\n"
              << "Youngs Modulus: " << Youngs.Get() << "\n"
              << "Poissons ratio: " << Poissons.Get() << "\n"
              << "Other: " << otherVal.Get() << "\n"
              << "log interval: " << printingCoutner.Get() << "\n"
              << "snapshot interval: " << snapshotCoutner.Get() << "\n";


    if (CreateDirectoryA(SaveFolderName.Get().c_str(), NULL)) {
        std::cout << "Folder created successfully." << std::endl;
    } else {
        std::cout << "Failed to create folder. (prob. exists)" << std::endl;
    }
 


    workingFolder = SaveFolderName.Get();
    snapshotEvery = snapshotCoutner.Get();
    printinEvery = printingCoutner.Get();


    std::string fileExtension = filepath.substr(filepath.length() - 4);

    


        // Load mesh
    //std::unique_ptr<SurfaceMesh> mesh;      /// ALRAFDY DECLARED
    //std::unique_ptr<VertexPositionGeometry> geometry;
    //std::unique_ptr<ElasticGeometry> EG;
    //std::unique_ptr<RichSurfaceMeshData> richData;

    //std::string datafile = "D:/code output/geometry/screenshots_raw/RichData_5.ply";
    //bool readfromfile = false;


   if (fileExtension == ".ply") {  //rich!
        std::unique_ptr<SurfaceMesh> sMesh;
        std::tie(sMesh, richData) = RichSurfaceMeshData::readMeshAndData(filepath); 
        mesh = std::move(sMesh->toManifoldMesh());
        geometry = richData->getGeometry();
        EG = std::move(std::unique_ptr<ElasticGeometry> (new ElasticGeometry(*mesh, geometry->vertexPositions)));
        EG->vertexPositions = geometry->vertexPositions;
        readRichaData(*richData, *EG);
        richData->outputFormat = happly::DataFormat::ASCII;
    
   }

   else { // fileExtension == ".obj" --> we have a regular start!

        std::tie(mesh, geometry) = readManifoldSurfaceMesh(filepath);
        // std::unique_ptr<ElasticGeometry> EG1(new ElasticGeometry(*mesh));
        // std::unique_ptr<VertexPositionGeometry> EG2(new VertexPositionGeometry(*mesh));
        // std::unique_ptr<ElasticGeometry> EG3(new ElasticGeometry(
        //*mesh, geometry->inputVertexPositions, EdgeData<double>(*mesh, 0), EdgeData<double>(*mesh, 0),
        // FaceData<double>(*mesh, 0), FaceData<Eigen::Matrix3f>(*mesh, Eigen::Matrix3f()), 0));

        geometry->requireEdgeCotanWeights();
        geometry->requireEdgeLengths();

        richData = std::move(std::unique_ptr<RichSurfaceMeshData>(new RichSurfaceMeshData(*mesh)));
        richData->addMeshConnectivity();
        richData->addGeometry(*geometry);
        richData->outputFormat = happly::DataFormat::ASCII;

        polyscope::init();
        polyscope::options::alwaysRedraw = true;

        // Set the callback function
        // polyscope::state::userCallback = functionCallback;

        // Add mesh to GUI
        /*psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->vertexPositions,
                                                mesh->getFaceVertexList(), polyscopePermutations(*mesh));
        psMesh->addEdgeScalarQuantity("Cotan", geometry->edgeCotanWeights);
        polyscope::show();
        */


        VertexData<Vector3> VP = geometry->vertexPositions;
        geometry->requireEdgeLengths();
        geometry->requireVertexTangentBasis();
        geometry->requireEdgeDihedralAngles();

        polyscope::init();
        polyscope::options::alwaysRedraw = true;







        double edgeLmin = geometry->edgeLengths.toVector().minCoeff(); //wiggle wigglwe
        double maxz = 0, maxy = 0;
        for (Vertex v : mesh->vertices()) {
            maxy = std::max(maxy, VP[v].y);
            maxz = std::max(maxz, VP[v].z);
        }
        for (Vertex v : mesh->vertices()) {
             VP[v] += 0.5 * edgeLmin * (randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][0].normalize() +
             randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][1].normalize());            
        }


        for (Vertex v : mesh->vertices()) { //rescale
            VP[v].y *= 1 / maxy * otherVal.Get() / 25;
            VP[v].x *= 1 / maxz; //  *25;
            VP[v].z *= 1 / maxz; // * 25;
        }

        // VP[10].x *= .99;
        // VP[10].y *= .99;
        // VP[10].z *= .99;
        geometry->vertexPositions = VP;
        geometry->refreshQuantities();

        // std::cout << atan2(0, -1) << "\n";
        // std::cout << atan2(0, 1) << "\n";
        // std::cout << atan2(1,0) << "\n";
        // std::cout << atan2(-1,0) << "\n";
        // std::cout << atan2(1, -1) << "\n";
        // std::cout << atan2(1, 1) << "\n";

        // psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath),
        // geometry->vertexPositions,
        //                                         mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        // psMesh->addEdgeScalarQuantity("PolyAngle", geometry->edgeDihedralAngles);
        // polyscope::show();

        
        Ttarget = thickness.Get()/25;
        
        Tinit = 5 * Ttarget;
        if (0.2 * maxz > Tinit) Tinit = 0.2 * maxz;
        std::unique_ptr<ElasticGeometry> BG(new ElasticGeometry(*mesh, geometry->vertexPositions,Tinit,Youngs.Get(),Poissons.Get(),pressure.Get() ));
        EG = std::move(BG);
        EG->requireElasticEnergy();
        EG->requireBendingEnergy();

      std:: cout << pressure.Get();


        /*std::cout << "\n  Reference Metric [0]: \n";
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


        std::cout << "\n  Inflating\n";*/

        //VertexData<Vector3> VP2 = EG->vertexPositions;
        //EdgeData<bool> tochange = EdgeData<bool>(EG->mesh, true);
        //double maxz = 0, maxy = 0, maxx = 0;
        //double minz = 0, miny = 0, minx = 0;
        //for (Vertex v : mesh->vertices()) {
        //    if (maxz < VP2[v].z) maxz = VP2[v].z;
        //    if (maxy < VP2[v].y) maxy = VP2[v].y;
        //    if (maxx < VP2[v].x) maxx = VP2[v].x;

        //    if (minz > VP2[v].z) minz = VP2[v].z;
        //    if (miny > VP2[v].y) miny = VP2[v].y;
        //    if (minx > VP2[v].x) minx = VP2[v].x;
        //}
        //Vector3 center1{(maxx + minx) / 2, maxy, (maxz + minz) / 2};
        //Vector3 center2{(maxx + minx) / 2, miny, (maxz + minz) / 2};
        //double sphereRad = (maxx - minx) / 2;
        //double rad = 1 * sphereRad;
        //for (Vertex v : mesh->vertices()) {
        //    // EG->vertexPositions[v].y *= 2;
        //    double dx1 = VP2[v].x - center1.x;
        //    double dy1 = VP2[v].y - center1.y;
        //    double dz1 = VP2[v].z - center1.z;

        //    double dx2 = VP2[v].x - center2.x;
        //    double dy2 = VP2[v].y - center2.y;
        //    double dz2 = VP2[v].z - center2.z;
        //    if (dx1 * dx1 + dy1 * dy1 + dz1 * dz1 < rad * rad || dx2 * dx2 + dy2 * dy2 + dz2 * dz2 < rad * rad) {
        //        for (Edge e : v.adjacentEdges()) {
        //            tochange[e] = false;
        //        }
        //    }
        //}
        //for (Vertex v : mesh->vertices()) {
        //    for (Edge e : v.adjacentEdges()) {
        //        if (tochange[e]) {
        //            EG->referenceLengths[e] *= 1;
        //            EG->referenceEdgeDihedralAngles[e] *= 1;
        //            tochange[e] = false;
        //        }
        //    }
        //}
        /*EG->vertexPositions[10].x += 0.1;
        EG->vertexPositions[10].y += .1;
        EG->vertexPositions[10].z += .1;*/


        for (Vertex v : mesh->vertices()) {
            // VP[v] += 0.5 * edgeLmin * (randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][0].normalize() +
            // randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][1].normalize());
            EG->vertexPositions[v].y *= 1.001;
            EG->vertexPositions[v].x *= 1.001;
            EG->vertexPositions[v].z *= 1.001;
        }



        EG->refreshQuantities();


       /* int randindex = 10;

        std::cout << "\n  Reference Metric[index]: \n";
        std::cout << EG->referenceMetric[randindex][0] << ", \t";
        std::cout << EG->referenceMetric[randindex][1] << ", \t";
        std::cout << EG->referenceMetric[randindex][2] << "\n";

        std::cout << "\n  Actual Metric[index]: \n";
        std::cout << EG->actualMetric[randindex][0] << ", \t";
        std::cout << EG->actualMetric[randindex][1] << ", \t";
        std::cout << EG->actualMetric[randindex][2] << "\n";


        std::cout << "\n  Reference Curature[index]: \n";
        std::cout << EG->referenceCurvature[randindex][0] << ", \t";
        std::cout << EG->referenceCurvature[randindex][1] << ", \t";
        std::cout << EG->referenceCurvature[randindex][2] << "\n";

        std::cout << "\n  Actual Curvature[index]: \n";
        std::cout << EG->actualCurvature[randindex][0] << ", \t";
        std::cout << EG->actualCurvature[randindex][1] << ", \t";
        std::cout << EG->actualCurvature[randindex][2] << "\n";


        std::cout << "\n  Cauchy Tensor[index]:\n";
        std::cout << EG->elasticCauchyTensor[randindex](0, 0) << ", \t";
        std::cout << EG->elasticCauchyTensor[randindex](0, 1) << ", \t";
        std::cout << EG->elasticCauchyTensor[randindex](0, 2);
        std::cout << "\n";
        std::cout << EG->elasticCauchyTensor[randindex](1, 0) << ", \t";
        std::cout << EG->elasticCauchyTensor[randindex](1, 1) << ", \t";
        std::cout << EG->elasticCauchyTensor[randindex](1, 2);
        std::cout << "\n";
        std::cout << EG->elasticCauchyTensor[randindex](2, 0) << ", \t";
        std::cout << EG->elasticCauchyTensor[randindex](2, 1) << ", \t";
        std::cout << EG->elasticCauchyTensor[randindex](2, 2) << "\n";

        std::cout << "\n  Energy[index]:\n";
        std::cout << EG->elasticEnergy[randindex] << "\n";
        std::cout << "\n  Stretching Energy[index]:\n";
        std::cout << EG->stretchingEnergy[randindex] << "\n";
        std::cout << "\n  Bending Energy[index]:\n";
        std::cout << EG->bendingEnergy[randindex] << "\n";


        std::cout << "\n  \n Number of tirangles:" << mesh->nFaces() << "\n";*/


        EG->computeGradient();
    }
    
    
   

    // Initialize polyscope
    polyscope::init();
    polyscope::options::alwaysRedraw = true;

    // Set the callback function
    //polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), EG->vertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    psMesh->addFaceScalarQuantity("Elastic Energy0", EG->elasticEnergy/EG->faceAreas);
    psMesh->addFaceScalarQuantity("stretch Energy0", EG->stretchingEnergy);
    psMesh->addFaceScalarQuantity("bend Energy0", EG->bendingEnergy);
    psMesh->addEdgeScalarQuantity("reference lengths0", EG->referenceLengths);
    psMesh->addEdgeScalarQuantity("reference angles0", EG->referenceEdgeDihedralAngles);
    psMesh->addVertexVectorQuantity("Grad0", EG->elasticGradient);
    psMesh->addEdgeScalarQuantity("dihedral angles difference0",
                                                          EG->edgeDihedralAngles - EG->referenceEdgeDihedralAngles);
    
    


    //polyscope::state::userCallback = myCallback;  
    
  // polyscope::show();

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
    */



    //polyscope::show();


    mySubroutine();


    return EXIT_SUCCESS;
}