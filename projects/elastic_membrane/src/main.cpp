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
#include "gsl/gsl_multimin.h"



using namespace geometrycentral;
using namespace geometrycentral::surface;
//using namespace gsl_multimin_fdfminimizer;


// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geometry;
std::unique_ptr<ElasticGeometry> EG;
std::unique_ptr<RichSurfaceMeshData> richData;

std::string workingFolder;
int snapshotEvery;
int printinEvery;

double Tinit, Ttarget;
bool restartQ = false;



// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;
polyscope::SurfaceMesh* psMesh0;

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
    geo.thickness = RD.getFaceProperty<double>("Thickness").reinterpretTo(geo.mesh);
    geo.pressure = RD.getFaceProperty<double>("Pressure").reinterpretTo(geo.mesh)[0];
    geo.youngsModulus =
        FaceData<double>(geo.mesh, RD.getFaceProperty<double>("Youngs_Modulus").reinterpretTo(geo.mesh).toVector());
    geo.poissonsRatio =
        FaceData<double>(geo.mesh, RD.getFaceProperty<double>("Poissons_Ratio").reinterpretTo(geo.mesh).toVector());
    // geo.elasticCauchyTensor = RD.getFaceProperty<Eigen::Matrix3f>("Elastic Tensor");

    geo.referenceLengths = RD.getEdgeProperty<double>("Reference_Lengths").reinterpretTo(geo.mesh);
    geo.referenceEdgeDihedralAngles = RD.getEdgeProperty<double>("Reference_Dihedral_Angles").reinterpretTo(geo.mesh);
    // geo.referenceMetric = RD.getFaceProperty<Eigen::Vector3f>("Reference Metric");
    // geo.referenceCurvature = RD.getFaceProperty<Eigen::Vector3f>("Reference Curvature");

    geo.edgeLengths = RD.getEdgeProperty<double>("Actual_Lengths").reinterpretTo(geo.mesh);
    geo.edgeDihedralAngles = RD.getEdgeProperty<double>("Actual_Dihedral_Angles").reinterpretTo(geo.mesh);
    // geo.actualMetric = RD.getFaceProperty<Eigen::Vector3f>("Actual Metric");
    // geo.actualCurvature = RD.getFaceProperty<Eigen::Vector3f>("Actual Curvature");

    //geo.elasticGradient = RD.getVertexProperty<Vector3>("Elastic Gradient");

    geo.bendingEnergy = RD.getFaceProperty<double>("Bending_Energy_density").reinterpretTo(geo.mesh);
    geo.stretchingEnergy = RD.getFaceProperty<double>("Stretching_Energy_density").reinterpretTo(geo.mesh);
    geo.elasticEnergy = RD.getFaceProperty<double>("Elastic_Energy").reinterpretTo(geo.mesh);


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
        if (diffnorm / (oldnorm + 1e-6) > 3) 
            stable = false;
    }
    return stable;
}

void mySubroutine2() {
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
    int max_steps = 300000;
    int count = 0;
    int reg_count = 0;
    int ss_count = 0;
    double last_ener = EG->elasticEnergy.toVector().sum();
    EG->requireVertexDualAreas();
    EG->requireVertexNormals();
    EG->requireVertexMeanCurvatures();
    EG->requireFaceAreas();
    EG->computeGradient();
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
    int press_reg_steps = 1000;
    stepsize = .0001;
    reg_count = 0;
    double pres_func = 0;
    int thickness_reg_steps = 1000;
    double thick_func;  
    double gradnorm;
    VertexData<Vector3> grad = VertexData<Vector3>(EG->mesh, {0, 0, 0});
    VertexData<Vector3> descent_dir=VertexData<Vector3>(EG->mesh, {0, 0, 0});
    double f,g1,g2,alphastep;
    double eps = 1e-6;
    FaceData<double> Ttarget = EG->thickness;
    FaceData<double> Tinit = FaceData<double>(EG->mesh, 1* EG->thickness.toVector().maxCoeff());
    FaceData<double> DeltaT = FaceData<double>(EG->mesh,0);
    for (Face f : EG->mesh.faces()) {
        DeltaT[f] = (Ttarget[f]-Tinit[f])/thickness_reg_steps;
        //std::cout << DeltaT[f]<< ",";
    }
    
    /*for (int reg_step = 1; reg_step <= press_reg_steps; reg_step++) {
     std::cout << "\n \n  Current reg_step: " << reg_step << "/" << press_reg_steps << "\n \n";*/
     do {
        count++;
        reg_count++;
        //if (count % 50 == 0) polyscope::show();
        if (stepsize <= 10 && count % 10 == 1 /* && count > press_reg_steps*/) stepsize *= 2;
        if (restartQ) {
            pres_func = EG->pressure;
            thick_func = EG->thickness[0];
        } else {  // A SHIDDY way to do it. MUST CLEAN THIS CODE
            pres_func = std::min(count * 1.0, press_reg_steps * 1.0) / press_reg_steps / 1.0 * EG->pressure;
        }
        // std::cout << pres_func << "\n";
        stepsize_updated_flag = false;
        // std::cout << count << "\n";
        
        /* Given the gradient at a point g , we define p= g/norm(g). define the line about x along p as  h(a)= f(x+a p).
        To leading order - h(a) = f(x)+ a*p.g + 0.5 * a^2 (p.g)^2 . Which implies a minimum at a= -pg/pg^2.  This is probably too unstable. 
        first pg migh be close to 0 and we need to regularize. so a = -pg /(pg^2 + eps).  Mor important, though, is that this seems very naive. 

        Maybe a better option would be to calculate p, and then g1 = g(x + eps p), g2= g(x - eps p) and from here g' = g1-g2/2eps. Then a = -pg/(g'+eps) 
        */
        gradnorm = 0;
        EG->computeGradient();
        for (Vertex v : mesh->vertices()) {
            grad[v] = -1.0 * EG->elasticGradient[v]; //grad is returned as steepest descent direction, it officaially the steepest ascent
            descent_dir[v] = -1.0 *  EG->elasticGradient[v];
            gradnorm += grad[v].norm2();
        }
        gradnorm = std::sqrt(gradnorm);
        //for (Vertex v : mesh->vertices()) {
        //    descent_dir[v] = descent_dir[v] / gradnorm;
        //}

        ////calculate g1,g2 and f2 =g1-g2/2eps
        //g1 = 0;
        //g2 = 0;
        //EG->vertexPositions += 2 * eps * descent_dir;//start from max   f(x+2eps p)
        //EG->refreshQuantities();
        //f = EG->totalEnergy.toVector().sum();
        //g2 += f;
        //EG->vertexPositions -= eps * descent_dir; //move back once f(x + eps p)
        //EG->refreshQuantities();
        //f = EG->totalEnergy.
        //g2 -= 2 * f;
        //g1 += f;
        //        
        //EG->vertexPositions -= 2 * eps * descent_dir; // f(x- eps p)
        //EG->refreshQuantities();
        //f = EG->totalEnergy.toVector().sum();
        //g2 -= 2 * f;
        //g1 -= f;
        //EG->vertexPositions -= eps * descent_dir; // f(x -2 eps p)
        //EG->refreshQuantities();
        //f = EG->totalEnergy.toVector().sum();
        //g2 += f;
        //EG->vertexPositions += 2* eps * descent_dir; // reset to f(x);
        //EG->refreshQuantities();
        //f = EG->totalEnergy.toVector().sum();
        //g2 += 2 * f;
        //g1 = g1 / eps / 2.0;
        //g2 = g2 / eps /eps / 2.0;
        //alphastep = - g1 / g2;        
         

            for (Vertex v : mesh->vertices()) {
                forces[v] = stepsize * (EG->elasticGradient[v] +
                                        0 * pres_func * EG->vertexNormals[v] * EG->vertexDualAreas[v] / 3);
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


        //update position according to best step;

        //EG->vertexPositions += alphastep * descent_dir;
        if (count <= thickness_reg_steps) for (Face f : mesh->faces()) EG->thickness[f] = std::max(EG->thickness[f]+ DeltaT[f], Ttarget[f]);

        EG->refreshQuantities();
        EG->computeGradient();
        

        /*psMesh->addVertexVectorQuantity("-Grad", 1. * EG->elasticGradient);
        polyscope::show();*/
        gradAmp = 0;
        height = 0;
        for (Vertex v : mesh->vertices()) gradAmp += EG->elasticGradient[v].norm2();
        for (Vertex v : mesh->vertices()) height = std::max(height, EG->vertexPositions[v].y);
        last_ener = tot_ener;
        tot_ener = EG->totalEnergy.toVector().sum();
        if (!stepsize_updated_flag) {
            d_ener = (tot_ener - last_ener)/stepsize;
            rel_ener = d_ener / last_ener;
        }
            
        if ((count - 1) % printinEvery == 0) {
            data.assign({std::to_string(count), 
                            std::to_string(std::time(stopwatch)),
                            std::to_string(gradnorm), // std::to_string(gradnorm),
                            std::to_string(tot_ener),
                            std::to_string(d_ener), 
                            std::to_string(std::abs(gradnorm/tot_ener)),//std::to_string(rel_ener),
                            std::to_string(mean_func(EG->vertexMeanCurvatures / EG->vertexDualAreas)),
                            std::to_string(pres_func),
                            std::to_string(alphastep),//std::to_string(stepsize),
                            std::to_string(mesh->nFaces())});
            printline(headers, data);
            saveLog(headers, data, Log_folder, log_file, true);
        }

        if ((count - 1) % snapshotEvery == 0) {
            //polyscope::view::resetCameraToHomeView();
            polyscope::refresh();
            std::string file1 =
                screen_folder + "energy_" + std::to_string(count) + "order_" + std::to_string(ss_count) + ".png";
            std::string file2 =
                screen_folder + "curvature_" + std::to_string(count) + "order_" + std::to_string(ss_count) + ".png";

            datafile_w =
                screen_folder + "RichData_" + std::to_string(count) + "order_" + std::to_string(ss_count) + ".ply";

            psMesh->updateVertexPositions(EG->vertexPositions);
            psMesh->setAllQuantitiesEnabled(false);
            auto energy_int = psMesh->addFaceScalarQuantity("energy integrated elastic", EG->elasticEnergy);
            /*auto energy_int_log =
                psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy.toVector().unaryExpr(&logfunc));*/
            auto energy_quantity = psMesh->addFaceScalarQuantity(
                "energy content eleastic", (EG->elasticEnergy / EG->faceAreas)); //.toVector().unaryExpr(&logfunc));
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

            auto stretch_ener = psMesh->addFaceScalarQuantity("energy content stretching", EG->stretchingEnergy);
            auto BEND_ener = psMesh->addFaceScalarQuantity("energy content bending", EG->bendingEnergy);
            auto gradient_quantity = psMesh->addVertexVectorQuantity("gradient", EG->elasticGradient);
            gradient_quantity->draw();
            auto normals_quantity = psMesh->addVertexVectorQuantity("face normal", EG->vertexNormals);
            normals_quantity->draw();
            auto angles_quant_vis = psMesh->addEdgeScalarQuantity(
                "dihedral angles difference", EG->edgeDihedralAngles - EG->referenceEdgeDihedralAngles);
            auto ref_angles_vis =
                psMesh->addEdgeScalarQuantity("dihedral angles reference", EG->referenceEdgeDihedralAngles);
            auto angles_vis = psMesh->addEdgeScalarQuantity("dihedral angles actual", EG->edgeDihedralAngles);
            // auto BEND_ener2 = psMesh->addFaceScalarQuantity("bending Energy", EG->bendingEnergy);

            // if(d_ener>0)  polyscope::show();
            auto curv_quantity =
                psMesh->addVertexScalarQuantity("cruvature vertex mean", EG->vertexMeanCurvatures / EG->vertexDualAreas);
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
            //polyscope::show();
        }
     } while ((std::abs(gradnorm / (tot_ener + 1e-4)) > 1e-6 && count <= max_steps && stepsize > 1e-9) ||
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

void assignVector(const gsl_vector* vertexpositions) {
    for (Vertex v: mesh->vertices()) {
        for (int direction = 0; direction < 3; direction++) {
            EG->vertexPositions[v][direction] = gsl_vector_get(vertexpositions, 3*v.getIndex() + direction);
        }
    }
    EG->refreshQuantities();
    EG->computeGradient();
}

void getVector(gsl_vector* vertexpositions) {
    for (Vertex v : mesh->vertices()) {
        for (int direction = 0; direction < 3; direction++) {
            gsl_vector_set(vertexpositions,3* v.getIndex() + direction,  EG->vertexPositions[v][direction]);
        }
    }
}
    //energy for minimizer
double my_f(const gsl_vector *vertexpositions, void* params) {
    assignVector(vertexpositions);
    EG->refreshQuantities();
    return EG->totalEnergy.toVector().sum();
    std::cout << "\n Called f \n";
}

void my_df(const gsl_vector* vertexpositions, void* params, gsl_vector* df) {
    gsl_vector* orig;
    orig = gsl_vector_alloc(3. * mesh->nVertices());
    getVector(orig);
    assignVector(vertexpositions); //can we get by without this? wouldn't the simulation always compute both f and df? or at least df only given we computed f?
    EG->computeGradient();
    for (Vertex v : mesh->vertices()) {
        for (int direction = 0; direction < 3; direction++) {
            gsl_vector_set(df, 3*v.getIndex() + direction, -1 * EG->elasticGradient[v][direction]);
        }
    }
    assignVector(orig);
    EG->computeGradient();
    gsl_vector_free(orig);
}

void my_fdf(const gsl_vector* vertexpositions, void* params, double* f, gsl_vector* df) {
    *f = my_f(vertexpositions, params);
    my_df(vertexpositions, params, df);

    std::cout << "\n Called fdf \n";
}

double get_step_estimate() {
    EG->requireVertexGaussianCurvatures();
    EG->requireVertexMaxPrincipalCurvatures();
    EG->requireVertexMinPrincipalCurvatures();
    EG->requireEdgeLengths();
    EG->computeGradient();
    double max_curv_metric = std::sqrt(EG->vertexGaussianCurvatures.toVector().cwiseAbs().maxCoeff());
    double max_curv1 = EG->vertexMaxPrincipalCurvatures.toVector().cwiseAbs().maxCoeff();
    double max_curv2 = EG->vertexMinPrincipalCurvatures.toVector().cwiseAbs().maxCoeff();
    double min_length = EG->edgeLengths.toVector().minCoeff();
    double max_grad_amp = 0; 
    for (Vertex v : mesh->vertices()) max_grad_amp = fmax(max_grad_amp, EG->elasticGradient[v].norm());
    double max_ratio =
        fmax(max_grad_amp / min_length,
             fmax(max_grad_amp * max_curv1, fmax(max_grad_amp * max_curv2, max_grad_amp * max_curv_metric)));
    EG->unrequireVertexGaussianCurvatures();
    EG->unrequireVertexMaxPrincipalCurvatures();
    EG->unrequireVertexMinPrincipalCurvatures();
    return 0.001/max_ratio;
}
int regulate_grad(VertexData<Vector3>* grad_vec) {
    EG->computeGradient();
    for (Vertex v : mesh->vertices()) (*grad_vec)[v] = EG->elasticGradient[v];
    VertexData<double> gradamps = VertexData<double>(*mesh, 0);
    for (Vertex v : mesh->vertices()) gradamps[v] = (*grad_vec)[v].norm();

    double amp_mean = gradamps.toVector().mean();
    double maxC = gradamps.toVector().maxCoeff();
    double amp_std = 0;
    Vertex maxV;
    for (Vertex v : mesh->vertices())
        if (gradamps[v] == maxC) maxV = v;
    bool is_large = gradamps.toVector().maxCoeff() > 2 * amp_mean;
    bool is_oscilating = false; 
    double oscilation_thresh = 0.3;
    for (Vertex v : maxV.adjacentVertices())
        is_oscilating = (is_oscilating || dot(EG->elasticGradient[maxV], EG->elasticGradient[v]) / gradamps[maxV] <
                                              -1. * oscilation_thresh * gradamps[maxV]);

    /*for (Vertex v : mesh->vertices()) amp_std += (gradamps[v] - am_mean) * (gradamps[v] - amp_mean);
    amp_std = amp_std / mesh->nVertices();
    for (Vertex v : mesh->vertices())
        if (gradamps[v] > amp_mean + 3 * amp_std) is_large = true;*/
    if (is_large) {
        /*for (Vertex v : mesh->vertices()) (*grad_vec)[v] = (*grad_vec)[v] / 2;
        std::cout << "\n reducing  grad \n";*/
        return 1;
    } else if (is_oscilating) {
        return 2;
    }    
    return 0;
}
double stability_check(double fac, int grad_stab, int* counter) {
    double rescale = 1;    
    if (grad_stab == 10) {
        rescale = 1; // 0.5;
        *counter = 0;
        std::cout << "\n uneven gradient, step size reduced by factor of " << rescale << "\n ";
    } else if (grad_stab == 2) {
        rescale = 0.5;
        *counter = 0;
        std::cout << "\n oscillating gradient, step size reduced by factor of " << rescale << "\n ";
    } else if (*counter > 20 && fac < .1) {
        *counter = 0;
        rescale = 2;
        std::cout << "\n stable, stepsize increased by factor of " << rescale << "\n ";
    }
    return rescale * fac;
}
int is_stop(double last_ener, double delta_ener) {
    double gradmean = 0;
    for (Vertex ver : mesh->vertices()) gradmean += EG->elasticGradient[ver].norm() / mesh->nVertices();
    if (gradmean < 1.e-6  || abs(delta_ener/last_ener) < 1.e-6) return 0;
    return -2;
}
int iterate_once() {
    return 0;
}

int ShowPolyscope(int snap = 0, std::string file = "") {
    EG->requireVertexDualAreas();
    EG->requireVertexMeanCurvatures();
    EG->requireElasticEnergy();
   
    polyscope::refresh();

    psMesh->updateVertexPositions(EG->vertexPositions);
    psMesh->setAllQuantitiesEnabled(false);
    psMesh->addFaceScalarQuantity("energy total integrated", EG->totalEnergy);
    auto energy_int = psMesh->addFaceScalarQuantity("energy elastic integrated", EG->elasticEnergy);
    auto energy_quantity = psMesh->addFaceScalarQuantity(
        "energy denisty elastic", (EG->elasticEnergy / EG->faceAreas)); //.toVector().unaryExpr(&logfunc));
    energy_quantity->setMapRange(std::make_pair(0, 0.005)); //(-5, -3)); // EG->elasticEnergy.toVector().maxCoeff()));
    energy_quantity->setColorMap("coolwarm");
    energy_quantity->setEnabled(false);
    energy_quantity->draw();

    auto stretch_ener = psMesh->addFaceScalarQuantity("energy density stretching", EG->stretchingEnergy);
    auto bending_ener = psMesh->addFaceScalarQuantity("energy density bending", EG->bendingEnergy);
    bending_ener->draw();
    bending_ener->setEnabled(true);
    auto gradient_quantity = psMesh->addVertexVectorQuantity("gradient", EG->elasticGradient);
    gradient_quantity->draw();
    gradient_quantity->setEnabled(true);
    auto normals_quantity = psMesh->addFaceVectorQuantity("normal faces", EG->faceNormals);
    normals_quantity->draw();
    auto angles_quant_vis = psMesh->addEdgeScalarQuantity("dihedral angles difference",
                                                          EG->edgeDihedralAngles - EG->referenceEdgeDihedralAngles);
    auto ref_angles_vis = psMesh->addEdgeScalarQuantity("dihedral angles reference ", EG->referenceEdgeDihedralAngles);
    auto angles_vis = psMesh->addEdgeScalarQuantity("dihedral angles actual", EG->edgeDihedralAngles);
    
    auto curv_quantity =
        psMesh->addVertexScalarQuantity("curvature mean vertex", EG->vertexMeanCurvatures / EG->vertexDualAreas);
    curv_quantity->setMapRange(std::make_pair(0, .5)); // EG->elasticEnergy.toVector().maxCoeff()));
    curv_quantity->setColorMap("viridis");
    curv_quantity->setEnabled(false);
    curv_quantity->draw();
    FaceData<double> mean_face_bar = FaceData<double>(*mesh, -100);
    FaceData<double> det_face_bar = FaceData<double>(*mesh, -100);
    FaceData<double> mean_face = FaceData<double>(*mesh, -100);
    FaceData<double> det_face = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_diff_b11 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_diff_b22 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_diff_b12 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_diff_a11 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_diff_a22 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_diff_a12 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_b11 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_b22 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_b12 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_a11 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_a22 = FaceData<double>(*mesh, -100);
    FaceData<double> elemnt_a12 = FaceData<double>(*mesh, -100);
    for (Face f : mesh->faces()) {
        mean_face_bar[f] = EG->referenceCurvature[f][0] + EG->referenceCurvature[f][1];
        mean_face[f] = EG->actualCurvature[f][0] + EG->actualCurvature[f][1];
        det_face_bar[f] = EG->referenceCurvature[f][0] * EG->referenceCurvature[f][1] -
                          EG->referenceCurvature[f][2] * EG->referenceCurvature[f][2];
        det_face[f] = EG->actualCurvature[f][0] * EG->actualCurvature[f][1] -
                      EG->actualCurvature[f][2] * EG->actualCurvature[f][2];
        elemnt_diff_b11[f] = EG->actualCurvature[f][0] - EG->referenceCurvature[f][0];
        elemnt_diff_b22[f] = EG->actualCurvature[f][1] - EG->referenceCurvature[f][1];
        elemnt_diff_b12[f] = EG->actualCurvature[f][2] - EG->referenceCurvature[f][2];
        elemnt_diff_a11[f] = EG->actualMetric[f][0] - EG->referenceMetric[f][0];
        elemnt_diff_a22[f] = EG->actualMetric[f][1] - EG->referenceMetric[f][1];
        elemnt_diff_a12[f] = EG->actualMetric[f][2] - EG->referenceMetric[f][2];

        elemnt_b11[f] = EG->actualCurvature[f][0];
        elemnt_b22[f] = EG->actualCurvature[f][1];
        elemnt_b12[f] = EG->actualCurvature[f][2];
        elemnt_a11[f] = EG->actualMetric[f][0];
        elemnt_a22[f] = EG->actualMetric[f][1];
        elemnt_a12[f] = EG->actualMetric[f][2];
    }


    auto curvefaceref = psMesh->addFaceScalarQuantity("curvature mean face reference", mean_face_bar);    
    auto curvefaceact = psMesh->addFaceScalarQuantity("curvature mean face actual", mean_face);    
    auto curvediff = psMesh->addFaceScalarQuantity("curvature mean face actual", mean_face);    
    auto detfaceref = psMesh->addFaceScalarQuantity("curvature det face reference", det_face_bar);    
    auto detfaceact = psMesh->addFaceScalarQuantity("curvature det face actual", det_face);  
    psMesh->addFaceScalarQuantity("curvature b11 diff", elemnt_diff_b11);
    psMesh->addFaceScalarQuantity("curvature b22 diff", elemnt_diff_b22);
    psMesh->addFaceScalarQuantity("curvature b12 diff", elemnt_diff_b12);
    psMesh->addFaceScalarQuantity("curvature a11 diff", elemnt_diff_a11);
    psMesh->addFaceScalarQuantity("curvature a22 diff", elemnt_diff_a22);
    psMesh->addFaceScalarQuantity("curvature a12 diff", elemnt_diff_a12);

    psMesh->addFaceScalarQuantity("curvature b11 elem", elemnt_b11);
    psMesh->addFaceScalarQuantity("curvature b22 elem", elemnt_b22);
    psMesh->addFaceScalarQuantity("curvature b12 elem", elemnt_b12);
    psMesh->addFaceScalarQuantity("curvature a11 elem", elemnt_a11);
    psMesh->addFaceScalarQuantity("curvature a22 elem", elemnt_a22);
    psMesh->addFaceScalarQuantity("curvature a12 elem", elemnt_a12);
    


    psMesh->addEdgeScalarQuantity("lengths reference", EG->referenceLengths);
    psMesh->addEdgeScalarQuantity("lengths actual", EG->edgeLengths);
    //    psMesh->addEdgeScalarQuantity("reference angles0", EG->referenceEdgeDihedralAngles);
    EG->requireVertexDualAreas();
    EG->requireVertexNormals();
    EG->requireFaceVolume();
    EG->requireTotalEnergy();
    EG->computeGradient();
    psMesh->addVertexScalarQuantity("area vertex", EG->vertexDualAreas);
    // psMesh->addVertexVectorQuantity("-Grad_Norm",-1 * EG->elasticGradient / EG->vertexDualAreas);
    psMesh->addVertexVectorQuantity("normal vertexes", EG->vertexNormals * EG->pressure);
    //  psMesh->addVertexVectorQuantity("Forces", EG->vertexNormals * EG->vertexDualAreas * EG->pressure);
    
    psMesh->addFaceScalarQuantity("area face", EG->faceAreas);
    EG->requireEdgeDihedralAngles();
    psMesh->addFaceScalarQuantity("volume faces", EG->faceVolume);
    VertexData<Vector3> deformation_vector = VertexData<Vector3>(*mesh, {0, 0, 0});
    for (Vertex v : mesh->vertices()) {
        deformation_vector[v] = -EG->vertexPositions[v] + geometry->vertexPositions[v];
    }
    psMesh->addVertexVectorQuantity("deformation", deformation_vector);

    polyscope::screenshot(file, true);
    

    
    if (snap > 0) return 0;
    else polyscope::show();
    return 0;
}

int mySubroutine3() { // REWRITE. WITHOUT GSL as this seems to fail after gfixing energ terms.
    // Either implement your own gradient descent algorithm, or do as follows (less eficient in terms of run time):
    //  1. For large gradient gradients (large changes of the gradianet along the surface), limit max gradient amp.
    //  (possibly through grad amp statistics)
    //  2. Max dx (u,v) step amp cant be larger than (somefraction) of minimal curvauture scale / position amp.
    //  3. Stopping condition (grad <= 1e-4  elastic energy scale) (relevant only for elastic problems under pressure)
    //  4. Stability test - several large gradient oscillations (spatial + temporal) -> step size reduced. Lack of ->
    //  step size increased up to max.

    // test if everything works do far (get and assignment of vectors and grads)
    size_t iter = 0;
    int status = -2;

    
    std::cout << "\n Initializing...";
    
    double pressuremax = EG->pressure;
    FaceData<double> StopThick = FaceData<double>(*EG->thickness.getMesh(), EG->thickness.toVector());
    FaceData<double> StartThick = FaceData<double>(*EG->thickness.getMesh(), 1 * EG->thickness.toVector());
    for (Face f : mesh->faces()) EG->thickness[f] = EG->thickness[f] = StartThick[f];
    EG->refreshQuantities();
    EG->computeGradient();
    
    std::cout << "\n Initial energy: " << EG->totalEnergy.toVector().sum();

    std::cout << "\n \n ********************************* \n Begin minimzation \n "
                 "********************************* \n \n";

    double fac = get_step_estimate();
    std::cout << "\n Initial step size: " << fac;
    
    int press_reg_iter = 0;
    int thick_reg_iter = 0;
    

    double last_ener = 0;
    double delta_ener = 0;

    int ss_count = 0;
    double normi = 0;

    std::time_t* stopwatch = 0;


    std::vector<std::string> headers = {"Iteration",
                                        "Epoch",
                                        "Gradient",
                                        "Step Size",
                                        "Total Energy",
                                        "Elastic Energy",
                                        "dEnergy(this step)",
                                        "Average Mean Curvature",
                                        "pressure",
                                        "Thickness",
                                        "Number of triangles",
                                        "Status"};
    std::string screen_folder = workingFolder;
    std::string Log_folder = workingFolder;
    std::string log_file = "log.log";
    std::string datafile_w;
    std::vector<std::string> data;
    saveLog(headers, data, Log_folder, log_file,true);

    std::cout << "\n pressure, thickness = " << pressuremax << "," << StopThick[0] << "  " << EG->pressure << ","
              << EG->thickness[12] << "\n";

    
    double test_Energy;
    double current_energy = EG->totalEnergy.toVector().sum();
    double delta_energy = 0.;   
    bool reduce_flag;
    VertexData<Vector3> grad_vec = VertexData<Vector3>(*mesh, {0., 0., 0.});
    VertexData<Vector3> last_pos = VertexData<Vector3>(*mesh, {0., 0., 0.});
    for (Vertex v : mesh->vertices()) last_pos[v] = EG->vertexPositions[v];
    fac = get_step_estimate();
    int grad_stab;
    int reduced_counter = 0;
    int succesive_non_reductions = 0;
    do {
        iter++;
        if (iter <= press_reg_iter) {
            EG->pressure = 1.0 * (1 + 0 * iter / press_reg_iter) * pressuremax;

        } else if (iter - press_reg_iter <= thick_reg_iter) {
            for (Face f : mesh->faces())
                EG->thickness[f] = StopThick[f] + 0 * (StartThick[f] - (StartThick[f] - StopThick[f]) *
                                                                           (iter - press_reg_iter) / (thick_reg_iter));
        }
        EG->refreshQuantities();
        EG->computeGradient();        
        grad_stab = regulate_grad(&grad_vec);
        fac = stability_check(fac, grad_stab, &succesive_non_reductions);
        for (Vertex v : mesh->vertices()) {
            EG->vertexPositions[v] += fac*grad_vec[v];
        }
        EG->refreshQuantities();
        test_Energy = EG->totalEnergy.toVector().sum();
        if (test_Energy <= current_energy) {
            for (Vertex v : mesh->vertices()) last_pos[v] = EG->vertexPositions[v];
            delta_energy = test_Energy - current_energy;
            current_energy = test_Energy;
            if (test_Energy == current_energy)
                reduced_counter = fmax(reduced_counter - 1, 0);
            succesive_non_reductions++;
            //ShowPolyscope();

        } else /*if(iter>100)*/ { // roleback
            iter--; 
            std::cout << "\n energy unstable, rollback and step reduced\n";
            for (Vertex v : mesh->vertices()) EG->vertexPositions[v] = last_pos[v];
            EG->refreshQuantities();
            EG->computeGradient();
            reduced_counter++;
            succesive_non_reductions = 0;
            fac *=0.5;
            if (reduced_counter == 100 || fac< 1e-6) status = 27;                          
           // ShowPolyscope();
        }
        if (status == 27) {
            std::cout << "\n Does not converge!"; 
            break;
        }
        status = is_stop(current_energy, delta_energy/fac);    
        
        if (iter % 100 == 0 && false) {

            EG->refreshQuantities();
            EG->computeGradient();
            psMesh->updateVertexPositions(EG->vertexPositions);
            auto gradient_quantity = psMesh->addVertexVectorQuantity("Gradient", EG->elasticGradient);
            gradient_quantity->draw();
            auto angles_vis = psMesh->addEdgeScalarQuantity("actual dihedral angles", EG->edgeDihedralAngles);
            auto angles_quant_vis = psMesh->addEdgeScalarQuantity(
                "dihedral angles difference", EG->edgeDihedralAngles - EG->referenceEdgeDihedralAngles);
            auto ref_angles_vis =
                psMesh->addEdgeScalarQuantity("reference dihedral angles", EG->referenceEdgeDihedralAngles);

            auto stretch_ener = psMesh->addFaceScalarQuantity("Stretching Energy density", EG->stretchingEnergy);
            auto BEND_ener = psMesh->addFaceScalarQuantity("bending Energy density", EG->bendingEnergy);
            // auto gradient_quantity = psMesh->addVertexVectorQuantity("Gradient", EG->elasticGradient);
            // gradient_quantity->draw();
            auto normals_quantity = psMesh->addVertexVectorQuantity("Normal (pressur)", EG->vertexNormals);
            normals_quantity->draw();

            // auto curv_quantity =
            //     psMesh->addVertexScalarQuantity("MeanCurvature", EG->vertexMeanCurvatures /
            //     EG->vertexDualAreas);
            // curv_quantity->setMapRange(std::make_pair(0, .5)); // EG->elasticEnergy.toVector().maxCoeff()));
            // curv_quantity->setColorMap("viridis");
            // curv_quantity->setEnabled(true);
            // curv_quantity->draw();


            psMesh->refresh();
            polyscope::show();
        }
        
        if (status == 0) {
            printf("\n Minimum found!\n");
            break;
        }

       

        if ((iter - 1) % printinEvery == 0) {
            EG->requireVertexDualAreas();
            EG->requireVertexMeanCurvatures();
            EG->requireElasticEnergy();
            for (Vertex ver : mesh->vertices()) normi += EG->elasticGradient[ver].norm() / mesh->nVertices();
            data.assign({std::to_string(iter), std::to_string(std::time(stopwatch)),
                         std::to_string(normi), // std::to_string(gradnorm),
                         std::to_string(fac),
                         std::to_string(/* s->f*/ EG->totalEnergy.toVector().sum()),
                         std::to_string(EG->elasticEnergy.toVector().sum()), std::to_string(delta_energy/fac),
                         std::to_string(mean_func(EG->vertexMeanCurvatures / EG->vertexDualAreas)),
                         std::to_string(EG->pressure),
                         std::to_string(EG->thickness[1]), // std::to_string(stepsize),
                         std::to_string(mesh->nFaces()), std::to_string(status)});
            printline(headers, data);
            saveLog(headers, data, Log_folder, log_file, true);
            normi = 0;
        }
        if ((iter - 1) % snapshotEvery == 0 ) {
            // polyscope::view::resetCameraToHomeView();
            polyscope::refresh();
            std::string file1 =
                screen_folder + "energy_" + std::to_string(iter) + "order_" + std::to_string(ss_count) + ".png";
            std::string file2 =
                screen_folder + "curvature_" + std::to_string(iter) + "order_" + std::to_string(ss_count) + ".png";

            datafile_w =
                screen_folder + "RichData_" + std::to_string(iter) + "order_" + std::to_string(ss_count) + ".ply";

            ShowPolyscope(1, file1);

            writeRichData(*richData, *EG, datafile_w);
            polyscope::screenshot(file2, true);
            ss_count += 1;
        }
    } while (status == -2  && iter < 5000);

    std::cout << "Status: " << status << "\n";

    


    std::cout << "\n  \n \t \t SIM COMPLETE! \n \n";
    if (iter >= 5000) std::cout << "max interation exceeded \n";
    if (status == 27) std::cout << "unstable \n";
    
    std::cout << "Total iterations: " << iter << "\n ";
    std::cout << "Final log: \n";
    EG->requireVertexDualAreas();
    EG->requireVertexMeanCurvatures();
    EG->requireElasticEnergy();
    for (Vertex ver : mesh->vertices()) normi += EG->elasticGradient[ver].norm() / mesh->nVertices();
    data.assign({
        std::to_string(iter), std::to_string(std::time(stopwatch)),
            std::to_string(normi), // std::to_string(gradnorm),
            std::to_string(fac),
            std::to_string(EG->totalEnergy.toVector().sum()), std::to_string(EG->elasticEnergy.toVector().sum()),
            std::to_string(delta_energy/fac),
                 std::to_string(mean_func(EG->vertexMeanCurvatures / EG->vertexDualAreas)),
                 std::to_string(EG->pressure),
                 std::to_string(EG->thickness[1]), // std::to_string(stepsize),
                 std::to_string(mesh->nFaces()), std::to_string(status)});
    printline(headers, data);
    saveLog(headers, data, Log_folder, log_file, true);

    datafile_w = screen_folder + "RichData_Final.ply";
    writeRichData(*richData, *EG, datafile_w);
    ShowPolyscope(1, screen_folder + "energy_" + std::to_string(iter) + "_Final.png");

    return 0;

    // psMesh->updateVertexPositions(EG->vertexPositions);
    // polyscope::view::resetCameraToHomeView();
    // polyscope::refresh();

    // psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy / EG->faceAreas);
    // psMesh->addFaceScalarQuantity("thickness", EG->thickness);
    // psMesh->addFaceScalarQuantity("stretch Energy", EG->stretchingEnergy);
    // psMesh->addFaceScalarQuantity("bend Energy", EG->bendingEnergy);
    // psMesh->addEdgeScalarQuantity("reference lengths0", EG->referenceLengths);
    // psMesh->addEdgeScalarQuantity("actual lengths", EG->edgeLengths);
    // psMesh->addEdgeScalarQuantity("lengths difference", EG->edgeLengths - EG->referenceLengths);
    ////    psMesh->addEdgeScalarQuantity("reference angles0", EG->referenceEdgeDihedralAngles);
    // EG->requireVertexDualAreas();
    // EG->requireVertexNormals();
    // EG->requireFaceVolume();
    // EG->requireTotalEnergy();
    // psMesh->addVertexVectorQuantity("-Grad", 1. * EG->elasticGradient);
    // psMesh->addVertexScalarQuantity("vertex Area", EG->vertexDualAreas);
    //// psMesh->addVertexVectorQuantity("-Grad_Norm",-1 * EG->elasticGradient / EG->vertexDualAreas);
    // psMesh->addVertexVectorQuantity("Vertex Normal", EG->vertexNormals * EG->pressure);
    ////  psMesh->addVertexVectorQuantity("Forces", EG->vertexNormals * EG->vertexDualAreas * EG->pressure);
    // psMesh->addFaceScalarQuantity("Total Energy", EG->totalEnergy);
    // psMesh->addFaceScalarQuantity("faceVolume", EG->faceVolume);
    // auto angles_quant_vis = psMesh->addEdgeScalarQuantity("dihedral angles difference",
    //                                                       EG->edgeDihedralAngles -
    //                                                       EG->referenceEdgeDihedralAngles);
    // auto ref_angles_vis = psMesh->addEdgeScalarQuantity("reference dihedral angles",
    // EG->referenceEdgeDihedralAngles); auto angles_vis = psMesh->addEdgeScalarQuantity("actual dihedral angles",
    // EG->edgeDihedralAngles);
}

int mySubroutine() {  // REWRITE. WITHOUT GSL as this seems to fail after gfixing energ terms.  
    //Either implement your own gradient descent algorithm, or do as follows (less eficient in terms of run time): 
    // 1. For large gradient gradients (large changes of the gradianet along the surface), limit max gradient amp. (possibly through grad amp statistics)
    // 2. Max dx (u,v) step amp cant be larger than (somefraction) of minimal curvauture scale / position amp.
    // 3. Stopping condition (grad <= 1e-4  elastic energy scale) (relevant only for elastic problems under pressure)
    // 4. Stability test - several large gradient oscillations (spatial + temporal) -> step size reduced. Lack of -> step size increased up to max.

    // test if everything works do far (get and assignment of vectors and grads)
    size_t iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type* T;
    gsl_multimin_fdfminimizer* s;

    gsl_vector* x;
    gsl_multimin_function_fdf my_func;

    int vec_size = 3 * mesh->nVertices(); // we are flattening the position vevtor.

    my_func.n = vec_size;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = 0; // no function params, real function is implemented in EG object

    x = gsl_vector_alloc(vec_size);
    getVector(x);

    std::cout << "\n (f,energy) = (" << my_func.f(x, 0) << "," << EG->totalEnergy.toVector().sum() << ")";
    //T = gsl_multimin_fdfminimizer_conjugate_fr;
    T = gsl_multimin_fdfminimizer_conjugate_pr;
   //  T = gsl_multimin_fdfminimizer_vector_bfgs2;
    // T = gsl_multimin_fdfminimizer_steepest_descent;
    s = gsl_multimin_fdfminimizer_alloc(T, vec_size);

    std::cout << "\n f before init: " << my_func.f(x, 0);
    double pressuremax = EG->pressure;
    FaceData<double> StopThick = FaceData<double>(*EG->thickness.getMesh(), EG->thickness.toVector());
    FaceData<double> StartThick = FaceData<double>(*EG->thickness.getMesh(),  3*EG->thickness.toVector());
    for (Face f : mesh->faces()) EG->thickness[f] = EG->thickness[f] = StartThick[f];
    EG->refreshQuantities();
    EG->computeGradient();

    std::cout << "\n \n ********************************* \n Begin minimzation \n "
                 "********************************* \n \n";
    gsl_vector* dir = gsl_vector_alloc(vec_size);
    gsl_vector* grad = gsl_vector_alloc(vec_size);
    my_func.df(x, 0, grad);

    int press_reg_iter = 100;
    int thick_reg_iter = 100;
   


    double last_ener = 0;
    double delta_ener = 0;

    int ss_count = 0;
    double normi = 0;

    std::time_t* stopwatch = 0;


    std::vector<std::string> headers = {"Iteration",
                                        "Epoch",
                                        "Gradient",
                                        "Total Energy",
                                        "Elastic Energy",
                                        "dEnergy(this step)",
                                        "Average Mean Curvature",
                                        "pressure",
                                        "Thickness",
                                        "Number of triangles",
                                        "Status"};
    std::string screen_folder = workingFolder;
    std::string Log_folder = workingFolder;
    std::string log_file = "log.log";
    std::string datafile_w;
    std::vector<std::string> data;
    saveLog(headers, data, Log_folder, log_file);

     std::cout << "\n pressure, thickness = " << pressuremax << "," << StopThick[0] << "  " << EG->pressure << ","
              << EG->thickness[12] << "\n";

     gsl_multimin_fdfminimizer_set(s, &my_func, x, 1e-6, 1e-6);
     int reset_count = 0;

    ShowPolyscope(1);

    do {
        iter++;
        if (iter % 10 == 0) {

            if (iter <= press_reg_iter) {
                EG->pressure = 1.0 * (0 * 1 + 1. * iter / press_reg_iter) * pressuremax;
                EG->refreshQuantities();
                // EG->computeGradient();

            } else if (iter - press_reg_iter <= thick_reg_iter) {
                for (Face f : mesh->faces())
                    EG->thickness[f] =
                        0. * StopThick[f] + 1. * (StartThick[f] - (StartThick[f] - StopThick[f]) *
                                                                      (iter - press_reg_iter) / (thick_reg_iter));
                EG->refreshQuantities();
                // EG->computeGradient();
            }
        }
        
        
    status = gsl_multimin_fdfminimizer_iterate(s);
    if (status && iter > (press_reg_iter + thick_reg_iter)) break;        
    status = gsl_multimin_test_gradient(s->gradient, 1e-4);
    if (status == GSL_SUCCESS) printf("Minimum found at:\n");
    delta_ener = s->f - last_ener;       
    last_ener = s->f;


        if ((iter - 1) % printinEvery == 0) {

            EG->requireVertexDualAreas();
            EG->requireVertexMeanCurvatures();
            EG->requireElasticEnergy();
            for (Vertex ver : mesh->vertices()) normi += EG->elasticGradient[ver].norm() / mesh->nVertices();
            data.assign({std::to_string(iter), std::to_string(std::time(stopwatch)),
                         std::to_string(normi), // std::to_string(gradnorm),
                         std::to_string(/* s->f*/ EG->totalEnergy.toVector().sum()), std::to_string(EG->elasticEnergy.toVector().sum()),
                         std::to_string(delta_ener),
                         std::to_string(mean_func(EG->vertexMeanCurvatures / EG->vertexDualAreas)),
                         std::to_string(EG->pressure),
                         std::to_string(EG->thickness[1]), // std::to_string(stepsize),
                         std::to_string(mesh->nFaces()), 
                         std::to_string(status)});
            printline(headers, data);
            saveLog(headers, data, Log_folder, log_file, true);
            normi = 0;
        }

        if ((iter - 1) % snapshotEvery == 0) {
            polyscope::refresh();
            std::string file1 =
                screen_folder + "energy_" + std::to_string(iter) + "order_" + std::to_string(ss_count) + ".png";
            std::string file2 =
                screen_folder + "curvature_" + std::to_string(iter) + "order_" + std::to_string(ss_count) + ".png";

            datafile_w =
                screen_folder + "RichData_" + std::to_string(iter) + "order_" + std::to_string(ss_count) + ".ply";

            ShowPolyscope(1,file1);


            writeRichData(*richData, *EG, datafile_w);           
            ss_count += 1;
        }
        if (status == GSL_ENOPROG && reset_count <10) {
            std::cout << "\n reset minimization \n";
            status = GSL_CONTINUE;
            gsl_multimin_fdfminimizer_restart(s);
            reset_count++;
        }

        //double gradmean = 0;
        //for (Vertex ver : mesh->vertices()) gradmean += EG->elasticGradient[ver].norm() / mesh->nVertices();
        //if (gradmean < 1e-4) status = GSL_SUCCESS;
    } while (status == GSL_CONTINUE && iter < 1000);

    

    std::cout << "Status: " << status << "\n";

    assignVector(s->x);

    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);


    std::cout << "\n  \n \t \t SIM COMPLETE! \n \n";
    if (iter > 10000) std::cout << "max interation exceeded \n";
    if (status == 27) std::cout << "Didn't converge \n";
    if (std::abs(delta_ener) <= 1e-6) std::cout << "minimum found! \n";
    std::cout << "Total iterations: " << iter << "\n ";
    std::cout << "Final log: \n";
    EG->requireVertexDualAreas();
    EG->requireVertexMeanCurvatures();
    EG->requireElasticEnergy();
    for (Vertex ver : mesh->vertices()) normi += EG->elasticGradient[ver].norm() / mesh->nVertices();
    data.assign({std::to_string(iter), 
                 std::to_string(std::time(stopwatch)),
                 std::to_string(normi), // std::to_string(gradnorm),
                 std::to_string(s->f),
                 std::to_string(EG->elasticEnergy.toVector().sum()), 
                 std::to_string(delta_ener),
                 std::to_string(mean_func(EG->vertexMeanCurvatures / EG->vertexDualAreas)),
                 std::to_string(EG->pressure),
                 std::to_string(EG->thickness[1]), // std::to_string(stepsize),
                 std::to_string(mesh->nFaces()), 
                 std::to_string(status)});
    printline(headers, data);
    saveLog(headers, data, Log_folder, log_file, true); 

    datafile_w = screen_folder + "RichData_before_Final.ply";
    writeRichData(*richData, *EG, datafile_w);
    ShowPolyscope(1, screen_folder + "energy_" + std::to_string(iter) + "_before_final.png");
    return 0;
}




int main(int argc, char** argv) {

    bool DEBUG_MODE = false;


    // Configure the argument parser
    args::ArgumentParser parser("15-458 HW2");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");
    args::ValueFlag<double> thickness(parser, "thickness", "thickness value (uniform)", {'t', "thickness"}, .1);
    args::ValueFlag<double> pressure(parser, "pressure", "pressure value", {'p', 'P', "pressure"}, 0.1);
    args::ValueFlag<double> Youngs(parser, "Young's Modulus", "Young's Modulus", {'Y', 'E', "youngs"}, 1.);
    args::ValueFlag<double> Poissons(parser, "Poisson's ratio", "Poisson's ratio", {'u', "poisson", "poissons"}, .5);
    args::ValueFlag<std::string> SaveFolderName(parser, "folder name", "Folder to save snapshots and final result",
                                                {'f', "folder"}, "D:/code_output/geometry/debug/");
    args::ValueFlag<int> printingCoutner(
        parser, "print counter", "print meta data after every [print counter] iterations", {'l', "log", "print"}, 10);
    args::ValueFlag<int> snapshotCoutner(parser, "snapshot counter",
                                         "snapshot after every [snapshot counter] iterations",
                                         {'s', "snap", "snapshot"}, 100);
    args::ValueFlag<double> otherVal(parser, "other", "other", {'o', "other"}, -100);
    args::ValueFlag<double> anotherVal(parser, "another", "another", {'k', "another"}, -100);
    /*args::MapFlag<std::string, bool> readfromfileflag(parser, "input", "an input file .ply",
                                                      {'r', "read", "in", "input", "readfile", "inputfile"},
                                 {{"true", true}, {"false", false}});*/


    ////// CALCS
    ////
    //polyscope::init();
    //polyscope::options::alwaysRedraw = true;
    //std::string file1 = "D:/code_output/geometry/Completed/height_2.8_pressure_0.01/RichData_1order_0.ply";
    //std::string file2 = "D:/code_output/geometry/Completed/height_2.8_pressure_0.01/RichData_Final.ply";
    //std::string csvname = "D:/code_output/geometry/Completed/height_2.8_pressure_0.01/justa.csv";
    //std::unique_ptr<SurfaceMesh> sMesh;
    //std::unique_ptr<SurfaceMesh> sMesh2;
    //std::tie(sMesh, richData) = RichSurfaceMeshData::readMeshAndData(file2);
    //geometry = richData->getGeometry();
    //// mesh = std::move(sMesh->toManifoldMesh());
    //// sMesh.reset();
    //EG = std::move(std::unique_ptr<ElasticGeometry>(new ElasticGeometry(*sMesh, geometry->vertexPositions)));
    //EG->vertexPositions = geometry->vertexPositions;
    //readRichaData(*richData, *EG);
    //mesh.reset();
    //EG->refreshQuantities();
    //EG->requireElasticEnergy();
    ////
    //std::tie(sMesh2, richData) = RichSurfaceMeshData::readMeshAndData(file1);
    //geometry = richData->getGeometry();
    //VertexData<Vector3> deformation = VertexData<Vector3>(EG->mesh, Vector3{0., 0., 0.});
    //VertexData<double> deformation_mag = VertexData<double>(EG->mesh, 0.);
    //for (Vertex v : EG->mesh.vertices()) {
    //    deformation[v] = EG->vertexPositions[v] - geometry->vertexPositions[v];
    //    deformation_mag[v] = deformation[v].norm();
    //}
    ////
    //std::ofstream csvfile;
    //csvfile.open(csvname);
    //for (Face f : EG->mesh.faces()) {
    //    double counter = 0;
    //    double ycor = 0;
    //    double ycororiginal = 0;
    //    double dformag = 0;
    //    Vector3 defor = Vector3{0., 0., 0.};
    //    for (Vertex v : f.adjacentVertices()) {
    //        ycor += EG->vertexPositions[v].y;
    //        ycororiginal += geometry->vertexPositions[v].y;
    //        defor.x += deformation[v].x;
    //        defor.y += deformation[v].y;
    //        defor.z += deformation[v].z;
    //        counter += 1;
    //        //
    //    }
    //    ycor = ycor / counter;
    //    ycororiginal = ycororiginal / counter;
    //    dformag = defor.norm() / counter;
    //    //
    //    csvfile << f.getIndex() << "," << ycororiginal << "," << ycor << "," << EG->elasticEnergy[f] / EG->faceArea(f)
    //            << "," << dformag << "\n";
    //}
    //csvfile.close();
    ////
    ////
    //psMesh = polyscope::registerSurfaceMesh("deformed", EG->vertexPositions, EG->mesh.getFaceVertexList(),
    //                                        polyscopePermutations(EG->mesh));
    //psMesh = polyscope::registerSurfaceMesh("base", geometry->vertexPositions, EG->mesh.getFaceVertexList(),
    //                                        polyscopePermutations(EG->mesh));
    //psMesh->addVertexVectorQuantity("deformation", deformation);
    //psMesh->addVertexScalarQuantity("deformation_mag", deformation_mag);
    //EG->requireFaceAreas();
    //psMesh->addFaceScalarQuantity("Elastic Energy", EG->elasticEnergy / EG->faceAreas);
    //polyscope::show();

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

    std::string filepath =
        "C:/Users/dgrossma/Documents/GitHub/ElasticSim/input/torus.obj";//"C:/Users/dgrossma/Documents/GitHub/ElasticSim/input/smooth_cilinder_proto.obj";//"D:/code_output/geometry/fucus_200A+_synth_thickness_1.0_pressure_0.01/RichData_Final.ply";
    ////"D:/code_output/geometry/Completed/height_2.8_pressure_0.01/final/RichData_Final.ply"; ////"D:/code_output/geometry/inputs/Fucus_Synth/fucus_rhiz_synth.obj";//sphere.obj"; //
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
              << "Another: " << anotherVal.Get() << "\n"
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
        restartQ = true;
        std::unique_ptr<SurfaceMesh> sMesh;
        std::unique_ptr<VertexPositionGeometry> tGEO;
        std::tie(mesh, richData) = RichSurfaceMeshData::readManifoldMeshAndData(filepath);            
        tGEO = richData->getGeometry();
        geometry = std::move(tGEO->reinterpretTo(*mesh));

        // centering the object
        Vector3 centerMass = {0, 0, 0};
        for (Vertex v : mesh->vertices()) {
            centerMass += geometry->vertexPositions[v];
        }
        centerMass = centerMass / mesh->nVertices();

        for (Vertex v : mesh->vertices()) {
            geometry->vertexPositions[v] -= centerMass;
        }

        EG = std::move(std::unique_ptr<ElasticGeometry> (new ElasticGeometry(*mesh, geometry->vertexPositions)));
       // EG->vertexPositions = geometry->vertexPositions;
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



         /*psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->vertexPositions,
                                                mesh->getFaceVertexList(), polyscopePermutations(*mesh));
        polyscope::show();*/



        //centering the object
        Vector3 centerMass = {0, 0, 0};
        for (Vertex v : mesh->vertices()) {
            centerMass += geometry->vertexPositions[v];
        }
        centerMass = centerMass / mesh->nVertices();

        for (Vertex v : mesh->vertices()) {
            VP[v].x -= centerMass.x;
            VP[v].z -= centerMass.y;
            VP[v].y -= centerMass.z;
        }



        double edgeLmin = geometry->edgeLengths.toVector().minCoeff(); //wiggle wigglwe
        double maxz = 0, maxy = 0,miny=1e10;
        for (Vertex v : mesh->vertices()) {
            maxy = std::max(maxy, VP[v].y);
            miny = std::min(miny, VP[v].y);
            maxz = std::max(maxz, VP[v].z);
            maxz = std::max(maxz, -VP[v].z);

            
        }
        std::cout << "\n maxz: " << maxz << "\n";

        for (Vertex v : mesh->vertices()) {
             VP[v] += 0 * edgeLmin * (randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][0].normalize() +
             randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][1].normalize());            
        }

        geometry->requireEdgeLengths();
        double geometry_rescale = maxz;
        //geometry->edgeLengths.toVector().mean();
        for (Vertex v : mesh->vertices()) { //rescale            
            VP[v].y *= 1. / geometry_rescale; //* otherVal.Get(); 
            VP[v].x *= 1. / geometry_rescale; // seting radius to 1 (widdth = 2r)
            VP[v].z *= 1. / geometry_rescale; 
        }

         for (Vertex v : mesh->vertices()) {
            VP[v].x += 0;
            VP[v].z += 0;
            VP[v].y += 0;
        }

        std::cout << "\n"
                  << maxy << "," << maxz<< "," << miny ;
        geometry->vertexPositions = VP;
        geometry->refreshQuantities();   


        for (Vertex v : mesh->vertices()) {
            maxy = std::max(maxy, VP[v].y);
            miny = std::min(miny, VP[v].y);
            maxz = std::max(maxz, VP[v].z);
            maxz = std::max(maxz, -VP[v].z);
        }
        std::cout << "\n maxz after rescale: " << maxz << "\n";

       // psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), geometry->vertexPositions,
         //                                       mesh->getFaceVertexList(), polyscopePermutations(*mesh));
      //  polyscope::show();
        

        double scale_factor = geometry_rescale;
        if (otherVal.Get() != otherVal.GetDefault()) scale_factor = otherVal.Get();
        Ttarget = thickness.Get() / scale_factor;
        if (Ttarget > 0.2 * fmin(maxz,maxy)) {
            std::cout << "\n thickness " << Ttarget <<" incompatible with curvature(too thick), forcing a limit.setting thickenss to :";
            //Ttarget =0.2 *fmin(maxz, maxy);
            std::cout << Ttarget << "\n";
        }
         // 12  is the distance normalization
        
        // OK for non-apical cells (everything is just same thickenss for init). We can create from base and run from there. However it is probably better to modify the resutls, save, and run after prep for all.
        Tinit = 5 * Ttarget; 


        if (0.2 * maxz > Tinit) Tinit = 0.2 * maxz;
        //std::unique_ptr<ElasticGeometry> BG(new ElasticGeometry(*mesh, geometry->vertexPositions,Tinit,Youngs.Get(),Poissons.Get(),pressure.Get() ));
        


        // PREP - create   .ply  and 
        std::unique_ptr<ElasticGeometry> BG(new ElasticGeometry(*mesh, geometry->vertexPositions, Ttarget, Youngs.Get(),
                                                                Poissons.Get(), pressure.Get()));
        FaceData<double> thickness = FaceData<double>(*mesh, Ttarget);
        double thickenning_stops = 4;
        double thickenning_begins = maxy-miny;        // 2.0 + thickenning_stops / 2.0;        
        //for(Face f: mesh->faces()) 
        //{
        //    Vector3 faceCenter = {0, 0, 0};
        //    for (Vertex v : f.adjacentVertices())
        //    {
        //        faceCenter.x += 1.0 / 3 * BG->vertexPositions[v].x;
        //        faceCenter.y += 1.0 / 3 * BG->vertexPositions[v].y;
        //        faceCenter.z += 1.0 / 3 * BG->vertexPositions[v].z;
        //    }
        //    double y_pos = faceCenter.y - miny;
        //    if ( y_pos > thickenning_stops) {
        //        thickness[f] = std::min((thickenning_begins - y_pos) / (thickenning_begins - thickenning_stops) * Ttarget +
        //                       0.3 / scale_factor,Ttarget);
        //       /* std::cout << "\n"
        //                  << y_pos << "," << thickness[f] << "," << thickenning_begins << "," << thickenning_stops
        //                  << "," << (thickenning_begins - y_pos) / (thickenning_begins - thickenning_stops) * Ttarget
        //                  << "," << Ttarget;*/
        //    }
        //}
        BG->thickness = thickness;       //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// FINISH IMPLEMENTATION: save to .PLY also change sim process! (non uniform relax)
        BG->refreshQuantities();
        EG = std::move(BG);
        EG->requireElasticEnergy();
        EG->requireBendingEnergy();      
        

       //writeRichData(*richData, *EG, "D:/code_output/geometry/inputs/apical_base.ply");


      // return 0;

        std::cout << "\n" << pressure.Get() << "\n";    


        for (Vertex v : mesh->vertices()) {
            // VP[v] += 0.5 * edgeLmin * (randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][0].normalize() +
            // randomReal(-0.5, 0.5) * geometry->vertexTangentBasis[v][1].normalize());
            EG->vertexPositions[v].y *= 1.0001 + 1.0* pressure.Get() / Youngs.Get();
            EG->vertexPositions[v].x *= 1.0001 + 1.0 * pressure.Get() / Youngs.Get();
            EG->vertexPositions[v].z *= 1.0001 + 1.0 *  pressure.Get() / Youngs.Get();            
            EG->vertexPositions[v].x += 0;
        }
        EG->requireVertexNormals();
        EG->requireFaceVolume();
        VertexData<bool> moved = VertexData<bool>(*mesh, false);
        for (Face f : mesh->faces()) {           
            if (EG->faceVolume[f] > 0) {
                for (Vertex v: f.adjacentVertices())
                    if (!moved[v]) {
                        EG->vertexPositions[v] -= 0.00 * EG->vertexNormals[v];
                    }
            }
        }
        EG->refreshQuantities();
        EG->computeGradient();



        // SET REFERENCE CURVATURES VIA PARAMETERS
        /*EG->requireReferenceCurvature();
        double curvfact = -1.;
        if (anotherVal.Get() != anotherVal.GetDefault()) curvfact = anotherVal.Get();
        for (Edge e : mesh->edges())
        {
           EG->referenceEdgeDihedralAngles[e] = curvfact * EG->referenceEdgeDihedralAngles[e];
        }

        EG->refreshQuantities();
        EG->computeGradient();*/

        double meangrad = EG->elasticGradient[10].norm();

        std::cout << "\n  Gradient Mean: " << meangrad;
    }

    
    ///delete this , just test
    //EG->requireReferenceCurvature();
    //double curvfact =-1.0;
    //if (anotherVal.Get() != anotherVal.GetDefault()) curvfact = anotherVal.Get();
    //for (Edge e : mesh->edges()) {
    //    //std::cout << EG->referenceEdgeDihedralAngles[e] << ",";
    //    EG->referenceEdgeDihedralAngles[e] = curvfact * EG->referenceEdgeDihedralAngles[e];
    //    //std::cout << EG->referenceEdgeDihedralAngles[e] << "!\n";
    //}
    //EG->refreshQuantities();
    //EG->computeGradient();


    
   std::cout << "\n" << "Finished angles" << "\n";   

    // Initialize polyscope
    polyscope::init();
    polyscope::options::alwaysRedraw = true;

    // Set the callback function
    //polyscope::state::userCallback = functionCallback;

    // Add mesh to GUI
    psMesh = polyscope::registerSurfaceMesh(polyscope::guessNiceNameFromPath(filepath), EG->vertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    psMesh->addFaceScalarQuantity("initial energy density elastic", EG->elasticEnergy/EG->faceAreas);
    psMesh->addFaceScalarQuantity("thickness", EG->thickness);
    psMesh->addFaceScalarQuantity("initial energy density elastic stretching", EG->stretchingEnergy);
    psMesh->addFaceScalarQuantity("initial energy density bending", EG->bendingEnergy);
    psMesh->addEdgeScalarQuantity("lengths reference", EG->referenceLengths);
//    psMesh->addEdgeScalarQuantity("reference angles0", EG->referenceEdgeDihedralAngles);
    EG->requireVertexDualAreas();
    EG->requireVertexNormals();
    EG->requireFaceVolume();
    EG->requireTotalEnergy();
    EG->computeGradient();
    psMesh->addVertexVectorQuantity("initial gradient", 1 * EG->elasticGradient);
    psMesh->addVertexScalarQuantity("area vertex", EG->vertexDualAreas);
    //psMesh->addVertexVectorQuantity("-Grad_Norm",-1 * EG->elasticGradient / EG->vertexDualAreas);
    psMesh->addVertexVectorQuantity("normal vertexes", EG->vertexNormals * EG->pressure);
  //  psMesh->addVertexVectorQuantity("Forces", EG->vertexNormals * EG->vertexDualAreas * EG->pressure);
    psMesh->addFaceScalarQuantity("intial energy total integrated", EG->totalEnergy);
    psMesh->addFaceScalarQuantity("intial area face", EG->faceAreas);
    EG->requireEdgeDihedralAngles();
    psMesh->addEdgeScalarQuantity("initial dihedral angles actual", EG->edgeDihedralAngles);
    psMesh->addFaceScalarQuantity("inital volume faces", EG->faceVolume);
    psMesh->addFaceVectorQuantity("faceNorml", EG->faceNormals);
    psMesh->addVertexVectorQuantity("positionNorml", EG->vertexPositions);


    psMesh0 = polyscope::registerSurfaceMesh("OLD", geometry->vertexPositions,
                                            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    psMesh0->addFaceScalarQuantity("Elastic Energy0", EG->elasticEnergy / EG->faceAreas);
    psMesh0->addFaceScalarQuantity("thickness", EG->thickness);
    psMesh0->addFaceScalarQuantity("stretch Energy0", EG->stretchingEnergy);
    psMesh0->addFaceScalarQuantity("bend Energy0", EG->bendingEnergy);
    psMesh0->addEdgeScalarQuantity("reference lengths0", EG->referenceLengths);
    //    psMesh->addEdgeScalarQuantity("reference angles0", EG->referenceEdgeDihedralAngles);
    EG->requireVertexDualAreas();
    EG->requireVertexNormals();
    EG->requireFaceVolume();
    EG->requireTotalEnergy();
    EG->computeGradient();
    psMesh0->addVertexVectorQuantity("-Grad0", 1 * EG->elasticGradient);
    psMesh0->addVertexScalarQuantity("vertex Area", EG->vertexDualAreas);
    // psMesh->addVertexVectorQuantity("-Grad_Norm",-1 * EG->elasticGradient / EG->vertexDualAreas);
    psMesh0->addVertexVectorQuantity("Vertex Normal", EG->vertexNormals * EG->pressure);
    //  psMesh->addVertexVectorQuantity("Forces", EG->vertexNormals * EG->vertexDualAreas * EG->pressure);
    psMesh0->addFaceScalarQuantity("Total Energy0", EG->totalEnergy);
    psMesh0->addFaceScalarQuantity("face Area", EG->faceAreas);
    EG->requireEdgeDihedralAngles();
    psMesh0->addEdgeScalarQuantity("DihedralAngle", EG->edgeDihedralAngles);
    psMesh0->addFaceScalarQuantity("faceVolume", EG->faceVolume);
    psMesh0->addFaceVectorQuantity("faceNorml", EG->faceNormals);
    psMesh0->addVertexVectorQuantity("positionNorml", EG->vertexPositions);
    


    std::cout << "\n \n Total Volume:" << EG->faceVolume.toVector().sum();
    std::cout << "\n \n Total Free Energy:" << EG->totalEnergy.toVector().sum();
    std::cout << "\n \n Total Elastic Energy:" << EG->elasticEnergy.toVector().sum() << "\n  \n";
    
    VertexData<double> force_dist = VertexData<double>(*mesh, -100);
    for (Vertex v : mesh->vertices()) {
        Vector3 force_diff = EG->elasticGradient[v] + EG->vertexNormals[v] * EG->vertexDualAreas[v] * EG->pressure;
        force_dist[v] = force_diff.norm() / EG->vertexDualAreas[v] / EG->pressure;    
    }
    psMesh->addVertexScalarQuantity("forces fit", force_dist);
       // psMesh->addEdgeScalarQuantity("dihedral angles difference0",
                                                              //EG->edgeDihedralAngles - EG->referenceEdgeDihedralAngles);
    
    


        //polyscope::state::userCallback = myCallback;  

    for (Face f : mesh->faces()) {
        if (EG->faceArea(f) < 0) {
            std::cout << "Negative face area detected (inline)";
        }
        if (EG->faceAreas[f] < 0) {
            std::cout << "Negative face area detected (data)";
        }
    }
    
  //polyscope::show();

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

    
    std::cout << "\n pressure, thickness = " << pressure.Get() << "," << thickness.Get() << "  " << EG->pressure << ","
              << EG->thickness[12] << "\n";

    //polyscope::show();
    EG->requireEdgeDihedralAngles();
    EG->requireActualMetric();    
    EG->requireElasticEnergy();
    EG->requireBendingEnergy();
    EG->requireStretchingEnergy();
    std::cout << "\n ********************************************************************************** \n"
              << " ****************  PRINTING STUFF FOR DEBUG **************************************** \n"
              << " *********************************************************************************** \n \n";
    std::cout << "Geometry scale: " << EG->coordinate_scale << "\n";
    std::cout << "Mean edge length: " << EG->referenceLengths.toVector().mean() << "\n";
    std::cout << "Min edge length: " << EG->referenceLengths.toVector().minCoeff() << "\n";
    std::cout << "Max edge length: " << EG->referenceLengths.toVector().maxCoeff() << "\n";
    std::cout << "Reference metric: " << EG->referenceMetric[0][1] << "\n";
    Vector3 avg_metric = {0, 0, 0};
    Vector3 avg_curv = {0, 0, 0};
    Eigen::Matrix3f avg_elastictensor = Eigen::Matrix3f().setZero();
    for (Face f : mesh->faces()) {
        avg_metric[0] += EG->referenceMetric[f][0] / mesh->nFaces();
        avg_metric[1] += EG->referenceMetric[f][1] / mesh->nFaces();
        avg_metric[2] += EG->referenceMetric[f][2] / mesh->nFaces();

        avg_curv[0] += EG->referenceCurvature[f][0] / mesh->nFaces();
        avg_curv[1] += EG->referenceCurvature[f][1] / mesh->nFaces();
        avg_curv[2] += EG->referenceCurvature[f][2] / mesh->nFaces();

        avg_elastictensor(0, 0) += (0+1*EG->elasticCauchyTensor[f](0, 0) )/ mesh->nFaces();
        avg_elastictensor(0, 1) += (0+1*EG->elasticCauchyTensor[f](0, 1) )/ mesh->nFaces();
        avg_elastictensor(0, 2) += (0+1*EG->elasticCauchyTensor[f](0, 2) )/ mesh->nFaces();
        avg_elastictensor(1, 0) += (0+1*EG->elasticCauchyTensor[f](1, 0) )/ mesh->nFaces();
        avg_elastictensor(1, 1) += (0+1*EG->elasticCauchyTensor[f](1, 1) )/ mesh->nFaces();
        avg_elastictensor(1, 2) += (0+1*EG->elasticCauchyTensor[f](1, 2) )/ mesh->nFaces();
        avg_elastictensor(2, 0) += (0+1*EG->elasticCauchyTensor[f](2, 0) )/ mesh->nFaces();
        avg_elastictensor(2, 1) += (0+1*EG->elasticCauchyTensor[f](2, 1) )/ mesh->nFaces();
        avg_elastictensor(2, 2) += (0+1*EG->elasticCauchyTensor[f](2, 2) )/ mesh->nFaces();
    }

    std::cout << "mean_ref_metric: " << "{" << avg_metric[0] << "," << avg_metric[1] << "," << avg_metric[2] << "} \n";
    std::cout << "mean_ref_curv: " << "{" << avg_curv[0] << "," << avg_curv[1] << "," << avg_curv[2] << "} \n";
    std::cout << "mean_elastic_tensor (line 1): " << "{" << avg_elastictensor(0, 0) << "," << avg_elastictensor(0, 1) << "," << avg_elastictensor(0, 2) << "} \n";
    std::cout << "mean_elastic_tensor (line 2): " << "{" << avg_elastictensor(1, 0) << "," << avg_elastictensor(1, 1) << "," << avg_elastictensor(1, 2) << "} \n";
    std::cout << "mean_elastic_tensor (line 3): " << "{" << avg_elastictensor(2, 0) << "," << avg_elastictensor(2, 1) << "," << avg_elastictensor(2, 2) << "} \n";

    EG->requireFaceGaussianCurvatures();
    EG->requireVertexGaussianCurvatures();
    EG->requireVertexMeanCurvatures();
    EG->requireVertexDualAreas();

    if (false) {
        for (Edge e: mesh->edges()) {
            std::cout << "\n edge index: " << e.getIndex() << "\n";
            std::cout << "edge length: " << EG->edgeLengths[e] << "\n";
            Halfedge he = e.halfedge();
            std::cout << "halfe edge index: " << he.getIndex() << "\n";
            std::cout << "twin index: " << he.twin().getIndex() << "\n";
            std::cout << "twin edge index: " << he.twin().edge().getIndex() << "\n";
            std::cout << "twin edge length: " << EG->edgeLength(he.twin().edge())<< "\n";

            
        }
    }

   /* for (Edge e : mesh->edges()) {
        EG->referenceLengths[e] = e.getIndex();
    }*/

    if (false) {

        for (Face f : mesh->faces()) {
            std::cout << "\n \n metrics: {(" << EG->actualMetric[f][0] << "," << EG->referenceMetric[f][0] << "),("
                      << EG->actualMetric[f][1] << "," << EG->referenceMetric[f][1] << "),(" << EG->actualMetric[f][2]
                      << "," << EG->referenceMetric[f][2] << ")} \n\n";
            std::cout << "edge length: " << EG->edgeLengths[f.halfedge().edge()] << " scale: " << EG->coordinate_scale;
            std::cout << "\n curvaures: {(" << EG->actualCurvature[f][0] << "," << EG->referenceCurvature[f][0] << "),("
                      << EG->actualCurvature[f][1] << "," << EG->referenceCurvature[f][1] << "),("
                      << EG->actualCurvature[f][2] << "," << EG->referenceCurvature[f][2] << ")} \n";
            std::cout << "semiangle: " << 0.5 * EG->edgeDihedralAngles[f.halfedge().edge()] << " scale: " << EG->coordinate_scale;
            double met_det =
                EG->actualMetric[f][0] * EG->actualMetric[f][1] - EG->actualMetric[f][2] * EG->actualMetric[f][2];
            Vector3 invmet = {EG->actualMetric[f][1] / met_det, EG->actualMetric[f][0] / met_det,
                              -EG->actualMetric[f][2] / met_det};
            double shape_operator_mean =
                0.5 * (invmet[0] * EG->actualCurvature[f][0] + 2 * EG->actualCurvature[f][2] * invmet[2] +
                       EG->actualCurvature[f][1] * invmet[1]);
            double shape_operator_det = 1 / met_det *
                                        (EG->actualCurvature[f][0] * EG->actualCurvature[f][1] -
                                         EG->actualCurvature[f][2] * EG->actualCurvature[f][2]);
            std::cout << "\n  Shape: { " << shape_operator_mean / EG->coordinate_scale << ","
                      << shape_operator_det / EG->coordinate_scale / EG->coordinate_scale << "} \n";
            double mean_mean = 0;
            double mean_gauss = 0;
            for (Vertex v : f.adjacentVertices()) {
                mean_mean += EG->vertexMeanCurvatures[v] / EG->vertexDualAreas[v] /3;
                mean_gauss += EG->vertexGaussianCurvatures[v] / EG->vertexDualAreas[v]/ 3;
            }
            std::cout << " geometry: {" << mean_mean << "," << mean_gauss << ","
                      << EG->faceGaussianCurvatures[f] / EG->faceAreas[f] << "} \n";
            // std::cout << "elasrtic_tensor (line 1): " << "{" << EG->elasticCauchyTensor[f](0, 0) << "," <<
            // EG->elasticCauchyTensor[f](0, 1) << "," << EG->elasticCauchyTensor[f](0, 2) << "} \n"; std::cout <<
            // "elasrtic_tensor (line 2): " << "{" << EG->elasticCauchyTensor[f](1, 0) << "," <<
            // EG->elasticCauchyTensor[f](1, 1) << "," << EG->elasticCauchyTensor[f](1, 2) << "} \n"; std::cout <<
            // "elasrtic_tensor (line 3): " << "{" << EG->elasticCauchyTensor[f](2, 0) << "," <<
            // EG->elasticCauchyTensor[f](2, 1) << "," << EG->elasticCauchyTensor[f](2, 2) << "} \n";
            std::cout << "\n stetching content: " << EG->stretchingEnergy[f] << "\n";
            std::cout << "\n bending content: " << EG->bendingEnergy[f] << "\n";
        }
    }

    //ShowPolyscope(0);
    //polyscope::show();
    mySubroutine();
    mySubroutine3();


    return EXIT_SUCCESS;
}