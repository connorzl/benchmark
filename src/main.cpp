#include "polyscope/polyscope.h"
//#include "polyscope/overlap.h"

#include <iostream>
#include <fstream> 
#include <iomanip>

#include "geometrycentral/geometry.h"
#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/polygon_soup_mesh.h"
#include "geometrycentral/direction_fields.h"
#include "geometrycentral/distortion.h"

// for testing purposes
#include "geometrycentral/tutte.h"
#include "geometrycentral/spectral_conformal.h"
#include "geometrycentral/least_squares_conformal.h"
#include "geometrycentral/arap.h"
#include "geometrycentral/quad_cover.h"
#include "geometrycentral/stripes.h"
#include "geometrycentral/harmonicbases.h"
#include "geometrycentral/BHalfedgemesh.h"
#include "polyscope/quad_mesh.h"

#include "imgui.h"
#include "args/args.hxx"

#include "polyscope/scatterplot.h"
#include "polyscope/gl/colormap_sets.h"

#include <map> 

using namespace geometrycentral;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

// == Program data
// (in general, such data should probably be stored in a class, or whatever makes sense for your situation -- these globals are just for the sake of a simple example app)
Geometry<Euclidean>* geom;
HalfedgeMesh* mesh;
polyscope::Scatterplot* scatter;

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  // Begin an ImGUI window
  static bool showGui = true;

  ImGui::Begin("Sample Scatterplot", &showGui); //ImGuiWindowFlags_AlwaysAutoResize);
  ImGui::PushItemWidth(100);
  (*scatter).buildUI();

  // Cleanup the ImGUI window
  ImGui::PopItemWidth();
  ImGui::End();
}

void writeToFile(std::ofstream &outfile, std::string path, Vector3 areaDistortion, Vector3 angleDistortion, 
                 size_t trianglesFlipped, bool globalOverlap, size_t seamLength, double numFaces) {
  std::size_t start = path.find_last_of('/');
  std::size_t end = path.find_last_of('.');
  std::string fileName = path.substr(start+1, end-start-1);

  std::ostringstream area;
  area << areaDistortion[0] << " " << areaDistortion[1] << " " << areaDistortion[2];

  std::ostringstream angle;
  angle << angleDistortion[0] << " " << angleDistortion[1] << " " << angleDistortion[2];
  
  outfile << std::setw(20) << std::left << fileName;
  outfile << std::setw(50) << std::left << path;
  outfile << std::setw(35) << std::left << area.str();
  outfile << std::setw(35) << std::left << angle.str();
  outfile << std::setw(25) << std::left << trianglesFlipped;
  outfile << std::setw(25) << std::left << globalOverlap;
  outfile << std::setw(25) << std::left << seamLength;
  outfile << std::setw(25) << std::left << numFaces << std::endl;
}

void generateAnalysisFile(std::vector<std::string> objFiles) {
    std::ofstream outfile ("analysis.txt");
    outfile << std::setw(20) << std::left << "File" 
            << std::setw(50) << std::left << "Path"
              << std::setw(35) << std::left << "Area Distortion (Min,Max,Avg)" 
              << std::setw(35) << std::left << "Angle Distortion (Min,Max,Avg)" 
              << std::setw(25) << std::left << "Triangles Flipped"
              << std::setw(25) << std::left << "Global Overlap"
              << std::setw(25) << std::left << "Seam Length"
              << std::setw(25) << std::left << "Number Faces"
              << std::endl;

    for (size_t i = 0; i < objFiles.size(); i++) {
      mesh = new HalfedgeMesh(PolygonSoupMesh(objFiles[i]), geom);
      if (geom->paramCoords.size() == 0) {
        std::cout<<"gg"<<std::endl;
        continue;
      }
      // Compute distortion metrics
      Distortion* d = new Distortion(mesh, geom);
      Vector3 areaDistortion = d->computeAreaScaling();
      Vector3 angleDistortion = d->computeQuasiConformalError();
      size_t trianglesFlipped = d->computeTriangleFlips();
      bool globalOverlap = d->computeGlobalOverlap();
      size_t seamLength = d->computeSeamLength();
      double numFaces = mesh->nFaces();
      writeToFile(outfile, objFiles[i], areaDistortion, angleDistortion, trianglesFlipped, globalOverlap, seamLength, numFaces);
    }
    outfile.close();
}

void parseAnalysisFile(std::string filename, std::vector<std::string> &meshnames, std::vector<std::string> &filepaths,
  std::vector<double> &minAreaDistortion, std::vector<double> &maxAreaDistortion, std::vector<double> &avgAreaDistortion,
  std::vector<double> &minAngleDistortion, std::vector<double> &maxAngleDistortion, std::vector<double> &avgAngleDistortion,
  std::vector<double> &trianglesFlipped, std::vector<double> &globalOverlap, std::vector<double> &seamLengths, 
  std::vector<double> &numFaces, const std::map<std::string, bool>& sharedMeshes = std::map<std::string, bool>(),
  const std::string suffix = "") {

  std::ifstream in(filename);
  if (!in) throw std::invalid_argument("Could not open analysis file " + filename);

  std::string line;
  getline(in,line); // skip first line, which only contains labels for readability
  while (getline(in, line)) {
    std::stringstream ss(line);
    std::string token;
    double minVal, maxVal, avgVal;
    double d;

    // parse name
    ss >> token;
    if (sharedMeshes.size() > 0) {
      if (sharedMeshes.find(token) == sharedMeshes.end()) {
        continue;
      } else {
        meshnames.push_back(token + "-" + suffix);
      }
    } else {
      meshnames.push_back(token);
    }

    // parse filepath
    ss >> token;
    filepaths.push_back(token);

    // parse area distortion
    ss >> minVal >> maxVal >> avgVal;
    minAreaDistortion.push_back(minVal);
    maxAreaDistortion.push_back(maxVal);
    avgAreaDistortion.push_back(avgVal);

    // parse angle distortion
    ss >> minVal >> maxVal >> avgVal;
    minAngleDistortion.push_back(minVal);
    maxAngleDistortion.push_back(maxVal);
    avgAngleDistortion.push_back(avgVal);

    // parse triangles flipped
    ss >> d;
    trianglesFlipped.push_back(d);
    
    // parse global overlap
    ss >> d;
    globalOverlap.push_back(d);

    // parse total seam length
    ss >> d;
    seamLengths.push_back(d);

    // parse number faces
    ss >> d;
    numFaces.push_back(d);
  }
}

std::map<std::string,bool> findSharedMeshes(std::vector<std::string> analysisFiles) {
  std::map<std::string,bool> sharedObjFiles;
  std::map<std::string,size_t> meshCounts;

  for (size_t i = 0; i < analysisFiles.size(); i++) {
    std::ifstream in(analysisFiles[i]);
    if (!in) throw std::invalid_argument("Could not open analysis file " + analysisFiles[i]);

    std::string line;
    getline(in,line); // skip first line, which only contains labels for readability

    while (getline(in, line)) {
      std::stringstream ss(line);
      std::string fileName;
      std::string filePath;

      // parse fileName and filePath
      ss >> fileName;
      if ( meshCounts.find(fileName) == meshCounts.end() ) {
        meshCounts[fileName] = 1;
      } else {
        meshCounts[fileName] += 1;
      }
    }
  }

  for(std::map<std::string, size_t>::const_iterator it = meshCounts.begin();
    it != meshCounts.end(); ++it) {
    // only take the meshes that all the analysis files share
    if (it->second == analysisFiles.size()) {
        sharedObjFiles[it->first] = true;
    }
  }
  return sharedObjFiles;
}

void parseAnalysisMeshes(std::map<std::string,bool> objFiles, std::vector<std::string> analysisFiles) {
  // information that we want to parse for each mesh
  std::vector<std::string> meshnames;
  std::vector<std::string> filepaths;
  std::vector<double> minAreaDistortion;
  std::vector<double> maxAreaDistortion;
  std::vector<double> avgAreaDistortion;
  std::vector<double> minAngleDistortion;
  std::vector<double> maxAngleDistortion;
  std::vector<double> avgAngleDistortion;
  std::vector<double> trianglesFlipped;
  std::vector<double> globalOverlap;
  std::vector<double> seamLengths;
  std::vector<double> numFaces;

 for (size_t i = 0; i < analysisFiles.size(); i++) {
    parseAnalysisFile(analysisFiles[i], meshnames, filepaths, minAreaDistortion, maxAreaDistortion, avgAreaDistortion,
    minAngleDistortion, maxAngleDistortion, avgAngleDistortion, trianglesFlipped, globalOverlap, seamLengths,
    numFaces, objFiles, std::to_string(i));
 }

  std::vector<std::vector<double>> data = {minAreaDistortion, maxAreaDistortion, avgAreaDistortion,
              minAngleDistortion,maxAngleDistortion,avgAngleDistortion};
  std::vector<char*> labels = {(char*)"min area distortion",(char*)"max area distortion",(char*)"avg area distortion",
  (char*)"min angle distortion",(char*)"max angle distortion",(char*)"avg angle distortion"};
  std::vector<std::vector<double>> additionalData = {trianglesFlipped, globalOverlap, seamLengths};
  std::vector<std::string> additionalDataLabels = {"Triangles Flipped", "Global Overlap", "Total Seam Length"};

  scatter = new polyscope::Scatterplot();
  (*scatter).buildScatterplot(data, labels, meshnames, filepaths, additionalData, additionalDataLabels, numFaces, analysisFiles.size());
  (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[0]);
  polyscope::show();
}

int main(int argc, char** argv) {
  if (strcmp(argv[1],"-a") == 0) {
    std::cout << "Generating analysis file" << std::endl;
    std::vector<std::string> objFiles;
    for (int i = 2; i < argc; i++) {
      objFiles.push_back(argv[i]);
    }
    generateAnalysisFile(objFiles);

  } else if (strcmp(argv[1],"-d") == 0) {
    std::cout << "Visualizing" << std::endl;
    polyscope::init();
    polyscope::state::userCallback = myCallback;

    // information that we want to parse for each mesh
    std::vector<std::string> meshnames;
    std::vector<std::string> filepaths;
    std::vector<double> minAreaDistortion;
    std::vector<double> maxAreaDistortion;
    std::vector<double> avgAreaDistortion;
    std::vector<double> minAngleDistortion;
    std::vector<double> maxAngleDistortion;
    std::vector<double> avgAngleDistortion;
    std::vector<double> trianglesFlipped;
    std::vector<double> globalOverlap;
    std::vector<double> seamLengths;
    std::vector<double> numFaces;

    parseAnalysisFile(argv[2], meshnames, filepaths, minAreaDistortion, maxAreaDistortion, avgAreaDistortion, 
      minAngleDistortion, maxAngleDistortion, avgAngleDistortion, trianglesFlipped, globalOverlap, seamLengths, numFaces);
    
    std::vector<std::vector<double>> data = {minAreaDistortion, maxAreaDistortion, avgAreaDistortion,
              minAngleDistortion,maxAngleDistortion,avgAngleDistortion};
    std::vector<char*> labels = {(char*)"min area distortion",(char*)"max area distortion",(char*)"avg area distortion",
    (char*)"min angle distortion",(char*)"max angle distortion",(char*)"avg angle distortion"};
    std::vector<std::vector<double>> additionalData = {trianglesFlipped, globalOverlap, seamLengths};
    std::vector<std::string> additionalDataLabels = {"Triangles Flipped", "Global Overlap", "Total Seam Length"};

    scatter = new polyscope::Scatterplot();
    (*scatter).buildScatterplot(data, labels, meshnames, filepaths, additionalData, additionalDataLabels, numFaces);
    (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[0]);
    polyscope::show();

  } else if (strcmp(argv[1],"-c") == 0) {
     std::cout << "Comparing" << std::endl;
     polyscope::init();
     polyscope::state::userCallback = myCallback;
     std::vector<std::string> analysisFiles;
     for (int i = 2; i < argc; i++) {
       analysisFiles.push_back(argv[i]);
     }

     // find the meshes that both files share using a map
     std::map<std::string,bool> sharedObjFiles = findSharedMeshes(analysisFiles);
     // construct a scatterplot using only the shared meshes
     parseAnalysisMeshes(sharedObjFiles, analysisFiles);

  } else {
    // Configure the argument parser
    args::ArgumentParser parser("Polyscope sample program. See github.com/nmwsharp/polyscope/examples.");
    args::Positional<string> inFileName(parser, "input_file", "An .obj file to visualize");

    // Parse args
    try {
      parser.ParseCLI(argc, argv);
    } catch (args::Help) {
      std::cout << parser;
      return 0;
    } catch (args::ParseError e) {
      std::cerr << e.what() << std::endl;
      std::cerr << parser;
      return 1;
    }

    // Make sure a mesh name was given
    if(args::get(inFileName) == "") {
      std::cerr << "Please specify .obj file as argument" << std::endl;
      return EXIT_FAILURE;
    }
  
    // == Build the mesh object from the input file
    mesh = new HalfedgeMesh(PolygonSoupMesh(args::get(inFileName)), geom);
    //Tutte t = Tutte(mesh,geom);
    //t.computeTutteEmbedding();
    //SpectralConformal s = SpectralConformal(mesh,geom);
    //s.computeSpectralConformal();
    //LSCM l = LSCM(mesh,geom);
    //l.computeLSCM();
    //ARAP a = ARAP(mesh,geom);
    //a.computeARAP();
    
    //QuadCover Q = QuadCover(mesh,geom);
    //VertexData<std::complex<double>> field = Q.computeCrossField();
    //FaceData<int> singularities = Q.computeSingularities();
    //Q.computeBranchCover();
    //VertexData<double> offsets = Q.computeOffset();
   
    //Stripes S = Stripes(mesh,geom);
    //VertexData<std::complex<double>> field = S.computeField();
    //FaceData<int> singularities = S.computeSingularities();
    //S.edgeData();
    //S.computeStripes();

    //VertexData<double> offsets = Q.computeOffset(); this is for quad_cover double cover
    polyscope::init();
    std::string meshNiceName = polyscope::utilities::guessNiceNameFromPath(args::get(inFileName));
    polyscope::registerSurfaceMesh(meshNiceName, geom);
    //polyscope::getSurfaceMesh()->enabled = true;

    QuadMesh M = QuadMesh(mesh,geom);
    M.computeCrossField();
    M.computeSingularities();
    M.uniformize();
    M.computeBranchCover();
    M.computeCrossFieldCMBranchCover();
    M.computeStripes();

    //M.optimizeHarmonic();
    //M.textureCoordinates();
    //M.visualize();
    //polyscope::show();
    
    //M.optimizeSimpleLocally();
    //M.optimizeSimpleGlobally();

    int iter = 0;
    while (!M.textureCoordinates()) {
      std::cout << "Optimization: " << iter << std::endl;
      M.optimizeHarmonic();  
      iter++;
    }
    std::cout << "TOTAL ITERS: " << iter << std::endl;
    M.visualize();
    polyscope::show();
  
    /*
    //std::ofstream outfile ("eigenvalues.txt");
    double scale = 100;//16 * PI;
    std::complex<double> init(1,0);
    // 32, 122, 212, 302 for bunny
    for (int t = 0; t < 90; t+=1) {
      double rad = t * (PI / 180.0);
      std::complex<double> rot = std::exp(IM_I * rad);
      M.computeCrossFieldCMBranchCover(rot * init, scale);
      double lambda = M.computeStripes();
      std::cout << "rotation and lambda: " << t << "," << lambda << std::endl;
      //outfile << lambda << std::endl;
      M.textureCoordinates();
      M.visualize();
      polyscope::screenshot();
      polyscope::show();
      
      //break;
    }
    //outfile.close();
    polyscope::show();
    */
    /*
    HarmonicBases HB = HarmonicBases(mesh,geom);
    std::vector<Eigen::MatrixXd> bases = HB.compute();
    
    std::vector<std::vector<HalfedgePtr>> generators = HB.generators;
    for (size_t i = 0; i < generators.size(); i++) {
      FaceData<double> generatorFaces(mesh);
      for (FacePtr f : mesh->faces()) {
        generatorFaces[f] = 0;
      }
      for (size_t j = 0; j < generators[i].size(); j++) {
        HalfedgePtr he = generators[i][j];
        generatorFaces[he.face()] = 1;
      }
      std::string name = "Generator " + std::to_string(i);
      polyscope::getSurfaceMesh()->addQuantity(name,generatorFaces);
    }

    EdgeData<size_t> edgeIndices = mesh->getEdgeIndices();
    EdgeData<double> data(mesh);
    EdgeData<double> data2(mesh);
    for (EdgePtr e : mesh->edges()) {
      size_t index = edgeIndices[e];
      data[e] = HB.bases[0](index,0);
      data2[e] = HB.bases[1](index,0);
    }
    polyscope::getSurfaceMesh()->addVectorQuantity("1-form on edges 0", data);
    polyscope::getSurfaceMesh()->addVectorQuantity("1-form on edges 1", data2);

    FaceData<std::complex<double>> X = HB.visualize();
    polyscope::getSurfaceMesh()->addVectorQuantity("Harmonic Field", X);
    polyscope::show();
    */
    // == Add Distortion Data to Mesh, if any
    if (geom->paramCoords.size() > 0) {
      Distortion* d = new Distortion(mesh, geom);
      
      // Area Distortion
      Vector3 areaDistortion = d->computeAreaScaling();
      std::cout << "AREA DISTORTION: min: " << areaDistortion[0] << " max: " << areaDistortion[1] << " avg: " << areaDistortion[2] << std::endl;
      polyscope::getSurfaceMesh()->addQuantity("Area Distortion", d->areaDistortion);

      // Angle Distortion
      Vector3 angleDistortion = d->computeQuasiConformalError();
      std::cout << "ANGLE DISTORTION: min: " << angleDistortion[0] << " max: " << angleDistortion[1] << " avg: " << angleDistortion[2] << std::endl;
      polyscope::getSurfaceMesh()->addQuantity("Angle Distortion", d->angleDistortion);

      // Triangles flipped
      size_t trianglesFlipped = d->computeTriangleFlips();
      std::cout << "TRIANGLES FLIPPED: " << trianglesFlipped << std::endl;
      polyscope::getSurfaceMesh()->addQuantity("Flipped Triangles", d->trianglesFlipped);
      
      // Global Overlap
      bool globalOverlap = d->computeGlobalOverlap();
      std::cout << "GLOBAL OVERLAP: " << globalOverlap << std::endl;

      // Total Seam Length
      size_t seamLength = d->computeSeamLength();
      std::cout << "TOTAL SEAM LENGTH: " << seamLength << std::endl;
    }
    
    // branch cover visualization
    //for (VertexPtr v : mesh->vertices()) {
    //  offsets[v] = std::abs(offsets[v]);
    //}
    //polyscope::getSurfaceMesh()->addDistanceQuantity("Offsets", offsets);
    
    // Register the user callback 
    /*
    polyscope::state::userCallback = myCallback;
    scatter = new polyscope::Scatterplot();
    std::vector<double> xs = {1, 2, 3, 4, 5, 6};
    std::vector<double> ys = {2, 4, 8, 16, 32, 64};
    std::vector<double> zs = {1, 2, 3, 4, 5, 6};
    (*scatter).buildScatterplot(xs, ys, zs);
    (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[0]);
    std::vector<std::vector<double>> data = {xs, ys, zs};
    std::vector<char*> labels = {(char*)"A", (char*)"B", (char*)"C"};
    
    (*scatter).buildScatterplot(data,labels);
    (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[0]);
    */
  }
  return EXIT_SUCCESS;
}
