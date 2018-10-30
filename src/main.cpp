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

#include "imgui.h"
#include "args/args.hxx"

#include "polyscope/scatterplot.h"
#include "polyscope/gl/colormap_sets.h"

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
                 size_t trianglesFlipped, bool globalOverlap) {
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
  outfile << std::setw(25) << std::left << globalOverlap << std::endl;
}

int main(int argc, char** argv) {
  if (strcmp(argv[1],"-a") == 0) {
    std::cout << "Generating analysis file" << std::endl;
    std::vector<std::string> objFiles;
    for (int i = 2; i < argc; i++) {
      objFiles.push_back(argv[i]);
    }

    std::ofstream outfile ("analysis.txt");
    outfile << std::setw(20) << std::left << "File" 
            << std::setw(50) << std::left << "Path"
              << std::setw(35) << std::left << "Area Distortion (Min,Max,Avg)" 
              << std::setw(35) << std::left << "Angle Distortion (Min,Max,Avg)" 
              << std::setw(25) << std::left << "Triangles Flipped"
              << std::setw(25) << std::left << "Global Overlap"
              << std::endl;

    for (size_t i = 0; i < objFiles.size(); i++) {
      mesh = new HalfedgeMesh(PolygonSoupMesh(objFiles[i]), geom);
      if (geom->paramCoords.size() == 0) continue;
    
      // Compute distortion metrics
      Distortion* d = new Distortion(mesh, geom);
      Vector3 areaDistortion = d->computeAreaScaling();
      Vector3 angleDistortion = d->computeQuasiConformalError();
      size_t trianglesFlipped = d->computeTriangleFlips();
      bool globalOverlap = d->computeGlobalOverlap();

      // Total Seam Length
      size_t seamLength = d->computeSeamLength();

      //std::cout << "AREA DISTORTION: min: " << areaDistortion[0] << " max: " << areaDistortion[1] << " avg: " << areaDistortion[2] << std::endl;
      //std::cout << "ANGLE DISTORTION: min: " << angleDistortion[0] << " max: " << angleDistortion[1] << " avg: " << angleDistortion[2] << std::endl;
      //std::cout << "TRIANGLES FLIPPED: " << trianglesFlipped << std::endl;
      //std::cout << "GLOBAL OVERLAP: " << globalOverlap << std::endl;
      //std::cout << "TOTAL SEAM LENGTH: " << seamLength << std::endl;
      writeToFile(outfile, objFiles[i], areaDistortion, angleDistortion, trianglesFlipped, globalOverlap);
    }
    outfile.close();
  } else if (strcmp(argv[1],"-d") == 0) {
    std::cout << "Visualizing" << std::endl;

    std::string filename = argv[2];
    std::ifstream in(filename);
    if (!in) throw std::invalid_argument("Could not open mesh file " + filename);
    
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

    std::string line;
    getline(in,line); // skip first line, which only contains labels for readability
    while (getline(in, line)) {
      std::stringstream ss(line);
      std::string token;
      double minVal, maxVal, avgVal;
      double d;

      // parse name
      ss >> token;
      meshnames.push_back(token);

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
    }

    polyscope::init();

     // Register the user callback 
    polyscope::state::userCallback = myCallback;
    
    scatter = new polyscope::Scatterplot();
    (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[0]);
    std::vector<std::vector<double>> data;
    data.push_back(minAreaDistortion);
    data.push_back(maxAreaDistortion);
    data.push_back(avgAreaDistortion);
    data.push_back(minAngleDistortion);
    data.push_back(maxAngleDistortion);
    data.push_back(avgAngleDistortion);

    std::vector<char*> labels;
    labels.push_back((char*)"min area distortion");
    labels.push_back((char*)"max area distortion");
    labels.push_back((char*)"avg area distortion");
    labels.push_back((char*)"min angle distortion");
    labels.push_back((char*)"max angle distortion");
    labels.push_back((char*)"avg angle distortion");

    (*scatter).buildScatterPlot(data,labels);
    
    polyscope::show();

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

    // Initialize polyscope
    polyscope::init();
  
    // == Build the mesh object from the input file
    mesh = new HalfedgeMesh(PolygonSoupMesh(args::get(inFileName)), geom);
    std::string meshNiceName = polyscope::utilities::guessNiceNameFromPath(args::get(inFileName));
    polyscope::registerSurfaceMesh(meshNiceName, geom);

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

    // Register the user callback 
    polyscope::state::userCallback = myCallback;
    
    scatter = new polyscope::Scatterplot();
    (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[0]);
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    for (int i = 0; i < 100; i++) {
      xs.push_back(3 * unitRand() - .5);
      ys.push_back(3 * unitRand() - .5);
      zs.push_back(3 * unitRand() - .5);
    }
    //(*scatter).buildScatterplot(xs, ys, zs);
    
    std::vector<std::vector<double>> data;
    data.push_back(xs);
    data.push_back(ys);
    data.push_back(zs);

    std::vector<char*> labels;
    labels.push_back((char*)"A");
    labels.push_back((char*)"B");
    labels.push_back((char*)"C");

    (*scatter).buildScatterPlot(data,labels);

    // Give control to the polyscope gui
    polyscope::show();
  }
  return EXIT_SUCCESS;

}
