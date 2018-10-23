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
int iColorMap = 0;

// Parameters 
size_t iGeneratedPoints = 0;
int nPts = 100;
float rangeLow = -5.0;
float rangeHigh = 5.0;

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {

  // Begin an ImGUI window
  static bool showGui = true;

  ImGui::Begin("Sample Scatterplot", &showGui); //ImGuiWindowFlags_AlwaysAutoResize);
  ImGui::PushItemWidth(100);
  /*

  // Generate a random function
  ImGui::TextUnformatted("Generate random function:");
  ImGui::DragFloatRange2("Data range", &rangeLow, &rangeHigh);
  if (ImGui::Button("Generate")) {
    VertexData<double> randF(mesh);
    for (VertexPtr v : mesh->vertices()) {
      randF[v] = randomReal(rangeLow, rangeHigh);
    }
    polyscope::getSurfaceMesh()->addQuantity("generated_function", randF);
  }
  ImGui::Separator();


  // Add points 
  ImGui::TextUnformatted("Add new points clouds:");
  ImGui::InputInt("# pts", &nPts, 0, 1000000);
  if (ImGui::Button("Add another")) {
    std::vector<Vector3> points;
    for (int i = 0; i < nPts; i++) {
      points.push_back(3 * Vector3{unitRand() - .5, unitRand() - .5, unitRand() - .5});
    }
    polyscope::registerPointCloud("generated_points_"+std::to_string(iGeneratedPoints), points);
    iGeneratedPoints++;
  }
  ImGui::Separator();

  if (ImGui::Button("Batman")) {

    polyscope::warning("Na na na na na na na na na na na na na Batman!");
  }
  ImGui::Separator();
  */

   // Set colormap
  ImGui::SameLine();
  ImGui::PushItemWidth(100);
  (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[iColorMap]);
  int iColormapBefore = iColorMap;
  ImGui::Combo("##colormap", &iColorMap, polyscope::gl::quantitativeColormapNames,
                IM_ARRAYSIZE(polyscope::gl::quantitativeColormapNames));
  ImGui::PopItemWidth();
  if (iColorMap != iColormapBefore) {
    (*scatter).updateColormap(polyscope::gl::quantitativeColormaps[iColorMap]);
  }

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
  area << areaDistortion[0] << "," << areaDistortion[1] << "," << areaDistortion[2];

  std::ostringstream angle;
  angle << angleDistortion[0] << "," << angleDistortion[1] << "," << angleDistortion[2];
  
  outfile << std::setw(20) << std::left << fileName;
  outfile << std::setw(35) << std::left << area.str();
  outfile << std::setw(35) << std::left << angle.str();
  outfile << std::setw(25) << std::left << trianglesFlipped;
  outfile << std::setw(25) << std::left << globalOverlap << std::endl;
}

int main(int argc, char** argv) {
  bool ANALYZE_DIR = false;
  std::vector<std::string> objFiles;
  if (strcmp(argv[1],"-a") == 0) {
    ANALYZE_DIR = true;
    for (int i = 2; i < argc; i++) {
      objFiles.push_back(argv[i]);
    }
  } 

  if (ANALYZE_DIR) {

    /*
    // Initialize polyscope
    polyscope::init();

    std::ofstream outfile ("analysis.txt");
    outfile << std::setw(20) << std::left << "File" 
              << std::setw(35) << std::left << "Area Distortion (Min,Max,Avg)" 
              << std::setw(35) << std::left << "Angle Distortion (Min,Max,Avg)" 
              << std::setw(25) << std::left << "Triangles Flipped"
              << std::setw(25) << std::left << "Global Overlap"
              << std::endl;

    std::vector<Vector3> areaDistortionPoints;
    std::vector<std::string> meshNames;
    for (size_t i = 0; i < objFiles.size(); i++) {
      mesh = new HalfedgeMesh(PolygonSoupMesh(objFiles[i]), geom);
      if (geom->paramCoords.size() == 0) continue;
      
      std::string meshNiceName = polyscope::utilities::guessNiceNameFromPath(objFiles[i]);
      polyscope::registerSurfaceMesh(meshNiceName, geom);
      Distortion* d = new Distortion(mesh, geom);
      meshNames.push_back(meshNiceName);

      // Area Distortion
      Vector3 areaDistortion = d->computeAreaScaling();
      std::cout << "AREA DISTORTION: min: " << areaDistortion[0] << " max: " << areaDistortion[1] << " avg: " << areaDistortion[2] << std::endl;
      polyscope::getSurfaceMesh(meshNiceName)->addQuantity("Area Distortion", d->areaDistortion);
      areaDistortionPoints.push_back(areaDistortion);

      // Angle Distortion
      Vector3 angleDistortion = d->computeQuasiConformalError();
      std::cout << "ANGLE DISTORTION: min: " << angleDistortion[0] << " max: " << angleDistortion[1] << " avg: " << angleDistortion[2] << std::endl;
      polyscope::getSurfaceMesh(meshNiceName)->addQuantity("Angle Distortion", d->angleDistortion);

      // Triangles flipped
      size_t trianglesFlipped = d->computeTriangleFlips();
      std::cout << "TRIANGLES FLIPPED: " << trianglesFlipped << std::endl;
      polyscope::getSurfaceMesh(meshNiceName)->addQuantity("Flipped Triangles", d->trianglesFlipped);
      
      // Global Overlap
      bool globalOverlap = d->computeGlobalOverlap();
      std::cout << "GLOBAL OVERLAP: " << globalOverlap << std::endl;

      // Total Seam Length
      size_t seamLength = d->computeSeamLength();
      std::cout << "TOTAL SEAM LENGTH: " << seamLength << std::endl;
    }
    outfile.close();

    // Register the user callback 
    polyscope::state::userCallback = myCallback;

    // Give control to the polyscope gui
    polyscope::show();
    */
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
/*
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
*/
    }

    // == Add some data to the mesh we just created
    // Note: Since the viewer only currently only has one mesh, we can omit the mesh name field
    //       from these commands -- otherwise a correct name must be specified to getSurfaceMesh()
    {
      /*
      // Two function on vertices (x coord and a random color)
      VertexData<double> valX(mesh);
      VertexData<Vector3> randColor(mesh);
      for (VertexPtr v : mesh->vertices()) {
        valX[v] = geom->position(v).x;
        randColor[v] = Vector3{unitRand(), unitRand(), unitRand()};
      }
      polyscope::getSurfaceMesh()->addQuantity("x coord", valX);
      polyscope::getSurfaceMesh()->addColorQuantity("random color", randColor);

      // Face area
      FaceData<double> fArea(mesh);
      for (FacePtr f : mesh->faces()) {
        fArea[f] = geom->area(f);
      }
      polyscope::getSurfaceMesh()->addQuantity("face area", fArea, polyscope::DataType::MAGNITUDE);

      // Edge cotan weights
      EdgeData<double> cWeight(mesh);
      geom->getEdgeCotanWeights(cWeight);
      polyscope::getSurfaceMesh()->addQuantity("cotan weight", cWeight, polyscope::DataType::SYMMETRIC);
  
      // Vertex normals
      VertexData<Vector3> normals(mesh);
      geom->getVertexNormals(normals);
      polyscope::getSurfaceMesh()->addVectorQuantity("vertex normals", normals);

      // Smoothest 4-symmetric direction field
      if(mesh->nBoundaryLoops() == 0) { // (haven't implemented for boundary yet...)
        FaceData<Complex> smoothestField = computeSmoothestFaceDirectionField(geom, 4, true);
        polyscope::getSurfaceMesh()->addVectorQuantity("smoothest 4-field", smoothestField, 4);
      }
      */
    }

    // Register the user callback 
    polyscope::state::userCallback = myCallback;
    
    scatter = new polyscope::Scatterplot();
    std::vector<double> xs;
    std::vector<double> ys;
    for (int i = 0; i < 100; i++) {
      xs.push_back(3 * unitRand() - .5);
      ys.push_back(3 * unitRand() - .5);
    }
    (*scatter).buildScatterplot(xs, ys);
    
    // Give control to the polyscope gui
    polyscope::show();
  }
  return EXIT_SUCCESS;

}
