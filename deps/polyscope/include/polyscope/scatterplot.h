#pragma once

#include "polyscope/gl/gl_utils.h"
#include "polyscope/gl/colormap_sets.h"
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include <vector>
#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/polygon_soup_mesh.h"
#include "geometrycentral/distortion.h"

namespace polyscope {

// A Scatterplot that shows up in ImGUI
class Scatterplot {
public:
  Scatterplot();
  ~Scatterplot();

  // regular scatterplot with x,y,z data
  void buildScatterplot(std::vector<double>& xs_, std::vector<double>& ys_, const std::vector<double>& zs_ = {});

  // general scatterplot that allows user to pick what to plot on each axis
  void buildScatterplot(std::vector<std::vector<double>>& data, std::vector<char*>& labels);

  // scatterplot where each point represents a mesh
  void buildScatterplot(std::vector<std::vector<double>>& data, std::vector<char*>& labels,
                        std::vector<std::string> &pointNames, std::vector<std::string> &filePaths,
                        std::vector<std::vector<double>> &additionalData, std::vector<std::string> &additionalDataLabels,
                        const std::vector<double>& pointSizes = {}, int numPartitions_ = 0);

  // scatterplot for comparing two algorithms
  void updateColormap(const gl::Colormap* newColormap);

  // Width = -1 means set automatically
  void buildUI(float width=-1.0);

private:
  std::vector<Vector3> computePointCoords();
  void prepareColorCoords();
  void computeDistortion(std::string meshName);
  void computeDataNorm(std::vector<std::vector<double>> &data, bool useLog = false);
  void computeLogData();
  void fillBuffers();

  float radius = 0.025;
  int numSides = 120;
  int selectedIndex = -1;
  bool hasZ = false;
  int iColorMap = 0;

  // current data index for each axis
  int xData = 0;
  int yData = 0;
  int zData = 0;

  // data, normalized data, and the mins and maxes
  std::vector<std::vector<double>> data;
  std::vector<std::vector<double>> dataNorm;
  std::vector<std::pair<double,double>> dataMinMax;

  // Render to texture
  void renderToTexture(int prog);
  void prepareBuffers(int prog);
  void prepare();
  bool prepared = false;
  void unprepare();
  GLuint texDim = 600;
  GLuint framebufferInd, textureInd;
  gl::GLProgram* scatterProgram = nullptr;
  GLuint framebufferInd2, textureInd2;
  gl::GLProgram* colormapProgram = nullptr;
  const gl::Colormap* colormap = &gl::CM_CONST_RED;

  // point coords
  std::vector<Vector3> pointCoords;

  // colormap coords
  std::vector<Vector3> colorCoordsGradient;
  std::vector<Vector3> colorCoordsSolid;

  // Generalized scatterplot
  bool isGeneral = false;
  std::vector<char*> labels;
  // z label needs its own vector because "None" should always be an option
  std::vector<char*> zLabels;

  // Additional information for points
  bool meshInfo = false;
  std::vector<std::string> meshNames;
  std::vector<std::string> filePaths;
  std::vector<std::vector<double>> additionalData;
  std::vector<std::string> additionalDataLabels;
  bool scalePoints = false;
  std::vector<double> pointScaling;
  Geometry<Euclidean>* geom;
  HalfedgeMesh* mesh;
  Distortion* distort;

  // log scale information
  bool logX = false;
  bool logY = false;
  std::vector<std::vector<double>> logData;
  std::vector<std::vector<double>> logDataNorm;
  std::vector<std::pair<double,double>> logDataMinMax;

  // partitions for multiple data
  int numPartitions = 0;

  // Helpers for building UI
  void buildScatter(float w, float h);
  void buildColormap(float w, float h);
};


}; // namespace polyscope
