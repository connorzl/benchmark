#pragma once

#include "polyscope/gl/gl_utils.h"
#include "polyscope/gl/colormap_sets.h"
#include <iomanip> // setprecision
#include <sstream> // stringstream

#include <vector>


namespace polyscope {

// A Scatterplot that shows up in ImGUI
class Scatterplot {
public:
  Scatterplot();
  ~Scatterplot();

  void buildScatterplot(std::vector<double>& xs_, std::vector<double>& ys_, const std::vector<double>& zs_ = {});
  void buildScatterPlot(std::vector<std::vector<double>>& data, std::vector<char*>& labels);
  void updateColormap(const gl::Colormap* newColormap);
  std::vector<Vector3> computePointCoords();
  void prepareColorCoords();

  // Width = -1 means set automatically
  void buildUI(float width=-1.0);

private:
  void fillBuffers();

  double xminVal = 0;
  double xmaxVal = 0;
  double yminVal = 0;
  double ymaxVal = 0;
  double zminVal = 0;
  double zmaxVal = 0;
  float radius = 0.025;
  int numSides = 120;
  
  bool toggleInfo = false;
  int selectedIndex = -1;

  bool hasZ = false;
  int iColorMap = 0;

  std::vector<double> xs;
  std::vector<double> xs_norm;
  std::vector<double> ys;
  std::vector<double> ys_norm;
  std::vector<double> zs;
  std::vector<double> zs_norm;

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

  // used to render the colormap
  std::vector<Vector3> colorCoordsGradient;
  std::vector<Vector3> colorCoordsSolid;

  // Generalized scatterplot
  bool isGeneral = false;
  int xData = 0;
  int yData = 0;
  int zData = 0;
  std::vector<std::vector<double>> generalData;
  std::vector<std::vector<double>> generalDataNorm;
  std::vector<char*> generalLabels;
  std::vector<char*> generalLabelsZ;
  std::vector<std::pair<double,double>> generalMinMax;

  // Helpers for building UI
  void buildScatter(float w, float h);
  void buildColormap(float w, float h);
};


}; // namespace polyscope
