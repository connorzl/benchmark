#pragma once

#include "polyscope/gl/gl_utils.h"

#include <vector>


namespace polyscope {

// A Scatterplot that shows up in ImGUI
class Scatterplot {
public:
  Scatterplot();
  Scatterplot(const std::vector<double>& xs_, const std::vector<double>& ys_);

  ~Scatterplot();

  void buildScatterplot(const std::vector<double>& xs_, const std::vector<double>& ys_);
  void updateColormap(const gl::Colormap* newColormap);
  std::vector<Vector2> computePointCoords(float radius, int numberOfSides);

  // Width = -1 means set automatically
  void buildUI(float width=-1.0);

private:
  // Manage the actual Scatterplot
  void fillBuffers();

  double xminVal;
  double xmaxVal;
  double yminVal;
  double ymaxVal;
  float radius = 0.025;
  bool radiusChanged = false;

  std::vector<double> xs;
  std::vector<double> xs_norm;
  std::vector<double> ys;
  std::vector<double> ys_norm;

  // Render to texture
  void renderToTexture();
  void prepare();
  bool prepared = false;
  void unprepare();
  GLuint texDim = 600;
  GLuint framebufferInd, textureInd;
  gl::GLProgram* program = nullptr;
  const gl::Colormap* colormap = &gl::CM_CONST_RED;
};


}; // namespace polyscope
