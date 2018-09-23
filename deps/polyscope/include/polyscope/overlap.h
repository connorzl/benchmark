#pragma once

#include "polyscope/gl/gl_utils.h"
#include "polyscope/polyscope.h"
#include <vector>

namespace polyscope {

// A histogram that shows up in ImGUI
class Overlap {
public:
  Overlap(HalfedgeMesh* mesh, Geometry<Euclidean>* geom);
  ~Overlap();

  // Width = -1 means set automatically
  void buildUI(float width=-1.0);

  float colormapRangeMin, colormapRangeMax; // in DATA values, not [0,1]

private:
  // Manage boundary vertices
  void fillBuffers();

  HalfedgeMesh* mesh;
  Geometry<Euclidean>* geom;
  
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