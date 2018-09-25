#pragma once

#include "polyscope/gl/gl_utils.h"
#include "polyscope/polyscope.h"
#include <vector>

namespace polyscope {

class Overlap {
public:
  Overlap(HalfedgeMesh* mesh, Geometry<Euclidean>* geom);
  ~Overlap();

  void renderToTexture();
  
private:
  // Manage boundary vertices
  void fillBuffers();

  HalfedgeMesh* mesh;
  Geometry<Euclidean>* geom;
  
  void prepare();
  bool prepared = false;
  void unprepare();
  GLuint texDim = 600;
  GLuint framebufferInd, textureInd;
  gl::GLProgram* program = nullptr;
  const gl::Colormap* colormap = &gl::CM_CONST_RED;
};


}; // namespace polyscope