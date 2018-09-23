#include "polyscope/overlap.h"

//#include "polyscope/affine_remapper.h"
#include "polyscope/gl/shaders/histogram_shaders.h"

//#include "imgui.h"

#include <algorithm>
#include <limits>

using std::cout;
using std::endl;

namespace polyscope {

Overlap::Overlap(HalfedgeMesh* mesh, Geometry<Euclidean>* geom) : mesh(mesh), geom(geom) {
    prepare();
    fillBuffers();
}

Overlap::~Overlap() {
    if (prepared) {
        unprepare();
    }
}

void Overlap::prepare() {

  if (prepared) {
    unprepare();
  }

  // Generate a framebuffer to hold our texture
  glGenFramebuffers(1, &framebufferInd);
  glBindFramebuffer(GL_FRAMEBUFFER, framebufferInd);

  // Create the texture
  glGenTextures(1, &textureInd);

  // Bind to the new texture so we can set things about it
  glBindTexture(GL_TEXTURE_2D, textureInd);

  // Configure setttings
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texDim, texDim, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  // Attach the texture to the framebuffer
  glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, textureInd, 0);
  GLenum drawBuffers[1] = {GL_COLOR_ATTACHMENT0};
  glDrawBuffers(1, drawBuffers);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    throw std::logic_error("Histogram framebuffer problem");
  }

  // Create the program, might have to change these shaders
  program = new gl::GLProgram(&HISTOGRAM_VERT_SHADER, &HISTORGRAM_FRAG_SHADER, gl::DrawMode::Triangles);

  prepared = true;
}

void Overlap::unprepare() {
  safeDelete(program);
  glDeleteTextures(1, &textureInd);
  glDeleteFramebuffers(1, &framebufferInd);
  prepared = false;
}

void Overlap::fillBuffers() {
  // Push to buffer
  std::vector<Vector2> coords;

  // loop through "imaginary" boundary edges
  for (size_t i = 0; i < mesh->nImaginaryHalfedges(); i++) {
    HalfedgePtr he = mesh->imaginaryHalfedge(i).twin();
    Vector2 A = geom->paramCoords[he.next()];
    coords.push_back(A);
  }

  program->setAttribute("a_coord", coords);
  program->setTextureFromColormap("t_colormap", *colormap, true);
}

void Overlap::renderToTexture() {

}


}; // namespace polyscope