#include "polyscope/Scatterplot.h"

#include "polyscope/affine_remapper.h"
#include "polyscope/gl/shaders/scatterplot_shaders.h"
#include "polyscope/polyscope.h"

#include "imgui.h"

#include <algorithm>
#include <limits>
#include <math.h>

using std::cout;
using std::endl;

namespace polyscope {

Scatterplot::Scatterplot() {
  prepare();
  fillBuffers();
}

Scatterplot::Scatterplot(const std::vector<double>& xs_, const std::vector<double>& ys_) {
  prepare();
  buildScatterplot(xs_, ys_);
}

Scatterplot::~Scatterplot() {
  if (prepared) {
    unprepare();
  }
}

void Scatterplot::buildScatterplot(const std::vector<double>& xs_, const std::vector<double>& ys_) {
    xs = xs_;
    ys = ys_;
    xs_norm = xs_;
    ys_norm = ys_;

    // normalize xs to [0, 1] range
    std::pair<double, double> xminmax = robustMinMax(xs_);
    xminVal = xminmax.first;
    xmaxVal = xminmax.second;
    double xrange = xmaxVal - xminVal;

    std::pair<double, double> yminmax = robustMinMax(ys_);
    yminVal = yminmax.first;
    ymaxVal = yminmax.second;
    double yrange = ymaxVal - yminVal;

    for (size_t i = 0; i < xs_norm.size(); i++) {
        xs_norm[i] = (xs_norm[i] - xminVal) / xrange * 0.9 + 0.05;
        ys_norm[i] = (ys_norm[i] - yminVal) / yrange * 0.9 + 0.05;
    }
  
    fillBuffers();
}

void Scatterplot::updateColormap(const gl::Colormap* newColormap) {
  colormap = newColormap;
  fillBuffers();
}

std::vector<Vector2> Scatterplot::computePointCoords(float radius, int numberOfSides) {
  std::vector<Vector2> coords;
  float doublePi = 2.0f * M_PI;

  for (size_t j = 0; j < xs_norm.size(); j++) {
    float x0 = xs_norm[j];
    float y0 = ys_norm[j];

    for (int i = 0; i < numberOfSides; i++) {
      GLfloat x1 = x0 + ( radius * cos( doublePi * i / float(numberOfSides)) );
      GLfloat y1 = y0 + ( radius * sin( doublePi * i / float(numberOfSides)) );

      int iNext = (i+1) % numberOfSides;
      GLfloat x2 = x0 + ( radius * cos( doublePi * iNext / float(numberOfSides)) );
      GLfloat y2 = y0 + ( radius * sin( doublePi * iNext / float(numberOfSides)) );

      coords.push_back(Vector2{x0,y0});
      coords.push_back(Vector2{x1,y1});
      coords.push_back(Vector2{x2,y2});
    }
  }

  return coords;
}

void Scatterplot::fillBuffers() {
  // draw circles using xs_norm and ys_norm...
  int numSides = 120;
  std::vector<Vector2> coords = computePointCoords(radius, numSides);

  program->setAttribute("a_coord", coords);
  program->setTextureFromColormap("t_colormap", *colormap, true);
}

void Scatterplot::prepare() {

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
    throw std::logic_error("Scatterplot framebuffer problem");
  }

  // Create the program
  program = new gl::GLProgram(&SCATTERPLOT_VERT_SHADER, &SCATTERPLOT_FRAG_SHADER, gl::DrawMode::Triangles);

  prepared = true;
}

void Scatterplot::unprepare() {
  safeDelete(program);
  glDeleteTextures(1, &textureInd);
  glDeleteFramebuffers(1, &framebufferInd);
  prepared = false;
}

void Scatterplot::renderToTexture() {
  if (radiusChanged) {
    fillBuffers();
  }

  // Bind to the texture buffer
  glBindFramebuffer(GL_FRAMEBUFFER, framebufferInd);

  // Bind to the new texture so we can do things
  glBindTexture(GL_TEXTURE_2D, textureInd);


  GLenum drawBuffers[1] = {GL_COLOR_ATTACHMENT0};
  glDrawBuffers(1, drawBuffers);

  // Make sure we render to the whole buffer
  glViewport(0, 0, texDim, texDim);
  glClearColor(0.0, 0.0, 0.0, 0.2);
  glClear(GL_COLOR_BUFFER_BIT);

  // Draw
  program->draw();

  bindDefaultBuffer();
}


void Scatterplot::buildUI(float width) {
  renderToTexture();

  // Compute size for image
  float aspect = 3.0;
  float w = width;
  if (w == -1.0) {
    w = .8 * ImGui::GetWindowWidth();
  }
  float h = w;// / aspect;

  // Render image
  ImGui::Image(reinterpret_cast<void*>((size_t)textureInd) /* yes, really. */, ImVec2(w, h), ImVec2(0, 1),
               ImVec2(1, 0));

  // Draw a cursor popup on mouseover
  if (ImGui::IsItemHovered()) {

    // Get mouse x coodinate within image
    float mouseX = ImGui::GetMousePos().x - ImGui::GetCursorScreenPos().x - ImGui::GetScrollX();
    double mouseT = mouseX / w;
    double valX = xminVal + mouseT * (xmaxVal - xminVal);

    // Get mouse y coordinate within image
    float mouseY = ImGui::GetCursorScreenPos().y - ImGui::GetMousePos().y - ImGui::GetScrollY();
    mouseT = mouseY / h;
    double valY = yminVal + mouseT * (ymaxVal - yminVal);
    
    float x_p = mouseX / w;
    float y_p = mouseY / h;
    float closestD = -1;
    size_t index = 0;

    for (size_t i = 0; i < xs.size(); i++) {
      float x_c = xs_norm[i];
      float y_c = ys_norm[i];
      float d = sqrt( (x_c - x_p) * (x_c - x_p) + (y_c - y_p) * (y_c - y_p) );
      if (d <= radius) {
        if (closestD == -1 || d < closestD) {
          closestD = d;
          index = i;
        } 
      }
    }

    if (closestD != -1) {
      ImGui::SetTooltip("%g, %g", xs[index], ys[index]);
    }
    
    /*
    // Draw line
    ImVec2 imageUpperLeft(ImGui::GetCursorScreenPos().x, ImGui::GetCursorScreenPos().y);
    ImVec2 lineStart(imageUpperLeft.x + mouseX, imageUpperLeft.y - h - 3);
    ImVec2 lineEnd(imageUpperLeft.x + mouseX, imageUpperLeft.y - 4);
    ImGui::GetWindowDrawList()->AddLine(lineStart, lineEnd,
                                        ImGui::ColorConvertFloat4ToU32(ImVec4(254 / 255., 221 / 255., 66 / 255., 1.0)));

    ImVec2 lineStart2(imageUpperLeft.x, imageUpperLeft.y - mouseY);
    ImVec2 lineEnd2(imageUpperLeft.x + w, imageUpperLeft.y - mouseY);
    ImGui::GetWindowDrawList()->AddLine(lineStart2, lineEnd2,
                                        ImGui::ColorConvertFloat4ToU32(ImVec4(254 / 255., 221 / 255., 66 / 255., 1.0)));
    */
  }

  // Slider for point radius
  float prevRadius = radius;
  ImGui::SliderFloat("Point Radius", &radius, 0.01, .075, "%.5f", 1.);
  radiusChanged = (radius != prevRadius);
}


} // namespace polyscope
