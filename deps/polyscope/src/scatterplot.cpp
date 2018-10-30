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

Scatterplot::~Scatterplot() {
  if (prepared) {
    unprepare();
  }
}

void Scatterplot::buildScatterplot(std::vector<double>& xs_, std::vector<double>& ys_, const std::vector<double>& zs_) {
    xs = xs_;
    ys = ys_;
    xs_norm = xs_;
    ys_norm = ys_;

    if (zs_.size() != 0) {
      hasZ = true;
      zs = zs_;
      zs_norm = zs_;
    }

    // normalize to [0, 1] range
    std::pair<double, double> xminmax = robustMinMax(xs_);
    xminVal = xminmax.first;
    xmaxVal = xminmax.second;
    double xrange = xmaxVal - xminVal;

    std::pair<double, double> yminmax = robustMinMax(ys_);
    yminVal = yminmax.first;
    ymaxVal = yminmax.second;
    double yrange = ymaxVal - yminVal;

    double zrange = 0;
    if (hasZ) {
      std::pair<double, double> zminmax = robustMinMax(zs_);
      zminVal = zminmax.first;
      zmaxVal = zminmax.second;
      zrange = zmaxVal - zminVal;
    }

    for (size_t i = 0; i < xs_norm.size(); i++) {
        xs_norm[i] = (xs_norm[i] - xminVal) / xrange * 0.9 + 0.05;
        ys_norm[i] = (ys_norm[i] - yminVal) / yrange * 0.9 + 0.05;

        if (hasZ) {
          zs_norm[i] = (zs_norm[i] - zminVal) / zrange * 0.9 + 0.05;
        }
    }
    
    fillBuffers();
}

void Scatterplot::buildScatterPlot(std::vector<std::vector<double>>& data_, std::vector<char*>& labels_) {
  isGeneral = true;
  generalData = data_;
  generalLabels = labels_;

  double minVal;
  double maxVal;
  double range;
  for (size_t i = 0; i < generalData.size(); i++) {
    std::vector<double> currData = generalData[i];
    std::pair<double, double> minmax = robustMinMax(currData);
    minVal = minmax.first;
    maxVal = minmax.second;
    range = maxVal - minVal;
    generalMinMax.push_back(minmax);

    std::vector<double> currDataNorm;
    for (size_t j = 0; j < currData.size(); j++) {
      currDataNorm.push_back( (currData[j] - minVal) / range * 0.9 + 0.05 );
    }
    generalDataNorm.push_back(currDataNorm);
  }

  // handle color coordinate, off by default
  hasZ = false;
  generalLabelsZ = generalLabels;
  generalLabelsZ.push_back((char*)"None");
  zData = generalLabels.size();

  xminVal = generalMinMax[xData].first;
  xmaxVal = generalMinMax[xData].second;
  yminVal = generalMinMax[yData].first;
  ymaxVal = generalMinMax[yData].second;

  fillBuffers();
}

void Scatterplot::updateColormap(const gl::Colormap* newColormap) {
  colormap = newColormap;
  fillBuffers();
}

std::vector<Vector3> Scatterplot::computePointCoords() {
  std::vector<Vector3> coords;
  float doublePi = 2.0f * M_PI;

  for (size_t j = 0; j < xs_norm.size(); j++) {
    float x0 = xs_norm[j];
    float y0 = ys_norm[j];

    for (int i = 0; i < numSides; i++) {
      GLfloat x1 = x0 + ( radius * cos( doublePi * i / float(numSides)) );
      GLfloat y1 = y0 + ( radius * sin( doublePi * i / float(numSides)) );

      int iNext = (i+1) % numSides;
      GLfloat x2 = x0 + ( radius * cos( doublePi * iNext / float(numSides)) );
      GLfloat y2 = y0 + ( radius * sin( doublePi * iNext / float(numSides)) );

      float zCoord = 1.0;
      if (hasZ) {
        zCoord = zs_norm[j];
      }
      coords.push_back(Vector3{x0,y0,zCoord});
      coords.push_back(Vector3{x1,y1,zCoord});
      coords.push_back(Vector3{x2,y2,zCoord});
    }
  }

  return coords;
}

void Scatterplot::prepareColorCoords() {
  float stepsize = 0.001;
  int numsteps = 1.0 / stepsize;
  float leftX, rightX, leftY, rightY;
  float color;

  for (int i = 0; i < numsteps; i++) {
    leftX = i * stepsize;
    rightX = i * (stepsize+1);
    leftY = 1.0;
    rightY = 1.0;
    
    color = i * stepsize;
    // = Lower triangle (lower left, lower right, upper left)
    colorCoordsGradient.push_back(Vector3{leftX, 0.0, color});
    colorCoordsGradient.push_back(Vector3{rightX, 0.0, color});
    colorCoordsGradient.push_back(Vector3{leftX, leftY, color});

    // = Upper triangle (lower right, upper right, upper left)
    colorCoordsGradient.push_back(Vector3{rightX, 0.0, color});
    colorCoordsGradient.push_back(Vector3{rightX, rightY, color});
    colorCoordsGradient.push_back(Vector3{leftX, leftY, color});

    color = 1.0;
    // = Lower triangle (lower left, lower right, upper left)
    colorCoordsSolid.push_back(Vector3{leftX, 0.0, color});
    colorCoordsSolid.push_back(Vector3{rightX, 0.0, color});
    colorCoordsSolid.push_back(Vector3{leftX, leftY, color});

    // = Upper triangle (lower right, upper right, upper left)
    colorCoordsSolid.push_back(Vector3{rightX, 0.0, color});
    colorCoordsSolid.push_back(Vector3{rightX, rightY, color});
    colorCoordsSolid.push_back(Vector3{leftX, leftY, color});
  }
}

void Scatterplot::fillBuffers() {
  if (isGeneral) {
    xs = generalData[xData];
    xs_norm = generalDataNorm[xData];
    ys = generalData[yData];
    ys_norm = generalDataNorm[yData];
    if (hasZ) {
      zs = generalData[zData];
      zs_norm = generalDataNorm[zData];
    }
  }

  pointCoords = computePointCoords();
  scatterProgram->setAttribute("a_coord", pointCoords);
  scatterProgram->setTextureFromColormap("t_colormap", *colormap, true);

  if (hasZ) {
    colormapProgram->setAttribute("a_coord", colorCoordsGradient);
  } else {
    colormapProgram->setAttribute("a_coord", colorCoordsSolid);
  }
  colormapProgram->setTextureFromColormap("t_colormap", *colormap, true);
}

void Scatterplot::prepareBuffers(int prog) {
  GLuint frameInd, texInd;
  // Generate a framebuffer to hold our texture
  glGenFramebuffers(1, &frameInd);
  glBindFramebuffer(GL_FRAMEBUFFER, frameInd);

  // Create the texture
  glGenTextures(1, &texInd);

  // Bind to the new texture so we can set things about it
  glBindTexture(GL_TEXTURE_2D, texInd);

  // Configure setttings
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texDim, texDim, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  // Attach the texture to the framebuffer
  glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texInd, 0);
  GLenum drawBuffers[1] = {GL_COLOR_ATTACHMENT0};
  glDrawBuffers(1, drawBuffers);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    throw std::logic_error("Scatterplot framebuffer problem");
  }

  if (prog == 0) {
    framebufferInd = frameInd;
    textureInd = texInd;
  } else {
    framebufferInd2 = frameInd;
    textureInd2 = texInd;
  }
}
void Scatterplot::prepare() {
  if (prepared) {
    unprepare();
  }
  
  prepareBuffers(0);
  prepareBuffers(1);
  prepareColorCoords();
  
  // Create the programs
  scatterProgram = new gl::GLProgram(&SCATTERPLOT_VERT_SHADER, &SCATTERPLOT_FRAG_SHADER, gl::DrawMode::Triangles);
  colormapProgram = new gl::GLProgram(&SCATTERPLOT_VERT_SHADER, &SCATTERPLOT_FRAG_SHADER, gl::DrawMode::Triangles);
  prepared = true;

}

void Scatterplot::unprepare() {
  safeDelete(scatterProgram);
  glDeleteTextures(1, &textureInd);
  glDeleteFramebuffers(1, &framebufferInd);

  safeDelete(colormapProgram);
  glDeleteTextures(1, &textureInd2);
  glDeleteFramebuffers(1, &framebufferInd2);
  prepared = false;
}

void Scatterplot::renderToTexture(int prognum) {
  GLuint frameInd, texInd;
  gl::GLProgram* prog;
  if (prognum == 0) {
    frameInd = framebufferInd;
    texInd = textureInd;
    prog = scatterProgram;
  } else {
    frameInd = framebufferInd2;
    texInd = textureInd2;
    prog = colormapProgram;
  }

  // Bind to the texture buffer
  glBindFramebuffer(GL_FRAMEBUFFER, frameInd);

  // Bind to the new texture so we can do things
  glBindTexture(GL_TEXTURE_2D, texInd);
 
  GLenum drawBuffers[1] = {GL_COLOR_ATTACHMENT0};
  glDrawBuffers(1, drawBuffers);

  // Make sure we render to the whole buffer
  glViewport(0, 0, texDim, texDim);
  glClearColor(0.0, 0.0, 0.0, 0.2);
  glClear(GL_COLOR_BUFFER_BIT);

  // Draw
  prog->draw();
  bindDefaultBuffer();
}

void Scatterplot::buildScatter(float w, float h) {
  float windowX = ImGui::GetCursorScreenPos().x;
  float windowY = ImGui::GetCursorScreenPos().y;

  // Render image
  renderToTexture(0);
  float x_offset = 50;
  
  // Shift the plot horizontally to leave room for y-axis labels
  ImGui::SetCursorScreenPos(ImVec2{windowX+x_offset, windowY});
  ImGui::Image(reinterpret_cast<void*>((size_t)textureInd) /* yes, really. */, ImVec2(w, h), ImVec2(0, 1),ImVec2(1, 0));
  ImGui::SetCursorScreenPos(ImVec2{ImGui::GetCursorScreenPos().x-x_offset, ImGui::GetCursorScreenPos().y});
  
  ImDrawList* L = ImGui::GetWindowDrawList();
  // Draw a cursor popup on mouseover
  if (ImGui::IsItemHovered()) {
    // Get mouse x coordinate within image
    float mouseX = ImGui::GetMousePos().x - windowX - ImGui::GetScrollX() - x_offset;
    double x_p = mouseX / w;

    // Get mouse y coordinate within image
    float mouseY = ImGui::GetMousePos().y - windowY; //- ImGui::GetScrollY();
    double y_p = 1.0 - mouseY / h;

    float closestD = -1;
    int index = 0;
    // Check if mouse is selecting over a scatterplot point
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

    // Show coordinates of selected point
    if (closestD != -1) { 
      float circle_x = xs_norm[index] * w + x_offset + windowX + ImGui::GetScrollX();
      float circle_y = (1.0 - ys_norm[index]) * h + windowY;
      L->AddCircle(ImVec2{circle_x,circle_y},radius*w,ImColor(255, 255, 255, 255), numSides, 3.0);
      
      ImGui::BeginTooltip();
      if (hasZ) {
        ImGui::Text("%g, %g, %g", xs[index], ys[index], zs[index]);
      } else {
        ImGui::Text("%g, %g", xs[index], ys[index]);
      }
      ImGui::EndTooltip();

      if (ImGui::IsMouseClicked(0)) {
        selectedIndex = index;
      }
    }
  }

  // Draw axis stuff
  ImGui::Text("");
  int intervals = 5;
  float stepsizeX = (xmaxVal - xminVal) / intervals;
  float stepsizeY = (ymaxVal - yminVal) / intervals;
  windowX = windowX + x_offset;
  float y_offset = 15;
  float y_offset2 = 5;
  float labelSize = 30;
  float labelHeight = 15;

  L->AddLine(ImVec2{windowX, windowY+h}, ImVec2{windowX+w, windowY+h}, ImColor(255, 255, 255, 255), 2.0);
  L->AddLine(ImVec2{windowX, windowY}, ImVec2{windowX, windowY+h}, ImColor(255, 255, 255, 255), 2.0);
  for (int i = 0; i <= intervals; i++) {
    float currT = (w-labelSize) / intervals * i;
    float f = xminVal + (xmaxVal - xminVal) * (currT / (w-labelSize));
    std::stringstream stream;
    stream << std::fixed << std::setprecision(3) << f;
    std::string s = stream.str();
    L->AddText(ImVec2{windowX + currT, windowY+h+y_offset2},ImColor(255, 255, 255, 255), s.c_str());

    currT = (h-labelHeight) / intervals * i;
    f = yminVal + (ymaxVal - yminVal) * (currT / (h-labelHeight));
    std::stringstream().swap(stream);
    stream << std::fixed << std::setprecision(3) << f;
    s = stream.str();
    L->AddText(ImVec2{windowX-x_offset, windowY + h - currT - y_offset},ImColor(255, 255, 255, 255), s.c_str());
  }  

  // Additional information about selected point
  ImGui::Separator();
  if (selectedIndex != -1) {
    ImGui::BulletText("X Coordinate: %g", xs[selectedIndex]);
    ImGui::BulletText("Y Coordinate: %g", ys[selectedIndex]);
    if (hasZ) {
      ImGui::BulletText("Z Coordinate: %g", zs[selectedIndex]);
    }

    // if there is additional information passed in...
    if (ImGui::Button("Load Mesh")) {
      std::cout<<"weijaers"<<std::endl;
    }
  } else {
    ImGui::Text("Click on a point for additional information!");
  }

  // if general scatterplot, allow user to choose which data to plot
  if(isGeneral) {
    ImGui::Separator();
    int xDataBefore = xData;
    int yDataBefore = yData;
    int zDataBefore = zData;
    ImGui::Combo("X-Axis Data", &xData, generalLabels.data(), generalLabels.size());
    ImGui::Combo("Y-Axis Data", &yData, generalLabels.data(), generalLabels.size());
    ImGui::Combo("Colors Data", &zData, generalLabelsZ.data(), generalLabelsZ.size());
    if (xDataBefore != xData || yDataBefore != yData || zDataBefore != zData) {
      if (strcmp(generalLabelsZ[zData], (char*)"None") == 0) {
        hasZ = false;
      } else {
        hasZ = true;
        zminVal = generalMinMax[zData].first;
        zmaxVal = generalMinMax[zData].second;
      }

      xminVal = generalMinMax[xData].first;
      xmaxVal = generalMinMax[xData].second;
      yminVal = generalMinMax[yData].first;
      ymaxVal = generalMinMax[yData].second;
      fillBuffers();
    }
  }
}

void Scatterplot::buildColormap(float w, float h) {
  // Set colormap
  ImGui::Separator();
  float widgetWidth = 100;
  float widgetPadding = 8.6;
  int iColormapBefore = iColorMap;
  ImGui::Combo("##Color Style", &iColorMap, polyscope::gl::quantitativeColormapNames,
                IM_ARRAYSIZE(polyscope::gl::quantitativeColormapNames));
  if (iColorMap != iColormapBefore) {
    updateColormap(polyscope::gl::quantitativeColormaps[iColorMap]);
  }
  
  // Draw colormap image
  ImGui::SameLine();
  renderToTexture(1);
  ImGui::Image(reinterpret_cast<void*>((size_t)textureInd2) /* yes, really. */, ImVec2(w/1.5, h/10.0), ImVec2(0, 1),
               ImVec2(1, 0));

  // Draw a cursor popup on mouseover
  if (ImGui::IsItemHovered() && hasZ) {
    // Get mouse x coodinate within image
    float mouseX = ImGui::GetMousePos().x - ImGui::GetCursorScreenPos().x - ImGui::GetScrollX();
    double mouseT = (mouseX - widgetWidth - widgetPadding) / ( w/1.5 );
    double val = zminVal + mouseT * (zmaxVal - zminVal);
    ImGui::SetTooltip("%g", val);

    // Draw line showing value of colormap at current mouse position
    ImVec2 imageUpperLeft(ImGui::GetCursorScreenPos().x, ImGui::GetCursorScreenPos().y);
    ImVec2 lineStart(imageUpperLeft.x + mouseX, imageUpperLeft.y - h/10.0 - 4);
    ImVec2 lineEnd(imageUpperLeft.x + mouseX, imageUpperLeft.y - 4);
    ImGui::GetWindowDrawList()->AddLine(lineStart, lineEnd,
                                        ImGui::ColorConvertFloat4ToU32(ImVec4(254 / 255., 221 / 255., 66 / 255., 1.0)));
  }
}

void Scatterplot::buildUI(float width) {
  // Compute size for image
  float w = width;
  if (w == -1.0) {
    w = .8 * ImGui::GetWindowWidth();
  }
  float h = w;// aspect;

  // build scatterplot
  buildScatter(w, h);

  // build colormap
  buildColormap(w, h);

  // Add slider for point radius
  float prevRadius = radius;
  ImGui::SliderFloat("Point Radius", &radius, 0.005, .075, "%.5f", 1.);
  if (radius != prevRadius) {
    fillBuffers();
  }
}


} // namespace polyscope
