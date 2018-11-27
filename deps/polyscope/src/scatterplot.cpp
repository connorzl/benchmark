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
}

Scatterplot::~Scatterplot() {
  if (prepared) {
    unprepare();
  }
}

void Scatterplot::computeLogData() {
  logData = data;
  for (size_t i = 0; i < logData.size(); i++) {
    for (size_t j = 0; j < logData[i].size(); j++) {
      logData[i][j] = log2(logData[i][j]);
    }
  }
  computeDataNorm(logData, true);
}

void Scatterplot::computeDataNorm(std::vector<std::vector<double>> &data, bool useLog) { 
    double minVal;
    double maxVal;
    double range;

    std::vector<std::pair<double,double>> dataMinMax_;
    std::vector<std::vector<double>> dataNorm_;
    for (size_t i = 0; i < data.size(); i++) {
      std::vector<double> currData = data[i];
      std::pair<double, double> minmax = robustMinMax(currData);
      minVal = minmax.first;
      maxVal = minmax.second;
      range = maxVal - minVal;
      dataMinMax_.push_back(minmax);

      std::vector<double> currDataNorm;
      for (size_t j = 0; j < currData.size(); j++) {
        currDataNorm.push_back( (currData[j] - minVal) / range * 0.9 + 0.05 );
      }
      dataNorm_.push_back(currDataNorm);
    }

    if (useLog) {
      logDataMinMax = dataMinMax_;
      logDataNorm = dataNorm_;
    } else {
      dataMinMax = dataMinMax_;
      dataNorm = dataNorm_;
    }
}

void Scatterplot::buildScatterplot(std::vector<double>& xs_, std::vector<double>& ys_, const std::vector<double>& zs_) {
    if (xs_.size() != ys_.size() || ys_.size() != xs_.size()) {
      polyscope::error("Data vectors must have same length!");
    }
    data.push_back(xs_);
    data.push_back(ys_);
    xData = 0;
    yData = 1;
    if (zs_.size() != 0) {
      hasZ = true;
      zData = 2;
      data.push_back(zs_);
    }

    computeDataNorm(data);
    computeLogData();
    fillBuffers();
}

void Scatterplot::buildScatterplot(std::vector<std::vector<double>>& data_, std::vector<char*>& labels_) {
  isGeneral = true;
  data = data_;
  labels = labels_;

  computeDataNorm(data);
  computeLogData();

  // handle color coordinate, off by default
  hasZ = false;
  zLabels = labels;
  zLabels.push_back((char*)"None");
  zData = labels.size();
  fillBuffers();
}

void Scatterplot::buildScatterplot(std::vector<std::vector<double>>& data_, std::vector<char*>& labels_,
                      std::vector<std::string> &meshNames_, std::vector<std::string> &filePaths_,
                      std::vector<std::vector<double>> &additionalData_, std::vector<std::string> &additionalDataLabels_,
                      const std::vector<double>& pointSizes_, int numParitions_) {
  // perform some error checking
  if (labels_.size() != data_.size()) {
    error("There must be the same number of data vectors and labels!");
  }
  if (meshNames_.size() != filePaths_.size()) {
    error("There must be the same number of mesh names and filepaths!");
  }
  if (additionalData_.size() != additionalDataLabels_.size()) {
    error("There must be the same number of additional data vectors and labels!");
  }
  for (size_t i = 0; i < data_.size(); i++) {
    if (data_[i].size() != meshNames_.size()) {
      error("There must be the same number of points per data vector as mesh names and filepaths!");
    }
  }
  for (size_t i = 0; i < additionalData_.size(); i++) {
    if (meshNames_.size() != additionalData_[i].size()) {
      error("There must be the same number of points and labels per additional data vector as mesh names and filepaths!");
    }
  }

  // Set properties
  meshInfo = true;
  meshNames = meshNames_;
  filePaths = filePaths_;
  additionalData = additionalData_;
  additionalDataLabels = additionalDataLabels_;
  numPartitions = numParitions_;

  if (pointSizes_.size() != 0) {
    std::pair<double, double> minmax = robustMinMax(pointSizes_);
    pointScaling = pointSizes_;
    for (size_t i = 0; i < pointScaling.size(); i++) {
      pointScaling[i] = (pointScaling[i] - minmax.first) / (minmax.second - minmax.first) * .5 + .5;
    }
  }
  buildScatterplot(data_, labels_);
}

void Scatterplot::updateColormap(const gl::Colormap* newColormap) {
  colormap = newColormap;
  fillBuffers();
}

std::vector<Vector3> Scatterplot::computePointCoords() {
  std::vector<Vector3> coords;
  float doublePi = 2.0f * M_PI;

  std::vector<double> xdataNorm;
  if (logX) {
    xdataNorm = logDataNorm[xData];
  } else {
    xdataNorm = dataNorm[xData];
  }

  std::vector<double> ydataNorm;
  if (logY) {
    ydataNorm = logDataNorm[yData];
  } else {
    ydataNorm = dataNorm[yData];
  }

  float pointRadius;
  for (size_t j = 0; j < xdataNorm.size(); j++) {
    pointRadius = radius;
    if (scalePoints) {
      pointRadius = radius * pointScaling[j];
    } 
    float x0 = xdataNorm[j];
    float y0 = ydataNorm[j];

    for (int i = 0; i < numSides; i++) {
      GLfloat x1 = x0 + ( pointRadius * cos( doublePi * i / float(numSides)) );
      GLfloat y1 = y0 + ( pointRadius * sin( doublePi * i / float(numSides)) );

      int iNext = (i+1) % numSides;
      GLfloat x2 = x0 + ( pointRadius * cos( doublePi * iNext / float(numSides)) );
      GLfloat y2 = y0 + ( pointRadius * sin( doublePi * iNext / float(numSides)) );

      float zCoord = 1.0;
      if (hasZ) {
        zCoord = dataNorm[zData][j];
      } else if (numPartitions > 0) {
        int elemsPerPartition = xdataNorm.size() / numPartitions;
        int currPartition = j / elemsPerPartition;
        zCoord = currPartition / (numPartitions - 1.0);
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
  
  // first, retrieve all the relevant data
  std::vector<double> xdataNorm;
  std::vector<double> xdata;
  std::pair<double,double> xMinMax;
  if (logX) {
    xdataNorm = logDataNorm[xData];
    xdata = logData[xData];
    xMinMax = logDataMinMax[xData];
  } else {
    xdataNorm = dataNorm[xData];
    xdata = data[xData];
    xMinMax = dataMinMax[xData];
  }

  std::vector<double> ydataNorm;
  std::vector<double> ydata;
  std::pair<double,double> yMinMax;
  if (logY) {
    ydataNorm = logDataNorm[yData];
    ydata = logData[yData];
    yMinMax = logDataMinMax[yData];
  } else {
    ydataNorm = dataNorm[yData];
    ydata = data[yData];
    yMinMax = dataMinMax[yData];
  }

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
    for (size_t i = 0; i < xdata.size(); i++) {
      float x_c = xdataNorm[i];
      float y_c = ydataNorm[i];
      float d = sqrt( (x_c - x_p) * (x_c - x_p) + (y_c - y_p) * (y_c - y_p) );
      float r = radius;
      if (scalePoints) {
        r *= pointScaling[i];
      }
      if (d <= radius) {
        if (closestD == -1 || d < closestD) {
          closestD = d;
          index = i;
        } 
      }
    }

    // Show coordinates of selected point
    if (closestD != -1) {     
      float circle_x = xdataNorm[index] * w + x_offset + windowX + ImGui::GetScrollX();
      float circle_y = (1.0 -ydataNorm[index]) * h + windowY;
      float r = radius;
      if (scalePoints) {
        r *= pointScaling[index];
      }
      L->AddCircle(ImVec2{circle_x,circle_y},r*w,ImColor(255, 255, 255, 255), numSides, 3.0);  
      
      ImGui::BeginTooltip();
      if (hasZ) {
        ImGui::Text("%g, %g, %g", xdata[index], ydata[index], data[zData][index]);
      } else {
        ImGui::Text("%g, %g", xdata[index], ydata[index]);
      }
      ImGui::EndTooltip();

      if (ImGui::IsMouseClicked(0)) {
        selectedIndex = index;
      }
    }
  }

  // Keep highlighting the selected point
  if (selectedIndex != -1) {
    float circle_x = xdataNorm[selectedIndex] * w + x_offset + windowX + ImGui::GetScrollX();
    float circle_y = (1.0 - ydataNorm[selectedIndex]) * h + windowY;
    float r = radius;
    if (scalePoints) {
      r *= pointScaling[selectedIndex];
    }
    L->AddCircle(ImVec2{circle_x,circle_y},r*w,ImColor(255, 255, 255, 255), numSides, 3.0);
  }

  // Draw axis stuff
  ImGui::Text("");
  int intervals = 5;
  double xminVal = xMinMax.first;
  double xmaxVal = xMinMax.second;
  double yminVal = yMinMax.first;
  double ymaxVal = yMinMax.second;
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

  // Option to use log axis
  ImGui::Separator();
  bool prevLogX = logX;
  ImGui::Checkbox("Log scale for X-axis", &logX);
  ImGui::SameLine();
  bool prevLogY = logY;
  ImGui::Checkbox("Log scale for Y-axis", &logY);
  if (logX != prevLogX || logY != prevLogY) {
    fillBuffers();
  }

  // Additional information about selected point
  ImGui::Separator();
  if (selectedIndex != -1) {
    ImGui::Text("Point Coordinates");
    ImGui::BulletText("X : %g", xdata[selectedIndex]);
    ImGui::BulletText("Y : %g", ydata[selectedIndex]);
    if (hasZ) {
      ImGui::BulletText("Z : %g", data[zData][selectedIndex]);
    }

    // if there is additional information passed in about the meshes
    if (meshInfo) {
      ImGui::Text("Mesh Name: %s", meshNames[selectedIndex].c_str());
      for (size_t i = 0; i < additionalData.size(); i++) {
        ImGui::BulletText("%s : %g", additionalDataLabels[i].c_str(), additionalData[i][selectedIndex]);
      }
      if (ImGui::Button("Load Mesh")) {
        mesh = new HalfedgeMesh(PolygonSoupMesh(filePaths[selectedIndex]), geom);
        polyscope::registerSurfaceMesh(meshNames[selectedIndex], geom);
        computeDistortion(meshNames[selectedIndex]);
      }
      ImGui::SameLine();
      if (ImGui::Button("Remove Mesh")) {
        polyscope::removeStructure("Surface Mesh", meshNames[selectedIndex]);
      }
    }
  } else {
    ImGui::Text("Click on a point for additional information!");
  }

  // if general scatterplot, allow user to choose which data to plot
  if(isGeneral) {
    ImGui::Separator();
    int xDataBefore = xData;
    int yDataBefore = yData;
    ImGui::Combo("X-Axis Data", &xData, labels.data(), labels.size());
    ImGui::Combo("Y-Axis Data", &yData, labels.data(), labels.size());

    int zDataBefore = zData;
    if (numPartitions == 0) {
      ImGui::Combo("Colors Data", &zData, zLabels.data(), zLabels.size());
    }
    if (xDataBefore != xData || yDataBefore != yData || zDataBefore != zData) {
      hasZ = (strcmp(zLabels[zData], (char*)"None") != 0);
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

  if (!hasZ) {
    return;
  }

  // Draw colormap image
  ImGui::SameLine();
  renderToTexture(1);
  float colorMapW = w / 1.5;
  float colorMapH = h / 10.0;
  ImGui::Image(reinterpret_cast<void*>((size_t)textureInd2) /* yes, really. */, ImVec2(colorMapW, colorMapH), 
              ImVec2(0, 1), ImVec2(1, 0));

  // Draw a cursor popup on mouseover
  if (ImGui::IsItemHovered() && hasZ) {
    double zminVal = dataMinMax[zData].first;
    double zmaxVal = dataMinMax[zData].second;
    // Get mouse x coordinate within image
    float mouseX = ImGui::GetMousePos().x - ImGui::GetCursorScreenPos().x - ImGui::GetScrollX();
    double mouseT = (mouseX - widgetWidth - widgetPadding) / ( colorMapW );
    double val = zminVal + mouseT * (zmaxVal - zminVal);
    ImGui::SetTooltip("%g", val);

    // Draw line showing value of colormap at current mouse position
    ImVec2 imageUpperLeft(ImGui::GetCursorScreenPos().x, ImGui::GetCursorScreenPos().y);
    ImVec2 lineStart(imageUpperLeft.x + mouseX, imageUpperLeft.y - colorMapH - 4);
    ImVec2 lineEnd(imageUpperLeft.x + mouseX, imageUpperLeft.y - 4);
    ImGui::GetWindowDrawList()->AddLine(lineStart, lineEnd,
                                        ImGui::ColorConvertFloat4ToU32(ImVec4(254 / 255., 221 / 255., 66 / 255., 1.0)));
  }
}

void Scatterplot::computeDistortion(std::string meshName) {
  distort = new Distortion(mesh, geom);
  distort->computeAreaScaling();
  polyscope::getSurfaceMesh(meshName)->addQuantity("Area Distortion", distort->areaDistortion);
  distort->computeQuasiConformalError();
  polyscope::getSurfaceMesh(meshName)->addQuantity("Angle Distortion", distort->angleDistortion);
  distort->computeTriangleFlips();
  polyscope::getSurfaceMesh(meshName)->addQuantity("Flipped Triangles", distort->trianglesFlipped);
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

  // Add option to scale points
  if (pointScaling.size() > 0) {
    bool prevScalePoints = scalePoints;
    ImGui::Checkbox("Scale Points by Mesh Sizes", &scalePoints);
    if (prevScalePoints != scalePoints) {
      fillBuffers();
    }
  }
}


} // namespace polyscope
