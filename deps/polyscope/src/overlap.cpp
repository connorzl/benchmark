#include "polyscope/overlap.h"
#include "polyscope/gl/shaders/histogram_shaders.h"

#include <algorithm>
#include <limits>

using std::cout;
using std::endl;

namespace polyscope {

static const VertShader OVERLAP_VERT_SHADER =  {
    
    // uniforms
    {
    },

    // attributes
    {
        {"a_coord", GLData::Vector2Float},
    },

    // source
    GLSL(150,
      in vec2 a_coord;
      out float t;

      void main()
      {
          t = 0.5f;
          vec2 scaledCoord = vec2(a_coord.x, a_coord.y * .85);
          gl_Position = vec4(2.*scaledCoord - vec2(1.0, 1.0),0.,1.);
      }
    )
};

static const FragShader OVERLAP_FRAG_SHADER = {
    
    // uniforms
    {
    }, 

    // attributes
    {
    },
    
    // textures 
    {
        {"t_colormap", 1}
    },
    
    // output location
    "outputF",
    
    // source 
    GLSL(330,
      in float t;
      uniform sampler1D t_colormap;
      layout(location = 0) out vec4 outputF;

      void main()
      {
        outputF = vec4(texture(t_colormap, t).rgb, 1.0);
      }
    )
};


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

  // Set the list of draw buffers
  GLenum drawBuffers[1] = {GL_COLOR_ATTACHMENT0};
  glDrawBuffers(1, drawBuffers);

  if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
    throw std::logic_error("Overlap framebuffer problem");
  }

  // Create the program
  program = new gl::GLProgram(&OVERLAP_VERT_SHADER, &OVERLAP_FRAG_SHADER, gl::DrawMode::Triangles);

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

  for (FacePtr f : mesh->faces()) {
      HalfedgePtr he = f.halfedge();
      coords.push_back(geom->paramCoords[he.prev()]);
      coords.push_back(geom->paramCoords[he]);
      coords.push_back(geom->paramCoords[he.next()]);
  }

  program->setAttribute("a_coord", coords);
  program->setTextureFromColormap("t_colormap", *colormap, true);
}

void Overlap::renderToTexture() {
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
/*
  unsigned char* pixels = (unsigned char*) malloc(texDim * texDim * 4);
  glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
  for (GLuint i = 0; i < texDim; i++) {
      for (GLuint j = 0; j < texDim; j++) {
          std::cout << pixels[4 * (i * texDim + j)] << std::endl;
      }
  }
*/
  bindDefaultBuffer();
}

}; // namespace polyscope