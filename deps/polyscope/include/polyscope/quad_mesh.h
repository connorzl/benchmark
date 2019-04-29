#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/polygon_soup_mesh.h"
#include "geometrycentral/operators.h"
#include "geometrycentral/BHalfedgemesh.h"

#include "polyscope/polyscope.h"
#include "polyscope/gl/gl_utils.h"
#include "polyscope/gl/colormap_sets.h"
#include "polyscope/gl/shaders/scatterplot_shaders.h"

using namespace geometrycentral;

class QuadMesh {
    public:
        QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g);
        void computeCrossField();
        void computeSingularities();
        void computeBranchCover(bool improve=true);
        
        // cone metric
        void uniformize();
        void computeCrossFieldCM();  // not really needed anymore 
        void computeCrossFieldCMBranchCover();
        
        // stripes
        void computeStripes();
        void textureCoordinates();

        // visualization
        void visualize();

    private:
        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;

        int nPowerIterations = 20;
        double eps = std::pow(10.,-8.);
        int n = 4;

        // smoothest field quantities
        HalfedgeData<std::complex<double>> theta;
        HalfedgeData<std::complex<double>> r;
        FaceData<std::complex<double>> field;

        // smoothest field helpers
        void setup();
        Eigen::SparseMatrix<std::complex<double>> assembleM();
        Eigen::SparseMatrix<std::complex<double>> assembleA();
        void computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A);

        VertexData<int> singularities;
        HalfedgeData<int> eta;
        int numSingularities = 0;

        // uniformization results
        EdgeData<double> edgeLengthsCM;
        HalfedgeData<std::complex<double>> thetaCM;
        HalfedgeData<std::complex<double>> rCM;
        HalfedgeData<double> cmAngles;
        VertexData<double> curvatures;
        
        // cross field
        FaceData<std::complex<double>> fieldCM;    // field computed on original surface
        std::vector<FaceData<std::complex<double>>> branchCoverFields; // field computed on branch cover
        BranchCoverTopology BC;
        double scale = 1;

        // stripes
        void computeOmega();
        Eigen::SparseMatrix<double> energyMatrix();
        Eigen::SparseMatrix<double> massMatrix();
        std::vector<EdgeData<double>> omega;
        std::vector<VertexData<size_t>> BVertexIndices;

        // texture coordinates
        std::vector<VertexData<std::complex<double>>> psi;
        std::vector<VertexData<double>> coords;

        // helpers
        void setupCM();
};