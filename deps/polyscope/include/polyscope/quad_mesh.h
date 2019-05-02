#pragma once

#include "polyscope/polyscope.h"
#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/polygon_soup_mesh.h"
#include "geometrycentral/operators.h"
#include "geometrycentral/BHalfedgemesh.h"

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

        // constants
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

        // branch cover quantities
        VertexData<int> singularities;
        HalfedgeData<int> eta;
        int numSingularities = 0;

        // uniformization quantities
        EdgeData<double> edgeLengthsCM;
        HalfedgeData<std::complex<double>> thetaCM;
        HalfedgeData<std::complex<double>> rCM;
        HalfedgeData<double> cmAngles;
        VertexData<double> curvatures;

        // uniformization helpers
        void setupCM();

        // cross field quantities
        FaceData<std::complex<double>> fieldCM;    // field computed on original surface
        std::vector<FaceData<std::complex<double>>> branchCoverFields; // field computed on branch cover
        BranchCoverTopology BC;
        double scale = 1;

        // stripes helpers
        void computeOmega();
        Eigen::SparseMatrix<double> energyMatrix();
        Eigen::SparseMatrix<double> massMatrix();

        // stripes quantities
        std::vector<EdgeData<double>> omega;
        std::vector<VertexData<size_t>> BVertexIndices;

        // texture coordinates quantities
        std::vector<VertexData<std::complex<double>>> psi;
        std::vector<VertexData<double>> coords;   // these are the direct coords 
        FaceData<std::vector<Vector2>> texCoords;    

};