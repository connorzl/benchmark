#pragma once

#include "polyscope/polyscope.h"
#include "polyscope/hodge_decomposition.h"
#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/polygon_soup_mesh.h"
#include "geometrycentral/operators.h"
#include "geometrycentral/BHalfedgemesh.h"
#include "geometrycentral/linear_solvers.h"

using namespace geometrycentral;

class QuadMesh {
    public:
        QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g);
        void computeCrossField(bool isCM = false);
        void computeSingularities();
        void computeBranchCover(bool improve=true);
        
        // cone metric
        void uniformize();
        void computeCrossFieldCMBranchCover(std::complex<double> init = std::complex<double>(1,0), double scale = 100);
        
        // stripes
        double computeStripes();
        void optimizeSimpleLocally();
        void optimizeSimpleGlobally();
        void optimizeHarmonic(bool visualize = false);
        bool textureCoordinates();

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
        Eigen::SparseMatrix<std::complex<double>> assembleM(bool isCM = false);
        Eigen::SparseMatrix<std::complex<double>> assembleA(bool isCM = false);
        void computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A, bool isCM = false);

        // branch cover quantities
        VertexData<int> singularities;
        HalfedgeData<int> eta;
        int numSingularities = 0;

        // uniformization quantities
        EdgeData<double> edgeLengthsCM;
        HalfedgeData<std::complex<double>> thetaCM;
        HalfedgeData<std::complex<double>> rCM;
        FaceData<std::complex<double>> fieldCM;
        HalfedgeData<double> cmAngles;
        FaceData<double> cmAreas;
        VertexData<double> curvatures;

        // uniformization helpers
        void setupCM();
        double updateAreas();
        void uniformizeBoundary();

        // cross field quantities
        std::vector<FaceData<std::complex<double>>> branchCoverFields; // field computed on branch cover
        BranchCoverTopology BC;
        std::vector<FaceData<std::complex<double>>> xBasis;

        // stripes helpers
        void computeOmega(double scale);
        Eigen::SparseMatrix<double> energyMatrix();
        Eigen::SparseMatrix<double> energyMatrix2();
        Eigen::SparseMatrix<double> massMatrix();

        // stripes quantities
        std::vector<EdgeData<double>> omega;
        std::vector<VertexData<size_t>> BVertexIndices;

        // texture coordinates helpers
        std::complex<double> getPsi(BVertex Bv);
        double getSigma(BEdge Be);
        void computeSigma();

        // texture coordinates quantities
        std::vector<VertexData<std::complex<double>>> psi;
        std::vector<VertexData<double>> coords;           // these are the direct coords 
        std::vector<EdgeData<double>> sigma;
        FaceData<std::vector<Vector2>> texCoords;    
        FaceData<int> zeros;

        EdgeData<double> errors; 
};