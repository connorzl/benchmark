#pragma once

#include "polyscope/polyscope.h"
#include "polyscope/hodge_decomposition.h"

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/operators.h"
#include "geometrycentral/BHalfedgemesh.h"
#include "geometrycentral/linear_solvers.h"

#include "geometrycentral/god_field.h"
#include "geometrycentral/uniformization.h"

using namespace geometrycentral;

class QuadMesh {
    public:
        QuadMesh() {};
        QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g, double rot=0, double scale=100, bool visualize=true);
        void generateQuadMesh();
        void updateQuadMesh(double rot, double scale);
        
    private:
        // ------------------ HELPER FUNCTIONS ------------------ // 
     
        // branch cover related
        //void computeBranchCover(bool improve=true);
        void computeCrossFieldCMBranchCover();

        // stripes
        void computeStripes();
        
        // optimization and texture coordinates
        //void optimizeSimpleLocally();
        //void optimizeSimpleGlobally();
        void optimizeHarmonic(bool visualize = false);
        bool textureCoordinates();

        // visualization
        void visualize();

        // ------------------ OTHER QUANTITIES ------------------ // 

        HalfedgeMesh* mesh;
        Geometry<Euclidean>* geom;
        bool vis;
        double rot;
        double scale;

        // globally optimal direction field
        GodField GF;
        VertexData<int> singularities;
        size_t numSingularities;

        // uniformization
        Uniformization U;
        
        // branch cover quantities
        BranchCoverTopology BC;
        //HalfedgeData<int> eta;
        std::vector<FaceData<std::complex<double>>> branchCoverFields; // cross field lifted to branch cover
        EdgeData<double> errors; 
        //std::vector<FaceData<std::complex<double>>> xBasis;

        // stripes helpers
        void computeOmega();
        void fixOmegaBranchPoints();
        Eigen::SparseMatrix<double> energyMatrix();
        //Eigen::SparseMatrix<double> energyMatrix2();
        Eigen::SparseMatrix<double> massMatrix();
        Eigen::SparseMatrix<double> linkMatrix();

        // stripes quantities
        std::vector<EdgeData<double>> omega;
        std::vector<VertexData<size_t>> BVertexIndices;

        // texture coordinates helpers
        std::complex<double> getPsi(BVertex Bv);
        double getSigma(BEdge Be);
        void computeSigma();
        //void computeSigma2();
        //Eigen::MatrixXd buildLocalEnergy(BFace Bf);
        //void computeIdealPsi(BVertex Bv_i, size_t i_re, size_t j_re, size_t k_re, Eigen::MatrixXd A, 
        //                     std::complex<double> &psi_j, std::complex<double> &psi_k);

        // texture coordinates quantities
        std::vector<VertexData<std::complex<double>>> psi;
        std::vector<VertexData<double>> coords;           // these are the direct coords 
        std::vector<EdgeData<double>> sigma;
        FaceData<std::vector<Vector2>> texCoords;    
        FaceData<int> zeros;
};