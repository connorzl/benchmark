#pragma once

#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/geometry.h"
#include "geometrycentral/BHalfedgemesh.h"
#include "geometrycentral/linear_solvers.h"

using namespace geometrycentral;

class HodgeDecomposition {
    public:
        HodgeDecomposition(std::vector<VertexData<size_t>> &vertexIndices_, 
                           std::vector<EdgeData<size_t>> &edgeIndices_, 
                           std::vector<FaceData<size_t>> &faceIndices_, 
                           BranchCoverTopology &BC_,
                           HalfedgeData<double> &cmAngles_, FaceData<double> &cmAreas_,
                           HalfedgeMesh* mesh_, int numSingularities_);

        Eigen::MatrixXd computeHarmonicComponent(Eigen::MatrixXd omega);
        Eigen::MatrixXd computeExactComponent(Eigen::MatrixXd omega);
        Eigen::MatrixXd computeCoExactComponent(Eigen::MatrixXd omega);

    private:
        std::vector<VertexData<size_t>> vertexIndices; 
        std::vector<EdgeData<size_t>> edgeIndices;
        std::vector<FaceData<size_t>> faceIndices; 
        BranchCoverTopology BC;
        HalfedgeData<double> cmAngles;
        FaceData<double> cmAreas;
        HalfedgeMesh* mesh;
        int numSingularities;

        Eigen::SparseMatrix<double> buildHodgeStar1Form();
        Eigen::SparseMatrix<double> buildExteriorDerivative0Form();
        Eigen::SparseMatrix<double> buildExteriorDerivative1Form();
        Eigen::SparseMatrix<double> buildLaplacian();
        Eigen::SparseMatrix<double> invertDiagonal(Eigen::SparseMatrix<double> M);

        Eigen::SparseMatrix<double> buildExteriorDerivative0FormInterior(size_t numInterior);

        Eigen::SparseMatrix<double> hodge1;
        Eigen::SparseMatrix<double> hodge1Inv;
        Eigen::SparseMatrix<double> d0;
        Eigen::SparseMatrix<double> d0T;
        Eigen::SparseMatrix<double> d1;
        Eigen::SparseMatrix<double> d1T;

        Eigen::SparseMatrix<double> A;
        Eigen::SparseMatrix<double> B;
        std::shared_ptr<PositiveDefiniteSolver<double>> PDS;
        std::shared_ptr<SquareSolver<double>> SS;

        // for computing exact component
        std::vector<VertexData<size_t>> interiorVertexIndices;
        size_t iInterior;
        Eigen::SparseMatrix<double> d0I;
};