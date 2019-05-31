#include "polyscope/hodge_decomposition.h"

HodgeDecomposition::HodgeDecomposition(std::vector<VertexData<size_t>> &v, 
    std::vector<EdgeData<size_t>> &e, std::vector<FaceData<size_t>> &f, 
    BranchCoverTopology &BC_, HalfedgeData<double> &cmCotans_, 
    FaceData<double> &cmAreas_, HalfedgeMesh* mesh_, int numSingularities_) : 
        vertexIndices(v), edgeIndices(e), faceIndices(f), 
        BC(BC_), cmCotans(cmCotans_), cmAreas(cmAreas_),
        mesh(mesh_), numSingularities(numSingularities_) {

            hodge1 = buildHodgeStar1Form();
            hodge1Inv = invertDiagonal(hodge1);
            d0 = buildExteriorDerivative0Form();
            d0T = d0.transpose();
            d1 = buildExteriorDerivative1Form();
            d1T = d1.transpose();

            A = buildLaplacian();
            B = d1 * hodge1Inv * d1T;

            // to support meshes with boundary, make sure A only contains interior vertices
            for (size_t i = 0; i < 4; i++) {
                VertexData<size_t> indexSheet(mesh);
                interiorVertexIndices.push_back(indexSheet);
            }

            Eigen::Array<bool,Eigen::Dynamic,1> interior(4 * mesh->nVertices(),1);
            std::vector<BVertex> allBVertices = BC.allVertices();
            iInterior = 0;
            for (BVertex Bv : allBVertices) {
                if (!Bv.v.isBoundary()) {
                    interiorVertexIndices[Bv.sheet][Bv.v] = iInterior++;
                } 
                interior(vertexIndices[Bv.sheet][Bv.v],0) = (!Bv.v.isBoundary());
            }
            BlockDecompositionResult<double> blocks = blockDecomposeSquare(A, interior);
            A = blocks.AA;

            d0I = buildExteriorDerivative0FormInterior(iInterior);

            PDS = std::make_shared<PositiveDefiniteSolver<double>>(A);
            SS = std::make_shared<SquareSolver<double>>(B);
}

Eigen::SparseMatrix<double> HodgeDecomposition::buildHodgeStar1Form() {
    size_t nBEdges = 4 * mesh->nEdges();
    Eigen::SparseMatrix<double> hodge1(nBEdges, nBEdges);
    std::vector<Eigen::Triplet<double>> triplets;

    // loop through all BEdges
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BVertex Bv_i = Be.halfedge().vertex();
        BVertex Bv_j = Be.halfedge().next().vertex();

        size_t i = edgeIndices[Be.sheet][Be.e];
        BHalfedge BHe = Be.halfedge();
        double cotA = (BHe.he.isReal()) ? cmCotans[BHe.he] : 0;
        double cotB = (BHe.he.twin().isReal()) ? cmCotans[BHe.he.twin()] : 0; 
        double w = (cotA + cotB) / 2.0;
        if (w < 1e-8) w = 1e-8;
        //std::cout << w << std::endl;

        triplets.push_back( Eigen::Triplet<double>(i, i, w) );
    }

    hodge1.setFromTriplets(triplets.begin(), triplets.end());
    return hodge1;
}

Eigen::SparseMatrix<double> HodgeDecomposition::buildExteriorDerivative0Form() {
    size_t nBEdges = 4 * mesh->nEdges();
    size_t nNonSingular = mesh->nVertices() - numSingularities;
    size_t nBVertices = 4 * nNonSingular + numSingularities;
    Eigen::SparseMatrix<double> d0(nBEdges, nBVertices);
    std::vector<Eigen::Triplet<double>> triplets;

    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BVertex Bv_i = Be.halfedge().vertex();
        BVertex Bv_j = Be.halfedge().next().vertex();

        size_t k = edgeIndices[Be.sheet][Be.e];
        size_t i = vertexIndices[Bv_i.sheet][Bv_i.v];
        size_t j = vertexIndices[Bv_j.sheet][Bv_j.v];
        
        triplets.push_back( Eigen::Triplet<double>(k, i, 1) );
        triplets.push_back( Eigen::Triplet<double>(k, j, -1) );
    }

    d0.setFromTriplets(triplets.begin(), triplets.end());
    return d0;
}

Eigen::SparseMatrix<double> HodgeDecomposition::buildExteriorDerivative0FormInterior(size_t numInterior) {
    size_t nBEdges = 4 * mesh->nEdges();
    size_t nNonSingular = mesh->nVertices() - numSingularities;
    size_t nBVertices = 4 * nNonSingular + numSingularities;
    Eigen::SparseMatrix<double> d0(nBEdges, numInterior);
    std::vector<Eigen::Triplet<double>> triplets;

    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BVertex Bv_i = Be.halfedge().vertex();
        BVertex Bv_j = Be.halfedge().next().vertex();

        size_t k = edgeIndices[Be.sheet][Be.e];

        if (!Bv_i.v.isBoundary()) {
            size_t i = interiorVertexIndices[Bv_i.sheet][Bv_i.v];
            triplets.push_back( Eigen::Triplet<double>(k, i, 1) );
        }

        if (!Bv_j.v.isBoundary()) {
            size_t j = interiorVertexIndices[Bv_j.sheet][Bv_j.v];
            triplets.push_back( Eigen::Triplet<double>(k, j, -1) );
        }
    }

    d0.setFromTriplets(triplets.begin(), triplets.end());
    return d0;
}

Eigen::SparseMatrix<double> HodgeDecomposition::buildExteriorDerivative1Form() { 
    size_t nBFaces = 4 * mesh->nFaces();
    size_t nBEdges = 4 * mesh->nEdges();
    Eigen::SparseMatrix<double> d1(nBFaces, nBEdges);
    std::vector<Eigen::Triplet<double>> triplets;

    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        size_t i = faceIndices[Bf.sheet][Bf.f];
        BHalfedge BHe = Bf.halfedge();
        do {
            BEdge Be = BHe.edge();
            size_t j = edgeIndices[Be.sheet][Be.e];
            double sign = (Be.halfedge().he == BHe.he) ? 1 : -1;
            triplets.push_back( Eigen::Triplet<double>(i,j,sign) );
            BHe = BHe.next();
        } while (BHe != Bf.halfedge());
    }

    d1.setFromTriplets(triplets.begin(),triplets.end());
    return d1;
}

Eigen::SparseMatrix<double> HodgeDecomposition::buildLaplacian() {
    Eigen::SparseMatrix<double> A = d0T * hodge1 * d0;    
    /*
    // build this the classic way, gets equivalent results with the above 
    size_t nNonSingular = mesh->nVertices() - numSingularities;
    size_t nBVertices = 4 * nNonSingular + numSingularities;
    Eigen::SparseMatrix<double> A(nBVertices, nBVertices);
    std::vector<Eigen::Triplet<double>> tripletsCot;

    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        BHalfedge BHe_ij = Bf.halfedge();
        BHalfedge BHe_jk = BHe_ij.next();
        BHalfedge BHe_ki = BHe_jk.next();

        BVertex Bv_i = BHe_ij.vertex();
        BVertex Bv_j = BHe_jk.vertex();
        BVertex Bv_k = BHe_ki.vertex();

        double cot_ij, cot_jk, cot_ki;

        // compute cotan
        cot_ij = cmCotans[BHe_ij.he];
        cot_jk = cmCotans[BHe_jk.he];
        cot_ki = cmCotans[BHe_ki.he];

        size_t i = vertexIndices[Bv_i.sheet][Bv_i.v];
        size_t j = vertexIndices[Bv_j.sheet][Bv_j.v];
        size_t k = vertexIndices[Bv_k.sheet][Bv_k.v];

        // add weights to entries
        // ii, ij, ik
        tripletsCot.push_back( Eigen::Triplet<double>(i,i, (cot_ij + cot_ki) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(i,j, -cot_ij / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(i,k, -cot_ki / 2.0) );

        // jj, ji, jk
        tripletsCot.push_back( Eigen::Triplet<double>(j,j, (cot_jk + cot_ij) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j,i, -cot_ij / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j,k, -cot_jk / 2.0) );

        // kk, ki, kj
        tripletsCot.push_back( Eigen::Triplet<double>(k,k, (cot_jk + cot_ki) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k,i, -cot_ki / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k,j, -cot_jk / 2.0) );
    }
    A.setFromTriplets(tripletsCot.begin(),tripletsCot.end());
    */

    Eigen::SparseMatrix<double> offset(A.rows(),A.cols());
    offset.setIdentity();
    offset = offset * 1e-8;
    return A + offset;
}

Eigen::SparseMatrix<double> HodgeDecomposition::invertDiagonal(Eigen::SparseMatrix<double> M) {
    Eigen::MatrixXd ones = Eigen::MatrixXd::Ones(M.rows(),1);
    if (M.rows() != M.cols()) throw std::logic_error ("not square matrix");
    Eigen::MatrixXd mult = M * ones;

    Eigen::SparseMatrix<double> result(M.rows(), M.cols());
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t i = 0; i < (size_t)M.rows(); i++) {
        if (mult(i,0) == 0) throw std::logic_error ("division by 0 rip");
        triplets.push_back( Eigen::Triplet<double>( i, i, 1.0 / mult(i,0)) );
    }

    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

// split the laplacian to be ii on interior vertices for meshes with boundary
Eigen::MatrixXd HodgeDecomposition::computeExactComponent(Eigen::MatrixXd omega) {
    Eigen::MatrixXd b = d0T * hodge1 * omega;

    Eigen::MatrixXd b_interior(iInterior,1);
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (!Bv.v.isBoundary()) {
            size_t i = interiorVertexIndices[Bv.sheet][Bv.v];
            size_t j = vertexIndices[Bv.sheet][Bv.v];
            b_interior(i,0) = b(j,0);
        }
    }
    Eigen::MatrixXd alpha = PDS->solve(b_interior); 
    Eigen::MatrixXd dAlpha = d0I * alpha;
    return dAlpha;
}

// TODO: have d1 take into account imaginary faces, and then split into II for meshes with boundary
Eigen::MatrixXd HodgeDecomposition::computeCoExactComponent(Eigen::MatrixXd omega) {
    Eigen::MatrixXd b = d1 * omega;

    Eigen::MatrixXd betaTilde = SS->solve(b);
    Eigen::MatrixXd dBeta = hodge1Inv * d1T * betaTilde;
    return dBeta;
}

Eigen::MatrixXd HodgeDecomposition::computeHarmonicComponent(Eigen::MatrixXd omega) {
    Eigen::MatrixXd dAlpha = computeExactComponent(omega);
    Eigen::MatrixXd dBeta = computeCoExactComponent(omega);

    std::cout << dAlpha.norm() << "," << dBeta.norm() << std::endl;
    // γ = ω - dα - 𝛿β
    return omega - dAlpha - dBeta;
}
