#include "polyscope/quad_mesh.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <queue>

QuadMesh::QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g, 
                   double rot_, double scale_, bool vis_) : mesh(m), geom(g), singularities(m), texCoords(m) {
        
    rot = rot_;
    scale = scale_;
    vis = vis_;

    for (int i = 0; i < 4; i++) {
        FaceData<std::complex<double>> sheetField(mesh,0);
        VertexData<size_t> sheetIndices(mesh);
        EdgeData<double> omegaSheet(mesh);
        branchCoverFields.push_back(sheetField);
        BVertexIndices.push_back(sheetIndices);
        omega.push_back(omegaSheet);
    }

    for (int i = 0; i < 2; i++) {
        VertexData<double> sheetCoords(mesh);
        VertexData<std::complex<double>> sheetPsi(mesh,0);
        EdgeData<double> sheetSigma(mesh);

        coords.push_back(sheetCoords);
        psi.push_back(sheetPsi);
        sigma.push_back(sheetSigma);
    }

    generateQuadMesh();
}

void QuadMesh::generateQuadMesh() {
    // first, compute cone points via globally optimal direction fields
    GF = GodField(mesh, geom);
    EdgeData<double> edgeLengths(mesh);
    geom->getEdgeLengths(edgeLengths);
    HalfedgeData<double> heAngles(mesh);
    geom->getHalfedgeAngles(heAngles);

    GF.computeCrossField(edgeLengths,heAngles);
    numSingularities = GF.computeSingularities(singularities);

    // next, uniformize the surface
    U = Uniformization(mesh, geom, singularities);
    U.uniformize();

    // update the cross field on the surface
    GF.computeCrossField(U.cmEdgeLengths, U.cmAngles);

    // build branch cover on uniformized surface
    BC = BranchCoverTopology(mesh, singularities);
    BC.validateEta(U.cmEdgeVector);
    BC.validateConnectivity();

    // compute a constant cross field on the surface
    computeCrossFieldCMBranchCover();
    // project to edges
    computeOmega();

    // minimize stripes energy + constraints along links
    computeStripes();

    // extract texture coordinates for visualization
    textureCoordinates();
    
    if (vis) {
        visualize();
    }
}

void QuadMesh::updateQuadMesh(double rot_, double scale_) {
    rot = rot_;
    scale = scale_;

    // compute a constant cross field on the surface
    computeCrossFieldCMBranchCover();
    // project to edges
    computeOmega();

    // minimize stripes energy + constraints along links
    computeStripes();

    // extract texture coordinates for visualization
    textureCoordinates();
    
    if (vis) {
        visualize();
    }
}

/* 
void QuadMesh::computeBranchCover(bool improve) {
    std::cout<< "Computing Branch Cover..." << std::endl;

    // OLD WAY OF COMPUTING IT BASED ON THE FIELD
    for (VertexPtr v : mesh->vertices()) {
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) {
                eta[he] = 0;
                continue;
            }
            std::complex<double> u_i = std::pow(GF.field[he.face()], 1.0 / 4);
            std::complex<double> u_j = std::pow(GF.field[he.twin().face()], 1.0 / 4);
            
            std::complex<double> theta_ij = U.cmEdgeVector[he];
            std::complex<double> theta_ji = U.cmEdgeVector[he.twin()];
            std::complex<double> r_ji = -theta_ij / theta_ji;
            r_ji = r_ji / std::abs(r_ji);
            double ang = std::arg(u_i / (r_ji * u_j));

            if (ang >= -M_PI_4 && ang < M_PI_4) {
                eta[he] = 0;
            } else if (ang >= M_PI_4 && ang < 3.0 * M_PI_4) {
                eta[he] = 1;
            } else if ((ang >= 3.0 * M_PI_4 && ang <= PI) || 
                       (ang < -3.0 * M_PI_4 && ang >= -PI)) {
                eta[he] = 2;
            } else {
                if (!(ang >= -3.0 * M_PI_4 && ang < -M_PI_4)) {
                    throw std::logic_error("angle not in right range for branch cover");
                }
                eta[he] = 3;
            }
        }
    }
    if (improve) {
        // perform BFS starting on some FacePtr
        FacePtr root = mesh->face(0);
        std::map<FacePtr,bool> visited;
        visited[root] = true;

        std::queue<FacePtr> Q;
        Q.push(root);
        size_t count = 1;
        
        while(!Q.empty()) {
            FacePtr currFace = Q.front();
            Q.pop();

            for (HalfedgePtr he : currFace.adjacentHalfedges()) {
                if (he.edge().isBoundary()) continue;
                FacePtr neighbor = he.twin().face();

                if (visited.find(neighbor) == visited.end()) {
                    Q.push(neighbor);
                    visited[neighbor] = true;
                    
                    int sheetDiff = eta[he.twin()];
                    int offset = (4 - sheetDiff) % 4;
                    for (HalfedgePtr Nhe : neighbor.adjacentHalfedges()) {
                        eta[Nhe] = (eta[Nhe] + offset + 4) % 4;
                        if (Nhe.twin().isReal()) {
                            eta[Nhe.twin()] = (eta[Nhe.twin()] - offset + 4) % 4;
                        }
                        if ( !Nhe.edge().isBoundary() && (eta[Nhe] + eta[Nhe.twin()]) % 4 != 0) {
                            throw std::logic_error("bad update to eta");
                        }
                    }
                    count++;
                }
            }
        }
        if (count != mesh->nFaces()) throw std::logic_error("BFS did not traverse all BFaces");
    }
    
    BC = BranchCoverTopology(mesh, eta, singularities);
    BC.validateConnectivity();

    // Sanity check the rotations here
    std::vector<BVertex> allVertices = BC.allVertices();
    for (BVertex Bv : allVertices) {
        if (Bv.v.isBoundary()) continue;
        std::complex<double> v0(((double) rand() / (RAND_MAX)), ((double) rand() / (RAND_MAX)));
        std::complex<double> vCurr = v0;

        // rotate vector around fan
        BHalfedge BHe = Bv.halfedge();
        BHalfedge firstHe = BHe;
        do {
            std::complex<double> theta_ij = U.cmEdgeVector[BHe.he];
            std::complex<double> theta_ji = U.cmEdgeVector[BHe.he.twin()];
            std::complex<double> r_ij = -theta_ji / theta_ij;
            r_ij = r_ij / std::abs(r_ij);
            vCurr = r_ij * vCurr;
            BHe = BHe.twin().next();
        } while (BHe.he != firstHe.he);

        // check that the resulting vector makes sense
        std::complex<double> i(0,1);
        std::complex<double> rotPos90 = std::exp(i * M_PI_2);
        std::complex<double> rotNeg90 = std::exp(i * -M_PI_2);
        if (singularities[Bv.v] == 0 && std::norm(vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around non-singular vertex");
        } else if (singularities[Bv.v] == 1 && std::norm(rotPos90 * vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around +1 singular vertex");
        } else if (singularities[Bv.v] == -1 && std::norm(rotNeg90 * vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around -1 singular vertex");
        }
    }

    std::cout << "Done!" << std::endl;
}
*/

void QuadMesh::computeCrossFieldCMBranchCover() {
    std::cout << "Computing Cross Field CM on Branch Cover..." << std::endl;
    
    // get faces for traversing branch cover
    std::vector<BFace> allBFaces = BC.allFaces();
    for (int i = 0; i < 4; i++) {
        FaceData<std::complex<double>> sheetField(mesh,0);
        branchCoverFields[i] = sheetField;

        //FaceData<std::complex<double>> sheetXBasis(mesh);
        //xBasis.push_back(sheetXBasis);
    }

    // perform BFS starting on some BFace
    BFace root = allBFaces[0];
    std::map<BFace,bool> visited;
    visited[root] = true;
    branchCoverFields[root.sheet][root.f] = std::exp(IM_I * rot);
    //xBasis[root.sheet][root.f] = std::complex<double> (1,0);

    std::queue<BFace> Q;
    Q.push(root);
    size_t count = 1;
    errors = EdgeData<double>(mesh, 0);

    while(!Q.empty()) {
        BFace currFace = Q.front();
        Q.pop();

        // get neighbors
        BHalfedge BHe = currFace.halfedge();
        do {
            if (BHe.he.edge().isBoundary()) {
                BHe = BHe.next();
                continue;
            }
            BFace neighbor = BHe.twin().face();
            if (visited.find(neighbor) == visited.end()) {
                Q.push(neighbor);
                visited[neighbor] = true;

                std::complex<double> theta_ij = U.cmEdgeVector[BHe.he];
                std::complex<double> theta_ji = U.cmEdgeVector[BHe.he.twin()];
                std::complex<double> r_ij = -theta_ji / theta_ij;
                r_ij = r_ij / std::abs(r_ij);

                branchCoverFields[neighbor.sheet][neighbor.f] = r_ij * branchCoverFields[currFace.sheet][currFace.f];
                //xBasis[neighbor.sheet][neighbor.f] = r_ij * xBasis[currFace.sheet][currFace.f];
                count++;
            } 
            BHe = BHe.next();
        } while (BHe != currFace.halfedge());
        
        // look at neighbor in the sheet above
        BFace neighbor = currFace;
        neighbor.sheet = (neighbor.sheet + 1) % 4;
        if (visited.find(neighbor) == visited.end()) {
            Q.push(neighbor);
            visited[neighbor] = true;
            branchCoverFields[neighbor.sheet][neighbor.f] = IM_I * branchCoverFields[currFace.sheet][currFace.f];
            //xBasis[neighbor.sheet][neighbor.f] = IM_I * xBasis[currFace.sheet][currFace.f];
            count++;
        }      
    }
    if (count != allBFaces.size()) throw std::logic_error("BFS did not traverse all BFaces");

    // compute errors here
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        if (Be.e.isBoundary()) continue;
        BHalfedge BHe = Be.halfedge();
        BFace Bf_ij = BHe.face();
        BFace Bf_ji = BHe.twin().face();

        std::complex<double> theta_ij = U.cmEdgeVector[BHe.he];
        std::complex<double> theta_ji = U.cmEdgeVector[BHe.he.twin()];
        std::complex<double> r_ji = -theta_ij / theta_ji;
        r_ji = r_ji / std::abs(r_ji);

        std::complex<double> f_ij = branchCoverFields[Bf_ij.sheet][Bf_ij.f];
        std::complex<double> f_ijt = r_ji * branchCoverFields[Bf_ji.sheet][Bf_ji.f];
        errors[Be.e] = std::max(errors[Be.e], std::abs(f_ijt - f_ij));
    }

    std::cout << "Done!" << std::endl;
}

void QuadMesh::computeOmega() {
    // compute omega on each edge of the branch cover
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BHalfedge he_ij = Be.halfedge();
        BHalfedge he_ji = he_ij.twin();

        std::complex<double> f_ij, f_ji, e_ij, e_ji;
        
        double boundaryScale = 1.0;
        if (he_ij.he.isReal()){
            f_ij = scale * branchCoverFields[he_ij.sheet][he_ij.face().f];
            e_ij = U.cmEdgeVector[he_ij.he];
        } else {
            f_ij = 0;
            e_ij = 0;
            boundaryScale = 2.0;
        }
        if (he_ji.he.isReal()) {
            f_ji = scale * branchCoverFields[he_ji.sheet][he_ji.face().f];
            e_ji = -U.cmEdgeVector[he_ji.he];
        } else {
            f_ji = 0;
            e_ji = 0;
            boundaryScale = 2.0;
        }

        // these 2 should equal each other, unless on boundary
        double prod_ij = dot(e_ij, f_ij);
        double prod_ji = dot(e_ji, f_ji);

        if (!Be.e.isBoundary() && std::abs(prod_ij - prod_ji) > 1e-6) {
            polyscope::warning("non-parallel vectors", std::to_string(std::abs(prod_ij - prod_ji)));
        }
        omega[Be.sheet][Be.e] = boundaryScale * 0.5 * (prod_ij + prod_ji);
    }

    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        BHalfedge Be_ij = Bf.halfedge();
        BHalfedge Be_jk = Bf.halfedge().next();
        BHalfedge Be_ki = Bf.halfedge().next().next();

        int c_ij = (Be_ij.edge().halfedge() == Be_ij) ? 1 : -1;
        int c_jk = (Be_jk.edge().halfedge() == Be_jk) ? 1 : -1;
        int c_ki = (Be_ki.edge().halfedge() == Be_ki) ? 1 : -1;

        double w_ij = c_ij * omega[Be_ij.edge().sheet][Be_ij.edge().e];
        double w_jk = c_jk * omega[Be_jk.edge().sheet][Be_jk.edge().e];
        double w_ki = c_ki * omega[Be_ki.edge().sheet][Be_ki.edge().e];

        double sum = w_ij + w_jk + w_ki;
        if (std::abs(sum) > 1e-6) {
            polyscope::warning("omega not locally integrable", std::to_string(sum));
        }
    }

    // boundary conditions for branch point holes
    fixOmegaBranchPoints();
}

void QuadMesh::fixOmegaBranchPoints() {
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] == 0) continue;

        BHalfedge BHe = Bv.halfedge();

        size_t numMaxima = 0;
        size_t numMinima = 0;

        size_t deg = Bv.v.degree();
        for (size_t i = 0; i < 2*deg; i++) {
            BHalfedge BHe_prev = BHe.twin().next().next();
            BHalfedge BHe_next = BHe.next();

            BEdge Be_prev = BHe_prev.edge();
            BEdge Be_next = BHe_next.edge();

            int sign_prev = (Be_prev.halfedge() == BHe_prev) ? 1 : -1;
            int sign_next = (Be_next.halfedge() == BHe_next) ? 1 : -1;

            double omega_prev = sign_prev * omega[Be_prev.sheet][Be_prev.e];
            double omega_next = sign_next * omega[Be_next.sheet][Be_next.e];

            if (omega_prev < 0 && omega_next > 0) {
                numMinima++;
            } else if (omega_prev > 0 && omega_next < 0) {
                numMaxima++;
            }
            BHe = BHe.next().next().twin();
        }

        if (singularities[Bv.v] == 1) {
            if ( !((numMaxima == 1 && numMinima == 2) || (numMaxima == 2 && numMinima == 1)) ) {
                std::cout << "numMax: " << numMaxima << ", numMin: " << numMinima << std::endl;
                throw std::logic_error("bad omega around a +1 singularity link");
            }
        } else if (singularities[Bv.v] == -1) {
            if ( !((numMaxima == 2 && numMinima == 3) || (numMaxima == 3 && numMinima == 2)) ) {
                throw std::logic_error("bad omega around a -1 singularity link");
            }
        } else {
            throw std::logic_error("omega check only implemented for +1/-1 singularities");
        }
 
        //std::cout << "singular index: " << singularities[Bv.v] << " num max: " << numMaxima << " num min: " << numMinima << std::endl;
    }
}

#if 0
Eigen::SparseMatrix<double> QuadMesh::energyMatrix2() {
    size_t iNonSingular = 0;
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] != 0) {
            continue;
        } else if (Bv.sheet == 0 || Bv.sheet == 1) {
            BVertexIndices[Bv.sheet][Bv.v] = iNonSingular++;
        }
    }

    // 1 copy of singular vertices (sheet 0) and 2 copies of nonsingular vertices (sheets 0 and 1)
    // multiply by 2 to get real DOFs rather than complex    
    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular); // + numSingularities);
    Eigen::SparseMatrix<double> A(numPsi, numPsi);
    std::vector<Eigen::Triplet<double>> triplets;
    
    Eigen::SparseMatrix<double> Acot(numPsi,numPsi);
    Eigen::SparseMatrix<double> AzArea(numPsi,numPsi);
    Eigen::SparseMatrix<double> Adrift(numPsi,numPsi);
    std::vector<Eigen::Triplet<double>> tripletsCot;
    std::vector<Eigen::Triplet<double>> tripletszArea;
    std::vector<Eigen::Triplet<double>> tripletsDrift;
    
    // For vertex on some sheet, returns indices to index real and imaginary parts of psi,
    // as well as sign for imaginary componenet to implement conjugate
    auto getIndexData = [&](BVertex Bv, size_t& realInd, size_t& imagInd, double& imagSign) {
        bool v_singular = (singularities[Bv.v] != 0);
        if (v_singular) {  
            if (Bv.sheet != BC.singularSheet) throw std::logic_error("wrong sheet for singular vertex");
            realInd = 2 * BVertexIndices[BC.singularSheet][Bv.v];    // all singularities live on a single sheet (0)
            imagSign = 1;
        } else if ((!v_singular && Bv.sheet == 0) || (!v_singular && Bv.sheet == 1)) {
            realInd = 2 * BVertexIndices[Bv.sheet][Bv.v];            // sheet 0 and 1 treated as normal
            imagSign = 1;
        } else if (!v_singular && Bv.sheet == 2) {
            realInd = 2 * BVertexIndices[0][Bv.v];                   // sheet 2 shares indices with sheet 0
            imagSign = -1;
        } else if (!v_singular && Bv.sheet == 3) {
            realInd = 2 * BVertexIndices[1][Bv.v];                   // sheet 3 shares indices with sheet 1
            imagSign = -1;
        }
        imagInd = realInd+1;
    };

    double scale = 100;
    double adjust = 1.0;
    std::vector<BFace> allBFaces = BC.allFaces();

    /* 
    for (BFace Bf : allBFaces) {
        BHalfedge BHe_ij = Bf.halfedge();
        BHalfedge BHe_jk = BHe_ij.next();
        BHalfedge BHe_ki = BHe_jk.next();

        BVertex Bv_i = BHe_ij.vertex();
        BVertex Bv_j = BHe_jk.vertex();
        BVertex Bv_k = BHe_ki.vertex();

        // skip over faces that contain a singular vertex
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || singularities[Bv_k.v] != 0) {
            continue;
        }

        // compute cotan, zArea, and drift
        double cot_ij = 1.0 / tan(U.cmAngles[BHe_ij.he]);
        double cot_jk = 1.0 / tan(U.cmAngles[BHe_jk.he]);
        double cot_ki = 1.0 / tan(U.cmAngles[BHe_ki.he]);
        double zArea = adjust * adjust * scale * scale * U.cmAreas[Bf.f];

        std::complex<double> z = adjust * scale * branchCoverFields[Bf.sheet][Bf.f];
        std::complex<double> e_i_perp = IM_I * thetaCM[BHe_jk.he];
        std::complex<double> e_j_perp = IM_I * thetaCM[BHe_ki.he];
        std::complex<double> e_k_perp = IM_I * thetaCM[BHe_ij.he];
       
        double drift_ij = (1.0 / 6.0) * dot(e_i_perp - e_j_perp, z); 
        double drift_jk = (1.0 / 6.0) * dot(e_j_perp - e_k_perp, z); 
        double drift_ki = (1.0 / 6.0) * dot(e_k_perp - e_i_perp, z); 
        double drift_ji = -drift_ij;
        double drift_kj = -drift_jk;
        double drift_ik = -drift_ki;

        // get index and sign information
        size_t i_re, i_im, j_re, j_im, k_re, k_im;
        double si_im, sj_im, sk_im;
        getIndexData(Bv_i, i_re, i_im, si_im);
        getIndexData(Bv_j, j_re, j_im, sj_im);
        getIndexData(Bv_k, k_re, k_im, sk_im);

        double A_ii_00, A_ii_01, A_ii_10, A_ii_11;
        double A_ij_00, A_ij_01, A_ij_10, A_ij_11;
        double A_ik_00, A_ik_01, A_ik_10, A_ik_11;
        double A_ji_00, A_ji_01, A_ji_10, A_ji_11;
        double A_jj_00, A_jj_01, A_jj_10, A_jj_11;
        double A_jk_00, A_jk_01, A_jk_10, A_jk_11;
        double A_ki_00, A_ki_01, A_ki_10, A_ki_11;
        double A_kj_00, A_kj_01, A_kj_10, A_kj_11;
        double A_kk_00, A_kk_01, A_kk_10, A_kk_11;

        A_ii_00 = (cot_ij + cot_ki) / 2.0 + zArea / 6.0;  A_ii_01 = 0;
        A_ii_10 = 0;                                      A_ii_11 = (cot_ij + cot_ki) / 2.0 + zArea / 6.0;

        A_ij_00 = -cot_ij / 2.0 + zArea / 12.0;           A_ij_01 = sj_im * drift_ij;
        A_ij_10 =  si_im * -drift_ij;                     A_ij_11 = si_im * sj_im * (-cot_ij / 2.0 + zArea / 12.0);   

        A_ik_00 = -cot_ki / 2.0 + zArea / 12.0;           A_ik_01 = sk_im * drift_ik;
        A_ik_10 = si_im * -drift_ik;                      A_ik_11 = si_im * sk_im * (-cot_ki / 2.0 + zArea / 12.0);

        //-----------------------------------------------------------------------------------------------------------//
        
        A_ji_00 = -cot_ij / 2.0 + zArea / 12.0;           A_ji_01 = si_im * drift_ji;
        A_ji_10 = sj_im * -drift_ji;                      A_ji_11 = sj_im * si_im * (-cot_ij / 2.0 + zArea / 12.0);

        A_jj_00 = (cot_jk + cot_ij) / 2.0 + zArea / 6.0;  A_jj_01 = 0;
        A_jj_10 = 0;                                      A_jj_11 = (cot_jk + cot_ij) / 2.0 + zArea / 6.0;

        A_jk_00 = -cot_jk / 2.0 + zArea / 12.0;           A_jk_01 = sk_im * drift_jk;
        A_jk_10 = sj_im * -drift_jk;                      A_jk_11 = sj_im * sk_im * (-cot_jk / 2.0 + zArea / 12.0);

        //-----------------------------------------------------------------------------------------------------------//
        
        A_ki_00 = -cot_ki / 2.0 + zArea / 12.0;           A_ki_01 = si_im * drift_ki;
        A_ki_10 = sk_im * -drift_ki;                      A_ki_11 = sk_im * si_im * (-cot_ki / 2.0 + zArea / 12.0);

        A_kj_00 = -cot_jk / 2.0 + zArea / 12.0;           A_kj_01 = sj_im * drift_kj;
        A_kj_10 = sk_im * -drift_kj;                      A_kj_11 = sk_im * sj_im * (-cot_jk / 2.0 + zArea / 12.0);

        A_kk_00 = (cot_jk + cot_ki) / 2.0 + zArea / 6.0;  A_kk_01 = 0;
        A_kk_10 = 0;                                      A_kk_11 = (cot_jk + cot_ki) / 2.0 + zArea / 6.0;

        triplets.push_back( Eigen::Triplet<double>(i_re,i_re, A_ii_00) ); triplets.push_back( Eigen::Triplet<double>(i_re,i_im, A_ii_01) );
        triplets.push_back( Eigen::Triplet<double>(i_im,i_re, A_ii_10) ); triplets.push_back( Eigen::Triplet<double>(i_im,i_im, A_ii_11) );
        triplets.push_back( Eigen::Triplet<double>(i_re,j_re, A_ij_00) ); triplets.push_back( Eigen::Triplet<double>(i_re,j_im, A_ij_01) );
        triplets.push_back( Eigen::Triplet<double>(i_im,j_re, A_ij_10) ); triplets.push_back( Eigen::Triplet<double>(i_im,j_im, A_ij_11) );
        triplets.push_back( Eigen::Triplet<double>(i_re,k_re, A_ik_00) ); triplets.push_back( Eigen::Triplet<double>(i_re,k_im, A_ik_01) );
        triplets.push_back( Eigen::Triplet<double>(i_im,k_re, A_ik_10) ); triplets.push_back( Eigen::Triplet<double>(i_im,k_im, A_ik_11) );

        triplets.push_back( Eigen::Triplet<double>(j_re,i_re, A_ji_00) ); triplets.push_back( Eigen::Triplet<double>(j_re,i_im, A_ji_01) );
        triplets.push_back( Eigen::Triplet<double>(j_im,i_re, A_ji_10) ); triplets.push_back( Eigen::Triplet<double>(j_im,i_im, A_ji_11) );
        triplets.push_back( Eigen::Triplet<double>(j_re,j_re, A_jj_00) ); triplets.push_back( Eigen::Triplet<double>(j_re,j_im, A_jj_01) );
        triplets.push_back( Eigen::Triplet<double>(j_im,j_re, A_jj_10) ); triplets.push_back( Eigen::Triplet<double>(j_im,j_im, A_jj_11) );
        triplets.push_back( Eigen::Triplet<double>(j_re,k_re, A_jk_00) ); triplets.push_back( Eigen::Triplet<double>(j_re,k_im, A_jk_01) );
        triplets.push_back( Eigen::Triplet<double>(j_im,k_re, A_jk_10) ); triplets.push_back( Eigen::Triplet<double>(j_im,k_im, A_jk_11) );

        triplets.push_back( Eigen::Triplet<double>(k_re,i_re, A_ki_00) ); triplets.push_back( Eigen::Triplet<double>(k_re,i_im, A_ki_01) );
        triplets.push_back( Eigen::Triplet<double>(k_im,i_re, A_ki_10) ); triplets.push_back( Eigen::Triplet<double>(k_im,i_im, A_ki_11) );
        triplets.push_back( Eigen::Triplet<double>(k_re,j_re, A_kj_00) ); triplets.push_back( Eigen::Triplet<double>(k_re,j_im, A_kj_01) );
        triplets.push_back( Eigen::Triplet<double>(k_im,j_re, A_kj_10) ); triplets.push_back( Eigen::Triplet<double>(k_im,j_im, A_kj_11) );
        triplets.push_back( Eigen::Triplet<double>(k_re,k_re, A_kk_00) ); triplets.push_back( Eigen::Triplet<double>(k_re,k_im, A_kk_01) );
        triplets.push_back( Eigen::Triplet<double>(k_im,k_re, A_kk_10) ); triplets.push_back( Eigen::Triplet<double>(k_im,k_im, A_kk_11) );
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    */

    // build cotan component
    for (BFace Bf : allBFaces) {
        BHalfedge BHe_ij = Bf.halfedge();
        BHalfedge BHe_jk = BHe_ij.next();
        BHalfedge BHe_ki = BHe_jk.next();

        BVertex Bv_i = BHe_ij.vertex();
        BVertex Bv_j = BHe_jk.vertex();
        BVertex Bv_k = BHe_ki.vertex();

        // skip over faces that contain a singular vertex
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || singularities[Bv_k.v] != 0) {
            continue;
        }

        // compute cotan
        double cot_ij = 1.0 / tan(U.cmAngles[BHe_ij.he]);
        double cot_jk = 1.0 / tan(U.cmAngles[BHe_jk.he]);
        double cot_ki = 1.0 / tan(U.cmAngles[BHe_ki.he]);

        // get index and sign information
        size_t i_re, i_im, j_re, j_im, k_re, k_im;
        double si_im, sj_im, sk_im;
        getIndexData(Bv_i, i_re, i_im, si_im);
        getIndexData(Bv_j, j_re, j_im, sj_im);
        getIndexData(Bv_k, k_re, k_im, sk_im);

        // add weights to entries
        // ii, ij, ik
        tripletsCot.push_back( Eigen::Triplet<double>(i_re,i_re, (cot_ij + cot_ki) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(i_im,i_im, (cot_ij + cot_ki) / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(i_re,j_re, -cot_ij / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(i_im,j_im, si_im * sj_im * -cot_ij / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(i_re,k_re, -cot_ki / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(i_im,k_im, si_im * sk_im * -cot_ki / 2.0) );

        // ji, jj, jk
        tripletsCot.push_back( Eigen::Triplet<double>(j_re,i_re, -cot_ij / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j_im,i_im, sj_im * si_im * -cot_ij / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(j_re,j_re, (cot_jk + cot_ij) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j_im,j_im, (cot_jk + cot_ij) / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(j_re,k_re, -cot_jk / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j_im,k_im, sj_im * sk_im * -cot_jk / 2.0) );

        // ki, kj, kk
        tripletsCot.push_back( Eigen::Triplet<double>(k_re,i_re, -cot_ki / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k_im,i_im, sk_im * si_im * -cot_ki / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(k_re,j_re, -cot_jk / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k_im,j_im, sk_im * sj_im * -cot_jk / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(k_re,k_re, (cot_jk + cot_ki) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k_im,k_im, (cot_jk + cot_ki) / 2.0) );
    }
    Acot.setFromTriplets(tripletsCot.begin(),tripletsCot.end());

    // build zArea component
    for (BFace Bf : allBFaces) {
        BHalfedge BHe_ij = Bf.halfedge();
        BHalfedge BHe_jk = BHe_ij.next();
        BHalfedge BHe_ki = BHe_jk.next();

        BVertex Bv_i = BHe_ij.vertex();
        BVertex Bv_j = BHe_jk.vertex();
        BVertex Bv_k = BHe_ki.vertex();

        // skip over faces that contain a singular vertex
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || singularities[Bv_k.v] != 0) {
            continue;
        }

        // compute zArea weight
        double zArea = adjust * adjust * scale * scale * U.cmAreas[Bf.f];

        // get index and sign information
        size_t i_re, i_im, j_re, j_im, k_re, k_im;
        double si_im, sj_im, sk_im;
        getIndexData(Bv_i, i_re, i_im, si_im);
        getIndexData(Bv_j, j_re, j_im, sj_im);
        getIndexData(Bv_k, k_re, k_im, sk_im);

        // add weights to entries
        // ii, ij, ik
        tripletszArea.push_back( Eigen::Triplet<double>(i_re,i_re, zArea / 6.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(i_im,i_im, zArea / 6.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(i_re,j_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(i_im,j_im, si_im * sj_im * zArea / 12.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(i_re,k_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(i_im,k_im, si_im * sk_im * zArea / 12.0) );

        // ji, jj, jk
        tripletszArea.push_back( Eigen::Triplet<double>(j_re,i_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(j_im,i_im, sj_im * si_im * zArea / 12.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(j_re,j_re, zArea / 6.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(j_im,j_im, zArea / 6.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(j_re,k_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(j_im,k_im, sj_im * sk_im * zArea / 12.0) );

        // ki, kj, kk
        tripletszArea.push_back( Eigen::Triplet<double>(k_re,i_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(k_im,i_im, sk_im * si_im * zArea / 12.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(k_re,j_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(k_im,j_im, sk_im * sj_im * zArea / 12.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(k_re,k_re, zArea / 6.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(k_im,k_im, zArea / 6.0) );
    }
    AzArea.setFromTriplets(tripletszArea.begin(),tripletszArea.end());

    // build drift component
    for (BFace Bf : allBFaces) {
        BHalfedge BHe_ij = Bf.halfedge();
        BHalfedge BHe_jk = BHe_ij.next();
        BHalfedge BHe_ki = BHe_jk.next();

        BVertex Bv_i = BHe_ij.vertex();
        BVertex Bv_j = BHe_jk.vertex();
        BVertex Bv_k = BHe_ki.vertex();

        // skip over faces that contain a singular vertex
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || singularities[Bv_k.v] != 0) {
            continue;
        }

        // compute drift weight
        std::complex<double> z = adjust * scale * branchCoverFields[Bf.sheet][Bf.f];
        std::complex<double> e_i_perp = IM_I * thetaCM[BHe_jk.he];
        std::complex<double> e_j_perp = IM_I * thetaCM[BHe_ki.he];
        std::complex<double> e_k_perp = IM_I * thetaCM[BHe_ij.he];
       
        double drift_ij = (1.0 / 6.0) * ( dot(e_i_perp, z) - dot(e_j_perp, z) ); 
        double drift_jk = (1.0 / 6.0) * ( dot(e_j_perp, z) - dot(e_k_perp, z) ); 
        double drift_ki = (1.0 / 6.0) * ( dot(e_k_perp, z) - dot(e_i_perp, z) ); 

        double drift_ji = -drift_ij;
        double drift_kj = -drift_jk;
        double drift_ik = -drift_ki;

        // get index and sign information
        size_t i_re, i_im, j_re, j_im, k_re, k_im;
        double si_im, sj_im, sk_im;
        getIndexData(Bv_i, i_re, i_im, si_im);
        getIndexData(Bv_j, j_re, j_im, sj_im);
        getIndexData(Bv_k, k_re, k_im, sk_im);

        // add weights to entries
        // ii, ij, ik
        tripletsDrift.push_back( Eigen::Triplet<double>(i_re,j_im, sj_im * drift_ij) );
        tripletsDrift.push_back( Eigen::Triplet<double>(i_im,j_re, si_im * -drift_ij) );

        tripletsDrift.push_back( Eigen::Triplet<double>(i_re,k_im, sk_im * drift_ik) );
        tripletsDrift.push_back( Eigen::Triplet<double>(i_im,k_re, si_im * -drift_ik) );

        // jj, ji, jk
        tripletsDrift.push_back( Eigen::Triplet<double>(j_re,i_im, si_im * drift_ji) );
        tripletsDrift.push_back( Eigen::Triplet<double>(j_im,i_re, sj_im * -drift_ji) );

        tripletsDrift.push_back( Eigen::Triplet<double>(j_re,k_im, sk_im * drift_jk) );
        tripletsDrift.push_back( Eigen::Triplet<double>(j_im,k_re, sj_im * -drift_jk) );

        // kk, ki, kj
        tripletsDrift.push_back( Eigen::Triplet<double>(k_re,i_im, si_im * drift_ki) );
        tripletsDrift.push_back( Eigen::Triplet<double>(k_im,i_re, sk_im * -drift_ki) );

        tripletsDrift.push_back( Eigen::Triplet<double>(k_re,j_im, sj_im * drift_kj) );
        tripletsDrift.push_back( Eigen::Triplet<double>(k_im,j_re, sk_im * -drift_kj) );
    }
    Adrift.setFromTriplets(tripletsDrift.begin(), tripletsDrift.end());
    
    checkHermitian(Acot);
    std::cout << "Acot symmetric" << std::endl;
    checkHermitian(AzArea);
    std::cout << "AzArea symmetric" << std::endl;
    checkHermitian(Adrift);
    std::cout << "Adrift symmetric" << std::endl;
    return Acot + AzArea + Adrift;
    
   //checkHermitian(A);
   //return A;
}
#endif

Eigen::SparseMatrix<double> QuadMesh::energyMatrix() {
    // 1 copy of singular vertices (sheet 0) and 2 copies of nonsingular vertices (sheets 0 and 1)
    // multiply by 2 to get real DOFs rather than complex    
    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular); // + numSingularities);
    Eigen::SparseMatrix<double> A(numPsi,numPsi);
    std::vector<Eigen::Triplet<double>> triplets;
    
    std::vector<BHalfedge> allBHalfedges = BC.allHalfedges();
    for (BHalfedge BHe : allBHalfedges) {
        // skip over halfedges that contain a singular vertex
        if (singularities[BHe.he.vertex()] != 0 || singularities[BHe.next().he.vertex()] != 0) continue;

        // compute cotan weights
        double cotA, cotB, cotWeight;
        VertexPtr A,B;

        // first check if boundary halfedges
        if (BHe.he.isReal()) {
            cotA = 1.0 / tan(U.cmAngles[BHe.he]);
            A = BHe.he.prev().vertex();
        } else {
            cotA = 0;
        }
        if (BHe.he.twin().isReal()) {
            cotB = 1.0 / tan(U.cmAngles[BHe.twin().he]);
            B = BHe.twin().he.prev().vertex();
        } else {
            cotB = 0;
        }

        // check for singularities
        if (BHe.he.isReal() && singularities[A] != 0 && 
            BHe.he.twin().isReal() && singularities[B] != 0) {
            throw std::logic_error("halfedge between 2 singularities");
        } else if (BHe.he.isReal() && singularities[A] != 0) {
            cotA = 0;                  
        } else if (BHe.he.twin().isReal() && singularities[B] != 0) {
            cotB = 0;
        }
        cotWeight = (cotA + cotB) / 2.0;

        // For vertex on some sheet, returns indices to index real and imaginary parts of psi,
        // as well as sign for imaginary componenet to implement conjugate
        auto getIndexData = [&](BVertex Bv, size_t& realInd, size_t& imagInd, double& imagSign) {
            bool v_singular = (singularities[Bv.v] != 0);
            if (v_singular) {  
                if (Bv.sheet != BC.singularSheet) throw std::logic_error("wrong sheet for singular vertex");
                realInd = 2 * BVertexIndices[BC.singularSheet][Bv.v];    // all singularities live on a single sheet (0)
                imagSign = 1;
            } else if ((!v_singular && Bv.sheet == 0) || (!v_singular && Bv.sheet == 1)) {
                realInd = 2 * BVertexIndices[Bv.sheet][Bv.v];            // sheet 0 and 1 treated as normal
                imagSign = 1;
            } else if (!v_singular && Bv.sheet == 2) {
                realInd = 2 * BVertexIndices[0][Bv.v];                   // sheet 2 shares indices with sheet 0
                imagSign = -1;
            } else if (!v_singular && Bv.sheet == 3) {
                realInd = 2 * BVertexIndices[1][Bv.v];                   // sheet 3 shares indices with sheet 1
                imagSign = -1;
            }
            imagInd = realInd+1;
        };
   
        // get vertex indices and imaginary sign
        BVertex Bv_A, Bv_B;
        if (BHe.he.isReal()) {
            Bv_A = BHe.vertex();
            Bv_B = BHe.next().vertex();
        } else {
            Bv_A = BHe.twin().next().vertex();
            Bv_B = BHe.twin().vertex();
        }

        size_t iA_re, iA_im, iB_re, iB_im;
        double sA_im, sB_im;
        getIndexData(Bv_A, iA_re, iA_im, sA_im);
        getIndexData(Bv_B, iB_re, iB_im, sB_im);

        // add diagonal entries
        triplets.push_back( Eigen::Triplet<double>(iA_re,iA_re,cotWeight) );
        triplets.push_back( Eigen::Triplet<double>(iA_im,iA_im,cotWeight) );

        // transport coefficient components
        BEdge Be = BHe.edge();
        double sign = (Be.halfedge() == BHe) ? 1 : -1;
        double omega_ij = sign * omega[Be.sheet][Be.e];
        double x = -cotWeight * cos(omega_ij);
        double yA = sA_im * -cotWeight * sin(omega_ij);
        double yB = sB_im * -cotWeight * sin(omega_ij);
        double xAB = sA_im * sB_im * -cotWeight * cos(omega_ij);

        // add non-diagonal entries
        triplets.push_back( Eigen::Triplet<double>(iA_re,iB_re,  x) ); triplets.push_back( Eigen::Triplet<double>(iA_re,iB_im, yB) );
        triplets.push_back( Eigen::Triplet<double>(iA_im,iB_re,-yA) ); triplets.push_back( Eigen::Triplet<double>(iA_im,iB_im,xAB) );
    }
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

Eigen::SparseMatrix<double> QuadMesh::massMatrix() {
    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular); // + numSingularities);
    Eigen::SparseMatrix<double> M(numPsi, numPsi);
    std::vector<Eigen::Triplet<double>> triplets;

    auto getIndexData = [&](BVertex Bv, size_t& realInd, size_t& imagInd) {
        bool v_singular = (singularities[Bv.v] != 0);
        if (v_singular) {  
            if (Bv.sheet != BC.singularSheet) throw std::logic_error("wrong sheet for singular vertex");
            realInd = 2 * BVertexIndices[BC.singularSheet][Bv.v];    // all singularities live on a single sheet (0)
        } else if ((!v_singular && Bv.sheet == 0) || (!v_singular && Bv.sheet == 1)) {
            realInd = 2 * BVertexIndices[Bv.sheet][Bv.v];            // sheet 0 and 1 treated as normal
        } else if (!v_singular && Bv.sheet == 2) {
            realInd = 2 * BVertexIndices[0][Bv.v];                   // sheet 2 shares indices with sheet 0
        } else if (!v_singular && Bv.sheet == 3) {
            realInd = 2 * BVertexIndices[1][Bv.v];                   // sheet 3 shares indices with sheet 1
        }
        imagInd = realInd+1;
    };
    
    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        BHalfedge BHe = Bf.halfedge();

        // skip faces that contain a singular vertex
        VertexPtr v_i = BHe.vertex().v;
        VertexPtr v_j = BHe.next().vertex().v;
        VertexPtr v_k = BHe.next().next().vertex().v;
        if (singularities[v_i] != 0 || singularities[v_j] != 0 || singularities[v_k] != 0) {
            continue;
        }

        size_t iA_re, iA_im, iB_re, iB_im, iC_re, iC_im;
        getIndexData(BHe.vertex(),               iA_re, iA_im);
        getIndexData(BHe.next().vertex(),        iB_re, iB_im);
        getIndexData(BHe.next().next().vertex(), iC_re, iC_im);

        double area = U.cmAreas[Bf.f]; 
        
        // golurkin mass matrix
        triplets.push_back(Eigen::Triplet<double>(iA_re, iA_re, area/6.));
        triplets.push_back(Eigen::Triplet<double>(iA_im, iA_im, area/6.));
        triplets.push_back(Eigen::Triplet<double>(iA_re, iB_re, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iA_im, iB_im, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iA_re, iC_re, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iA_im, iC_im, area/12.));

        triplets.push_back(Eigen::Triplet<double>(iB_re, iB_re, area/6.));
        triplets.push_back(Eigen::Triplet<double>(iB_im, iB_im, area/6.));
        triplets.push_back(Eigen::Triplet<double>(iB_re, iA_re, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iB_im, iA_im, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iB_re, iC_re, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iB_im, iC_im, area/12.));

        triplets.push_back(Eigen::Triplet<double>(iC_re, iC_re, area/6.));
        triplets.push_back(Eigen::Triplet<double>(iC_im, iC_im, area/6.));
        triplets.push_back(Eigen::Triplet<double>(iC_re, iB_re, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iC_im, iB_im, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iC_re, iA_re, area/12.));
        triplets.push_back(Eigen::Triplet<double>(iC_im, iA_im, area/12.));
        
        /*
        // lumped mass matrix
        triplets.push_back(Eigen::Triplet<double>(iA_re, iA_re, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iA_im, iA_im, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iB_re, iB_re, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iB_im, iB_im, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iC_re, iC_re, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iC_im, iC_im, area/3.));
        */
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

// builds the constrained least squares matrix
Eigen::SparseMatrix<double> QuadMesh::linkMatrix() {
    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular);
    size_t numBEdges = 0;
    for (VertexPtr v : mesh->vertices()) {
        if (singularities[v] != 0) {
            numBEdges += 2 * 4 * v.degree();
        }
    }
    Eigen::SparseMatrix<double> C(numPsi+numBEdges, numPsi+numBEdges);
    std::vector<Eigen::Triplet<double>> triplets;

    // For vertex on some sheet, returns indices to index real and imaginary parts of psi,
    // as well as sign for imaginary componenet to implement conjugate
    auto getIndexData = [&](BVertex Bv, size_t& realInd, size_t& imagInd, double& imagSign) {
        bool v_singular = (singularities[Bv.v] != 0);
        if (v_singular) {  
            if (Bv.sheet != BC.singularSheet) throw std::logic_error("wrong sheet for singular vertex");
            realInd = 2 * BVertexIndices[BC.singularSheet][Bv.v];    // all singularities live on a single sheet (0)
            imagSign = 1;
        } else if ((!v_singular && Bv.sheet == 0) || (!v_singular && Bv.sheet == 1)) {
            realInd = 2 * BVertexIndices[Bv.sheet][Bv.v];            // sheet 0 and 1 treated as normal
            imagSign = 1;
        } else if (!v_singular && Bv.sheet == 2) {
            realInd = 2 * BVertexIndices[0][Bv.v];                   // sheet 2 shares indices with sheet 0
            imagSign = -1;
        } else if (!v_singular && Bv.sheet == 3) {
            realInd = 2 * BVertexIndices[1][Bv.v];                   // sheet 3 shares indices with sheet 1
            imagSign = -1;
        }
        imagInd = realInd+1;
    };

    // upper left corner = 2*identity matrix
    for (size_t i = 0; i < numPsi; i++) {
        triplets.push_back( Eigen::Triplet<double>(i,i,2) );
    }

    // bottom left and upper right corners
    size_t iBEdge = 0;
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] == 0) continue;

        BHalfedge BHe = Bv.halfedge();
        size_t deg = Bv.v.degree();

        for (size_t i = 0; i < 4*deg; i++) {
            BHalfedge BHe_link = BHe.next();
            BEdge Be_link = BHe_link.edge();

            // For vertex on some sheet, returns indices to index real and imaginary parts of psi,
            // as well as sign for imaginary componenet to implement conjugate
            auto getIndexData = [&](BVertex Bv, size_t& realInd, size_t& imagInd, double& imagSign) {
                bool v_singular = (singularities[Bv.v] != 0);
                if (v_singular) {  
                    if (Bv.sheet != BC.singularSheet) throw std::logic_error("wrong sheet for singular vertex");
                    realInd = 2 * BVertexIndices[BC.singularSheet][Bv.v];    // all singularities live on a single sheet (0)
                    imagSign = 1;
                } else if ((!v_singular && Bv.sheet == 0) || (!v_singular && Bv.sheet == 1)) {
                    realInd = 2 * BVertexIndices[Bv.sheet][Bv.v];            // sheet 0 and 1 treated as normal
                    imagSign = 1;
                } else if (!v_singular && Bv.sheet == 2) {
                    realInd = 2 * BVertexIndices[0][Bv.v];                   // sheet 2 shares indices with sheet 0
                    imagSign = -1;
                } else if (!v_singular && Bv.sheet == 3) {
                    realInd = 2 * BVertexIndices[1][Bv.v];                   // sheet 3 shares indices with sheet 1
                    imagSign = -1;
                }
                imagInd = realInd+1;
            };
    
            // get vertex indices and imaginary sign
            BVertex Bv_A = BHe_link.vertex();
            BVertex Bv_B = BHe_link.next().vertex();

            size_t iA_re, iA_im, iB_re, iB_im;
            double sA_im, sB_im;
            getIndexData(Bv_A, iA_re, iA_im, sA_im);
            getIndexData(Bv_B, iB_re, iB_im, sB_im);

            // transport coefficient components
            double sign = (Be_link.halfedge() == BHe_link) ? 1 : -1;
            double omega_ij = sign * omega[Be_link.sheet][Be_link.e];
            double x = -cos(omega_ij);
            double yA = -sA_im * sin(omega_ij);
            double yB = -sB_im * sin(omega_ij);
            double xAB = -sA_im * sB_im * cos(omega_ij);

            size_t row_re = iBEdge + numPsi;
            size_t row_im = row_re + 1;

            // build C
            triplets.push_back( Eigen::Triplet<double>(row_re,iA_re,1) );
            triplets.push_back( Eigen::Triplet<double>(row_im,iA_im,1) );

            triplets.push_back( Eigen::Triplet<double>(row_re,iB_re,  x) ); triplets.push_back( Eigen::Triplet<double>(row_re,iB_im, yB) );
            triplets.push_back( Eigen::Triplet<double>(row_im,iB_re,-yA) ); triplets.push_back( Eigen::Triplet<double>(row_im,iB_im,xAB) );

            // build C^T
            triplets.push_back( Eigen::Triplet<double>(iA_re,row_re,1) );
            triplets.push_back( Eigen::Triplet<double>(iA_im,row_im,1) );

            triplets.push_back( Eigen::Triplet<double>(iB_re,row_re,  x) ); triplets.push_back( Eigen::Triplet<double>(iB_im,row_re, yB) );
            triplets.push_back( Eigen::Triplet<double>(iB_re,row_im,-yA) ); triplets.push_back( Eigen::Triplet<double>(iB_im,row_im,xAB) );

            iBEdge += 2;
            BHe = BHe.next().next().twin();
        }
    }
    if (iBEdge != numBEdges) {
        throw std::logic_error("built link constraint matrix wrong");
    }
    C.setFromTriplets(triplets.begin(),triplets.end());
    return C;
}

/* 
double QuadMesh::computeStripes() {
    std::cout << "Computing Stripes..." << std::endl;
    int nPowerIterations = 20;
    
    // build matrices for inverse power method
    Eigen::SparseMatrix<double> A = energyMatrix();
    Eigen::SparseMatrix<double> B = massMatrix();

    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular); // + numSingularities);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(numPsi,1);
    Eigen::MatrixXd prevX, num, denom;
    double lambda;

    Eigen::SparseMatrix<double> C = linkMatrix();

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    PositiveDefiniteSolver<double> s(A);
    Solver<double> s2(C);
    for (int i = 0; i < nPowerIterations; i++) {
        prevX = x;
        
        Eigen::MatrixXd rhs = B * x;
        x = s.solve(rhs);

        // project onto the Cx = 0 constraint
        Eigen::MatrixXd b = Eigen::MatrixXd::Zero(C.rows(),1);
        for (size_t j = 0; j < numPsi; j++) {
            b(j,0) = 2 * x(j,0);
        }
        Eigen::MatrixXd y = s2.solve(b);
        for (size_t j = 0; j < numPsi; j++) {
            x(j,0) = y(j,0);
        }

        // normalize result and check residual
        double norm2 = (x.transpose() * B * x)(0,0);
        x = x / std::sqrt(norm2);

        num = x.transpose() * A * x;
        denom = x.transpose() * B * x;
        lambda = num(0,0) / denom(0,0);
        double resid = (A*x - lambda * B * x).norm();
        std::cout << "Resid: " << i << "," << resid << std::endl;
    }

    // store resulting field and naive coords
    for (int i = 0; i < 2; i++) {
        for (VertexPtr v : mesh->vertices()) {
            if (singularities[v] != 0) continue;
            
            size_t index = BVertexIndices[i][v];
            double real = x(2*index,  0);
            double imag = x(2*index+1,0);
            
            std::complex<double> p(real,imag);
            psi[i][v] = p;
            coords[i][v] = std::arg(p);
        }
    }

    std::cout << "Done!" << std::endl;
    return lambda;
}
*/
  
void QuadMesh::computeStripes() {
    // index BVertices with singular vertices coming first, then sheet 0 and sheet 1 vertices
    //size_t iSingular = 0;
    //size_t iNonSingular = numSingularities;
    size_t iNonSingular = 0;
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] != 0) {
            //assert(Bv.sheet == BC.singularSheet);
            //BVertexIndices[BC.singularSheet][Bv.v] = iSingular++;
            continue;
        } else if (Bv.sheet == 0 || Bv.sheet == 1) {
            BVertexIndices[Bv.sheet][Bv.v] = iNonSingular++;
        }
    }

    // first, fix the psi values along each link of a branch point
    std::vector<VertexData<int>> allLinkVertices;
    for (size_t i = 0; i < 4; i++) {
        VertexData<int>linkSheet(mesh,0);
        allLinkVertices.push_back(linkSheet);
    }
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] == 0) continue;

        size_t deg = Bv.v.degree();
        BHalfedge BHe = Bv.halfedge();

        // compute a psi for each vertex along the link
        for (size_t i = 0; i < 4*deg; i++) {
            BHalfedge currBHe = BHe;

            double W = 0;
            // compute psi by integrating omega along the halflink of the vertex
            for (size_t j = 0; j < 2*deg; j++) {
                BHalfedge BHe_link = currBHe.next();
                BEdge Be_link = BHe_link.edge();

                int sign_ij = (Be_link.halfedge() == BHe_link) ? 1 : -1;
                double omega_ij = sign_ij * omega[Be_link.sheet][Be_link.e];
                W += omega_ij;

                currBHe = currBHe.next().next().twin();
            }
            
            BVertex Bv_i = BHe.next().vertex();
            double theta_i = -W / 2.0;
            std::complex<double> psi_i = std::exp(IM_I * theta_i);
            if (Bv_i.sheet == 0 || Bv_i.sheet == 1) {
                psi[Bv_i.sheet][Bv_i.v] = psi_i;
                coords[Bv_i.sheet][Bv_i.v] = std::arg(psi_i);
            } else if (Bv_i.sheet == 2) {
                psi[0][Bv_i.v] = std::conj(psi_i);
                coords[0][Bv_i.v] = std::arg(std::conj(psi_i));
            } else if (Bv_i.sheet == 3) {
                psi[1][Bv_i.v] = std::conj(psi_i);
                coords[1][Bv_i.v] = std::arg(std::conj(psi_i));
            } else {
                throw std::logic_error("impossible sheet");
            }
            allLinkVertices[Bv_i.sheet][Bv_i.v] = 1;
            BHe = BHe.next().next().twin();
        }
    }

    // now just minimize the energy by computing a solve with Dirichlet boundary conditions
    size_t iN = 0;
    size_t iB = 0;
    std::vector<VertexData<size_t>> interiorVertexIndices;
    std::vector<VertexData<size_t>> boundaryVertexIndices;
    for (size_t i = 0; i < 4; i++) {
        VertexData<size_t> indexSheet(mesh);
        VertexData<size_t> indexSheet2(mesh);
        interiorVertexIndices.push_back(indexSheet);
        boundaryVertexIndices.push_back(indexSheet2);
    }
    Eigen::SparseMatrix<double> A = energyMatrix();
    Eigen::Array<bool,Eigen::Dynamic,1> interior(A.rows(),1);

    for (BVertex Bv : allBVertices) {
        if(Bv.sheet == 2 || Bv.sheet == 3) continue;
        if (singularities[Bv.v] != 0) continue;

        if (allLinkVertices[Bv.sheet][Bv.v] == 0) {
            interiorVertexIndices[Bv.sheet][Bv.v] = iN++;
        } else {
            boundaryVertexIndices[Bv.sheet][Bv.v] = iB++;
        }
        interior(2*BVertexIndices[Bv.sheet][Bv.v],  0) = (allLinkVertices[Bv.sheet][Bv.v] == 0);
        interior(2*BVertexIndices[Bv.sheet][Bv.v]+1,0) = (allLinkVertices[Bv.sheet][Bv.v] == 0);
    }

    BlockDecompositionResult<double> A_all = blockDecomposeSquare(A, interior);
    Eigen::SparseMatrix<double> A_II = A_all.AA;
    Eigen::SparseMatrix<double> A_IB = A_all.AB;

    // build the RHS
    Eigen::MatrixXd x_b(2*iB,1);
    for (BVertex Bv : allBVertices) {
        if (Bv.sheet == 2 || Bv.sheet == 3) continue;
        if (singularities[Bv.v] != 0) continue;
        if (allLinkVertices[Bv.sheet][Bv.v] == 1) {
            std::complex<double> psi_v = psi[Bv.sheet][Bv.v];
            size_t index = boundaryVertexIndices[Bv.sheet][Bv.v];
            x_b(2*index,  0) = psi_v.real();
            x_b(2*index+1,0) = psi_v.imag();
        }
    }
    Eigen::MatrixXd rhs = -A_IB * x_b;

    // solve and extract result
    Eigen::MatrixXd x_i = solvePositiveDefinite<double>(A_II,rhs);
    for (BVertex Bv : allBVertices) {
        if (Bv.sheet == 2 || Bv.sheet == 3) continue;
        if (singularities[Bv.v] != 0) continue;
        if (allLinkVertices[Bv.sheet][Bv.v] == 0) {
            size_t index = interiorVertexIndices[Bv.sheet][Bv.v];
            double real = x_i(2*index,  0);
            double imag = x_i(2*index+1,0);
            std::complex<double> psi_v(real,imag);

            psi[Bv.sheet][Bv.v] = psi_v;
            coords[Bv.sheet][Bv.v] = std::arg(psi_v);
        }
    }
}

/*
// local field per face
void QuadMesh::optimizeSimpleLocally() {
    std::cout << "optimizing quad mesh output..." << std::endl;
    // psi = ae^itheta
    // theta = arg(psi)
    // compute updated omega = dtheta = difference between vertices

    // set up an edgedata per sheet
    std::vector<EdgeData<double>> newOmega;
    for (int i = 0; i < n; i++) {
        EdgeData<double> omegaSheet(mesh);
        newOmega.push_back(omegaSheet);
    }

    computeSigma();
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BVertex Bv_i = Be.halfedge().vertex();
        BVertex Bv_j = Be.halfedge().next().vertex();

        // if contains a singular vertex, copy over old value of omega
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0) {
            newOmega[Be.sheet][Be.e] = omega[Be.sheet][Be.e];
        } else {
            newOmega[Be.sheet][Be.e] = getSigma(Be);
        }
    }

    // compute new face vectors
    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        BHalfedge BHe_ij = Bf.halfedge();
        BHalfedge BHe_jk = BHe_ij.next();
        BHalfedge BHe_ki = BHe_jk.next();

        BEdge Be_ij = BHe_ij.edge();
        BEdge Be_jk = BHe_jk.edge();
        BEdge Be_ki = BHe_ki.edge();

        double sign_ij = (Be_ij.halfedge() == BHe_ij) ? 1 : -1;
        double sign_jk = (Be_jk.halfedge() == BHe_jk) ? 1 : -1;
        double sign_ki = (Be_ki.halfedge() == BHe_ki) ? 1 : -1;

        double omega_ij = sign_ij * newOmega[Be_ij.sheet][Be_ij.e];
        double omega_jk = sign_jk * newOmega[Be_jk.sheet][Be_jk.e];
        double omega_ki = sign_ki * newOmega[Be_ki.sheet][Be_ki.e];

        std::complex<double> e_ij = thetaCM[BHe_ij.he];
        std::complex<double> e_jk = thetaCM[BHe_jk.he];
        std::complex<double> e_ki = thetaCM[BHe_ki.he];

        // build 3x2 linear system
        Eigen::MatrixXd A(3,2);
        Eigen::MatrixXd b(3,1);
        Eigen::MatrixXd x;

        A(0,0) = e_ij.real();
        A(0,1) = e_ij.imag();
        b(0,0) = omega_ij;

        A(1,0) = e_jk.real();
        A(1,1) = e_jk.imag();
        b(1,0) = omega_jk;

        A(2,0) = e_ki.real();
        A(2,1) = e_ki.imag();
        b(2,0) = omega_ki;

        //Eigen::MatrixXd oldX(2,1);
        //oldX(0,0) = 100 * branchCoverFields[Bf.sheet][Bf.f].real();
        //oldX(1,0) = 100 * branchCoverFields[Bf.sheet][Bf.f].imag();

        //Eigen::MatrixXd AoldX = A * oldX;
        //std::cout << "AoldX: " << AoldX << std::endl;
        //std::cout << "--------------" << std::endl;
        //std::cout<< "new: " << b << std::endl;

        x = A.colPivHouseholderQr().solve(b);
        std::complex<double> z(x(0,0),x(1,0));
        //std::cout << "mag: " << std::abs(z) << std::endl;
        //std::cout << std::abs(100.0 * branchCoverFields[Bf.sheet][Bf.f] - z) << std::endl;
        branchCoverFields[Bf.sheet][Bf.f] = z;
    }
    computeOmega(1.0);
    computeStripes();
}

// global field 
void QuadMesh::optimizeSimpleGlobally() {
 std::cout << "optimizing quad mesh output..." << std::endl;
    // set up an edgedata per sheet
    std::vector<EdgeData<double>> newOmega;
    std::vector<EdgeData<size_t>> edgeIndices;
    for (int i = 0; i < n; i++) {
        EdgeData<double> omegaSheet(mesh);
        EdgeData<size_t> indexSheet(mesh);
        newOmega.push_back(omegaSheet);
        edgeIndices.push_back(indexSheet);
    }

    // compute 1-form on edges and index them
    size_t iEdge = 0;
    computeSigma();
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BVertex Bv_i = Be.halfedge().vertex();
        BVertex Bv_j = Be.halfedge().next().vertex();

        // if contains a singular vertex, copy over old value of omega
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0) {
            newOmega[Be.sheet][Be.e] = omega[Be.sheet][Be.e];
        } else {
            newOmega[Be.sheet][Be.e] = getSigma(Be);
        }

        edgeIndices[Be.sheet][Be.e] = iEdge++; 
    }

    // compute new face vectors
    size_t nBEdges = 4 * mesh->nEdges();
    Eigen::MatrixXd A(nBEdges,2);
    Eigen::MatrixXd b(nBEdges,1);
    Eigen::MatrixXd x;
    for (BEdge Be : allBEdges) {
        BHalfedge BHe_ij, BHe_ji;
        BFace Bf_ij, Bf_ji;
        std::complex<double> e_ij, e_ji, xAxis_ij, xAxis_ji, yAxis_ij, yAxis_ji;
        double boundaryScale = 1.0;
        
        if (Be.halfedge().he.isReal()) {
            BHe_ij = Be.halfedge();
            Bf_ij = BHe_ij.face();
            e_ij = thetaCM[BHe_ij.he];
            xAxis_ij = xBasis[Bf_ij.sheet][Bf_ij.f];
            yAxis_ij = IM_I * xAxis_ij;
        } else {
            e_ij = 0;
            xAxis_ij = 0;
            yAxis_ij = 0;
            boundaryScale = 2.0;
        }
        if (Be.halfedge().he.twin().isReal()) {
            BHe_ji = Be.halfedge().twin();
            Bf_ji = BHe_ji.face();
            e_ji = -thetaCM[BHe_ji.he];
            xAxis_ji = xBasis[Bf_ji.sheet][Bf_ji.f];
            yAxis_ji = IM_I * xAxis_ji;
        } else {
            e_ji = 0;
            xAxis_ji = 0;
            yAxis_ji = 0;
            boundaryScale = 2.0;
        }

        // build a 3x2 linear system locally
        size_t i = edgeIndices[Be.sheet][Be.e];
        A(i,0) = boundaryScale * 0.5 * (dot(xAxis_ij, e_ij) + dot(xAxis_ji, e_ji));
        A(i,1) = boundaryScale * 0.5 * (dot(yAxis_ij, e_ij) + dot(yAxis_ji, e_ji));
        b(i,0) = newOmega[Be.sheet][Be.e];
    }

    // gets 2x1 in basis of the root face
    x = A.colPivHouseholderQr().solve(b);
    std::cout << x << std::endl;
    
    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        std::complex<double> xAxis = xBasis[Bf.sheet][Bf.f];
        if (xAxis != branchCoverFields[Bf.sheet][Bf.f]) {
            std::cout << xAxis << "," <<  branchCoverFields[Bf.sheet][Bf.f] << std::endl;
            throw std::logic_error("bad basis value");
        }
        std::complex<double> yAxis = IM_I * xAxis;
        std::complex<double> z = x(0,0) * xAxis + x(1,0) * yAxis;
        branchCoverFields[Bf.sheet][Bf.f] = z;
    }
    
    computeOmega(1.0);
    computeStripes();
}
*/

// need to do some extra work for meshes with boundary, this is broken for the new interpolant / constraints
void QuadMesh::optimizeHarmonic(bool visualize) {
    std::cout << "Helmholz hodge decomposition..." << std::endl;

    // index each element
    std::vector<VertexData<size_t>> vertexIndices;
    std::vector<EdgeData<size_t>> edgeIndices;
    std::vector<FaceData<size_t>> faceIndices;
    for (size_t i = 0; i < 4; i++) {
        VertexData<size_t> vSheet(mesh);
        EdgeData<size_t> eSheet(mesh);
        FaceData<size_t> fSheet(mesh);
        vertexIndices.push_back(vSheet);
        edgeIndices.push_back(eSheet);
        faceIndices.push_back(fSheet);
    }

    size_t iVertex = 0;
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        vertexIndices[Bv.sheet][Bv.v] = iVertex++;
    }

    size_t iEdge = 0;
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        edgeIndices[Be.sheet][Be.e] = iEdge++;
    }

    size_t iFace = 0;
    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        faceIndices[Bf.sheet][Bf.f] = iFace++;
    }

    // compute newOmega in vector form
    computeSigma();
    Eigen::MatrixXd newOmega(4 * mesh->nEdges(),1);
    std::vector<EdgeData<double>> dTheta;
    for (size_t i = 0; i < 4; i++) {
        EdgeData<double> dThetaSheet(mesh);
        dTheta.push_back(dThetaSheet);
    }

    for (BEdge Be : allBEdges) {
        size_t i = edgeIndices[Be.sheet][Be.e];
        //newOmega(i,0) = omega[Be.sheet][Be.e];
        //continue;
        
        BVertex Bv_i = Be.halfedge().vertex();
        BVertex Bv_j = Be.halfedge().next().vertex();

        // if contains a singular vertex, copy over old value of omega
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || Be.e.isBoundary()) {
            newOmega(i,0) = getSigma(Be);
        } else {
            newOmega(i,0) = getSigma(Be);
        }

        dTheta[Be.sheet][Be.e] = newOmega(i,0);
    }

    if (visualize) {
        polyscope::registerSurfaceMesh("post-uniformization",geom);
        polyscope::getSurfaceMesh("post-uniformization")->addQuantity("dTheta 0", dTheta[0]);
        polyscope::getSurfaceMesh("post-uniformization")->addQuantity("dTheta 1", dTheta[1]);
        polyscope::getSurfaceMesh("post-uniformization")->addQuantity("dTheta 2", dTheta[2]);
        polyscope::getSurfaceMesh("post-uniformization")->addQuantity("dTheta 3", dTheta[3]);
        polyscope::show();
    }

    // compute harmonic component
    HodgeDecomposition H = HodgeDecomposition(vertexIndices, edgeIndices, faceIndices, 
                                              BC, U.cmAngles, U.cmAreas, mesh, numSingularities);
    Eigen::MatrixXd gamma = H.computeHarmonicComponent(newOmega);
    
    // now update omega
    for (BEdge Be : allBEdges) {
        size_t i = edgeIndices[Be.sheet][Be.e];
        omega[Be.sheet][Be.e] = gamma(i,0);
    }
    std::cout << "Done!" << std::endl;
    computeStripes();
}

std::complex<double> QuadMesh::getPsi(BVertex Bv) {
    if (Bv.sheet == 0 || Bv.sheet == 1) {
        return psi[Bv.sheet][Bv.v];
    } else if (Bv.sheet == 2) {
        return std::conj(psi[0][Bv.v]);
    } else {
        return std::conj(psi[1][Bv.v]);
    }
}

double QuadMesh::getSigma(BEdge Be) {
    if (Be.sheet == 0 || Be.sheet == 1) {
        return sigma[Be.sheet][Be.e];
    } else if (Be.sheet == 2) {
        return -sigma[0][Be.e];
    } else {
        return -sigma[1][Be.e];
    }
}

void QuadMesh::computeSigma() {
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        if (Be.sheet != 0 && Be.sheet != 1) continue;

        BHalfedge BHe = Be.halfedge();
        std::complex<double> psi_i, psi_j;
        psi_i = getPsi(BHe.vertex());       
        psi_j = getPsi(BHe.next().vertex());
        if (singularities[BHe.vertex().v] != 0) {
            psi_i = 0;         
        } 
        if (singularities[BHe.next().vertex().v] != 0) {
            psi_j = 0;
        } 
        double w_ij = omega[Be.sheet][Be.e];

        //double y = (w_ij + std::arg(psi_i) - std::arg(psi_j) ) / (2 * PI);
        //int n_ij = std::round(y);
        //sigma[Be.sheet][Be.e] = std::arg(psi_j) - std::arg(psi_i) + 2 * PI * n_ij;
        
        double delta_ij = std::arg( (std::exp(IM_I * w_ij) * psi_i) / psi_j);
        sigma[Be.sheet][Be.e] = w_ij - delta_ij; 
    }

    /* 
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] == 0) continue;
        BHalfedge BHe = Bv.halfedge();
        size_t deg = Bv.v.degree();
        double sigma_sum = 0;
        for (size_t i = 0; i < 4*deg; i++) {
            BHalfedge BHe_link = BHe.next();
            BEdge Be_link = BHe_link.edge();
            std::cout << "diff: " << getSigma(Be_link) - omega[Be_link.sheet][Be_link.e] << std::endl; 
            BHe = BHe.next().next().twin();
        }
    }
    */
    
}

/* 
Eigen::MatrixXd QuadMesh::buildLocalEnergy(BFace Bf) {
    // For vertex on some sheet, returns indices to index real and imaginary parts of psi,
    // as well as sign for imaginary componenet to implement conjugate
    auto getIndexData = [&](BVertex Bv, double& imagSign) {
        bool v_singular = (singularities[Bv.v] != 0);
        if (v_singular) {  
            if (Bv.sheet != BC.singularSheet) throw std::logic_error("wrong sheet for singular vertex");
            imagSign = 1;
        } else if ((!v_singular && Bv.sheet == 0) || (!v_singular && Bv.sheet == 1)) {
            imagSign = 1;
        } else if (!v_singular && Bv.sheet == 2) {
            imagSign = -1;
        } else if (!v_singular && Bv.sheet == 3) {
            imagSign = -1;
        }
    };

    BHalfedge BHe_ij = Bf.halfedge();
    BHalfedge BHe_jk = BHe_ij.next();
    BHalfedge BHe_ki = BHe_jk.next();

    BVertex Bv_i = BHe_ij.vertex();
    BVertex Bv_j = BHe_jk.vertex();
    BVertex Bv_k = BHe_ki.vertex();

    // compute cotan, zArea, and drift
    double scale = 100;
    double cot_ij = 1.0 / tan(U.cmAngles[BHe_ij.he]);
    double cot_jk = 1.0 / tan(U.cmAngles[BHe_jk.he]);
    double cot_ki = 1.0 / tan(U.cmAngles[BHe_ki.he]);
    double zArea = scale * scale * U.cmAreas[Bf.f];

    std::complex<double> z = scale * branchCoverFields[Bf.sheet][Bf.f];
    std::complex<double> e_i_perp = IM_I * thetaCM[BHe_jk.he];
    std::complex<double> e_j_perp = IM_I * thetaCM[BHe_ki.he];
    std::complex<double> e_k_perp = IM_I * thetaCM[BHe_ij.he];
    
    double drift_ij = (1.0 / 6.0) * dot(e_i_perp - e_j_perp, z); 
    double drift_jk = (1.0 / 6.0) * dot(e_j_perp - e_k_perp, z); 
    double drift_ki = (1.0 / 6.0) * dot(e_k_perp - e_i_perp, z); 
    double drift_ji = -drift_ij;
    double drift_kj = -drift_jk;
    double drift_ik = -drift_ki;

    // get sign information
    double si_im, sj_im, sk_im;
    getIndexData(Bv_i, si_im);
    getIndexData(Bv_j, sj_im);
    getIndexData(Bv_k, sk_im);

    size_t i_re = 0;
    size_t i_im = 1;
    size_t j_re = 2;
    size_t j_im = 3;
    size_t k_re = 4;
    size_t k_im = 5;

    // ---------------------------------------------------------------------------------- //

    Eigen::MatrixXd A(6,6);

    // A_ii, A_ij, A_ik
    A(i_re,i_re) = (cot_ij + cot_ki) / 2.0 + zArea / 6.0;  A(i_re,i_im) = 0;
    A(i_im,i_re) = 0;                                      A(i_im,i_im) = (cot_ij + cot_ki) / 2.0 + zArea / 6.0;

    A(i_re,j_re) = -cot_ij / 2.0 + zArea / 12.0;           A(i_re,j_im) = sj_im * drift_ij;
    A(i_im,j_re) = si_im * -drift_ij;                      A(i_im,j_im) = si_im * sj_im * (-cot_ij / 2.0 + zArea / 12.0);   

    A(i_re,k_re) = -cot_ki / 2.0 + zArea / 12.0;           A(i_re,k_im) = sk_im * drift_ik;
    A(i_im,k_re) = si_im * -drift_ik;                      A(i_im,k_im) = si_im * sk_im * (-cot_ki / 2.0 + zArea / 12.0);

    // A_ji, A_jj, A_jk
    A(j_re,i_re) = -cot_ij / 2.0 + zArea / 12.0;           A(j_re,i_im) = si_im * drift_ji;
    A(j_im,i_re) = sj_im * -drift_ji;                      A(j_im,i_im) = sj_im * si_im * (-cot_ij / 2.0 + zArea / 12.0);

    A(j_re,j_re) = (cot_jk + cot_ij) / 2.0 + zArea / 6.0;  A(j_re,j_im) = 0;
    A(j_im,j_re) = 0;                                      A(j_im,j_im) = (cot_jk + cot_ij) / 2.0 + zArea / 6.0;

    A(j_re,k_re) = -cot_jk / 2.0 + zArea / 12.0;           A(j_re,k_im) = sk_im * drift_jk;
    A(j_im,k_re) = sj_im * -drift_jk;                      A(j_im,k_im) = sj_im * sk_im * (-cot_jk / 2.0 + zArea / 12.0);

    // A_ki, A_kj, A_kk
    A(k_re,i_re) = -cot_ki / 2.0 + zArea / 12.0;           A(k_re,i_im) = si_im * drift_ki;
    A(k_im,i_re) = sk_im * -drift_ki;                      A(k_im,i_im) = sk_im * si_im * (-cot_ki / 2.0 + zArea / 12.0);

    A(k_re,j_re) = -cot_jk / 2.0 + zArea / 12.0;           A(k_re,j_im) = sj_im * drift_kj;
    A(k_im,j_re) = sk_im * -drift_kj;                      A(k_im,j_im) = sk_im * sj_im * (-cot_jk / 2.0 + zArea / 12.0);

    A(k_re,k_re) = (cot_jk + cot_ki) / 2.0 + zArea / 6.0;  A(k_re,k_im) = 0;
    A(k_im,k_re) = 0;                                      A(k_im,k_im) = (cot_jk + cot_ki) / 2.0 + zArea / 6.0;

    return A;
}

// this function assumes the psi value at Bv_i is fixed, and computes the ideal psi_j, psi_k
void QuadMesh::computeIdealPsi(BVertex Bv_i, size_t i_re, size_t j_re, size_t k_re, Eigen::MatrixXd A, 
                            std::complex<double> &psi_j, std::complex<double> &psi_k) {

    std::complex<double> psi_i = getPsi(Bv_i);
    Eigen::MatrixXd psi_i_mat(2,1);
    psi_i_mat(0,0) = psi_i.real();
    psi_i_mat(1,0) = psi_i.imag();

    Eigen::MatrixXd A_ji = A.block(j_re,i_re,2,2);
    Eigen::MatrixXd A_ki = A.block(k_re,i_re,2,2);
    Eigen::MatrixXd b_ji = A_ji * psi_i_mat;
    Eigen::MatrixXd b_ki = A_ki * psi_i_mat;

    Eigen::MatrixXd b(4,1);
    b(0,0) = -b_ji(0,0);
    b(1,0) = -b_ji(1,0);
    b(2,0) = -b_ki(0,0);
    b(3,0) = -b_ki(1,0);

    Eigen::MatrixXd A_jk(4,4);
    A_jk.block<2,2>(0,0) = A.block(j_re,j_re,2,2);
    A_jk.block<2,2>(0,2) = A.block(j_re,k_re,2,2);
    A_jk.block<2,2>(2,0) = A.block(k_re,j_re,2,2);
    A_jk.block<2,2>(2,2) = A.block(k_re,k_re,2,2);

    Eigen::MatrixXd x = A_jk.colPivHouseholderQr().solve(b);
    psi_j = std::complex<double>(x(0,0),x(1,0));
    psi_k = std::complex<double>(x(2,0),x(3,0));
}

void QuadMesh::computeSigma2() {
    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        BHalfedge BHe_ij = Bf.halfedge();
        BHalfedge BHe_jk = BHe_ij.next();
        BHalfedge BHe_ki = BHe_jk.next();
        BEdge Be_ij = BHe_ij.edge();
        BEdge Be_jk = BHe_jk.edge();
        BEdge Be_ki = BHe_ki.edge();
        BVertex Bv_i = BHe_ij.vertex();
        BVertex Bv_j = BHe_jk.vertex();
        BVertex Bv_k = BHe_ki.vertex();

        // skip over faces that contain a singular vertex
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || singularities[Bv_k.v] != 0) {
            continue;
        }

        Eigen::MatrixXd A = buildLocalEnergy(Bf);

        int sign_ij = (Be_ij.e.halfedge() == BHe_ij.he) ? 1 : -1;
        int sign_jk = (Be_jk.e.halfedge() == BHe_jk.he) ? 1 : -1;
        int sign_ki = (Be_ki.e.halfedge() == BHe_ki.he) ? 1 : -1;
        double omega_ij = sign_ij * omega[Be_ij.sheet][Be_ij.e];
        double omega_jk = sign_jk * omega[Be_jk.sheet][Be_jk.e];
        double omega_ki = sign_ki * omega[Be_ki.sheet][Be_ki.e];
        double omega_ji = -omega_ij;
        double omega_kj = -omega_jk;
        double omega_ik = -omega_ki;

        size_t i_re = 0;
        size_t j_re = 2;
        size_t k_re = 4;

        std::complex<double> psi_i = getPsi(Bv_i);
        std::complex<double> psi_j = getPsi(Bv_j);
        std::complex<double> psi_k = getPsi(Bv_k);

        // assuming psi_i is fixed, compute the ideal psi_j, psi_k values
        std::complex<double> psi_ij, psi_ik;
        computeIdealPsi(Bv_i, i_re, j_re, k_re, A, psi_ij, psi_ik);
        double sigma_ij = omega_ij - std::arg(psi_ij / psi_i);
        double sigma_ik = omega_ik - std::arg(psi_ik / psi_k);

        // assuming psi_j is fixed, compute the ideal psi_i, psi_k values
        std::complex<double> psi_ji, psi_jk;
        computeIdealPsi(Bv_j, j_re, i_re, k_re, A, psi_ji, psi_jk);
        double sigma_ji = omega_ji - std::arg(psi_ji / psi_i);
        double sigma_jk = omega_jk - std::arg(psi_jk / psi_k);

        // assuming psi_k is fixed, compute the ideal psi_i, psi_j values
        std::complex<double> psi_ki, psi_kj;
        computeIdealPsi(Bv_k, k_re, i_re, j_re, A, psi_ki, psi_kj);
        double sigma_ki = omega_ki - std::arg(psi_ki / psi_i);
        double sigma_kj = omega_kj - std::arg(psi_kj / psi_j);

        sigmaHe[BHe_ij.sheet][BHe_ij.he] = sigma_ij;
        sigmaHe[BHe_jk.sheet][BHe_jk.he] = sigma_jk;
        sigmaHe[BHe_ki.sheet][BHe_ki.he] = sigma_ki;
        //std::cout << sigma_ij << "," << sigma_ji << std::endl;
    }

    // check that sigma_ij = -sigma_ji
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BHalfedge BHe_ij = Be.halfedge();
        BHalfedge BHe_ji = Be.halfedge().twin();
        double sigma_ij = sigmaHe[BHe_ij.sheet][BHe_ij.he];
        double sigma_ji = sigmaHe[BHe_ji.sheet][BHe_ji.he];

        if (std::abs(sigma_ij + sigma_ji) > 1e-6) {
            std::cout << "diff: " << std::abs(sigma_ij + sigma_ji) << std::endl;
            throw std::logic_error("sigma's don't agree");
        }
    }
}
*/

bool QuadMesh::textureCoordinates() {
    std::cout << "Computing texture coordinates..." << std::endl;
    FaceData<int> singularFace(mesh,0);

    computeSigma();
    zeros = FaceData<int>(mesh,0);
    int numNewSingularities = 0;

    // COMPUTE SIGMAS ON EDGES POINTING TOWARDS THE BRANCH POINT
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] == 0) continue;

        BHalfedge BHe = Bv.halfedge();
        BEdge Be = BHe.edge();

        // fix an initial sigma
        if (Be.sheet == 0 || Be.sheet == 1) {
            sigma[Be.sheet][Be.e] = omega[Be.sheet][Be.e];
        } else if (Be.sheet == 2) {
            sigma[0][Be.e] = -omega[Be.sheet][Be.e];
        } else if (Be.sheet == 3) {
            sigma[1][Be.e] = -omega[Be.sheet][Be.e];
        }

        size_t deg = Bv.v.degree();
        for (size_t i = 0; i < 4*deg; i++) {
            BHalfedge BHe_ij = BHe;
            BHalfedge BHe_jk = BHe_ij.next();
            BHalfedge BHe_ki = BHe_jk.next();

            BEdge Be_ij = BHe_ij.edge();
            BEdge Be_jk = BHe_jk.edge();
            BEdge Be_ki = BHe_ki.edge();

            int sign_ij = (Be_ij.halfedge() == BHe_ij) ? 1 : -1;
            int sign_jk = (Be_jk.halfedge() == BHe_jk) ? 1 : -1;
            int sign_ki = (Be_ki.halfedge() == BHe_ki) ? 1 : -1;

            double sigma_ij = sign_ij * getSigma(Be_ij);
            double sigma_jk = sign_jk * getSigma(Be_jk);
            double sigma_ki = -sigma_ij - sigma_jk;
            if (Be_ki.sheet == 0 || Be_ki.sheet == 1) {
                sigma[Be_ki.sheet][Be_ki.e] = sign_ki * sigma_ki; 
            } else if (Be_ki.sheet == 2) {
                sigma[0][Be_ki.e] = -(sign_ki * sigma_ki);
            } else if (Be_ki.sheet == 3) {
                sigma[1][Be_ki.e] = -(sign_ki * sigma_ki);
            }
            BHe = BHe.next().next().twin();
        }
    }

    HalfedgeData<double> xCoords(mesh,0);
    HalfedgeData<double> yCoords(mesh,0);

    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        // sheet 0 = x coordinate, sheet 1 = y coordinate
        if (Bf.sheet != 0 && Bf.sheet != 1) continue;

        // grab halfedges, vertices, and edges
        BHalfedge Bhe_ij = Bf.halfedge();
        BHalfedge Bhe_jk = Bhe_ij.next();
        BHalfedge Bhe_ki = Bhe_jk.next();

        // skip singularities for now, since they have no psi associated with them
        if (singularities[Bhe_ij.vertex().v] != 0 || singularities[Bhe_jk.vertex().v] != 0 || singularities[Bhe_ki.vertex().v] != 0) {
            singularFace[Bf.f] = 1;

            if (singularities[Bhe_ij.vertex().v] != 0) {
                Bhe_ij = Bf.halfedge().next();
                Bhe_jk = Bhe_ij.next();
                Bhe_ki = Bhe_jk.next();
            }
            //continue;
        }

        BVertex Bv_i = Bhe_ij.vertex();
        BVertex Bv_j = Bhe_jk.vertex();
        BVertex Bv_k = Bhe_ki.vertex();

        BEdge Be_ij = Bhe_ij.edge();
        BEdge Be_jk = Bhe_jk.edge();
        BEdge Be_ki = Bhe_ki.edge();

        // get orientations of edges
        int c_ij = (Be_ij.halfedge() == Bhe_ij) ? 1 : -1;
        int c_jk = (Be_jk.halfedge() == Bhe_jk) ? 1 : -1;
        int c_ki = (Be_ki.halfedge() == Bhe_ki) ? 1 : -1;

        // retrieve sigmas
        double sigma_ij = c_ij * getSigma(Be_ij);
        double sigma_jk = c_jk * getSigma(Be_jk);
        double sigma_ki = c_ki * getSigma(Be_ki);

        std::complex<double> z_i = getPsi(Bv_i);
        
        // compute coordinates
        double c_i = std::arg(z_i);
        double c_j = c_i + sigma_ij; 
        double c_k = c_j + sigma_jk;

        double c_l = c_k + sigma_ki;
        double n = (c_l - c_i) / (2.0 * PI);
        
        if (std::round(n) != 0 && singularFace[Bf.f] == 0) {
            zeros[Bf.f] = std::round(n);
            numNewSingularities++;        
        } else if (std::abs(n) > 1e-6) {
            std::cout << "sigmas don't close: " << std::abs(n) << std::endl;
        }
        
        std::vector<double> currCoords = {c_i, c_j, c_k};
        if (Bf.sheet == 0) {
            xCoords[Bhe_ij.he] = c_i;
            xCoords[Bhe_jk.he] = c_j;
            xCoords[Bhe_ki.he] = c_k;
        } else {
            yCoords[Bhe_ij.he] = c_i;
            yCoords[Bhe_jk.he] = c_j;
            yCoords[Bhe_ki.he] = c_k;
        }
    }

    /* 
    // now take care of singular vertices and their neighboorhood of faces
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] == 0) continue;
        size_t numIters = 2 * Bv.v.degree();

        for (int i = 0; i < 2; i++) {
            BHalfedge BHe_orig = Bv.halfedge();
            BHe_orig.sheet = i;
            BHalfedge BHe = BHe_orig;

            std::complex<double> psi_i = getPsi(BHe.next().vertex());
            double Bi = std::arg(psi_i);
            double Bj;
            double Bm = Bi;

            // compute midpoint coordinate
            double sigma_sum = 0;
            for (size_t j = 0; j < numIters; j++) {
                int sign_ij = (BHe.next().edge().halfedge() == BHe.next()) ? 1 : -1;
                double sigma_ij = sign_ij * getSigma(BHe.next().edge());

                Bm += sigma_ij / 2.0;
                BHe = BHe.next().next().twin();

                sigma_sum += sigma_ij;
                //std::cout << "he: " << currHe << " sheet: " << BHe.sheet << " sigma_sum: " << sigma_sum << std::endl;
                //std::cout << "--------------------------------" << std::endl;
            } 

            Bm = (2.*M_PI) * std::round( Bm / (2.*M_PI) );

            // compute other coordinates
            BHe = BHe_orig;
            for (size_t j = 0; j < numIters; j++) {
                int sign_ij = (BHe.next().edge().halfedge() == BHe.next()) ? 1 : -1;
                double sigma_ij = sign_ij * getSigma(BHe.next().edge());
                Bj = Bi + sigma_ij;

                if (i == 0) {
                    xCoords[BHe.he] = Bm;
                    xCoords[BHe.next().he] = Bi;
                    xCoords[BHe.next().next().he] = Bj;
                } else {
                    yCoords[BHe.he] = Bm;
                    yCoords[BHe.next().he] = Bi;
                    yCoords[BHe.next().next().he] = Bj;
                }

                Bi = Bj;
                BHe = BHe.next().next().twin();
            } 
        }
    }
    */

    // prepare data for shader
    for (FacePtr f : mesh->faces()) {
        std::vector<Vector2> coords; 
        HalfedgePtr he = f.halfedge();
        do {
            Vector2 vertCoords;
            if (singularFace[f] != 0) {   
                vertCoords = Vector2{ xCoords[he], yCoords[he] };
                //vertCoords = {0, 0}; // default value for singularities
            } else {
                vertCoords = Vector2{ xCoords[he], yCoords[he] };
            }
            coords.push_back(vertCoords);
            he = he.next();
        } while (he != f.halfedge());
        texCoords[f] = coords;
    }

    std::cout << "Done! Num new singularities: " << numNewSingularities << std::endl;
    return (numNewSingularities == 0);
}

void QuadMesh::visualize() {
    // smoothest field used to compute singularities and branch cover
    polyscope::getSurfaceMesh("pre-uniformization")->addVectorQuantity("Smoothest Cross Field", GF.field, 4);

    // visualize singularities
    std::vector<Vector3> singularCloud_pos;
    std::vector<Vector3> singularCloud_neg;
    std::vector<Vector3> singularCloud_other;
    for (VertexPtr v : mesh->vertices()) {
        if (singularities[v] == 1) {
            singularCloud_pos.push_back(geom->position(v));
        } else if (singularities[v] == -1) {
            singularCloud_neg.push_back(geom->position(v));
        } else if (singularities[v] != 0) {
            singularCloud_other.push_back(geom->position(v));
        }
    }
    polyscope::registerPointCloud("+1 singularities", singularCloud_pos);
    polyscope::registerPointCloud("-1 singularities", singularCloud_neg);
    polyscope::registerPointCloud("other singularities", singularCloud_other);

    polyscope::registerSurfaceMesh("post-uniformization", geom);

    // curvatures post-uniformization
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("Curvatures", U.cmCurvatures);
 
    // updated cross field on cone metric
    polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("Smoothest Cross Field CM", GF.field, 4);

    // cross frame on branch cover
    for (int i = 0; i < 4; i++) {
      polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("Branch Cover Field " + std::to_string(i), branchCoverFields[i], 1);
    }

    // errors after bfs
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("edge error", errors);

    polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("omega 0", omega[0]);
    polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("omega 1", omega[1]);
    polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("omega 2", omega[2]);
    polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("omega 3", omega[3]);
    
    // direct stripe coords
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("X Coords", coords[0]);
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("Y Coords", coords[1]);

    FaceData<std::pair<int,std::vector<Vector2>>> shaderInfo(mesh);
    for (FacePtr f : mesh->faces()) {
        std::vector<Vector2> coords = texCoords[f];
        int n = zeros[f];
        shaderInfo[f] = std::make_pair(n, coords);
    }
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("stripes shader", shaderInfo);

    polyscope::getSurfaceMesh("post-uniformization")->getSurfaceQuantity("stripes shader")->enabled = true;  
    polyscope::getSurfaceMesh("post-uniformization")->setActiveSurfaceQuantity((polyscope::SurfaceQuantityThatDrawsFaces*)polyscope::getSurfaceMesh("post-uniformization")->getSurfaceQuantity("stripes shader"));

    VertexData<double> psi_real(mesh,0);
    VertexData<double> psi_imag(mesh,0);
    for (VertexPtr v : mesh->vertices()){
        psi_real[v] = psi[0][v].real();
        psi_imag[v] = psi[0][v].imag();
    }
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("Psi real", psi_real);
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("Psi imag", psi_imag);

    VertexData<double> sigma_sums(mesh,0);
    for (VertexPtr v : mesh->vertices()) {  
        if (singularities[v] == 0) continue;
        BHalfedge BHe_first = BHalfedge{v.halfedge(), 0, &BC};
        BHalfedge BHe = BHe_first;
        double sigma_sum = 0;
        size_t iter = 0;
        size_t d = v.degree();
        for (size_t iter = 0; iter < 2*d; iter++) {
            int sign_ij = (BHe.next().edge().halfedge() == BHe.next()) ? 1 : -1;
            double sigma_ij = sign_ij * getSigma(BHe.next().edge());
            sigma_sum += sigma_ij;
            //std::cout << "sigma_ij: " << sigma_ij << std::endl; 
            BHe = BHe.next().next().twin();
        }
        //std::cout << "sigma_sum: " << sigma_sum << std::endl;
        //std::cout << "-----------" << std::endl;
        sigma_sums[v] = sigma_sum;
    }
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("sigma sum 0", sigma_sums);

}