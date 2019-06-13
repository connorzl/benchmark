#include "polyscope/quad_mesh.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <queue>

QuadMesh::QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), theta(m), r(m), field(m),
                                                                singularities(m), eta(m),
                                                                edgeLengthsCM(m), thetaCM(m), rCM(m), fieldCM(m),
                                                                cmAngles(m), cmAreas(m), curvatures(m), texCoords(m) {
    for (int i = 0; i < n; i++) {
        FaceData<std::complex<double>> sheetField(mesh,0);
        VertexData<size_t> sheetIndices(mesh);
        branchCoverFields.push_back(sheetField);
        BVertexIndices.push_back(sheetIndices);
    }

    for (int i = 0; i < 2; i++) {
        VertexData<double> sheetCoords(mesh);
        VertexData<std::complex<double>> sheetPsi(mesh,0);
        EdgeData<double> sheetSigma(mesh);

        coords.push_back(sheetCoords);
        psi.push_back(sheetPsi);
        sigma.push_back(sheetSigma);
    }
}

void QuadMesh::setup() {
    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_jl = geom->length(he.edge());
        double l_ij = geom->length(he.prev().edge());
        double theta_ijl = geom->angle(he.next()); // radians

        // Place first vertex at (0,0)
        Vector2 j = Vector2{0,0};
        Vector2 l = Vector2{l_jl,0};
        Vector2 i = Vector2{cos(theta_ijl) * l_ij, sin(theta_ijl) * l_ij}; 

        Vector2 jl = l - j;
        Vector2 li = i - l;
        Vector2 ij = j - i;

        theta[he] = std::complex<double>(jl.x,jl.y);
        theta[he.next()] = std::complex<double>(li.x,li.y);
        theta[he.prev()] = std::complex<double>(ij.x,ij.y);
    }
    
    // Compute d_ij * r_ji = -d_ji
    for (VertexPtr v : mesh->vertices()) {
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) {
                continue;
            }
            std::complex<double> theta_ij = theta[he];
            std::complex<double> theta_ji = theta[he.twin()];
            r[he] = std::pow((-theta_ij / theta_ji), n);
            r[he] = r[he] / std::abs(r[he]);
        }
    }   
}

void QuadMesh::setupCM() {
    // make sure that the CM edge lengths satisfy what we want for curvature
    Eigen::MatrixXd K  = Operators::intrinsicCurvature(mesh, edgeLengthsCM);
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    for (VertexPtr v : mesh->vertices()) {
        if (v.isBoundary()) continue;
        double C = K(vertexIndices[v],0);
        if (singularities[v] == 0 && std::abs(C) > 1e-8) {
            std::cout << "wrong curvature on non-singular vertex" << std::endl;
        } else if (singularities[v] == 1 && std::abs(C - M_PI_2) > 1e-8) {
            std::cout << "wrong curvature on +1 singular vertex" << std::endl;
        } else if (singularities[v] == -1 && std::abs(C + M_PI_2) > 1e-8) {
            std::cout << "wrong curvature on -1 singular vertex" << std::endl;
        }
    }

    for (FacePtr f : mesh->faces()) {
        // Gather elements
        HalfedgePtr he = f.halfedge();
        double l_ij = edgeLengthsCM[he.edge()];
        double l_ki = edgeLengthsCM[he.prev().edge()];
        double theta_ijk = cmAngles[he.next()];

        // Place first vertex at (0,0)
        Vector2 i = Vector2{0,0};
        Vector2 j = Vector2{l_ij,0};
        Vector2 k = Vector2{cos(theta_ijk) * l_ki, sin(theta_ijk) * l_ki}; 

        Vector2 ij = j - i;
        Vector2 jk = k - j;
        Vector2 ki = i - k;

        thetaCM[he] = std::complex<double>(ij.x,ij.y);
        thetaCM[he.next()] = std::complex<double>(jk.x,jk.y);
        thetaCM[he.prev()] = std::complex<double>(ki.x,ki.y);
    }
    
    // Compute d_ji * r_ij = -d_ij
    for (VertexPtr v : mesh->vertices()) {
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) continue;
            std::complex<double> theta_ij = thetaCM[he];
            std::complex<double> theta_ji = thetaCM[he.twin()];
            rCM[he] = std::pow((-theta_ij / theta_ji), n);
            rCM[he] /= std::abs(rCM[he]);
        }
    }
}

Eigen::SparseMatrix<std::complex<double>> QuadMesh::assembleM(bool isCM) {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> M(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];
        if (isCM) {
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, cmAreas[f]));
        } else {
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, geom->area(f)));
        }
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

Eigen::SparseMatrix<std::complex<double>> QuadMesh::assembleA(bool isCM) {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];
        HalfedgePtr he_ij = f.halfedge();

        std::complex<double> r_ij, r_jk, r_ki;
        if (isCM) {
            r_ij = rCM[he_ij];
            r_jk = rCM[he_ij.next()];
            r_ki = rCM[he_ij.prev()];
        } else {
            r_ij = r[he_ij];
            r_jk = r[he_ij.next()];
            r_ki = r[he_ij.prev()];
        }

        double numReal = 0;
        HalfedgePtr he = f.halfedge();
        do {
            if (!he.edge().isBoundary()) {
                numReal++;
            }
            he = he.next();
        } while (he != f.halfedge());

        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, numReal));
        if (!he_ij.edge().isBoundary()) {
            size_t i_ij = faceIndices[he_ij.twin().face()];
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ij,-r_ij));
        }
        if (!he_ij.next().edge().isBoundary()) {
            size_t i_jk = faceIndices[he_ij.next().twin().face()];
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_jk,-r_jk));
        }
        if (!he_ij.prev().edge().isBoundary()) {
            size_t i_ki = faceIndices[he_ij.prev().twin().face()];
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ki,-r_ki));   
        }
    }
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

void QuadMesh::computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A, bool isCM) {
    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd u = Eigen::MatrixXcd::Random(mesh->nFaces(),1);
    Eigen::MatrixXcd x;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    PositiveDefiniteSolver<std::complex<double>> s(A);
    for (int i = 0; i < nPowerIterations; i++) {
        Eigen::MatrixXcd rhs = M * u;
        x = s.solve(rhs);

        std::complex<double> norm2 = (x.transpose() * M * x)(0,0);
        u = x / sqrt(norm2);
    }

    // map resulting vector to VertexData
    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        std::complex<double> c = u(faceIndices[f],0);
        std::complex<double> val;
        if (std::abs(c) == 0) {
            val = 0;
        } else {
            val = c / std::abs(c);   
        }
        if (isCM) {
            fieldCM[f] = val;
        } else {
            field[f] = val;
        }
    }
} 

void QuadMesh::computeCrossField(bool isCM) {
    std::cout << "Computing Smoothest Cross Field" << std::endl;
    if (isCM) {
        setupCM();
    } else {
        setup();
    }

    // Algorithm 2 : Smoothest Field
    Eigen::SparseMatrix<std::complex<double>> M = assembleM(isCM);
    Eigen::SparseMatrix<std::complex<double>> A = assembleA(isCM);
    A = A + eps * M;
    computeSmoothestField(M,A,isCM);
    std::cout << "Done!" << std::endl;
}

void QuadMesh::computeSingularities() {
    std::cout << "Computing Singularities...";

    // compute index for each vertex v
    int total = 0;
    for (VertexPtr v : mesh->vertices()) {
        double angleSum = 0;

        if (v.isBoundary()) {
            singularities[v] = 0;
            continue;
        } 

        for (HalfedgePtr he : v.outgoingHalfedges()) {
            std::complex<double> u_i = field[he.face()];
            std::complex<double> u_j = field[he.twin().face()];
            std::complex<double> r_ji = r[he];
            angleSum += std::arg(u_i / (r_ji * u_j));
        }

        double phi = (angleSum + n*geom->angleDefect(v)) / (2.0 * M_PI);
        singularities[v] = std::round(phi);
        if (singularities[v] != 0) {
            numSingularities++;
        }
        total += singularities[v];
    }
    std::cout << "Done! Singularities Index Sum: " << total << std::endl;
}

double QuadMesh::updateAreas() {
    double totalArea = 0;
    for (FacePtr f : mesh->faces()) {
        double l_ij = edgeLengthsCM[f.halfedge().edge()       ];
        double l_jk = edgeLengthsCM[f.halfedge().next().edge()];
        double l_ki = edgeLengthsCM[f.halfedge().prev().edge()];

        double s = (l_ij + l_jk + l_ki) / 2.0;
        cmAreas[f] = sqrt( s * (s - l_ij) * (s - l_jk) * (s - l_ki) );
        totalArea += cmAreas[f];
    }
    return totalArea;
}

// this needs to include hyperbolic edge flips
void QuadMesh::uniformizeBoundary() {
    std::cout << "Boundary Uniformization..." << std::endl;
    // re-index the interior vertices starting from 0
    // we want to omit boundary vertices from the energy and pin u to 0 
    size_t iN = 0;
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    VertexData<size_t> interiorVertexIndices(mesh);
    Eigen::Array<bool,Eigen::Dynamic,1> interior(mesh->nVertices(),1);

    for (VertexPtr v : mesh->vertices()) {
        if (!v.isBoundary()) {
            interiorVertexIndices[v] = iN++;
        }
        interior(vertexIndices[v],0) = (!v.isBoundary());
    }
    if (iN != mesh->nInteriorVertices()) {
        throw std::logic_error("error indexing interior vertices!");
    }
    
    // Ax = b, where A = |V| x |V| laplacian, x = |V| x 1 vertex scaling, b = |V| x 1 curvature diff K - K*
    Eigen::MatrixXd KTarg(mesh->nInteriorVertices(),1);
    Eigen::MatrixXd K(mesh->nInteriorVertices(),1);
    for (VertexPtr v : mesh->vertices()) {
        if (v.isBoundary()) continue;
        size_t index = interiorVertexIndices[v];
        if (singularities[v] == 0) {
            KTarg(index,0) = 0;
        } else if (singularities[v] == 1) {
            KTarg(index,0) = M_PI_2;
        } else if (singularities[v] == -1) {
            KTarg(index,0) = -M_PI_2;
        } else if (singularities[v] == -2) {
            KTarg(index,0) = -M_PI;
        } else {
            throw std::logic_error("non +1/-1/-2 singularity");
        }
    }

    // make mesh delaunay
    geom->getEdgeLengths(edgeLengthsCM);
    Operators::hyperbolicEdgeFlips(mesh, edgeLengthsCM);
    
    double resid = 0;
    do {
        Eigen::SparseMatrix<double> A_all = Operators::intrinsicLaplaceMatrix(mesh,edgeLengthsCM);
        BlockDecompositionResult<double> B = blockDecomposeSquare(A_all, interior);

        Eigen::MatrixXd K_all = Operators::intrinsicCurvature(mesh,edgeLengthsCM);
        for (VertexPtr v : mesh->vertices()) {
            if (v.isBoundary()) continue;
            K(interiorVertexIndices[v],0) = K_all(vertexIndices[v],0);
        }

        Eigen::MatrixXd rhs = KTarg - K;
        Eigen::MatrixXd x = solveSquare<double>(B.AA, rhs);

        // update edge lengths
        for (EdgePtr e : mesh->edges()) {
            VertexPtr vi = e.halfedge().vertex();
            VertexPtr vj = e.halfedge().twin().vertex();

            double ui, uj;
            if (vi.isBoundary()) {
                ui = 0;
            } else {
                ui = x(interiorVertexIndices[vi],0);
            }
            if (vj.isBoundary()) {
                uj = 0;
            } else {
                uj = x(interiorVertexIndices[vj],0);
            }
            double s = std::exp( (ui + uj) / 2 );
            edgeLengthsCM[e] *= s;
        }
        
        Operators::hyperbolicEdgeFlips(mesh, edgeLengthsCM);

        // check to see if triangle inequality still holds
        for (FacePtr f : mesh->faces()) {
            HalfedgePtr he = f.halfedge();
            double a = edgeLengthsCM[he.edge()];
            double b = edgeLengthsCM[he.next().edge()];
            double c = edgeLengthsCM[he.prev().edge()];
            if (a > b + c || b > a + c || c > a + b) {
                throw std::logic_error("Triangle Inequality Violated during Uniformization!");
            }
        }

        K_all = Operators::intrinsicCurvature(mesh,edgeLengthsCM);
        for (VertexPtr v : mesh->vertices()) {
            if (v.isBoundary()) continue;
            K(interiorVertexIndices[v],0) = K_all(vertexIndices[v],0);
        }
        resid = (KTarg - K).array().abs().maxCoeff();
        std::cout << "Resid: " << resid << std::endl;  
    } while (resid > 1e-12);

    K = Operators::intrinsicCurvature(mesh, edgeLengthsCM);
    // store curvatures for visualization
    for (VertexPtr v : mesh->vertices()) {
        curvatures[v] = K(vertexIndices[v],0);        
    }

    // update CM areas and angles
    updateAreas();
    cmAngles = Operators::computeAngles(mesh, edgeLengthsCM);
    std::cout << "Done!" << std::endl;
}

void QuadMesh::uniformize() {
    if (mesh->nBoundaryLoops() > 0) {
        //throw std::logic_error ("not implemented for boundary meshes yet");
        uniformizeBoundary();
        return;
    }
    if (mesh->eulerCharacteristic() != 2) {
        throw std::logic_error("non genus-0 shape");
    }
    std::cout << "Boundary-less Uniformization..." << std::endl;
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();

    // Ax = b, where A = |V| x |V| laplacian, x = |V| x 1 vertex scaling, b = |V| x 1 curvature diff K - K*
    Eigen::MatrixXd KTarg(mesh->nVertices(),1);
    for (VertexPtr v : mesh->vertices()) {
        size_t index = vertexIndices[v];
        if (singularities[v] == 0) {
            KTarg(index,0) = 0;
        } else if (singularities[v] == 1) {
            KTarg(index,0) = M_PI_2;
        } else if (singularities[v] == -1) {
            KTarg(index,0) = -M_PI_2;
        } else if (singularities[v] == -2) {
            KTarg(index,0) = -M_PI;
        } else {
            polyscope::warning("WEIRD SINGULARITY INDEX");
            std::cout << "warning: singularity index " << singularities[v] << std::endl;
            if (singularities[v] == 2) {
                KTarg(index,0) = M_PI;
                continue;
            }
            std::cout << singularities[v] << std::endl;
            throw std::logic_error("non +1/-1/-2 singular index");
        }
    }
    if (std::abs(KTarg.sum() - mesh->eulerCharacteristic() * 2.0 * M_PI) > 1e-8) {
        throw std::logic_error("Target curvatures do not satisfy Gauss Bonnet");
    }

    // make mesh delaunay
    geom->getEdgeLengths(edgeLengthsCM);    
    Operators::hyperbolicEdgeFlips(mesh, edgeLengthsCM);

    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd x, K;
    double resid = 0;
    do {
        A = Operators::intrinsicLaplaceMatrix(mesh,edgeLengthsCM);
        K = Operators::intrinsicCurvature(mesh,edgeLengthsCM);

        Eigen::MatrixXd rhs = KTarg - K;
        x = solveSquare<double>(A, rhs);
        x = x.array() - x.mean();

        // update edge lengths
        for (EdgePtr e : mesh->edges()) {
            VertexPtr vi = e.halfedge().vertex();
            VertexPtr vj = e.halfedge().twin().vertex();

            double ui = x(vertexIndices[vi],0);
            double uj = x(vertexIndices[vj],0);
            double s = std::exp( (ui + uj) / 2 );
            edgeLengthsCM[e] *= s;
        }

        // perform hyperbolic edge flips to make sure things are okay
        Operators::hyperbolicEdgeFlips(mesh, edgeLengthsCM);

        // check to see if triangle inequality still holds
        for (FacePtr f : mesh->faces()) {
            HalfedgePtr he = f.halfedge();
            double a = edgeLengthsCM[he.edge()];
            double b = edgeLengthsCM[he.next().edge()];
            double c = edgeLengthsCM[he.prev().edge()];
            if (a > b + c || b > a + c || c > a + b) {
                throw std::logic_error("Triangle Inequality Violated during Uniformization!");
            }
        }

        K = Operators::intrinsicCurvature(mesh,edgeLengthsCM);
        resid = (KTarg - K).array().abs().maxCoeff();
        std::cout << "Resid: " << resid << std::endl;  
    } while (resid > 1e-12);

    // store curvatures for visualization
    for (VertexPtr v : mesh->vertices()) {
        curvatures[v] = K(vertexIndices[v],0);        
    }

    // update CM areas and angles
    updateAreas();
    cmAngles = Operators::computeAngles(mesh, edgeLengthsCM);
    std::cout << "Done!" << std::endl;
}

void QuadMesh::computeBranchCover(bool improve) {
    computeCrossField(true);

    std::cout<< "Computing Branch Cover...";
    std::complex<double> i(0, 1);
    for (VertexPtr v : mesh->vertices()) {
        int total = 0;
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) {
                eta[he] = 0;
                continue;
            }
            std::complex<double> u_i = std::pow(fieldCM[he.face()], 1.0 / n);
            std::complex<double> u_j = std::pow(fieldCM[he.twin().face()], 1.0 / n);
            
            std::complex<double> theta_ij = thetaCM[he];
            std::complex<double> theta_ji = thetaCM[he.twin()];
            std::complex<double> r_ji = -theta_ij / theta_ji;
            r_ji = r_ji / std::abs(r_ji);
            double ang = std::arg(u_i / (r_ji * u_j));

            if (ang >= -M_PI_4 && ang < M_PI_4) {
                eta[he] = 0;
            } else if (ang >= M_PI_4 && ang < 3.0 * M_PI_4) {
                eta[he] = 1;
                total = (total + 1) % n;
            } else if ((ang >= 3.0 * M_PI_4 && ang <= PI) || 
                       (ang < -3.0 * M_PI_4 && ang >= -PI)) {
                eta[he] = 2;
                total = (total + 2) % n;
            } else {
                if (!(ang >= -3.0 * M_PI_4 && ang < -M_PI_4)) {
                    throw std::logic_error("angle not in right range for branch cover");
                }
                eta[he] = 3;
                total = (total + 3) % n;
            }
        }
        if (v.isBoundary()) continue;
        if (singularities[v] != 0 && total == 0) {
            std::cout << "difference at singularity: " << total << std::endl;
        } else if (singularities[v] == 0 && total != 0) {
            std::cout << "difference at non-singularity: " << total << std::endl; 
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
                    int offset;
                    if (sheetDiff == 0) {  
                        offset = 0;
                    } else if (sheetDiff == 1) {
                        offset = -1;
                    } else if (sheetDiff == 2) {
                        offset = 2;
                    } else if (sheetDiff == 3) {
                        offset = 1;
                    } else {
                        throw std::logic_error("impossible offset");
                    }
                    for (HalfedgePtr Nhe : neighbor.adjacentHalfedges()) {
                        int oldVal = eta[Nhe];
                        eta[Nhe] = (eta[Nhe] + offset) % n;
                        if (Nhe.twin().isReal()) {
                            eta[Nhe.twin()] = (eta[Nhe.twin()] - offset) % n;
                        }
                        if (eta[Nhe] == -1) eta[Nhe] = 3;
                        if (eta[Nhe.twin()] == -1) eta[Nhe.twin()] = 3;
                        if (eta[Nhe] == -3) eta[Nhe] = 1;
                        if (eta[Nhe.twin()] == -3) eta[Nhe.twin()] = 1;
                        if (eta[Nhe] == -2) eta[Nhe] = 2;
                        if (eta[Nhe.twin()] == -2) eta[Nhe.twin()] = 2;
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
            std::complex<double> theta_ij = thetaCM[BHe.he];
            std::complex<double> theta_ji = thetaCM[BHe.he.twin()];
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

void QuadMesh::computeCrossFieldCMBranchCover(std::complex<double> init, double scale) {
    std::cout << "Computing Cross Field CM on Branch Cover..." << std::endl;
    
    // get faces for traversing branch cover
    std::vector<BFace> allBFaces = BC.allFaces();
    for (int i = 0; i < n; i++) {
        FaceData<std::complex<double>> sheetField(mesh,0);
        branchCoverFields[i] = sheetField;

        FaceData<std::complex<double>> sheetXBasis(mesh);
        xBasis.push_back(sheetXBasis);
    }

    // perform BFS starting on some BFace
    BFace root = allBFaces[0];
    std::map<BFace,bool> visited;
    visited[root] = true;
    branchCoverFields[root.sheet][root.f] = init;
    xBasis[root.sheet][root.f] = std::complex<double> (1,0);

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

                std::complex<double> theta_ij = thetaCM[BHe.he];
                std::complex<double> theta_ji = thetaCM[BHe.he.twin()];
                std::complex<double> r_ij = -theta_ji / theta_ij;
                r_ij = r_ij / std::abs(r_ij);

                branchCoverFields[neighbor.sheet][neighbor.f] = r_ij * branchCoverFields[currFace.sheet][currFace.f];
                xBasis[neighbor.sheet][neighbor.f] = r_ij * xBasis[currFace.sheet][currFace.f];
                count++;
            } 
            BHe = BHe.next();
        } while (BHe != currFace.halfedge());
        
        // for the torus
        BFace neighbor = currFace;
        neighbor.sheet = (neighbor.sheet + 1) % 4;
        if (visited.find(neighbor) == visited.end()) {
            Q.push(neighbor);
            visited[neighbor] = true;
            branchCoverFields[neighbor.sheet][neighbor.f] = IM_I * branchCoverFields[currFace.sheet][currFace.f];
            xBasis[neighbor.sheet][neighbor.f] = IM_I * xBasis[currFace.sheet][currFace.f];
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

        std::complex<double> theta_ij = thetaCM[BHe.he];
        std::complex<double> theta_ji = thetaCM[BHe.he.twin()];
        std::complex<double> r_ji = -theta_ij / theta_ji;
        r_ji = r_ji / std::abs(r_ji);

        std::complex<double> f_ij = branchCoverFields[Bf_ij.sheet][Bf_ij.f];
        std::complex<double> f_ijt = r_ji * branchCoverFields[Bf_ji.sheet][Bf_ji.f];
        errors[Be.e] = std::max(errors[Be.e], std::abs(f_ijt - f_ij));
    }
    // compute omega using field
    computeOmega(scale);
    std::cout << "Done!" << std::endl;
}

void QuadMesh::computeOmega(double scale) {
    // set up an edgedata per sheet
    for (int i = 0; i < n; i++) {
        EdgeData<double> omegaSheet(mesh);
        omega.push_back(omegaSheet);
    }

    // compute omega on each edge of the branch cover
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BHalfedge he_ij = Be.halfedge();
        BHalfedge he_ji = he_ij.twin();

        std::complex<double> f_ij, f_ji, e_ij, e_ji;
        
        double boundaryScale = 1.0;
        if (he_ij.he.isReal()){
            f_ij = scale * branchCoverFields[he_ij.sheet][he_ij.face().f];
            e_ij = thetaCM[he_ij.he];
        } else {
            f_ij = 0;
            e_ij = 0;
            boundaryScale = 2.0;
        }
        if (he_ji.he.isReal()) {
            f_ji = scale * branchCoverFields[he_ji.sheet][he_ji.face().f];
            e_ji = -thetaCM[he_ji.he];
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
}

Eigen::SparseMatrix<double> QuadMesh::energyMatrix() {
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

    // build cotan component
    std::vector<BFace> allBFaces = BC.allFaces();
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
        double cot_ij = 1.0 / tan(cmAngles[BHe_ij.he]);
        double cot_jk = 1.0 / tan(cmAngles[BHe_jk.he]);
        double cot_ki = 1.0 / tan(cmAngles[BHe_ki.he]);

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

        // jj, ji, jk
        tripletsCot.push_back( Eigen::Triplet<double>(j_re,j_re, (cot_jk + cot_ij) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j_im,j_im, (cot_jk + cot_ij) / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(j_re,i_re, -cot_ij / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j_im,i_im, sj_im * si_im * -cot_ij / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(j_re,k_re, -cot_jk / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(j_im,k_im, sj_im * sk_im * -cot_jk / 2.0) );

        // kk, ki, kj
        tripletsCot.push_back( Eigen::Triplet<double>(k_re,k_re, (cot_jk + cot_ki) / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k_im,k_im, (cot_jk + cot_ki) / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(k_re,i_re, -cot_ki / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k_im,i_im, sk_im * si_im * -cot_ki / 2.0) );

        tripletsCot.push_back( Eigen::Triplet<double>(k_re,j_re, -cot_jk / 2.0) );
        tripletsCot.push_back( Eigen::Triplet<double>(k_im,j_im, sk_im * sj_im * -cot_jk / 2.0) );
    }
    Acot.setFromTriplets(tripletsCot.begin(),tripletsCot.end());

    // build zArea component
    double scale = 100;
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
        double zArea = scale * scale * cmAreas[Bf.f];

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

        // jj, ji, jk
        tripletszArea.push_back( Eigen::Triplet<double>(j_re,j_re, zArea / 6.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(j_im,j_im, zArea / 6.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(j_re,i_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(j_im,i_im, sj_im * si_im * zArea / 12.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(j_re,k_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(j_im,k_im, sj_im * sk_im * zArea / 12.0) );

        // kk, ki, kj
        tripletszArea.push_back( Eigen::Triplet<double>(k_re,k_re, zArea / 6.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(k_im,k_im, zArea / 6.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(k_re,i_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(k_im,i_im, sk_im * si_im * zArea / 12.0) );

        tripletszArea.push_back( Eigen::Triplet<double>(k_re,j_re, zArea / 12.0) );
        tripletszArea.push_back( Eigen::Triplet<double>(k_im,j_im, sk_im * sj_im * zArea / 12.0) );
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
        std::complex<double> z = scale * branchCoverFields[Bf.sheet][Bf.f];
        std::complex<double> e_i_perp = IM_I * thetaCM[BHe_jk.he];
        std::complex<double> e_j_perp = IM_I * thetaCM[BHe_ki.he];
        std::complex<double> e_k_perp = IM_I * thetaCM[BHe_ij.he];
       
        double drift_ij = -(1.0 / 6.0) * dot(e_i_perp - e_j_perp, z); 
        double drift_jk = -(1.0 / 6.0) * dot(e_j_perp - e_k_perp, z); 
        double drift_ki = -(1.0 / 6.0) * dot(e_k_perp - e_i_perp, z); 

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
}

Eigen::SparseMatrix<double> QuadMesh::energyMatrix2() {
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
            cotA = 1.0 / tan(cmAngles[BHe.he]);
            A = BHe.he.prev().vertex();
        } else {
            cotA = 0;
        }
        if (BHe.he.twin().isReal()) {
            cotB = 1.0 / tan(cmAngles[BHe.twin().he]);
            B = BHe.twin().he.prev().vertex();
        } else {
            cotB = 0;
        }

        // check for singularities
        if (BHe.he.isReal() && singularities[A] != 0 && 
            BHe.he.twin().isReal() && singularities[B] != 0) {
            continue;
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

        double area = cmAreas[Bf.f]; 
        
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

double QuadMesh::computeStripes() {
    std::cout << "Computing Stripes..." << std::endl;
    
    // build matrices for inverse power method
    Eigen::SparseMatrix<double> A = energyMatrix();
    Eigen::SparseMatrix<double> B = massMatrix();

    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular); // + numSingularities);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(numPsi,1);
    Eigen::MatrixXd prevX;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    PositiveDefiniteSolver<double> s(A);
    for (int i = 0; i < nPowerIterations; i++) {
        prevX = x;
        
        Eigen::MatrixXd rhs = B * x;
        x = s.solve(rhs);

        double norm2 = (x.transpose() * B * x)(0,0);
        x = x / sqrt(norm2);

        Eigen::MatrixXd resid = x - prevX;
        std::cout << "Resid: " << resid.transpose() * B * resid  << std::endl;
    }
    std::cout << "Done!" << std::endl;

    Eigen::MatrixXd num = x.transpose() * A * x;
    Eigen::MatrixXd denom = x.transpose() * B * x;
    double lambda = num(0,0) / denom(0,0);

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
    return lambda;
}

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

// need to do some extra work for meshes with boundary
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
            newOmega(i,0) = omega[Be.sheet][Be.e];
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
                                              BC, cmAngles, cmAreas, mesh, numSingularities);
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

        double y = ( w_ij + std::arg(psi_i) - std::arg(psi_j) ) / (2 * PI);
        int n_ij = std::round(y);
        sigma[Be.sheet][Be.e] = std::arg(psi_j) - std::arg(psi_i) + 2 * PI * n_ij;
    }
}

bool QuadMesh::textureCoordinates() {
    std::cout << "Computing texture coordinates...";
    FaceData<std::vector<double>> xCoords(mesh);
    FaceData<std::vector<double>> yCoords(mesh);
    FaceData<int> singularFace(mesh,0);

    computeSigma();
    zeros = FaceData<int>(mesh,0);
    int numNewSingularities = 0;

    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        // sheet 0 = x coordinate, sheet 1 = y coordinate
        if (Bf.sheet != 0 && Bf.sheet != 1) continue;

        // grab halfedges, vertices, and edges
        BHalfedge Bhe_ij = Bf.halfedge();
        BHalfedge Bhe_jk = Bhe_ij.next();
        BHalfedge Bhe_ki = Bhe_jk.next();

        BVertex Bv_i = Bhe_ij.vertex();
        BVertex Bv_j = Bhe_jk.vertex();
        BVertex Bv_k = Bhe_ki.vertex();

        BEdge Be_ij = Bhe_ij.edge();
        BEdge Be_jk = Bhe_jk.edge();
        BEdge Be_ki = Bhe_ki.edge();

        // skip singularities for now, since they have no psi associated with them
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || singularities[Bv_k.v] != 0) {
            singularFace[Bf.f] = 1;
            continue;
        }

        // get orientations of edges
        int c_ij = (Be_ij.halfedge() == Bhe_ij) ? 1 : -1;
        int c_jk = (Be_jk.halfedge() == Bhe_jk) ? 1 : -1;
        int c_ki = (Be_ki.halfedge() == Bhe_ki) ? 1 : -1;

        // retrieve sigmas
        double sigma_ij = c_ij * getSigma(Be_ij);
        double sigma_jk = c_jk * getSigma(Be_jk);
        double sigma_ki = c_ki * getSigma(Be_ki);

        // get the vectors at each vertex
        std::complex<double> z_i = getPsi(Bv_i);
        std::complex<double> z_j = getPsi(Bv_j);
        std::complex<double> z_k = getPsi(Bv_k);
        
        // compute coordinates
        std::complex<double> i(0,1);
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
            std::cout << std::arg(getPsi(Bv_i)) << "," << std::arg(getPsi(Bv_j)) << "," << std::arg(getPsi(Bv_k)) << std::endl;
            std::cout << omega[Be_ij.sheet][Be_ij.e] << "," << omega[Be_jk.sheet][Be_jk.e] << "," << omega[Be_ki.sheet][Be_ki.e] << std::endl; 
            std::cout << sigma_ij << "," << sigma_jk << "," << sigma_ki << std::endl;
            double y = ( omega[Be_ij.sheet][Be_ij.e] + std::arg(getPsi(Bv_i)) - std::arg(getPsi(Bv_j)) ) / (2 * PI);
            int n_ij = std::round(y);
            std::cout << c_ij * (std::arg(getPsi(Bv_j)) - std::arg(getPsi(Bv_i)) + 2 * PI * n_ij) << std::endl;
        }
        
        std::vector<double> currCoords = {c_i, c_j, c_k};
        if (Bf.sheet == 0) {
            xCoords[Bf.f] = currCoords;
        } else {
            yCoords[Bf.f] = currCoords;
        }
    }

    // prepare data for shader
    for (FacePtr f : mesh->faces()) {
        std::vector<double> xs = xCoords[f];
        std::vector<double> ys = yCoords[f];
        std::vector<Vector2> coords; 
        for (size_t i = 0; i < 3; i++) {
            Vector2 vertCoords;
            if (singularFace[f] != 0) {   
                //vertCoords = {.x = xs[i], .y = ys[i]};
                vertCoords = {0, 0}; // default value for singularities
            } else {
                vertCoords = {.x = xs[i], .y = ys[i]};
            }
            coords.push_back(vertCoords);
        } 
        texCoords[f] = coords;
    }
    std::cout << "Done!" << std::endl;
    std::cout << "Num new singularities: " << numNewSingularities << std::endl;
    return (numNewSingularities == 0);
}

void QuadMesh::visualize() {
    polyscope::registerSurfaceMesh("post-uniformization", geom);

    // smoothest field used to compute singularities and branch cover
    polyscope::getSurfaceMesh("pre-uniformization")->addVectorQuantity("Smoothest Cross Field", field, 4);
    polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("Smoothest Cross Field", field, 4);
    
    // singularities
    VertexData<Vector3> singularityColors(mesh);
    for (VertexPtr v : mesh->vertices()) {
      if (singularities[v] == 1) {
        singularityColors[v] = Vector3{1,0,0};
      } else if (singularities[v] == -1) {
        singularityColors[v] = Vector3{0,0,1};
      } else if (singularities[v] != 0) {
        singularityColors[v] = Vector3{0,1,0};
      } else {
        singularityColors[v] = Vector3{0.75,0.75,0.75};
      }
    }
    polyscope::getSurfaceMesh("post-uniformization")->addColorQuantity("Singularities", singularityColors);

    // curvatures post-uniformization
    polyscope::getSurfaceMesh("post-uniformization")->addQuantity("Curvatures", curvatures);

    // updated cross field on cone metric
    polyscope::getSurfaceMesh("pre-uniformization")->addVectorQuantity("Smoothest Cross Field CM", fieldCM, 4);
    polyscope::getSurfaceMesh("post-uniformization")->addVectorQuantity("Smoothest Cross Field CM", fieldCM, 4);
    
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

    FaceData<Vector3> zerosColors(mesh);
    for (FacePtr f : mesh->faces()) {
      if (zeros[f] == 1) {
        zerosColors[f] = Vector3{1,0,0};
      } else if (zeros[f] == -1) {
        zerosColors[f] = Vector3{0,0,1};
      } else {
        zerosColors[f] = Vector3{0,1.0,0};
      }
    }
    polyscope::getSurfaceMesh("post-uniformization")->addColorQuantity("New Zeros", zerosColors);  

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

}