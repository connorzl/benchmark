#include "polyscope/quad_mesh.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <queue>

QuadMesh::QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), theta(m), r(m), field(m), 
                                                                singularities(m), eta(m),
                                                                edgeLengthsCM(m), thetaCM(m), rCM(m), cmAngles(m),
                                                                cmAreas(m), curvatures(m), fieldCM(m), texCoords(m) {
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

Eigen::SparseMatrix<std::complex<double>> QuadMesh::assembleM() {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> M(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, geom->area(f)));
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

Eigen::SparseMatrix<std::complex<double>> QuadMesh::assembleA() {
    size_t n = mesh->nFaces();
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    FaceData<size_t> faceIndices = mesh->getFaceIndices();
    for (FacePtr f : mesh->faces()) {
        size_t i = faceIndices[f];
        HalfedgePtr he_ij = f.halfedge();

        std::complex<double> r_ij = r[he_ij];
        std::complex<double> r_jk = r[he_ij.next()];
        std::complex<double> r_ki = r[he_ij.prev()];

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

void QuadMesh::computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A) {
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
        /*
        // for testing the torus
        VertexPtr v_i = f.halfedge().vertex();
        VertexPtr v_j = f.halfedge().next().vertex();
        VertexPtr v_k = f.halfedge().prev().vertex();
        
        Vector3 p_i = geom->position(v_i);
        Vector3 p_j = geom->position(v_j);
        Vector3 p_k = geom->position(v_k);

        Vector3 centroid = (p_i + p_j + p_k) / 3.0;
        centroid = Vector3{-centroid.z, 0, centroid.x};

        Vector3 x = p_j - p_i;
        Vector3 normal = geom->normal(f);
        Vector3 y = cross(normal, x);

        double X1 = dot(x, centroid);
        double X2 = dot(y, centroid);
        std::complex<double> c(X1,X2);

        c = std::pow(c,4);
        */
        if (std::abs(c) == 0) {
            field[f] = 0;
        } else {
            field[f] = c / std::abs(c);   
        }
    }
} 

void QuadMesh::computeCrossField() {
    std::cout << "Computing Smoothest Cross Field..." << std::endl;
    // Algorithm 1 : Setup
    setup();

    // Algorithm 2 : Smoothest Field
    Eigen::SparseMatrix<std::complex<double>> M = assembleM();
    Eigen::SparseMatrix<std::complex<double>> A = assembleA();
    A = A + eps * M;
    computeSmoothestField(M,A);

    std::cout << "Done!" << std::endl;
}

void QuadMesh::computeSingularities() {
    std::cout << "Computing Singularities...";

    // finally, compute index for each vertex v
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

void QuadMesh::computeBranchCover(bool improve) {
    std::cout<< "Computing Branch Cover...";
    std::complex<double> i(0, 1);
    for (VertexPtr v : mesh->vertices()) {
        int total = 0;
        for (HalfedgePtr he : v.outgoingHalfedges()) {
            if (he.edge().isBoundary()) {
                eta[he] = 0;
                continue;
            }
            std::complex<double> u_i = std::pow(field[he.face()], 1.0 / n);
            std::complex<double> u_j = std::pow(field[he.twin().face()], 1.0 / n);
            
            std::complex<double> theta_ij = theta[he];
            std::complex<double> theta_ji = theta[he.twin()];
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
    std::cout << "Done!" << std::endl;
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
        } else {
            if (singularities[v] != -1) {
                throw std::logic_error("non +1/-1 singular index");
            }
            KTarg(index,0) = -M_PI_2;
        }
    }

    EdgeData<double> l0(mesh);
    geom->getEdgeLengths(l0);
    geom->getEdgeLengths(edgeLengthsCM);
    
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(mesh->nInteriorVertices(),1);
    Eigen::MatrixXd u_prev, x;
    int iter = 0;
    double resid = 0;
    
    do {
        u_prev = u;
        Eigen::SparseMatrix<double> A_all = Operators::intrinsicLaplaceMatrix(mesh,edgeLengthsCM);
        BlockDecompositionResult<double> B = blockDecomposeSquare(A_all, interior);

        Eigen::MatrixXd K_all = Operators::intrinsicCurvature(mesh,edgeLengthsCM);
        for (VertexPtr v : mesh->vertices()) {
            if (v.isBoundary()) continue;
            K(interiorVertexIndices[v],0) = K_all(vertexIndices[v],0);
        }

        Eigen::MatrixXd rhs = KTarg - K;
        x = solveSquare<double>(B.AA, rhs);
        u = u + x;

        // update edge lengths
        for (EdgePtr e : mesh->edges()) {
            VertexPtr vi = e.halfedge().vertex();
            VertexPtr vj = e.halfedge().twin().vertex();

            double ui, uj;
            if (vi.isBoundary()) {
                ui = 0;
            } else {
                ui = u(interiorVertexIndices[vi],0);
            }
            if (vj.isBoundary()) {
                uj = 0;
            } else {
                uj = u(interiorVertexIndices[vj],0);
            }
            double s = std::exp( (ui + uj) / 2 );
            edgeLengthsCM[e] = l0[e] * s;
        }
        
        resid = (u - u_prev).array().abs().maxCoeff();
        std::cout << "Norm of change: " << resid << std::endl;    
        
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

        iter++;
    } while (resid > 1e-12);

    K = Operators::intrinsicCurvature(mesh, edgeLengthsCM);
    // store curvatures for visualization
    for (VertexPtr v : mesh->vertices()) {
        curvatures[v] = K(vertexIndices[v],0);        
    }

    // update CM areas
    updateAreas();
    std::cout << "Done!" << std::endl;
}

void QuadMesh::uniformize() {
    if (mesh->nBoundaryLoops() > 0) {
        uniformizeBoundary();
        return;
    }

    std::cout << "Boundary-less Uniformization..." << std::endl;
    size_t n = mesh->nVertices();
    VertexData<size_t> vertexIndices = mesh->getVertexIndices();
    
    // Ax = b, where A = |V| x |V| laplacian, x = |V| x 1 vertex scaling, b = |V| x 1 curvature diff K - K*
    Eigen::MatrixXd KTarg(n,1);
    for (VertexPtr v : mesh->vertices()) {
        size_t index = vertexIndices[v];
        if (singularities[v] == 0) {
            KTarg(index,0) = 0;
        } else if (singularities[v] == 1) {
            KTarg(index,0) = M_PI_2;
        } else {
            if (singularities[v] != -1) {
                throw std::logic_error("non +1/-1 singular index");
            }
            KTarg(index,0) = -M_PI_2;
        }
    } 

    EdgeData<double> l0(mesh);
    geom->getEdgeLengths(l0);
    geom->getEdgeLengths(edgeLengthsCM);
    
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(n,1);
    Eigen::MatrixXd u_prev;
    Eigen::SparseMatrix<double> A;
    Eigen::MatrixXd K,x;
    int iter = 0;
    double resid = 0;

    do {
        u_prev = u;
        A = Operators::intrinsicLaplaceMatrix(mesh,edgeLengthsCM);
        K = Operators::intrinsicCurvature(mesh,edgeLengthsCM);

        Eigen::MatrixXd rhs = KTarg - K;
        x = solveSquare<double>(A, rhs);
        u = u + x;

        // store new areas for CM mass matrix
        double totalArea = updateAreas();
        VertexData<double> vertexAreas(mesh);
        for (VertexPtr v : mesh->vertices()) {
            double area = 0;
            for (FacePtr f : v.adjacentFaces()) {    // iterates around real faces
                area += cmAreas[f] / 3.0;
            }  
            vertexAreas[v] = area;
        }

        double u_avg = 0;
        for (VertexPtr v : mesh->vertices()){
            size_t index = vertexIndices[v];
            u_avg += vertexAreas[v] * u(index,0);
        }
        u_avg /= totalArea;
        u = u.array() - u_avg;

        // update edge lengths
        for (EdgePtr e : mesh->edges()) {
            VertexPtr vi = e.halfedge().vertex();
            VertexPtr vj = e.halfedge().twin().vertex();

            double ui = u(vertexIndices[vi],0);
            double uj = u(vertexIndices[vj],0);

            double s = std::exp( (ui + uj) / 2 );
            edgeLengthsCM[e] = l0[e] * s;
        }

        resid = (u - u_prev).array().abs().maxCoeff();
        std::cout << "Norm of change: " << resid << std::endl;    

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

        iter++;
    } while (resid > 1e-12);

    // store curvatures for visualization
    for (VertexPtr v : mesh->vertices()) {
        curvatures[v] = K(vertexIndices[v],0);        
    }

    // update CM areas
    updateAreas();
    std::cout << "Done!" << std::endl;
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

    cmAngles = Operators::computeAngles(mesh, edgeLengthsCM); 
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
            rCM[he] = -theta_ij / theta_ji;
            rCM[he] /= std::abs(rCM[he]);
        }
    }   

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
            vCurr = vCurr * rCM[BHe.he];
            BHe = BHe.twin().next();
        } while (BHe.he != firstHe.he);

        // check that the resulting vector makes sense
        std::complex<double> i(0,1);
        std::complex<double> rotPos90 = std::exp(i * M_PI_2);
        std::complex<double> rotNeg90 = std::exp(i * -M_PI_2);
        if (singularities[Bv.v] == 0 && std::norm(vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around non-singular vertex");
        } else if (singularities[Bv.v] == 1 && std::norm(rotNeg90 * vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around +1 singular vertex");
        } else if (singularities[Bv.v] == -1 && std::norm(rotPos90 * vCurr - v0) > 1e-8) {
            throw std::logic_error("wrong rotation around -1 singular vertex");
        }
    }
}

void QuadMesh::computeCrossFieldCMBranchCover(std::complex<double> init) {
    std::cout << "Computing Cross Field CM on Branch Cover..." << std::endl;
    
    // compute change of basis rotations
    setupCM();

    // get faces for traversing branch cover
    std::vector<BFace> allBFaces = BC.allFaces();

    for (int i = 0; i < n; i++) {
        FaceData<std::complex<double>> sheetField(mesh,0);
        branchCoverFields[i] = sheetField;
    }

    // perform BFS starting on some BFace
    BFace root = allBFaces[0];
    std::map<BFace,bool> visited;
    visited[root] = true;
    branchCoverFields[root.sheet][root.f] = init;  // this is a degree of freedom

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
                branchCoverFields[neighbor.sheet][neighbor.f] = rCM[BHe.twin().he] * branchCoverFields[currFace.sheet][currFace.f];
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
        std::complex<double> f_ij = branchCoverFields[Bf_ij.sheet][Bf_ij.f];
        std::complex<double> f_ijt = rCM[BHe.he] * branchCoverFields[Bf_ji.sheet][Bf_ji.f];
        errors[Be.e] = std::max(errors[Be.e], std::abs(f_ijt - f_ij));
    }
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

        // get cotan weights
        double cotA, cotB, w;
        VertexPtr A,B;
        // first check if boundary halfedges
        if (BHe.he.isReal()) {
            cotA = 1.0 / tan(cmAngles[BHe.he]);
            A = BHe.he.prev().vertex();
        } else {
            cotA = 0;
        }
        if (BHe.twin().he.isReal()) {
            cotB = 1.0 / tan(cmAngles[BHe.twin().he]);
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
        w = (cotA + cotB) / 2.0;

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
        triplets.push_back( Eigen::Triplet<double>(iA_re,iA_re,w) );
        triplets.push_back( Eigen::Triplet<double>(iA_im,iA_im,w) );

        // transport coefficient components
        BEdge Be = BHe.edge();
        double sign = (Be.halfedge() == BHe) ? 1 : -1;
        double omega_ij = sign * omega[Be.sheet][Be.e];
        double x = w * cos(omega_ij);
        double yA = sA_im * w * sin(omega_ij);
        double yB = sB_im * w * sin(omega_ij);
        double xAB = sA_im * sB_im * w * cos(omega_ij);

        // add non-diagonal entries
        triplets.push_back( Eigen::Triplet<double>(iA_re,iB_re,-x) ); triplets.push_back( Eigen::Triplet<double>(iA_re,iB_im,-yB ) );
        triplets.push_back( Eigen::Triplet<double>(iA_im,iB_re,yA) ); triplets.push_back( Eigen::Triplet<double>(iA_im,iB_im,-xAB) );
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
        triplets.push_back(Eigen::Triplet<double>(iA_re, iA_re, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iA_im, iA_im, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iB_re, iB_re, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iB_im, iB_im, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iC_re, iC_re, area/3.));
        triplets.push_back(Eigen::Triplet<double>(iC_im, iC_im, area/3.));
    }
    M.setFromTriplets(triplets.begin(),triplets.end());
    return M;
}

double QuadMesh::computeStripes(double scale) {
    std::cout << "Computing Stripes..." << std::endl;
    
    // compute 1-form and build matrices
    computeOmega(scale);
    Eigen::SparseMatrix<double> A = energyMatrix();
    Eigen::SparseMatrix<double> B = massMatrix();

    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular); // + numSingularities);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(numPsi,1);
    Eigen::MatrixXd prevX;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    SquareSolver<double> s(A);
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
            std::cout << "New Singularity: " << n << std::endl;
            //return false;
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
    return true;
}

void QuadMesh::visualize() {
    // smoothest field used to compute singularities and branch cover
    polyscope::getSurfaceMesh()->addVectorQuantity("Smoothest Cross Field", field, 4);
    
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
        singularityColors[v] = Vector3{0.9,0.9,0.9};
      }
    }
    polyscope::getSurfaceMesh()->addColorQuantity("Singularities", singularityColors);

    // curvatures post-uniformization
    polyscope::getSurfaceMesh()->addQuantity("Curvatures", curvatures);
    
    // cross frame on branch cover
    for (int i = 0; i < 4; i++) {
      polyscope::getSurfaceMesh()->addVectorQuantity("Branch Cover Field " + std::to_string(i), branchCoverFields[i], 1);
    }

    // errors after bfs
    polyscope::getSurfaceMesh()->addQuantity("edge error", errors);

    // direct stripe coords
    polyscope::getSurfaceMesh()->addQuantity("X Coords", coords[0]);
    polyscope::getSurfaceMesh()->addQuantity("Y Coords", coords[1]);

    FaceData<std::pair<int,std::vector<Vector2>>> shaderInfo(mesh);
    for (FacePtr f : mesh->faces()) {
        std::vector<Vector2> coords = texCoords[f];
        int n = zeros[f];
        shaderInfo[f] = std::make_pair(n, coords);
    }
    polyscope::getSurfaceMesh()->addQuantity("stripes shader", shaderInfo);

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
    polyscope::getSurfaceMesh()->addColorQuantity("New Zeros", zerosColors);  

    polyscope::getSurfaceMesh()->getSurfaceQuantity("stripes shader")->enabled = true;  
    polyscope::getSurfaceMesh()->setActiveSurfaceQuantity((polyscope::SurfaceQuantityThatDrawsFaces*)polyscope::getSurfaceMesh()->getSurfaceQuantity("stripes shader"));
    
}