#include "polyscope/quad_mesh.h"
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <queue>

// ONLY WORKS FOR MESHES WITHOUT BOUNDARY
QuadMesh::QuadMesh(HalfedgeMesh* m, Geometry<Euclidean>* g) : mesh(m), geom(g), theta(m), r(m), field(m), 
                                                                singularities(m), eta(m),
                                                                edgeLengthsCM(m), thetaCM(m), rCM(m), cmAngles(m),
                                                                cmAreas(m), curvatures(m), fieldCM(m), texCoords(m) {
    assert(mesh->nBoundaryLoops() == 0);

    for (int i = 0; i < n; i++) {
        FaceData<std::complex<double>> sheetField(mesh);
        VertexData<size_t> sheetIndices(mesh);
        branchCoverFields.push_back(sheetField);
        BVertexIndices.push_back(sheetIndices);
    }

    for (int i = 0; i < 2; i++) {
        VertexData<double> sheetCoords(mesh);
        VertexData<std::complex<double>> sheetPsi(mesh);
        EdgeData<double> sheetSigma(mesh);

        coords.push_back(sheetCoords);
        psi.push_back(sheetPsi);
        sigma.push_back(sheetSigma);
    }

    zeros = FaceData<int>(mesh,0);
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

        std::complex<double> r_ij = r[f.halfedge()];
        std::complex<double> r_jk = r[f.halfedge().next()];
        std::complex<double> r_ki = r[f.halfedge().prev()];

        size_t i_ij = faceIndices[f.halfedge().twin().face()];
        size_t i_jk = faceIndices[f.halfedge().next().twin().face()];
        size_t i_ki = faceIndices[f.halfedge().prev().twin().face()];

        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, i, 3));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ij,-r_ij));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_jk,-r_jk));
        triplets.push_back(Eigen::Triplet<std::complex<double>>(i,i_ki,-r_ki));
    }
    
    A.setFromTriplets(triplets.begin(),triplets.end());
    return A;
}

void QuadMesh::computeSmoothestField(Eigen::SparseMatrix<std::complex<double>> M, Eigen::SparseMatrix<std::complex<double>> A) {
    // u <- UniformRand(-1,1)
    Eigen::MatrixXcd u = Eigen::MatrixXcd::Random(mesh->nFaces(),1);
    Eigen::MatrixXcd x;

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    for (int i = 0; i < nPowerIterations; i++) {
        Eigen::MatrixXcd rhs = M * u;
        x = solvePositiveDefinite<std::complex<double>>(A, rhs);
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
    std::cout << "Computing Smoothest Cross Field...";
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
                assert(ang >= -3.0 * M_PI_4 && ang < -M_PI_4);
                eta[he] = 3;
                total = (total + 3) % n;
            }
        }
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
                        eta[Nhe.twin()] = (eta[Nhe.twin()] - offset) % n;
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

void QuadMesh::uniformize() {
    std::cout << "Uniformization..." << std::endl;
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
            assert(singularities[v] == -1);
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
            for (FacePtr f : v.adjacentFaces()) {
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

        resid = (u - u_prev).array().abs().maxCoeff();
        std::cout << "Norm of change: " << resid << std::endl;
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
            std::complex<double> theta_ij = thetaCM[he];
            std::complex<double> theta_ji = thetaCM[he.twin()];
            rCM[he] = -theta_ij / theta_ji;
            rCM[he] /= std::abs(rCM[he]);
        }
    }   

    // Sanity check the rotations here
    std::vector<BVertex> allVertices = BC.allVertices();
    for (BVertex Bv : allVertices) {
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

void QuadMesh::computeCrossFieldCMBranchCover() {
    std::cout << "Computing Cross Field CM on Branch Cover..." << std::endl;
    
    // compute change of basis rotations
    setupCM();

    // get faces for traversing branch cover
    std::vector<BFace> allBFaces = BC.allFaces();

    // perform BFS starting on some BFace
    BFace root = allBFaces[0];
    std::map<BFace,bool> visited;
    visited[root] = true;
    branchCoverFields[root.sheet][root.f] = std::complex<double>(0,1);  // this is a degree of freedom

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
        BHalfedge BHe = Be.halfedge();
        BFace Bf_ij = BHe.face();
        BFace Bf_ji = BHe.twin().face();
        std::complex<double> f_ij = branchCoverFields[Bf_ij.sheet][Bf_ij.f];
        std::complex<double> f_ijt = rCM[BHe.he] * branchCoverFields[Bf_ji.sheet][Bf_ji.f];
        errors[Be.e] = std::max(errors[Be.e], std::abs(f_ijt - f_ij));
    }
    std::cout << "Done!" << std::endl;
}

void QuadMesh::computeOmega() {
    // set up an edgedata per sheet
    for (int i = 0; i < n; i++) {
        EdgeData<double> omegaSheet(mesh);
        omega.push_back(omegaSheet);
    }

    // compute omega on each edge of the branch cover
    // visualize this
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        BHalfedge he_ij = Be.halfedge();
        BHalfedge he_ji = he_ij.twin();

        std::complex<double> f_ij = scale * branchCoverFields[he_ij.sheet][he_ij.face().f];
        std::complex<double> f_ji = scale * branchCoverFields[he_ji.sheet][he_ji.face().f];

        std::complex<double> e_ij = thetaCM[he_ij.he];
        std::complex<double> e_ji = -thetaCM[he_ji.he];
        assert( std::abs(rCM[he_ji.he] * e_ij - e_ji) < 1e-6 );

        // these 2 should equal each other
        double prod_ij = dot(e_ij, f_ij);
        double prod_ji = dot(e_ji, f_ji);

        if (std::abs(prod_ij - prod_ji) > 1e-6) {
            polyscope::warning("non-parallel vectors", std::to_string(std::abs(prod_ij - prod_ji)));
            //std::cout << prod_ij << "," << prod_ji << std::endl;
        }
        omega[Be.sheet][Be.e] = 0.5 * (prod_ij + prod_ji);
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
    size_t iSingular = 0;
    size_t iNonSingular = numSingularities;
    std::vector<BVertex> allBVertices = BC.allVertices();
    for (BVertex Bv : allBVertices) {
        if (singularities[Bv.v] != 0) {
            assert(Bv.sheet == BC.singularSheet);
            BVertexIndices[BC.singularSheet][Bv.v] = iSingular++;
        } else if (Bv.sheet == 0 || Bv.sheet == 1) {
            BVertexIndices[Bv.sheet][Bv.v] = iNonSingular++;
        }
    }
    assert(iSingular == (size_t)numSingularities);

    // 1 copy of singular vertices (sheet 0) and 2 copies of nonsingular vertices (sheets 0 and 1)
    // multiply by 2 to get real DOFs rather than complex    
    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular + numSingularities);
    Eigen::SparseMatrix<double> A(numPsi,numPsi);
    std::vector<Eigen::Triplet<double>> triplets;
    
    std::vector<BHalfedge> allBHalfedges = BC.allHalfedges();
    for (BHalfedge BHe : allBHalfedges) {
        // get cotan weights
        double cotA = 1.0 / tan(cmAngles[BHe.he]);
        double cotB = 1.0 / tan(cmAngles[BHe.twin().he]);
        double w = (cotA + cotB) / 2.0;

        // For vertex on some sheet, returns indices to index real and imaginary parts of psi,
        // as well as sign for imaginary componenet to implement conjugate
        auto getIndexData = [&](BVertex Bv, size_t& realInd, size_t& imagInd, double& imagSign) {
            bool v_singular = (singularities[Bv.v] != 0);
            if (v_singular) {  
                assert(Bv.sheet == BC.singularSheet);
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
        BVertex Bv_A = BHe.vertex();
        BVertex Bv_B = BHe.twin().vertex();
        size_t iA_re, iA_im, iB_re, iB_im;
        double sA_im, sB_im;
        getIndexData(Bv_A, iA_re, iA_im, sA_im);
        getIndexData(Bv_B, iB_re, iB_im, sB_im);

        // add diagonal entries
        triplets.push_back( Eigen::Triplet<double>(iA_re,iA_re,w) );
        triplets.push_back( Eigen::Triplet<double>(iA_im,iA_im,w) ); // multiply by sA_im and sB_im?

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
    size_t numPsi = 2 * (2 * numNonSingular + numSingularities);
    Eigen::SparseMatrix<double> M(numPsi, numPsi);
    std::vector<Eigen::Triplet<double>> triplets;

    auto getIndexData = [&](BVertex Bv, size_t& realInd, size_t& imagInd) {
        bool v_singular = (singularities[Bv.v] != 0);
        if (v_singular) {  
            assert(Bv.sheet == BC.singularSheet);
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

void QuadMesh::computeStripes() {
    std::cout << "Computing Stripes..." << std::endl;
    // compute 1-form and build matrices
    computeOmega();

    Eigen::SparseMatrix<double> A = energyMatrix();
    Eigen::SparseMatrix<double> B = massMatrix();

    size_t numNonSingular = mesh->nVertices() - numSingularities;
    size_t numPsi = 2 * (2 * numNonSingular + numSingularities);
    Eigen::MatrixXd x = Eigen::MatrixXd::Random(numPsi,1);

    // inverse power iteration to find eigenvector belonging to the smallest eigenvalue
    Eigen::MatrixXd prevX = x;
    for (int i = 0; i < nPowerIterations/2; i++) {
        Eigen::MatrixXd rhs = B * x;
        x = solveSquare<double>(A, B * x);
        double norm2 = (x.transpose() * B * x)(0,0);
        x = x / sqrt(norm2);
        Eigen::MatrixXd resid = x - prevX;
        std::cout << "Resid: " << resid.transpose() * B * resid  << std::endl;
        prevX = x;
    }
    std::cout << "Done!" << std::endl;

    // store resulting field and naive coords
    for (int i = 0; i < 2; i++) {
        for (VertexPtr v : mesh->vertices()) {
            int sheet = (singularities[v] == 0) ? i : BC.singularSheet;
            size_t index = BVertexIndices[sheet][v];
            double real = x(2*index,  0);
            double imag = x(2*index+1,0);
            std::complex<double> p(real,imag);
            //if (singularities[v] != 0) std::cout << p << std::endl;
            psi[i][v] = p;
            coords[i][v] = std::arg(p);
        }
    }
}

std::complex<double> QuadMesh::getPsi(BVertex Bv) {
    if (Bv.sheet == 0 || Bv.sheet == 1) {
        return psi[Bv.sheet][Bv.v];
    } else if (Bv.sheet == 2) {
        return std::conj(psi[0][Bv.v]);
    } else {
        assert(Bv.sheet == 3);
        return std::conj(psi[1][Bv.v]);
    }
}

double QuadMesh::getSigma(BEdge Be) {
    if (Be.sheet == 0 || Be.sheet == 1) {
        return sigma[Be.sheet][Be.e];
    } else if (Be.sheet == 2) {
        return -sigma[0][Be.e];
    } else {
        assert(Be.sheet == 3);
        return -sigma[1][Be.e];
    }
}

void QuadMesh::computeSigma() {
    std::vector<BEdge> allBEdges = BC.allEdges();
    for (BEdge Be : allBEdges) {
        if (Be.sheet != 0 && Be.sheet != 1) continue;

        BHalfedge BHe = Be.halfedge();
        std::complex<double> psi_i = getPsi(BHe.vertex());
        std::complex<double> psi_j = getPsi(BHe.twin().vertex());
        double w_ij = omega[Be.sheet][Be.e];

        double y = ( w_ij + std::arg(psi_i) - std::arg(psi_j) ) / (2 * PI);
        int n_ij = std::round(y);

        sigma[Be.sheet][Be.e] = std::arg(psi_j) - std::arg(psi_i) + 2 * PI * n_ij;
    }
}

void QuadMesh::textureCoordinates() {
    std::cout << "Computing texture coordinates...";
    FaceData<std::vector<double>> xCoords(mesh);
    FaceData<std::vector<double>> yCoords(mesh);
    FaceData<int> singularFace(mesh);

    computeSigma();

    std::vector<BFace> allBFaces = BC.allFaces();
    for (BFace Bf : allBFaces) {
        // sheet 0 = x coordinate, sheet 1 = y coordinate
        if (Bf.sheet != 0 && Bf.sheet != 1) continue; // need to fill in default quantities here i think

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

        // skip singularities for now
        if (singularities[Bv_i.v] != 0 || singularities[Bv_j.v] != 0 || singularities[Bv_k.v] != 0) {
            singularFace[Bf.f] = 1;
            continue;
        } else {
            singularFace[Bf.f] = 0;
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
        int zero_ind = (c_l - c_i) / (2.0 * PI);
        if (zero_ind != 0) {
            zeros[Bf.f] = zero_ind;
            std::cout << "New Singularity: " << zero_ind << std::endl;
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
                vertCoords = {0, 0}; // default value for singularities
            } else {
                vertCoords = {.x = xs[i], .y = ys[i]};
            }
            coords.push_back(vertCoords);
        } 
        texCoords[f] = coords;
    }
    std::cout << "Done!" << std::endl;
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

    // direct stripe coords
    polyscope::getSurfaceMesh()->addQuantity("X Coords", coords[0]);
    polyscope::getSurfaceMesh()->addQuantity("Y Coords", coords[1]);

    // omega
    polyscope::getSurfaceMesh()->addVectorQuantity("omega 0", omega[0]);
    polyscope::getSurfaceMesh()->addVectorQuantity("omega 1", omega[1]);
    polyscope::getSurfaceMesh()->addVectorQuantity("omega 2", omega[2]);
    polyscope::getSurfaceMesh()->addVectorQuantity("omega 3", omega[3]);

    polyscope::getSurfaceMesh()->addQuantity("edge error", errors);

    polyscope::getSurfaceMesh()->addQuantity("stripes shader", texCoords);

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
}