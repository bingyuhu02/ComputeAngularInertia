#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/constants.hpp>

//EIGEN - SVD 
#include <Eigen/Dense>
#include <Eigen\SVD>
#include <iostream>
#include <fstream>
#include <vector>
#include "TetrahedralMesh.h"

using namespace Eigen;

void convertGlmMatToEigen(glm::dmat3 &glmMat, Eigen::MatrixXd &eigenMatOut) {
	//eigenMatOut(row, col)
	eigenMatOut(0, 0) = glmMat[0][0];
	eigenMatOut(0, 1) = glmMat[0][1];
	eigenMatOut(0, 2) = glmMat[0][2];
	eigenMatOut(1, 0) = glmMat[1][0];
	eigenMatOut(1, 1) = glmMat[1][1];
	eigenMatOut(1, 2) = glmMat[1][2];
	eigenMatOut(2, 0) = glmMat[2][0];
	eigenMatOut(2, 1) = glmMat[2][1];
	eigenMatOut(2, 2) = glmMat[2][2];
}

void printVec(glm::dvec3 vec) {
	std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
}

void printMat(glm::dmat3 mat) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << mat[i][j] << " ";
		}
		std::cout << std::endl;
	}
	
}

int main() {
	//T2etra0802  //the problem tetra
	//New810BallPoints3dMeters2
	TetrahedralMesh mesh = TetrahedralMesh::BuildFromFile("Tetra1215.txt");

	if (mesh.TetrahedraCount == 0) {
		std::cout << "Error building from file" << std::endl;
	}

	//1. Center of mass:
	mesh.ComputeCentroid();
	std::cout << "Original center of mass: " << std::endl;
	printVec(mesh.MeshCentroid);
	std::cout << std::endl;

	//2. Translate to origin and set the mass of the mesh
	mesh.Translate(-mesh.MeshCentroid);

	//3. Compute angular inertia (recomputes centroid)
	mesh.ComputeMeshAttributes();
	std::cout << "Angular inertia: " << std::endl;
	printMat(mesh.AngularInertia);
	std::cout << std::endl;

	//Getting eigenvectors..
	MatrixXd angularInertia(3, 3);
	convertGlmMatToEigen(mesh.AngularInertia, angularInertia);

	SelfAdjointEigenSolver<MatrixXd> es(3);
	es.compute(angularInertia);
	typedef Eigen::SelfAdjointEigenSolver<MatrixXd>::EigenvectorsType EigenVec;

	EigenVec eigenVec1 = es.eigenvectors().col(0);
	EigenVec eigenVec2 = es.eigenvectors().col(1);
	EigenVec eigenVec3 = es.eigenvectors().col(2);

	double eigenVal1 = es.eigenvalues().col(0).row(0).real().value();
	double eigenVal2 = es.eigenvalues().col(0).row(1).real().value();
	double eigenVal3 = es.eigenvalues().col(0).row(2).real().value();

	glm::dvec3 eigenVecA(eigenVec1.row(0).real().value(), eigenVec1.row(1).real().value(), eigenVec1.row(2).real().value());
	glm::dvec3 eigenVecB(eigenVec2.row(0).real().value(), eigenVec2.row(1).real().value(), eigenVec2.row(2).real().value());
	glm::dvec3 eigenVecC(eigenVec3.row(0).real().value(), eigenVec3.row(1).real().value(), eigenVec3.row(2).real().value());

	//Eigenvectors are columns of matrix P:
	glm::dmat3 P(0.0);
	P[0][0] = eigenVecA.x;
	P[0][1] = eigenVecA.y;
	P[0][2] = eigenVecA.z;

	P[1][0] = eigenVecB.x;
	P[1][1] = eigenVecB.y;
	P[1][2] = eigenVecB.z;

	P[2][0] = eigenVecC.x;
	P[2][1] = eigenVecC.y;
	P[2][2] = eigenVecC.z;

	std::cout << "P matrix: " << std::endl;
	printMat(P);
	std::cout << std::endl;

	// A = P D trans(P)
	// => D = trans(P) A P
	glm::dmat3 pTrans = glm::transpose(P);
	glm::dmat3 D = pTrans * mesh.AngularInertia * P;

	std::cout << "D matrix: " << std::endl;
	printMat(D);
	std::cout << std::endl;
	getchar();
	return 0;
}