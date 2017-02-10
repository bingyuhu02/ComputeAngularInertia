#pragma once
#include <glm\common.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <string>
#include <fstream>

class Tetrahedron {
public:
	//Properties set by caller:
	glm::dvec3 a;
	glm::dvec3 b;
	glm::dvec3 c;
	glm::dvec3 d;
	//Normal vectors for lighting
	glm::dvec3 normalACD;
	glm::dvec3 normalABC;
	glm::dvec3 normalCDB;
	glm::dvec3 normalABD;
	double volume;
	double mass;
	double density;

	//Set by TetrahedralMesh class:
	glm::dvec3 centroid;
	glm::dmat3 angularInertia;

	Tetrahedron() {
		volume = 1.0;
		mass = 1.0;
		density = 1.0;
	}

	Tetrahedron(glm::dvec3 pa, glm::dvec3 pb, glm::dvec3 pc, glm::dvec3 pd) : a(pa), b(pb), c(pc), d(pd) {
		computeVolume();
		mass = 1.0;
		density = 1.0;
		normalACD = getNormalOfTriangle(a, c, d);
		normalABC = getNormalOfTriangle(a, b, c);
		normalCDB = getNormalOfTriangle(c, d, b);
		normalABD = getNormalOfTriangle(a, b, d);
	}

	Tetrahedron(glm::dvec3 pa, glm::dvec3 pb, glm::dvec3 pc, glm::dvec3 pd, double pVolume, double pMass, double pDensity) :
		a(pa), b(pb), c(pc), d(pd), volume(pVolume), mass(pMass), density(pDensity) {

		mass = 1.0;
		density = 1.0;
		normalACD = getNormalOfTriangle(a, c, d);
		normalABC = getNormalOfTriangle(a, b, c);
		normalCDB = getNormalOfTriangle(c, d, b);
		normalABD = getNormalOfTriangle(a, b, d);
	}

	Tetrahedron(glm::dvec3 pa, glm::dvec3 pb, glm::dvec3 pc, glm::dvec3 pd, glm::dvec3 pNormalACD, glm::dvec3 pNormalABC,
		glm::dvec3 pNormalCDB, glm::dvec3 pNormalABD) : a(pa), b(pb), c(pc), d(pd), normalACD(pNormalACD),
		normalABC(pNormalABC), normalCDB(pNormalCDB), normalABD(pNormalABD) {
		computeVolume();
		mass = 1.0;
		density = 1.0;
	}

	Tetrahedron(glm::dvec3 pa, glm::dvec3 pb, glm::dvec3 pc, glm::dvec3 pd, glm::dvec3 pNormalACD, glm::dvec3 pNormalABC,
		glm::dvec3 pNormalCDB, glm::dvec3 pNormalABD, double pVolume, double pMass, double pDensity) :
		a(pa), b(pb), c(pc), d(pd), normalACD(pNormalACD), normalABC(pNormalABC), normalCDB(pNormalCDB),
		normalABD(pNormalABD), volume(pVolume), mass(pMass), density(pDensity) {

	}

	void computeVolume() {
		glm::mat4 volMat;
		volMat[0] = glm::vec4(a.x, b.x, c.x, d.x);
		volMat[1] = glm::vec4(a.y, b.y, c.y, d.y);
		volMat[2] = glm::vec4(a.z, b.z, c.z, d.z);
		volMat[3] = glm::vec4(1.0f, 1.0f, 1.0f, 1.0f);


		double det = abs(glm::determinant(volMat));
		this->volume = det / (6.0f); // volume formula divides by 3!

		this->mass = volume * density;
	}

	void Translate(glm::dvec3 trans) {
		a += trans;
		b += trans;
		c += trans;
		d += trans;
	}

	glm::dvec3 getNormalOfTriangle(glm::dvec3 &p1, glm::dvec3 &p2, glm::dvec3 &p3) {
		glm::dvec3 normal;
		glm::dvec3 V = p2 - p1;
		glm::dvec3 W = p3 - p1;
		normal.x = (V.y * W.z) - (V.z * W.y);
		normal.y = (V.z * W.z) - (V.x * W.z);
		normal.z = (V.x * W.y) - (V.y * W.x);
		return normal;
	}

};

class TetrahedralMesh {
public:
	int TetrahedraCount;
	std::vector<Tetrahedron> Tetrahedra;
	glm::dvec3 MeshCentroid;
	glm::dmat3 AngularInertia;
	double Volume;
	double Mass;
	double Density;

	TetrahedralMesh() : TetrahedraCount(0) {}

	void AddTetrahedron(Tetrahedron &th) {
		Tetrahedra.push_back(th);
		TetrahedraCount++;
	}

	//First aggregates mass & volumes of tetrahedra, while computing their 
	// centroids & individual angular inertia matrices
	//Then averages centroids for mesh centroid and 
	// sums the tetrahedra inertia w.r.t. the mesh centroid
	// to get the mesh's angular inertia matrix
	void ComputeMeshAttributes() {
		typedef std::vector<Tetrahedron>::iterator iter;

		glm::dvec3 centroid(0.0, 0.0, 0.0);
		double vol = 0.0f;
		double mass = 0.0f;
		for (iter it = Tetrahedra.begin(); it != Tetrahedra.end(); it++) {
			this->computeTetrahedronAngularInertia(*it);
			centroid += it->centroid * it->volume;
			vol += it->volume;
			mass += it->mass;
		}
		//Take average of tetrahedra centroids as mesh centroid
		this->MeshCentroid = centroid / vol;
		this->Volume = vol;
		this->Mass = mass;

		//Translate to origin:
		//this->Translate(-this->MeshCentroid);
		//this->MeshCentroid = glm::dvec3(0.0, 0.0, 0.0);

		//Now have inertia matrix of each tetrahedra wrt its own centroid, compute
		// with respect to mesh centroid:
		glm::mat3 meshInertia(0.0);
		for (iter it = Tetrahedra.begin(); it != Tetrahedra.end(); it++) {
			//glm::mat3 inertiaWrtMesh;
			//this->angularInertiaWithRespectToPoint(this->MeshCentroid, *it, inertiaWrtMesh);
			//meshInertia += inertiaWrtMesh;

			computeTetrahedronAngularInertia2(*it, this->MeshCentroid);
			meshInertia += it->angularInertia;
		}
		this->AngularInertia = meshInertia;
	}

	void Translate(glm::dvec3 trans) {
		for (std::vector<Tetrahedron>::iterator it = Tetrahedra.begin(); it != Tetrahedra.end(); it++) {
			it->Translate(trans);
		}
	}

	int GetTetrahedraCount() {
		return this->TetrahedraCount;
	}

	static TetrahedralMesh BuildFromFile(std::string filename) {
		TetrahedralMesh mesh;
		std::string line;
		std::ifstream meshfile(filename);
		if (meshfile.is_open()) {
			std::vector<glm::dvec3> vertices;
			bool vertex = false;
			bool tetrahedra = false;
			while (getline(meshfile, line)) {
				if (!vertex) {
					if (line.compare("Vertices") == 0) {
						vertex = true;
					}
				}
				else if (vertex && !tetrahedra) {
					if (line.compare("Tetrahedra") == 0) {
						tetrahedra = true;
						continue;
					}

					std::vector<std::string> splitLine = TetrahedralMesh::split(line, ' ');
					vertices.push_back(glm::dvec3(stof(splitLine[0]), stof(splitLine[1]), stof(splitLine[2]))); //drop the homogeneous coord
				}
				else if (tetrahedra) {
					//Read in tetrahedra indices:
					if (line.compare("End") == 0)
						break;
					std::vector<std::string> splitLine = TetrahedralMesh::split(line, ' ');
					int indexA = stoi(splitLine[0]) - 1;
					int indexB = stoi(splitLine[1]) - 1;
					int indexC = stoi(splitLine[2]) - 1;
					int indexD = stoi(splitLine[3]) - 1;
					Tetrahedron t(vertices[indexA], vertices[indexB], vertices[indexC], vertices[indexD]);
					mesh.AddTetrahedron(t);
				}
			}
			meshfile.close();
		}
		else {
			std::cout << "Failed opening " << filename << std::endl;
		}

		return mesh;
	}


	static void computeTetrahedronCentroid(Tetrahedron &th) {
		double centroidX = (th.a.x + th.b.x + th.c.x + th.d.x) / 4.0;
		double centroidY = (th.a.y + th.b.y + th.c.y + th.d.y) / 4.0;
		double centroidZ = (th.a.z + th.b.z + th.c.z + th.d.z) / 4.0;
		th.centroid = glm::dvec3(centroidX, centroidY, centroidZ);
	}

	//Computes A.I. wrt the mesh's centroid instead of tetrahedron's centroid,
	// then there's no need for the angularInertiaWithRespectToPoint conversion
	void computeTetrahedronAngularInertia2(Tetrahedron &th, glm::dvec3 &meshCentroid) {

		computeTetrahedronCentroid(th);

		//use vertices wrt centroid:
		glm::dvec3 a = th.a - meshCentroid; //basically translates centroid to origin
		glm::dvec3 b = th.b - meshCentroid;
		glm::dvec3 c = th.c - meshCentroid;
		glm::dvec3 d = th.d - meshCentroid;

		th.computeVolume();

		//DET of Jacobian, see http://docsdrive.com/pdfs/sciencepublications/jmssp/2005/8-11.pdf
		double DETJ = 6 * th.volume;
		//inertia tensor matrix:
		// | a  -b' -c' |
		// |-b'  b  -a' |
		// |-c' -a'  c  |

		double density = th.density;
		if (density == 0.0)
			density = 1.0;
		double aa = (density * DETJ * (a.y*a.y + a.y*b.y + b.y*b.y + a.y*c.y + b.y*c.y + c.y*c.y + a.y * d.y
			+ b.y*d.y + c.y*d.y + d.y*d.y + a.z*a.z + a.z*b.z + b.z*b.z + a.z*c.z + b.z*c.z + c.z*c.z + a.z * d.z
			+ b.z*d.z + c.z*d.z + d.z*d.z)) / 60.0;

		double bb = (density * DETJ * (a.x*a.x + a.x*b.x + b.x*b.x + a.x*c.x + b.x*c.x + c.x*c.x + a.x * d.x
			+ b.x*d.x + c.x*d.x + d.x*d.x + a.z*a.z + a.z*b.z + b.z*b.z + a.z*c.z + b.z*c.z + c.z*c.z + a.z * d.z
			+ b.z*d.z + c.z*d.z + d.z*d.z)) / 60.0;

		double cc = (density * DETJ * (a.x*a.x + a.x*b.x + b.x*b.x + a.x*c.x + b.x*c.x + c.x*c.x + a.x * d.x
			+ b.x*d.x + c.x*d.x + d.x*d.x + a.y*a.y + a.y*b.y + b.y*b.y + a.y*c.y + b.y*c.y + c.y*c.y + a.y * d.y
			+ b.y*d.y + c.y*d.y + d.y*d.y)) / 60.0;

		double aPrime = (density * DETJ * (2.0*a.y*a.z + b.y*a.z + c.y*a.z + d.y*a.z + a.y*b.z + 2.0*b.y*b.z
			+ c.y*b.z + d.y*b.z + a.y*c.z + b.y*c.z + 2.0*c.y*c.z + d.y*c.z + a.y*d.z + b.y*d.z
			+ c.y*d.z + 2.0*d.y*d.z)) / 120.0;

		double bPrime = (density * DETJ * (2.0*a.x*a.z + b.x*a.z + c.x*a.z + d.x*a.z + a.x*b.z + 2.0*b.x*b.z
			+ c.x*b.z + d.x*b.z + a.x*c.z + b.x*c.z + 2.0*c.x*c.z + d.x*c.z + a.x*d.z + b.x*d.z
			+ c.x*d.z + 2.0*d.x*d.z)) / 120.0;

		double cPrime = (density * DETJ * (2.0*a.x*a.y + b.x*a.y + c.x*a.y + d.x*a.y + a.x*b.y + 2.0*b.x*b.y
			+ c.x*b.y + d.x*b.y + a.x*c.y + b.x*c.y + 2.0*c.x*c.y + d.x*c.y + a.x*d.y + b.x*d.y
			+ c.x*d.y + 2.0*d.x*d.y)) / 120.0;

		th.angularInertia = glm::mat3(aa, -bPrime, -cPrime, -bPrime, bb, -aPrime, -cPrime, -aPrime, cc);
	}


private:

	void angularInertiaWithRespectToPoint(glm::dvec3 &point, Tetrahedron &th, glm::mat3 &inertiaWrtPoint) {
		// I is th's angular inertia matrix, I_3 is 3x3 identitiy, r is vector from point to th's centroid , r^t = r transpose,
		// I + m (r^t * r * I_3 - r * r^t)    

		glm::dmat3 I3; // 3x3 identity
		glm::dvec3 ptToCentroid = th.centroid - point; // vector from point to centroid
		glm::dmat3 ptToCtimesTranspose = glm::outerProduct(ptToCentroid, ptToCentroid); // ptToC * (ptToC)^T

		inertiaWrtPoint = th.angularInertia + th.mass * (glm::dot(ptToCentroid, ptToCentroid) * I3 - ptToCtimesTranspose);
	}

	void computeTetrahedronAngularInertia(Tetrahedron &th) {

		computeTetrahedronCentroid(th);

		//use vertices wrt centroid:
		glm::dvec3 a = th.a - th.centroid; //basically translates centroid to origin
		glm::dvec3 b = th.b - th.centroid;
		glm::dvec3 c = th.c - th.centroid;
		glm::dvec3 d = th.d - th.centroid;

		th.computeVolume();

		//DET of Jacobian, see http://docsdrive.com/pdfs/sciencepublications/jmssp/2005/8-11.pdf
		double DETJ = 6 * th.volume;
		//inertia tensor matrix:
		// | a  -b' -c' |
		// |-b'  b  -a' |
		// |-c' -a'  c  |

		double density = th.density;
		if (density == 0.0)
			density = 1.0;
		double aa = (density * DETJ * (a.y*a.y + a.y*b.y + b.y*b.y + a.y*c.y + b.y*c.y + c.y*c.y + a.y * d.y
			+ b.y*d.y + c.y*d.y + d.y*d.y + a.z*a.z + a.z*b.z + b.z*b.z + a.z*c.z + b.z*c.z + c.z*c.z + a.z * d.z
			+ b.z*d.z + c.z*d.z + d.z*d.z)) / 60.0;

		double bb = (density * DETJ * (a.x*a.x + a.x*b.x + b.x*b.x + a.x*c.x + b.x*c.x + c.x*c.x + a.x * d.x
			+ b.x*d.x + c.x*d.x + d.x*d.x + a.z*a.z + a.z*b.z + b.z*b.z + a.z*c.z + b.z*c.z + c.z*c.z + a.z * d.z
			+ b.z*d.z + c.z*d.z + d.z*d.z)) / 60.0;

		double cc = (density * DETJ * (a.x*a.x + a.x*b.x + b.x*b.x + a.x*c.x + b.x*c.x + c.x*c.x + a.x * d.x
			+ b.x*d.x + c.x*d.x + d.x*d.x + a.y*a.y + a.y*b.y + b.y*b.y + a.y*c.y + b.y*c.y + c.y*c.y + a.y * d.y
			+ b.y*d.y + c.y*d.y + d.y*d.y)) / 60.0;

		double aPrime = (density * DETJ * (2.0*a.y*a.z + b.y*a.z + c.y*a.z + d.y*a.z + a.y*b.z + 2.0*b.y*b.z
			+ c.y*b.z + d.y*b.z + a.y*c.z + b.y*c.z + 2.0*c.y*c.z + d.y*c.z + a.y*d.z + b.y*d.z
			+ c.y*d.z + 2.0*d.y*d.z)) / 120.0;

		double bPrime = (density * DETJ * (2.0*a.x*a.z + b.x*a.z + c.x*a.z + d.x*a.z + a.x*b.z + 2.0*b.x*b.z
			+ c.x*b.z + d.x*b.z + a.x*c.z + b.x*c.z + 2.0*c.x*c.z + d.x*c.z + a.x*d.z + b.x*d.z
			+ c.x*d.z + 2.0*d.x*d.z)) / 120.0;

		double cPrime = (density * DETJ * (2.0*a.x*a.y + b.x*a.y + c.x*a.y + d.x*a.y + a.x*b.y + 2.0*b.x*b.y
			+ c.x*b.y + d.x*b.y + a.x*c.y + b.x*c.y + 2.0*c.x*c.y + d.x*c.y + a.x*d.y + b.x*d.y
			+ c.x*d.y + 2.0*d.x*d.y)) / 120.0;

		th.angularInertia = glm::mat3(aa, -bPrime, -cPrime, -bPrime, bb, -aPrime, -cPrime, -aPrime, cc);
	}

	static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			elems.push_back(item);
		}
		return elems;
	}


	static std::vector<std::string> split(const std::string &s, char delim) {
		std::vector<std::string> elems;
		split(s, delim, elems);
		return elems;
	}

};