#ifndef BVH_HPP
#define BVH_HPP

#include <Eigen/Sparse>
#include <Eigen/Geometry> 
#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>

using namespace std;

typedef Eigen::Vector3d Vert;
typedef Eigen::Vector3i Face;

// double NegativeInf = numeric_limits<double>::min();

struct Node
{
	Eigen::Vector3d center;
	
	// [xmin, xmax, ymin, ymax, zmin, zmax]
	vector<double> box;

	// face contained => index of face in faces variable
	vector<int> node_faces;

	// if it is leaf node
	bool is_leaf;

	bool is_null;

	int left;  //left child
	int right; //right child

	int dimensionIdx; // 0:x, 1:y, 2:z
};

class bvh
{
public:
	bvh(int _max_faces_per_leaf, double _min_thickness);
	~bvh();

	void assignObj(vector<Vert> _vertices, vector<Face> _faces);

	void selfCollision(vector<Eigen::Vector2i>& intersectFaces);
	void collisonDetection(vector<Node>& T, int idx, Face f, vector<int>& intersectFaces);
	void collisonDetection(vector<Node>& T1, vector<Node>& T2, int idx1, int idx2, vector<Eigen::Vector2i>& intersectFaces);
	bool checkBoxIntersect(Node& N1, Node& N2);
	bool checkBoxIntersect(vector<double>& box1, vector<double>& box2, Vert& v1, Vert& v2);

	bool checkFaceIntersect(Face f1, Face f2, Vert& forceVec, double& outMinDist, int& f2Index);
	bool checkFaceIntersect(Face f1, Face f2);
	bool checkFaceBoxIntersect(Face f, Node& N);
	void facePlane(Vert v1, Vert v2, Vert v3, vector<double>& outPlane);
	void facePlane(Face f, std::vector<double>& outPlane);
	double vertex2plane(Vert v, std::vector<double> plane);
	bool vertexInTriangle(Vert v, Face f, Vert faceNorm);
	bool checkedBefore(Face f1, Face f2, vector<Eigen::Vector2i>& intersectFaces);
	bool sameFace(Face f1, Face f2);

	void buildTree();
	void build(int p, int r, int nodeIdx);
	int splitFaces(int p, int r, vector<double> box, int& dimensionIdx);
	Vert findCenter(int p, int r);
	vector<double> findBoundBox(int p, int r);

	Vert faceCenter(Face f);

	Face getFace(int idx);
	bool shareSide(Face f1, Face f2);

public:
	vector<Vert> vertices;
	vector<Face> faces;
	vector<Vert> faceCenters;

	vector<Node> treeNodes;


	int n_vertices;
	int n_faces;
	int max_faces_per_leaf;
	double min_thickness;
};

#endif