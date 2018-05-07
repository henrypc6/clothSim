#include "bvh.hpp"

bvh::bvh(int _max_faces_per_leaf, double _min_thickness):max_faces_per_leaf(_max_faces_per_leaf), min_thickness(_min_thickness)
{
}

bvh::~bvh()
{

}

void bvh::assignObj(vector<Vert> _vertices, vector<Face> _faces)
{
	vertices.clear();
	faces.clear();

	vertices = _vertices;
	faces = _faces;

	// assign centers of face in prior
	faceCenters.clear();
	for(int i = 0 ; i < faces.size(); ++i)
	{
		faceCenters.push_back(faceCenter(faces[i]));
	}

	n_vertices = vertices.size();
	n_faces = faces.size();
}

void bvh::selfCollision(vector<Eigen::Vector2i>& intersectFaces)
{
	collisonDetection(treeNodes, treeNodes, 0, 0, intersectFaces);

	cout<<"intersecting faces "<<intersectFaces.size()<<endl;
}

// check collison of a face to a set of faces 
void bvh::collisonDetection(vector<Node>& T, int idx, Face f, vector<int>& intersectFaces)
{
	if(idx<0 || idx >= T.size())
		return;

	Node N = T[idx];

	if(N.is_null)
		return;

	if(!checkFaceBoxIntersect(f,N))
		return;

	if(N.is_leaf)
	{
		for(int i = 0 ; i < N.node_faces.size(); ++i)
		{

			Face f2 = faces[N.node_faces[i]];
			if(checkFaceIntersect(f,f2))
			{
				intersectFaces.push_back(N.node_faces[i]);
			}
		}
		return;
	}

	if(!N.is_leaf)
	{
		collisonDetection(T, N.left, f, intersectFaces);
		collisonDetection(T, N.right, f, intersectFaces);
	}


}



// check self collison
void bvh::collisonDetection(vector<Node>& T1, vector<Node>& T2, int idx1, int idx2, vector<Eigen::Vector2i>& intersectFaces)
{

	// check if the node exists
	if(idx1 < 0 || idx2 < 0)
		return;

	if(idx1 >= T1.size() || idx2 >= T2.size())
		return;

	Node N1 = T1[idx1];
	Node N2 = T2[idx2];

	if(N1.is_null || N2.is_null)
		return;

	// check if the box is intersecting
	bool boxInter = checkBoxIntersect(N1,N2);


	if(!boxInter)
		return;
	
	// both are at the leaf node
	if(N1.is_leaf && N2.is_leaf)
	{
		for(int i = 0; i < N1.node_faces.size(); ++i)
		{
			Face f1 = faces[N1.node_faces[i]];
			for(int j = 0; j < N2.node_faces.size(); ++j)
			{
				Face f2 = faces[N2.node_faces[j]];

				if(f1 == f2)
					continue;

				if(shareSide(f1,f2))
					continue;

				if(checkedBefore(f1,f2, intersectFaces))
					continue;

				if(checkFaceIntersect(f1,f2))
				{
					// cout<<"intersect "<<N1.node_faces[i]<<" "<<N2.node_faces[j]<<" "<<endl;
					intersectFaces.push_back(Eigen::Vector2i(N1.node_faces[i],N2.node_faces[j]));
				}
			}
		}
		return;
	}
	
	if(!N1.is_leaf || !N2.is_leaf)
	{
		collisonDetection(T1,T2,N1.left, N2.left, intersectFaces);
		collisonDetection(T1,T2,N1.right, N2.right, intersectFaces);
		collisonDetection(T1,T2,N1.left, N2.right, intersectFaces);
		collisonDetection(T1,T2,N1.right, N2.left, intersectFaces);
	}
}

bool bvh::checkedBefore(Face f1, Face f2, vector<Eigen::Vector2i>& intersectFaces)
{
	for(int i =0; i < intersectFaces.size(); ++i)
	{
		if(f2 == faces[intersectFaces[i][0]] && f1 == faces[intersectFaces[i][1]])
			return true;
	}
	return false;
}

bool bvh::sameFace(Face f1, Face f2)
{
	return (f1[0]==f2[0] && f1[1]==f2[1] && f1[2]==f2[2]);
}

bool bvh::checkBoxIntersect(Node& N1, Node& N2)
{
	return checkBoxIntersect(N1.box, N2.box, N1.center, N2.center);
}

bool bvh::checkBoxIntersect(vector<double>& box1, vector<double>& box2, Vert& v1, Vert& v2)
{
	// check x direction
	double dx = (box1[1] - box1[0])/2.0 +(box2[1]-box2[0])/2.0;
	if(fabs(v1[0] - v2[0]) < dx )
		return true;

	// check y direction
	double dy = (box1[3] - box1[2])/2.0 +(box2[3]-box2[2])/2.0;
	if(fabs(v1[1] - v2[1]) < dy )
		return true;

	// check z direction
	double dz = (box1[5] - box1[4])/2.0 +(box2[5]-box2[4])/2.0;
	if(fabs(v1[2] - v2[2]) < dz )
		return true;

	return false;
}

bool bvh::checkFaceIntersect(Face f1, Face f2, Vert& forceVec, double& outMinDist, int& f2Index)
{
	/////////////////segment f2 into 3 line seg and check them with triangle f1//////////////////////
	// find face plane
	Vert f2_v1 = vertices[f2[0]];
	Vert f2_v2 = vertices[f2[1]];
	Vert f2_v3 = vertices[f2[2]];

	vector<double> p1(4,0);
	facePlane(f1,p1);

	// check if the vertices are all on the same side 
	double v1Dist = vertex2plane(f2_v1, p1);
	double v2Dist = vertex2plane(f2_v2, p1);
	double v3Dist = vertex2plane(f2_v3, p1);


	if( (v1Dist * v2Dist > 0) && (v1Dist * v3Dist > 0) )
		return false;

	// check if parallel
	double tol = 0.0001;
	if( fabs(v1Dist-v2Dist) < tol && fabs(v1Dist - v3Dist) < tol)
	{
		// if( fabs(v1Dist) > min_thickness)
		return false;
		// else
			// return true;
	}

	// cout<<"checking face 1 "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
	// cout<<"checking face 2 "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;

	//////////////////////////Check if each line seg in f2 intersect with triangle f1////////////////
	Vert N(p1[0],p1[1],p1[2]); // plane normal
	double d = p1[3];

	// first ray starting from v2
	bool seg1B = false;
	Vert seg1 = f2_v1 - f2_v2;
	double t1 = -(d + N.dot(f2_v2))/N.dot(seg1);
	if(t1 >=0 && t1<=1.0) // only check if the line seg (seg1) intersect the plane
	{
		Vert v = f2_v2 + t1 * seg1;
		seg1B = vertexInTriangle(v, f1, N);
	}

	//second  
	bool seg2B = false;
	Vert seg2 = f2_v1 - f2_v3;
	double t2 = -(d + N.dot(f2_v3))/N.dot(seg2);
	if(t2 >=0 && t2<=1.0)
	{
		Vert v = f2_v3 + t2 * seg2;
		seg2B = vertexInTriangle(v, f1, N);
	}

	//third  
	bool seg3B = false;
	Vert seg3 = f2_v2 - f2_v3;
	double t3 = -(d + N.dot(f2_v3))/N.dot(seg3);
	if(t3 >=0 && t3<=1.0)
	{
		Vert v = f2_v3 + t3 * seg3;
		seg3B = vertexInTriangle(v, f1, N);
	}

	bool ret = (seg1B || seg2B || seg3B);

	if(ret)
	{
		double minDist = min(min(fabs(v1Dist),fabs(v2Dist)),fabs(v3Dist));
		if(minDist == fabs(v1Dist))
		{
			f2Index = 0;
			outMinDist = v1Dist;
		}
		else if(minDist == fabs(v2Dist))
		{
			f2Index = 1;
			outMinDist = v2Dist;
		}
		else
		{
			f2Index = 2;
			outMinDist = v3Dist;
		}

		forceVec = N;

	}

	return ret;
}

bool bvh::checkFaceIntersect(Face f1, Face f2)
{
	/////////////////segment f2 into 3 line seg and check them with triangle f1//////////////////////
	// find face plane
	Vert f2_v1 = vertices[f2[0]];
	Vert f2_v2 = vertices[f2[1]];
	Vert f2_v3 = vertices[f2[2]];

	vector<double> p1(4,0);
	facePlane(f1,p1);

	// check if the vertices are all on the same side 
	double v1Dist = vertex2plane(f2_v1, p1);
	double v2Dist = vertex2plane(f2_v2, p1);
	double v3Dist = vertex2plane(f2_v3, p1);

	double tol = 0.0001;
	if( (v1Dist * v2Dist > 0) && (v1Dist * v3Dist > 0) )
		return false;

	// check if parallel
	if( fabs(v1Dist-v2Dist) < tol && fabs(v1Dist - v3Dist) < tol)
	{
		// if( fabs(v1Dist) > min_thickness)
		return false;
		// else
			// return true;
	}

	// cout<<"checking face 1 "<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
	// cout<<"checking face 2 "<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;

	//////////////////////////Check if each line seg in f2 intersect with triangle f1////////////////
	Vert N(p1[0],p1[1],p1[2]); // plane normal
	double d = p1[3];

	// first ray starting from v2
	bool seg1B = false;
	Vert seg1 = f2_v1 - f2_v2;
	double t1 = -(d + N.dot(f2_v2))/N.dot(seg1);
	if(t1 >=0 && t1<=1.0) // only check if the line seg (seg1) intersect the plane
	{
		Vert v = f2_v2 + t1 * seg1;
		seg1B = vertexInTriangle(v, f1, N);
	}

	//second  
	bool seg2B = false;
	Vert seg2 = f2_v1 - f2_v3;
	double t2 = -(d + N.dot(f2_v3))/N.dot(seg2);
	if(t2 >=0 && t2<=1.0)
	{
		Vert v = f2_v3 + t2 * seg2;
		seg2B = vertexInTriangle(v, f1, N);
	}

	//third  
	bool seg3B = false;
	Vert seg3 = f2_v2 - f2_v3;
	double t3 = -(d + N.dot(f2_v3))/N.dot(seg3);
	if(t3 >=0 && t3<=1.0)
	{
		Vert v = f2_v3 + t3 * seg3;
		seg3B = vertexInTriangle(v, f1, N);
	}
	// if(seg1B || seg2B || seg3B)
	// 	cout<<seg1B<<" "<<seg2B<<" "<<seg3B<<endl;

	return (seg1B || seg2B || seg3B);
}

bool bvh::checkFaceBoxIntersect(Face f, Node& N)
{
	Vert fc = faceCenter(f);
	Vert Nc = N.center;

	// find face bounding box 
	Vert v1 = vertices[f[0]];
	Vert v2 = vertices[f[1]];
	Vert v3 = vertices[f[2]];

	vector<double> faceBox(6,0);
	faceBox[0] = min(min(v1[0],v2[0]),v3[0]);
	faceBox[1] = max(max(v1[0],v2[0]),v3[0]);
	faceBox[2] = min(min(v1[1],v2[1]),v3[1]);
	faceBox[3] = max(max(v1[1],v2[1]),v3[1]);
	faceBox[4] = min(min(v1[2],v2[2]),v3[2]);
	faceBox[5] = max(max(v1[2],v2[2]),v3[2]);

	return checkBoxIntersect(faceBox, N.box, fc, Nc);
}

bool bvh::vertexInTriangle(Vert v, Face f, Vert faceNorm)
{
	Vert v1 = vertices[f[0]];
	Vert v2 = vertices[f[1]];
	Vert v3 = vertices[f[2]];

	Vert l1 = v1-v2;
	Vert l2 = v2-v3;
	Vert l3 = v3-v1;

	bool l1b = faceNorm.dot(l1.cross((v-v2))) <=0;
	bool l2b = faceNorm.dot(l2.cross((v-v3))) <=0;
	bool l3b = faceNorm.dot(l3.cross((v-v1))) <=0;
	
	// cout<<"in triangle check "<<l1b<<" "<<l2b<<" "<<l3b<<endl;

	return (l1b && l2b && l3b);
}

double bvh::vertex2plane(Vert v, std::vector<double> plane)
{
	Eigen::Vector4d P(plane[0],plane[1],plane[2],plane[3]);
	Eigen::Vector4d V(v[0],v[1],v[2],1);

	return V.dot(P);
}

void bvh::facePlane(Face f, std::vector<double>& outPlane)
{
	facePlane(vertices[f[0]],vertices[f[1]],vertices[f[2]],outPlane);
}

void bvh::facePlane(Vert v1, Vert v2, Vert v3, vector<double>& outPlane)
{
	Vert v21 = v2-v1;
	Vert v31 = v3-v1;
	Vert N = v21.cross(v31);

	N = N/N.norm();

	double d = -N.dot(v1);

	outPlane[0] = N[0];
	outPlane[1] = N[1];
	outPlane[2] = N[2];
	outPlane[3] = d;
}


void bvh::buildTree()
{
	// cout<<"n faces "<<n_faces<<endl;
	// cout<<"treeNodes "<<treeNodes.size()<<endl;
	double maxLevel = ceil(log2((double)n_faces));
	int treeSize = pow(2,maxLevel)-1;
	// cout<<"tree Size "<<treeSize<<endl;

	treeNodes.clear();
	for(int i = 0 ; i < treeSize; ++i)
	{
		treeNodes.push_back(Node());
		treeNodes[i].is_leaf = false;
		treeNodes[i].center = Eigen::Vector3d(0,0,0);
		treeNodes[i].box = vector<double>(6,0);
		treeNodes[i].left = -1;
		treeNodes[i].right = -1;
		treeNodes[i].dimensionIdx = 0;
		treeNodes[i].is_null = true;
	}
	build(0, faces.size()-1, 0);
}

void bvh::build(int p, int r, int nodeIdx)
{
	// add to current node
	treeNodes[nodeIdx].center = findCenter(p, r);
	treeNodes[nodeIdx].box = findBoundBox(p, r);
	treeNodes[nodeIdx].is_leaf = false;
	treeNodes[nodeIdx].is_null = false;

	// check if it is leaf node 
	if((r-p) < max_faces_per_leaf)
	{
		treeNodes[nodeIdx].is_leaf = true;

		//store the faces in the leaf node 
		for(int i = p; i<=r; ++i)
		{
			treeNodes[nodeIdx].node_faces.push_back(i);
		}
		return;
	}

	// find bounding box 
	vector<double> box = findBoundBox(p, r);

	// split it in two parts 
	int outDim = 0;
	int q = splitFaces(p,r, box, outDim); // range from [p,r+1]

	// for bad splits 
	if((q == p || q == r+1) && (r-q) >= max_faces_per_leaf)
	{
		// cout<<" it is a bad split "<<endl;
		q = ceil((p+r)/2.0);
	}

	// split ratio
	double badnessRatio = (double)(q-p)/(double)(r-p+1);
	if (badnessRatio > 0.8 || badnessRatio < 0.2)
	{
		q = ceil((p+r)/2.0);
	}

	// assign split dimension
	treeNodes[nodeIdx].dimensionIdx = outDim;

	// left parts
	treeNodes[nodeIdx].left = nodeIdx*2+1;
	build(p,q-1, nodeIdx*2+1);
	

	// right parts
	treeNodes[nodeIdx].right = nodeIdx*2+2;
	build(q,r, nodeIdx*2+2);
}

int bvh::splitFaces(int p, int r, vector<double> box, int& rangeIdx)
{
	// find the largest ranges
	double maxRange = -1;
	for (int i=0; i<3; ++i)
	{
		if(maxRange < (box[2*i+1]-box[2*i]))
		{
			maxRange = box[2*i+1]-box[2*i];
			rangeIdx = i;
		}
	}

	// 
	vector<double> splitRange;
	double pivot = 0;
	for(int i = p ; i<=r ; ++i)
	{
		Vert center = faceCenters[i];
		pivot += center[rangeIdx];
	}
	pivot = pivot/(double)(r-p+1);

	// split based on the means
	int i = p;
	for(int j = p; j<=r; ++j)
	{
		if(faceCenters[j][rangeIdx] <= pivot)
		{
			// swap the faces and centers 
			Vert tmp1 = faceCenters[i];
			faceCenters[i] = faceCenters[j];
			faceCenters[j] = tmp1;

			Face tmp2 = faces[i];
			faces[i] = faces[j];
			faces[j] = tmp2;
			++i;

		}
	}
	return i;
}


Vert bvh::faceCenter(Face f)
{
	return (vertices[f[0]] + vertices[f[1]] + vertices[f[2]])/3.0;
}

vector<double> bvh::findBoundBox(int p, int r)
{
	vector<double> box(6);

	vector<double> xvec,yvec,zvec;
	for(int i = p ; i <= r; ++i)
	{
		xvec.push_back(vertices[faces[i][0]][0]);
		xvec.push_back(vertices[faces[i][1]][0]);
		xvec.push_back(vertices[faces[i][2]][0]);

		yvec.push_back(vertices[faces[i][0]][1]);
		yvec.push_back(vertices[faces[i][1]][1]);
		yvec.push_back(vertices[faces[i][2]][1]);

		zvec.push_back(vertices[faces[i][0]][2]);
		zvec.push_back(vertices[faces[i][1]][2]);
		zvec.push_back(vertices[faces[i][2]][2]);
	}

	// sort the vector
	sort(xvec.begin(), xvec.end());
	sort(yvec.begin(), yvec.end());
	sort(zvec.begin(), zvec.end());

	box[0] = xvec[0];
	box[1] = xvec[xvec.size()-1];


	box[2] = yvec[0];
	box[3] = yvec[yvec.size()-1];
	

	box[4] = zvec[0];
	box[5] = zvec[zvec.size()-1];
	
	// make sure the bounding box is not too small
	if((box[1]-box[0]) < min_thickness)
		box[1] = box[0] + min_thickness;
	if((box[3]-box[2]) < min_thickness)
		box[3] = box[2] + min_thickness;
	if((box[5]-box[4]) < min_thickness)
		box[5] = box[4] + min_thickness;

	return box;
}

// check if the faces share a side 
bool bvh::shareSide(Face f1, Face f2)
{
	bool ret = false;
	for(int i = 0 ; i < 3; ++i)
	{
		for(int j = 0; j <3; ++j)
		{
			ret = (ret || (f1[i]==f2[j]));
		}
	}
	return ret;
}

Vert bvh::findCenter(int p, int r)
{
	Eigen::Vector3d center(0,0,0);
	for(int i = p ; i <= r; ++i)
	{
		center += faceCenter(faces[i]);
	}
	return center/(double)(r-p+1);
}

Face bvh::getFace(int idx)
{
	return faces[idx];
}

