#ifndef PRIMITIVES_TRIMESH_HPP
#define PRIMITIVES_TRIMESH_HPP

#include <vector>
#include "shapes/particle.hpp"

typedef Particle* Vertex;
typedef Eigen::Vector2i Edge;
typedef Eigen::Vector3i Face;

class TriMesh
{
public:
	std::vector<Vertex> vertices;
	std::vector<Edge> edges;
	std::vector<Face> faces;
};

#endif