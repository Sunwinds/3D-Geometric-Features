#pragma once
#include <QWidget>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include "Particle.h"
#include "NanoKdTree.h"
#include "RenderObjectExt.h"

#include "voxelization.h"

class ParticleMesh
{
public:
	ParticleMesh(SurfaceMeshModel * mesh, int gridsize = 64, double particle_raidus = 0.1);
	~ParticleMesh();

	std::vector<Eigen::Vector3d> extractSurface();

	void process();

	void drawParticles();
	void drawDebug(QGLWidget & widget);
	static QVector<QColor> rndcolors;

	typedef Eigen::Vector3f VoxelVector;
	VoxelContainer<VoxelVector> grid;
	std::map<uint64_t,size_t> mortonToParticleID;

	std::vector< std::vector<float> > desc;

    std::vector<Particle> particles;
    double raidus;
	Eigen::Vector3d tranlsation;
	Eigen::AlignedBox3d bbox;

	SurfaceMeshModel * surface_mesh;
	NanoKdTree * kdtree;

	QVector<RenderObject::Base*> debug;
};
