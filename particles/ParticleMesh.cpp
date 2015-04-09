#include "ParticleMesh.h"

#include "voxelization.h"

#include "RenderObjectExt.h"
#include <QGLWidget>

inline QVector<QColor> rndColors(int count){
	QVector<QColor> c;
	for(int i = 0; i < count; i++) c << starlab::qRandomColor3();
	return c;
}

QVector<QColor> ParticleMesh::rndcolors = rndColors(512);

ParticleMesh::ParticleMesh(SurfaceMeshModel * mesh, int gridsize, double particle_raidus) : surface_mesh(NULL),
	raidus(particle_raidus), tranlsation(Eigen::Vector3d(0,0,0))
{
	// Voxelization
	grid = ComputeVoxelization<VoxelVector>(mesh, gridsize, true, true);

	// Build voxelization mesh
	QString objectName = mesh->name;
	surface_mesh = new SurfaceMeshModel(objectName + ".obj", objectName);
	{
		int voffset = 0;
		for(auto q : grid.quads){
			std::vector<Vertex> quad_verts;
			for(int i = 0; i < 4; i++){
				surface_mesh->add_vertex( q[i].cast<double>() );
				quad_verts.push_back( Vertex(voffset + i) );
			}
			surface_mesh->add_face( quad_verts );
			voffset += 4;
		}

		surface_mesh->triangulate();
		meregeVertices( surface_mesh );

		for(auto v : surface_mesh->vertices())
			if(surface_mesh->is_isolated(v))
				surface_mesh->remove_vertex(v);

		surface_mesh->garbage_collection();
	}

	// Remove outer most voxel
	//container.data.erase(std::remove_if(container.data.begin(), container.data.end(), 
	//	[](const VoxelData<Eigen::Vector3f> & vd) { return !vd.isOuter; }), container.data.end());

	// Insert particles
    for( auto voxel : grid.data )
    {
		Eigen::Vector3f point = grid.voxelPos(voxel.morton);

        Particle particle( point.cast<double>() );

        particle.id = particles.size();
		particle.morton = voxel.morton;
		mortonToParticleID[voxel.morton] = particle.id;

        particles.push_back( particle );
		bbox.extend( point.cast<double>() );
    }

	grid.findOccupied();

	// KD-tree
	{
		kdtree = new NanoKdTree;
		Vector3 sizes = bbox.sizes();
		for( auto & particle : particles )
		{
			Vector3 mapped = (particle.pos - bbox.min());
			for(int i = 0; i < 3; i++) mapped[i] /= sizes[i];

			particle.relativePos = mapped;
			kdtree->addPoint( particle.relativePos );
		}
		kdtree->build();
	}

	process();
}

void ParticleMesh::process()
{
	double bmin = bbox.min().z(), bmax = bbox.max().z();

	for(auto & particle : particles)
	{
		particle.measure = (particle.pos.z() - bmin) / (bmax - bmin);
	}
}

void ParticleMesh::drawParticles()
{
	glDisable( GL_LIGHTING );
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask(GL_FALSE);

	glPointSize(4);
	glBegin(GL_POINTS);

	for(auto particle : particles)
	{
		//QColor c = starlab::qtJetColor(particle.measure);
		//Eigen::Vector4d color( c.redF(), c.greenF(), c.blueF(), particle.alpha );
		//Eigen::Vector4d color(abs(particle.direction[0]), abs(particle.direction[1]), abs(particle.direction[2]), particle.alpha);
		//color[0] *= color[0];color[1] *= color[1];color[2] *= color[2];
		//glColor4dv( color.data() );
		
		QColor c = rndcolors[ particle.flag ];
		Eigen::Vector4d color(c.redF(),c.greenF(),c.blueF(), particle.alpha);
		glColor4dv( color.data() );
		glVertex3dv( particle.pos.data() );
	}

	glEnd();

	glDepthMask(GL_TRUE);
}

void ParticleMesh::drawDebug(QGLWidget & widget)
{
	for(auto d : debug) d->draw( widget );
}

ParticleMesh::~ParticleMesh()
{
	if(surface_mesh) delete surface_mesh;
	if(kdtree) delete kdtree;
}
