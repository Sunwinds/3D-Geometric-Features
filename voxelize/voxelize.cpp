#include <QElapsedTimer>
#include <QFileDialog>

#include "RenderObjectExt.h"
#include "SurfaceMeshHelper.h"

#include "voxelize.h"
#include "voxelization.h"

#include "mc.h"

#include "NanoKdTree.h"
inline void snapCloseVertices( std::vector<SurfaceMesh::Vector3> & vertices, double threshold ){
	NanoKdTree tree;
	for(auto p : vertices) tree.addPoint(p);
	tree.build();
	for(auto & p : vertices){
		KDResults matches;
		tree.ball_search(p, threshold, matches);
		for(auto m : matches) vertices[m.first] = p;
	}
} 

inline bool CompareVector3(const Vector3& p, const Vector3& q){ 
	if(p.x() == q.x()){
		if(p.y() == q.y()) return p.z() < q.z();	
		return p.y() < q.y();
	}
	return p.x() < q.x();
}

inline void meregeVertices(SurfaceMesh::SurfaceMeshModel * m)
{
	// Collect original vertices
	std::vector<SurfaceMesh::Vector3> original;
	SurfaceMesh::Vector3VertexProperty points = m->vertex_coordinates();
	for( auto v : m->vertices() ) original.push_back( points[v] );

	// Snap close vertices, sort them, then remove duplicates
	snapCloseVertices( original, 1e-12 );	
	std::vector<SurfaceMesh::Vector3> clean = original;
	std::sort( clean.begin(), clean.end(), CompareVector3);
	clean.erase( std::unique(clean.begin(), clean.end()), clean.end() );

	// Find new ids
	std::vector<size_t> xrefs( original.size(), 0 );
	for (size_t i = 0; i != original.size(); i += 1)
		xrefs[i] = std::lower_bound(clean.begin(), clean.end(), original[i], CompareVector3) - clean.begin();
	
	// Replace face vertices
	std::vector< std::vector<SurfaceMesh::Vertex> > faces;
	for(auto f: m->faces()){
		std::vector<SurfaceMesh::Vertex> faceverts;
		for(auto v: m->vertices(f)) faceverts.push_back( SurfaceMesh::Vertex(xrefs[ v.idx() ]) );
		faceverts.erase( std::unique(faceverts.begin(), faceverts.end()), faceverts.end() );
		if(faceverts.size() == 3) faces.push_back(faceverts); // skip degenerate faces
	}

	// Rebuild
	m->clear();
	for(auto v: clean) m->add_vertex(v);
	for(auto face: faces) m->add_face(face);
}

inline Vector3 findNearestPoint(SurfaceMeshModel * mesh, SurfaceMeshModel * original_mesh, 
								NanoKdTree * kdtree, const Vector3& _point, SurfaceMeshModel::Face& _fh, double* _dbest)
{
	Vector3VertexProperty orig_points = original_mesh->vertex_property<Vector3>( VPOINT );

	double fc = original_mesh->bbox().diagonal().norm() * 2;
	Vector3  p_best = Vector3(fc,fc,fc) + Vector3(original_mesh->bbox().center());
	SurfaceMeshModel::Scalar d_best = (_point - p_best).squaredNorm();
	SurfaceMeshModel::Face fh_best;

	if( true )
	{
		foreach(Face f, original_mesh->faces())
		{
			Surface_mesh::Vertex_around_face_circulator cfv_it = original_mesh->vertices(f);

			// Assume triangular
			const Vector3& pt0 = orig_points[   cfv_it];
			const Vector3& pt1 = orig_points[ ++cfv_it];
			const Vector3& pt2 = orig_points[ ++cfv_it];

			Vector3 ptn = _point;

			//SurfaceMeshModel::Scalar d = distPointTriangleSquared( _point, pt0, pt1, pt2, ptn );
			SurfaceMeshModel::Scalar d = ClosestPointTriangle( _point, pt0, pt1, pt2, ptn );

			if( d < d_best)
			{
				d_best = d;
				p_best = ptn;

				fh_best = f;
			}
		}
	} else {
		KDResults matches;
		kdtree->k_closest(_point, 32, matches);

		foreach(KDResultPair match, matches){
			foreach(Halfedge h, original_mesh->onering_hedges(Vertex((int)match.first))){
				Face f = original_mesh->face(h);
				if(!mesh->is_valid(f)) continue;
				SurfaceMeshModel::Vertex_around_face_circulator cfv_it = original_mesh->vertices( f );

				// Assume triangular
				const Vector3& pt0 = orig_points[   cfv_it],pt1 = orig_points[ ++cfv_it],pt2 = orig_points[ ++cfv_it];

				Vector3 ptn = _point;
				SurfaceMeshModel::Scalar d = ClosestPointTriangle( _point, pt0, pt1, pt2, ptn );
				if( d < d_best){
					d_best = d;
					p_best = ptn;
					fh_best = f;
				}
			}
		}
	}

	_fh = fh_best;
	if(_dbest) *_dbest = sqrt(d_best);
	return p_best;
}

inline void projectToSurface( SurfaceMeshModel * mesh, SurfaceMeshModel * orginal_mesh, double maxDistance )
{
	SurfaceMeshModel::Vertex_iterator v_it;
	SurfaceMeshModel::Vertex_iterator v_end = mesh->vertices_end();
	Vector3VertexProperty points = mesh->vertex_coordinates();
	Vector3VertexProperty original_points = orginal_mesh->vertex_coordinates();

	// Build KD-tree
	NanoKdTree kdtree;
	foreach(Vertex v, mesh->vertices()) kdtree.addPoint(original_points[v]);
	kdtree.build();

	for (v_it = mesh->vertices_begin(); v_it != v_end; ++v_it){
		SurfaceMeshModel::Face fhNear;
		double distance;
		Vector3 pNear = findNearestPoint(mesh, orginal_mesh, &kdtree, points[v_it], fhNear, &distance);
		if( distance < maxDistance ) points[v_it] = pNear;
	}
}

void voxelize::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichInt("Resolution", 64, "Resolution"));
	pars->addParam(new RichBool("Solid", false, "Solid"));
	pars->addParam(new RichBool("Manifold", false, "Manifold"));
	pars->addParam(new RichBool("Project", false, "Project"));
	pars->addParam(new RichBool("Marching", false, "Marching"));
	pars->addParam(new RichInt("Smooth", 0, "Smooth"));
	pars->addParam(new RichBool("Intersection", false, "Intersection"));
	pars->addParam(new RichInt("Operation", 0, "Operation"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

void voxelize::applyFilter(RichParameterSet *pars)
{
	drawArea()->clear();

	size_t gridsize = pars->getInt("Resolution");

	if(pars->getBool("Intersection"))
	{
		QList<Starlab::Model*> models = mainWindow()->document()->models();

		if(models.size() < 2){
			foreach( QString filename, QFileDialog::getOpenFileNames(0, "Open Meshes", 
				mainWindow()->settings()->getString("lastUsedDirectory"), "Mesh Files (*.obj)") )
			{
				SurfaceMeshModel * m = new SurfaceMeshModel(filename, QFileInfo(filename).baseName());
				m->read( filename.toStdString() );
				m->updateBoundingBox();
				m->update_face_normals();
				m->update_vertex_normals();

				document()->addModel(m);
			}

			models = mainWindow()->document()->models();
		}

		SurfaceMeshModel * modelA = dynamic_cast<SurfaceMeshModel*>(models.at(models.size() - 2));
		SurfaceMeshModel * modelB = dynamic_cast<SurfaceMeshModel*>(models.at(models.size() - 1));

		double unitlength = 0;
		vector<VoxelData> data = ComputeVoxelizationCSG(modelA, modelB, unitlength, gridsize, BooleanOperation(pars->getInt("Operation")));

		starlab::CubeSoup * cs = new starlab::CubeSoup;
		starlab::PointSoup * ps = new starlab::PointSoup(1);

		for(auto voxel : data)
		{
			unsigned int x,y,z;
			mortonDecode(voxel.morton, z, y, x);
			Vector3 p = Vector3(Vector3(x,y,z) * unitlength) + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );

			if(pars->getBool("Visualize"))
			{
				if(gridsize < 64)
					cs->addCube(p, unitlength);
				else
					ps->addPoint(p, Qt::yellow);
			}
		}

		drawArea()->addRenderObject( cs );
		drawArea()->addRenderObject( ps );
	}
	else
	{
		QElapsedTimer timer; timer.start();

		VoxelContainer voxels = ComputeVoxelization( mesh(), gridsize, pars->getBool("Solid"), pars->getBool("Manifold") );
		double unitlength = voxels.unitlength;

		mainWindow()->setStatusBarMessage( QString("Time (%1 ms) / Count (%2 voxels)").arg(timer.elapsed()).arg(voxels.data.size()));

		// Build mesh
		if( pars->getBool("Manifold") )
		{
			SurfaceMeshModel * newMesh = new SurfaceMeshModel("voxelized_mesh.obj", "voxelized");
			
			int voffset = 0;
			for(auto q : voxels.quads)
			{
				std::vector<Vertex> quad_verts;
				for(int i = 0; i < 4; i++){
					newMesh->add_vertex( q[i] );
					quad_verts.push_back( Vertex(voffset + i) );
				}
				newMesh->add_face( quad_verts );
				voffset += 4;
			}

			newMesh->garbage_collection();
			newMesh->triangulate();
			meregeVertices( newMesh );

			// Smooth if requested
			SurfaceMeshHelper h(newMesh);
			h.smoothVertexProperty<Vector3>(VPOINT, pars->getInt("Smooth"), Vector3(0,0,0));

			// remove exciting voxelization
			SurfaceMeshModel * oldMesh = (SurfaceMeshModel *) document()->getModel( newMesh->name );
			if( oldMesh ) document()->deleteModel( oldMesh );

			// Prepare then add to document
			newMesh->updateBoundingBox();
			newMesh->update_face_normals();
			newMesh->update_vertex_normals();

			document()->addModel( newMesh );

			// Report any holes
			bool hasHoles = false;
			Vector3VertexProperty points = newMesh->vertex_coordinates();
			for( auto v : newMesh->vertices() ){
				if(newMesh->is_boundary(v)) {
					hasHoles = true;
					drawArea()->drawPoint(points[v], 10);
				}
			}
			if( hasHoles ) mainWindow()->setStatusBarMessage( QString("Holes are found in mesh") );

			if(pars->getBool("Project"))
			{
				projectToSurface( newMesh, mesh(), voxels.unitlength * 5 );
			}
		}

		if(pars->getBool("Marching"))
		{
			ScalarVolume volume = initScalarVolume( gridsize, 1.0 );

			for(auto voxel : voxels.data)
			{
				unsigned int x,y,z;
				mortonDecode(voxel.morton, x, y, z);

				volume[x+MC_VOLUME_PADDING][y+MC_VOLUME_PADDING][z+MC_VOLUME_PADDING] = -1;
			}

			// Smooth volume
			int smoothIterations = pars->getInt("Smooth");
			ScalarVolume smoothVolume = volume;
			for(int i = 0; i < smoothIterations; i++)
			{
				for(int z = 1; z < volume.size()-1; z++){
					for(int y = 1; y < volume[z].size()-1; y++){
						for(int x = 1; x < volume[z][y].size()-1; x++)
						{
							double sum = 0;
							int count = 0;

							for(int u = -1; u <= 1; u++){
								for(int v = -1; v <= 1; v++){
									for(int w = -1; w <= 1; w++){
										sum += volume[x+u][y+v][z+w];
										count++;
									}
								}
							}

							smoothVolume[x][y][z] = sum / count;
						}
					}
				}
			}
			volume = smoothVolume;

			SurfaceMeshModel * newMesh = new SurfaceMeshModel("marching_cubes.obj", "marching_cubes");

			int voffset = 0;
			for( auto tri : march( volume, 0.5 ) ){
				std::vector<Vertex> quad_verts;
				for(int i = 0; i < 3; i++){
					Vector3 p = Vector3(tri[i].x,tri[i].y,tri[i].z);

					p = (p * voxels.unitlength) + voxels.translation + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );

					newMesh->add_vertex( p );
					quad_verts.push_back( Vertex(voffset + i) );
				}
				newMesh->add_face( quad_verts );
				voffset += 3;
			}

			meregeVertices( newMesh );

			// Smooth if requested
			//SurfaceMeshHelper h(newMesh);
			//h.smoothVertexProperty<Vector3>(VPOINT, pars->getInt("Smooth"), Vector3(0,0,0));

			// remove any exciting
			SurfaceMeshModel * oldMesh = (SurfaceMeshModel *) document()->getModel( newMesh->name );
			if( oldMesh ) document()->deleteModel( oldMesh );

			// Prepare then add to document
			newMesh->updateBoundingBox();
			newMesh->update_face_normals();
			newMesh->update_vertex_normals();

			document()->addModel( newMesh );
		}

		if(pars->getBool("Visualize"))
		{
			if( false )
			{
				starlab::CubeSoup * cs = new starlab::CubeSoup(1.0, true);
				std::vector<Vector3> centers;

				for(auto voxel : voxels.data){
					unsigned int x,y,z;
					mortonDecode(voxel.morton, z, y, x);
					Vector3 p = voxels.translation + Vector3(Vector3(x,y,z) * unitlength) + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );
					centers.push_back( p );
				}

				for(auto p : centers) cs->addCube( p, unitlength );

				drawArea()->addRenderObject( cs );
			}
			
			if( false )
			{
				starlab::PolygonSoup * ps = new starlab::PolygonSoup;
				for(auto quad : voxels.quads){
					QVector<starlab::QVector3> pnts;
					for(auto p : quad) pnts.push_back(p);
					ps->addPoly(pnts);
				}

				drawArea()->addRenderObject( ps );
			}

			if( true )
			{
				// Aux. voxels
				{
					starlab::CubeSoup * cs = new starlab::CubeSoup(1.0, false);
					std::vector<Vector3> centers;
					for(auto voxel : voxels.aux){
						unsigned int x,y,z;
						mortonDecode(voxel.morton, z, y, x);
						Vector3 p = voxels.translation + Vector3(Vector3(x,y,z) * unitlength) + ( 0.5 * Vector3(unitlength,unitlength,unitlength) );
						
						cs->addCube( p, unitlength * 0.5, QColor::fromRgbF(voxel.normal.x(),voxel.normal.y(),voxel.normal.z(), 0.2) );
					}
					drawArea()->addRenderObject( cs );
				}
			}
		}
	}
}
