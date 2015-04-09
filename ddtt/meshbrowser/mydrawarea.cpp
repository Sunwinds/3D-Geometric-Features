#include "mydrawarea.h"
#include "surface_mesh/gl_wrappers.h"
#include <QMouseEvent>
#include <QDir>
#include <QFileInfo>
#include <stack>
#include "Octree.h"

using namespace SurfaceMesh;
#include "SurfaceMeshHelper.h"

QString curLabel = "";
QStringList AllLabels;
QVector<QColor> UniqueColors;
QVector<QString> labelNames; 

using namespace SurfaceMesh;

Q_DECLARE_METATYPE(std::vector< std::vector<double> >)


#pragma warning(disable:4309)
#pragma warning(disable:4267)

#include "weld.h"
inline void meregeVertices(SurfaceMesh::SurfaceMeshModel * m){
		std::vector<SurfaceMesh::Vector3> vertices;
		SurfaceMesh::Vector3VertexProperty points = m->vertex_coordinates();
		for(auto v: m->vertices()) vertices.push_back(points[v]);

		std::vector<size_t> xrefs;
		weld(vertices, xrefs, std::hash_Vector3d(), std::equal_to<Eigen::Vector3d>());

		std::vector< std::vector<SurfaceMesh::Vertex> > faces;
		for(auto f: m->faces()){
			std::vector<SurfaceMesh::Vertex> face;
			for(auto v: m->vertices(f)) face.push_back( SurfaceMesh::Vertex(xrefs[v.idx()]) );
			faces.push_back(face);
		}

		m->clear();

		for(auto v: vertices) m->add_vertex(v);
		for(auto face: faces) m->add_face(face);
}

MyDrawArea::MyDrawArea(SurfaceMesh::SurfaceMeshModel * mesh, QString filename) : m(mesh), filename(filename)
{
	isDeleted = false; 
	isDrawWireframe = true; 
	isDoubleLight = false;

	bg = QColor(51,51,51);
	fg = QColor(255,255,255);
}


void MyDrawArea::draw()
{
	SurfaceMesh::SurfaceMeshModel * mesh = m;

	if(!mesh) return;

	if(!mesh->property("hasNormals").toBool())
	{
		mesh->update_face_normals();
		mesh->update_vertex_normals();
		mesh->setProperty("hasNormals",true);
	}

	if(isDoubleLight){
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	}

	glClearColor( bg.redF(), bg.greenF(), bg.blueF(), bg.alphaF() );
	glClear(GL_COLOR_BUFFER_BIT);

	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LIGHTING);

	Surface_mesh::Vertex_property<Vector3> points = mesh->vertex_property<Vector3>("v:point");
	Surface_mesh::Face_property<Vector3> fnormals = mesh->face_property<Vector3>("f:normal");
	ScalarFaceProperty fsegment = mesh->face_property<Scalar>("f:segment", -1);

	Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;

	std::vector< std::vector<double> > colors = mesh->property("colors").value< std::vector< std::vector<double> > >();

	{
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
		glEnable( GL_POLYGON_OFFSET_FILL );
		glPolygonOffset( 1, 1 );

		glColor3f(fg.redF(),fg.greenF(),fg.blueF());

		glBegin(GL_TRIANGLES);
		for (fit=mesh->faces_begin(); fit!=fend; ++fit){
			glNormal3d( fnormals[fit][0], fnormals[fit][1], fnormals[fit][2] );
			fvit = fvend = mesh->vertices(fit);
			if(fsegment[fit] != -1){
				int i = fsegment[fit];
				glColor3f( colors[i][0], colors[i][1], colors[i][2] );
			}
			do{ glVertex3d(points[fvit][0],points[fvit][1],points[fvit][2]); } while (++fvit != fvend);
		}
		glEnd();

		if( isDrawWireframe )
		{
			glDisable( GL_POLYGON_OFFSET_FILL );
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
			{
				/// Render wire-frame
				glLineWidth(1);
				glColor4f(0,0,0,0.5);
				glDisable(GL_LIGHTING);
				glBegin(GL_TRIANGLES);
				for (fit=mesh->faces_begin(); fit!=fend; ++fit){
					glNormal3d( fnormals[fit][0], fnormals[fit][1], fnormals[fit][2] );
					fvit = fvend = mesh->vertices(fit);
					do{ glVertex3d(points[fvit][0],points[fvit][1],points[fvit][2]); } while (++fvit != fvend);
				}
				glEnd();
				glEnable(GL_LIGHTING);
			}
		}

        glDisable( GL_POLYGON_OFFSET_FILL );
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	}

	glDisable(GL_DEPTH_TEST);
	for(auto item: debugItems) item->draw(*this);
	glEnable(GL_DEPTH_TEST);

	if( isDeleted )
	{
		this->startScreenCoordinatesSystem();
		glLineWidth(5);

		glColor3f(1,0,0); //red
		glBegin(GL_LINES);
		gl::glVertex(Vector3(0,0,0));
		gl::glVertex(Vector3(this->width(), this->height(),0));
		glEnd();

		glLineWidth(2);
		this->stopScreenCoordinatesSystem();
    }

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    startScreenCoordinatesSystem();
    glColor3d(0,0,0);
	renderText(11,21, QFileInfo(filename).baseName());
	glColor3d(1,1,1);
	renderText(10,20, QFileInfo(filename).baseName());
    stopScreenCoordinatesSystem();
    glPopAttrib();

    glEnable( GL_MULTISAMPLE );
}

std::vector< SurfaceMesh::SurfaceMeshModel* > connectedPieces(SurfaceMesh::SurfaceMeshModel * mesh) {
	std::vector< SurfaceMesh::SurfaceMeshModel* > pieces;
	std::vector< std::vector< SurfaceMesh::Face > > pieceFaces;

	Vector3VertexProperty points = mesh->vertex_coordinates();
	BoolFaceProperty fvisisted = mesh->face_property<bool>("f:visisted", false);
	ScalarFaceProperty fsegment = mesh->face_property<Scalar>("f:segment", -1);
	for(auto f : mesh->faces()){
		fvisisted[f] = false;
		fsegment[f] = -1;
	}

	std::vector<int> sizes;

	for(Face f : mesh->faces())
	{
		if( fvisisted[f] ) continue;

		pieceFaces.push_back( std::vector< SurfaceMesh::Face >() );
		
		std::stack<Face> unprocessed;
		unprocessed.push(f);

		while( !unprocessed.empty() ) 
		{
			Face curFace = unprocessed.top();
			unprocessed.pop();
			if(fvisisted[curFace]) continue;

			fvisisted[curFace] = true;

			for(auto v : mesh->vertices(curFace))
			{ 
				for(auto fj : mesh->faces(v))
				{
					if(!mesh->is_valid(fj) || fvisisted[fj]) continue;
					unprocessed.push(fj);
				}
			}

			fsegment[curFace] = (pieceFaces.size() - 1);

			pieceFaces.back().push_back( curFace );
		}

		SurfaceMesh::SurfaceMeshModel * piece = new SurfaceMesh::SurfaceMeshModel;
		piece->reserve( pieceFaces.back().size() * 3, pieceFaces.back().size() * 6, pieceFaces.back().size() );

		for(auto v : mesh->vertices()) 
		{
			piece->add_vertex( points[v] );
		}

		for(auto f : pieceFaces.back()){
			std::vector<Vertex> face;
			for(auto v : mesh->vertices(f)) face.push_back(v);
			piece->add_face(face);
		}

		for(auto v : piece->vertices()) if(piece->is_isolated(v)) piece->remove_vertex(v);
		piece->garbage_collection();

		sizes.push_back(piece->n_vertices());

		pieces.push_back( piece );
	}

	// Assign colors
	{
		std::set<int> segments;
		for(auto f : mesh->faces()) segments.insert(fsegment[f]);
		segments.erase(-1);
		std::vector< std::vector<double> > colors = starlab::randomColors(segments.size());
		QVariant colorsVar;
		colorsVar.setValue(colors);
		mesh->setProperty("colors", colorsVar);
	}

	return pieces;
}

struct IsectData{
	double u,v,t;
};

bool rayTriangleIntersectionTest( const Ray &ray, Vector3 v0, Vector3 v1, Vector3 v2, HitResult & res, bool allowBack = true )
{
	res.hit = false;
	res.distance = DBL_MAX;

	double EPS = 1e-7;

	Eigen::Vector3d vertex1 = v0;
	Eigen::Vector3d vertex2 = v1;
	Eigen::Vector3d vertex3 = v2;

	// Compute vectors along two edges of the triangle.
	Eigen::Vector3d edge1 = vertex2 - vertex1;
	Eigen::Vector3d edge2 = vertex3 - vertex1;

	// Compute the determinant.
	Eigen::Vector3d directionCrossEdge2 = cross(ray.direction, edge2);

	double determinant = dot(edge1, directionCrossEdge2);

	// If the ray is parallel to the triangle plane, there is no collision.
	if (fabs(determinant) < EPS)
		return false;

	double inverseDeterminant = 1.0 / determinant;

	// Calculate the U parameter of the intersection point.
	Eigen::Vector3d distVector = ray.origin - vertex1;
	double triangleU = dot(distVector, directionCrossEdge2);
	triangleU *= inverseDeterminant;

	// Make sure it is inside the triangle.
	if (triangleU < 0 - EPS || triangleU > 1 + EPS)
		return false;

	// Calculate the V parameter of the intersection point.
	Eigen::Vector3d distanceCrossEdge1 = cross(distVector, edge1);
	double triangleV = dot(ray.direction, distanceCrossEdge1);
	triangleV *= inverseDeterminant;

	// Make sure it is inside the triangle.
	if (triangleV < 0 - EPS || triangleU + triangleV > 1 + EPS)
		return false;

	// Compute the distance along the ray to the triangle.
	double rayDistance = dot(edge2, distanceCrossEdge1);
	rayDistance *= inverseDeterminant;

	if(!allowBack){
		// Is the triangle behind the ray origin?
		if (rayDistance < 0)
			return false;
	}

	res.hit = true;
	res.distance = rayDistance;

	res.u = triangleU;
	res.v = triangleV;
	return true;
}


void MyDrawArea::mousePressEvent(QMouseEvent * event)
{
	QGLViewer::mousePressEvent(event);

    if(!(event->buttons() & Qt::RightButton)) return;

    Surface_mesh::Vertex_property<Vector3> points = m->vertex_property<Vector3>("v:point");

	if(curOp == MeshOperation::UNDO_ALL){
		m->clear();
		for(auto p: originalVertices) m->add_vertex(p);
		for(auto face: originalFaces) m->add_face(face);
	}

	if(curOp == MeshOperation::UNDO){
		m->clear();
		for(auto p: oldVertices) m->add_vertex(p);
		for(auto face: oldFaces) m->add_face(face);
	}else{
		oldVertices.clear();
		oldFaces.clear();

		for(auto v : m->vertices()) oldVertices.push_back(points[v]);
		for(auto f: m->faces()){
			std::vector<SurfaceMesh::Vertex> face;
			for(auto v: m->vertices(f)) face.push_back(v);
			oldFaces.push_back(face);
		}
	}

    if(curOp == MeshOperation::ROTATE_LEFT){
        Eigen::AngleAxisd rot(0.5 * M_PI, Vector3::UnitZ());
        for(auto v : m->vertices())
            points[v] = rot * points[v];
    }

    if(curOp == MeshOperation::ROTATE_RIGHT){
        Eigen::AngleAxisd rot(-0.5 * M_PI, Vector3::UnitY());
        for(auto v : m->vertices())
            points[v] = rot * points[v];
    }

	if(curOp == MeshOperation::ROTATE_UP){
		Eigen::AngleAxisd rot(-0.5 * M_PI, Vector3::UnitX());
		for(auto v : m->vertices())
			points[v] = rot * points[v];
	}

	if(curOp == MeshOperation::FLIP){
		std::vector<Vector3> verts;

		// Store coordinates
		Vector3VertexProperty points = m->vertex_coordinates();
		for(auto v: m->vertices()) verts.push_back( points[v] );

		// Store reversed faces
		std::vector< std::vector<SurfaceMesh::Vertex> > faces;
		for(auto f: m->faces()){
			std::vector<SurfaceMesh::Vertex> face;
			for(auto v: m->vertices(f)) face.push_back(v);
			std::reverse(face.begin(), face.end());
			faces.push_back(face);
		}

		// Remove everything
		m->clear();
		m->garbage_collection();

		// Reconstruct
		for(auto v: verts) m->add_vertex(v);
		for(auto f: faces) m->add_face(f);
	}

	if(curOp == MeshOperation::INVERT_PART)
	{

	}

	if(curOp == MeshOperation::REMOVE_SELECTED)
	{
		qglviewer::Vec _orig, _dir;
		camera()->convertClickToLine( event->pos(), _orig, _dir );
		Vector3 orig(_orig[0],_orig[1],_orig[2]);
		Vector3 dir(_dir[0],_dir[1],_dir[2]);
		Ray r(orig,dir);
		std::vector< SurfaceMesh::SurfaceMeshModel* > allpieces = connectedPieces(m);

		QSet<int> used;
		for(size_t i = 0; i < allpieces.size(); i++)
		{
			Vector3VertexProperty points = allpieces[i]->vertex_coordinates();

			for(auto f : allpieces[i]->faces()){
				std::vector<Vector3> v;
				for( auto vert: allpieces[i]->vertices(f) ) 
					v.push_back(points[vert]);

				HitResult hit;
				if( rayTriangleIntersectionTest(r, v[0], v[1], v[2], hit) ) 
					break;

				used.insert(i);
			}
		}

		m->clear();

		std::vector< SurfaceMesh::SurfaceMeshModel* > pieces;
		for(auto key : used) pieces.push_back( allpieces[key] );

		for(auto piece: pieces)
		{
			Vector3VertexProperty points = piece->vertex_coordinates();

			int offset = m->n_vertices();

			for(auto v: piece->vertices()){
				m->add_vertex( points[v] );
			}

			for(auto f: piece->faces()){
				std::vector<SurfaceMesh::Vertex> face;
				for( auto v: piece->vertices(f) ) 
					face.push_back( Vertex(offset + v.idx()) );
				m->add_face( face );
			}
		}

		for(auto piece: allpieces) delete piece;
	}

	if(curOp == MeshOperation::COMPONENTS)
	{
		std::vector< SurfaceMesh::SurfaceMeshModel* > pieces = connectedPieces(m);
		for(auto piece: pieces) delete piece;
	}

	if(curOp == MeshOperation::REMOVE_SMALL || curOp == MeshOperation::REMOVE_LARGE)
	{
		std::vector< SurfaceMesh::SurfaceMeshModel* > allpieces = connectedPieces(m);

		QMap<int, double> size;
		for(size_t i = 0; i < allpieces.size(); i++){
			ScalarFaceProperty fareas = SurfaceMeshHelper(allpieces[i]).computeFaceAreas();

			double area = 0;
			for(auto f: allpieces[i]->faces()) area += fareas[f];
			size[i] = area;
		}

		QList<double> values = size.values();
		int largest = std::max_element(values.begin(), values.end()) - values.begin();
		int smallest = std::min_element(values.begin(), values.end()) - values.begin();
		if(curOp == MeshOperation::REMOVE_SMALL) size.remove(smallest);
		if(curOp == MeshOperation::REMOVE_LARGE) size.remove(largest);

		m->clear();

		std::vector< SurfaceMesh::SurfaceMeshModel* > pieces;
		for(auto key : size.keys()) pieces.push_back( allpieces[key] );

		std::vector<int> sizes;
		for(auto piece: pieces)	sizes.push_back( piece->n_vertices() );

		for(auto piece: pieces)
		{
			Vector3VertexProperty points = piece->vertex_coordinates();

			int offset = m->n_vertices();

			for(auto v: piece->vertices()){
				m->add_vertex( points[v] );
			}

			for(auto f: piece->faces()){
				std::vector<SurfaceMesh::Vertex> face;
				for( auto v: piece->vertices(f) ) 
					face.push_back( Vertex(offset + v.idx()) );
				m->add_face( face );
			}
		}

		for(auto piece: allpieces) delete piece;
	}

	if(curOp == MeshOperation::CLOSE_HOLES)
	{
		Surface_mesh::Halfedge_property<bool> hvisisted = m->halfedge_property<bool>("h:visisted", false);

		std::vector< std::vector<Halfedge> > holes;

		for(Halfedge h : m->halfedges()){
			if( !m->is_boundary(h) || hvisisted[h] ) continue;

			std::vector<Halfedge> hole;

			Halfedge hinit = h;
			h = m->next_halfedge(h);

			while( h != hinit ) {
				hole.push_back( h );
				hvisisted[h] = true;
				h = m->next_halfedge(h);
			}
			hole.push_back( h );

			holes.push_back(hole);
		}

		Vector3VertexProperty points = m->vertex_coordinates();

		for(std::vector<Halfedge> hole : holes){
			std::vector<Vertex> holeVerts;
			for(auto h : hole){
				Vertex v = m->to_vertex(h);
				holeVerts.push_back(v);
			}

			// Compute hole center
			Vector3 center(0,0,0);
			for(auto v: holeVerts) center += points[v];
			center /= holeVerts.size();

			Vertex c = m->add_vertex(center);

			for(size_t vi = 0; vi < holeVerts.size(); vi++){
				std::vector<Vertex> face;

				face.push_back( holeVerts[vi] );
				face.push_back( holeVerts[(vi+1) % holeVerts.size()] );
				face.push_back( c );

				m->add_face(face);
			}
		}
	}

	if(curOp == MeshOperation::MERGE_VERTICES)
	{
		meregeVertices(m);
	}

	if(curOp == MeshOperation::UNIFY)
	{
		Vector3VertexProperty points = m->vertex_coordinates();
		Vector3FaceProperty fnormals = m->face_normals();

		Octree * octree = m->property("octree").value<Octree*>();
		if(!octree){
			octree = new Octree( m );
			QVariant octreeVal;
			octreeVal.setValue(octree);
			m->setProperty("octree", octreeVal);
		}

		Vector3 cameraPos, cameraDirection;
		qglviewer::Vec cpos = camera()->position();
		qglviewer::Vec cdir = camera()->viewDirection();
		for(int i = 0; i < 3; i++){
			cameraPos[i] = cpos[i];
			cameraDirection[i] = cdir[i];
		}

		Ray r(cameraPos, cameraDirection);

		/*double tiny = m->bbox().sizes().norm() * 1e-6;

		for(Face f : m->faces())
		{
			// Test triangles
			double dot = fnormals[f].dot(cameraDirection.normalized());
			if(dot < -0.4){
				QVector<starlab::QVector3> verts;
				int c = 0;
				for(auto v: m->vertices(f)){
					verts.push_back(points[v]);
					Ray r( points[v] + (-cameraDirection * tiny), -cameraDirection);
					QSet<int> tris = octree->intersectRay(r, r.thickness, true);
					if(tris.size() < 1)	c++;
				}
				if(c != 3) continue;
			}

			// Add to collection to try to flip

		}*/
		starlab::PolygonSoup * ps = new starlab::PolygonSoup;

		int step = 5;
		int w = width();
		int h = height();
		int nrows = std::floor( double(h) / step );
		int ncols = std::floor( double(w) / step );

		std::set<int> tris;

		for(int winX=0,i=0; i<ncols; winX+=step, i++){
			for(int winY=0,j=0; j<nrows; winY+=step, j++){
				qglviewer::Vec _orig, _dir;
				camera()->convertClickToLine( QPoint(winX, winY), _orig, _dir );
				Vector3 orig(_orig[0],_orig[1],_orig[2]);
				Vector3 dir(_dir[0],_dir[1],_dir[2]);
				dir.normalize(); ///< just to be sure
				int isectHit = -1;
				Eigen::Vector3d ipoint = octree->closestIntersectionPoint( Ray(orig, dir), &isectHit, false );
				if( isectHit >= 0 ) tris.insert(isectHit);
			}
		}

		// Collect inverted triangles
		std::set<int> invertedTris;
		for(int fidx: tris){
			Face f(fidx);
			double dot = fnormals[f].dot(cameraDirection.normalized());

			if(dot < 0) continue;

			QVector<starlab::QVector3> verts;
			for(auto v: m->vertices(f)) verts.push_back(points[v]);
			//ps->addPoly(verts);
			invertedTris.insert( fidx );
		}
		debugItems.clear();
		debugItems.push_back(ps);

		// Classify regular and bad faces
		std::vector< std::vector<Vertex> > regularFaces, badFaces;
		for(auto f : m->faces()){
			if( invertedTris.count(f.idx()) != 0){
				std::vector<Vertex> face;
				for(auto v: m->vertices(f)) face.push_back(v);
				std::reverse(face.begin(), face.end());
				badFaces.push_back(face);
				continue;
			}
			std::vector<Vertex> face;
			for(auto v: m->vertices(f)) face.push_back(v);
			regularFaces.push_back(face);
		}

		// Add vertices
		SurfaceMeshModel newMesh;
		for(auto v : m->vertices()) newMesh.add_vertex(points[v]);
		int n_vertices = newMesh.n_vertices();

		// Duplicate vertices
		for(auto v : m->vertices()) newMesh.add_vertex(points[v]);

		// Add regular faces
		for(auto f : regularFaces) newMesh.add_face(f);

		// Add inverted bad faces
		for(auto f : badFaces){
			for(int i = 0; i < 3; i++) f[i] = Vertex(n_vertices + f[i].idx());
			newMesh.add_face(f);
		}
	
		// Remove isolated vertices
		for(auto v : newMesh.vertices()) if(newMesh.is_isolated(v))	newMesh.remove_vertex(v);
		newMesh.garbage_collection();

		m->clear();

		Vector3VertexProperty newPoints = newMesh.vertex_coordinates();
		for(auto v : newMesh.vertices()) m->add_vertex( newPoints[v] );
		for(auto f : newMesh.faces()){
			std::vector<Vertex> face;
			for(auto v : newMesh.vertices(f)) face.push_back(v);
			m->add_face(face);
		}

		m->setProperty("octree", NULL);
	}

    cleanUp(m);
}

void MyDrawArea::postSelection(const QPoint&)
{

}

static void saveOBJ(SurfaceMesh::Model * mesh, QString filename)
{
	QFile file(filename);

	// Create folder
	QFileInfo fileInfo(file.fileName());
	QDir d(""); d.mkpath(fileInfo.absolutePath());

	// Open for writing
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;

	QTextStream out(&file);
	out << "# Exported from MeshBrowser" << "\n";
	out << "# NV = " << mesh->n_vertices() << " NF = " << mesh->n_faces() << "\n";
	SurfaceMesh::Vector3VertexProperty points = mesh->vertex_property<Eigen::Vector3d>("v:point");
	foreach( SurfaceMesh::Vertex v, mesh->vertices() )
		out << "v " << points[v][0] << " " << points[v][1] << " " << points[v][2] << "\n";
	foreach( SurfaceMesh::Face f, mesh->faces() ){
		out << "f ";
		Surface_mesh::Vertex_around_face_circulator fvit=mesh->vertices(f), fvend=fvit;
		do{     out << (((Surface_mesh::Vertex)fvit).idx()+1) << " ";} while (++fvit != fvend);
		out << "\n";
	}
	file.close();
}

MyDrawArea::~MyDrawArea()
{
	if(isDeleted) return;

	saveOBJ(m, filename);
}

void MyDrawArea::focusInEvent(QFocusEvent * event)
{
	QGLViewer::focusInEvent(event);

	emit( gotFocus(this) );
}
