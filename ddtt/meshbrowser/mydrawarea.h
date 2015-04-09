#ifndef MYDRAWAREA_H
#define MYDRAWAREA_H

#include <qglviewer/qglviewer.h>
#include "SurfaceMeshModel.h"
#include "RenderObjectExt.h"

extern QString curLabel;
extern QStringList AllLabels;
extern QVector<QColor> UniqueColors;
extern QVector<QString> labelNames;

static QString meshOps[] = {"none", "rotateleft", "rotateright", 
	"rotateup", "flip", "invertpart", "unifyorientation", "removesmallest", "removelargest", 
	"removeselected", "closeholes", "mergevertices", "UNDO", "undoAll", "COMPONENTS"};

enum MeshOperation{ NONE_OP, ROTATE_LEFT, ROTATE_RIGHT, ROTATE_UP, FLIP, INVERT_PART, UNIFY,
	REMOVE_SMALL, REMOVE_LARGE, REMOVE_SELECTED, CLOSE_HOLES, MERGE_VERTICES, UNDO, UNDO_ALL, COMPONENTS };

extern MeshOperation curOp;

class MyDrawArea : public QGLViewer
{
	Q_OBJECT
public:
    MyDrawArea(SurfaceMesh::SurfaceMeshModel * mesh, QString filename);
	~MyDrawArea();

	void draw();

	void postSelection(const QPoint& point);
	void mousePressEvent(QMouseEvent * event);
	void focusInEvent(QFocusEvent * event);

    SurfaceMesh::SurfaceMeshModel * m;
    QString filename;
	bool isDeleted;

	bool isDrawWireframe;
	bool isDoubleLight;
	QColor bg, fg;

	// Undo
	std::vector< SurfaceMesh::Vector3 > oldVertices;
	std::vector< std::vector<SurfaceMesh::Vertex> > oldFaces;

	// Original
	std::vector< SurfaceMesh::Vector3 > originalVertices;
	std::vector< std::vector<SurfaceMesh::Vertex> > originalFaces;

	// Debug:
	QVector<RenderObject::Base*> debugItems;

signals:
	void gotFocus(MyDrawArea *);
};

extern MyDrawArea * lastSelected;

inline void cleanUp(SurfaceMesh::SurfaceMeshModel * m){
	SurfaceMesh::Vector3VertexProperty points = m->vertex_coordinates();

	// Center to zero point
	SurfaceMesh::Vector3 mean(0,0,0);
	foreach(SurfaceMesh::Vertex v, m->vertices()) mean += points[v];
	mean /= m->n_vertices();
	foreach(SurfaceMesh::Vertex v, m->vertices()) points[v] -= mean;

	// Scale maximum dimension to 1.0
	m->updateBoundingBox();
	Eigen::AlignedBox3d orig_bbox = m->bbox();
	SurfaceMesh::Vector3 d = orig_bbox.sizes();
	double s = (d.x() > d.y())? d.x():d.y();
	s = (s>d.z())? s : d.z();
	foreach(SurfaceMesh::Vertex v, m->vertices()) points[v] /= s;

	// Move above floor
	double minZ = DBL_MAX;
	foreach(SurfaceMesh::Vertex v, m->vertices()) minZ = qMin(minZ, points[v].z());
	foreach(SurfaceMesh::Vertex v, m->vertices()) points[v] -= SurfaceMesh::Vector3(0,0,minZ);

	m->updateBoundingBox();
	m->update_face_normals();
	m->update_vertex_normals();
}

#endif // MYDRAWAREA_H
