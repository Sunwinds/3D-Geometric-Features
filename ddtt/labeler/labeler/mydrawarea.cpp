#include "mydrawarea.h"
#include <QMouseEvent>

QString curLabel = "";
QStringList AllLabels;
QVector<QColor> UniqueColors;
QVector<QString> labelNames; 

void MyDrawArea::draw()
{
	g->draw();
}

void MyDrawArea::mousePressEvent(QMouseEvent * event)
{
	QGLViewer::mousePressEvent(event);

	if(event->modifiers() & Qt::SHIFT)
		select(event->pos());

	// Clear all labels
	if(event->buttons() & Qt::RightButton)
	{
		for(Structure::Node * node : g->nodes){
			node->vis_property["meshColor"].setValue( QColor(180, 180, 180) );
			node->meta.remove("label");
		}
	}
}

void MyDrawArea::postSelection(const QPoint& point)
{
	QMap<Structure::Node*,double> dists;

	for(Structure::Node* n : g->nodes)
	{
		if(!n->property.contains("mesh")) continue;
		QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();

		Octree * octree = n->property["octree"].value<Octree*>();
		if(!octree){
			Octree * octree = new Octree( nodeMesh.data() );
			n->property["octree"].setValue( octree );
		}
		
		qglviewer::Vec orig, dir;
		this->camera()->convertClickToLine(point, orig, dir);

		Vector3 o(orig[0],orig[1],orig[2]);
		Vector3 d(dir[0],dir[1],dir[2]);

		int faceIndex = -1;
		Vector3 isect = octree->closestIntersectionPoint( Ray(o,d), &faceIndex );
		if(faceIndex < 0) continue;

		dists[n] = (isect - o).norm();
	}

	if(!dists.size()) return;

	Structure::Node* n = sortQMapByValue(dists).first().second;

	int idx = AllLabels.indexOf(curLabel);
	if(idx < 0) return;

	n->vis_property["meshColor"].setValue( UniqueColors[idx] );
	n->meta["label"] = AllLabels[idx];
}

MyDrawArea::~MyDrawArea()
{
	QString filename = g->property["name"].toString();

	g->saveToFile(filename);

	delete g;
}
