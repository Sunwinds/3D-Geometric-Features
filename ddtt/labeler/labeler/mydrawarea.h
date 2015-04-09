#ifndef MYDRAWAREA_H
#define MYDRAWAREA_H

#include <qglviewer/qglviewer.h>

#include "StructureGraph.h"

#include "Octree.h"

extern QString curLabel;
extern QStringList AllLabels;
extern QVector<QColor> UniqueColors;
extern QVector<QString> labelNames;

class MyDrawArea : public QGLViewer
{
public:
	MyDrawArea(Structure::Graph * graph) : g(graph){}
	~MyDrawArea();

	Structure::Graph * g;

	void draw();

	void postSelection(const QPoint& point);
	void mousePressEvent(QMouseEvent * event);
};

#endif // MYDRAWAREA_H
