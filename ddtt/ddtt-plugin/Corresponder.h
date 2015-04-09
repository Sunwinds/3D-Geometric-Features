#pragma once
#include <QVector>
#include "RenderObject.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Structure{ struct Graph; }

typedef QMap<QString,QVariant> PropertyMap;

class Corresponder{
public:
	enum Metric{ ALL, POSITIONAL, HEIGHT, ROTATIONAL, SCALING, STRUCTURE, DISTORTION };

private:
	Structure::Graph * source;
	Structure::Graph * target;

	Eigen::MatrixXd metric;

	Eigen::MatrixXd geometricDistance();
	Eigen::MatrixXd structuralDistance();

	int n, m;
	int selected_nidx;
	Metric selected_metric;
	int num_splits;

public:
	Corresponder(Structure::Graph * g1, Structure::Graph * g2, int nidx = 0, Corresponder::Metric metric = ALL, int splitTo = 1);

	void compute( PropertyMap prop = PropertyMap() );

	QVector<RenderObject::Base *> debug;
};
