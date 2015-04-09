#pragma warning(disable:4267)

#pragma once

#include <QThread>
#include <QVector>
#include <QVariant>
#include "RenderObjectExt.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Structure{ struct Graph; }

#include "DeformationPath.h"

class ShapeCorresponder : public QThread{
	Q_OBJECT
public:
	Structure::Graph * source;
	Structure::Graph * target;
	PropertyMap property;

	std::vector<DeformationPath> paths;

	ShapeCorresponder(Structure::Graph * g1, Structure::Graph * g2, QString knowledge = QString());

	QVector<RenderObject::Base *> debug;

public slots:
	void run();
	void increaseProgress();

signals:
	void done();
	void pathComputed();
};
