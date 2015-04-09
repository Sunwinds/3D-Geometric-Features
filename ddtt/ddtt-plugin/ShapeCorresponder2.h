#pragma once

#include <QThread>
#include <QVector>
#include <QVariant>
#include <QProgressDialog>
#include "RenderObjectExt.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Structure{ struct Graph; }
class BofSearchManager;

#include "DeformationPath.h"

class ShapeCorresponder2 : public QThread{
    Q_OBJECT
public:
    Structure::Graph * source;
    Structure::Graph * target;
    PropertyMap property;

    std::vector<DeformationPath> paths;
    BofSearchManager * bmanager;

    ShapeCorresponder2(Structure::Graph * g1, Structure::Graph * g2);

    QVector<RenderObject::Base *> debug;
    QProgressDialog * pd;

public slots:
    void run();
    void increaseProgress();

signals:
    void done();
    void pathComputed();
};

