#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

#include <iterator>
#include <math.h>
#include <vector>

#include "Eigen/LU"
#include "Eigen/Geometry"
#include "Eigen/Eigenvalues"

typedef double											FT;
typedef Eigen::Vector3d									Point_3;
typedef Point_3											Vector_3;
typedef std::vector<Point_3>::iterator					InputIterator;
typedef Eigen::Transform<FT, 3, Eigen::AffineCompact>	Aff_transformation;
typedef Eigen::VectorXd									LAVector;
typedef Eigen::MatrixXd									LAMatrix;
typedef Surface_mesh::Vertex_iterator                   VertexIterator;
typedef Surface_mesh::Vertex_property<Surface_mesh::Point> Points;

class PCAShapeEst: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "PCAShapeEst.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Compute PCAScalar"; }
    QString description() { return "Shape Estimation computation"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

    Vector_3 computePCAScalar(VertexIterator begin, VertexIterator end, Points pts);
};
