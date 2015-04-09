#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

#include <Eigen/Sparse>

class bdf: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "bdf.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Compute BDF"; }
    QString description() { return "Surface BDF computation"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

private:
	Eigen::SparseMatrix<double> Lc();
	bool has_halfedge(Vertex start, Vertex end);
};
