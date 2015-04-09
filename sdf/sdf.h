#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

#include <Eigen/Sparse>

#include "Octree.h"

class sdf: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "sdf.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Compute SDF"; }
    QString description() { return "The Shape-Diameter Function"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

	static double compute_volume(const Vector3 & p, const Vector3 & n, 
		int numCones, double coneSeparation, int raysInCone, bool gaussianWeights, 
		Octree * octree, QList<Ray> * allRays = 0);
};
