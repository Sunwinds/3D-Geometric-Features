#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class curvature: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "curvature.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Compute curvature"; }
    QString description() { return "Surface curvature computation"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

	void collectEnoughRings(Vertex v, const size_t min_nb, std::vector<int> & all);
};
