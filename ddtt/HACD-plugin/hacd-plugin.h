#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class hacd: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "hacd.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "HACD"; }
    QString description() { return "Hierarchical Approximate Convex Decomposition by Khaled Mamou"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};
