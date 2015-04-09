#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class experiment: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "experiment.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Deformation Correspondence"; }
    QString description() { return "Deformation Correspondence"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};
