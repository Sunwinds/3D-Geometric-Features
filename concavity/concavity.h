#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class concavity: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "sdf.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Compute harmonic fields"; }
    QString description() { return "Compute a concavity-aware segmentation field"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};
