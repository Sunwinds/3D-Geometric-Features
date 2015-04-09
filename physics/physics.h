#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

// Compute physical quantities

class physics: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "physics.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Compute physical quantities"; }
    QString description() { return "Measures different physical quantities."; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};
