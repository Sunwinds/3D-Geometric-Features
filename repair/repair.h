#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class repair: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "repair.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Repair mesh"; }
    QString description() { return "Repair mesh by performing mesh fixing operations"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};
