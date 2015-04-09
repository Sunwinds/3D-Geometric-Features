#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

#include "StructureGraph.h"

class analyze: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "analyze.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Analyse structure"; }
    QString description() { return "Analyse structure"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

	// Run without mesh
	bool isApplicable(Starlab::Model*) { return true; }
};
