#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

#include "StructureGraph.h"

class functional: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "functional.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Compute functional metric"; }
    QString description() { return "Find the functional properties"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

	// Run without mesh
	bool isApplicable(Starlab::Model*) { return true; }

	// Measures:
	void measure_direction(Structure::Graph * graph);
	void measure_ground(Structure::Graph * graph);
	void measure_ground_parts(Structure::Graph * graph);
	void measure_height(Structure::Graph * graph);
	void measure_curvyness(Structure::Graph * graph);
	void measure_area_volume_ratio(Structure::Graph * graph);
	void measure_access(Structure::Graph * graph);
	void measure_graph_agd(Structure::Graph * graph);
	void measure_physical_system(Structure::Graph * graph);
	void measure_position(Structure::Graph * graph);

	// Other:
	void simplified_graph();
};
