#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

// Over segment using skeletonization

class segmentation: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "segmentation.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Segment using skeleton"; }
    QString description() { return "Over segment using skeletonization"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};

// Utility:
static inline QColor qRandomColor2(double saturation = 0.5, double val = 0.95){
	double golden_ratio_conjugate = 0.618033988749895;
	double h = ((double)rand() / RAND_MAX);
	h += golden_ratio_conjugate;
	h = fmod(h, 1.0);
	return QColor::fromHsvF(h, saturation, val);
}
