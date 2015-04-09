#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class symmetry: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "symmetry.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Symmetry analysis"; }
    QString description() { return "Symmetry analysis"; }
	QKeySequence shortcut(){ return QKeySequence(Qt::CTRL + Qt::Key_A); }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};
