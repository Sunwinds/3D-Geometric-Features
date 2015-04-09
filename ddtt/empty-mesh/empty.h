#pragma once
#include "SurfaceMeshPlugins.h"
#include "SurfaceMeshHelper.h"

class empty : public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "empty.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Create empty mesh"; }
    QString description() { return "Creates an empty surface mesh."; }
    bool isApplicable(Starlab::Model*) { return true; }
    void applyFilter(RichParameterSet*){
        SurfaceMeshModel * m = new SurfaceMeshModel("empty.obj", "empty");
        m->add_vertex(Vector3(0,0,0));
        document()->addModel(m);
    }
};
