#pragma once

#include "SurfaceMeshPlugins.h"
#include "SurfaceMeshHelper.h"
#include "RichParameterSet.h"

namespace Structure{ struct Graph; }

class ddtt : public SurfaceMeshModePlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "ddtt.plugin.starlab")
    Q_INTERFACES(ModePlugin)

    QIcon icon(){ return QIcon(":/images/icon.png"); }

public:
	ddtt() { widget = NULL; }

    /// Functions part of the EditPlugin system
    void create();
    void destroy(){}

    void decorate();
	void setSceneBounds();

	bool keyPressEvent( QKeyEvent* event );

    QWidget * widget;

	QVector<Structure::Graph*> graphs;
	Eigen::AlignedBox3d bigbox;

public slots:
	void execute();
	void loadModels(QStringList fileNames);
	void loadJob(QString filename);
	void clear();
	void loadGraphs();
	void loadKnowledge();
	void correspond();
	void postCorrespond();
	void experiment();
};
