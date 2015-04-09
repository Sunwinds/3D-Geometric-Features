#pragma once
#include "GraphCorresponder.h"
#include "Task.h"
#include "Scheduler.h"
#include "TopoBlender.h"
#include "SynthesisManager.h"
class ProjectedStructureGraph;

typedef QMap< QString,QVariant > PropertyMap;
typedef QPair< QVector<QString>, QVector<QString> > Pairing;
typedef QVector< Pairing > VectorPairings;
typedef QVector< VectorPairings > Assignments;
typedef QVector< QPair<QString,QString> > VectorPairStrings;

class DeformationPath
{
public:
    DeformationPath();

    double weight;
    VectorPairStrings pairsDebug;
    QVector<Pairing> pairs;
    QVector<double> errors;

    QMap<QString, QColor> scolors, tcolors;

    int idx, i, si;
    PropertyMap property;

    GraphCorresponder * gcorr;
    QSharedPointer<Scheduler> scheduler;
    QSharedPointer<TopoBlender> blender;
	QSharedPointer<SynthesisManager> synthman;

	QVector<Structure::Graph*> samples;

	QVector<ProjectedStructureGraph*> projected;

	void execute();

	// Visualization
	void renderSamples();
	void renderProxies();
	void renderProjected();

	// Experimental
	void badMorphing();
};

static inline bool DeformationPathCompare (const DeformationPath & i, const DeformationPath & j) { return (i.weight < j.weight); }

