#include "experiment.h"

#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#include "GraphCorresponder.h"
#include "StructureGraph.h"
#include "Scheduler.h"
#include "TopoBlender.h"

void experiment::initParameters(RichParameterSet *pars)
{

	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

void experiment::applyFilter(RichParameterSet *pars)
{
	QVector<Structure::Graph *> graphs;

	graphs << new Structure::Graph("C:/Temp/dataset/ChairBasic1/SimpleChair1.xml");
	graphs << new Structure::Graph("C:/Temp/dataset/kchair2bs/kchair2bs.xml");

	GraphCorresponder gcorr( graphs.front(), graphs.back() );
	


	QSharedPointer<Scheduler> scheduler ( new Scheduler );
	QSharedPointer<TopoBlender> blender = QSharedPointer<TopoBlender>( new TopoBlender(&gcorr, scheduler.data()) );

	scheduler->executeAll();
}
