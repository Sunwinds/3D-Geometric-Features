#include <QApplication>
#include "DeformationPath.h"
#include "SynthesisManager.h"
#include "ProjectedStructureGraph.h"

SynthesisManager * smanager = NULL;
Q_DECLARE_METATYPE(SynthesisManager *)

DeformationPath::DeformationPath(){
    gcorr = NULL;
    weight = 0.0;
    idx = i = si = 0;
}

void DeformationPath::execute()
{
	scheduler = QSharedPointer<Scheduler>( new Scheduler );
	blender = QSharedPointer<TopoBlender>( new TopoBlender( gcorr, scheduler.data() ) );
	
	scheduler->executeAll();
	for(auto & g : scheduler->allGraphs) g->moveBottomCenterToOrigin( true );
	property["isReady"].setValue( true );
}

void DeformationPath::renderSamples()
{
	if(!gcorr) return;

	qApp->setOverrideCursor(Qt::WaitCursor);

	scheduler = QSharedPointer<Scheduler>( new Scheduler );
	blender = QSharedPointer<TopoBlender>( new TopoBlender( gcorr, scheduler.data() ) );

	// Surface sampling
	smanager = new SynthesisManager(gcorr, scheduler.data(), blender.data(), 4000);
	smanager->genSynData();

	// Execute path
	scheduler->executeAll();
	for(auto & g : scheduler->allGraphs) g->moveBottomCenterToOrigin( true );

	if( scheduler->allGraphs.size() )
	{
		int N = 6;

		for(int i = 0; i < N; i++)
		{
			double t = double(i) / (N-1);
			QString path = QFileInfo(gcorr->sg->property["name"].toString()).absolutePath() + "/";
			QString name = "output";
			QString filename = path + QString("%1-%2-%3").arg(name).arg(i).arg(t);

			//int idx = t * (scheduler->allGraphs.size()-1);
			//smanager->renderGraph(*scheduler->allGraphs[idx], filename, false, 6);
		}
	}

	property["synthManager"].setValue( smanager );

	qApp->restoreOverrideCursor();
	QCursor::setPos(QCursor::pos());

	property["isReady"].setValue( true );
}

void DeformationPath::renderProxies()
{
	if(!gcorr) return;

	qApp->setOverrideCursor(Qt::WaitCursor);

	scheduler = QSharedPointer<Scheduler>( new Scheduler );
	blender = QSharedPointer<TopoBlender>( new TopoBlender( gcorr, scheduler.data() ) );

	// Surface sampling via proxies
	smanager = new SynthesisManager(gcorr, scheduler.data(), blender.data());
	smanager->makeProxies(60,20);

	// Execute path
	scheduler->executeAll();

	//this->badMorphing();

	for(auto & g : scheduler->allGraphs) g->moveBottomCenterToOrigin( true );

	property["synthManager"].setValue( smanager );

	qApp->restoreOverrideCursor();
	QCursor::setPos(QCursor::pos());

	property["isReady"].setValue( true );
}

void DeformationPath::renderProjected()
{	
	if(!gcorr) return;
	qApp->setOverrideCursor(Qt::WaitCursor);

	projected.push_back( new ProjectedStructureGraph(gcorr->sg, 256, true) );
	projected.push_back( new ProjectedStructureGraph(gcorr->tg, 256, true) );

	qApp->restoreOverrideCursor();
	QCursor::setPos(QCursor::pos());

	property["isReady"].setValue( true );
}

void DeformationPath::badMorphing()
{
	for(double t = 0; t <= 1.0; t += scheduler->timeStep)
	{
		// Current nodes geometry
		QMap<Structure::Node*, Array1D_Vector3> startGeometry;
		for(auto n : scheduler->activeGraph->nodes)
		{
			if(n->property["taskTypeReal"].toInt() == Task::GROW)
			{
				Array1D_Vector3 pnts = n->controlPoints();
				Vector3 p = pnts.front();
				Array1D_Vector3 newpnts(pnts.size(), p);
				for(auto & p: newpnts) p += Eigen::Vector3d(starlab::uniformRand(), starlab::uniformRand(), starlab::uniformRand()) * 1e-5;
				n->setControlPoints(newpnts);
			}

			if(n->property["taskTypeReal"].toInt() == Task::SHRINK)
			{
				auto tn = scheduler->targetGraph->getNode( n->property["correspond"].toString() );
				Array1D_Vector3 pnts = tn->controlPoints();
				Vector3 p = pnts[pnts.size() / 2];
				Array1D_Vector3 newpnts(pnts.size(), p);
				for(auto & p: newpnts) p += Eigen::Vector3d(starlab::uniformRand(), starlab::uniformRand(), starlab::uniformRand()) * 1e-5;
				tn->setControlPoints(newpnts);
			}

			startGeometry[n] = n->controlPoints();
		}

		// Morph nodes
		for(auto n : scheduler->activeGraph->nodes)
		{
			auto tn = scheduler->targetGraph->getNode( n->property["correspond"].toString() );
			if(!tn) continue;

			Array1D_Vector3 finalGeometry = tn->controlPoints();
			Array1D_Vector3 newGeometry;

			for(int i = 0; i < (int) finalGeometry.size(); i++)
				newGeometry.push_back( AlphaBlend(t, startGeometry[n][i], Vector3(finalGeometry[i])) );

			n->setControlPoints( newGeometry );
			n->property["t"].setValue( t );

			if(t > 0.5 && n->property["taskTypeReal"].toInt() == Task::SHRINK) 
				n->property["shrunk"].setValue( true );
		}

		scheduler->allGraphs.push_back(  new Structure::Graph( *scheduler->activeGraph )  );
	}

	property["progressDone"] = true;
}
