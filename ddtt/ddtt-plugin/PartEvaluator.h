#pragma once
#include "StructureGraph.h"
#include "GraphCorresponder.h"
#include "TopoBlender.h"
#include "Scheduler.h"
#include "SchedulerWidget.h"
#include "Task.h"
#include "SynthesisManager.h"
#include "SoftwareRenderer.h"
#include "ImageCompare.h"
#include "MarchingSquares.h"
#include "ContourMappingDistance.h"

struct PartEvaluator{
	static Eigen::MatrixXd renderGraphBinary( Structure::Graph * g, SynthesisManager * synthman )
	{
		int width = 128, height = 128;

		// Setup camera
		Eigen::AlignedBox3d graphBBox = g->bbox();
		double distance = graphBBox.sizes().maxCoeff() * 4.5;
		Vector3 direction (-2,-2,0.8);
		direction.normalize();
		Vector3 target = graphBBox.center();
		Vector3 eye = (direction * distance) + target;
		Vector3 up(0,0,1);
		Eigen::MatrixXd camera = SoftwareRenderer::CreateViewMatrix(eye, target, up);

		// Draw geometry into a buffer
		std::vector<SimplePolygon> polys = synthman->drawWithProxies(g);
		QVector< QVector<Vector3> > tris;

		for(auto p : polys){
			if(p.vertices.size() < 4)
				tris.push_back( QVector<Vector3>() << p.vertices[0] << p.vertices[1] << p.vertices[2] );
			else{
				tris.push_back( QVector<Vector3>() << p.vertices[0] << p.vertices[1] << p.vertices[2] );
				tris.push_back( QVector<Vector3>() << p.vertices[2] << p.vertices[3] << p.vertices[0] );
			}
		}

		if( !tris.size() ) return Eigen::MatrixXd::Zero(0,0);

		return SoftwareRenderer::render(tris, width, height, camera);
	}

	static void evaluate( Structure::Graph * g )
	{
		ImageCompare im;

		// Add source and target to knowledge
		{
			Structure::Graph *source = new Structure::Graph( *g ), *target = new Structure::Graph( *g );
			QString newName = target->property["name"].toString().replace(target->name(), "Target");
			target->property["name"].setValue( newName );
			DeformationPath path;
			path.gcorr = new GraphCorresponder(source, target);
			path.scheduler = QSharedPointer<Scheduler>( new Scheduler );
			path.blender = QSharedPointer<TopoBlender>( new TopoBlender( path.gcorr, path.scheduler.data() ) );
			path.synthman = QSharedPointer<SynthesisManager>( new SynthesisManager(path.gcorr, path.scheduler.data(), path.blender.data()) );
			path.synthman->makeProxies(60, 20);
			path.scheduler->executeAll();
			QVector<Eigen::MatrixXd> buffers;
			buffers << renderGraphBinary( path.scheduler->allGraphs.front(), path.synthman.data() );
			buffers << renderGraphBinary( path.scheduler->allGraphs.back(), path.synthman.data() );
			for(auto buffer : buffers){
				std::vector< std::pair<double,double> > contour;
				for(auto p : MarchingSquares::march(buffer, 1.0)) contour.push_back( std::make_pair(p.x(), p.y()) );
				ImageCompare::Instance inst( contour );
				inst.cachedImage = SoftwareRenderer::matrixToImage(buffer);
				im.addInstance("chairs", inst);
			}
		}

		double maxError = -DBL_MAX;
		double minError = DBL_MAX;

		for(auto n : g->nodes)
		{
			Structure::Graph * source = new Structure::Graph( *g );
			Structure::Graph * target = new Structure::Graph( *g );

			// Rename target
			QString newName = target->property["name"].toString().replace(target->name(), "Target");
			target->property["name"].setValue( newName );

			GraphCorresponder * gcorr = new GraphCorresponder( source, target );
			gcorr->setNonCorresSource( n->id );
			gcorr->computeCorrespondences();

			Scheduler* scheduler = new Scheduler;
			TopoBlender * blender = new TopoBlender( gcorr, scheduler );
			SynthesisManager * synthman = new SynthesisManager(gcorr, scheduler, blender);

			scheduler->startAllSameTime( Task::DEFAULT_LENGTH );
			int idx = scheduler->tasks.indexOf( scheduler->getTaskFromNodeID(n->id) );
			scheduler->tasks[idx]->setStart( 0 );
			int endTime = scheduler->tasks[idx]->endTime();

			if( false )
			{
				SchedulerWidget * widget = new SchedulerWidget( scheduler );
				widget->setAttribute(Qt::WA_DeleteOnClose);
				widget->show();
				break;
			}

			// Deform
			synthman->makeProxies(60, 20);
			scheduler->executeAll();

			// Compare graph after node removed
			int midx = (double(endTime) / scheduler->totalExecutionTime()) * (scheduler->allGraphs.size()-1);
			Structure::Graph * modified = scheduler->allGraphs[midx];

			Eigen::MatrixXd buffer = renderGraphBinary( modified, synthman );

			// Find outer most contour using marching squares
			std::vector< std::pair<double,double> > contour;
			for(auto p : MarchingSquares::march(buffer, 1.0)) contour.push_back( std::make_pair(p.x(), p.y()) );

			// Compute signature
			ImageCompare::Instance sampleInstance( contour );

			// Compare with respect to knowledge
			ImageCompare::InstanceMatches neighbours = im.kNearest( sampleInstance );

			std::vector<Eigen::Vector2d> cA, cB;
			for(auto p : sampleInstance.contour) cA.push_back(Eigen::Vector2d(p.first, p.second));
			for(auto p : neighbours.front().second.contour) cB.push_back(Eigen::Vector2d(p.first, p.second));

			double cmm = optMapMae(cA,cB);

			n->property["cmm"].setValue( cmm );

			maxError = std::max(maxError, cmm);
			minError = std::min(minError, cmm);
		}

		g->setColorAll( QColor(0,0,255) );

		for(auto n : g->nodes)
		{
			double cmm = n->property["cmm"].toDouble();
			cmm = (cmm - minError) / (maxError - minError);
			n->property["cmm"].setValue( cmm );

			g->setColorFor( n->id, starlab::qtJetColor(cmm) );
			n->vis_property["meshSolid"].setValue( true );
		}
	}
};
