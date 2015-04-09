#pragma warning(disable:4267)

#include <QFileDialog>

#include "ddtt-plugin.h"
#include "interfaces/ModePluginDockWidget.h"
#include "ui_ddtt_widget.h"
#include "ddtt_widget.h"

#include "StructureGraph.h"
#include "SynthesisManager.h"
#include "GraphExplorer.h"

#include "Corresponder.h"

#include "ShapeCorresponder.h"
#include "ImageCompare.h"

#include "ShapeCorresponder2.h"

#include "DeformScene.h"

#define BBOX_WIDTH(box) (box.max().x()-box.min().x())
#define PADDING_FACTOR 1.0

ddtt_widget * w = NULL;
//ShapeCorresponder * sc = NULL;

ShapeCorresponder2 * sc2 = NULL;

void ddtt::create()
{
	if(!widget)
	{
		ModePluginDockWidget * dockwidget = new ModePluginDockWidget("TopoBlender", mainWindow());

		// Setup widget and signals
        w = new ddtt_widget();

		connect(w->ui->testButton, SIGNAL(clicked()), SLOT(execute()));
		connect(w->ui->clearButton, SIGNAL(clicked()), SLOT(clear()));
		connect(w->ui->experimentButton, SIGNAL(clicked()), SLOT(experiment()));

		connect(w->ui->loadGraphs, SIGNAL(clicked()), SLOT(loadGraphs()));
		connect(w->ui->correspondButton, SIGNAL(clicked()), SLOT(correspond()));

		dockwidget->setWidget(w);
		dockwidget->setWindowTitle(w->windowTitle());
		mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);
		this->widget = w;

		drawArea()->setShortcut(QGLViewer::DRAW_AXIS, Qt::Key_A);
		drawArea()->setShortcut(QGLViewer::DRAW_GRID, Qt::Key_G);

		//w->ui->loadGraphs->click();
		//w->ui->correspondButton->click();

		//loadKnowledge();
	}
}

#include "PartEvaluator.h"
void ddtt::experiment()
{
	PartEvaluator::evaluate( graphs.front() );
	if(graphs.size()>1) PartEvaluator::evaluate( graphs.back() );
}

void ddtt::loadKnowledge()
{
	QElapsedTimer loadTimer; loadTimer.start();

	ImageCompare im;
	
	//im.loadKnowledge("C:/Temp/_imageSearch/test_data", "chairs");
	im.loadKnowledge("C:/Temp/_imageSearch/all_images_apcluster_data", "chairs");

	mainWindow()->setStatusBarMessage(QString("Loaded knowledge dataset [%1] with %2 shapes (%3 ms)").arg(
		"chairs").arg(im.datasetSize("chairs")).arg(loadTimer.elapsed()));

	// Find duplicates:
	if( false )
	{	
		ImageCompare::Instance inst = im.getInstance("chairs");

		// Timing:
		QElapsedTimer queryTimer; queryTimer.start();
		ImageCompare::InstanceMatches results = im.kNearest( inst, 7 );
		mainWindow()->setStatusBarMessage( QString("Query time (%1 ms)").arg(queryTimer.elapsed()) );

		// Debug:
		ImageCompare::showInstances( ImageCompare::InstanceMatches() << qMakePair(0,inst) );
		ImageCompare::showInstances( im.kNearest( inst, 7 ) );

		// Clean up:
		//im.removeDuplicateSets(0.2, "chairs");

		for(auto dupset : im.duplicateSets( 0.3 ))
		{
			QVector<QImage> setimgs;
			for(auto instance : dupset)	setimgs << instance.image();
			showImages( setimgs );
		}
	}
}

void ddtt::decorate()
{
	double startX = bigbox.min().x();

	for(int g = 0; g < (int) graphs.size(); g++)
	{
		// Place and draw graph
		glPushMatrix();

		Eigen::AlignedBox3d curbox = graphs[g]->bbox();

		double curwidth = (curbox.max().x() - curbox.min().x());
		double deltaX = curwidth * 0.5;

		double padding = 0;
		if(g > 0) padding = curwidth * PADDING_FACTOR;

		double posX = startX + deltaX + padding;

		if(graphs.size() < 2) posX = 0;

		glTranslated(posX, 0, 0);

		// store for later use
		graphs[g]->property["posX"] = posX;
		graphs[g]->draw();
		//drawBBox( curbox );

		glPopMatrix();

		startX += curwidth + padding;
	}
}

void ddtt::setSceneBounds()
{
	if(!graphs.size())
	{
		drawArea()->setSceneRadius(2);
		drawArea()->setSceneCenter(qglviewer::Vec(0,0,0));
		drawArea()->setSceneBoundingBox(qglviewer::Vec(-1,-1,-1), qglviewer::Vec(1,1,1));
		drawArea()->camera()->setPosition(qglviewer::Vec(-1,-3,2));
		drawArea()->showEntireScene();
		drawArea()->updateGL();
		return;
	}

	// Set scene bounds
	bigbox = graphs.front()->bbox();
	double deltaX = BBOX_WIDTH(bigbox);
	bigbox.translate( Vector3(deltaX * 0.5, 0, 0) ); // start from zero

	for(int i = 1; i < (int)graphs.size(); i++)
	{
		Eigen::AlignedBox3d curbox = graphs[i]->bbox();

		double curWidth = BBOX_WIDTH(curbox);
		double padding = curWidth * PADDING_FACTOR;

		curbox.translate( Vector3(deltaX + (0.5 * curWidth) + padding, 0, 0) );
		bigbox = bigbox.merged( Eigen::AlignedBox3d(curbox) );

		deltaX += BBOX_WIDTH(curbox) + padding; 
	}

	// Move to center
	bigbox.translate( Vector3(-bigbox.center().x(), 0, 0) );

	// Setup camera
	{
		Vector3 a = bigbox.min();
		Vector3 b = bigbox.max();

		qglviewer::Vec vecA(a.x(), a.y(), a.z());
		qglviewer::Vec vecB(b.x(), b.y(), b.z());

		drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
		drawArea()->camera()->setPosition(qglviewer::Vec(-2,-3,1));
		drawArea()->camera()->lookAt(qglviewer::Vec(0,0,0));

		drawArea()->setSceneCenter((vecA + vecB) * 0.5);
		drawArea()->setSceneBoundingBox(vecA, vecB);
		drawArea()->showEntireScene();
		drawArea()->updateGL();
	}
}

void ddtt::loadModels(QStringList fileNames)
{
	if(fileNames.isEmpty()) return;

	foreach(QString file, fileNames){
		graphs.push_back( new Structure::Graph(file) );
	}

	QFileInfo fileInfo(fileNames.back());
	mainWindow()->settings()->set( "lastUsedDirectory", fileInfo.absolutePath() );

	setSceneBounds();
}

void ddtt::loadJob(QString filename)
{
	QFile job_file( filename );

	if (!job_file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
	QFileInfo jobFileInfo(job_file.fileName());
	QTextStream in(&job_file);

	// Job file content
	QString sgFileName, tgFileName, correspondenceFileName, scheduleFileName;
	int samplesCount;
	double gdResolution, timeStep;
	int reconLevel, renderCount;

	QString curPath = jobFileInfo.absolutePath() + "/";

	sgFileName = in.readLine();
	tgFileName = in.readLine();
	correspondenceFileName = in.readLine();
	scheduleFileName = in.readLine();

	in >> samplesCount;
	in >> gdResolution >> timeStep;
	in >> reconLevel >> renderCount;

	/*
	// Load graphs
	tb->graphs.push_back( new Structure::Graph ( curPath + sgFileName ) );
	tb->graphs.push_back( new Structure::Graph ( curPath + tgFileName ) );
	tb->setSceneBounds();
	tb->updateDrawArea();

	// Load correspondence
	tb->gcoor = tb->c_manager->makeCorresponder();
	tb->gcoor->loadCorrespondences( curPath + correspondenceFileName );
	tb->gcoor->isReady = true;
	*/

	graphs.push_back( new Structure::Graph ( curPath + sgFileName ) );
	graphs.push_back( new Structure::Graph ( curPath + tgFileName ) );

	graphs.front()->property["correspondenceFile"].setValue( curPath + correspondenceFileName );
	graphs.back()->property["correspondenceFile"].setValue( curPath + correspondenceFileName );

	mainWindow()->settings()->set( "lastUsedDirectory", jobFileInfo.absolutePath() );

	setSceneBounds();
}

void ddtt::loadGraphs()
{
	if( graphs.size() < 2 )
	{
		QStringList files = QFileDialog::getOpenFileNames(0, "Open Model", 
			mainWindow()->settings()->getString("lastUsedDirectory"), "All Supported (*.job *.xml);;Model Files (*.xml);;Job files (*.job)");
		if(files.size() < 1) return;

		if(files.front().contains(".job"))
		{
			loadJob( files.front() );
		}
		else
			loadModels( files );
	}
}

void ddtt::clear()
{
	drawArea()->clear();
	graphs.clear();
	drawArea()->updateGL();
}

void ddtt::execute()
{
	loadGraphs();

	int nidx = w->ui->nodeIndex->value();
	int metric = w->ui->metricBox->currentIndex();
	int splits = w->ui->splits->value();

	PropertyMap prop;
	prop["forceBest"] = w->ui->bestAssign->isChecked();

	Corresponder c(graphs.front(), graphs.back(), nidx, (Corresponder::Metric)metric, splits);
	c.compute( prop );

	drawArea()->clear();

	for(auto r : c.debug) drawArea()->addRenderObject(r);

	drawArea()->updateGL();
}

void ddtt::correspond()
{
	loadGraphs();
	if(graphs.empty()) return;

	// Make first has smaller number of nodes
	if(graphs.front()->nodes.size() > graphs.back()->nodes.size())
		std::swap( graphs.front(), graphs.back() );

    //sc = new ShapeCorresponder( graphs.front(), graphs.back(), "chairs" );
    sc2 = new ShapeCorresponder2(graphs.front(), graphs.back());

    //connect(sc, SIGNAL(done()), SLOT(postCorrespond()));
    connect(sc2, SIGNAL(done()), SLOT(postCorrespond()));

    //sc->start();
    sc2->start();
}

void ddtt::postCorrespond()
{
	drawArea()->clear();
    //for(auto r : sc->debug) drawArea()->addRenderObject(r);
    for(auto r : sc2->debug) drawArea()->addRenderObject(r);
	drawArea()->update();

	// Stats
	QStringList message;
    message << QString("Number of paths ( %1 ).").arg( sc2->property["pathsCount"].toInt() );
    message << QString("Prepare time ( %1 ms ).").arg( sc2->property["prepareTime"].toInt() );
    message << QString("Compute time ( %1 ms ).").arg( sc2->property["computeTime"].toInt() );
    message << QString("Evaluate time ( %1 ms ).").arg( sc2->property["evaluateTime"].toInt() );
	mainWindow()->setStatusBarMessage( message.join("\n") );

	// View correspondences
	DeformScene * ds = new DeformScene;

	int limit = 100;

	int i = 0;
	
    for(auto & path : sc2->paths)
	{
		path.idx = i++;
		ds->addDeformationPath( &path );

		if(limit > 0 && i > limit) break; 
	}

	ds->pack();
}

bool ddtt::keyPressEvent(QKeyEvent* event)
{
	if(event->key() == Qt::Key_I)
	{
		GraphExplorer * ge = new GraphExplorer;
		ge->update( graphs.back() );
		ge->show();

		return true;
	}

	return false;
}
