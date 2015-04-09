#include "labeler.h"
#include "ui_labeler.h"
#include "StarlabDrawArea.h"
#include <QPushButton>
#include <QFileDialog>
#include "mydrawarea.h"
#include "StructureGraph.h"

Labeler::Labeler(QWidget *parent) : QWidget(parent), ui(new Ui::Labeler)
{
	// Anti-aliasing when using QGLWidget or subclasses
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSamples(8);
	QGLFormat::setDefaultFormat(glf);

    ui->setupUi(this);

	UniqueColors << QColor(255,97,121) <<  QColor(255,219,88) << QColor(107,255,135) << QColor(255,165,107) << QColor(104,126,255) <<
		QColor(242,5,135) << QColor(138,0,242) << QColor(3,166,60) << QColor(242,203,5) << QColor(255,0,0) << 
		QColor(0,255,0) << QColor(0,0,255) << QColor(0,255,255) << QColor(255,255,0) << QColor(200,200,200) << QColor(60,60,60);

	labelNames << "leg"<< "seat"<< "back"<< "armrest"<< "back_top"<< "leg_support" 
		<< "seat_support" << "back_support"<< "bottom" << "tire" << "handle";

	// Create label buttons
	{
		int numColors = UniqueColors.size();
		int numLabels = labelNames.size();

		QButtonGroup * group = new QButtonGroup(this);
		group->setExclusive(true);

		for(int i = 0; i < numColors; i++){
			QString lable = QString("part%1").arg(i+1);
			if(i < numLabels) lable = labelNames[i];

			QLayout * label_layout = ui->labelsFrame->layout();

			QPushButton * button = new QPushButton(lable);
			QString style = QString("background:%1").arg(UniqueColors[i].name());
			button->setStyleSheet( style );
			button->setAutoExclusive( true );
			button->setCheckable( true );

			AllLabels << lable;

			label_layout->addWidget( button );
			group->addButton( button );

			connect(group, SIGNAL(buttonClicked(QAbstractButton *)), this, SLOT(onGroupButtonClicked(QAbstractButton *)));
		}
	}

	{
		QString datasetPath = QFileDialog::getExistingDirectory(0, "Dataset");

		QDir datasetDir(datasetPath);
		QStringList subdirs = datasetDir.entryList(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

		foreach(QString subdir, subdirs)
		{
			// Special folders
			if(subdir == "corr") continue;
			QDir d(datasetPath + "/" + subdir);

			// Check if no graph is in this folder
			if( d.entryList(QStringList() << "*.xml", QDir::Files).isEmpty() ) continue;
			QString filename = d.absolutePath() + "/" + d.entryList(QStringList() << "*.xml", QDir::Files).join("");
			database << filename;
		}
	}

	// Create viewers
	{
		nU = 2;
		nV = 3;

		// Scroll bars
		int numPages = ceil(double(database.size()) / double(nU * nV)) - 1;
		ui->viewersScroll->setValue(0);
		ui->viewersScroll->setRange(0, numPages);
		connect(ui->viewersScroll, SIGNAL(valueChanged(int)), SLOT(loadMeshes()));

		// Load meshes
		loadMeshes();
	}

	curLabel = AllLabels.front();

}

void Labeler::onGroupButtonClicked(QAbstractButton * pbutton)
{
	QPushButton * button = (QPushButton *)pbutton;
	ui->curLabel->setText( "Current lable: " + button->text() );

	curLabel = button->text();
}

void Labeler::loadMeshes()
{
	QGridLayout * layout = (QGridLayout *)ui->viewersFrame->layout();
	QLayoutItem* item;
	while ( ( item = layout->takeAt( 0 ) ) != NULL ){
		delete item->widget(); delete item;
	}

	int offset = ui->viewersScroll->value() * (nU * nV);

	for(int i = 0; i < nU; i++){
		for(int j = 0; j < nV; j++){

			int idx = offset + (i*nV) + j;
			if(idx + 1 > database.size()) continue;
			QString filename = database[idx];

			Structure::Graph * g = new Structure::Graph( filename );
			g->property["showNodes"] = false;

			MyDrawArea * viewer = new MyDrawArea(g);
			layout->addWidget(viewer, i, j, 1, 1);
            //viewer->setGridIsDrawn();
			viewer->camera()->setSceneRadius(2);
			viewer->camera()->setUpVector(qglviewer::Vec(0,0,1));
			viewer->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
			viewer->camera()->lookAt(qglviewer::Vec());
			viewer->camera()->showEntireScene();

			Eigen::AlignedBox3d bbox = g->bbox();
			double s = bbox.sizes().norm() * 0.1;
			bbox.extend(bbox.min() + (Vector3(-1,-1,-1) * s));
			bbox.extend(bbox.max() + (Vector3(1,1,1) * s));

			viewer->camera()->fitBoundingBox(qglviewer::Vec(bbox.min().data()),qglviewer::Vec(bbox.max().data()));
			viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(bbox.center()));

			for(Structure::Node * node : g->nodes){
				node->vis_property["meshColor"].setValue( QColor(180, 180, 180) );
				node->vis_property["meshSolid"].setValue( true );

				if(node->meta.contains("label")){
					int idx = AllLabels.indexOf( node->meta["label"].toString() );
					node->vis_property["meshColor"].setValue( UniqueColors[idx] );
				}

				if(!node->property.contains("mesh")) continue;
				QSharedPointer<SurfaceMeshModel> nodeMesh = node->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
#ifndef QT_DEBUG
				Octree * octree = new Octree( nodeMesh.data() );
				node->property["octree"].setValue( octree );
#endif
			}

			g->saveToFile("C:/temp/graph.xml");
		}
	}

	this->activateWindow();
	this->setFocus();
	this->raise();
}

Labeler::~Labeler()
{
    delete ui;
}
