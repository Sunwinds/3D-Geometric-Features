#include "meshbrowser.h"
#include "ui_meshbrowser.h"
#include "StarlabDrawArea.h"
#include <QPushButton>
#include <QFileDialog>

#include "mydrawarea.h"
MeshOperation curOp = NONE_OP;

using namespace SurfaceMesh;
MyDrawArea * lastSelected = NULL;

QVector<MyDrawArea*> activeViewers;

// Binary images rendering
#include "SoftwareRenderer.h"

MeshBrowser::MeshBrowser(QWidget *parent) : QWidget(parent), ui(new Ui::MeshBrowser)
{
	// Anti-aliasing when using QGLWidget or subclasses
	QGLFormat glf = QGLFormat::defaultFormat();
	glf.setSamples(8);
	QGLFormat::setDefaultFormat(glf);

    ui->setupUi(this);

    {
        QString datasetPath = QFileDialog::getExistingDirectory(0, "Dataset");

        // Get list of files
        QStringList filters;
        filters << "*.obj" << "*.off";
        QStringList files = QDir(datasetPath).entryList(filters);

        foreach(QString mesh, files)
        {
            database << (datasetPath + "/" + mesh);
		}
	}

	// Settings
	{
		// Use local file
		QString path = QDir(qApp->applicationDirPath()).absoluteFilePath("meshbrowser.ini");
		QSettings qsettings(path, QSettings::IniFormat);

		// Defaults  
		if(!qsettings.allKeys().contains("nV")){
			qsettings.setValue("nV", 3);
			qsettings.setValue("nU", 2);
			qsettings.sync();
		}
 
		nV = qsettings.value("nV", 3).toInt();
		nU = qsettings.value("nU", 2).toInt();
	}

	// Create viewers
	{
		// Scroll bars
		int numPages = ceil(double(database.size()) / double(nU * nV)) - 1;
		ui->viewersScroll->setValue(0);
		ui->viewersScroll->setRange(0, numPages);
		connect(ui->viewersScroll, SIGNAL(valueChanged(int)), SLOT(loadMeshes()));

		// Load meshes
		loadMeshes();
	}

    // Connections
	connect(ui->noneButton, &QPushButton::released, [=]() {curOp = MeshOperation::NONE_OP; ui->curLabel->setText(meshOps[curOp]);});
    connect(ui->rotateLeft, &QPushButton::released, [=]() {curOp = MeshOperation::ROTATE_LEFT; ui->curLabel->setText(meshOps[curOp]);});
    connect(ui->rotateRight, &QPushButton::released, [=]() {curOp = MeshOperation::ROTATE_RIGHT; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->rotateUp, &QPushButton::released, [=]() {curOp = MeshOperation::ROTATE_UP; ui->curLabel->setText(meshOps[curOp]);});
    connect(ui->flipNormals, &QPushButton::released, [=]() {curOp = MeshOperation::FLIP; ui->curLabel->setText(meshOps[curOp]);});
    connect(ui->invertPart, &QPushButton::released, [=]() {curOp = MeshOperation::INVERT_PART; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->removeSmallest, &QPushButton::released, [=]() {curOp = MeshOperation::REMOVE_SMALL; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->removeLargest, &QPushButton::released, [=]() {curOp = MeshOperation::REMOVE_LARGE; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->removeSelected, &QPushButton::released, [=]() {curOp = MeshOperation::REMOVE_SELECTED; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->closeHoles, &QPushButton::released, [=]() {curOp = MeshOperation::CLOSE_HOLES; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->mergeVertices, &QPushButton::released, [=]() {curOp = MeshOperation::MERGE_VERTICES; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->undoButton, &QPushButton::released, [=]() {curOp = MeshOperation::UNDO; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->undoAll, &QPushButton::released, [=]() {curOp = MeshOperation::UNDO_ALL; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->unifyOrientation, &QPushButton::released, [=]() {curOp = MeshOperation::UNIFY; ui->curLabel->setText(meshOps[curOp]);});
	connect(ui->componentsButton, &QPushButton::released, [=]() {curOp = MeshOperation::COMPONENTS; ui->curLabel->setText(meshOps[curOp]);});

	connect(ui->addDelete, &QPushButton::released, [=](){ 
		if(lastSelected && deletedItems().contains(lastSelected->filename)){
			for(auto i : ui->deleteList->findItems(lastSelected->filename, Qt::MatchExactly))
				delete ui->deleteList->takeItem( ui->deleteList->row(i) );
			refreshViewers();
			return;
		}
		if(lastSelected) ui->deleteList->insertItem(0, lastSelected->filename); 
	});

	connect(ui->deleteButton, &QPushButton::released, [=](){
		for(auto filename : deletedItems()){
			QFile::remove(filename);
			for(int i = 0; i < ui->deleteList->count(); ++i){
				for(auto viewer: activeViewers){
					if(viewer->filename == filename){
						viewer->isDeleted = true;
						break;
					}
				}
			}
		}
		ui->deleteList->clear();
	});

	connect(ui->generateThumbnails, &QPushButton::released, [=](){
		foreach(QString filename, database)
		{
			SurfaceMesh::SurfaceMeshModel m(filename, QFileInfo(filename).baseName());
			m.read( qPrintable(filename) );

			Vector3VertexProperty points = m.vertex_coordinates();

			/// Normalize, center, and move to base
			cleanUp(&m);

			MyDrawArea * viewer = new MyDrawArea(&m, filename);

			if( !QFileInfo(filename).exists() ){
				viewer->isDeleted = true;
				continue;
			}

			viewer->setMinimumSize(256,256);
			viewer->setMaximumSize(256,256);
			
			viewer->bg = QColor(255,255,255);
			viewer->fg = QColor(160,160,200);
			viewer->isDrawWireframe = false;
			viewer->isDoubleLight = true;

			viewer->show();

			{
				Eigen::AlignedBox3d bbox = m.bbox();
				viewer->camera()->setSceneRadius(m.bbox().sizes().norm() * 2);
				viewer->camera()->setUpVector(qglviewer::Vec(0,0,1));
				viewer->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
				viewer->camera()->lookAt(qglviewer::Vec());
				viewer->camera()->showEntireScene();
				double S = bbox.sizes().norm() * 0.2;
				bbox.extend(bbox.min() + (Vector3(-1,-1,-1) * S));
				bbox.extend(bbox.max() + (Vector3(1,1,1) * S));
				viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(bbox.center()));
				viewer->camera()->setSceneCenter(qglviewer::Vec(bbox.center()));
				viewer->camera()->fitBoundingBox(qglviewer::Vec(bbox.min().data()),qglviewer::Vec(bbox.max().data()));
			}

			QFileInfo meshFileInfo( filename );
			QString meshname = meshFileInfo.baseName();
			QString path = meshFileInfo.absolutePath();
			QString thumbImgFilename = QString("%1.png").arg(path + "/" + meshname);

			qApp->processEvents();

			viewer->raise();
			viewer->makeCurrent();
			viewer->grabFrameBuffer().save( thumbImgFilename );

			qApp->processEvents();

			viewer->hide();
			delete viewer;

			qApp->processEvents();
		}

	});

	connect(ui->genBinaryImgsButton, &QPushButton::released, [=](){
		if(database.empty()) return;

		#pragma omp parallel for
		for(int i = 0; i < (int)database.size(); i++)
		{
			QString meshfile = database.at(i);

			SurfaceMesh::SurfaceMeshModel mesh;
			mesh.read( meshfile.toStdString() );
			mesh.updateBoundingBox();

			QVector< QVector<Eigen::Vector3d> > faces;

			Vector3VertexProperty points = mesh.vertex_coordinates();
			for(auto f : mesh.faces()){
				QVector<Eigen::Vector3d> face;
				for(auto v : mesh.vertices(f)){
					face.push_back( points[v] );
				}
				faces.push_back(face);
			}

			QFileInfo meshFileInfo(meshfile);
			QString meshname = meshFileInfo.baseName();
			QString path = meshFileInfo.absolutePath();
			QString binaryImgFilename = QString("%1.png").arg(path + "/" + meshname);

			// Setup camera
			Eigen::AlignedBox3d bbox = mesh.bbox();
			double distance = bbox.sizes().maxCoeff() * 3.5;
			Vector3 direction (-1.25,-2,0.7);
			direction.normalize();
			Vector3 target = bbox.center();
			Vector3 eye = (direction * distance) + target;
			Vector3 up(0,0,1);
			Eigen::MatrixXd camera = SoftwareRenderer::CreateViewMatrix(eye, target, up);

			SoftwareRenderer::matrixToImage( SoftwareRenderer::render(faces, 128, 128, camera) ).save( binaryImgFilename );
		}
	});
}

QStringList MeshBrowser::deletedItems(){
	QStringList items;
	for(int i = 0; i < ui->deleteList->count(); ++i){
		QListWidgetItem* item = ui->deleteList->item(i);
		items << item->text();
	}
	return items;
}

void MeshBrowser::loadMeshes()
{
	lastSelected = NULL;

	QGridLayout * layout = (QGridLayout *)ui->viewersFrame->layout();
	QLayoutItem* item;
	while ( ( item = layout->takeAt( 0 ) ) != NULL ){
		delete item->widget(); delete item;
	}
	layout->setMargin(0);
	layout->setSpacing(0);

	int offset = ui->viewersScroll->value() * (nU * nV);
	int numPages = (database.size()-1) / (nU * nV);
	int pageID = (double(offset) / (database.size()-1)) * numPages;
	ui->curPage->setText(QString("Page: %1 / %2").arg(pageID).arg(numPages));
	
	activeViewers.clear();

	for(int i = 0; i < nU; i++){
		for(int j = 0; j < nV; j++){

			int idx = offset + (i*nV) + j;
			if(idx + 1 > database.size()) continue;
			QString filename = database[idx];

            //Structure::Graph * g = new Structure::Graph( filename );
            //g->property["showNodes"] = false;


            SurfaceMesh::SurfaceMeshModel * m = new SurfaceMeshModel(filename, QFileInfo(filename).baseName());
            m->read( qPrintable(filename) );

            Vector3VertexProperty points = m->vertex_coordinates();

            /// Normalize, center, and move to base
			cleanUp(m);

            MyDrawArea * viewer = new MyDrawArea(m, filename);

			viewer->setForegroundColor(QColor(255,255,255));

			if( !QFileInfo(filename).exists() ){
				viewer->isDeleted = true;
				continue;
			}

			// Save original
			Vector3VertexProperty origPnts = m->vertex_coordinates();
			for(auto v: m->vertices()) viewer->originalVertices.push_back(origPnts[v]);
			for(auto f: m->faces()){
				std::vector<Vertex> face;
				for(auto v: m->vertices(f)) face.push_back(v);
				viewer->originalFaces.push_back(face);
			}

			QFrame * frame = new QFrame;
			QHBoxLayout * ll = new QHBoxLayout;
			frame->setLayout(ll);
			ll->addWidget(viewer);
			layout->setMargin(2);
			layout->setSpacing(2);
			layout->addWidget(frame, i, j, 1, 1);

			activeViewers.push_back(viewer);

			viewer->setGridIsDrawn(ui->showGrids->isChecked());
			viewer->isDrawWireframe = ui->drawWireframe->isChecked();

			Eigen::AlignedBox3d bbox = m->bbox();

            //viewer->setGridIsDrawn();
			viewer->camera()->setSceneRadius(m->bbox().sizes().norm() * 2);
			viewer->camera()->setUpVector(qglviewer::Vec(0,0,1));
			viewer->camera()->setPosition(qglviewer::Vec(2,-2,1.5));
			viewer->camera()->lookAt(qglviewer::Vec());
			viewer->camera()->showEntireScene();

            double S = bbox.sizes().norm() * 0.1;
            bbox.extend(bbox.min() + (Vector3(-1,-1,-1) * S));
            bbox.extend(bbox.max() + (Vector3(1,1,1) * S));

            viewer->camera()->setRevolveAroundPoint(qglviewer::Vec(bbox.center()));
            viewer->camera()->setSceneCenter(qglviewer::Vec(bbox.center()));
			viewer->camera()->fitBoundingBox(qglviewer::Vec(bbox.min().data()),qglviewer::Vec(bbox.max().data()));

			connect(viewer, &MyDrawArea::gotFocus, [=](MyDrawArea * caller) {
				lastSelected = caller; 
				refreshViewers();
			});
		}
	}

	refreshViewers();

	this->activateWindow();
	this->setFocus();
	this->raise();
}

void MeshBrowser::refreshViewers(){
	for(auto viewer : activeViewers)
	{
		viewer->parentWidget()->setStyleSheet("border: 2px solid #232323;");
		if(deletedItems().contains(viewer->filename))
			viewer->parentWidget()->setStyleSheet("border: 2px solid red;");
	}

	if(lastSelected)
	{
		lastSelected->parentWidget()->setStyleSheet("border: 2px solid rgb(255,100,100);");
		ui->selectedMesh->setText( QFileInfo(lastSelected->filename).baseName() );
	}
}

MeshBrowser::~MeshBrowser()
{
    delete ui;
}
