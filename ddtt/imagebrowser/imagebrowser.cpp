#include "imagebrowser.h"
#include "ui_imagebrowser.h"
#include <QWidget>
#include <QPushButton>
#include <QFileDialog>
#include <QSettings>

#include "myimagearea.h"

ImageOperation curOp = NONE_OP;
MyImageArea * lastSelected = NULL;
QVector<MyImageArea*> activeViewers;

ImageBrowser::ImageBrowser(QWidget *parent) : QWidget(parent), ui(new Ui::ImageBrowser)
{
	// Anti-aliasing when using QGLWidget or subclasses
    //QGLFormat glf = QGLFormat::defaultFormat();
    //glf.setSamples(8);
    //QGLFormat::setDefaultFormat(glf);

    ui->setupUi(this);

    {
        QString datasetPath = QFileDialog::getExistingDirectory(0, "Dataset");

        // Get list of files
        QStringList filters;
        filters << "*.png" << "*.jpg" << "*.gif" << "*.jpeg" << "*.bmp";
        QStringList files = QDir(datasetPath).entryList(filters);

        foreach(QString mesh, files)
        {
            database << (datasetPath + "/" + mesh);
		}
	}

	// Settings
	{
		// Use local file
		QString path = QDir(qApp->applicationDirPath()).absoluteFilePath("imagebrowser.ini");
		QSettings qsettings(path, QSettings::IniFormat);

		// Defaults  
		if(!qsettings.allKeys().contains("nV")){
			qsettings.setValue("nV", 6);
			qsettings.setValue("nU", 3);
			qsettings.sync();
		}
 
		nV = qsettings.value("nV", 6).toInt();
		nU = qsettings.value("nU", 3).toInt();
	}

	// Create viewers
	{
		// Scroll bars
		int numPages = ceil(double(database.size()) / double(nU * nV)) - 1;
		ui->viewersScroll->setValue(0);
		ui->viewersScroll->setRange(0, numPages);
        connect(ui->viewersScroll, SIGNAL(valueChanged(int)), SLOT(loadImages()));

		// Load meshes
        loadImages();
	}

    // Connections
	connect(ui->flipImage, &QPushButton::released, [=]() {curOp = ImageOperation::FLIP; ui->curLabel->setText(imageOps[curOp]);});
	connect(ui->cropImg, &QPushButton::released, [=]() {curOp = ImageOperation::CROP; ui->curLabel->setText(imageOps[curOp]);});
	connect(ui->undoAll, &QPushButton::released, [=]() {curOp = ImageOperation::UNDO_ALL; ui->curLabel->setText(imageOps[curOp]);});
	connect(ui->thresholdButton, &QPushButton::released, [=]() {curOp = ImageOperation::THRESHOLD; ui->curLabel->setText(imageOps[curOp]);});
	connect(ui->clickDelete, &QPushButton::released, [=]() {curOp = ImageOperation::CLICK_DELETE; ui->curLabel->setText(imageOps[curOp]);});

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

	// Background removal
	{
		connect(ui->beginBackground, &QPushButton::released, [=](){curOp = ImageOperation::BACK_REMOVE; ui->curLabel->setText(imageOps[curOp]);  
			if(!lastSelected) return; lastSelected->isBackground = true; } );
		connect(ui->beginForeground, &QPushButton::released, [=](){curOp = ImageOperation::BACK_REMOVE; ui->curLabel->setText(imageOps[curOp]);  
			if(!lastSelected) return; lastSelected->isBackground = false; } );
		connect(ui->clearStrokes, &QPushButton::released, [=](){ if(!lastSelected) return; lastSelected->fg.clear(); lastSelected->bg.clear(); lastSelected->update(); } );

		connect(ui->removeBackground, &QPushButton::released, [=](){
			lastSelected->removeBackground( ui->paramSigma->value() );
			update();
		});

		connect(ui->thresholdToStrokesButton, &QPushButton::released, [=](){
			lastSelected->strokesFromThreshold( ui->thresholdVal->value() );
			update();
		});

		connect(ui->autoRemove1, &QPushButton::released, [=](){
			lastSelected->autoRemoveBack1();
			update();
		});

		connect(ui->autoRemove2, &QPushButton::released, [=](){
			lastSelected->autoRemoveBack2();
			update();
		});
	}
}

QStringList ImageBrowser::deletedItems(){
	QStringList items;
	for(int i = 0; i < ui->deleteList->count(); ++i){
		QListWidgetItem* item = ui->deleteList->item(i);
		items << item->text();
	}
	return items;
}

void ImageBrowser::loadImages()
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

            MyImageArea * viewer = new MyImageArea( filename, ui );

            if( !QFileInfo(filename).exists() ){
                viewer->isDeleted = true;
                continue;
            }

            QFrame * frame = new QFrame;
            QHBoxLayout * ll = new QHBoxLayout;
            frame->setLayout(ll);
            ll->addWidget(viewer);
            layout->setMargin(2);
            layout->setSpacing(2);
            layout->addWidget(frame, i, j, 1, 1);

            activeViewers.push_back(viewer);

            connect(viewer, &MyImageArea::gotFocus, [=](MyImageArea * caller) {
                lastSelected = caller;
                refreshViewers();
            });

			connect(viewer, &MyImageArea::deleteMe, [=](QString filename) {
				ui->addDelete->click();
			});
		}
	}

	refreshViewers();

	this->activateWindow();
	this->setFocus();
	this->raise();
}

void ImageBrowser::refreshViewers(){
	for(auto viewer : activeViewers)
	{
		viewer->parentWidget()->setStyleSheet("border: 2px solid #232323;");
		if(deletedItems().contains(viewer->filename))
			viewer->parentWidget()->setStyleSheet("border: 2px solid red;");
	}

	if(lastSelected)
	{
		lastSelected->parentWidget()->setStyleSheet("border: 2px solid rgb(255,100,100);");
        ui->selectedImage->setText( QFileInfo(lastSelected->filename).baseName() );
	}
}

ImageBrowser::~ImageBrowser()
{
    delete ui;
}
