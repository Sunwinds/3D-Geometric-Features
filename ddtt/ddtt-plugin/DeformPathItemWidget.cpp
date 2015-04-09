#pragma warning(disable:4267)

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QSlider>
#include <QLabel>
#include <QLineEdit>
#include <QFileDialog>
#include "DeformPathItemWidget.h"
#include "ShapeCorresponder.h"

DeformPathItemWidget::DeformPathItemWidget(int width, int height, DeformationPath * usedPath) : path(usedPath)
{
	// Create a widget with a slider and a progress bar
	QWidget * w = new QWidget;

	width = width - (height * 0.5);
	height = (1.0/5.0) * height;

	w->setMinimumSize(width, height);
	w->setMaximumSize(width, height);

	//w->setStyleSheet("*{ border: 1px solid red; }");

	// Set widget
	setWidget(w);

	connect(this, SIGNAL(widgetCreated()), SLOT(init()));
	emit( widgetCreated() );
}

void DeformPathItemWidget::init()
{
	QVBoxLayout * layout = new QVBoxLayout;

	layout->setMargin(0);
	layout->setSpacing(0);

	// Sub-layout
	{
		QHBoxLayout * sublayout = new QHBoxLayout;
		QHBoxLayout * sliderlayout = new QHBoxLayout;

		sublayout->setMargin(0);
		sublayout->setSpacing(0);

		sliderlayout->setMargin(0);
		sliderlayout->setSpacing(0);

		// Slider
		slider = new QSlider(Qt::Horizontal);
		slider->setTickPosition(QSlider::TicksBothSides);
		slider->setMinimum(0);
		slider->setMaximum(150);
		QString style = "QSlider::groove:horizontal {background: blue;height: 4px;}";
		style += "QSlider::handle:horizontal {background: #ff0000;width: 20px;margin: -16px 0px -16px 0px;}";
		slider->setStyleSheet( style );
		sliderlayout->addWidget(slider);

		// Save Correspondence button
		saveCorrButton = new QPushButton("Save..");
		sublayout->addWidget(saveCorrButton);
		connect(saveCorrButton, &QPushButton::clicked, [=](){ 
			QString dpath = QFileInfo(path->gcorr->sg->property["name"].toString()).absolutePath() + "/";
			QString g1n = path->gcorr->sg->name(), g2n = path->gcorr->tg->name();
			QString defaultFilename = g1n+"_"+g2n;
			QString filename = QFileDialog::getSaveFileName(0, "Save Correspondence", dpath + defaultFilename, "Correspondece files (*.txt)");
			path->gcorr->saveCorrespondences( filename );
		});

		// Render proxies
		proxyButton = new QPushButton("Proxies..");
		sublayout->addWidget(proxyButton);
		connect(proxyButton, &QPushButton::clicked, [=](){
			this->setEnabled(false);
			path->renderProxies();
			this->setEnabled(true);
		});

		// Render samples
		renderButton = new QPushButton("Render..");
		sublayout->addWidget(renderButton);
		connect(renderButton, &QPushButton::clicked, [=](){
			this->setEnabled(false);
			path->renderSamples();
			this->setEnabled(true);
		});

		// Execute button
		executeButton = new QPushButton("Execute..");
		sublayout->addWidget(executeButton);
		connect(executeButton, &QPushButton::clicked, [=](){ 
			this->setEnabled(false);
			//path->execute();
			path->renderProjected();
			this->setEnabled(true);
		});

		layout->addLayout( sliderlayout );
		layout->addLayout( sublayout );
	}
    
	// Messages
	layout->addWidget(label = new QLabel( QString("Path %1 (%2)").arg(path->idx).arg(path->i) ));
	//QLineEdit *numberEdit = new QLineEdit;layout->addWidget(numberEdit);

	widget()->setLayout(layout);

    // Connections
    connect(slider,SIGNAL(valueChanged(int)),this,SLOT(sliderValueChanged(int)));
}

void DeformPathItemWidget::sliderValueChanged(int val)
{
    label->setText( QString("Path %1 (%2) - %3").arg(path->idx).arg(path->i).arg(QString("inbetween: %1").arg(val)) );

	if( !path->scheduler.isNull() && path->scheduler->allGraphs.size() )
		path->si = (double(val) / slider->maximum()) * (path->scheduler->allGraphs.size()-1);
	else
		path->si = val;
}
