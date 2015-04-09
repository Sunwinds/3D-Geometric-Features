#include <QPushButton>
#include <QFileDialog>
#include "particles-widget.h"
#include "ui_particles-widget.h"

ParticlesWidget::ParticlesWidget(QWidget *parent) : QWidget(parent), ui(new Ui::ParticlesWidget), isReady(false)
{
    ui->setupUi(this);
}

ParticlesWidget::~ParticlesWidget()
{
	qDeleteAll(pmeshes);

    delete ui;
}
