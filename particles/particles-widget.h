#pragma once

#include <QWidget>

#include "ParticleMesh.h"

namespace Ui { class ParticlesWidget; }

class ParticlesWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ParticlesWidget(QWidget *parent = 0);
    ~ParticlesWidget();

	std::vector<ParticleMesh*> pmeshes;
	bool isReady;

//private:
    Ui::ParticlesWidget *ui;

signals:
	void shapesLoaded();
};
