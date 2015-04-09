#pragma once

#include <QGraphicsProxyWidget>
#include <QSlider>
#include <QLabel>
#include <QPushButton>

class DeformationPath;

class DeformPathItemWidget : public QGraphicsProxyWidget
{
	Q_OBJECT
public:
	DeformPathItemWidget(int width, int height, DeformationPath * usedPath = NULL);

	DeformationPath * path;
	QWidget * w;

    QSlider * slider;
    QLabel * label;
	QPushButton * saveCorrButton;
	QPushButton * proxyButton;
	QPushButton * renderButton;
	QPushButton * executeButton;

public slots:
	void init();
    void sliderValueChanged(int val);

signals:
	void widgetCreated();
};
