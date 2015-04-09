#pragma once

#include <QWidget>
#include <QAbstractButton>

namespace Ui {
class ImageBrowser;
}

class ImageBrowser : public QWidget
{
    Q_OBJECT

public:
    explicit ImageBrowser(QWidget *parent = 0);
    ~ImageBrowser();
	Ui::ImageBrowser *ui;

	QStringList database;
	int nU, nV;

	QStringList deletedItems();

public slots:
	void loadImages();
	void refreshViewers();
};
