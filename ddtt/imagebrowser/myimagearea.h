#pragma once
#include <QStringList>
#include <QColor>
#include <QWidget>
#include <QImage>

extern QString curLabel;
extern QStringList AllLabels;
extern QVector<QColor> UniqueColors;
extern QVector<QString> labelNames;

static QString imageOps[] = { "none", "undo-all", "flip", "crop", "threshold", "back-remove", "click-delete" };
enum ImageOperation{ NONE_OP, UNDO_ALL, FLIP, CROP, THRESHOLD, BACK_REMOVE, CLICK_DELETE };
extern ImageOperation curOp;

namespace Ui {class ImageBrowser;}

class MyImageArea : public QWidget
{
	Q_OBJECT
public:
    MyImageArea(QString filename, Ui::ImageBrowser *ui);
	~MyImageArea();

	void paintEvent(QPaintEvent *);
	void mouseMoveEvent(QMouseEvent * event);
	void mousePressEvent(QMouseEvent * event);
	void mouseReleaseEvent(QMouseEvent * event);

    QImage img;
    QString filename;
	bool isDeleted;
	bool isMouseClicked;
	QPoint startPos;
	Ui::ImageBrowser *ui;

	// Background removal
	std::vector< std::vector<bool> > fg, bg;
	bool isBackground;
	bool isScribbling;
	QPoint lastPoint;

    // Backup
    QImage orig_img;

	QPixmap checkerboard;
	QStringList debugTxt;

	void removeBackground( double sigma );
	void saveBackRemove();
	void runBackRemove( double sigma );
	void cleanupRemove();
	void strokesFromThreshold(double t);
	void autoRemoveBack1();
	void autoRemoveBack2();

signals:
	void gotFocus(MyImageArea *);
	void deleteMe(QString filename);
};

extern MyImageArea * lastSelected;

