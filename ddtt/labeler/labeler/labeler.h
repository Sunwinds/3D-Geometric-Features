#ifndef LABELER_H
#define LABELER_H

#include <QWidget>
#include <QAbstractButton>

namespace Ui {
class Labeler;
}

class Labeler : public QWidget
{
    Q_OBJECT

public:
    explicit Labeler(QWidget *parent = 0);
    ~Labeler();
	Ui::Labeler *ui;

	QStringList database;
	int nU, nV;

public slots:
	void onGroupButtonClicked(QAbstractButton * button);
	void loadMeshes();
};

#endif // LABELER_H
