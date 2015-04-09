#pragma once

#include <QString>
#include <QVector>
#include <QMap>
#include <QVariant>
#include <QImage>
#include <QPainter>
#include <QDir>
#include <QMessageBox>

inline void showImages(QVector<QImage> imgs, QStringList labels = QStringList()){ 
	QMessageBox msg;
	QImage img(imgs.size() * imgs.front().width(), imgs.front().height(), QImage::Format_RGBA8888_Premultiplied);
	QPainter painter(&img); QFont font("Monospace",7); font.setStyleHint(QFont::Monospace); painter.setFont(font);
	for(int i = 0; i < imgs.size(); i++){
		painter.drawImage(0,0,imgs[i]); 
		if(!labels.isEmpty()) {
			painter.setPen(QPen(Qt::white)); painter.drawText(QPoint(5,10), labels[i]);
			painter.setPen(QPen(Qt::red)); painter.drawText(QPoint(4,9), labels[i]);
		}
		painter.translate(imgs[i].width(),0); 
	}
	msg.setIconPixmap(QPixmap::fromImage(img)); msg.exec(); 
}

class ImageCompare
{
public:
    ImageCompare();

	struct Instance{
		int index;
		QString id, filename;
		std::vector< std::pair<double,double> > contour;
		std::vector<double> signature;
		QImage cachedImage;
		QImage image() { if(cachedImage.isNull()) cachedImage.load(filename); return cachedImage; }
		QMap<QString,QVariant> property;
		Instance(){ index = -1; }
		Instance(std::vector< std::pair<double,double> > contour);
	};

	struct DataSet{
		QString name, path;
		QVector<Instance> data;
		QMap<QString,QVariant> property;
	};

    QMap<QString,DataSet> datasets;
    QMap<QString,QVariant> property;

	void loadKnowledge( QString folderPath, QString datasetName = QString() );
	void addMoreKnowledge(QString datasetName, QString folderPath);
	void addInstance( QString datasetName, ImageCompare::Instance instance );
	static QVector<ImageCompare::Instance> loadDatafiles(QDir d, QStringList imageFiles);
	
	/// Signatures for fast look up
	static bool isClockwise( const std::vector< std::pair<double,double> > & contour );
	static std::vector<double> generateSignature( Instance instance, bool isSaveToFile = true );
	static std::pair<double,double> getCentroid( const std::vector< std::pair<double,double> > & contour );
	static std::vector<double> centroidDistanceSignature( const std::vector< std::pair<double,double> > & contour );
	static std::vector<double> fourierDescriptor( std::vector<double> cds );

	/// Distance measure
	static double distance(const ImageCompare::Instance &instanceA, const ImageCompare::Instance &instanceB);

	/// k-nearest neighbors
	typedef QVector< QPair<double, Instance > > InstanceMatches;
	InstanceMatches kNearest(const ImageCompare::Instance &instance, int k = 1, bool isReversed = false) const;

	/// Duplicated items
	QVector< QVector<Instance> > duplicateSets( double threshold, QString datasetName = QString() );
	void removeDuplicateSets( double threshold, QString datasetName = QString() );

	/// Utility
	Instance getInstance(QString datasetName, int idx = 0);
	size_t datasetSize(QString datasetName);

	/// Visualization
	static void showInstances( InstanceMatches instances );
	static QImage visualizeInstance( Instance instance, QString lable = QString() );
};
