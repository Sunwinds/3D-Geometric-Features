#include <omp.h>

#include "ImageCompare.h"
#include <QDir>
#include <QFileInfo>
#include <QDebug>

#include "fft.h"
#include "QC.h"

ImageCompare::ImageCompare()
{
}

void ImageCompare::loadKnowledge(QString folderPath, QString datasetName)
{
	QDir d(folderPath);

	if(!d.exists()) {
		qDebug() << "Warning: Cannot find folder " << folderPath;
		return;
	}

	if(!datasetName.length()) datasetName = d.dirName();

	DataSet dataset;
	dataset.name = datasetName;
	dataset.path = folderPath;

	QStringList imageFiles = d.entryList( QStringList() << "*.png" );

	QVector<ImageCompare::Instance> data = ImageCompare::loadDatafiles( d, imageFiles );
	dataset.data = data;

	datasets[datasetName] = dataset;
}

void ImageCompare::addMoreKnowledge(QString datasetName, QString folderPath)
{
	if( !datasets.contains(datasetName) ){
		loadKnowledge(folderPath, datasetName);
		return;
	}

	QDir d(folderPath);
	datasets[datasetName].data << loadDatafiles( d, d.entryList( QStringList() << "*.png" ) );
}

void ImageCompare::addInstance(QString datasetName, Instance instance)
{
	if( instance.index < 0 ) instance.index = datasetSize(datasetName);
	if( instance.id.isNull() || instance.id.isEmpty() ) instance.id = QString::number(instance.index);
	if( instance.signature.empty() ) instance.signature = ImageCompare::generateSignature( instance.contour );

	datasets[datasetName].data.push_back( instance );
}

QVector<ImageCompare::Instance> ImageCompare::loadDatafiles(QDir d, QStringList imageFiles)
{
	QVector<ImageCompare::Instance> data( imageFiles.size() );

	// Go over image files in the folder
	#pragma omp parallel for
	for(int i = 0; i < imageFiles.size(); i++)
	{
		QString file = imageFiles[i];
		QString fullFilename = d.absolutePath() + "/" + file;

		ImageCompare::Instance inst;
		inst.index = i;
		inst.filename = fullFilename;
		inst.id = QFileInfo(inst.filename).baseName();

		// Check for pre computed contour file
		QString contourFilename = inst.id + ".txt";
		if( d.exists(contourFilename) )
		{
			// Load from disk
			QFile contourFile( d.absolutePath() + "/" + contourFilename );
			QTextStream in( &contourFile );
			if( contourFile.open(QIODevice::ReadOnly | QIODevice::Text) )
			{
				while( !in.atEnd() )
				{
					QStringList data = in.readLine().split(" ", QString::SkipEmptyParts);
					if(data.isEmpty() || data.size() < 2) continue;

					inst.contour.push_back( std::make_pair(data[0].toDouble(), data[1].toDouble()) );
				}
			}

			// Check orientation
			if( isClockwise(inst.contour) ) std::reverse( inst.contour.begin(), inst.contour.end() );
		}
		else
		{
			// Compute contour..
		}

		// Check for signatures
		QString sigFilename = inst.id + ".sig";
		if( d.exists(sigFilename) )
		{
			// Load from disk
			QFile sigFile( d.absolutePath() + "/" + sigFilename );
			QTextStream in( &sigFile );
			if( sigFile.open(QIODevice::ReadOnly | QIODevice::Text) )
			{
				QStringList signature = in.readLine().split(" ", QString::SkipEmptyParts);
				for(auto s : signature) inst.signature.push_back( s.toDouble() );
			}
		}
		else
		{
			// Compute signature..
			inst.signature = generateSignature(inst, true);
		}

		data[i] = inst;
	}

	QVector<ImageCompare::Instance> result = data;
	return result;
}

bool ImageCompare::isClockwise(const std::vector< std::pair<double,double> > & contour)
{
	double sum = 0.0;
	for (size_t i = 0; i < contour.size(); i++) {
		std::pair<double,double> v1 = contour[i];
		std::pair<double,double> v2 = contour[(i + 1) % contour.size()];
		sum += (v2.first - v1.first) * (v2.second + v1.second);
	}
	return sum > 0.0;
}

std::pair<double,double> ImageCompare::getCentroid(const std::vector< std::pair<double,double> > & contour)
{
	double accumulatedArea = 0.0;
	double centerx = 0.0;
	double centery = 0.0;

	for (size_t i = 0, j = contour.size() - 1; i < contour.size(); j = i++)
	{
		double temp = contour[i].first * contour[j].second - contour[j].first * contour[i].second;
		accumulatedArea += temp;
		centerx += (contour[i].first + contour[j].first) * temp;
		centery += (contour[i].second + contour[j].second) * temp;
	}

	// return midpoint if area is almost zero..
	if (abs(accumulatedArea) < 1e-8){
		std::pair<double,double> c(0,0);
		for(auto p : contour) { c.first += p.first ; c.second += p.second; }
		c.first /= contour.size(); c.second /= contour.size();
		return c;
	}

	accumulatedArea *= 3.0;
	return std::pair<double,double>(centerx / accumulatedArea, centery / accumulatedArea);
}

std::vector<double> ImageCompare::centroidDistanceSignature( const std::vector< std::pair<double,double> > & contour )
{
	// Find centroid
	std::pair<double,double> c = getCentroid( contour );

	std::vector<double> dists;
	for(auto p : contour) 
	{
		double dist = std::sqrt( pow(p.first - c.first,2) + pow(p.second - c.second,2) );
		dists.push_back( dist );
	}

	return dists;
}

std::vector<double> ImageCompare::fourierDescriptor( std::vector<double> cds )
{
	if(!cds.size()) return cds;

	int sig_length = 256;

	// FFT on Centroid Distance signature
	std::vector<double> real = cds, imag(real.size(), 0);
	Fft::transform(real, imag);

	bool isNormalizeFD = true;
	if( isNormalizeFD )
	{
		// Translation invariance
		real[0] = 0;
		imag[1] = 0;

		// Scale invariance
		double si = std::abs( std::complex<double>(real[1], imag[1]) );
		for(size_t i = 0; i < real.size(); i++)
		{
			real[i] /= si;
			imag[i] /= si;
		}

		// Rotation and changes in starting point
		for(size_t i = 0; i < real.size(); i++)
		{
			real[i] = std::abs( std::complex<double>(real[i], imag[i]) );
			imag[i] = 0;
		}
	}

	// Check against signature length
	if( real.size() < sig_length)
	{
		size_t padding = sig_length - real.size();
		size_t h = real.size() / 2;
		real.insert(real.begin() + h, padding, 0);
		imag.insert(imag.begin() + h, padding, 0);
	}

	// Flip
	Fft::fftshift(real, imag);

	// Shorten to assigned size 
	Fft::shrink(sig_length, real, imag);

	return real;
}

std::vector<double> ImageCompare::generateSignature( Instance inst, bool isSaveToFile )
{		
	if(!inst.contour.size()) return std::vector<double>();

	std::vector<double> sig = fourierDescriptor( centroidDistanceSignature(inst.contour) );

	if( isSaveToFile )
	{
		// Save to disk
		QString sigFilename = inst.id + ".sig";
		QFile sigFile( QFileInfo(inst.filename).absolutePath() + "/" + sigFilename );
		QTextStream out( &sigFile );

		if( sigFile.open(QIODevice::WriteOnly | QIODevice::Text) )
			for(auto s : sig) out << QString::number(s) << " ";
	}

	return sig;
}

double ImageCompare::distance(const Instance &instanceA, const Instance &instanceB)
{
	/* Quadratic-Chi Histogram Distance */
	//double dist = QC::distance( j.signature, instance.signature );

	/* L2 Euclidean distance */
	double dist = QC::L2distance( instanceA.signature, instanceB.signature );	

	return dist;
}

ImageCompare::InstanceMatches ImageCompare::kNearest(const Instance &instance, int k, bool isReversed) const
{
	ImageCompare::InstanceMatches result;

	typedef QPair<double, int> ScoreInstance;
	QVector< ScoreInstance > candidates;

	for(auto key : datasets.keys())
	{
		const DataSet & dataset = datasets[key];

		for(auto & j : dataset.data)
		{
			candidates << ScoreInstance( ImageCompare::distance(j, instance), j.index );
		}

		// Sort...
		std::sort( candidates.begin(), candidates.end(), [](const ScoreInstance& a, const ScoreInstance& b){ return a.first < b.first; } );

		// To get 'k' furthest, useful for debugging
		if(isReversed) std::reverse(candidates.begin(), candidates.end());

		// bound check
		int K = std::min(k, datasets[key].data.size());

		// Return first 'k'
		for(int i = 0; i < K; i++) result.push_back( qMakePair(candidates[i].first, datasets[key].data.at( candidates[i].second )) );
	}

	return result;
}

ImageCompare::Instance ImageCompare::getInstance( QString datasetName, int idx )
{
	return datasets[datasetName].data.at(idx);
}

size_t ImageCompare::datasetSize(QString datasetName)
{
	return datasets.value(datasetName).data.size();
}

QVector< QVector<ImageCompare::Instance> > ImageCompare::duplicateSets( double threshold, QString datasetName )
{
	if(!datasetName.length() || !datasets.keys().contains(datasetName))	datasetName = datasets.keys().front();
	DataSet & dataset = datasets[datasetName];

	int N = dataset.data.size();

	typedef QPair<double, int> ScoreInstance;
	QVector< QVector< ScoreInstance > > candidates(N);

	for(int i = 0; i < N; i++)
	{
		Instance & instance = dataset.data[i];

		for(int j = 0; j < N; j++)
		{
			if(i==j) continue;

			candidates[i].push_back( qMakePair(distance(instance, dataset.data.at(j)), j) );
		}

		// Sort
		std::sort( candidates[i].begin(), candidates[i].end(), 
			[](const ScoreInstance& a, const ScoreInstance& b){ return a.first < b.first; } );
	}

	QVector< QVector<ImageCompare::Instance> > result;

	for(int i = 0; i < N; i++)
	{
		std::vector<int> similar;
		for(auto p : candidates[i]) if(p.first < threshold) similar.push_back( p.second );

		QVector<ImageCompare::Instance> mysimilar;
		for(auto s : similar) mysimilar.push_back(dataset.data.at(s));

		if( !mysimilar.empty() ) 
		{
			mysimilar.push_front( dataset.data.at(i) );
			result.push_back( mysimilar );
		}
	}

	return result;
}

void ImageCompare::removeDuplicateSets( double threshold, QString datasetName )
{
	if(!datasetName.length() || !datasets.keys().contains(datasetName))	datasetName = datasets.keys().front();
	QString datasetPath = datasets[datasetName].path;

	QStringList filenamesToRemove;

	for(auto set : duplicateSets(threshold, datasetName)){
		set.removeFirst();
		for(auto instance : set) filenamesToRemove << instance.filename;	
	}
	if( filenamesToRemove.empty() ) return;

	if(QMessageBox::question(0, "Delete duplicates", QString("Are you sure you want to delete (%1) files?").arg(
		filenamesToRemove.size()), QMessageBox::Yes|QMessageBox::No) == QMessageBox::Yes){
			for(auto filename : filenamesToRemove)
				QFile::remove( filename );
	}

	// Reload
	loadKnowledge(datasetPath, datasetName);
}

void ImageCompare::showInstances(InstanceMatches instances)
{
	QVector<QImage> imgs;
	QStringList labels;

	for(auto match : instances){
		Instance & i = match.second;
		imgs.push_back( i.image() );
		labels.push_back( QString("%1 : %2").arg(i.index).arg(i.id) );
	}

	showImages(imgs, labels);
}

QImage ImageCompare::visualizeInstance( ImageCompare::Instance instance, QString lable )
{
	QImage img = instance.image().convertToFormat(QImage::Format_ARGB32_Premultiplied);

	QPainter painter( &img ); 
	QPainterPath qpath;	int i = 0;
	for(auto p : instance.contour) (i++ == 0) ? qpath.moveTo(p.first, p.second) : qpath.lineTo(p.first, p.second);

	painter.setPen(QPen(Qt::blue, 2)); painter.drawPath(qpath);

	painter.fillRect(img.rect().translated(QPoint(0,-80)), QColor(255,255,255,180));

	QString title = QString("%1:%2").arg(instance.index).arg(instance.id);
	painter.setFont(QFont("Courier", 9));
	painter.drawText(QPoint(10,15), title);

	QStringList lables = lable.split("\n");
	int y = 25;

	for(auto lable : lables){
		painter.drawText(QPoint(10,y), lable);
		y += 10;
	}

	return img;
}

ImageCompare::Instance::Instance(std::vector< std::pair<double,double> > contour) : contour(contour), index(-1)
{
	// Check orientation
	if( isClockwise(this->contour) ) std::reverse( this->contour.begin(), this->contour.end() );

	signature = ImageCompare::fourierDescriptor( ImageCompare::centroidDistanceSignature( contour ) );
}
