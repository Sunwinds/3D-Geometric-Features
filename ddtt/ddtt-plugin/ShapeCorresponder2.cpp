#include "ShapeCorresponder2.h"
#include "PathsGenerator.h"

#include "BofSearchManager.h"
QString knowledgeDataset = "C:/Shared/all_images_clean02";

#include "ProjectedStructureGraph.h"

ShapeCorresponder2::ShapeCorresponder2(Structure::Graph *g1, Structure::Graph *g2) : source(g1), target(g2)
{
    QElapsedTimer prepareTimer; prepareTimer.start();

    // Load knowledge
    bmanager = new BofSearchManager( knowledgeDataset, true );
	 
    // Generate possible paths
    int k = 2;
    paths = PathsGenerator( source, target, k ).paths;

    // Timing & Stats
    property["pathsCount"].setValue( (int)paths.size() );
    property["prepareTime"].setValue( (int)prepareTimer.elapsed() );

    // Progress bar
    pd = new QProgressDialog( "Evaluating Paths..", "Cancel", 0, (int)paths.size() );
    pd->setValue(0); pd->show();
    this->connect( this, SIGNAL(pathComputed()), SLOT(increaseProgress()) );
}

void ShapeCorresponder2::run()
{
    QElapsedTimer computeTimer; computeTimer.start();

    bool abort = false;
	bool isVisualizeProcess = true;

	/// Parameters
	int numSamples = 5;
	int k_neighbours = 6;

	/// Preparation / clean up
	QString jobPath = QFileInfo(source->property["correspondenceFile"].toString()).absolutePath() + "/";
	if(jobPath.length() < 2) jobPath = "";
	QDir dir(jobPath);
	foreach(QString pngfile, dir.entryList( QStringList() << "*.*" )){
		if(pngfile.endsWith("p.png"))
			dir.remove(pngfile);
	}

	// Using projected source and target graphs
	ProjectedStructureGraph pgSource(source, 128);
	ProjectedStructureGraph pgTarget(target, 128);

	// Weighted comparisons based on neighbours of source and target
	std::vector<dist_idx_t> sneighbours = bmanager->search( QImageToCvMat(pgSource.drawBoundaryImage()), true );
	std::vector<dist_idx_t> tneighbours = bmanager->search( QImageToCvMat(pgTarget.drawBoundaryImage()), true );
	std::vector<double> sweights(sneighbours.size(),1.0), tweights(sneighbours.size(),1.0);
	for(size_t i = 0; i < sneighbours.size(); i++){
		sweights[ sneighbours[i].second ] = std::max(1e-12, sneighbours[i].first);
		tweights[ tneighbours[i].second ] = std::max(1e-12, tneighbours[i].first);
	}

    // Evaluate deformations
    #pragma omp parallel for
    for(int pi = 0; pi < (int)paths.size(); pi++)
    {
        #pragma omp flush (abort)
        if(!abort)
        {
            if(pd->wasCanceled()){
                abort = true;
                #pragma omp flush (abort)
            }

            auto & path = paths[pi];

			// DEBUG:
			QImage pathImage;
			int pathImageWidth = 1280;
			if(isVisualizeProcess) pathImage = QImage( pathImageWidth, pathImageWidth * 0.4, QImage::Format_ARGB32_Premultiplied );

			// Sample the path
			for(int s = 0; s < numSamples; s++)
			{
				double t = double(s)/(numSamples-1);
				t = ((1.0 - 0.6) / 2.0) + (t * 0.6); // middle 60%

				// Draw blended shape
				QImage inbetween = pgSource.drawBlendedImage(&pgTarget, path.gcorr, t);

				// Compare with knowledge
				cv::Mat cvinbetween = QImageToCvMat(inbetween);
				std::vector<dist_idx_t> neighbours = bmanager->search( cvinbetween );
				neighbours.resize( k_neighbours );

				double error = 0.0;

				for( auto & against : neighbours )
				{
					index_t idx = against.second;
					double similarity = against.first * std::max(sweights[idx], tweights[idx]);
					error += 1.0 / similarity;
				}

				path.weight = error;

				/// DEBUG:
				if( isVisualizeProcess )
				{
					QPainter painter(&pathImage);

					int w = pathImageWidth / numSamples;
					painter.drawImage(s * w, 0, cvMatToQImage(cvinbetween));
					painter.drawRect(inbetween.rect().translated(QPoint(s*w,0)));

					int top = inbetween.height(), left = s * w, padding = 5;
					int x = left, y = top + padding;

					for( auto & against : neighbours )
					{
						index_t idx = against.second;
						double similarity = against.first * std::max(sweights[idx], tweights[idx]);
						
						QString cost = QString::number(similarity, 'g', 3);

						QImage img = bmanager->cachedImages[against.second];
						img = img.scaledToWidth(w * 0.47, Qt::SmoothTransformation);

						QPainter pimg(&img);
						pimg.drawText(10,10, cost);
						pimg.fillRect( QRect(0,img.height()-10,img.width(), 5), starlab::qtJetColor(similarity, 0.0, 0.5));
	
						painter.drawImage(x, y, img);

						x += img.width(); 
						if(x > (w + left) - img.width()) { x = left; y += img.height() + padding; }
					}
				}
			}

			if( isVisualizeProcess && !pathImage.isNull() )
			{
				QString filename = QString("%2-path-%1p.png").arg(pi).arg( QString::number(path.weight, 'g', 4) );

				// Special mark ground truth correspondence
				if(path.property.contains("isGroundTruth") && path.property.value("isGroundTruth").toBool()){
					QPainter painter(&pathImage);
					painter.setPen(QPen(Qt::red, 10));
					painter.drawRect(pathImage.rect());
				}

				pathImage.save( jobPath + filename );
			}

            emit( pathComputed() );
        }
    }

    // Timing
    property["computeTime"].setValue( (int)computeTimer.elapsed() );

    if( !paths.size() ){
        emit( done() );
        return;
    }

    // Find best
    {
        // Best = lowest error
        std::sort(paths.begin(), paths.end(), DeformationPathCompare);
    }

    emit( done() );
}

void ShapeCorresponder2::increaseProgress()
{
    if(!pd->wasCanceled()) pd->setValue( pd->value() + 1 );
}
