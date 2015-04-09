#include "myimagearea.h"
#include "ui_imagebrowser.h"

#include <QPainter>
#include <QImageReader>
#include <QMouseEvent>
#include <QtConcurrent/QtConcurrent>
#include <QMessageBox>

QImage centerWithWhiteBackground( const QImage & orig_img )
{
	int w = 256;
	int maxDim = qMax(orig_img.width(), orig_img.height());
	if(maxDim < 200) w = 128;

	QImage centered(w, w, QImage::Format_ARGB32_Premultiplied);
	centered.fill(QColor(255,255,255));
	QPainter pimg(&centered);

	QImage scaledImg = orig_img.scaled(w,w,Qt::KeepAspectRatio, Qt::SmoothTransformation);

	QRect rect(scaledImg.rect());
	QRect devRect(0, 0, pimg.device()->width(), pimg.device()->height());
	rect.moveCenter(devRect.center());
	pimg.drawImage(rect.topLeft(), scaledImg);

	return centered;
}

QImage noBackgroundImage( const QImage & orig_img ){
	int w = orig_img.width();
	int h = orig_img.height();

	QImage replaced(w, h, QImage::Format_ARGB32_Premultiplied);
	replaced.fill(QColor(0,0,0,0));

	int emptyPixels = 0;

	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			QRgb pixel = orig_img.pixel(x,y);
			if(pixel == qRgb(255,0,255)) {
				emptyPixels++;
				continue;
			}
			replaced.setPixel(x,y, pixel);
		}
	}

	if(emptyPixels > 0.99 * w * h) return orig_img;

	return replaced;
}

enum COLOR_CHANNEL{ CH_RED, CH_GREEN, CH_BLUE, CH_VALUE, CH_SAT, CH_HUE, CH_LIGHT };
int NUM_CHANNELS = 7;

std::vector<int> channels(QColor c){
	std::vector<int> ch(NUM_CHANNELS,0);
	ch[CH_RED	]	= c.red();
	ch[CH_GREEN	]	= c.green();
	ch[CH_BLUE	]	= c.blue();
	ch[CH_VALUE	]	= c.value();
	ch[CH_SAT	]	= c.saturation();
	ch[CH_HUE	]	= qMax(0,c.hue());
	ch[CH_LIGHT	]	= c.lightness();
	return ch;
}

struct ImageRef {
	QRgb *data;
	int width;
	std::vector< std::vector<int> > channels;
	ImageRef(QRgb *data, int width) : data(data), width(width) {}
};

std::vector< std::vector<unsigned int> > getHistograms( const QImage & orig_img )
{
	std::vector< std::vector<unsigned int> > histograms(NUM_CHANNELS, std::vector<unsigned int>(256, 0));
	histograms[CH_HUE].resize(360, 0);
	for(int x=0;x<orig_img.width();x++){
		for(int y=0;y<orig_img.height();y++){
			std::vector<int> ch = channels(QColor(orig_img.pixel(x,y)));
			for(int i = 0; i < (int)ch.size(); i++)
				histograms[i][ch[i]]++;
		}
	}
	return histograms;
}

std::vector< std::vector<unsigned int> > colorLimits( const QImage & orig_img )
{
	std::vector< std::vector<unsigned int> > limits(NUM_CHANNELS, std::vector<unsigned int>(2, 0));

	for(int i = 0; i < NUM_CHANNELS; i++){
		limits[i][0] = 360;
		limits[i][1] = 0;
	}

	for(int x=0;x<orig_img.width();x++){
		for(int y=0;y<orig_img.height();y++){
			QColor c(orig_img.pixel(x,y));
			std::vector<int> ch = channels(c);
			for(int i = 0; i < (int)ch.size(); i++){
				limits[i][0] = qMin((int)limits[i][0], ch[i]);
				limits[i][1] = qMax((int)limits[i][1], ch[i]);
			}
		}
	}

	return limits;
}

#include <numeric>
std::vector< std::vector<unsigned int> > histogramCDF( const std::vector< std::vector<unsigned int> > & histogram )
{
	std::vector< std::vector<unsigned int> > CDF(NUM_CHANNELS, std::vector<unsigned int>(256, 0));
	for(int c = 0; c < NUM_CHANNELS; c++)
		for(int i=0;i<256;i++)
			for(int j=0;j<=i;j++)
				CDF[c][i] = histogram[c][j] + CDF[c][i];
	return CDF;
}

QImage calcHistImage( const QImage & orig_img )
{
	int width = qMax(300, int(orig_img.width() * 1.5));
	int height = orig_img.height() * 0.2;
	QImage histImage(width,height,QImage::Format_ARGB32_Premultiplied);
	histImage.fill(QColor(255,255,255));
	QPainter hpainter(&histImage);

	int channel = CH_VALUE;

	std::vector<unsigned int> histogram = getHistograms(orig_img)[channel];
	double maxHist = *std::max_element(histogram.begin(), histogram.end());

	int binWidth = qMax(2, int(histImage.width() / histogram.size()));

	for(size_t h = 0; h < histogram.size(); h++){
		double value = histogram[h] == 0.0 ? 1 : histogram[h];
		int y = height * (log(value) / log(maxHist));
		int x = int(h) * binWidth;
		int binHeight = y;
		QRect binRect(x, 0, binWidth, binHeight);
		binRect.moveBottom(histImage.height());
		hpainter.fillRect(binRect, QBrush(QColor(0,0,255)));
	}
	return histImage;
}

QImage levelImageValue( const QImage & orig_img, double low = 0.0, double high = 1.0 )
{
	QImage tempImage = orig_img;

	std::vector< std::vector<unsigned int> > histograms = getHistograms(orig_img);

	double newRange = high - low;

	for(int x=0;x<orig_img.width();x++){
		for(int y=0;y<orig_img.height();y++){
			QColor c(orig_img.pixel(x,y));
			std::vector<int> ch = channels(c);

			double v = (double)ch[CH_LIGHT] / 255;

			if(v > high) v = 1.0;
			else if(v < low) v = 0;
			else v = (v - low) / newRange;

			int newVal = 255 * v;

			tempImage.setPixel( x, y, QColor::fromHsl(ch[CH_HUE], ch[CH_SAT], newVal).rgb() );
		}
	}

	return tempImage;
}

QImage thresholdImage(const QImage& orig_img, double threshold = 0.5)
{
	QImage tempImage = orig_img;

	for(int x=0;x<orig_img.width();x++){
		for(int y=0;y<orig_img.height();y++){
			QColor c(orig_img.pixel(x,y));
			std::vector<int> ch = channels(c);

			double v = (double)ch[CH_LIGHT] / 255;

			if(v > threshold) 
				tempImage.setPixel( x, y, qRgb(255,255,255) );
			else
				tempImage.setPixel( x, y, qRgb(0,0,0) );
		}
	}

	return tempImage;
}

MyImageArea::MyImageArea(QString filename, Ui::ImageBrowser *ui) : filename(filename), ui(ui)
{
	isDeleted = false; 
	isMouseClicked = false;
	isBackground = true;
	isScribbling = false;

	orig_img.load(filename);

	// Saved with wrong extension
	if(orig_img.isNull())
	{
		orig_img.load(filename, "jpg");

		// fall back:
		QList<QByteArray> formats = QImageReader::supportedImageFormats();
		size_t f = 0;
		while(orig_img.isNull() && f < formats.size())
			orig_img.load(filename, formats.at(int(f++)));
	}

	img = centerWithWhiteBackground(orig_img);

	// Background
	checkerboard = QPixmap(20, 20);
	QPainter pmp(&checkerboard);

	QColor lightColor(35,40,35);
	QColor darkColor(40,45,40);

	pmp.fillRect(0, 0, 10, 10, darkColor);
	pmp.fillRect(10, 10, 10, 10, darkColor);
	pmp.fillRect(0, 10, 10, 10, lightColor);
	pmp.fillRect(10, 0, 10, 10, lightColor);
	pmp.end();

	setMouseTracking(true);
}

MyImageArea::~MyImageArea()
{
    // Save changed image
    img.save( filename );
}

void MyImageArea::paintEvent(QPaintEvent *)
{
    QPainter painter(this);

	painter.fillRect(rect(), QBrush(checkerboard));

    painter.drawImage(0, 0, img);

	// Border
	//painter.setPen(QPen(Qt::green, 1));
	//painter.drawRect(img.rect());

	/// Background removal stuff
	{
		// Draw strokes
		int w = img.width();
		int h = img.height();

		// Prepare overlays for background and foreground strokes
		QImage fgOverlay(w,h, QImage::Format_ARGB32_Premultiplied);
		fgOverlay.fill( QColor(0,0,0,0) );
		QImage bgOverlay = fgOverlay;

		// Draw the strokes
		for(int y = 0; y < h; y++){
			for(int x = 0; x < w; x++){
				if(bg.size() && bg[y][x]) bgOverlay.setPixel(x,y, qRgb(255,0,0));
				if(fg.size() && fg[y][x]) fgOverlay.setPixel(x,y, qRgb(0,0,225));
			}
		}

		painter.drawImage(0, 0, bgOverlay);
		painter.drawImage(0, 0, fgOverlay);
	}

	/// DEBUG:
	if( false )
	{
		painter.setBrush(Qt::black);
		int y = 10, lineHeight = 12;

		for(auto s : debugTxt)
		{
			painter.drawText(10, (y += lineHeight), s);
		}
	}
}

void MyImageArea::mousePressEvent(QMouseEvent *event)
{
	QWidget::mousePressEvent(event);

	if (event->button() == Qt::LeftButton) {
		lastPoint = event->pos();
		isScribbling = true;
	}

	emit( gotFocus(this) );

	isMouseClicked = true;
	startPos = event->pos();

	if(curOp == ImageOperation::CROP)
	{
		int w = img.width(), h = img.height();

		if(event->buttons() & Qt::LeftButton)
		{
			double dx = double(event->pos().x()) / w;
			img = img.copy(0, 0, w * dx, h);
		}

		if(event->buttons() & Qt::RightButton)
		{
			double dx = double(event->pos().x()) / w;
			int wx = w * (1.0 - dx);

			img = img.copy(w - wx, 0, wx, h);
		}

		img = centerWithWhiteBackground(img);
	}

	if(!(event->buttons() & Qt::RightButton)){
		update();
		return;
	}

	if(curOp == ImageOperation::FLIP)		img = img.mirrored(true,false);
	if(curOp == ImageOperation::THRESHOLD)	img = thresholdImage(img, ui->thresholdVal->value());
	if(curOp == ImageOperation::UNDO_ALL)	img = centerWithWhiteBackground(orig_img);
	if(curOp == ImageOperation::CLICK_DELETE) emit( deleteMe(filename) );

	update();
}

void MyImageArea::mouseReleaseEvent(QMouseEvent *event)
{
	QWidget::mouseReleaseEvent(event);

	if (event->button() == Qt::LeftButton && isScribbling) {
		isScribbling = false;
	}

	isMouseClicked = false;

	update();
}


void horizontalLine(std::vector< std::vector<bool> > & m, bool color,  int xpos, int ypos, int x1){
	for(int x = xpos; x <= x1; ++x){
		if(x < 0 || x > m.front().size() - 1 || ypos < 0 || ypos > m.size() - 1) continue;
		m[ypos][x] = color;
	}
}

void plot4points(std::vector< std::vector<bool> > & buffer, bool color, int cx, int cy, int x, int y){
	horizontalLine(buffer, color, cx - x, cy + y, cx + x);
	if (x != 0 && y != 0) horizontalLine(buffer, color, cx - x, cy - y, cx + x);
}

void circle(std::vector< std::vector<bool> > & buffer, bool color, int cx, int cy, int radius){
	int error = -radius;
	int x = radius;
	int y = 0;

	while (x >= y){
		int lastY = y;
		error += y;
		++y;
		error += y;
		plot4points(buffer, color, cx, cy, x, lastY);
		if (error >= 0){
			if (x != lastY)
				plot4points(buffer, color, cx, cy, lastY, x);
			error -= x;
			--x;
			error -= x;
		}
	}
}

void DrawLine( float x1, float y1, float x2, float y2, std::vector< std::vector<bool> > & img, int radius = 10, bool color = true )
{
	// Bresenham's line algorithm
	const bool steep = (std::abs(y2 - y1) > std::abs(x2 - x1));
	if(steep){
		std::swap(x1, y1);
		std::swap(x2, y2);
	}
	if(x1 > x2){
		std::swap(x1, x2);
		std::swap(y1, y2);
	}

	const float dx = x2 - x1;
	const float dy = fabs(y2 - y1);

	float error = dx / 2.0f;
	const int ystep = (y1 < y2) ? 1 : -1;
	int y = (int)y1;

	const int maxX = (int)x2;

	if(!img.size()) return;

	for(int x=(int)x1; x<maxX; x++){
		if(steep)
			circle(img, color, y, x, radius);
		else
			circle(img, color, x, y, radius);
		error -= dy;
		if(error < 0){
			y += ystep;
			error += dx;
		}
	}
}

void MyImageArea::mouseMoveEvent(QMouseEvent *event)
{
	if ((event->buttons() & Qt::LeftButton) && isScribbling)
	{
		int prev_y = lastPoint.y();
		int prev_x = lastPoint.x();

		int y = event->pos().y();
		int x = event->pos().x();

		int w = img.width();
		int h = img.height();

		if(!bg.size()) bg.resize(h, std::vector<bool>(w, false));
		if(!fg.size()) fg.resize(h, std::vector<bool>(w, false));

		if(isBackground) DrawLine(prev_x,prev_y, x, y, bg, ui->brushRadius->value(), true);
		else DrawLine(prev_x,prev_y, x, y, fg, ui->brushRadius->value(), true);
	}

	lastPoint = event->pos();

	// DEBUG:
	debugTxt.clear();
	debugTxt << QString("Mouse pos: %1, %2 (click %3)").arg(event->pos().x()).arg(event->pos().y()).arg(isMouseClicked);

	update();
}

void MyImageArea::saveBackRemove()
{
	if(!bg.size() || !fg.size()) {qDebug() << "empty strokes.."; return;}

	QString path = QFileInfo( QCoreApplication::applicationFilePath() ).absolutePath() + "/";
	QString f = QFileInfo(filename).baseName();
	QStringList filenames;
	filenames << (path + QString("%1.bg.dat").arg(f)) << (path + QString("%1.fg.dat").arg(f));

	QString copyFile = path + f + filename.right(4);
	img.save(copyFile);

	std::vector<QFile*> file;
	std::vector<QTextStream*> out;
	for(int i = 0; i < 2; i++){
		file.push_back(new QFile(filenames[i]));
		if (!file.back()->open(QIODevice::WriteOnly | QIODevice::Text)) return;
		out.push_back(new QTextStream(file.back()) );
	}

	int w = img.width(), h = img.height();
	(*out.front()) << QString("%1 %2\n").arg(h).arg(w);
	(*out.back()) << QString("%1 %2\n").arg(h).arg(w);
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			if(bg[y][x]) (*out.front()) << QString("%1 %2\n").arg(y).arg(x);
			if(fg[y][x]) (*out.back()) << QString("%1 %2\n").arg(y).arg(x);
		}
	}

	qDeleteAll(out);
	qDeleteAll(file);
}

void MyImageArea::runBackRemove(double sigma)
{
	QString path = QFileInfo( QCoreApplication::applicationFilePath() ).absolutePath() + "/";
	QString f = QFileInfo(filename).baseName();
	QStringList filenames;
	filenames << QString("%1.bg.dat").arg(f) << QString("%1.fg.dat").arg(f);
	
	QStringList options;

	// default
	{
		options << QString("-v %1").arg( 0 ); // UseLUV
		options << QString("-x %1").arg( 1 ); // Down-sample
		options << QString("-d %1").arg( 5 ); // Dimensions
		options << QString("-D %1").arg( "EUC" ); // Dist metric
		options << QString("-n %1").arg( ui->paramNystrom->value() ); // Nystrom sampling rate
		options << QString("--crf_sigma %1").arg( sigma );
		options << QString("--crf_weight %1").arg( 20 );
	}

	options << QString("-p %1").arg( path ); // Path
	options << QString("-f %1").arg( f + filename.right(4) ); // File
	options << QString("-C %1").arg( filenames.back() ); // Foreground
	options << QString("-C %1").arg( filenames.front() ); // Background
	options << QString("-C %1").arg( path + "output" ); // Output
	qDebug() << "Running : " << "RobustSelect " + options.join(" ").toLatin1();

	/// Call using system():
	//options.push_front("RobustSelect");
	//system( options.join(" ").toLatin1() );
	
	// Using process
	QProcess robustselect;
	QString params = options.join(" ").toLatin1();
	robustselect.start("RobustSelect", params.split(" "));
	robustselect.waitForFinished();
	QByteArray result = robustselect.readAll();

	// Get the output image
	QDir d(path);
	QStringList outputs = d.entryList().filter(QRegExp(".*transparent.png$"));

	if( outputs.size() ){
		QString outputFilename = outputs.front();
		QImage segmented(path + outputFilename);

		img = noBackgroundImage(segmented);
	}

	fg.clear();
	bg.clear();
}

void MyImageArea::cleanupRemove()
{
	QString path = QFileInfo( QCoreApplication::applicationFilePath() ).absolutePath() + "/";
	QString f = QFileInfo(filename).baseName();

	QDir d(path);

	// Remove image and stroke files
	foreach(QString file, d.entryList().filter(QRegExp( QString("^%1").arg(f) )))
		d.remove(file);
	
	// Remove output
	foreach(QString file, d.entryList().filter(QRegExp( QString("^output.sig") )))
		d.remove(file);
}

void MyImageArea::removeBackground( double sigma )
{
	img = centerWithWhiteBackground(orig_img);

	saveBackRemove();
	runBackRemove(sigma);
	cleanupRemove();
}

void MyImageArea::strokesFromThreshold(double t)
{
	int w = img.width();
	int h = img.height();

	// Clear previous strokes
	bg.clear(); fg.clear();
	bg.resize(h, std::vector<bool>(w, false));
	fg.resize(h, std::vector<bool>(w, false));

	// Set foreground
	for(int y = 0; y < h; y++){
		for(int x = 0; x < w; x++){
			QColor c(img.pixel(x,y));
			std::vector<int> ch = channels(c);
			double v = (double)ch[CH_LIGHT] / 255;
			if(v < t) fg[y][x] = true;
		}
	}

	// Set default background strokes
	int k = 10;
	w -= k;
	h -= k;

	DrawLine(k, k, w, k, bg, 5, true);
	DrawLine(k, k, k, h, bg, 5, true);
	DrawLine(w, k, w, h, bg, 5, true);
	DrawLine(k, h, w, h, bg, 5, true);
}

void MyImageArea::autoRemoveBack1()
{
	strokesFromThreshold( ui->thresholdVal->value() );
	removeBackground( ui->paramSigma->value() );
}

void MyImageArea::autoRemoveBack2()
{
	QList<double> sigmas;
	sigmas << 0.2 << 0.1 << 4;

	for(auto sigma : sigmas){
		strokesFromThreshold( ui->thresholdVal->value() );
		removeBackground( sigma );

		// Test quality of output
		QColor p0 = img.pixel(0,0);
		if( p0 != qRgb(255,0,255) ) 
			break;
	}
}
