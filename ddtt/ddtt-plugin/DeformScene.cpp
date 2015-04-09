#include "DeformScene.h"
#include "ShapeCorresponder.h"
#include <QGraphicsView>
#include <QHBoxLayout>
#include <QDesktopWidget>
#include <QApplication>
#include <QSlider>
#include <QLabel>

#include "DeformPathItem.h"
#include "DeformPathItemWidget.h"

DeformScene::DeformScene()
{
	gview = new QGraphicsView(this);
    gview->setViewport(new QGLWidget);
    gview->setViewportUpdateMode(QGraphicsView::FullViewportUpdate);
	gview->setAlignment(Qt::AlignLeft | Qt::AlignTop);
	gview->setRenderHint(QPainter::HighQualityAntialiasing, true);
	gview->setRenderHint(QPainter::SmoothPixmapTransform, true);
	gview->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
	gview->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    gview->setDragMode(QGraphicsView::RubberBandDrag);

	connect(this, SIGNAL(viewPortReady()), SLOT(init()));
	emit( viewPortReady() );
}

void DeformScene::init()
{
	gview->resize(1580,820);

	// Center to screen
	QDesktopWidget* m = QApplication::desktop();
	QRect desk_rect = m->screenGeometry(m->screenNumber(QCursor::pos()));
	int desk_x = desk_rect.width();
	int desk_y = desk_rect.height();
	int x = gview->width();
	int y = gview->height();

	gview->move(desk_x / 2 - x / 2 + desk_rect.left(), desk_y / 2 - y / 2 + desk_rect.top());
	gview->show();
	gview->raise();
}

void DeformScene::drawBackground(QPainter *painter, const QRectF &rect)
{
    QGraphicsScene::drawBackground(painter,rect);
	painter->drawRect(this->sceneRect());
}

void DeformScene::drawForeground(QPainter *painter, const QRectF &rect)
{
    QGraphicsScene::drawForeground(painter, rect);
}

void DeformScene::addDeformationPath(DeformationPath * path)
{
	int w = 550;
	int h = w * 0.75;

	// Create a DeformationItem and its widget
	DeformPathItem * ditem = new DeformPathItem(w,h, path);
	DeformPathItemWidget * dwidget = new DeformPathItemWidget(w,h, path);
	
	QGraphicsItemGroup * group = new QGraphicsItemGroup;

	group->addToGroup(ditem);
	group->addToGroup(dwidget);
	group->setPos(0, path->idx * h);
	group->setAcceptHoverEvents(true);
	group->setHandlesChildEvents(false); // items in group get user inputs
	//group->setFlags( QGraphicsItem::ItemIsMovable | QGraphicsItem::ItemIsSelectable );

	addItem(group); 

	groups.push_back(group);
}

void DeformScene::wheelEvent( QGraphicsSceneWheelEvent *e )
{
    // zoom only when CTRL key pressed
    if (e->modifiers().testFlag(Qt::ControlModifier))
    {
        int numSteps = e->delta() / 15 / 8;
        if (numSteps == 0) {
            e->ignore();
            return;
        }

        qreal sc = pow(1.25, numSteps); // I use scale factor 1.25
        QGraphicsView * view = views().front();
        view->scale(sc, sc);
        view->centerOn(view->mapToScene(e->pos().toPoint()));
        e->accept();
    }
}

void DeformScene::pack()
{
	int itemWidth = groups.front()->boundingRect().width();
	int itemHeight = groups.front()->boundingRect().height();
	int viewWidth = this->views().front()->width();

	int xi = 0, row = 0;

	for(auto item : groups)
	{
		if((xi * itemWidth) > viewWidth) {xi = 0; row++;}
		
		int x = xi++ * itemWidth;
		int y = row * itemHeight;

		item->setPos(x,y);
	}
}
