#include "ShapeCorresponder.h"
#include "DeformPathItem.h"
#include <QGraphicsScene>
#include <QGraphicsView>
#include "SynthesisManager.h"
#include "ProjectedStructureGraph.h"

static inline QString shortName(QString name){
	if(name.length() < 3) return name;
	return QString("%1%2%3").arg(name.at(0)).arg(name.at(1)).arg(name.at(name.length()-1));
}

DeformPathItem::DeformPathItem(int width, int height, DeformationPath *usedPath) : path(usedPath){
	m_rect = QRectF(0,0,width,height);
}

QRectF DeformPathItem::boundingRect() const{
	return m_rect;
}

void qgluPerspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
	const GLdouble ymax = zNear * tan(fovy * M_PI / 360.0);
	const GLdouble ymin = -ymax;
	const GLdouble xmin = ymin * aspect;
	const GLdouble xmax = ymax * aspect;
	glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
}

QRectF setCamera(QRectF boundingRect, Structure::Graph * g)
{
	qglviewer::Camera * camera = new qglviewer::Camera;

	camera->setUpVector(qglviewer::Vec(0,0,1));
	camera->setPosition(qglviewer::Vec(-2,-2,0.8));
	camera->lookAt(qglviewer::Vec());
	camera->setSceneRadius( 10 );
	camera->showEntireScene();

	if(camera->type() != qglviewer::Camera::PERSPECTIVE) camera->setType(qglviewer::Camera::PERSPECTIVE);

	QRectF r = boundingRect;

	qglviewer::Vec viewDir = camera->viewDirection();
	Eigen::AlignedBox3d graphBBox = g->bbox();
	double distance = graphBBox.sizes().maxCoeff() * 2.25;
	Vector3 center = graphBBox.center();
	Vector3 newPos = center - (distance * Vector3(viewDir[0], viewDir[1], viewDir[2]));
	camera->setRevolveAroundPoint( qglviewer::Vec(center) );
	qglviewer::Vec new_pos(newPos);
	camera->frame()->setPositionWithConstraint(new_pos);
	camera->setScreenWidthAndHeight(r.width(), r.height());
	camera->loadProjectionMatrix();
	camera->loadModelViewMatrix();

	return r;
}

void setLights()
{
	// Setup lights and material
	GLfloat ambientLightColor[] = {0.06f,0.06f,0.06f,1};
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLightColor);

	GLfloat diffuseLightColor[] = {0.8f,0.8f,0.8f,1};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLightColor);

	GLfloat specularLightColor[] = {0.95f,0.95f,0.95f,1};
	glLightfv(GL_LIGHT0, GL_SPECULAR, specularLightColor);

	float posLight0[] = { 3, -3, 3, 0 };
	glLightfv(GL_LIGHT0, GL_POSITION, posLight0);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);

	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

	// Specular lighting
	float specReflection[] = { 0.8f, 0.8f, 0.8f, 1.0f };
	glMaterialfv(GL_FRONT, GL_SPECULAR, specReflection);
	glMateriali(GL_FRONT, GL_SHININESS, 40);
}

void DeformPathItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *, QWidget *)
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	QGraphicsView * view = scene()->views().front();
	QPointF sceneP = mapToScene( QPointF(0,0) );
	QPoint viewP = view->mapFromScene(sceneP);

	double dy = -viewP.y();
	double dx = viewP.x();

	if( dy > viewport[3] || dy < -viewport[3] || dx > viewport[2] || dx < -viewport[2] ) return;

	QRectF m_rect_scaled = view->transform().mapRect(m_rect);
	int width = m_rect_scaled.width();
	int height = m_rect_scaled.height();
	int smallWidth = height * 0.5;
	int tinyWidth = height * (1.0 / 5.0);

	QRectF srect(width - smallWidth ,0,smallWidth,smallWidth);
	QRectF trect(width - smallWidth,smallWidth,smallWidth,smallWidth);
	QRectF inbetween(0, tinyWidth, width - smallWidth, height - tinyWidth);

	// DEBUG:
	if( false ){
		painter->drawText(QPoint(10,200), QString("IDX %3  = %1, %2, dy = %4").arg(viewP.x()).arg(viewP.y()).arg(path->idx).arg(dy));
		painter->setPen( QPen(Qt::green, 2) );	painter->drawRect( srect );
		painter->setPen( QPen(Qt::blue, 2) );	painter->drawRect( trect );
		painter->setPen( QPen(Qt::black, 2) );	painter->drawRect( inbetween );
	}

	// Draw projected in between if any
	{
		bool isTopoBlend = !path->scheduler.isNull() && path->scheduler->allGraphs.size() && path->property["isReady"].toBool();
		if( !isTopoBlend && path->projected.size() )
		{
			double t = double(path->si) / 99;
			painter->drawImage(inbetween, path->projected.front()->drawBlendedImage(path->projected.back(), path->gcorr, t));
		}
	}

	// Draw score
	painter->drawText( QPoint(10,inbetween.top() + 20), QString("%1").arg(path->weight) );

	// Draw assignments
	painter->save();
	{
		QFont font("Monospace", 7); font.setStyleHint(QFont::TypeWriter);
		painter->setFont(font);
		//painter->drawText( QPoint(10,inbetween.top() + 30), "TestingAssignment");

		QVector< QPair<double, QString> > ssorted, tsorted;
		for(auto n : path->gcorr->sg->nodes) ssorted.push_back( qMakePair(n->bbox().center().z(), n->id) ); qSort(ssorted);
		for(auto n : path->gcorr->tg->nodes) tsorted.push_back( qMakePair(n->bbox().center().z(), n->id) ); qSort(tsorted);

		QMap<QString,int> sorder, torder;
		QRectF titleRect(0,0,30,10);
		QTextOption aligned(Qt::AlignAbsolute|Qt::AlignHCenter|Qt::AlignVCenter);
		int padding = 3; // px;

		painter->setPen(QPen(Qt::black, 1));

		for(auto n : ssorted){
			int row = sorder.size();
			sorder[n.second] = row;
			titleRect.moveTop( inbetween.top() + (row * (titleRect.height() + padding)) );
			titleRect.moveRight( srect.left() - (titleRect.width() * 2) );
			if(path->scolors.contains(n.second)) painter->setBrush(QBrush(path->scolors[n.second])); else painter->setBrush(Qt::NoBrush);
			painter->drawRoundedRect(titleRect, 3, 3);
			painter->drawText( titleRect, shortName(n.second), aligned );
		}

		for(auto n : tsorted){
			int row = torder.size();
			torder[n.second] = row;
			titleRect.moveTop( inbetween.top() + (row * (titleRect.height() + padding)) );
			titleRect.moveRight( srect.left() );
			if(path->tcolors.contains(n.second)) painter->setBrush(QBrush(path->tcolors[n.second])); else painter->setBrush(Qt::NoBrush);
			painter->drawRoundedRect(titleRect, 3, 3);
			painter->drawText( titleRect, shortName(n.second), aligned );
		}

		// Draw assignments
		QRectF stitleRect = titleRect, ttitleRect = titleRect;
		for(auto pair : path->pairs)
		{
			if(pair.first.contains("NOTHING") || pair.second.contains("NOTHING")) continue; // don't show to 'null'

			for(auto s : pair.first){
				int srow = sorder.contains(s) ? sorder[s] : sorder.size();

				for(auto t : pair.second){
					int trow = torder.contains(t) ? torder[t] : torder.size();

					stitleRect.moveTop( inbetween.top() + (srow * (titleRect.height() + padding)) );
					stitleRect.moveRight( srect.left() - (titleRect.width() * 2) );

					ttitleRect.moveTop( inbetween.top() + (trow * (titleRect.height() + padding)) );
					ttitleRect.moveRight( srect.left() );

					int hheight = stitleRect.height() * 0.5;
					painter->drawLine(stitleRect.topRight() + QPointF(0, hheight), ttitleRect.topLeft() + QPointF(0, hheight));
				}
			}
		}
	}
	painter->restore();

	// Draw border
	painter->setPen(QPen(Qt::gray, 1));

	// Special border for ground truth
	if(path->property.contains("isGroundTruth") && path->property.value("isGroundTruth").toBool())
		painter->setPen( QPen(Qt::green, 8) );

	if(this->parentItem()->isSelected()) painter->setPen(QPen(Qt::blue, 5));
	painter->drawRect(boundingRect());

	// Draw 3D parts
	painter->beginNativePainting();

	setCamera(srect, path->gcorr->sg);
	setLights();
	glEnable(GL_DEPTH_TEST);
	glClear(GL_DEPTH_BUFFER_BIT);

	bool isDrawMeshes = true;

	// Draw source
	{
		Structure::Graph * g = path->gcorr->sg;

		glViewport(dx + srect.x(), dy + viewport[3] - srect.height() - srect.top(), srect.width(), srect.height());

		for(auto n : g->nodes)
		{
			QColor c = QColor::fromRgbF(0.8,0.8,0.8);
			glColor3d(c.redF(), c.greenF(), c.blueF());
			if(path->scolors.contains(n->id)) glColorQt(path->scolors[n->id]);
			glDisable(GL_LIGHTING);
			n->draw( false, true );

			if( isDrawMeshes )
			{
				glEnable(GL_LIGHTING);
				if(path->scolors.contains(n->id)) c = path->scolors[n->id];
				g->drawNodeMesh(n->id, c);
			}
		}
	}

	// Draw target
	{
		Structure::Graph * g = path->gcorr->tg;

		glViewport(dx + trect.x(), dy + viewport[3] - trect.height() - trect.top(), trect.width(), trect.height());

		for(auto n : g->nodes)
		{
			QColor c = QColor::fromRgbF(0.8,0.8,0.8);
			glColor3d(c.redF(), c.greenF(), c.blueF());
			if(path->tcolors.contains(n->id)) glColorQt(path->tcolors[n->id]);
			glDisable(GL_LIGHTING);
			n->draw( false, true );

			if( isDrawMeshes )
			{
				glEnable(GL_LIGHTING);
				if(path->tcolors.contains(n->id)) c = path->tcolors[n->id];
				g->drawNodeMesh(n->id, c);
			}
		}
	}

	// Draw in between
	if( !path->scheduler.isNull() && path->scheduler->allGraphs.size() && path->property["isReady"].toBool() )
	{
		Structure::Graph * g = path->scheduler->allGraphs[path->si];

		glViewport(dx + inbetween.x(), dy + viewport[3] - inbetween.height() - inbetween.top(), inbetween.width(), inbetween.height());

		for(auto n : g->nodes)
		{
			if( Scheduler::shouldNotExist(n) ) continue;

			glColor3d(0.8,0.8,0.8);

			QString sid = n->property["original_ID"].toString();
			QString tid = n->property["correspond"].toString();

			if(path->scolors.contains(sid))
				glColorQt(path->scolors[sid]);
			else if(path->tcolors.contains(tid))
				glColorQt(path->tcolors[tid]);

			glDisable(GL_LIGHTING);
			n->draw( false, true );
		}

		// Draw with surface
		if(path->property.contains("synthManager"))
		{
			SynthesisManager* sm = path->property["synthManager"].value<SynthesisManager*>();
			if( sm )
			{
				glEnable(GL_LIGHTING);
				sm->color = Qt::red;

				if(sm->proxies.size())
				{
					for(auto poly : sm->drawWithProxies( g ))
					{
						glColorQt( sm->color );
						glBegin(GL_POLYGON);
						Vector3 normal = (poly.vertices[1] - poly.vertices[0]).cross(poly.vertices[2] - poly.vertices[0]).normalized();
						glNormal3( normal );
						for(auto v : poly.vertices) glVector3(v);
						glEnd();
					}
				}
				else
					sm->drawSynthesis( g );
			}
		}
	}

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	glViewport(viewport[0],viewport[1],viewport[2],viewport[3]);
	painter->endNativePainting();

	// Draw graph
	{
		QRectF graphRect = inbetween;
		graphRect.setHeight( inbetween.height() * 0.10 );

		QPainterPath graph;

		painter->save();
		painter->translate( graphRect.topLeft() );

		double sum = 0;

		for(int i = 0; i < path->errors.size(); i++)
		{
			double dx = (double(i) / (path->errors.size()-1));
			double dy = (sum += path->errors[i]) / path->weight;

			double x = dx * graphRect.width();
			double y = graphRect.height() - (dy * graphRect.height());

			if( i == 0 ) graph.moveTo(x,y);
			else graph.lineTo(x,y);
		}

		painter->setPen( QPen(QColor(0,128,0), 2) );
		painter->drawPath(graph);

		painter->restore();
	}
}