#pragma once

#include <QGraphicsObject>

class DeformationPath;
typedef QMap<QString,QVariant> PropertyMap;

class DeformPathItem : public QGraphicsObject
{
	Q_OBJECT
public:
    DeformPathItem(int width, int height, DeformationPath * usedPath);
	PropertyMap property;

	DeformationPath * path;
	QRectF m_rect;

public:
	QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
};
