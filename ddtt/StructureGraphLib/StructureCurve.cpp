#include "StructureCurve.h"
#include "LineSegment.h"
using namespace Structure;

#include "qglviewer/quaternion.h"

#if defined(Q_OS_MAC)
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include <Eigen/Geometry>

Curve::Curve(const NURBS::NURBSCurved & newCurve, QString newID, QColor color)
{
    this->curve = newCurve;
    this->id = newID;
    this->vis_property["color"] = color;
    this->vis_property["showControl"] = false;
}

Node * Curve::clone()
{
	if(!this) return NULL;

	Curve * cloneCurve = new Curve( this->curve, this->id );
	cloneCurve->property = this->property;
	cloneCurve->vis_property = this->vis_property;
	return cloneCurve;
}

QString Curve::type()
{
    return CURVE;
}

Eigen::AlignedBox3d Curve::bbox(double scaling)
{
    Eigen::AlignedBox3d box;

    foreach(Vector3d cp, curve.getControlPoints())
        box = box.merged( Eigen::AlignedBox3d(cp, cp) );

	// Scaling
	Eigen::Vector3d diagonal = box.diagonal() * 0.5;

	Eigen::Vector3d a = box.center() + (diagonal * scaling);
	Eigen::Vector3d b = box.center() - (diagonal * scaling);

	box = box.merged( Eigen::AlignedBox3d(a,a) );
	box = box.merged( Eigen::AlignedBox3d(b,b) );

    return box;
}

std::vector<int> Curve::controlCount()
{
	return std::vector<int>( 1, curve.GetNumCtrlPoints() );
}

std::vector<Vector3> Curve::controlPoints()
{
    return curve.mCtrlPoint;
}

std::vector<Scalar> Curve::controlWeights()
{
    return curve.mCtrlWeight;
}

void Curve::get( const Vector4d& coordinates, Vector3 & pos, std::vector<Vector3> & frame )
{
	double u = coordinates[0];
	Vector3 der1(0,0,0);

	frame.resize(3, Vector3(0,0,0));

	curve.GetFrame(u, pos, frame[0], frame[1], frame[2]);
}

SurfaceMesh::Vector3 Curve::position( const Vector4d& coordinates )
{
    std::vector<Vector3d> nf = noFrame();
    Vector3 p(0,0,0); 
	get(coordinates, p, nf);
	return p;
}

Vector4d Curve::approxCoordinates( const Vector3 & pos )
{
	Scalar t = curve.timeAt( pos );
	return Vector4d( t, 0, 0, 0 );
}

SurfaceMesh::Vector3 Curve::approxProjection( const Vector3 & point )
{
	Vector4d coords = approxCoordinates(point);
	return curve.GetPosition(coords[0]);
}

Array1D_Vector3 Curve::discretizedAsCurve(Scalar resolution)
{
    Array1D_Vector3 result;
    Scalar curveLength = curve.GetLength(0,1);

    // Degenerate cases
    bool validControlPoints = std::isfinite(curve.mCtrlPoint.front().x());
    if( !validControlPoints ){
        result.clear();
        for(int i = 0; i < 4; i++) result.push_back(Vector3(0,0,0));
        return result;
    }
    if( !std::isfinite(curveLength) || curveLength < resolution )
        return curve.mCtrlPoint;

    int np = qMax(2.0, 1.0 + (curveLength / resolution));
    curve.SubdivideByLength(np, result);

    return result;
}

Array2D_Vector3 Curve::discretized(Scalar resolution)
{
	return curve.toSegments( resolution );
}

Array2D_Vector4d Curve::discretizedPoints( Scalar resolution )
{
	Array2D_Vector4d result;

	Scalar curveLength = curve.GetLength(0,1);

	double firstTwoDist = (curve.mCtrlPoint[0] - curve.mCtrlPoint[1]).norm();
	if(firstTwoDist < resolution * 0.001){
		result.push_back( Array1D_Vector4d( 1, Vector4d(0,0,0,0) ) );
		return result;
	}

	// For singular cases
	if(curveLength < resolution){
		result.push_back( Array1D_Vector4d( 1, Vector4d(0,0,0,0) ) );
		return result;
	}

	if(!IsNumber(curveLength) || !IsFiniteNumber(curveLength)){
		result.push_back( Array1D_Vector4d( 1, Vector4d(0,0,0,0) ) );
		return result;
	}

    int np = qMax(2.0, 1.0 + (curveLength / resolution));
	std::vector<Scalar> ptsTimes;

	curve.SubdivideByLengthTime(np, ptsTimes);

	Array1D_Vector4d resultTimes;
	result.push_back( resultTimes );

	foreach(Scalar t, ptsTimes)
		result.back().push_back(Vector4d(t,0,0,0));

	return result;
}

void Curve::laplacianSmoothControls( int num_iterations, std::set<int> anchored )
{
	std::vector<Vector3> & cpnts = curve.mCtrlPoint;

	// Special case anchoring
	if(anchored.count(-1)){
		anchored.clear();		
		anchored.insert(0);
		anchored.insert(cpnts.size() - 1);
	}

	// Laplacian smoothing
	for(int itr = 0; itr < num_iterations; itr++)
	{
		std::vector<Vector3> newCtrlPnts = cpnts;

		for (int j = 0; j < (int)cpnts.size(); j++)
		{
			if(anchored.count(j) == 0)
				newCtrlPnts[j] = (cpnts[j-1] + cpnts[j+1]) / 2.0;
		}

		for (int j = 0; j < (int)cpnts.size(); j++)
			cpnts[j] = newCtrlPnts[j];
	}
}

void Curve::moveBy( const Vector3d & delta )
{
	curve.translate( delta );
}

std::vector<Vector3d> Curve::foldTo( Vector4d & foldPoint, bool isApply)
{
	int cpIDX = controlPointIndexFromCoord(foldPoint);
	Vector3 cp = curve.mCtrlPoint[cpIDX];

	std::vector<Vector3d> deltas;

	for(size_t i = 0; i < curve.mNumCtrlPoints; i++)
	{
		deltas.push_back(curve.mCtrlPoint[i] - cp);
		if( isApply ) curve.mCtrlPoint[i] -= deltas.back();
	}

	return deltas;
}

void Curve::setControlPoints( const std::vector<Vector3> & newPositions )
{
	assert(curve.mCtrlPoint.size() == newPositions.size());

	curve.mCtrlPoint = newPositions;
}

void Curve::scale( Scalar scaleFactor )
{
	this->curve.scale(scaleFactor);
}

void Curve::rotate( double angle, Vector3 axis )
{
	angle *= 3.14159265358979 /180; 

	std::vector<Vector3> &mCtrlPoint = this->curve.mCtrlPoint;

	for(int i = 0; i < (int)mCtrlPoint.size(); i++)
		mCtrlPoint[i] = rotatedVec(mCtrlPoint[i], angle, axis);
}


Vector3 & Curve::controlPoint( int idx )
{
	return curve.mCtrlPoint[idx];
}

int Curve::controlPointIndexFromCoord( Vector4d& coord )
{
	return (curve.mCtrlPoint.size() - 1) * coord[0];
}

Vector3 & Curve::controlPointFromCoord( Vector4d& coord )
{
	return curve.mCtrlPoint[controlPointIndexFromCoord( coord )];
}

SurfaceMesh::Scalar Curve::area()
{
	return curve.GetTotalLength();
}

SurfaceMesh::Vector3 Curve::center()
{
	Vector3 pos(0,0,0);
    std::vector<Vector3d> nf = noFrame();
    get(Vector4d(0.5,0.5,0,0), pos, nf);
	return pos;
}

void Curve::draw(bool isShowCtrlPts, bool isColorOverride)
{
	if(vis_property["glow"].toBool())
	{
        NURBS::CurveDraw::draw( &curve, Qt::yellow, isShowCtrlPts, isColorOverride, 2.5 );
	}
	else
	{
        NURBS::CurveDraw::draw( &curve, vis_property["color"].value<QColor>(), isShowCtrlPts, isColorOverride );
	}

    if(isColorOverride) return;

	glPointSize(10);
	glColor3d(1,0,1);
	glBegin(GL_POINTS);
	foreach(Vector3 p, curve.misc_points) glVector3(p);
	glEnd();

	// Draw selections
	GLUquadricObj *quadObj = gluNewQuadric();

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);

	foreach (int pID, selections.keys())
	{
		QColor color = selections[pID];
		glColor3d(color.red(), color.green(), color.blue());

		Vector3 p = curve.GetControlPoint(pID);

		glPushMatrix();
		glTranslatef(p.x(), p.y(), p.z());
		gluSphere(quadObj, 0.02, 16, 16);
		glPopMatrix();
	}

	gluDeleteQuadric(quadObj);
}

void Curve::drawWithNames( int nID, int pointIDRange )
{
	int pID = nID * pointIDRange;

	glPointSize(20.0f);
	for(size_t i = 0; i < curve.mNumCtrlPoints; i++)
	{
		glPushName(pID++);

		Vector3d p = curve.GetControlPoint(i);

		glBegin(GL_POINTS);
		glVertex3d(p.x(), p.y(), p.z());
		glEnd();

		glPopName();
	}
}

void Curve::equalizeControlPoints( Node * _other )
{
	Curve *other = (Curve *)_other;

	int targetNumber = other->curve.mNumCtrlPoints;
	int diff = targetNumber - curve.mNumCtrlPoints;
	if(diff == 0.0) return;

	std::vector<Vector3d> newPnts = this->curve.simpleRefine( diff );

    this->curve = NURBS::NURBSCurved(newPnts, std::vector<double>( newPnts.size(), 1.0 ));
}

int Curve::numCtrlPnts()
{
	return curve.mNumCtrlPoints;
}

void Curve::refineControlPoints( int nU, int nV /*= 0*/ )
{
    nV = nV;

	int diff = nU - curve.mNumCtrlPoints;
	if(diff == 0) return;

	std::vector<Vector3d> newPnts = this->curve.simpleRefine( diff );

	this->curve = NURBS::NURBSCurved(newPnts, std::vector<double>( newPnts.size(), 1.0 ));
}

Vector3d Curve::direction()
{
	Vector3d dir = curve.mCtrlPoint.back() - curve.mCtrlPoint.front();
	return dir.normalized();
}



/* Curve encoding, to decode you need two points A,B and a frame XYZ */
CurveEncoding Curve::encodeCurve( Array1D_Vector3 points, Vector3 start, Vector3 end, bool isFlip )
{
	CurveEncoding cpCoords;

	NURBS::Line segment(start, end);

	// compute frame
	qglviewer::Vec x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
	Vector3d Z = (end - start).normalized();
	qglviewer::Quaternion q(z, qglviewer::Vec(Z));
	x = q * x; y = q * y;
	Vector3d X(x[0], x[1], x[2]);
	Vector3d Y(y[0], y[1], y[2]);

	for(int i = 0; i < (int)points.size(); i++)
	{
		double t;
		Vector3 p = points[i], proj;
		segment.ClosestPoint(p, t, proj);

		Vector3 dir = p - proj;

		// Parameters: t, offset, theta, psi
		Array1D_Real params(4,0);
		params[0] = t;
		if(dir.norm() > 0){
			params[1] = dir.norm() / segment.length;
			dir = dir.normalized();
		}
		globalToLocalSpherical(X,Y,Z, params[2], params[3], dir);

		// Flipping case
		int idx = i;
		if(isFlip) idx = (points.size()-1) - i; 

		cpCoords[idx] = params;
	}

	return cpCoords;
}

CurveEncoding Curve::encodeCurve( Curve * curve, Vector3 start, Vector3 end, bool isFlip )
{
	return encodeCurve(curve->controlPoints(),start,end,isFlip);
}

Array1D_Vector3 Curve::decodeCurve(CurveEncoding cpCoords, Vector3 start, Vector3 end, double T)
{
	Array1D_Vector3 controlPoints (cpCoords.size(), Vector3(0,0,0));

	NURBS::Line segment(start, end);

	if(segment.length < DECODE_ZERO_THRESHOLD){
		return Array1D_Vector3( cpCoords.size(), start );
	}

	// compute frame
	qglviewer::Vec x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
	Vector3d Z = (end - start).normalized();
	qglviewer::Quaternion q(z, qglviewer::Vec(Z));
	x = q * x; y = q * y;
	Vector3d X(x[0], x[1], x[2]);
	Vector3d Y(y[0], y[1], y[2]);

	for(int i = 0; i < (int)controlPoints.size(); i++)
	{
		double t = cpCoords[i][0];
		double offset = cpCoords[i][1];
		double theta = cpCoords[i][2];
		double psi = cpCoords[i][3];

		Vector3 dir(0,0,0);
		localSphericalToGlobal(X,Y,Z,theta,psi,dir);

		controlPoints[i] = segment.pointAt(t) + (dir * ((offset * T) * segment.length));
	}

	return controlPoints;
}


void Curve::deformTo( const Vector4d & handle, const Vector3 & to, bool isRigid )
{
	Vector4d otherEndCoord = Vector4d((handle[0] > 0.5) ? 0 : 1);

	Vector3 p = position( handle );

	double diff = (p-to).norm();
	if(diff < DECODE_ZERO_THRESHOLD) return;

	if( isRigid )
	{
		Vector3d delta = to - p;
		this->moveBy( delta );
		return;
	}

	Eigen::Vector3d otherEnd = position( otherEndCoord );

	// Find new scale
	double d = (p - otherEnd).norm();
	if(d < DECODE_ZERO_THRESHOLD) return;

	double scale = (to - otherEnd).norm() / d;

	// Find minimum rotation
	Eigen::Vector3d dirFrom (p - otherEnd);
	Eigen::Vector3d dirTo (to - otherEnd);

	dirFrom.normalize();
	dirTo.normalize();

	Eigen::Quaterniond rotation = Eigen::Quaterniond::FromTwoVectors(dirFrom, dirTo).normalized();

	Array1D_Vector3 ctrlPnts = controlPoints();

	// Compute deltas and apply transformations
	for(int i = 0; i < (int)ctrlPnts.size(); i++)
	{
		Eigen::Vector3d delta = (ctrlPnts[i] - otherEnd);
		double length = delta.norm() * scale;
		delta.normalize();

		Eigen::Vector3d rotated = rotation * delta;

		ctrlPnts[i] = Vector3( (rotated * length) + otherEnd );
	}

	setControlPoints( ctrlPnts );
}

void Curve::deformTwoHandles( Vector4d& handleA, Vector3 newPosA, Vector4d& handleB, Vector3 newPosB )
{
	Vector3d oldA = position(handleA);
	Vector3d oldB = position(handleB);

	// Numerical checks
	double oldDiff = (oldA-oldB).norm(); if(oldDiff < DECODE_ZERO_THRESHOLD) return;
	double diff = (newPosA-newPosB).norm();	if(diff < DECODE_ZERO_THRESHOLD) return;

	setControlPoints( decodeCurve( encodeCurve(this, oldA, oldB), newPosA, newPosB ) );
}
