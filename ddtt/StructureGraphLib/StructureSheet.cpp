// Used for sheet folding
#include "NanoKdTree.h"
#include "GraphDistance.h"
#include "PCA.h"
using namespace qglviewer;

#if defined(Q_OS_MAC)
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif

#include "StructureSheet.h"
using namespace Structure;

Sheet::Sheet(const NURBS::NURBSRectangled & sheet, QString sheetID, QColor color)
{
    this->surface = sheet;
	this->surface.quads.clear();

    this->id = sheetID;
    this->vis_property["color"] = color;
    this->vis_property["showControl"] = false;
}

Node * Sheet::clone()
{
	if(!this) return NULL;

	Sheet * cloneSheet = new Sheet( this->surface, this->id );
	cloneSheet->surface.quads.clear();
	cloneSheet->property = this->property;
	cloneSheet->vis_property = this->vis_property;
	return cloneSheet;
}

QString Sheet::type()
{
    return SHEET;
}

Eigen::AlignedBox3d Sheet::bbox(double scaling)
{
    Eigen::AlignedBox3d box;

    foreach(std::vector<Vector3d> cps, surface.mCtrlPoint)
        foreach(Vector3d cp, cps)
            box = box.merged( Eigen::AlignedBox3d(cp, cp) );

	// Scaling
	Eigen::Vector3d diagonal = box.diagonal() * 0.5;
	Eigen::Vector3d a = box.center() + (diagonal * scaling);
	Eigen::Vector3d b = box.center() - (diagonal * scaling);

	box = box.merged( Eigen::AlignedBox3d(a,a) );
	box = box.merged( Eigen::AlignedBox3d(b,b) );

    return box;
}

std::vector<int> Sheet::controlCount()
{
	std::vector<int> count;

	count.push_back(surface.mNumUCtrlPoints);
	count.push_back(surface.mNumVCtrlPoints);

	return count;
}

std::vector<Vector3> Sheet::controlPoints()
{
	std::vector<Vector3> cpoints;

	for(size_t u = 0; u < surface.mNumUCtrlPoints; u++)
	{
		for(size_t v = 0; v < surface.mNumVCtrlPoints; v++)
		{
			cpoints.push_back( surface.mCtrlPoint[u][v] );
		}
	}

	return cpoints;
}

void Sheet::setControlPoints( const std::vector<Vector3> & newPositions )
{
	int i = 0;

	assert(newPositions.size() == surface.mCtrlPoint.size() * surface.mCtrlPoint.front().size());

	for(size_t u = 0; u < surface.mNumUCtrlPoints; u++)
		for(size_t v = 0; v < surface.mNumVCtrlPoints; v++)
			surface.mCtrlPoint[u][v] = newPositions[i++];
}

std::vector<Scalar> Sheet::controlWeights()
{
	std::vector<Scalar> cpoints;

	for(size_t v = 0; v < surface.mNumVCtrlPoints; v++)
		for(size_t u = 0; u < surface.mNumUCtrlPoints; u++)
			cpoints.push_back( surface.mCtrlWeight[u][v] );

	return cpoints;
}

void Sheet::get( const Vector4d& coordinates, Vector3 & pos, std::vector<Vector3> & frame )
{
	double u = coordinates[0];
	double v = coordinates[1];

	frame.resize(3, Vector3(0,0,0));

	surface.GetFrame(u, v, pos, frame[0], frame[1], frame[2]);
}

SurfaceMesh::Vector3 Sheet::position( const Vector4d& coordinates )
{
    std::vector<Vector3> nf = noFrame();
    Vector3 p(0,0,0); get(coordinates,p, nf);
	return p;
}

Vector4d Sheet::approxCoordinates( const Vector3 & pos )
{
	return surface.timeAt( pos );
}

SurfaceMesh::Vector3 Sheet::approxProjection( const Vector3 & point )
{
	Vector3 pos(0,0,0);
	Vector4d coords = approxCoordinates(point);
    surface.Get(coords[0], coords[1], &pos);
	return pos;
}

Array1D_Vector3 Sheet::discretizedAsCurve(Scalar resolution)
{
    // Return a representative curve along the longest side at the middle of the sheet
    NURBS::NURBSCurved curve;

    Vector3 uDir = surface.mCtrlPoint.front().back() - surface.mCtrlPoint.front().front();
    Vector3 vDir = surface.mCtrlPoint.back().front() - surface.mCtrlPoint.front().front();

    double projux = abs(dot(Vector3d(1,0,0), uDir.normalized()));
    double projvx = abs(dot(Vector3d(1,0,0), vDir.normalized()));

    //bool isU = uDir.norm() > vDir.norm();
	bool isU = projux > projvx;

    if( isU )
    {
        int v = (surface.mNumUCtrlPoints - 1) * 0.5;
        curve = NURBS::NURBSCurved(surface.GetControlPointsV(v), surface.GetControlWeightsV(v));
    }
    else
    {
        int u = (surface.mNumVCtrlPoints - 1) * 0.5;
        curve = NURBS::NURBSCurved(surface.GetControlPointsU(u), surface.GetControlWeightsU(u));
    }

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

Array2D_Vector3 Sheet::discretized(Scalar resolution)
{
	return surface.generateSurfaceTris( resolution );
}

Array2D_Vector4d Sheet::discretizedPoints( Scalar resolution )
{
	Array2D_Vector4d coords;
	surface.generateSurfacePointsCoords(resolution, coords);
	return coords;
}

void Sheet::moveBy( const Vector3d & delta )
{
	surface.translate( delta );
}

Vector3 & Sheet::controlPoint( int idx )
{
	int u = idx / surface.mNumVCtrlPoints;
	int v = idx - (u * surface.mNumVCtrlPoints);
	return surface.mCtrlPoint[u][v];
}

SurfaceMesh::Scalar Sheet::area()
{
	double factor = 0.5; // distance of first control points

	Scalar r = factor * (surface.mCtrlPoint[0][0] - surface.mCtrlPoint[0][1]).norm();
	Array2D_Vector3 tris = discretized(r);

	double a = 0;

	foreach(std::vector<Vector3> v, tris)
	{
		double triArea = 0.5 * cross(Vector3(v[1] - v[0]), Vector3(v[2] - v[0])).norm();
		a += triArea;
	}

	return a;
}

SurfaceMesh::Vector3 Sheet::center()
{
	Vector3 pos(0,0,0);
    std::vector<Vector3> nf = noFrame();
    get(Vector4d(0.5,0.5,0,0), pos, nf);
	return pos;
}

void Sheet::draw(bool isShowCtrlPts, bool isColorOverride)
{
	if(vis_property["glow"].toBool())
	{
        NURBS::SurfaceDraw::draw( &surface, Qt::yellow, isShowCtrlPts, isColorOverride, 2.5, Qt::yellow );
	}
	else
	{
        NURBS::SurfaceDraw::draw( &surface, vis_property["color"].value<QColor>(), isShowCtrlPts, isColorOverride );
	}

    if(isColorOverride) return;

	// Draw selections
	GLUquadricObj *quadObj = gluNewQuadric();

	gluQuadricDrawStyle(quadObj, GLU_FILL);
	gluQuadricNormals(quadObj, GLU_SMOOTH);

	foreach (int pID, selections.keys())
	{
		QColor color = selections[pID];
		glColor3d(color.red(), color.green(), color.blue());

		int u = pID / surface.mNumVCtrlPoints;
		int v = pID % surface.mNumVCtrlPoints;

		Vector3 p = surface.GetControlPoint(u, v);

		glPushMatrix();
		glTranslatef(p.x(), p.y(), p.z());
		gluSphere(quadObj, 0.02, 16, 16);
		glPopMatrix();
	}

	gluDeleteQuadric(quadObj);
}

void Sheet::drawWithNames( int nID, int pointIDRange )
{
	int pID = nID * pointIDRange;

    //float radius = 1.0f;
	glPointSize(20.0f);

	for(size_t u = 0; u < surface.mNumUCtrlPoints; u++)
	{
		for(size_t v = 0; v < surface.mNumVCtrlPoints; v++)
		{
			glPushName(pID++);

			Vector3d p = surface.GetControlPoint(u, v);

			glBegin(GL_POINTS);
			glVertex3d(p.x(), p.y(), p.z());
			glEnd();

			glPopName();
		}
	}
}

std::vector< std::vector<Vector3d> > Sheet::foldTo( const Array1D_Vector4d & curve, bool isApply )
{
	std::vector<Vector3d> projections, coords;
	NanoKdTree tree;
	foreach(Vector4d coord, curve) 
	{
		coords.push_back(Vector3d(coord[0], coord[1], 0));
		projections.push_back( surface.projectOnControl( coord[0], coord[1] ) );
		tree.addPoint( projections.back() );
	}
	tree.build();

	std::vector< std::vector<Vector3d> > deltas = std::vector< std::vector<Vector3d> > (surface.mNumUCtrlPoints, 
		std::vector<Vector3d>(surface.mNumVCtrlPoints, Vector3(0,0,0)));

	// Project onto itself (PCA then project on closest main axis)
	Vector3d mp = PCA::mainAxis(coords);
	//int midIdx = curve.size() * 0.5;

	for(size_t u = 0; u < surface.mNumUCtrlPoints; u++)
	{
		for(size_t v = 0; v < surface.mNumVCtrlPoints; v++)
		{
			Vector4d currCoord(double(u) / (surface.mNumUCtrlPoints-1), double(v) / (surface.mNumVCtrlPoints-1),0,0);

			double minDist = DBL_MAX;
			int minIDX = -1;

			for(int i = 0; i < (int)curve.size(); i++)
			{
				Vector4d coord = curve[i];

				if(qMax(mp[0],mp[1]) == mp[0])
					currCoord[1] = coord[1];
				else
					currCoord[0] = coord[0];

				double dist = (coord - currCoord).norm();

				if(dist < minDist)
				{
					minIDX = i;
					minDist = dist;
				}
			}

			deltas[u][v] = surface.mCtrlPoint[u][v] - projections[minIDX];
		}
	}

	if(isApply)
	{
		for(size_t u = 0; u < surface.mNumUCtrlPoints; u++)
		{
			for(size_t v = 0; v < surface.mNumVCtrlPoints; v++)
			{
				surface.mCtrlPoint[u][v] += -deltas[u][v];
			}
		}
	}

	return deltas;
}

void Sheet::scale( Scalar scaleFactor )
{
	this->surface.scale(scaleFactor);

	this->surface.quads.clear();
}

void Sheet::rotate( double angle, Vector3 axis )
{
	angle *= 3.14159265358979 /180; 

	Array2D_Vector3 &mCtrlPoint = this->surface.mCtrlPoint;

	for(int y = 0; y < (int)mCtrlPoint.size(); y++)
		for(int x = 0; x < (int)mCtrlPoint[0].size(); x++)
			mCtrlPoint[y][x]  = rotatedVec(mCtrlPoint[y][x], angle, axis);

	this->surface.quads.clear();
}

void Sheet::equalizeControlPoints( Node * _other )
{
	Sheet *other = (Sheet*) _other;

	int nU = qMax(other->surface.mNumUCtrlPoints,surface.mNumUCtrlPoints);
	int nV = qMax(other->surface.mNumVCtrlPoints,surface.mNumVCtrlPoints);

	Array2D_Vector3 my_cp, other_cp;

	NURBS::NURBSRectangled mySurface = surface;
	NURBS::NURBSRectangled otherSurface = other->surface;

	int k = 0;

	// Refine U
	k = nU - mySurface.mNumUCtrlPoints;
	my_cp = mySurface.simpleRefine(k, 0);

	k = nU - otherSurface.mNumUCtrlPoints;
	other_cp = otherSurface.simpleRefine(k, 0);

	// Both now have same number of U control points
	mySurface = NURBS::NURBSRectangled::createSheetFromPoints(my_cp);
	otherSurface = NURBS::NURBSRectangled::createSheetFromPoints(other_cp);

	// Refine V
	k = nV - mySurface.mNumVCtrlPoints;
	my_cp = mySurface.simpleRefine(k, 1);

	k = nV - otherSurface.mNumVCtrlPoints;
	other_cp = otherSurface.simpleRefine(k, 1);

	// Both have exactly same number of control points
	surface = NURBS::NURBSRectangled::createSheetFromPoints(my_cp);
	other->surface = NURBS::NURBSRectangled::createSheetFromPoints(other_cp);
}

int Sheet::numCtrlPnts()
{
	return surface.mNumUCtrlPoints * surface.mNumVCtrlPoints;
}

void Sheet::refineControlPoints( int nU, int nV /*= 0*/ )
{
	Array2D_Vector3 my_cp;

	NURBS::NURBSRectangled mySurface = surface;

	int k = 0;

	// Refine U
	k = nU - mySurface.mNumUCtrlPoints;
	my_cp = mySurface.simpleRefine(k, 0);
	mySurface = NURBS::NURBSRectangled::createSheetFromPoints(my_cp);

	// Refine V
	k = nV - mySurface.mNumVCtrlPoints;
	my_cp = mySurface.simpleRefine(k, 1);
	surface = NURBS::NURBSRectangled::createSheetFromPoints(my_cp);
}

int Sheet::numUCtrlPnts()
{
	return surface.mNumUCtrlPoints;
}

int Sheet::numVCtrlPnts()
{
	return surface.mNumVCtrlPoints;

}


NURBS::NURBSCurved Sheet::convertToNURBSCurve( Vector3d p, Vector3d dir )
{
	Vector3d uDir = surface.mCtrlPoint.front().back() - surface.mCtrlPoint.front().front();
	Vector3d vDir = surface.mCtrlPoint.back().front() - surface.mCtrlPoint.front().front();
	uDir.normalize(); vDir.normalize();

	double uDot = dot(uDir, dir);
	double vDot = dot(vDir, dir);

	std::vector<Vector3d> curvePoints;
	std::vector<double> curveWeights;
	QString curveDirection;

	// along U
	if (abs(uDot) > abs(vDot))
	{
		// get the row that closet to p
		double minDis = DBL_MAX;
		double minID = -1;
		std::vector<Vector3d> uPoints0 = surface.GetControlPointsU(0);
		for (int i = 0; i < (int)uPoints0.size(); i++)
		{
			double dis = (p - uPoints0[i]).norm();
			if (dis < minDis)	{minDis = dis;	minID = i;}
		}

		curvePoints = surface.GetControlPointsV(minID);
		curveWeights = surface.GetControlWeightsV(minID);
		curveDirection = "positiveU";

		if (uDot < 0)
		{
			std::reverse(curvePoints.begin(), curvePoints.end());
			std::reverse(curveWeights.begin(), curveWeights.end());
			curveDirection = "negativeU";
		}
	}
	// along V
	else
	{
		// get the row that closet to p
		double minDis = DBL_MAX;
		double minID = -1;
		std::vector<Vector3d> vPoints0 = surface.GetControlPointsV(0);
		for (int i = 0; i < (int)vPoints0.size(); i++)
		{
			double dis = (p - vPoints0[i]).norm();
			if (dis < minDis)	{minDis = dis;	minID = i;}
		}

		curvePoints = surface.GetControlPointsU(minID);
		curveWeights = surface.GetControlWeightsU(minID);
		curveDirection = "positiveV";

		if (vDot < 0)
		{
			std::reverse(curvePoints.begin(), curvePoints.end());
			std::reverse(curveWeights.begin(), curveWeights.end());
			curveDirection = "negativeV";
		}
	}

	this->property["curveDirection"] = curveDirection;
    return NURBS::NURBSCurved(curvePoints, curveWeights);
}

SheetEncoding Sheet::encodeSheet( Sheet * sheet, Vector3 origin, Vector3 X, Vector3 Y, Vector3 Z )
{
	SheetEncoding cpCoords;
	Array1D_Vector3 controlPoints = sheet->controlPoints();

	for(int i = 0; i < (int)controlPoints.size(); i++)
	{
		Vector3 p = controlPoints[i];
		Vector3 dir = p - origin;

		// Parameters: offset, theta, psi
		Array1D_Real params(3,0);
		if(dir.norm() > 0){
			params[0] = dir.norm();
			dir = dir.normalized();
		}
		globalToLocalSpherical(X,Y,Z, params[1], params[2], dir);

		cpCoords[i] = params;
	}

	return cpCoords;
}

Array1D_Vector3 Sheet::decodeSheet( SheetEncoding cpCoords, Vector3 origin, Vector3 X, Vector3 Y, Vector3 Z )
{
	Array1D_Vector3 pnts( cpCoords.size(), Vector3(0,0,0) );

	for(int i = 0; i < (int)pnts.size(); i++)
	{
		double offset = cpCoords[i][0];
		double theta = cpCoords[i][1];
		double psi = cpCoords[i][2];

		Vector3 dir(0,0,0);
		localSphericalToGlobal(X,Y,Z,theta,psi,dir);

		pnts[i] = origin + (dir * offset);
	}

	return pnts;
}

void Sheet::deformTo( const Vector4d & handle, const Vector3 & to, bool isRigid )
{
	Q_UNUSED(isRigid)

	Vector3 p = position( handle );

	double diff = (p - to).norm();
	if(diff < DECODE_ZERO_THRESHOLD) return;

	//if( isRigid )
	{
		Vector3d delta = to - p;
		this->moveBy( delta );
		return;
	}
}


void Sheet::deformTwoHandles( Vector4d& handleA, Vector3 newPosA, Vector4d& handleB, Vector3 newPosB )
{
	Vector3d oldA = position(handleA);
	Vector3d oldB = position(handleB);

	// Numerical checks
	double diffOld = (oldA - oldB).norm();	if(diffOld < DECODE_ZERO_THRESHOLD) return;
	double diff = (newPosA - newPosB).norm();	if(diff < DECODE_ZERO_THRESHOLD) return;

	this->moveBy(newPosA - oldA);

	Vector3d oldX = (oldB - oldA).normalized();
	Vector3d oldY = orthogonalVector(oldX);
	Vector3d oldZ = cross(oldX, oldY);
	SheetEncoding SE = encodeSheet(this, newPosA, oldX, oldY, oldZ);

	double scale = (newPosB - newPosA).norm() / (oldB - oldA).norm();
	foreach(int i, SE.keys()) SE[i][0] *= scale;

	Vector3d newX = (newPosB - newPosA).normalized();
	qglviewer::Vec tempOldX(oldX), tempNewX(newX);
	const qglviewer::Quaternion q(tempOldX, tempNewX);

	qglviewer::Vec Y = q * Vec(oldY);
	qglviewer::Vec Z = q * Vec(oldZ);

	this->setControlPoints( decodeSheet( SE, newPosA, newX, Vector3d(Y[0], Y[1], Y[2]), Vector3d(Z[0], Z[1], Z[2]) ) );
}
