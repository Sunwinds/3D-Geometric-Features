#pragma once

#include <assert.h>
#include <QVariant>
#include <QMap>
#include <QSet>
#include <QFileInfo>
#include <QDir>

#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#undef max
#define qRanged(min, v, max) ( qMax(min, qMin(v, max)) )

typedef QMap< QString, QVariant > PropertyMap;

/** Sorts a vector and returns index of the sorted values
 * \param Index Contains the index of sorted values in the original vector
 * \param data The vector to be sorted
 */
template<class T>
void paired_sort(std::vector<unsigned int> & Index, std::vector<T> & data, bool isReverse = false)
{
    // A vector of a pair which will contain the sorted value and its index in the original array
    std::vector< std::pair<T,unsigned int> > IndexedPair;
    IndexedPair.resize(data.size());
    for(unsigned int i=0;i<IndexedPair.size();++i)
    {
        IndexedPair[i].first = data[i];
        IndexedPair[i].second = i;
    }
    std::sort(IndexedPair.begin(),IndexedPair.end());
    Index.resize(data.size());
    for(size_t i = 0; i < Index.size(); ++i) Index[i] = IndexedPair[i].second;
	for(size_t i = 0; i < Index.size(); ++i) data[i] = IndexedPair[i].first;

	if(isReverse){
		std::reverse(Index.begin(), Index.end());
		std::reverse(data.begin(), data.end());
	}
}

// Sort by QMap second value
template<class F, class S>
bool sortByFirst(const QPair<F,S>& e1, const QPair<F,S>& e2) {
	return e1.first < e2.first;
}

template<class F, class S>
QList< QPair<S, F> > sortQMapByValue(const QMap<F,S> & map)
{
	QList< QPair<S, F> > result;

	// Append items to a list
	QMapIterator<F, S> i(map);
	while (i.hasNext()) {
		i.next();
		result.push_back(qMakePair(i.value(), i.key()));
	}

	// Sort that list
	qSort(result.begin(), result.end(), sortByFirst<S,F>);

	return result;
}

static inline std::vector<Eigen::Vector3d> noFrame(){
	return std::vector<Eigen::Vector3d>();
}

static inline Eigen::Vector3d orthogonalVector(const Eigen::Vector3d& n) {
	if ((abs(n.y()) >= 0.9 * abs(n.x())) &&
		abs(n.z()) >= 0.9 * abs(n.x())) return Eigen::Vector3d(0.0, -n.z(), n.y());
	else if ( abs(n.x()) >= 0.9 * abs(n.y()) &&
		abs(n.z()) >= 0.9 * abs(n.y()) ) return Eigen::Vector3d(-n.z(), 0.0, n.x());
	else return Eigen::Vector3d(-n.y(), n.x(), 0.0);
}

// Coordinates utility functions
static inline Eigen::Vector4d coord(double u = 0, double v = 0)	{ return Eigen::Vector4d(u, v, 0, 0); }
static inline Eigen::Vector4d inverseCoord(const Eigen::Vector4d& c)		{ return Eigen::Vector4d(1 - c.x(), 1 - c.y(), 0,0); }
static inline std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > inverseCoords(const std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> >& fromCoords)	{ 
	std::vector<Eigen::Vector4d, Eigen::aligned_allocator<Eigen::Vector4d> > invertedCoords;

	for(int i = 0; i < (int) fromCoords.size(); i++)
	{
		Eigen::Vector4d coord = fromCoords[i];
		invertedCoords.push_back( inverseCoord(coord) );
	}
	return invertedCoords; 
}
static inline bool isSameHalf(const Eigen::Vector4d& A, const Eigen::Vector4d& B){
	bool result = false;
	if( (A[0] <= 0.5 && B[0] <= 0.5) || (A[0] > 0.5 && B[0] > 0.5) )
		result = true;
	return result;
}

template<class VectorType>
static inline double signedAngle(const VectorType &a, const VectorType &b, const VectorType &axis)
{
	if(axis.norm() == 0) qDebug() << "warning: zero axis";
	double cosAngle = a.normalized().dot(b.normalized());
	double angle = acos( qRanged(-1.0, cosAngle, 1.0) );
	VectorType c = a.cross(b);
	if (c.dot(axis) < 0) return -angle;
	return angle;
}

template<class VectorType>
static inline VectorType rotatedVec(const VectorType & v, double theta, const VectorType & axis)
{
	if(theta == 0) return v;
	return (v * cos(theta) + axis.cross(v) * sin(theta) + axis * axis.dot(v) * (1 - cos(theta)));
}

/// Spherical Coordinates:
/* From Spherical-coordinates to Vector 'v' */
template<class VectorType, class ParamType>
static inline void localSphericalToGlobal( VectorType X, VectorType Y, VectorType Z, ParamType theta, ParamType psi, VectorType &v )
{
	Q_UNUSED(X);
	v = rotatedVec(Z, theta, Y);
	v = rotatedVec(v, psi, Z);
}

/* Encode vector 'v' to theta-psi in spherical coordinates of XYZ  */
template<class VectorType, class ParamType>
static inline void globalToLocalSpherical( VectorType X, VectorType Y, VectorType Z, ParamType &theta, ParamType &psi, VectorType v )
{
	// Theta: angle from Z [0, PI]
	// Psi: angle from X on XY plane [0, 2*PI)
	double dotZ = v.dot(Z);
	theta = acos( qRanged(-1.0, dotZ, 1.0) );

	double dotX = v.dot(X);
	double dotY = v.dot(Y);
	VectorType proj_v = (dotX * X + dotY * Y).normalized();
	psi = signedAngle(X, proj_v, Z);
}

static inline Eigen::Vector3d pointOnPlane(Eigen::Vector3d p, Eigen::Vector3d plane_normal, double plane_d = 0)
{
    double t = plane_normal.dot(p) - plane_d;
	return p - (t * plane_normal);
}

#define	POINT_ID_RANGE 1000
#define NODE_ID_RANGE	100

static void getIndicesFromSelectedName(int selectedName, int &gID, int &nID, int &pID)
{
	pID = selectedName % POINT_ID_RANGE;
	selectedName /= POINT_ID_RANGE;

	nID = selectedName % NODE_ID_RANGE;
	gID = selectedName / NODE_ID_RANGE;
}

typedef std::pair< QVector<QString>, QVector<QString> > PART_LANDMARK;
typedef std::pair< int, int > POINT_ID;
typedef std::pair< QVector<POINT_ID>, QVector<POINT_ID> > POINT_LANDMARK;

static bool operator<(const PART_LANDMARK & a, const PART_LANDMARK & b){
	QString keyA = a.first.front() + a.second.front();
	QString keyB = b.first.front() + b.second.front();
	return keyA < keyB;
}

#define AlphaBlend(alpha, start, end) ( ((1-alpha) * start) + (alpha * end) )

// Spatial Hausdorff distance
static double supInfDistance( const std::vector<Eigen::Vector3d> &A, const std::vector<Eigen::Vector3d> &B )
{
	double supinfDis = -1;
	for (int i = 0; i < (int)A.size(); i++)
	{
		double infDis = DBL_MAX;
		for (int j = 0; j < (int)B.size(); j++)
		{
			double dis = (A[i] - B[j]).norm();

			if (dis < infDis)
				infDis = dis;
		}

		if (infDis > supinfDis)
			supinfDis = infDis;
	}

	return supinfDis;
}

static double HausdorffDistance( const std::vector<Eigen::Vector3d> &A, const std::vector<Eigen::Vector3d> &B )
{
	double ABDis = supInfDistance(A, B);
	double BADis = supInfDistance(B, A);

	return std::max(ABDis, BADis);
}

static std::vector<Eigen::Vector3d> refineByNumber(const std::vector<Eigen::Vector3d> & fromPnts, int targetNumber)
{
	// Refine until targetNumber is achieved
	std::vector<Eigen::Vector3d> newPnts = fromPnts;
	while(targetNumber > (int)newPnts.size())
	{
		// Find index of largest edge and split
		double maxDist = -DBL_MAX;
		int idx = -1;
		for(int i = 1; i < (int) newPnts.size(); i++){
			double dist = (newPnts[i] - newPnts[i-1]).norm();
			if(dist > maxDist){
				maxDist = dist;
				idx = i;
			}
		}

		Eigen::Vector3d midPoint = (newPnts[idx] + newPnts[idx-1]) / 2.0;

		// Insert new point
		newPnts.insert( newPnts.begin() + (idx), midPoint );
	}

	// Contract to get to target
	while((int)newPnts.size() > targetNumber)
	{
		// Find index of shortest edge
		double minDist = DBL_MAX;
		int idx = -1;
		for(int i = 1; i < (int) newPnts.size(); i++){
			double dist = (newPnts[i] - newPnts[i-1]).norm();
			if(dist < minDist){
				minDist = dist;
				idx = i;
			}
		}

		// Remove not at end neighbor
		int toRemoveIdx = idx;
		if(idx == newPnts.size()-1) toRemoveIdx--;

		newPnts.erase( newPnts.begin() + toRemoveIdx );
	}

	return newPnts;
}

static std::vector<Eigen::Vector3d> refineByResolution(const std::vector<Eigen::Vector3d> & fromPnts, double resolution)
{
	std::vector<Eigen::Vector3d> resampled_polyline;

	// Parameterize line by finding segments
	double totalLineLength = 0.0;

	typedef std::pair<double,int> LengthSegment;
	std::vector< LengthSegment > segments;

	// Add a start segment
	segments.push_back(LengthSegment(totalLineLength, 1));

	for(int i = 1; i < (int)fromPnts.size(); i++)
	{
		double curLen = (fromPnts[i] - fromPnts[i-1]).norm();
		totalLineLength += curLen;
		segments.push_back( LengthSegment(totalLineLength, i) );
	}

	if(totalLineLength == 0.0)
		return resampled_polyline;

	double numSegments = totalLineLength / resolution;

	// Re-sample line
	for(int i = 0; i < numSegments; i++)
	{
		double t = totalLineLength * (double(i) / (numSegments-1));

		std::vector< LengthSegment >::iterator it = lower_bound(segments.begin(), 
			segments.end(), LengthSegment(qMin(t, totalLineLength), -1));

		// Find start
		int idx = it->second;

		double seg_range = segments[idx].first - segments[idx-1].first;
		double seg_start = segments[idx-1].first;

		double alpha = (t - seg_start) / seg_range;

		resampled_polyline.push_back( ( (1-alpha) * fromPnts[idx-1] ) + (alpha * fromPnts[idx]) );
	}

	return resampled_polyline;
}

static std::vector<Eigen::Vector3d> smoothPolyline(const std::vector<Eigen::Vector3d> & fromPnts, int num_iterations)
{
	std::vector<Eigen::Vector3d> pnts = fromPnts;

	// Laplacian smoothing - fixed ends
	for(int itr = 0; itr < num_iterations; itr++)
	{
		std::vector<Eigen::Vector3d> newPos(pnts.size(), Eigen::Vector3d(0,0,0));
		newPos[0] = pnts[0];

		for(int i = 1; i < (int)pnts.size() - 1; i++)
			newPos[i] = (pnts[i-1] + pnts[i+1]) * 0.5;

		newPos.back() = pnts.back();
		pnts = newPos;
	}

	return pnts;
}

/* Generic vector operations */
template<typename T>
static inline std::vector<T> concat(const std::vector<T> & A, const std::vector<T> & B){
	std::vector<T> dest;
	dest.insert(dest.end(), A.begin(), A.end());
	dest.insert(dest.end(), B.begin(), B.end());
	return dest;
}
template<typename T>
static inline std::vector<T> concat(const std::vector<T> & A, const std::vector<T> & B, const std::vector<T> & C){
	return concat(concat(A,B), C);
}

template<typename T>
static inline std::vector<T> subvec(const std::vector<T> & V, int start, int length = -1){
	if(length < 0) length = V.size() - start;
    typename std::vector<T>::const_iterator first = V.begin() + start;
    typename std::vector<T>::const_iterator last = V.begin() + qMin(start + length, (int)V.size());
	return std::vector<T> (first, last);
}

template<typename T>
static inline std::vector<T> leftvec(const std::vector<T> & V, int length){
	return subvec(V, 0, qMin((int)V.size(), length));
}
template<typename T>
static inline std::vector<T> rightvec(const std::vector<T> & V, int length){
	return subvec(V, qMax(0, (int)V.size() - length));
}
template<typename T>
static inline std::vector<T> chopvec(const std::vector<T> & V, int length){
	return subvec(V, 0, V.size() - length);
}

template<typename T>
static inline std::vector<T> reversedvec(const std::vector<T> & V){
	std::vector<T> reversed = V;
	std::reverse(reversed.begin(), reversed.end());
	return reversed;
}

template <typename T>
QList<T> reversed( const QList<T> & in ) {
	QList<T> result;
	result.reserve( in.size() );
	std::reverse_copy( in.begin(), in.end(), std::back_inserter( result ) );
	return result;
}

template <typename T>
QVector<T> reversed( const QVector<T> & in ) {
	QVector<T> result;
	result.reserve( in.size() );
	std::reverse_copy( in.begin(), in.end(), std::back_inserter( result ) );
	return result;
}

template<typename T>
static inline QVector<T> sumQVec( const QVector<T> & V ){
	T vector_sum;
	for(int i = 0; i < (int)V.size(); i++) vector_sum += V[i];
	return vector_sum;
}
template<typename T>
static inline T sumvec( const std::vector<T> & V ){
	T vector_sum;
	for(int i = 0; i < (int)V.size(); i++) vector_sum += V[i];
	return vector_sum;
}

template <class Type>
static inline std::vector<std::vector<Type> > transpose(const std::vector<std::vector<Type> > data)
{
	// this assumes that all inner vectors have the same size and
	// allocates space for the complete result in advance
	std::vector<std::vector<Type> > result(data[0].size(), std::vector<Type>(data.size()));
	for (int i = 0; i < (int)data[0].size(); i++){
		for (int j = 0; j < (int)data.size(); j++)
			result[i][j] = data[j][i];
	}

	return result;
}

static void saveOBJ(SurfaceMesh::Model * mesh, QString filename)
{
	QFile file(filename);

	// Create folder
	QFileInfo fileInfo(file.fileName());
	QDir d(""); d.mkpath(fileInfo.absolutePath());

	// Open for writing
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	
	QTextStream out(&file);
	out << "# NV = " << mesh->n_vertices() << " NF = " << mesh->n_faces() << "\n";
	SurfaceMesh::Vector3VertexProperty points = mesh->vertex_property<Eigen::Vector3d>("v:point");
	foreach( SurfaceMesh::Vertex v, mesh->vertices() )
		out << "v " << points[v][0] << " " << points[v][1] << " " << points[v][2] << "\n";
	foreach( SurfaceMesh::Face f, mesh->faces() ){
		out << "f ";
		Surface_mesh::Vertex_around_face_circulator fvit=mesh->vertices(f), fvend=fvit;
		do{	out << (((Surface_mesh::Vertex)fvit).idx()+1) << " ";} while (++fvit != fvend);
		out << "\n";
	}
	file.close();
}

static void combineMeshes( QStringList filenames, QString outputFilename )
{
	QFile file(outputFilename);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QTextStream out(&file);

	int v_offset = 0;

	foreach(QString filename, filenames)
	{
		QFileInfo fileInfo(filename);

		SurfaceMesh::SurfaceMeshModel * m = new SurfaceMesh::SurfaceMeshModel;
		m->read( qPrintable(filename) );
		SurfaceMesh::Vector3VertexProperty points = m->vertex_property<Eigen::Vector3d>("v:point");

		out << "# Start of mesh " << fileInfo.baseName() << "\n";

		// Vertices
		foreach( SurfaceMesh::Vertex v, m->vertices() )
			out << "v " << points[v][0] << " " << points[v][1] << " " << points[v][2] << "\n";

		// Triangles
		out << "g " << fileInfo.baseName() << "\n";
		foreach( SurfaceMesh::Face f, m->faces() ){
			out << "f ";
			Surface_mesh::Vertex_around_face_circulator fvit = m->vertices(f), fvend = fvit;
			do{	out << (((Surface_mesh::Vertex)fvit).idx() + 1 + v_offset) << " ";} while (++fvit != fvend);
			out << "\n";
		}

		v_offset += m->n_vertices();

		out << "# End of mesh " << fileInfo.baseName() ;
	}
}

/*
General Bezier curve
Number of control points is n+1
0 <= mu < 1    IMPORTANT, the last point is not computed
*/
template<class VectorType>
VectorType BezierCurve(const std::vector<VectorType> & p, double mu){
	int k,kn,nn,nkn;
	double blend,muk,munk;
	VectorType b(0,0,0);

	int n = p.size() - 1;
	muk = 1;
	munk = pow(1-mu,(double)n);

	for (k=0;k<=n;k++) {
		nn = n;
		kn = k;
		nkn = n - k;
		blend = muk * munk;
		muk *= mu;
		munk /= (1-mu);
		while (nn >= 1) {
			blend *= nn;
			nn--;
			if (kn > 1) {
				blend /= (double)kn;
				kn--;
			}
			if (nkn > 1) {
				blend /= (double)nkn;
				nkn--;
			}
		}
		b[0] += p[k][0] * blend;
		b[1] += p[k][1] * blend;
		b[2] += p[k][2] * blend;
	}

	return(b);
}

template<class VectorType>
std::vector<VectorType> BezierArc(const VectorType& pointA, const VectorType& pointB, int numSteps = 20, double bend = 1.0)
{
	std::vector<VectorType> arc, cpoints;

	VectorType midPoint = 0.5 * (pointA + pointB);
	VectorType axis(0,0,1);
	VectorType d(pointB - midPoint);
	VectorType delta = rotatedVec(d, M_PI_2, axis);

	// Set control points for Bezier curve
	cpoints.push_back(pointA);
	cpoints.push_back(midPoint + (delta * bend));
	cpoints.push_back(pointB);

	// Bound
	numSteps = qMax(3, numSteps);

	for(int i = 0; i < numSteps; i++){
		double t = double(i) / numSteps;
		arc.push_back( BezierCurve(cpoints, t) );
	}

	// Last point on arc
	arc.push_back(pointB);

	return arc;
}

static inline std::vector<bool> subsampleMask(int samples, int entire){
	std::vector<bool> mask = std::vector<bool>(entire, false);
	int difference = entire - samples;
	if(difference <= 0) return std::vector<bool>(entire, true);

	// Subtract the difference and shuffle
	for(int i = 0; i < samples; i++) mask[i] = true;
	std::random_shuffle(mask.begin(), mask.end());
	return mask;
}
