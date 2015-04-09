#define _USE_MATH_DEFINES

#include <QStack>
#include "SurfaceMeshHelper.h"
#include "segmentation.h"

#include "combinatorics.h"

#include "CurveskelModel.h"
#include "CurveskelHelper.h"

#include "disjointset.h"

#define qRanged(_min, v, _max) ( qMax(_min, qMin(v, _max)) )
#define RADIANS(deg)    ((deg)/180.0 * M_PI)

typedef QPair<Vector3,Vector3> LineSegment;

template<class VectorType>
static inline double signedAngle(const VectorType &a, const VectorType &b, const VectorType &axis){
	if(axis.norm() == 0) qDebug() << "warning: zero axis";
	double cosAngle = a.normalized().dot(b.normalized());
	double angle = acos( qRanged(-1.0, cosAngle, 1.0) );
	VectorType c = a.cross(b);
	if (c.dot(axis) < 0) return -angle;
	return angle;
}

template<class T1, class T2>
bool pairCompare(const std::pair<T1,T2> & x,const std::pair<T1,T2> & y) {
	return x.second < y.second; 
}

template<class T>
typename T::const_iterator map_max_element(const T & A){
	typedef typename T::value_type pair_type;
	typedef typename pair_type::first_type K;
	typedef typename pair_type::second_type V;
	return std::max_element(A.begin(), A.end(), pairCompare<K,V>);
}

class IntCombinations{ 
public: 
	std::vector< std::vector<int> > * result;
	explicit IntCombinations(std::vector< std::vector<int> > * result) : result(result) {}
	template <class It>	bool operator()(It first, It last){ 
		std::vector<int> s;
		for (++first; first != last; ++first){ CurveskelTypes::Edge e = *first; int x = e.idx(); s.push_back(x); };
		result->push_back(s);
		return false; 
	} 
};

double dist_Point_to_Segment( Vector3 P, Vector3 SP0, Vector3 SP1){
	Vector3 v = SP1 - SP0;
	Vector3 w = P - SP0;
	double c1 = dot(w,v);
	if ( c1 <= 0 ) return (P - SP0).norm();
	double c2 = dot(v,v);
	if ( c2 <= c1 )	return (P - SP1).norm();
	double b = c1 / c2;
	Point Pb = SP0 + b * v;
	return (P - Pb).norm();
}

inline double minAngle(Face f, SurfaceMeshModel * ofMesh)
{
	double minAngle(DBL_MAX);
	Vector3VertexProperty pts = ofMesh->vertex_property<Vector3>(VPOINT);

	SurfaceMesh::Model::Halfedge_around_face_circulator h(ofMesh, f), eend = h;
	do{ 
		Vector3 a = pts[ofMesh->to_vertex(h)];
		Vector3 b = pts[ofMesh->from_vertex(h)];
		Vector3 c = pts[ofMesh->to_vertex(ofMesh->next_halfedge(h))];

		double d = dot((b-a).normalized(), (c-a).normalized());
		double angle = acos(qRanged(-1.0, d, 1.0));

		minAngle = qMin(angle, minAngle);
	} while(++h != eend);

	return minAngle;
}

inline double signedVolumeOfTriangle(Vector3 p1, Vector3 p2, Vector3 p3) {
	double v321 = p3[0]*p2[1]*p1[2];
	double v231 = p2[0]*p3[1]*p1[2];
	double v312 = p3[0]*p1[1]*p2[2];
	double v132 = p1[0]*p3[1]*p2[2];
	double v213 = p2[0]*p1[1]*p3[2];
	double v123 = p1[0]*p2[1]*p3[2];
	return (1.0/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

inline double meshVolume(SurfaceMeshModel * m){
	double vol = 0;
	Vector3VertexProperty pts = m->vertex_property<Vector3>(VPOINT);
	for(Face f : m->faces()){
		std::vector<Vector3> p;
		for(Vertex v : m->vertices(f)) p.push_back(pts[v]);
		vol += signedVolumeOfTriangle(p[0],p[1],p[2]);
	}
	return vol;
}

void segmentation::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichBool("Visualize", true, "Visualize"));

	pars->addParam(new RichInt("PreSmoothing", 0, "Smoothing"));
	pars->addParam(new RichFloat("H", 0.4f, "H"));
	pars->addParam(new RichFloat("P", 0.2f, "P"));
}

void segmentation::applyFilter(RichParameterSet *pars)
{
	// Make a copy
	SurfaceMeshModel* copy = mesh()->clone();
	copy->name = "copy";
	copy->isVisible = false;

	mesh()->isVisible = false;

	Vector3VertexProperty points = copy->vertex_coordinates();
	ScalarVertexProperty segment = copy->vertex_property<Scalar>("v:segment", 0);

	document()->addModel(copy);

	// pre-processing
	{
		Vector3VertexProperty points = mesh()->vertex_coordinates();

		// Laplacian smoothing
		SurfaceMeshHelper h(mesh());
		h.smoothVertexProperty<Vector3>(VPOINT, pars->getInt("PreSmoothing"), Vector3(0,0,0));
	}

	// Medial axis
	FilterPlugin * voromatPlugin = pluginManager()->getFilter("Voronoi based MAT");
	RichParameterSet * voromat_params = new RichParameterSet;
	voromatPlugin->initParameters(voromat_params);
	voromatPlugin->applyFilter(voromat_params);

	// Skeletonize
	FilterPlugin * mcfPlugin = pluginManager()->getFilter("MCF Skeletonization");
	RichParameterSet * mcf_params = new RichParameterSet;
	mcfPlugin->initParameters( mcf_params );
	mcf_params->setValue("omega_H_0", pars->getFloat("H"));
	mcf_params->setValue("omega_P_0", pars->getFloat("P"));

	bool almostCurve = false;
	while(!almostCurve){
		almostCurve = true;
		double volume = meshVolume(mesh());

		mcfPlugin->applyFilter( mcf_params );

		if(volume - meshVolume(mesh()) > 1e-6) 
			almostCurve = false;
	}

	// Convert to curve
	FilterPlugin * makeskeletonPlugin = pluginManager()->getFilter("SurfaceMesh to Skeleton");
	RichParameterSet * makeskeleton_params = new RichParameterSet;
	makeskeletonPlugin->initParameters(makeskeleton_params);
	makeskeletonPlugin->applyFilter(makeskeleton_params);

	// Visualize joints and endpoints
	CurveskelTypes::CurveskelModel* skeleton = (CurveskelTypes::CurveskelModel*)document()->getModel("skeleton");
	CurveskelTypes::CurveskelModel::Vertex_property<CurveskelTypes::Point> skel_pnts = skeleton->vertex_property<CurveskelTypes::Point>("v:point");
	
	CurveskelTypes::CurveskelModel::Edge_property<int> skel_seg_id = skeleton->edge_property<int>("e:seg_id", 0);
	std::vector<QColor> rndColors;

	std::set<CurveskelTypes::Vertex> junctions = skeleton->junctions(), endpoints = skeleton->endpoints();

	foreach(CurveskelTypes::Vertex v, junctions){
		Vector3 p(skel_pnts[v][0], skel_pnts[v][1], skel_pnts[v][2]);
		drawArea()->drawPoint( p, 20, Qt::red );
	}

	foreach(CurveskelTypes::Vertex v, endpoints){
		Vector3 p(skel_pnts[v][0], skel_pnts[v][1], skel_pnts[v][2]);
		drawArea()->drawPoint( p, 20, Qt::blue );
	}

	/// Extract and visualize curves:
	
	CurveskelTypes::CurveskelModel::Edge_property<bool> evisisted = skeleton->edge_property<bool>("e:visisted", false);

	int sid = 0;

	foreach( CurveskelTypes::Edge e, skeleton->edges() )
	{
		if(evisisted[e]) continue;

		QStack<CurveskelTypes::Edge> tovisit;
		tovisit.push(e);

		while(!tovisit.isEmpty())
		{
			CurveskelTypes::Edge cur = tovisit.pop();
			evisisted[cur] = true;
			CurveskelTypes::Vertex v0 = skeleton->vertex(cur, 0);
			CurveskelTypes::Vertex v1 = skeleton->vertex(cur, 1);

			bool v0_stop = (junctions.find(v0) != junctions.end());
			bool v1_stop = (junctions.find(v1) != junctions.end());

			if(!v0_stop) for(CurveskelTypes::Edge j : skeleton->edges_of(v0)) { if(!evisisted[j]) tovisit.push(j); }
			if(!v1_stop) for(CurveskelTypes::Edge j : skeleton->edges_of(v1)) { if(!evisisted[j]) tovisit.push(j); }

			skel_seg_id[cur] = sid;
		}

		sid++;
	}

	// Generate nice random colors
	for(int i = 0; i < sid * 3; i++) rndColors.push_back(qRandomColor2());
	std::random_shuffle(rndColors.begin(), rndColors.end());

	// Collect segments
	QMap< int, std::set<CurveskelTypes::Edge> > segments;
	foreach( CurveskelTypes::Edge e, skeleton->edges() )
	{
		int sid = skel_seg_id[e];
		segments[ sid ].insert( e );

	}
		
	// Convert to line segments
	std::map<int, LineSegment> lines;
	foreach(int s, segments.keys())
	{
		std::vector<Vector3> two_ends;

		foreach( CurveskelTypes::Edge e, segments[s] )
		{
			CurveskelTypes::Vertex v0 = skeleton->vertex(e, 0);
			CurveskelTypes::Vertex v1 = skeleton->vertex(e, 1);

			bool v0_special = (junctions.find(v0) != junctions.end()) || (endpoints.find(v0) != endpoints.end());
			bool v1_special = (junctions.find(v1) != junctions.end()) || (endpoints.find(v1) != endpoints.end());

			Vector3 p0(skel_pnts[skeleton->vertex(e,0)][0],skel_pnts[skeleton->vertex(e,0)][1],skel_pnts[skeleton->vertex(e,0)][2]);
			Vector3 p1(skel_pnts[skeleton->vertex(e,1)][0],skel_pnts[skeleton->vertex(e,1)][1],skel_pnts[skeleton->vertex(e,1)][2]);

			if(v0_special) two_ends.push_back(p0);
			if(v1_special) two_ends.push_back(p1);
		}

		assert(two_ends.size() == 2);

		lines[s] = LineSegment(two_ends.front(), two_ends.back());
	}

	// Merge not so different segments
	DisjointSet set( segments.size() );

	foreach(CurveskelTypes::Vertex v, junctions)
	{
		std::vector<CurveskelTypes::Edge> edges = skeleton->edges_of(v);

		for(int j = 0; j < edges.size(); j++){
			for(int k = j + 1; k < edges.size(); k++){
				CurveskelTypes::Edge e1 = edges[j];
				CurveskelTypes::Edge e2 = edges[k];

				#define furthest(p, a, b) (( (p-a).norm() > (p-b).norm() ) ? a : b)
					
				Vector3 p0 = Vector3(skel_pnts[v][0], skel_pnts[v][1], skel_pnts[v][2]);
				Vector3 p1 = furthest(p0, lines[ skel_seg_id[e1] ].first, lines[ skel_seg_id[e1] ].second);
				Vector3 p2 = furthest(p0, lines[ skel_seg_id[e2] ].first, lines[ skel_seg_id[e2] ].second);

				Vector3 v1 = (p1-p0).normalized();
				Vector3 v2 = (p2-p0).normalized();

				if(acos( dot(v1,v2) ) > M_PI - (M_PI * 0.27)){
					//drawArea()->drawSegment( p2, p1, 5, rndColors[ skel_seg_id[e1] ] );
					set.Union( skel_seg_id[e1], skel_seg_id[e2] );
				}
			}
		}
	}

	// Merge segments
	foreach(int s, segments.keys()){
		int replaceTo = set.Parent[s];

		foreach( CurveskelTypes::Edge e, segments[s] )
			skel_seg_id[e] = replaceTo;
	}

	// Segment mesh copy
	{
		// Assign vertices to the closest segment
		foreach(Vertex v, copy->vertices()){
			CurveskelTypes::Edge bestEdge;
			double minDist = DBL_MAX;

			foreach(CurveskelTypes::Edge e, skeleton->edges()){
				Vector3 v0(skel_pnts[skeleton->vertex(e,0)][0],skel_pnts[skeleton->vertex(e,0)][1],skel_pnts[skeleton->vertex(e,0)][2]);
				Vector3 v1(skel_pnts[skeleton->vertex(e,1)][0],skel_pnts[skeleton->vertex(e,1)][1],skel_pnts[skeleton->vertex(e,1)][2]);

				double dist = dist_Point_to_Segment(points[v], v0, v1);

				if(dist < minDist){
					minDist = dist;
					bestEdge = e;
				}
			}

			segment[v] = skel_seg_id[bestEdge];
		}

		// Find faces of each segment
		QMap< int, QVector<Face> > seg;
		foreach(Face f, copy->faces()){
			std::map< int,int > counter;
			foreach(Vertex v, copy->vertices(f)) counter[segment[v]]++;
			seg[ (*map_max_element(counter)).first ].push_back(f);
		}

		// Extract components
		foreach(int s, seg.keys()){
			QVector<Face> faces = seg[s];

			SurfaceMeshModel * m = copy->clone(faces.toStdVector());
			m->name = QString("seg_%1").arg(s);
			m->color = rndColors[s];

			document()->addModel(m);
		}
	}

	// Colorize the segments
	foreach(CurveskelTypes::Edge e, skeleton->edges()){
		Vector3 v0(skel_pnts[skeleton->vertex(e,0)][0],skel_pnts[skeleton->vertex(e,0)][1],skel_pnts[skeleton->vertex(e,0)][2]);
		Vector3 v1(skel_pnts[skeleton->vertex(e,1)][0],skel_pnts[skeleton->vertex(e,1)][1],skel_pnts[skeleton->vertex(e,1)][2]);
		drawArea()->drawSegment(v0,v1, 6, rndColors[ skel_seg_id[e] ]);
	}

	foreach(Vertex v, copy->vertices())
		drawArea()->drawPoint(points[v], 8, rndColors[ segment[v] ]);
}
