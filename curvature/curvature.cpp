#include <numeric> 

#include "SurfaceMeshHelper.h"

#include "principal_curvature.h"
#include "curvature.h"
#include "RenderObjectExt.h"

#include "Monge_via_jet_fitting.h"

#ifdef WIN32
namespace std{  bool isnan(double x){ return _isnan(x); } }
#endif

void curvature::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichInt("Radius", 3, "Radius"));
	pars->addParam(new RichInt("Smoothing", 0, "Smoothing"));
	pars->addParam(new RichBool("Fitting", false, "By fitting"));
	pars->addParam(new RichBool("Clamp", true, "Clamp"));
	pars->addParam(new RichInt("ActiveItem", 0, "Item (mean, ..etc)"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
	pars->addParam(new RichBool("VisualizePoints", false, "Show colored points"));
	pars->addParam(new RichBool("VisualizeVectors", false, "Show directions"));
	pars->addParam(new RichFloat("VisualizeScale", 1.0, "Vectors scale"));
}

static inline Vector3 orthogonalVector(const Vector3& n) {
	if ((abs(n.y()) >= 0.9 * abs(n.x())) &&
		abs(n.z()) >= 0.9 * abs(n.x())) return Vector3(0.0, -n.z(), n.y());
	else if ( abs(n.x()) >= 0.9 * abs(n.y()) &&
		abs(n.z()) >= 0.9 * abs(n.y()) ) return Vector3(-n.z(), 0.0, n.x());
	else return Vector3(-n.y(), n.x(), 0.0);
}

// get the median of an unordered set of numbers of arbitrary type without modifying the
// underlying dataset
template <typename InputIterator>
double median(InputIterator const cbegin, InputIterator const cend){
	typedef std::iterator_traits<InputIterator>::value_type T;
	std::vector<T> dataSorted(cbegin, cend);
	std::sort(dataSorted.begin(), dataSorted.end());
	return dataSorted[dataSorted.size() / 2U];
}

void curvature::collectEnoughRings(Vertex v, const size_t min_nb, std::vector<int> & all)
{
	std::vector<Vertex> current_ring, next_ring;
	SurfaceMeshModel::Vertex_property<int> visited_map = mesh()->vertex_property<int>("v:visit_map", -1);

	foreach(Vertex v, mesh()->vertices()) visited_map[v] = -1;

	//initialize
	visited_map[v] = 0;
	current_ring.push_back(v);
	all.push_back(v.idx());

	int i = 1;

	while ( (all.size() < min_nb) &&  (current_ring.size() != 0) )
	{
		// collect ith ring
		std::vector<Vertex>::iterator it = current_ring.begin(), ite = current_ring.end();

		for(;it != ite; it++){
			SurfaceMeshModel::Halfedge_around_vertex_circulator hedgeb = mesh()->halfedges(*it), hedgee = hedgeb;

			do{
				Vertex vj = mesh()->to_vertex(hedgeb);

				if (visited_map[vj] == -1){
					visited_map[vj] = i;
					next_ring.push_back(vj);
					all.push_back(vj.idx());

					if(all.size() > min_nb) return;
				}

				++hedgeb;
			} while(hedgeb != hedgee);
		}

		//next round must be launched from p_next_ring...
		current_ring = next_ring;
		next_ring.clear();

		i++;
	}
}

void curvature::applyFilter(RichParameterSet *pars)
{
	Eigen::MatrixXd V(mesh()->n_vertices(), 3);
	Eigen::MatrixXi F(mesh()->n_faces(), 3);

	Eigen::MatrixXd NV(mesh()->n_vertices(), 3);
	Eigen::MatrixXd NF(mesh()->n_faces(), 3);

	Vector3VertexProperty points = mesh()->vertex_coordinates();
	Vector3VertexProperty normals = mesh()->vertex_normals(true);
	Vector3FaceProperty fnormals = mesh()->face_normals(true);

	// Fill in vertices
	int vi = 0;
	foreach(Vertex v, mesh()->vertices()) { V.row(vi) = points[v]; NV.row(vi) = normals[v]; vi++; }

	// Fill in triangles
	int fi = 0;
	foreach(Face f, mesh()->faces()){
		int i = 0;
		foreach(Vertex v, mesh()->vertices(f)) F(fi,i++) = v.idx();
		NF.row(fi) = fnormals[f];
		fi++;
	};

	// Save results to mesh
	ScalarVertexProperty k1 = mesh()->vertex_property("v:k1", 0.0);
	ScalarVertexProperty k2 = mesh()->vertex_property("v:k2", 0.0);
	Vector3VertexProperty pd1 = mesh()->vertex_property("v:pd1", Vector3(0,0,0));
	Vector3VertexProperty pd2 = mesh()->vertex_property("v:pd2", Vector3(0,0,0));

	ScalarVertexProperty cmax = mesh()->vertex_property("v:cmax", 0.0);
	ScalarVertexProperty cmin = mesh()->vertex_property("v:cmin", 0.0);
	ScalarVertexProperty cmean = mesh()->vertex_property("v:cmean", 0.0);
    ScalarVertexProperty cgauss = mesh()->vertex_property("v:cgauss", 0.0);

	// Compute average edge length
	double avgEdge = 0;
	SurfaceMeshHelper h(mesh());
	ScalarEdgeProperty elengths = h.computeEdgeLengths();
	foreach(Edge e, mesh()->edges()) avgEdge += elengths[e];
	avgEdge /= mesh()->n_edges();

	if( !pars->getBool("Fitting") )
	{
		// Temporary results
		Eigen::MatrixXd K1(mesh()->n_vertices(), 1);
		Eigen::MatrixXd K2(mesh()->n_vertices(), 1);

		Eigen::MatrixXd PD1(mesh()->n_vertices(), 3);
		Eigen::MatrixXd PD2(mesh()->n_vertices(), 3);

		igl::principal_curvature( V, F, NV, NF, K1, K2, PD1, PD2, pars->getInt("Radius"), 0 );

		foreach(Vertex v, mesh()->vertices()){
			int vidx = v.idx();

			k1[v] = K1(vidx, 0);
			k2[v] = K2(vidx, 0);

			pd1[v] = PD1.row(vidx);
			pd2[v] = PD2.row(vidx);
		}
	}
	else
	{
		std::vector< Vector3 > pnts;
		std::vector< std::vector<int> > neighbours(mesh()->n_vertices());
		std::vector< std::vector<int> > vertex_to_vertices;
		//double scaledRadius = avgEdge * pars->getInt("Radius");
		//igl::adjacency_list(F, vertex_to_vertices);

		foreach(Vertex v, mesh()->vertices())
		{
			// Do a sphere search to collect neighbours
			//igl::CurvatureCalculator::getSphere(V, vertex_to_vertices, v.idx(), scaledRadius, neighbours[v.idx()], 6);
			collectEnoughRings(v, qMax(7, 2 * pars->getInt("Radius")), neighbours[v.idx()]);

			pnts.push_back(points[v]);
		}

		// Fit around each point
		std::vector<Monge_via_jet_fitting::Monge_form> fit = Monge_via_jet_fitting::fit_points(pnts, neighbours);

		// Extract results
		foreach(Vertex v, mesh()->vertices())
		{
			Monge_via_jet_fitting::Monge_form m = fit[v.idx()];
			m.comply_wrt_given_normal(normals[v]);

			k1[v] = m.principal_curvatures(0);
			k2[v] = m.principal_curvatures(1);

			pd1[v] = m.maximal_principal_direction() * k1[v];
			pd2[v] = m.minimal_principal_direction() * k2[v];
		}
	}

	std::vector<double> maxVals(4, -DBL_MAX);
	std::vector<double> minVals(4,  DBL_MAX);


	foreach(Vertex v, mesh()->vertices()){
		if(std::isnan(pd1[v][0])) pd1[v] = Vector3(0,0,0);
		if(std::isnan(pd2[v][0])) pd2[v] = Vector3(0,0,0);

		cmax[v] = k1[v];
		cmin[v] = k2[v];
		cmean[v] = (cmax[v] + cmin[v]) / 2.0;
		cgauss[v] = cmin[v] * cmax[v];
	}

	QStringList nms;
    nms << "v:cmean" << "v:cgauss" << "v:cmax" << "v:cmin";

	if(pars->getInt("Smoothing"))
	{
		foreach(QString name, nms)
		{
			std::string prop_name = name.toStdString();
			ScalarVertexProperty prop = mesh()->vertex_property<Scalar>( prop_name );
			SurfaceMeshHelper(mesh()).smoothVertexProperty<Scalar>( prop_name, pars->getInt("Smoothing") );
		}
	}

	if(pars->getBool("Clamp"))
	{
		foreach(QString name, nms)
		{
			ScalarVertexProperty prop = mesh()->vertex_property<Scalar>( name.toStdString() );

			SurfaceMeshHelper(mesh()).smoothVertexProperty<Scalar>(name.toStdString(), pars->getInt("Smoothing"));

			std::vector<double> values(mesh()->n_vertices(), 0.0);
            foreach(Vertex v, mesh()->vertices()){
				values[v.idx()] = prop[v];
			}

			std::string clamped_name = QString(name + "_clamped").toStdString();
			ScalarVertexProperty clamped = mesh()->vertex_property(clamped_name, 0.0);

			double sum = std::accumulate(values.begin(), values.end(), 0.0);
			double mean = sum / values.size();
			double mad = median(values.begin(), values.end());
			Eigen::VectorXd V(mesh()->n_vertices());
			foreach(Vertex v, mesh()->vertices()) V(v.idx()) = (values[v.idx()] - mean) / (2.0 * mad);
			foreach(Vertex v, mesh()->vertices()) clamped[v] = qBound(-2.0, V(v.idx()), 2.0);
		}
	}

	if(pars->getBool("Visualize")){
		starlab::VectorSoup * vs1 = new starlab::VectorSoup(Qt::red);
		starlab::VectorSoup * vs2 = new starlab::VectorSoup(Qt::blue);
		starlab::PointSoup * ps = new starlab::PointSoup();

		// Change this arg(..) to visualize different values
		std::string prop_name = QString(pars->getBool("Clamp") ? "%1_clamped" : "%1").arg( nms[pars->getInt("ActiveItem")] ).toStdString();

		ScalarVertexProperty tovisualize = mesh()->vertex_property<Scalar>( prop_name );

		double minVal = DBL_MAX;
		double maxVal = -DBL_MAX;
		foreach(Vertex v, mesh()->vertices()){
			minVal = qMin(minVal, tovisualize[v]);
			maxVal = qMax(maxVal, tovisualize[v]);
		}

		Vector3VertexProperty vcolor = mesh()->vertex_property<Vector3>(VCOLOR);

		foreach(Vertex v, mesh()->vertices())
		{
			//QColor color = starlab::qtJetColorMap( clamped[v], -2.0, 2.0 );
			QColor color = pars->getBool("Clamp") ? starlab::qtJetColor( tovisualize[v], -2.0, 2.0 ) : 
													starlab::qtJetColor( tovisualize[v], minVal, maxVal );

			vcolor[v] = Vector3(color.redF(), color.greenF(), color.blueF());

			if(pars->getBool("VisualizePoints"))
				ps->addPoint(points[v], color);

			if(pars->getBool("VisualizeVectors"))
			{
				// Draw principle directions
				Vector3 v1 = Vector3(pd1[v].normalized() * avgEdge * 0.5 * pars->getFloat("VisualizeScale"));
				Vector3 v2 = Vector3(pd2[v].normalized() * avgEdge * 0.5 * pars->getFloat("VisualizeScale"));

				vs1->addVector(Vector3(points[v] - (v1 * 0.5)), v1);
				vs2->addVector(Vector3(points[v] - (v2 * 0.5)), v2);
			}
		}

		drawArea()->clear();

		drawArea()->addRenderObject(vs1);
		drawArea()->addRenderObject(vs2);
		drawArea()->addRenderObject(ps);
	}
}
