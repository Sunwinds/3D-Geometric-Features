#include "concavity.h"
#include "RenderObjectExt.h"

#include "SurfaceMeshHelper.h"
#include "CurvatureEstimationHelper.h"
#include "StatisticsHelper.h"
#include "ColorizeHelper.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/CholmodSupport>

typedef Eigen::SimplicialLLT< Eigen::SparseMatrix<double> > SparseSolver;
//typedef Eigen::CholmodSupernodalLLT< Eigen::SparseMatrix<double> > SparseSolver;

void concavity::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichFloat("ConcaveWeight", 0.01f, "Concave weight"));
	pars->addParam(new RichInt("Iterations", 10, "Iterations"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

ScalarVertexProperty computeHarmonicField( std::pair<int,int> vertex_pair, SurfaceMeshModel * m_mesh, 
						  const ScalarVertexProperty & concaveVertices, 
						  double avgEdgeLength, double ConcaveWeight)
{
	// Vertex properties
	Vector3VertexProperty points = m_mesh->vertex_coordinates();
	Vector3VertexProperty normals = m_mesh->vertex_normals();

	// Field
	ScalarVertexProperty harmonicField = m_mesh->vertex_property<Scalar>("v:harmonicField", 0.0);

	// Efficient sparse matrix construction
	typedef Eigen::Triplet<double> T;
	std::vector< T > A;

	for(Vertex v : m_mesh->vertices())
	{
		double tot = 0;

		Vector3 v1 = points[v];
		Vector3 n1 = normals[v];

		std::vector<double> wei(m_mesh->valence(v), 0);
		int count = 0;

		/// Part 1: change weight here and remember to change Part 2 of weight
		for(Halfedge h : m_mesh->onering_hedges(v))
		{
			Vertex vj = m_mesh->to_vertex(h);
			Vector3 v2 = points[vj];

			double w = 1.0;

			if (concaveVertices[vj] || concaveVertices[v])
			{
				w = ((v1-v2).norm() / avgEdgeLength) * ConcaveWeight;
			}

			wei[count++] = w;
			tot += w;
		}

		count = 0;
		for(Halfedge h : m_mesh->onering_hedges(v))
			A.push_back(T(v.idx(), m_mesh->to_vertex(h).idx(), wei[count++] / tot));
		A.push_back(T(v.idx(), v.idx(), -1.0));
	}

	int n = m_mesh->n_vertices();
	int m = n + 2;

	// add positional weights
	int I = vertex_pair.first;
	int J = vertex_pair.second;

	double ConstantWeight = 1000;
	A.push_back(T(n, I, 1.0 * ConstantWeight));
	A.push_back(T(n + 1, J, 1.0 * ConstantWeight));

	Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
	b[I] = 1.0 * ConstantWeight * ConstantWeight;

	// Build sparse matrix
	Eigen::SparseMatrix<Scalar> L(m, n);
	L.setFromTriplets(A.begin(), A.end());
	Eigen::SparseMatrix<Scalar> ATA = L.transpose() * L;

	// Solve
	SparseSolver sparseSolver;
	sparseSolver.compute(ATA); 
	Eigen::VectorXd x = sparseSolver.solve(b);

	// Store solution
	for(Vertex v : m_mesh->vertices())
		harmonicField[v] = x[v.idx()];

	return harmonicField;
}

void concavity::applyFilter(RichParameterSet *pars)
{
	drawArea()->clear();

	// Vertex properties
	Vector3VertexProperty points = mesh()->vertex_coordinates();
	Vector3VertexProperty normals = mesh()->vertex_normals();

	// Find concave vertices
	ScalarVertexProperty concaveVertices = mesh()->vertex_property<Scalar>("v:concave", 0);
	for(Vertex vi : mesh()->vertices())
	{
		Vector3 v = points[vi];
		Vector3 normal = normals[vi];
		Vector3 center(0,0,0);

		for(Halfedge h : mesh()->onering_hedges(vi))
		{
			Vertex vj = mesh()->to_vertex(h);
			center += points[vj];
		}

		center /= mesh()->valence(vi);

		double conv = (v - center).dot(normal);

		if (conv < 0) concaveVertices[vi] = 1.0;
	}

	// Edge properties
	ScalarEdgeProperty elengths = SurfaceMeshHelper(mesh()).computeEdgeLengths();
	double avgEdgeLength = StatisticsHelper<Scalar>(mesh()).mean(ELENGTH);

	// Parameters
	double ConcaveWeight = pars->getFloat("ConcaveWeight");

	// Get 'X' most distant points
	{
		FilterPlugin * gd_plugin = pluginManager()->getFilter("Geodesic distance: heat kernel");
		BoolVertexProperty src_points = mesh()->vertex_property<bool>("v:geo_src_points", false);
		ScalarVertexProperty geo_dist = mesh()->vertex_property<Scalar>("v:uniformDistance", 0);

		// Consider a subset of vertices
		src_points[Vertex(0)] = true;
		int numIter = pars->getInt("Iterations");

		for(int itr = 0; itr < numIter; itr++){
			gd_plugin->applyFilter( NULL );

			// Find furthest point, add it to source set
			double maxDist = -DBL_MAX;
			int maxIDX = -1;
			foreach(Vertex v, mesh()->vertices()){
				if(geo_dist[v] > maxDist){
					maxDist = geo_dist[v];
					maxIDX = v.idx();
				}
			}
			src_points[Vertex(maxIDX)] = true;
		}
		
		// Collect the used ones
		std::vector<Vertex> used;
		for(Vertex v : mesh()->vertices()) 
		{
			if(src_points[v]) 
			{
				used.push_back(v);
				drawArea()->drawPoint(points[v],5);
			}
		}
		
		// Compute pair-wise
		ScalarVertexProperty sum = mesh()->vertex_property<Scalar>("v:sumHarmonicField", 0.0);
		for(Vertex v : mesh()->vertices()) sum[v] = 0.0; // clear sum

		for(size_t i = 0; i < used.size(); i++){
			for(size_t j = i + 1; j < used.size(); j++)
			{
				std::pair<int,int> vertex_pair(used[i].idx(), used[j].idx());
				ScalarVertexProperty harmonicField = computeHarmonicField(vertex_pair, mesh(), concaveVertices, avgEdgeLength, ConcaveWeight);

				for(Vertex v : mesh()->vertices()) sum[v] += harmonicField[v];
			}
		}
	}

	// Visualize
	if(pars->getBool("Visualize"))
	{	
		//ScalarVertexProperty field = mesh()->vertex_property<Scalar>("v:harmonicField", 0.0);
		ScalarVertexProperty field = mesh()->vertex_property<Scalar>("v:sumHarmonicField", 0.0);

		// Find limits
		double minVal = DBL_MAX, maxVal = -DBL_MAX;
		for(Vertex v : mesh()->vertices()){
			minVal = qMin(minVal, field[v]);
			maxVal = qMax(maxVal, field[v]);
		}

		// Assign colors
		Vector3VertexProperty vcolor = mesh()->vertex_property<Vector3>(VCOLOR);
		for(Vertex v : mesh()->vertices()){
			QColor color = starlab::qtJetColor( field[v], minVal, maxVal );
			vcolor[v] = Vector3(color.redF(), color.greenF(), color.blueF());
		}
	}
}
