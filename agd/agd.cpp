#include "agd.h"

void agd::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichBool("IsSubset", true, "Subset of vertices"));
	pars->addParam(new RichInt("NumSubset", 30, "Count"));
	pars->addParam(new RichBool("VisualizeSubset", true, "Visualize subset"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

void agd::applyFilter(RichParameterSet *pars)
{
	// Clear any past visualizations
	drawArea()->deleteAllRenderObjects();

    FilterPlugin * gd_plugin = pluginManager()->getFilter("Geodesic distance: heat kernel");

	BoolVertexProperty src_points = mesh()->vertex_property<bool>("v:geo_src_points", false);
	ScalarVertexProperty geo_dist = mesh()->vertex_property<Scalar>("v:uniformDistance", 0);

	// Consider a subset of vertices
	if(pars && pars->getBool("IsSubset"))
	{
		// Starting point
		src_points[Vertex(0)] = true;

		int numIter = pars->getInt("NumSubset");

		for(int itr = 0; itr < numIter; itr++)
		{
            gd_plugin->applyFilter( NULL );

			// Find furthest point, add it to source set
			double maxDist = -DBL_MAX;
			int maxIDX = -1;
			foreach(Vertex v, mesh()->vertices())
			{
				if(geo_dist[v] > maxDist){
					maxDist = geo_dist[v];
					maxIDX = v.idx();
				}
			}

			src_points[Vertex(maxIDX)] = true;

			mainWindow()->setStatusBarMessage( QString("iter (%1) out of (%2)").arg(itr).arg(numIter)  );

			qApp->processEvents();
		}

		// Visualize selected subset
		if(pars->getBool("VisualizeSubset"))
		{
			Vector3VertexProperty points = mesh()->vertex_coordinates();
			foreach(Vertex v, mesh()->vertices()){
				if(src_points[v]) drawArea()->drawPoint( points[v], 16.0);
			}
		}
	}
	else
	{
		foreach(Vertex v, mesh()->vertices()) src_points[v] = true;
	}

	drawArea()->updateGL();
	qApp->processEvents();

	// Compute vertex-to-all distances
	ScalarVertexProperty all_dist = mesh()->vertex_property<Scalar>("v:agd", 0);

    RichParameterSet * gd_params = new RichParameterSet;
    gd_plugin->initParameters( gd_params );
    gd_params->setValue("visualizeDistance", false);

	foreach(Vertex v, mesh()->vertices())
	{
		if(!src_points[v]) continue;

        gd_params->setValue("sourceVertex", int( v.idx() ));
        gd_plugin->applyFilter( gd_params );
		ScalarVertexProperty geo_dist = mesh()->vertex_property<Scalar>("v:uniformDistance", 0);

		foreach(Vertex v, mesh()->vertices())
		{
			all_dist[v] += geo_dist[v];
		}
	}

	double minVal = DBL_MAX;
	double maxVal = -DBL_MAX;

	foreach(Vertex v, mesh()->vertices())
	{
		minVal = qMin( minVal, all_dist[v] );
		maxVal = qMax( maxVal, all_dist[v] );
	}

	foreach(Vertex v, mesh()->vertices())
	{
		all_dist[v] = (all_dist[v] - minVal) / (maxVal - minVal);
	}

	// Visualize
	if(pars && pars->getBool("Visualize"))
	{
		Vector3VertexProperty points = mesh()->vertex_coordinates();

		foreach(Vertex v, mesh()->vertices())
		{
			drawArea()->drawPoint( points[v], 8.0, qtJetColorMap(1.0 - all_dist[v]) );
		}
	}
}
