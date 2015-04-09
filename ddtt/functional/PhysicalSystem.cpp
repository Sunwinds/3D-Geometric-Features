#include "PhysicalSystem.h"
#include "PhysicsHelper.h"

#include "RenderObjectExt.h"
using namespace starlab;

#include "lemon_graphviz.h"
#include <lemon/dijkstra.h>
using namespace lemon;

Q_DECLARE_METATYPE(Vector3)
Q_DECLARE_METATYPE(PropertyMap)
Q_DECLARE_METATYPE(ListGraph::Node)

PropertyMap physicalSystem(Structure::Graph * graph, QVector<RenderObject::Base*> & debug)
{
	PropertyMap result;

	/// Find set of ground nodes
	std::vector<Structure::Node*> floorNodes;
	{
		Eigen::AlignedBox3d box = graph->bbox();
		double bottom = box.min().z();
		double height = box.max().z() - bottom;

		double threshold = height * 0.05;

		for(Structure::Node * n : graph->nodes){
			double minz = n->bbox().min().z();
			double dist = abs(minz - bottom);
			if( dist > threshold ) continue;
			floorNodes.push_back(n);
		}
	}

	// Compute shortest paths
	ListDigraph order_graph;
	ListDigraph::NodeMap<Structure::Node*> orderNode(order_graph);
	ListDigraph::NodeMap<double> totalForce(order_graph);
	{
		ListGraph g;
		ListGraph::EdgeMap<double> length(g);
		ListGraph::NodeMap<double> dist(g);
		ListGraph::NodeMap<ListDigraph::Node> nodeMap(g);
		ListDigraph::NodeMap<QString> node_names(order_graph);

		// Nodes
		for(Structure::Node * n : graph->nodes){
			ListGraph::Node node = g.addNode();
			n->property["node"].setValue( node );
			
			ListDigraph::Node ordered_node = order_graph.addNode();
			nodeMap[node] = ordered_node;
			node_names[ordered_node] = n->id;
			orderNode[ordered_node] = n;
			totalForce[ordered_node] = 0.0;
		}

		// Edges
		for(Structure::Link * l : graph->edges){
			ListGraph::Edge a = g.addEdge( l->n1->property["node"].value<ListGraph::Node>(), l->n2->property["node"].value<ListGraph::Node>() );
			length[a] = 1.0;
		}

		// Prepare path finder
		Dijkstra< ListGraph, ListGraph::EdgeMap<double> > dijkstra(g, length);
		dijkstra.distMap(dist);
		dijkstra.init();

		// Start nodes
		for(Structure::Node * n : floorNodes) dijkstra.addSource( n->property["node"].value<ListGraph::Node>() );

		// Compute paths
		dijkstra.start();

		// Store nodes as layers based on distance
		QMap< int, std::vector<ListGraph::Node> > layers;
		for (ListGraph::NodeIt n(g); n != INVALID; ++n){
			layers[int(dist[n])].push_back(n);
		}

		// Create directed graph based on distance
		// Arcs only go from higher to lower
		size_t max_layer = layers.keys().back();
		for(size_t i = 0; i < max_layer; i++)
		{
			for(size_t j = 0; j < layers[i].size(); j++)
			{
				for(size_t k = 0; k < layers[i+1].size(); k++)
				{
					ListGraph::Node u = layers[i][j];
					ListGraph::Node v = layers[i+1][k];

					if( findEdge(g, u, v) == INVALID ) continue;

					order_graph.addArc( nodeMap[v], nodeMap[u] );
				}
			}
		}

		saveToGraphVizFile(digraphToGraphViz(order_graph, node_names));
	}

	// Propagate forces
	InDegMap<ListDigraph> indeg(order_graph);

	while( countNodes(order_graph) )
	{
		std::vector<ListDigraph::Node> sources;
		for (ListDigraph::NodeIt n(order_graph); n != INVALID; ++n) 
			if (indeg[n] == 0) sources.push_back(n);

		for (ListDigraph::Node n : sources)
		{
			Structure::Node * node = orderNode[n];

			double myMass = PhysicsHelper(node->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >().data()).volume();

			totalForce[n] += myMass;

			// Distribute force
			double eachForce = totalForce[n] / countOutArcs(order_graph, n);

			for (ListDigraph::OutArcIt a(order_graph, n); a != INVALID; ++a)
			{
				totalForce[ order_graph.target(a) ] += eachForce;
			}

			node->property["totalForce"] = totalForce[n];
			order_graph.erase(n);
		}
	}

	LineSegments * ls = new LineSegments (2);
	double r = graph->bbox().sizes().minCoeff() * 0.1;

	for(Structure::Node * n : graph->nodes)
	{
		if(!n->property.contains("mesh")) continue;
		QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();

		// Get physical properties
		Eigen::Matrix3d J;
		double mass;
		Vector3 com = PhysicsHelper(nodeMesh.data()).centerOfMass(J, 1.0, &mass);

		// Visualize
		//spheres->addSphere(com, r * 0.5);
		//ls->addLine(com, Vector3(com + (Vector3(0,0,-1) * r * 1)));

		// Store results
		PropertyMap nprop;
		nprop["com"].setValue(com);
		result[n->id] = nprop;

		n->property["mass"] = mass;
	}

	/// Assign based on forces toward ground
	{
		QMap<Structure::Node*, double> forces;

		foreach(Structure::Node * n, graph->nodes)	
			forces[n] = n->property["totalForce"].toDouble();

		QList<double> forces_sorted = forces.values();
		qSort(forces_sorted);

		foreach(Structure::Node * n, graph->nodes){
			std::vector<Vector3> curve = n->discretizedAsCurve( 0.01 );
			ls->addLines(curve, starlab::qtJetColor(forces[n], forces_sorted.front(), forces_sorted.back()));
		}
	}

	debug.push_back(ls);

	return result;
}
