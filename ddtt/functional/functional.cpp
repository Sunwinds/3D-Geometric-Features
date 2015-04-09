#include <QFileDialog>
#include "functional.h"
#include "RenderObjectExt.h"

#include "SurfaceMeshHelper.h"
#include "SamplingHelper.h"
#include "psimpl.h"
#include "MinOBB.h"
#include "GraphDistance.h"
#include "area3D_Polygon.h"
#include "Octree.h"
#include "GenericGraph.h"
#include "PhysicalSystem.h"
#include "PhysicsHelper.h"

using namespace Structure;

double stdev(const std::vector<double> & v)
{
	double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
	double m = sum / v.size();

	double accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
		accum += (d - m) * (d - m);
	});

	return sqrt(accum / (v.size()-1));
}

Vector3 lowest_point( std::vector<Vector3> pnts )
{
	Vector3 result(DBL_MAX, DBL_MAX, DBL_MAX);
	for(Vector3 p : pnts) if(p.z() < result.z()) result = p;
	return result;
}

template<typename T>
static inline T sumvec3( const std::vector<T> & V ){
	T vector_sum(0,0,0);
	for(int i = 0; i < (int)V.size(); i++) vector_sum += V[i];
	return vector_sum;
}

inline double meshArea( SurfaceMeshModel * mesh ){
	double area = 0;
	ScalarFaceProperty fareas = SurfaceMeshHelper(mesh).computeFaceAreas();
	for(auto f : mesh->faces()) area += fareas[f];
	return area;
}

inline double meshVolume( SurfaceMeshModel * mesh ){
	double volume = 0;

	Vector3VertexProperty points = mesh->vertex_coordinates();

	for(auto f: mesh->faces())
	{
		std::vector<Vector3> p;
		for(auto v : mesh->vertices(f)) p.push_back(points[v]);
		double signedVolumeTri = p[0].dot(p[1].cross(p[2])) / 6.0;
		volume += signedVolumeTri;
	}

	return volume;
}

SurfaceMeshModel * singleMesh(Structure::Graph * graph){
	SurfaceMeshModel * m = new SurfaceMeshModel;
	SurfaceMeshModel::Face_property<Node*> nodes = m->face_property<Node*>("f:node", NULL);

	{
		int offset = 0;
		foreach(Node * n, graph->nodes)	{
			if(!n->property.contains("mesh")) continue;
			QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
			Vector3VertexProperty nodePoints = nodeMesh->vertex_coordinates();
			for(auto v : nodeMesh->vertices()) m->add_vertex(nodePoints[v]);
			for(auto f : nodeMesh->faces()){
				std::vector<SurfaceMeshModel::Vertex> verts;
				for(Vertex v : nodeMesh->vertices(f)) verts.push_back( Vertex(v.idx() + offset) );
				Face ff = m->add_face(verts);

				nodes[ff] = n;
			}
			offset += nodeMesh->n_vertices();
		}
		m->updateBoundingBox();
		m->update_face_normals();
		m->update_vertex_normals();

	}

	return m;
}

void functional::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

void functional::simplified_graph()
{
	Structure::Graph * graph = new Structure::Graph( QFileDialog::getOpenFileName(0, "Open Graph") );

	starlab::LineSegments * ls = new starlab::LineSegments( 3 );

	foreach(Node * n, graph->nodes)
	{
		std::vector<Vector3> pnts = n->controlPoints();

		std::vector<Vector3> simplified_curve;

		if(n->type() == Structure::CURVE){
			std::vector<Vector3> original_curve = n->discretizedAsCurve( 0.005 );

			std::vector<double> curve;
			std::vector<double> simple_curve;
			foreach(Vector3 p, original_curve) for(int i = 0 ; i < 3; i++) curve.push_back( p[i] );

			psimpl::simplify_douglas_peucker_n<3>( curve.begin(), curve.end(), 6, std::back_inserter(simple_curve) );

			for(size_t i = 0; i < simple_curve.size() / 3; i++)
				simplified_curve.push_back( Vector3(simple_curve[(i*3)+0], simple_curve[(i*3)+1], simple_curve[(i*3)+2]) );
		}

		if(n->type() == Structure::SHEET){
			Structure::Sheet * sheet = (Structure::Sheet*) n;

			simplified_curve.push_back( sheet->surface.mCtrlPoint.front().front() );
			simplified_curve.push_back( sheet->surface.mCtrlPoint.front().back() );
			simplified_curve.push_back( sheet->surface.mCtrlPoint.back().back() );
			simplified_curve.push_back( sheet->surface.mCtrlPoint.back().front() );
			simplified_curve.push_back( sheet->surface.mCtrlPoint.front().front() );
		}

		ls->addLines( simplified_curve );
	}

	drawArea()->addRenderObject( ls );
}

void functional::measure_direction(Structure::Graph * graph)
{
	double threshold = 0.7;

	foreach(Node * n, graph->nodes)
	{
		if(!n->property.contains("mesh")) continue;
		QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();

		Vector3FaceProperty fnormals = nodeMesh->face_normals();
		Vector3VertexProperty points = nodeMesh->vertex_coordinates();
		ScalarFaceProperty farea = SurfaceMeshHelper( nodeMesh.data() ).computeFaceAreas();

		std::vector<Vector3> directions;
		directions.push_back( Vector3(1,0,0) );
		directions.push_back( Vector3(0,1,0) );
		directions.push_back( Vector3(0,0,1) );

		std::vector<double> classes(4, 0.0);

		foreach(Face f, nodeMesh->faces()){
			std::vector<Vector3> pts;
			foreach(Vertex v, nodeMesh->vertices(f)) pts.push_back(points[v]);

			for(size_t i = 0; i < directions.size(); i++)
				if( abs(dot( fnormals[f], directions[i] )) > threshold ) classes[i] += farea[f];
		}

		if( n->type() == Structure::CURVE )
		{
			std::vector<Vector3> curve = n->discretizedAsCurve( 0.005 );
			Vector3 dir = (curve.back() - curve.front()).normalized();

			double vertical = abs(dot(dir, Vector3(0,0,1)));
			double horizontal = abs(dot(dir, Vector3(1,0,0)));
			double sideways = abs(dot(dir, Vector3(0,1,0)));

			if(classes[2] == qMax(classes[0],qMax(classes[1],classes[2])))
			{
				//classes = std::vector<double>(4, 0.0);
				//classes[2] = 1.0;
			}

			if( vertical > threshold )			
			{ 
				classes[1] = 0; classes[2] = 0; 
			}

			if( vertical < (1.0 - threshold) )	
			{ 
				classes[0] = 0; classes[2] = 0; 
			}

			if( sideways > threshold )
			{
				classes = std::vector<double>(4, 0.0);
				classes[3] = 1.0;
			}
		}

		int best_idx = (std::max_element(classes.begin(), classes.end()) - classes.begin());

		std::vector<Qt::GlobalColor> color;
		color.push_back(Qt::blue);
		color.push_back(Qt::green);
		color.push_back(Qt::red);
		color.push_back(Qt::yellow);

		foreach(Face f, nodeMesh->faces()){
			std::vector<Vector3> pts;
			foreach(Vertex v, nodeMesh->vertices(f)) pts.push_back(points[v]);
			drawArea()->drawTriangle(pts[0], pts[1], pts[2], color[best_idx]);
		}

		QStringList dnames;
		dnames << "vertical" << "horizantal" << "flat" << "horizantal2";

		n->meta["geo_direction"] = dnames[best_idx];
	}
}

void functional::measure_ground(Structure::Graph * graph)
{
	Eigen::AlignedBox3d box = graph->bbox();
	double bottom = box.min().z();
	double height = box.max().z() - bottom;

	starlab::PointSoup * ps = new starlab::PointSoup;
		
	std::vector<Vector3> floorPoints;

	double threshold = height * 0.1;

	// Find set of ground nodes
	foreach(Node * n, graph->nodes)
	{
		double minz = n->bbox().min().z();
		double dist = abs(minz - bottom);

		std::vector<Vector3> curve = n->discretizedAsCurve( 0.01 );

		if( dist > threshold )
			continue;
		else
			floorPoints.push_back( lowest_point(curve) );
	}

	double DIST_RESOLUTION = 0.03;

	GraphDistance gd(graph);
	gd.computeDistances(floorPoints, DIST_RESOLUTION);
	
	double max_dist = *std::max_element(gd.min_distance.begin(), gd.min_distance.end());
	double min_dist = *std::min_element(gd.min_distance.begin(), gd.min_distance.end());

	for(size_t i = 0; i < gd.allPoints.size(); i++)
	{
		double val = (gd.min_distance[i] - min_dist) / (max_dist - min_dist);
		ps->addPoint( gd.allPoints[i], starlab::qtJetColor(val) );
	}

	drawArea()->addRenderObject( ps );
}

void functional::measure_ground_parts(Structure::Graph * graph)
{
	/// Find set of ground nodes
	std::vector<Structure::Node*> floorNodes;
	{
		Eigen::AlignedBox3d box = graph->bbox();
		double bottom = box.min().z();
		double height = box.max().z() - bottom;

		double threshold = height * 0.1;

		foreach(Node * n, graph->nodes)
		{
			double minz = n->bbox().min().z();
			double dist = abs(minz - bottom);

			std::vector<Vector3> curve = n->discretizedAsCurve( 0.01 );

			if( dist > threshold )
				continue;
			else
				floorNodes.push_back(n);
		}
	}

	/// Build a graph
	GenericGraphs::Graph<size_t, double> g;
	int floorNode = graph->nodes.size() * 10;
	{
		// Add regular edges
		for(Link * e : graph->edges){
			g.AddEdge(e->n1->property["index"].toUInt(), e->n2->property["index"].toUInt(), 1, e->property["uid"].toInt());
		}

		for(Node * n : floorNodes){
			g.AddEdge(floorNode, n->property["index"].toUInt(), 0);
		}
	}

	starlab::LineSegments * ls = new starlab::LineSegments( 3 );

	/// Assign based on distance to ground
	{
		QMap<Node*, double> dists;

		foreach(Node * n, graph->nodes)	
			dists[n] = g.NodeDistance(floorNode, n->property["index"].toUInt());

		QList<double> dists_sorted = dists.values();
		qSort(dists_sorted);

		foreach(Node * n, graph->nodes){
			std::vector<Vector3> curve = n->discretizedAsCurve( 0.01 );

			double val = (dists[n] - dists_sorted.front()) / (dists_sorted.back() - dists_sorted.front());
			ls->addLines( curve, starlab::qtJetColor(val) );

			n->meta["structure_ground"] = QString::number(val);
		}
	}

	drawArea()->addRenderObject( ls );
}

void functional::measure_height(Structure::Graph * graph)
{
	starlab::PointSoup * ps = new starlab::PointSoup;

	Eigen::AlignedBox3d box = graph->bbox();
	double min_z = box.min().z(), max_z = box.max().z();

	foreach(Node * n, graph->nodes){
		std::vector<Vector3> curve = n->discretizedAsCurve( 0.01 );
		for(auto p : curve){
			ps->addPoint( p, starlab::qtJetColor(p.z(), min_z, max_z) );

		}

		double val = (n->bbox().center().z() - min_z) / (max_z - min_z);
		n->meta["geo_height"] = QString::number(val);
	}

	drawArea()->addRenderObject( ps );
}

void functional::measure_curvyness(Structure::Graph * graph)
{
	starlab::LineSegments * ls = new starlab::LineSegments( 3 );

	Eigen::AlignedBox3d graph_bbox = graph->bbox();
	double res = graph_bbox.sizes().norm() * 0.002;

	foreach(Node * n, graph->nodes)
	{
		std::vector<Vector3> curve = n->discretizedAsCurve( res );

		Vector3 center = sumvec3(curve) / curve.size();

		if(curve.size() < 2) continue;
		int idx = (curve.size()-1) * 0.25;

		Vector3 a = (curve[idx] - center).normalized();
		Vector3 b = (curve[0] - center).normalized();
		Vector3 axis = a.cross(b).normalized();
		
		double sum_angle = signedAngle( (curve.front() - center).normalized(), (curve.back() - center).normalized(), axis );

		double val = abs(sum_angle) / M_PI;
		ls->addLines( curve, starlab::qtJetColor(val) );

		n->meta["geo_curvy"] = val;
	}

	drawArea()->addRenderObject(ls);
}

void functional::measure_area_volume_ratio(Structure::Graph * graph)
{
	QMap<Node*, double> vals;

	double totalArea = 0, totalVolume = 0;

	foreach(Node * n, graph->nodes)
	{
		if(!n->property.contains("mesh")) continue;
		QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();

		totalArea += meshArea(nodeMesh.data());
		totalVolume += meshVolume(nodeMesh.data());
	}

	foreach(Node * n, graph->nodes)
	{
		if(!n->property.contains("mesh")) continue;
		QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();

		double mesh_area = meshArea(nodeMesh.data());
		double mesh_volume = meshVolume(nodeMesh.data());

		vals[n] = mesh_area / mesh_volume;
	}

	QList<double> all_vals = vals.values();
	qSort(all_vals);

	foreach(Node * n, graph->nodes)
	{
		if(!n->property.contains("mesh")) continue;
		QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
		Vector3VertexProperty points = nodeMesh->vertex_coordinates();

		vals[n] = (vals[n] - all_vals.front()) / ( all_vals.back() - all_vals.front() );

		foreach(Face f, nodeMesh->faces()){
			std::vector<Vector3> pts;
			foreach(Vertex v, nodeMesh->vertices(f)) pts.push_back(points[v]);
			drawArea()->drawTriangle(pts[0], pts[1], pts[2], starlab::qtJetColor( vals[n] ));
		}

		n->meta["geo_areaVolume"] = QString::number(vals[n]);
	}
}

void functional::measure_access(Structure::Graph * graph)
{
	starlab::LineSegments * ls = new starlab::LineSegments;
	starlab::PointSoup * ps = new starlab::PointSoup;

	/// Combine into single mesh
	SurfaceMeshModel * m = singleMesh(graph);

	/// Collect samples
	QVector<SamplePoint> samples;
	QVector<Vector3> normals;
	double avgSpacing = 0;

	SamplingHelper::SimilarSampler::FaceSamples(m, 50000, normals, &avgSpacing, &samples);

	std::map<int,int> samplesPerFace;
	for(SamplePoint sp : samples) samplesPerFace[sp.findex]++;
	
	// Force minimum number of samples (for bad triangles)
	Vector3VertexProperty m_points = m->vertex_coordinates();
	Vector3FaceProperty fnormals = m->face_normals();
	for(Face f : m->faces()){
		std::vector<Vector3> verts;
		for(Vertex v : m->vertices(f)) verts.push_back(m_points[v]);

		int N = 8;
		if(samplesPerFace[f.idx()] >= N * N) continue;

		for(int i = 1; i + 1 < N; i++){
			for(int j = 1; j + 1 < N; j++){
				double w = (1 - (double(i) / (N-1)));
				
				double t = double(j) / (N-1);
				double u = (1 - w) * (t);
				double v = (1 - w) * (1 - t);

				Vector3 uvw(u, v, w);

				Vector3 p(0,0,0);
				for(int vi = 0; vi < 3; vi++) p = p + (verts[vi] * uvw[vi]);

				samples.push_back( SamplePoint( p, fnormals[f], f.idx(), uvw[0], uvw[1] ) );
			}
		}
	}

	/// Build octree
	Octree octree(m);

	Vector3 direction(0,0,1);
	double max_dist = m->bbox().sizes().z();

	SurfaceMeshModel::Face_property<Node*> nodes = m->face_property<Node*>("f:node", NULL);
	QMap<Node*,double> nodeAccess;

	/// Shoot rays
	double dot_threshold = 0.7;
	for(SamplePoint sp : samples)
	{
		if( dot(sp.n, direction) < dot_threshold ) continue;
		
		int faceIndex = -1;
		Vector3 isect = octree.closestIntersectionPoint(Ray(sp.pos + direction * (max_dist * 0.0001), direction), &faceIndex);
		if(faceIndex < 0) 
		{
			isect = sp.pos + direction * max_dist;
		}
		else
		{
			// Ignore back face intersections
			if(dot(fnormals[Face(faceIndex)], direction) > -dot_threshold) continue;
		}

		double dist = (isect - sp.pos).norm();

		ls->addLine(sp.pos, Vector3(sp.pos + (direction * (max_dist * 0.04))), starlab::qtJetColor(dist, 0, max_dist) );

		// Record result
		if( faceIndex < 0 ) continue;
		
		Node * n = nodes[Face(faceIndex)];
		nodeAccess[n] += 1;
	}

	double min_val = nodeAccess.values().front();
	double max_val = nodeAccess.values().back();

	for(Node * n : graph->nodes){
		double range = (max_val - min_val);
		double val = (range == 0) ? 0 : (nodeAccess[n] - min_val) / range;
		n->meta["geo_access"] = QString::number(val);
	}

	drawArea()->addRenderObject( ls );

	document()->addModel(m);

	drawArea()->setRenderer(m, "Wireframe");
	drawArea()->updateGL();
}

void functional::measure_graph_agd(Structure::Graph * graph)
{
	GenericGraphs::Graph<size_t, double> g;
	for(Link * e : graph->edges){
		g.AddEdge(e->n1->property["index"].toUInt(), e->n2->property["index"].toUInt(), 1, e->property["uid"].toInt());
	}

	QMap<Node*,double> node_dists;

	for(Node * ni : graph->nodes)
	{
		double sum_dist = 0;

		for(Node * nj : graph->nodes)
			sum_dist += g.NodeDistance(ni->property["index"].toUInt(), nj->property["index"].toUInt());

		node_dists[ni] = sum_dist;
	}

	starlab::LineSegments * ls = new starlab::LineSegments( 3 );

	QList<double> dists_sorted = node_dists.values();
	qSort(dists_sorted);

	foreach(Node * n, graph->nodes){
		std::vector<Vector3> curve = n->discretizedAsCurve( 0.01 );

		double val = (node_dists[n] - dists_sorted.front()) / (dists_sorted.back() - dists_sorted.front());
		ls->addLines( curve, starlab::qtJetColor(val) );

		n->meta["structure_agd"] = val;
	}

	drawArea()->addRenderObject( ls );
}

void functional::measure_physical_system(Structure::Graph * graph)
{
	QVector<RenderObject::Base*> debug;
	PropertyMap result = physicalSystem( graph, debug );

	QMap<Structure::Node*, double> forces;

	foreach(Structure::Node * n, graph->nodes)	
		forces[n] = n->property["totalForce"].toDouble();

	QList<double> forces_sorted = forces.values();
	qSort(forces_sorted);

	foreach(Structure::Node * n, graph->nodes){
		double val = (forces[n] - forces_sorted.front()) / (forces_sorted.back() - forces_sorted.front());
		n->meta["structure_physical"] = val;
	}

	drawArea()->addRenderObject(debug.front());
}

void functional::measure_position(Structure::Graph * graph)
{
	starlab::LineSegments * ls = new starlab::LineSegments;

	Eigen::AlignedBox3d bbox = graph->bbox();
	Vector3 bbox_min = bbox.min();
	Vector3 bbox_max = bbox.max();

	foreach(Node * n, graph->nodes)
	{
		Vector3 p = n->bbox().center();
		std::vector<Vector3> curve = n->discretizedAsCurve( 0.01 );

		double r = (p[0] - bbox_min[0]) / (bbox_max[0] - bbox_min[0]);
		double g = (p[1] - bbox_min[1]) / (bbox_max[1] - bbox_min[1]);
		double b = (p[2] - bbox_min[2]) / (bbox_max[2] - bbox_min[2]);

		QColor color = QColor::fromRgbF(r,g,b);
		ls->addLines(curve, color);

		n->meta["geo_x"] = r;
		n->meta["geo_y"] = g;
		n->meta["geo_z"] = b;
	}

	drawArea()->addRenderObject( ls );
}

void functional::applyFilter(RichParameterSet *pars)
{
	// Default view
	{
		drawArea()->camera()->setSceneCenter(qglviewer::Vec(0,0,0.5));
		drawArea()->camera()->setUpVector(qglviewer::Vec(0,0,1));
		drawArea()->camera()->setPosition(qglviewer::Vec(-1.5,-2,0.5));
		drawArea()->camera()->lookAt(qglviewer::Vec(0,0,0));
		drawArea()->setSceneRadius(1);
		drawArea()->showEntireScene();
		drawArea()->updateGL();
	}

	QString datasetPath = QFileDialog::getExistingDirectory(0, "Dataset", settings()->getString("lastUsedDirectory"));

	settings()->set("lastUsedDirectory", datasetPath);

	QDir datasetDir(datasetPath);
	QStringList subdirs = datasetDir.entryList(QDir::Dirs | QDir::NoSymLinks | QDir::NoDotAndDotDot);

	QVector< QMap<QString,QString> > graphMetas;

	foreach(QString subdir, subdirs)
	{
		// Special folders
		if(subdir == "corr") continue;
		QDir d(datasetPath + "/" + subdir);
		// Check if no graph is in this folder
		if( d.entryList(QStringList() << "*.xml", QDir::Files).isEmpty() ) continue;
		QString filename = d.absolutePath() + "/" + d.entryList(QStringList() << "*.xml", QDir::Files).join("");
		Structure::Graph * graph = new Structure::Graph( filename );

		document()->clear();

		/* Apply measure */
		measure_direction( graph );
		measure_ground( graph );
		measure_ground_parts( graph );
		//measure_height( graph );
		measure_curvyness( graph );
		measure_area_volume_ratio( graph );
		measure_access( graph );
		measure_graph_agd( graph );
		measure_physical_system( graph );
		measure_position( graph );

		bool isDrawMesh = true;
		if( isDrawMesh ){
			SurfaceMeshModel * m = singleMesh(graph);
			document()->addModel(m);
			drawArea()->setRenderer(m, "Transparent");
		}

		foreach(Node * n, graph->nodes)
		{
			QMap<QString,QString> nodeMeta;
			foreach(QString key, n->meta.keys()) 
				nodeMeta[key] = n->meta[key].toString();

			nodeMeta["id"] = n->id;
			nodeMeta["graph"] = QFileInfo(graph->property["name"].toString()).baseName();

			graphMetas.push_back(nodeMeta);
		}

		drawArea()->updateGL();
		qApp->processEvents();

		if(subdirs.size() > 1) 
		{
			QString thumb_filename = d.absolutePath() + "/../" + d.entryList(QStringList() << "*.png", QDir::Files).join("");
			drawArea()->saveSnapshot(thumb_filename);

			drawArea()->clear();
		}
	}

	QString filename = "dataset_functional.csv";
	QFile file(filename); file.open(QIODevice::WriteOnly | QIODevice::Text);
	QTextStream out(&file);
	
	typedef QMap<QString,QString> StringMap;
	QList<QString> keys = graphMetas.front().keys();
	
	QStringList allkeys;
	foreach(QString key, keys) 
	{
		allkeys << key;
	}
	out << allkeys.join(",") << "\n";

	foreach(StringMap meta, graphMetas){
		QStringList data;
		foreach(QString key, keys) 
		{
			QString value = meta[key];
			if(key == "label" && value.trimmed().isEmpty()){
				value = "none";
			}
			data << value;
		}
		out << data.join(",") << "\n";
	}

	file.close();
}
