#include <stack>
#include <omp.h>
#include <QApplication>
#include "analyze.h"

#include "MinOBB.h"

#include "StructureGraph.h"

std::vector< SurfaceMesh::SurfaceMeshModel* > connectedPieces(SurfaceMesh::SurfaceMeshModel * mesh) {
	std::vector< SurfaceMesh::SurfaceMeshModel* > pieces;
	std::vector< std::vector< SurfaceMesh::Face > > pieceFaces;
	Vector3VertexProperty points = mesh->vertex_coordinates();
	BoolFaceProperty fvisisted = mesh->face_property<bool>("f:visisted", false);
	ScalarFaceProperty fsegment = mesh->face_property<Scalar>("f:segment", -1);
	for(auto f : mesh->faces()){
		fvisisted[f] = false;
		fsegment[f] = -1;
	}
	std::vector<int> sizes;
	for(Face f : mesh->faces()){
		if( fvisisted[f] ) continue;
		pieceFaces.push_back( std::vector< SurfaceMesh::Face >() );
		std::stack<Face> unprocessed;
		unprocessed.push(f);
		while( !unprocessed.empty() ){
			Face curFace = unprocessed.top();
			unprocessed.pop();
			if(fvisisted[curFace]) continue;
			fvisisted[curFace] = true;
			for(auto v : mesh->vertices(curFace)){ 
				for(auto fj : mesh->faces(v)){
					if(!mesh->is_valid(fj) || fvisisted[fj]) continue;
					unprocessed.push(fj);
				}
			}
			fsegment[curFace] = (pieceFaces.size() - 1);
			pieceFaces.back().push_back( curFace );
		}

		SurfaceMesh::SurfaceMeshModel * piece = new SurfaceMesh::SurfaceMeshModel;
		piece->reserve( pieceFaces.back().size() * 3, pieceFaces.back().size() * 6, pieceFaces.back().size() );

		for(auto v : mesh->vertices()) 
			piece->add_vertex( points[v] );

		for(auto f : pieceFaces.back()){
			std::vector<Vertex> face;
			for(auto v : mesh->vertices(f)) face.push_back(v);
			piece->add_face(face);
		}

		for(auto v : piece->vertices()) if(piece->is_isolated(v)) piece->remove_vertex(v);
		piece->garbage_collection();
		sizes.push_back(piece->n_vertices());
		pieces.push_back( piece );
	}

	// Assign colors
	{
		std::set<int> segments;
		for(auto f : mesh->faces()) segments.insert(fsegment[f]);
		segments.erase(-1);
		std::vector< std::vector<double> > colors = starlab::randomColors(segments.size());
		QVariant colorsVar;
		colorsVar.setValue(colors);
		mesh->setProperty("colors", colorsVar);
	}

	return pieces;
}

void analyze::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichBool("IsParallel", false, "IsParallel"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
	pars->addParam(new RichBool("VisualizeBoxes", false, "VisualizeBoxes"));
	pars->addParam(new RichBool("VisualizeGraph", true, "VisualizeGraph"));
}

void analyze::applyFilter(RichParameterSet* pars)
{
	starlab::LineSegments * bbox_lines = new starlab::LineSegments(1);
	std::vector<starlab::LineSegments*> dir_lines(3, new starlab::LineSegments(2));
	starlab::PolygonSoup * ps = new starlab::PolygonSoup;
	std::vector< std::vector<double> > colors = starlab::randomColors(12);

	mainWindow()->setStatusBarMessage("Finding connected pieces..");
	QApplication::processEvents();
	std::vector< SurfaceMesh::SurfaceMeshModel* > pieces = connectedPieces( mesh() );
	std::vector<Vector3> all_points;

	std::vector<MinOBB::OBB> obb( pieces.size() );

	mainWindow()->setStatusBarMessage("Finding OBBs..");
	QApplication::processEvents();
	QElapsedTimer timer; timer.start();

	if( pars->getBool("IsParallel") )
	{
		int N = pieces.size();

		#pragma omp parallel for
		for(int i = 0; i < N; i++)
		{
			Vector3VertexProperty points = pieces[i]->vertex_coordinates();
			std::vector<Vector3> all_points;
			for(auto v: pieces[i]->vertices()) all_points.push_back( points[v] );

			obb[i] = MinOBB::OBB( all_points, true );
		}
	}
	else
	{
		for(size_t pi = 0; pi < pieces.size(); pi++)
		{
			auto piece = pieces[pi];

			Vector3VertexProperty points = piece->vertex_coordinates();
			std::vector<Vector3> all_points;
			for(auto v: piece->vertices()) all_points.push_back( points[v] );

			bool isAddJitter = true;
			obb[pi] = MinOBB::OBB( all_points, isAddJitter );

			std::vector<Vector3> dir(3);
			for(int i = 0; i < 3; i++){
				dir[0][i] = obb[pi].bb.dir_1[i];
				dir[1][i] = obb[pi].bb.dir_2[i];
				dir[2][i] = obb[pi].bb.dir_3[i];
			}

			std::vector<Vector3> c = obb[pi].corners<Vector3>();

			bbox_lines->addLines( (QVector<Vector3>() << c[0] << c[1] << c[5] << c[4] << c[0]).toStdVector(), Qt::black );
			bbox_lines->addLines( (QVector<Vector3>() << c[0] << c[1] << c[3] << c[2] << c[0]).toStdVector(), Qt::black );
			bbox_lines->addLines( (QVector<Vector3>() << c[6] << c[7] << c[3] << c[2] << c[6]).toStdVector(), Qt::black );
			bbox_lines->addLines( (QVector<Vector3>() << c[6] << c[7] << c[5] << c[4] << c[6]).toStdVector(), Qt::black );

			Vector3 mean(0,0,0);
			for(auto p: c) mean += p;
			mean /= c.size();

			dir_lines[0]->addLine(mean, Vector3(mean + dir[0] * 0.5 * obb[pi].bb.length(0)), Qt::red);
			dir_lines[1]->addLine(mean, Vector3(mean + dir[1] * 0.5 * obb[pi].bb.length(1)), Qt::green);
			dir_lines[2]->addLine(mean, Vector3(mean + dir[2] * 0.5 * obb[pi].bb.length(2)), Qt::blue);

			if( pars->getBool("VisualizeBoxes") ){
				std::vector< std::vector<Vector3> > face = obb[pi].faces<Vector3>();
				for(size_t i = 0; i < face.size(); i++){
					ps->addPoly((QVector<starlab::QVector3>() << face[i][2] << face[i][1] << face[i][0]), QColor::fromRgbF(colors[i][0],colors[i][1],colors[i][2], 0.1));
				}
			}
		}
	}

	if( pars->getBool("Visualize") )
	{
		drawArea()->addRenderObject(bbox_lines);
		for(size_t i = 0; i < dir_lines.size(); i++) drawArea()->addRenderObject( dir_lines[i] );
		drawArea()->addRenderObject(ps);
	}

	mainWindow()->setStatusBarMessage( QString("Done finding OBBs. (%1 ms)").arg(timer.elapsed()));
	QApplication::processEvents();

	Structure::Graph * graph = new Structure::Graph();

	for(size_t bi = 0; bi < obb.size(); bi++)
	{
		MinOBB::OBB::OBBOX<Vector3> oi( obb[bi].bb );
		
		std::vector<Vector3> dir(3);
		for(int i = 0; i < 3; i++){
			dir[0][i] = obb[bi].bb.dir_1[i];
			dir[1][i] = obb[bi].bb.dir_2[i];
			dir[2][i] = obb[bi].bb.dir_3[i];
		}

		Vector3 delta = dir[2] * oi.HalfSize[2];

		Vector3 from = oi.Center + delta;
		Vector3 to = oi.Center - delta;

		NURBS::NURBSCurved curve = NURBS::NURBSCurved::createCurve(from, to);
		Structure::Curve * n = new Structure::Curve(curve, QString("node%1").arg(bi));

		graph->addNode( n );
	}

	mainWindow()->setStatusBarMessage("Now finding edges..");
	QApplication::processEvents();

	starlab::LineSegments * graphLines = new starlab::LineSegments(4);
	for(size_t i = 0; i < obb.size(); i++){
		MinOBB::OBB::OBBOX<Vector3> oi( obb[i].bb );
		for(size_t j = i + 1; j < obb.size(); j++){
			MinOBB::OBB::OBBOX<Vector3> oj( obb[j].bb );

			bool isIsect = MinOBB::OBB::testObbObb<Vector3>( oi, oj );
			if( isIsect ) 
			{
				graph->addEdge(QString("node%1").arg(i), QString("node%1").arg(j));

				if( pars->getBool("VisualizeGraph") )
					graphLines->addLine( oi.Center, oj.Center, Qt::magenta );
			}
		}
	}

	drawArea()->addRenderObject( graphLines );

	graph->saveToFile("graph.xml");
}
