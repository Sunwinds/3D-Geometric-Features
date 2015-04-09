#include "repair.h"
#include "SurfaceMeshHelper.h"

void repair::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichFloat("MaxAspectRatio", 20, "Max aspect ratio"));
	pars->addParam(new RichBool("CloseHoles", true, "Close Holes"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

Scalar triangleAspectRatio( const Vector3& _v0, const Vector3& _v1, const Vector3& _v2 ){
	Vector3 d0 = _v0 - _v1;
	Vector3 d1 = _v1 - _v2;

	// finds the max squared edge length
	Scalar l2, maxl2 = d0.squaredNorm();
	if ((l2=d1.squaredNorm()) > maxl2)
		maxl2 = l2;
	// keep searching for the max squared edge length
	d1 = _v2 - _v0;
	if ((l2=d1.squaredNorm()) > maxl2)
		maxl2 = l2;

	// squared area of the parallelogram spanned by d0 and d1
	Scalar a2 = cross(d0, d1).squaredNorm();

	// the area of the triangle would be
	// sqrt(a2)/2 or length * height / 2
	// aspect ratio = length / height
	//              = length * length / (2*area)
	//              = length * length / sqrt(a2)

	// returns the length of the longest edge
	//         divided by its corresponding height
	return sqrt( (maxl2 * maxl2) / a2 );
}

void repair::applyFilter(RichParameterSet *pars)
{
	Vector3VertexProperty points = mesh()->vertex_coordinates();

	if( pars->getBool("CloseHoles") )
	{
		Surface_mesh::Halfedge_property<bool> hvisisted = mesh()->halfedge_property<bool>("h:visisted", false);

		std::vector< std::vector<Halfedge> > holes;

        foreach(Halfedge h, mesh()->halfedges()){
			if( !mesh()->is_boundary(h) || hvisisted[h] ) continue;

			std::vector<Halfedge> hole;

			Halfedge hinit = h;
			h = mesh()->next_halfedge(h);

			while( h != hinit ) {
				hole.push_back( h );
				hvisisted[h] = true;
				h = mesh()->next_halfedge(h);
			}
			hole.push_back( h );

			holes.push_back(hole);
		}

		Vector3VertexProperty points = mesh()->vertex_coordinates();

        foreach(std::vector<Halfedge> hole, holes){
			std::vector<Vertex> holeVerts;
            foreach(auto h, hole){
				Vertex v = mesh()->to_vertex(h);
				holeVerts.push_back(v);
			}

			// Compute hole center
			Vector3 center(0,0,0);
            foreach(auto v, holeVerts) center += points[v];
			center /= holeVerts.size();

			Vertex c = mesh()->add_vertex(center);

			for(size_t vi = 0; vi < holeVerts.size(); vi++){
				std::vector<Vertex> face;

				face.push_back( holeVerts[vi] );
				face.push_back( holeVerts[(vi+1) % holeVerts.size()] );
				face.push_back( c );

				mesh()->add_face(face);
			}
		}
	}

	ScalarFaceProperty faspect = mesh()->face_property("f:aspect", 1.0);

	// Remove degenerate faces
	bool isGoodTriangle = false;

	while( !isGoodTriangle )
	{
		isGoodTriangle = true;

        foreach(Face f, mesh()->faces())
		{
			std::vector<Vector3> p;
            foreach(Vertex v, mesh()->vertices(f)) p.push_back(points[v]);

			double aspectRatio = triangleAspectRatio(p[0], p[1], p[2]);
			faspect[f] = aspectRatio;

			if(aspectRatio > pars->getFloat("MaxAspectRatio"))
			{
				Halfedge h = mesh()->halfedge(f);

				QMap<Halfedge,double> edgeLengths;
				Surface_mesh::Halfedge_around_face_circulator adjE(mesh(), f), eend = adjE;
				do { 
					Edge e = mesh()->edge(adjE);
					Vector3 p1 = points[mesh()->vertex(e,0)];
					Vector3 p2 = points[mesh()->vertex(e,1)];
					edgeLengths[adjE] = (p2-p1).norm();
				} while(++adjE != eend);

				bool isCollapsed = false;
				foreach(Halfedge e, edgeLengths.keys()){
					if(mesh()->is_collapse_ok(e)){
						mesh()->collapse(e);
						isCollapsed = true;
						break;
					}
				}
				
				if(!isCollapsed){
					qDebug() << "[edge collapse] something went wrong.";
					return;
				}

				isGoodTriangle = false;

				break;
			}
		}
	}

	// Clean up
	mesh()->remove_face_property(faspect);
    foreach(Vertex v, mesh()->vertices()) if(mesh()->is_isolated(v)) mesh()->remove_vertex(v);
	mesh()->garbage_collection();

	// Update properties
	mesh()->update_face_normals();
	mesh()->update_vertex_normals();
	mesh()->updateBoundingBox();
}
