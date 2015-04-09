#include <QApplication> // For mouse icon changing

#include <QFile>
#include <QFileInfo>
#include <QDir>
#include <QElapsedTimer>
#include <QDomElement>
#include <QStack>
#include <QMatrix4x4>

#include "StructureGraph.h"
using namespace Structure;

#include "GraphDistance.h"

#include "Synthesizer.h"

#include "GraphEmbed.h"
#include "GraphDraw2D.h"
#include "GenericGraph.h"

#include "QuickMeshDraw.h"

#include "Task.h"

Q_DECLARE_METATYPE( Eigen::AlignedBox3d )

Graph::Graph()
{
	init();
}

Graph::Graph( QString fileName )
{
	init();

	loadFromFile( fileName );
	property["name"] = fileName;
}

void Graph::init()
{
	property["showNodes"]	= true;
	property["showNames"]	= false;
	property["embeded2D"]	= false;
	property["showAABB"]	= false;
	property["showMeshes"]	= true;
	property["showEdges"]	= false;
	property["showTasks"]	= false;
	property["showCtrlPts"] = false;
	property["isSplatsHQ"]	= false;

	property["isBusy"] = false;

    vs = starlab::VectorSoup(Qt::blue);
    vs2 = starlab::VectorSoup(Qt::red);
    vs3 = starlab::VectorSoup(Qt::yellow);

    spheres = starlab::SphereSoup(Qt::yellow);
    spheres2 = starlab::SphereSoup(Qt::blue);

	ueid = 1001;
}

Graph::Graph( const Graph & other )
{
	foreach(Node * n, other.nodes)
	{
		this->addNode( n->clone() );
	}

	foreach(Link * e, other.edges)
	{
		Node * n1 = this->getNode(e->n1->id);
		Node * n2 = this->getNode(e->n2->id);
		
		Structure::Link * newEdge = this->addEdge(n1, n2, e->coord[0], e->coord[1], e->id);
		newEdge->property = e->property;
	}

	groups = other.groups;
	property = other.property;
	misc = other.misc;

	debugPoints = other.debugPoints;
	debugPoints2 = other.debugPoints2;
	debugPoints3 = other.debugPoints3;

	vs = other.vs; vs2 = other.vs2; vs3 = other.vs3;
	ps = other.ps; ps2 = other.ps2; ps3 = other.ps3;
	spheres = other.spheres; spheres2 = other.spheres2;

	ueid = other.ueid;
}

Graph::~Graph()
{
    qDeleteAll( nodes );
	nodes.clear();

	qDeleteAll( edges );
	edges.clear();

    //nodes.clear();
}

Eigen::AlignedBox3d Graph::bbox(bool isSkipUnready)
{
    Eigen::AlignedBox3d box;

    foreach(Node * n, nodes)
	{
		if( isSkipUnready )
		{
			if(getEdges(n->id).isEmpty()) continue;
			
			Task * t = n->property["task"].value<Task*>();
			if( t )
			{
				if(t->type == Task::GROW && !t->isReady) continue;
				if(t->type == Task::SHRINK && !t->isDone) continue;
			}
		}

        box = box.merged( n->bbox() );
    }

	//box.transform(QMatrix4x4() * 1.25);

    return box;
}

Node *Graph::addNode(Node * n)
{
    assert(getNode( n->id ) == NULL);

    nodes.push_back(n);

	// Add property : id
	n->property["index"] = nodes.size() - 1;
    n->property["groupID"] = -1;
    n->property["original_ID"] = n->id;

    return n;
}

Link * Graph::addEdge(QString n1_id, QString n2_id)
{
	return addEdge(getNode(n1_id), getNode(n2_id));
}

Link * Graph::addEdge(Node *n1, Node *n2)
{
	Vector3 intersectPoint = nodeIntersection(n1, n2);

	Array1D_Vector4d c1,c2;

	QString edgeType = POINT_EDGE;

	if(n1->type() == SHEET && n2->type() == SHEET)
	{
		Sheet *s1 = (Sheet*) n1, *s2 = (Sheet*) n2;
		double length = qMin(s1->bbox().diagonal().norm(), s2->bbox().diagonal().norm());
		std::vector<Vector3d> pnts = s1->surface.intersect(s2->surface, 0.005 * length, c1, c2);
		
		// DEBUG:
		//foreach(Vector3d p, pnts) debugPoints3.push_back(p);

		// Fall back
		double scaleFactor = 1.1;
		while(c1.size() < 1)
		{
			Vector3d closetPoint1 = s1->position(s1->approxCoordinates(intersectPoint));
			Vector3d closetPoint2 = s2->position(s2->approxCoordinates(intersectPoint));
			
			Vector3d delta = closetPoint2 - closetPoint1;
			Vector3d moveDelta = (scaleFactor * delta.norm()) * delta.normalized();
			
			s1->moveBy(moveDelta);

			std::vector<Vector3d> pnts = s1->surface.intersect(s2->surface, 0.005 * length, c1, c2);

			// DEBUG:
			//foreach(Vector3d p, pnts) debugPoints3.push_back(p);

			s1->moveBy(-moveDelta);

			scaleFactor += 0.2;
		}

		edgeType = LINE_EDGE;
	}
	else
	{
		c1.push_back(n1->approxCoordinates(intersectPoint));
		c2.push_back(n2->approxCoordinates(intersectPoint));
	}

    return addEdge(n1,n2,c1,c2,linkName(n1, n2));
}

Link * Graph::addEdge(Node *n1, Node *n2, Array1D_Vector4d coord1, Array1D_Vector4d coord2, QString linkName)
{
	if(linkName.isEmpty()) linkName = this->linkName(n1,n2);

	QString edgeType = POINT_EDGE;
	if(n1->type() == SHEET && n2->type() == SHEET) edgeType = LINE_EDGE;

	Link * e = new Link( n1, n2, coord1, coord2, edgeType, linkName );
	edges.push_back( e );

	e->property["uid"] = ueid++;

	return e;
}

void Graph::removeEdge( int uid )
{
	int edge_idx = -1;

	for(int i = 0; i < (int)edges.size(); i++){
		Link * e = edges[i];
		if(e->property["uid"].toInt() == uid){
			edge_idx = i;
			break;
		}
	}

	if(edge_idx < 0) return;

	Node *n1 = edges[edge_idx]->n1, *n2 = edges[edge_idx]->n2;

	delete edges[edge_idx];
	edges[edge_idx] = NULL;
	edges.remove(edge_idx);

	qDebug() << QString("[%1]->removeEdge( %2, %3 )").arg(name()).arg(n1->id).arg(n2->id);
}

void Graph::removeEdge( Node * n1, Node * n2 )
{
	int edge_idx = -1;

	for(int i = 0; i < (int)edges.size(); i++)
	{
		Link * e = edges[i];

		if((e->n1->id == n1->id && e->n2->id == n2->id) || (e->n2->id == n1->id && e->n1->id == n2->id)){
			edge_idx = i;
			break;
		}
	}
	
	if(edge_idx < 0) return;

	removeEdge( edges[edge_idx]->property["uid"].toInt() );
}

void Graph::removeEdge( QString n1_id, QString n2_id )
{
	removeEdge( getNode(n1_id), getNode(n2_id) );
}

void Graph::removeEdges( QString nodeID )
{
	QVector<Link*> edgs = this->getEdges(nodeID);
	foreach(Link * l, edgs) 
		this->removeEdge( l->n1->id, l->n2->id );
}

int Graph::indexOfNode( Node * node )
{
	for(int i = 0; i < (int) nodes.size(); i++)
		if(nodes[i]->id == node->id) return i;
	return -1;
}

QString Graph::linkName(Node *n1, Node *n2)
{
    return QString("%1 : %2").arg(n1->id).arg(n2->id);
}

QString Graph::linkName( QString n1_id, QString n2_id )
{
	return linkName(getNode(n1_id),getNode(n2_id));
}

Node *Graph::getNode(QString nodeID)
{
    foreach(Node* n, nodes)
	{
        if(n->id == nodeID) 
			return n;
	}

    return NULL;
}

Link *Graph::getEdge(QString id1, QString id2)
{
	for(int i = 0; i < (int)edges.size(); i++)
	{
		Link * e = edges[i];

		QString nid1 = e->n1->id;
		QString nid2 = e->n2->id;

		if( (nid1 == id1 && nid2 == id2) || (nid1 == id2 && nid2 == id1))
			return e;
	}
	
	return NULL;
}

Link* Graph::getEdge( int edgeUID )
{
	for(int i = 0; i < (int)edges.size(); i++)
	{
		Link * e = edges[i];

		if( e->property["uid"].toInt() == edgeUID) return e;
	}

	return NULL;
}

QVector<Link*> Graph::getEdges( QVector<int> edgeUIDs )
{
	QVector<Link*> result;
	foreach(int eid, edgeUIDs){
		Link * l = getEdge(eid);
		if(l) result.push_back( l );
	}
	return result;
}

QVector<int> Structure::Graph::getEdgeIDs( QVector<Link*> forEdges )
{
	QVector<int> result;
	foreach(Link * e, forEdges)	result.push_back( e->property["uid"].toInt() );
	return result;
}

Eigen::AlignedBox3d Graph::cached_bbox()
{
	if(!property.contains("bbox"))
		property["bbox"].setValue( bbox() );

	return property["bbox"].value<Eigen::AlignedBox3d>();
}

void Graph::setPropertyAll( QString prop_name, QVariant value )
{
	foreach(Node* n, nodes) 
		n->property[prop_name].setValue(value);
}

void Graph::setVisPropertyAll( QString prop_name, QVariant value )
{
	foreach(Node* n, nodes) 
		n->vis_property[prop_name].setValue(value);
}

void Graph::setPropertyFor( QVector<QString> nodeIDs, QString prop_name, QVariant value )
{
	foreach(QString nodeID, nodeIDs) 
		(getNode(nodeID))->property[prop_name].setValue(value);
}

QVector<Node*> Graph::nodesWithProperty( QString propertyName )
{
	QVector<Node*> result;
	
	foreach(Node* n, nodes) {
		if( n->hasProperty(propertyName) ) 
			result.push_back(n);
	}

	return result;
}

QVector<Node*> Graph::nodesWithProperty( QString propertyName, QVariant value )
{
	QVector<Node*> result;
	QVector<Node*> nodesWithProp = nodesWithProperty(propertyName);
	
	foreach(Node* n, nodesWithProp){
        if( n->hasProperty(propertyName) && n->property.value(propertyName) == value )
			result.push_back(n);
	}

	return result;
}

QVector<Node*> Graph::path( Node * from, Node * to )
{
	QVector<Node*> result;

	// Add nodes to generic graph
	GenericGraphs::Graph<int,int> g;
	QMap<Node*,int> nodeIdxMap;
	QMap<int,Node*> idxNodeMap;
	foreach(Node * n, nodes){
		int idx = nodeIdxMap.size();
		idxNodeMap[idx] = n;
		nodeIdxMap[n] = idx;
		g.AddVertex(idx);
	}

	// Add edges
	foreach(Link * e, edges)
		g.AddEdge(nodeIdxMap[e->n1],nodeIdxMap[e->n2],1);

	// Get path
	std::list<int> shortestPath = g.DijkstraShortestPath(nodeIdxMap[from], nodeIdxMap[to]);

	foreach(int idx, shortestPath)
		result.push_back(idxNodeMap[idx]);

	return result;
}

bool Graph::shareEdge( Node * n1, Node * n2 )
{
	foreach( Node * adj, this->adjNodes(n1) )
		if( adj->id == n2->id ) return true;
	return false;
}

bool Graph::shareEdge( QString nid1, QString nid2 )
{
	return shareEdge(getNode(nid1), getNode(nid2));
}

void Graph::clearDebug()
{
	debugPoints.clear(); debugPoints2.clear(); debugPoints3.clear();
	vs.clear(); vs2.clear(); vs3.clear();
	ps.clear(); ps2.clear(); ps3.clear();
	spheres.clear(); spheres2.clear();
}

void Graph::draw( QGLViewer * drawArea )
{
	if( property["isBusy"].toBool() ) return;

	if( property["showTasks"].toBool() )
	{
		vs.draw();	vs2.draw(5); vs3.draw();
		ps.draw();	ps2.draw();	ps3.draw();
		spheres.draw(); spheres2.draw();
	}

    foreach(Node * n, nodes)
    {
		// Debug points for node
		{
			glDisable(GL_LIGHTING);
			glPointSize(10);	glColor3d(1,0,1);
			glBegin(GL_POINTS);
			foreach(Vector3 p, n->debugPoints) glVector3(p);
			glEnd();
			glPointSize(12);	glColor3d(0,1,1);
			glBegin(GL_POINTS);
			foreach(Vector3 p, n->debugPoints2) glVector3(p);
			glEnd();
			glPointSize(14);	glColor3d(1,1,0);
			glBegin(GL_POINTS);
			foreach(Vector3 p, n->debugPoints3) glVector3(p);
			glEnd();
			glEnable(GL_LIGHTING);
		}

		if (n->property.contains("isReady") && !n->property["isReady"].toBool())
			continue;

		if ( n->property["shrunk"].toBool() || n->property["zeroGeometry"].toBool() )
			continue;

		if (property["showCurveFrames"].toBool())
		{
			if (n->property.contains("rmf_frames"))
			{
				std::vector<RMF::Frame> frames = n->property["rmf_frames"].value< std::vector<RMF::Frame> >();
                starlab::FrameSoup fs(0.05f);
				foreach (RMF::Frame f, frames)
					fs.addFrame(f.r, f.s, f.t, f.center);
				fs.draw();
			}
		}

		if( !property.contains("reconMesh") )
		{
			// Draw node mesh
			if( property["showMeshes"].toBool() )
			{
				if(n->property.contains("mesh"))
				{
					QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();

					QColor meshColor(255,255,255,8);
					if(!n->vis_property.contains("meshColor"))
						n->vis_property["meshColor"] = meshColor;
					else
						meshColor = n->vis_property["meshColor"].value<QColor>();

					bool meshSolid = n->vis_property["meshSolid"].toBool();

					if(!meshSolid) glDisable(GL_DEPTH_TEST);
					QuickMeshDraw::drawMeshSolid( nodeMesh.data(), meshColor );
					glEnable(GL_DEPTH_TEST);
				}
			}
		}

		// Task visualization:
		if( property["showTasks"].toBool() )
		{
			// Morph-time coloring
			if( n->property.contains("isActive")  )
			{
				if(n->property["isActive"].toBool())
					n->vis_property["color"] = QColor(150,150,255);
				else
					n->vis_property["color"] = QColor(150,255,150);
			}

			glDisable(GL_LIGHTING);

			if( n->property["isActive"].toBool() )
			{
                if( n->property.contains("consistentFrame") )
                {
                    RMF consistentFrame = n->property["consistentFrame"].value<RMF>();
                    starlab::FrameSoup fs(0.02f);
                    foreach(RMF::Frame f, consistentFrame.U) fs.addFrame(f.r, f.s, f.t, f.center);
                    fs.draw();
                }

				if( n->property.contains("rmf") )
				{
					RMF rmf = n->property["rmf"].value<RMF>();
                    starlab::FrameSoup fs(0.01f);
					foreach(RMF::Frame f, rmf.U) fs.addFrame(f.r, f.s, f.t, f.center);
					fs.draw();
				}

				if( n->property.contains("rmf2") )
				{
					RMF rmf = n->property["rmf2"].value<RMF>();
                    starlab::FrameSoup fs(0.01f);
					foreach(RMF::Frame f, rmf.U) fs.addFrame(f.r, f.s, f.t, f.center);
					fs.draw();
				}

				if( n->property.contains("frame") )
				{
					RMF::Frame f = n->property["frame"].value<RMF::Frame>();
                    starlab::FrameSoup fs(0.01f);
					fs.addFrame(f.r, f.s, f.t, Vector3(n->bbox().center()) );
					fs.draw();
				}

				if( n->property.contains("frame2") )
				{
					RMF::Frame f = n->property["frame2"].value<RMF::Frame>();
                    starlab::FrameSoup fs(0.01f);
					fs.addFrame(f.r, f.s, f.t, Vector3(n->bbox().center()) );
					fs.draw();
				}


				// For edges
				foreach(Structure::Link * l, edges)
				{
					if( l->property.contains("path") )
					{
						Array1D_Vector3 path;
						if (!l->property.contains("CachedPath"))
						{
							QVector< GraphDistance::PathPointPair > pathRelative = l->property["path"].value< QVector< GraphDistance::PathPointPair > >();
							if(!pathRelative.size()) continue;

							path = GraphDistance::positionalPath(this, pathRelative);
							l->property["CachedPath"].setValue(path);
						}
						else
							path = l->property["CachedPath"].value<Array1D_Vector3>();

                        starlab::PointSoup ps(10);
                        starlab::LineSegments ls;
						Vector3d lastP = path.front();
						foreach(Vector3 p, path)
						{
							ps.addPoint(p, Qt::green);
							ls.addLine(lastP, p, QColor(255,255,255,100));
							lastP = p;
						}
						ps.draw();
						ls.draw();
					}
				}

				glEnable(GL_LIGHTING);

			}
		}

		// Default node draw
		if(property["showNodes"].toBool())
		{
			glDisable(GL_CULL_FACE);
			n->draw( property["showCtrlPts"].toBool() );
		}

		if(n->type() == SHEET)
		{
			//Sheet * s = (Sheet*) n;

			//glBegin(GL_TRIANGLES);
			//foreach(std::vector<Vector3> tri, s->surface.generateSurfaceTris(0.4))
			//{
			//	glColor3d(1,0,0); glVector3(tri[0]);
			//	glColor3d(0,1,0); glVector3(tri[1]);
			//	glColor3d(0,0,1); glVector3(tri[2]);
			//}
			//glEnd();
		}
    }

	// Draw node names
	if( property["showNames"].toBool() && drawArea )
	{
		glColor3d(1,1,1);

		foreach(Node * n, nodes)
		{
			if (n->property.contains("isReady") && !n->property["isReady"].toBool()) continue;
			Vector3d position = n->bbox().center();

			if(property.contains("posX")) position[0] += property["posX"].toDouble();

			Vec proj = drawArea->camera()->projectedCoordinatesOf(Vec(position.x(), position.y(), position.z()));
			drawArea->renderText(proj.x,proj.y,n->id);
		}
	}

	if(property["showEdges"].toBool())
	{
		foreach(Link * e, edges)
		{
			if (e->n1->property.contains("isReady") && !e->n1->property["isReady"].toBool()) continue;
			if (e->n2->property.contains("isReady") && !e->n2->property["isReady"].toBool()) continue;

			e->draw();
		}
	}

	if(property.contains("selectedSample"))
	{
		glDisable(GL_LIGHTING);
		glColor3d(1,1,0); glPointSize(20);
		glBegin(GL_POINTS); glVector3(property["selectedSample"].value<Vector3d>()); glEnd();
	}

	glDisable(GL_LIGHTING);
	glColor3d(1,1,0); glPointSize(20);
	glBegin(GL_POINTS); foreach(Vector3 p, debugPoints) glVector3(p); glEnd();

	glColor3d(0,1,1); glPointSize(30);
	glBegin(GL_POINTS); foreach(Vector3 p, debugPoints2) glVector3(p); glEnd();

	glColor3d(1,0,1); glPointSize(40);
	glBegin(GL_POINTS); 
	int idx = 0;
	foreach(Vector3 p, debugPoints3){ 
		double t = double(idx) / debugPoints3.size(); 
		glColor3d(1,1,t*1);
		glVector3(p);
		idx++;
	} 
	glEnd();
	glEnable(GL_LIGHTING);

	// AABB
	if (property["showAABB"].toBool()) drawAABB();

    // LEGACY:
    /*
    // Materialized
    glBegin(GL_QUADS);
    foreach(QuadFace f, cached_mesh.faces)
    {
        Vector3 aa = cached_mesh.points[f[1]]-cached_mesh.points[f[0]];
        Vector3 bb = cached_mesh.points[f[2]]-cached_mesh.points[f[0]];
        Vector3 nn = cross(aa,bb).normalized();

        glNormal3(nn);
        glVector3(cached_mesh.points[f[0]]);
        glVector3(cached_mesh.points[f[1]]);
        glVector3(cached_mesh.points[f[2]]);
        glVector3(cached_mesh.points[f[3]]);
    }
    glEnd();

    glDisable(GL_LIGHTING);
    glPointSize(15);
    glColor3d(1,1,0);
    glBegin(GL_POINTS);
    foreach(Vector3 p, cached_mesh.debug) glVector3(p);
    glEnd();
    glEnable(GL_LIGHTING);
    */
}

void Graph::drawNodeMesh(QString nid, QColor meshColor)
{
    Node * n = getNode(nid);
    if(n->property.contains("mesh"))
    {
        QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
        QuickMeshDraw::drawMeshSolid( nodeMesh.data(), meshColor );
    }
}

void Graph::drawAABB()
{
	if (!property.contains("AABB"))
		property["AABB"].setValue( bbox() );

	Eigen::AlignedBox3d aabb = property["AABB"].value<Eigen::AlignedBox3d>();
	Eigen::Vector3d qc (aabb.center());
	Eigen::Vector3d diagonal (aabb.max() - aabb.min());

	double width = diagonal.x()/2;
	double length = diagonal.y()/2;
	double height = diagonal.z()/2;

	Vector3d center(qc.x(), qc.y(), qc.z());
	Vector3d  c1, c2, c3, c4;
	Vector3d  bc1, bc2, bc3, bc4;

	c1 = Vector3d (width, length, height) + center;
	c2 = Vector3d (-width, length, height) + center;
	c3 = Vector3d (-width, -length, height) + center;
	c4 = Vector3d (width, -length, height) + center;

	bc1 = Vector3d (width, length, -height) + center;
	bc2 = Vector3d (-width, length, -height) + center;
	bc3 = Vector3d (-width, -length, -height) + center;
	bc4 = Vector3d (width, -length, -height) + center;

	glDisable(GL_LIGHTING);
	glLineWidth(1.0f);

	glBegin(GL_LINES);
	glVertex3dv(c1.data());glVertex3dv(bc1.data());
	glVertex3dv(c2.data());glVertex3dv(bc2.data());
	glVertex3dv(c3.data());glVertex3dv(bc3.data());
	glVertex3dv(c4.data());glVertex3dv(bc4.data());
	glVertex3dv(c1.data());glVertex3dv(c2.data());
	glVertex3dv(c3.data());glVertex3dv(c4.data());
	glVertex3dv(c1.data());glVertex3dv(c4.data());
	glVertex3dv(c2.data());glVertex3dv(c3.data());
	glVertex3dv(bc1.data());glVertex3dv(bc2.data());
	glVertex3dv(bc3.data());glVertex3dv(bc4.data());
	glVertex3dv(bc1.data());glVertex3dv(bc4.data());
	glVertex3dv(bc2.data());glVertex3dv(bc3.data());
	glEnd();

	glEnable(GL_LIGHTING);

}

void Graph::draw2D(int width, int height)
{
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

    // Embed structure graph to 2D screen plane
    if(!property["embeded2D"].toBool())
	{
        //GraphEmbed::embed( this );
		GraphEmbed::circleEmbed( this );

		fontImage = QImage(":/images/font.png");
	}

	// Find 2D bounds
	Eigen::AlignedBox3d bound2D;
	foreach(Node * n, nodes){
		QVector3D center = n->vis_property[NODE_CENTER].value<QVector3D>();
		Eigen::Vector3d c(center.x(),center.y(),center.z());
		bound2D = bound2D.merged( Eigen::AlignedBox3d(c) );
	}
	Scalar boundRadius = 0.5 * bound2D.diagonal().norm();

	Eigen::Vector3d b2dmin = bound2D.min();
	QVector3D minBound(b2dmin.x(),b2dmin.y(),b2dmin.z());
    QVector3D corner = QVector3D(qMin(0.0f, minBound.x()), qMin(0.0f, minBound.y()), qMin(0.0f, minBound.z()));
	QVector3D graph_center = 0.5 * QVector3D(width, height,0);
	QVector3D scaling = QVector3D(width, height, 0) / boundRadius;

	// Compute 2D node centers
	QMap<int, Node*> nmap;
	std::vector<Vector3> node_centers;

	foreach(Node * n, nodes){
		QVector3D center = graph_center + (scaling * (n->vis_property[NODE_CENTER].value<QVector3D>() - corner));
		node_centers.push_back(Vector3(center.x(),center.y(),center.z()));
		nmap[nmap.size()] = n;
	}

	int N = nodes.size();

	// Draw edges:
	glColorQt(QColor(0,0,0));
	foreach(Link * e, edges)
	{
		Vector3 c1 = node_centers[nmap.key(e->n1)];
		Vector3 c2 = node_centers[nmap.key(e->n2)];
		drawLine(Vector2i(c1.x(), c1.y()), Vector2i(c2.x(), c2.y()), 1);
	}

	// Draw nodes:
	for(int i = 0; i < N; i++)
	{
		glColorQt( (nmap[i]->type() == SHEET ? QColor(50,0,0) : QColor(0,50,0)) );
		Vector3 c = node_centers[i];

		int w = stringWidthGL(qPrintable(nmap[i]->id)) * 1.1;
		int h = stringHeightGL();

		drawRect(Vector2i(c.x(), c.y()), w, h);
	}

	// Draw titles:
	glColorQt(QColor(255,255,255));
	beginTextDraw(fontImage);
	for(int i = 0; i < N; i++){
		Vector3 c = node_centers[i];

		int w = stringWidthGL(qPrintable(nmap[i]->id)) * 1.1;
		int h = stringHeightGL();

		drawStringGL( c.x() - (0.5 * w), c.y() - (0.5 * h), qPrintable(nmap[i]->id) );
	}
	endTextDraw();

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
}

void Graph::drawNodeMeshNames( int & offSet )
{
	foreach(Node * n, nodes){
		if(!n->property.contains("mesh")) continue;

		QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();

		n->property["meshSelectID"] = offSet;

		QuickMeshDraw::drawMeshName( nodeMesh.data(), offSet );

		offSet++;
	}
}

void Graph::saveToFile( QString fileName, bool isOutParts /*= true*/ ) const
{
	if(nodes.size() < 1) return;

	QFile file(fileName);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QFileInfo fileInfo(file.fileName());
	QDir graphDir( fileInfo.absolutePath() );

	QTextStream out(&file);
	out << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n<document>\n";

	// Save nodes
	foreach(Node * n, nodes)
	{
		out << "<node>\n";

		// Type and ID
		out << QString("\t<id>%1</id>\n").arg(n->id);
		out << QString("\t<type>%1</type>\n").arg(n->type());

		if(n->property.contains("mesh_filename"))
		{
			QString meshesFolder = "meshes";

			graphDir.mkdir( meshesFolder );

			QString mesh_filename = n->property["mesh_filename"].toString();
			QFileInfo meshFileInfo( mesh_filename );

			QString relativeFileName = meshesFolder + "/" + meshFileInfo.baseName() + ".obj";

			// Save mesh from memory
			QSharedPointer<SurfaceMeshModel> mesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
			if(mesh_filename.size() && mesh)
			{
				if( isOutParts )
					saveOBJ( mesh.data(), graphDir.path() + "/" + relativeFileName );
			}

			out << QString("\t<mesh>%1</mesh>\n\n").arg( relativeFileName ); // relative
		}

		// Control count
		out << "\t<controls>\n";
		std::vector<int> controlCount = n->controlCount();
		foreach(int c, controlCount)
			out << QString("\t\t<c>%1</c>\n").arg(c);
		out << "\t</controls>\n\n";

		// Control Points
		foreach(Vector3d p, n->controlPoints())
			out << QString("\t<point>%1 %2 %3</point>\n").arg(p.x()).arg(p.y()).arg(p.z());

		// Weights
		out << "\n\t<weights>";
		foreach(Scalar w, n->controlWeights())
			out << QString::number(w) + " ";
		out << "</weights>\n";

		// Meta data
        foreach(QString key, n->meta.keys()){
            QString value = n->meta[key].toString();
            out << "\n\t<meta><key>" << key << "</key><value>" << value << "</value></meta>\n";
		}

		out << "</node>\n\n";
	}

	// Save edges
	foreach(Link * e, edges)
	{
		out << "<edge>\n";
		out << QString("\t<id>%1</id>\n").arg(e->id);
		out << QString("\t<type>%1</type>\n\n").arg(e->type);
		out << QString("\t<n>%1</n>\n").arg(e->n1->id);
		out << QString("\t<n>%1</n>\n").arg(e->n2->id);

		for(int k = 0; k < 2; k++)
		{
			out << "\t<coord>\n";
			foreach(Vector4d c, e->coord[k])
				out << "\t\t<uv>" << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << "</uv>\n";
			out << "\t</coord>\n";
		}

		out << "</edge>\n\n";
	}

	// Save groups
	foreach(QVector<QString> group, groups){
		out << "<group>";
		foreach(QString nid, group){
			out << QString("\n\t<n>%1</n>").arg(nid);
		}
		out << "\n</group>\n\n";
	}

	out << "\n</document>\n";
	file.close();
}

void Graph::loadFromFile( QString fileName )
{
	// Clear data
	nodes.clear();
	edges.clear();

	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;
	QFileInfo fileInfo(file.fileName());

	QDomDocument mDocument;
	mDocument.setContent(&file, false);    
	
	int degree = 3;

	// For each node
	QDomNodeList node_list = mDocument.firstChildElement("document").elementsByTagName("node");
	int num_nodes = node_list.count();

	bool hasMeshes = false;

	for(int i = 0; i < num_nodes; i++)
	{
		QDomNode node = node_list.at(i);

		QString id = node.firstChildElement("id").text();
		QString node_type = node.firstChildElement("type").text();
		QString mesh_filename = node.firstChildElement("mesh").text();
		
		// Find number of control elements
		QVector<int> control_count;
		QDomNodeList controls_list = node.firstChildElement("controls").toElement().elementsByTagName("c");
		for(int c = 0; c < (int)controls_list.size(); c++)
			control_count.push_back( controls_list.at(c).toElement().text().toInt() );

		// Load all control points
		std::vector<Vector3d> ctrlPoints;
		QDomNode n = node.firstChildElement("point");
		while (!n.isNull()) {
			if (n.isElement()) {
				QStringList point = n.toElement().text().split(" ");
				ctrlPoints.push_back( Vector3d(point[0].toDouble(),point[1].toDouble(),point[2].toDouble()) );
			}
			n = n.nextSiblingElement("point");
		}

		// Load all control weights
		std::vector<Scalar> ctrlWeights;
		QStringList weights = node.firstChildElement("weights").toElement().text().trimmed().split(" ");
		foreach(QString w, weights) ctrlWeights.push_back(w.toDouble());

		// Add node
		Node * new_node = NULL;

		if(node_type == CURVE)
		{
            new_node = addNode( new Curve( NURBS::NURBSCurved(ctrlPoints, ctrlWeights, degree, false, true), id) );
		}
		else if(node_type == SHEET)
		{
			if(control_count.size() < 2) continue;

			std::vector< std::vector<Vector3d> > cp = std::vector< std::vector<Vector3d> > (control_count.first(), std::vector<Vector3d>(control_count.last(), Vector3(0,0,0)));
			std::vector< std::vector<Scalar> > cw = std::vector< std::vector<Scalar> > (control_count.first(), std::vector<Scalar>(control_count.last(), 1.0));

			for(int u = 0; u < control_count.first(); u++)
			{
				for(int v = 0; v < control_count.last(); v++)
				{
					int idx = (u * control_count.last()) + v;
					cp[u][v] = ctrlPoints[idx];
					cw[u][v] = ctrlWeights[idx];
				}
			}

            new_node = addNode( new Sheet( NURBS::NURBSRectangled(cp, cw, degree, degree, false, false, true, true), id ) );
		}

		// Load meta data
        PropertyMap meta_values;
		QDomNodeList meta_list = node.toElement().elementsByTagName("meta");
		int meta_nodes = meta_list.count();
		for(int m = 0; m < meta_nodes; m++)
		{
			QDomNode meta = meta_list.at(m);
			QDomNode key = meta.firstChildElement("key");
            QDomNode v = meta.firstChildElement("value");

            meta_values[ key.toElement().text() ] = v.toElement().text();
		}
        new_node->meta = meta_values;

		// Mesh file path
		new_node->property["mesh_filename"].setValue( mesh_filename );
		QString fullMeshPath = fileInfo.dir().path() + "/" + mesh_filename;
		QFile mfile( fullMeshPath );

		// Load node's mesh and make it ready
		if (mfile.exists())
		{
			QSharedPointer<SurfaceMeshModel> nodeMesh(new SurfaceMeshModel(mesh_filename, id));

			nodeMesh->read( qPrintable(fullMeshPath) );
			nodeMesh->update_face_normals();
			nodeMesh->update_vertex_normals();
			nodeMesh->updateBoundingBox();
			new_node->property["mesh"].setValue( nodeMesh );

			hasMeshes = true;
		}
	}

	// Original shape bounding box
	if( hasMeshes ){
		Eigen::AlignedBox3d shapeBox;
		for(int i = 0; i < nodes.size(); i++){
			QSharedPointer<SurfaceMeshModel> nodeMesh = nodes[i]->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
			shapeBox = shapeBox.merged( nodeMesh->bbox() );
		}
		property["shapeBox"].setValue( shapeBox );
		property["hasMeshes"].setValue( hasMeshes );
	}

	// For each edge
	QDomNodeList edge_list = mDocument.firstChildElement("document").elementsByTagName("edge");
	int num_edges = edge_list.count();

	for(int i = 0; i < num_edges; i++)
	{
		QDomNode edge = edge_list.at(i);

		QString id = edge.firstChildElement("id").text();
		QString edge_type = edge.firstChildElement("type").text();

		QDomNodeList n = edge.toElement().elementsByTagName("n");
		
		QString n1_id = n.at(0).toElement().text();
		QString n2_id = n.at(1).toElement().text();

		QDomNodeList coordList = edge.toElement().elementsByTagName("coord");
		Array2D_Vector4d coords((int)coordList.size());
		for(int j = 0; j < (int) coords.size(); j++)
		{
			QDomNodeList uv = coordList.at(j).toElement().elementsByTagName("uv");
			int uv_count = uv.count();

			for(int k = 0; k < uv_count; k++)
			{
				QStringList c = uv.at(k).toElement().text().split(" ");
				coords[j].push_back(Vector4d(c[0].toDouble(), c[1].toDouble(), c[2].toDouble(), c[3].toDouble()));
			}
		}

		addEdge(getNode(n1_id), getNode(n2_id), coords.front(), coords.back(), id);
	}

    // For each group
    QDomNodeList group_list = mDocument.firstChildElement("document").elementsByTagName("group");
    int num_groups = group_list.count();
    for(int i = 0; i < num_groups; i++)
    {
        QDomNode group = group_list.at(i);

        // Find elements of group
        QVector<QString> elements_ids;
        QDomNodeList elements_list = group.toElement().elementsByTagName("n");
        for(int c = 0; c < (int)elements_list.size(); c++)
            elements_ids.push_back( elements_list.at(c).toElement().text() );

        QColor groupColor = starlab::qRandomColor2();

        QVector<QString> element_nodes;
        foreach(QString nid, elements_ids){
			Node * n = getNode(nid);
			if(n){
				element_nodes.push_back( nid );
				n->vis_property["color"].setValue( groupColor );
			}
        }
        addGroup(element_nodes);
    }

	file.close();
}

void Graph::exportAsOBJ( QString filename )
{
	QFile file(filename);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) return;
	QFileInfo fileInfo(file.fileName());
	QTextStream out(&file);
	out << "# Exported graph model [" << fileInfo.baseName() << "] by TopoBlender\n\n";

	int v_offset = 0;

	foreach(Structure::Node * n, nodes)
	{
		if(!n->property.contains("mesh")) continue;

		out << "# Starting mesh " << n->id << "\n";

		SurfaceMesh::Model* m = getMesh(n->id);

		// Write out vertices
		Vector3VertexProperty points = m->vertex_property<Vector3>(VPOINT);
		foreach( Vertex v, m->vertices() )
			out << "v " << points[v][0] << " " << points[v][1] << " " << points[v][2] << "\n";

		// Write out triangles
		out << "g " << n->id << "\n";
		foreach( Face f, m->faces() ){
			out << "f ";
			Surface_mesh::Vertex_around_face_circulator fvit = m->vertices(f), fvend = fvit;
			do{	out << (((Surface_mesh::Vertex)fvit).idx() + 1 + v_offset) << " ";} while (++fvit != fvend);
			out << "\n";
		}

		v_offset += m->n_vertices();

		out << "# End of mesh " << n->id << "\n\n";
	}

	file.close();
}

Node *Graph::rootBySize()
{
	int idx = 0;
	Scalar maxArea = 0;

	for(int i = 0; i < nodes.size(); i++){
		Vector3 diagonal = nodes[i]->bbox().diagonal();

		double area = (diagonal[0] == 0 ? 1 : diagonal[0]) * 
					  (diagonal[1] == 0 ? 1 : diagonal[1]) * 
					  (diagonal[2] == 0 ? 1 : diagonal[2]);

		if( area > maxArea){
			idx = i;
			maxArea = area;
		}
	}

    return nodes[idx];
}

Node *Graph::rootByValence()
{
	int maxIdx = 0;
	int maxValence = valence(nodes[maxIdx]);

	for(int i = 0; i < (int)nodes.size(); i++)
	{
		int curVal = valence(nodes[i]);
		if(curVal > maxValence)
		{
			maxValence = curVal;
			maxIdx = i;
		}
	}

    return nodes[maxIdx];
}

SurfaceMesh::Vector3 Graph::nodeIntersection( Node * n1, Node * n2 )
{
	double s1 = n1->bbox().diagonal().norm();
	double s2 = n2->bbox().diagonal().norm();

	Scalar r = 0.04 * qMin(s1, s2);

	//if(n1->type() == SHEET && n2->type() == SHEET)
	//	r *= 10;

	Array2D_Vector3 parts1 = n1->discretized(r * (n1->type() == CURVE ? 1 : 1));
	Array2D_Vector3 parts2 = n2->discretized(r * (n1->type() == CURVE ? 1 : 1));

	// Fall back
	if( !parts1.size() )
	{
		if(n1->type() == SHEET && n2->type() == SHEET)
		{
			Sheet * s1 = (Sheet *)n1;
			Sheet * s2 = (Sheet *)n2;

			parts1 = s1->surface.mCtrlPoint;
			parts2 = s2->surface.mCtrlPoint;
		}
	}

	Scalar minDist = DBL_MAX;
	int minI = 0, minJ = 0;
	int minIm = 0, minJn = 0;

	for(int i = 0; i < (int)parts1.size(); i++)
	{
		for(int j = 0; j < (int)parts2.size(); j++)
		{
			Scalar local_min_dis = DBL_MAX;
			int min_m = 0;
			int min_n = 0;

			for (int m = 0; m < (int)parts1[i].size(); m++)
			{
				for (int n = 0;  n < (int)parts2[j].size(); n++ )
				{
					double dis = (parts1[i][m] - parts2[j][n]).norm();
					if (dis < local_min_dis)
					{
						local_min_dis = dis;
						min_m = m;
						min_n = n;
					}
				}
			}

			if (local_min_dis < minDist)
			{
				minDist = local_min_dis;
				minI = i;
				minJ = j;

				minIm = min_m;
				minJn = min_n;
			}

		}
	}

	Vector3 p1 = parts1[minI][minIm];
	Vector3 p2 = parts2[minJ][minJn];
	return (p1 + p2) / 2;

	//// Compare parts bounding boxes
	//for(int i = 0; i < (int)parts1.size(); i++)
	//{
	//	Vector3 mean1(0);
	//	Scalar r1 = 0;

	//	foreach(Vector3 p, parts1[i])	mean1 += p;
	//	mean1 /= parts1[i].size();
	//	foreach(Vector3 p, parts1[i])	r1 = qMax((p - mean1).norm(), r1);

	//	for(int j = 0; j < (int)parts2.size(); j++)
	//	{
	//		Vector3 mean2(0);
	//		Scalar r2 = 0;

	//		foreach(Vector3 p, parts2[j])	mean2 += p;
	//		mean2 /= parts2[j].size();
	//		foreach(Vector3 p, parts2[j])	r2 = qMax((p - mean2).norm(), r2);

	//		Vector3 diff = mean1 - mean2;
	//		Scalar dist = diff.norm();

	//		Scalar sphereDist = dist - r1 - r2;

	//		if(sphereDist <= 0)
	//		{
	//			std::vector<Vector3> p1 = parts1[ i ];
	//			std::vector<Vector3> p2 = parts2[ j ];
	//			int g1 = p1.size(), g2 = p2.size();
	//			Vector3 c1,c2;
	//			Scalar s,t;

	//			if(g1 == 2 && g2 == 2) ClosestPointSegments		(p1[0],p1[1],  p2[0],p2[1],       s, t, c1, c2);	// curve -- curve
	//			if(g1 == 2 && g2 == 3) ClosestSegmentTriangle	(p1[0],p1[1],  p2[0],p2[1],p2[2],       c1, c2);	// curve -- sheet
	//			if(g1 == 3 && g2 == 2) ClosestSegmentTriangle	(p2[0],p2[1],  p1[0],p1[1],p1[2],       c1, c2);	// sheet -- curve
	//			if(g1 == 3 && g2 == 3) TriTriIntersect			(p1[0],p1[1],p1[2],  p2[0],p2[1],p2[2], c1, c2);	// sheet -- sheet
	//			
	//			Scalar dist = (c1 - c2).norm();

	//			if(dist < minDist)
	//			{
	//				minI = i;
	//				minJ = j;
	//				minDist = dist;
	//			}

	//			//foreach(Vector3 w, parts1[i]) debugPoints2.push_back(w);
	//			//foreach(Vector3 w, parts2[j]) debugPoints2.push_back(w);
	//		}
	//	}
	//}

	//// Find the exact intersection point
	//std::vector<Vector3> p1 = parts1[ minI ];
	//std::vector<Vector3> p2 = parts2[ minJ ];
	//int g1 = p1.size(), g2 = p2.size();

	//Vector3 c1,c2;
	//Scalar s,t;

	//if(g1 == 2 && g2 == 2) ClosestPointSegments		(p1[0],p1[1],  p2[0],p2[1],       s, t, c1, c2);	// curve -- curve
	//if(g1 == 2 && g2 == 3) ClosestSegmentTriangle	(p1[0],p1[1],  p2[0],p2[1],p2[2],       c1, c2);	// curve -- sheet
	//if(g1 == 3 && g2 == 2) ClosestSegmentTriangle	(p2[0],p2[1],  p1[0],p1[1],p1[2],       c1, c2);	// sheet -- curve
	//if(g1 == 3 && g2 == 3) TriTriIntersect			(p1[0],p1[1],p1[2],  p2[0],p2[1],p2[2], c1, c2);	// sheet -- sheet

	//debugPoints.push_back(c1);
	//debugPoints.push_back(c2);
	//foreach(Vector3 w, p1) debugPoints2.push_back(w);
	//foreach(Vector3 w, p2) debugPoints3.push_back(w);

	//return (c1 + c2) / 2.0;
}

int Graph::valence( Node * n )
{
	// Collect set of neighbours
	QSet<QString> neighbours;
	foreach(Link * l, this->getEdges(n->id))
		neighbours.insert( l->otherNode( n->id )->id );
	
	// remove any self edges
	neighbours.remove(n->id);

	return neighbours.size();
}

Curve* Graph::getCurve( Link * l )
{
	Node *n1 = l->n1, *n2 = l->n2;
	return (Curve *) ((n1->type() == CURVE) ? n1: n2);
}

QVector<Link*> Graph::getEdges( QString nodeID )
{
	QVector<Link*> mylinks;

	for(int i = 0; i < edges.size(); i++)
	{
		Link * l = edges[i];
		if( l->hasNode(nodeID) ) mylinks.push_back(l);
	}

	return mylinks;
}

QVector<Node*> Graph::adjNodes( Node * node )
{
	QVector<Node*> adj;
	foreach(Link* edge, getEdges(node->id))
		adj.push_back( edge->otherNode(node->id) );
	return adj;
}

QMap<Link*, Array1D_Vector4d > Graph::linksCoords( QString nodeID )
{
    QMap< Link*, Array1D_Vector4d > coords;

	for(int i = 0; i < edges.size(); i++)
	{
		Link * l = edges[i];
		if(!l->hasNode(nodeID)) continue;

		coords[l] = l->getCoord(nodeID);
	}

	return coords;
}

QVector<Link*> Graph::nodeEdges( QString nodeID )
{
	QVector<Link*> nodeLinks;
	
	foreach(Link * e, edges) if(e->hasNode(nodeID)) 
		nodeLinks.push_back(e);

	return nodeLinks;
}

void Graph::removeNode( QString nodeID )
{
	Node * n = getNode(nodeID);

	foreach(Link * e, getEdges(nodeID)){
		int edge_idx = edges.indexOf(e);

		delete edges[edge_idx];
		edges[edge_idx] = NULL;
		edges.remove(edge_idx);
	}

	int node_idx = nodes.indexOf(n);

	if( node_idx < 0) return;

	delete nodes[node_idx];
	nodes[node_idx] = NULL;
	nodes.remove(node_idx);
}

void Graph::removeIsolatedNodes()
{
	foreach(Node * n, nodes){
		if( getEdges(n->id).size() == 0 )
			removeNode(n->id);
	}
}

QList<Link*> Graph::furthermostEdges( QString nodeID )
{
	QMap<double, Link*> sortedLinks;
	foreach(Link * e, nodeEdges(nodeID))
	{
		foreach(Vector4d c, e->getCoord(nodeID))
			sortedLinks[c.norm()] = e;
	}

	return sortedLinks.values();
}

QSet<QString> Graph::nodesCanVisit( Node * node )
{
	QMap<QString, bool> visitedNodes;
	QStack<QString> nodesToVisit;

	nodesToVisit.push( node->id );

	while(! nodesToVisit.isEmpty() )
	{
		QString id = nodesToVisit.pop();
		visitedNodes[id] = true;

		foreach(Link * l, getEdges(id))
		{
			QString adjID = l->otherNode(id)->id;

			if( !visitedNodes.contains(adjID) && !nodesToVisit.contains(adjID))
				nodesToVisit.push(adjID);
		}
	}

	return visitedNodes.keys().toSet();
}

int Graph::numCanVisit( Node * node )
{
	return nodesCanVisit( node ).size();
}

QVector< QSet<QString> > Graph::connectedComponents()
{
	QVector< QSet<QString> > connected;

	foreach(Node * n, nodes)
	{
		QSet<QString> curSet = nodesCanVisit( n );

		if(!connected.contains(curSet))
			connected.push_back( curSet );
	}

	return connected;
}

bool Graph::isConnected()
{
	return numCanVisit( nodes.front() ) == nodes.size();
}

bool Graph::isCutNode( QString nodeID )
{
	QSet<Link*> isDeleted;

	QMap<QString, bool> visitedNodes;
	QStack<QString> nodesToVisit;

	foreach(Link * l, getEdges(nodeID))
	{
		if(nodesToVisit.size() == 0) 
			nodesToVisit.push( l->otherNode(nodeID)->id );
		isDeleted.insert(l);
	}

	// Visit all nodes without going through a deleted edge
	while( ! nodesToVisit.isEmpty() )
	{
		QString id = nodesToVisit.pop();
		visitedNodes[id] = true;

		foreach(Link * l, getEdges(id))
		{
			if(isDeleted.contains(l)) continue;

			QString adjID = l->otherNode(id)->id;

			if( !visitedNodes.contains(adjID) && !nodesToVisit.contains(adjID))
				nodesToVisit.push(adjID);
		}
	}

	// Skip flying nodes
	int countConnectNodes = 0;
	foreach(Node * n, nodes)
	{
		if(this->getEdges(n->id).size() > 0)
		{
			countConnectNodes++;
		}
	}

	int countVisitNodes = visitedNodes.size();

	if( countVisitNodes != countConnectNodes - 1 )
		return true;
	else
		return false;
}

bool Graph::isInCutGroup( QString nodeID )
{
	Structure::Graph copy(*this);

	foreach(QVector<QString> nids, groupsOf(nodeID)){
		foreach(QString nid, nids) 
			if(copy.getNode(nid))
				copy.removeNode(nid);
	}

	return !copy.isConnected();
}

bool Graph::isBridgeEdge( Link * link )
{
	Graph g = *this;
	
	int beforeCount = g.numCanVisit( link->n1 );

	g.removeEdge(link->n1, link->n2);

	int afterCount = g.numCanVisit( link->n1 );

	return beforeCount != afterCount;
}

QVector< QVector<Node*> > Graph::split( QString nodeID )
{
	QVector< QVector<Node*> > parts;

	Node * n = getNode(nodeID);
	if(!n) return parts;

	QMap<QString, bool> visitedNodes;
	visitedNodes[nodeID] = true;

	foreach( Link * link, this->getEdges(nodeID) )
	{
		Node * other = link->otherNode(nodeID);

		QStack<QString> nodesToVisit;
		nodesToVisit.push(other->id);

		QSet<Node*> curSet;

		while( ! nodesToVisit.isEmpty() )
		{
			QString id = nodesToVisit.pop();
			if(visitedNodes[id]) continue;

			visitedNodes[id] = true;

			curSet.insert(getNode(id));

			foreach( Link * l, getEdges(id) )
			{
				if(l == link) continue;

				QString adjID = l->otherNode(id)->id;

				if( !visitedNodes.contains(adjID) && !nodesToVisit.contains(adjID))
					nodesToVisit.push(adjID);
			}
		}

		QVector<Node*> part;
		foreach(Node* n, curSet) part.push_back(n);
		if(part.size())	parts.push_back( part );
	}

	return parts;
}

void Graph::replaceCoords( QString nodeA, QString nodeB, Array1D_Vector4d coordA, Array1D_Vector4d coordB )
{
	// Get shared link
	Link * sharedLink = getEdge(nodeA, nodeB);

	if(!sharedLink) return;

	sharedLink->setCoord(nodeA, coordA);
	sharedLink->setCoord(nodeB, coordB);
}

Vector3 Graph::position( QString nodeID, Vector4d& coord )
{
	Node * node = getNode(nodeID);
	if(!node) return bbox().center(); // Something went wrong..
	return node->position(coord);
}

SurfaceMesh::Model* Graph::getMesh( QString nodeID )
{
	Node * node = getNode(nodeID);

    if(node && node->property.contains("mesh")) return node->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >().data();

	return NULL;
}

void Graph::translate( Vector3 delta, bool isKeepMeshes )
{
	foreach (Node * node, nodes)
	{
		node->moveBy( delta );

		// Update needed for Sheets
		if(node->type() == SHEET)
			((Sheet*)node)->surface.quads.clear();

		// Apply to actual geometry
		if( !isKeepMeshes )
		{
			if(!node->property.contains("mesh")) continue;
			SurfaceMesh::Model* model = getMesh(node->id);

            if(model)
            {
                Vector3VertexProperty points = model->vertex_property<Vector3d>("v:point");
                foreach(Vertex v, model->vertices()) points[v] += delta;
            }
		}
	}

	// Update the bounding box
	property["AABB"].setValue(bbox());
}

void Graph::rotate( double angle, Vector3 axis )
{
	foreach (Node * node, nodes)
	{
		node->rotate(angle, axis);

		if(!node->property.contains("mesh")) continue;

		// Update needed for Sheets
		if(node->type() == SHEET)
			((Sheet*)node)->surface.quads.clear();

		// Move actual geometry
        SurfaceMesh::Model* model = getMesh(node->id);
		Vector3VertexProperty points = model->vertex_property<Vector3d>("v:point");
		model->updateBoundingBox();

		foreach(Vertex v, model->vertices())
		{
			points[v] = rotatedVec(points[v],  (3.14159265358979 /180) * angle, axis);
		}
	}
}

void Graph::scale( double scaleFactor )
{
	double current_scale = 1.0;
	if (property.contains("scale"))
		current_scale = property["scale"].toDouble();

	double relative_scale = scaleFactor / current_scale;

	foreach (Node * node, nodes)
	{
		node->scale(relative_scale);

		// Update needed for Sheets
		if(node->type() == SHEET)
			((Sheet*)node)->surface.quads.clear();

		if(!node->property.contains("mesh")) continue;

		// Move actual geometry
        SurfaceMesh::Model* model = getMesh(node->id);
		Vector3VertexProperty points = model->vertex_property<Vector3d>("v:point");
		foreach(Vertex v, model->vertices())
			points[v] *= relative_scale;
	}

	property["scale"] = scaleFactor;

	property["AABB"].setValue(bbox());
}

void Graph::transform( QMatrix4x4 mat )
{
	Vector3 c = bbox().center();

	// Selective transformation
	bool isSelective = false;
	foreach (Node * node, nodes) if(node->selections.size()) isSelective = true;

	foreach (Node * node, nodes)
	{
		if(isSelective && !node->selections.size()) continue;
		
		Array1D_Vector3 controlPoints = node->controlPoints();
		for(int i = 0; i < (int)controlPoints.size(); i++)
            controlPoints[i] = starlab::QVector3(mat * starlab::QVector3(controlPoints[i]));
		node->setControlPoints(controlPoints);
		
		// Update needed for Sheets
		if(node->type() == SHEET)
			((Sheet*)node)->surface.quads.clear();
		
		// Transform actual geometry
		if(!node->property.contains("mesh")) continue;
		SurfaceMesh::Model* model = getMesh( node->id );
		Vector3VertexProperty points = model->vertex_property<Vector3d>("v:point");
        foreach(Vertex v, model->vertices()) points[v] = starlab::QVector3(mat * starlab::QVector3(points[v]));
	}
}

void Graph::moveBottomCenterToOrigin(bool isKeepMeshes)
{
	Eigen::AlignedBox3d aabb = bbox();
	double height = aabb.max().z() - aabb.min().z();
	Vector3d bottom_center(aabb.center().x(), aabb.center().y(), aabb.center().z() - height/2);

    translate( -bottom_center, isKeepMeshes );

	// Update the bounding box
	property["AABB"].setValue(bbox());
}

void Graph::normalize()
{
	// Normalize the height to be 1
	Eigen::AlignedBox3d aabb = bbox();
	double height = aabb.max().z() - aabb.min().z();

	if(height == 0.0) 
		height = aabb.max().x() - aabb.min().x();

	double scaleFactor = 1.0 / height;

	foreach (Node * node, nodes)
	{
		node->scale(scaleFactor);

		// Apply to actual geometry
		if(!node->property.contains("mesh")) continue;
		SurfaceMesh::Model* model = getMesh( node->id );
		Vector3VertexProperty points = model->vertex_property<Vector3d>("v:point");
		foreach(Vertex v, model->vertices())
			points[v] *= scaleFactor;
	}

	// Rebuild visualization geometry
	foreach (Node * node, nodes){
		if(node->type() == SHEET)
		{
			Sheet * sheet = (Sheet *)node;
			sheet->surface.quads.clear();
		}
	}

	// Update the bounding box
	property["AABB"].setValue(bbox());
}

void Graph::addGroup(QVector<QString> newGroup)
{
    if(!newGroup.size()) return;

    int gid = groups.size();
    groups.push_back( newGroup );

    // Tell nodes they belong to group
    foreach (QString nid, newGroup){
        getNode(nid)->property["groupID"] = gid;
    }
}

void Graph::removeGroup( QVector<QString> groupElements )
{
	int idx = -1;
	for(int i = 0; i < groups.size(); i++) if(groups[i] == groupElements) idx = i;

	removeGroup(idx);
}

void Graph::removeGroup(int groupIDX)
{
	if(groupIDX < 0) return;
    if(!groups.size()) return;

    // Remove
    groups.remove(groupIDX);
}

QVector< QVector<QString> > Graph::groupsOf( QString nodeID )
{
	QVector< QVector<QString> > result;

	foreach(QVector<QString> group, groups)
	{
		foreach(QString nid, group)
		{
			if(nid == nodeID){
				result.push_back(group);
				break;
			}
		}
	}

	if(result.isEmpty()) result.push_back( QVector<QString>() );

    return result;
}

QVector< QVector<QString> > Graph::nodesAsGroups()
{
    QVector< QVector<QString> > asgroups;

    QSet<QString> nodeSet;
    foreach(Node * n, nodes) nodeSet.insert(n->id);

    foreach(QVector<QString> group, groups){
        foreach(QString nid, group) nodeSet.remove(nid);
        asgroups.push_back(group);
    }

    foreach(QString nid, nodeSet)
        asgroups.push_back( QVector<QString>() << nid );

    return asgroups;
}

QVector<POINT_ID> Graph::selectedControlPointsByColor(QColor color)
{
	QVector<POINT_ID> result;
	for (int nID = 0; nID < (int)nodes.size(); nID++)
	{
		Node *node = nodes[nID];
		foreach (int pID, node->selections.keys())
		{
			if (node->selections[pID] == color)
				result.push_back(std::make_pair(nID, pID));
		}
	}

	return result;
}

void Graph::clearSelections()
{
	foreach(Node *node, nodes)
		node->selections.clear();
}

QString Graph::name()
{
    return property["name"].toString().section('\\', -1).section('/', -1).section('.', 0, 0);
}

void Graph::setColorAll( QColor newNodesColor )
{
	foreach(Node * n, nodes)
		n->vis_property["color"] = newNodesColor;
}

void Graph::setColorFor( QString nodeID, QColor newColor )
{
	Node * n = getNode(nodeID);
    if(!n) return;

	n->vis_property["color"] = newColor;

	QSharedPointer<SurfaceMeshModel> nodeMesh = n->property["mesh"].value< QSharedPointer<SurfaceMeshModel> >();
	if(!nodeMesh.isNull())	n->vis_property["meshColor"] = newColor;
}

void Graph::renameNode( QString oldNodeID, QString newNodeID )
{
	Structure::Node * n = getNode(oldNodeID);
	if(!n) return;

	n->id = n->id.replace(oldNodeID, newNodeID);

	foreach(Link * l, getEdges(oldNodeID))
	{
		if(l->id.contains(oldNodeID))
			l->id = l->id.replace(oldNodeID, newNodeID);
	}
}

void Graph::clearAll()
{
	foreach(Node * n, nodes){
		n->property.clear();
	}
}

QVector<Node*> Graph::articulationPoints()
{
	int N = nodes.size();
	setPropertyAll("articulation", false);

	QVector<int> low, pre;
	low.fill(-1,N);
	pre.fill(-1,N);

	int cnt = 0;

	for(int v = 0; v < nodes.size(); v++)
		if (pre[v] == -1)
			articulationDFS(cnt, v, v, low, pre);

	return nodesWithProperty("articulation", true);
}

void Graph::articulationDFS( int & cnt, int u, int v, QVector<int> & low, QVector<int> & pre )
{
	int children = 0;
	pre[v] = cnt++;
	low[v] = pre[v];
	
	foreach(Link * l, getEdges(nodes[v]->id))
	{
		int w = indexOfNode(l->otherNode(nodes[v]->id));

		if (pre[w] == -1) {
			children++;
			articulationDFS(cnt, v, w, low, pre);

			// update low number
			low[v] = qMin(low[v], low[w]);

			// non-root of DFS is an articulation point if low[w] >= pre[v]
			if (low[w] >= pre[v] && u != v) 
				nodes[v]->property["articulation"] = true;
		}

		// update low number - ignore reverse of edge leading to v
		else if (w != u)
			low[v] = qMin(low[v], pre[w]);
	}

	// root of DFS is an articulation point if it has more than 1 child
	if (u == v && children > 1)
		nodes[v]->property["articulation"] = true;
}

QVector<Node*> Graph::leaves()
{
	QVector<Node*> result;

	foreach(Node* n, nodes)
		if(valence(n) == 1)
			result.push_back(n);

	return result;
}

void Graph::moveCenterTo( Vector3 newCenter, bool isKeepMeshes )
{
	Eigen::AlignedBox3d curBox = bbox(true);
	Vector3 delta = newCenter - curBox.center();
	
	if(!isKeepMeshes) 
		translate( delta );
	else
	{
		foreach (Node * node, nodes)
		{
			node->moveBy( delta );
			if(node->type() == SHEET) ((Sheet*)node)->surface.quads.clear();
		}
	}
}

Structure::Graph * Graph::actualGraph(Structure::Graph * fromGraph)
{
	Graph * actual = new Graph(*fromGraph);

	bool stillCleaning = true;

	while( stillCleaning )
	{
		stillCleaning = false;

		foreach(Node* n, actual->nodes)
		{
			QVector<Link*> edges = actual->getEdges(n->id);

			// Skip disconnected
			if(edges.size() == 0) continue;

			// Get type of task on this node
			if(!n->property.contains("taskTypeReal")) continue;
			int taskType = n->property["taskTypeReal"].toInt();

			// Real morph tasks are not topology altering
			if(taskType == Task::MORPH) continue;

			// Cases:
			//	Nodes that finished shrinking or 
			//	have not yet grown or
			//	Going to split
			bool isShrinking = taskType == Task::SHRINK && n->property["taskIsDone"].toBool();
			bool isGrowing = taskType == Task::GROW && !n->property["taskIsReady"].toBool();
			bool isSpliting = taskType == Task::SPLIT && !n->property["taskIsReady"].toBool();

			if( isShrinking || isGrowing || isSpliting )
			{
				/// Transfer edges to someone else:
				if( !isSpliting )
				{
					// Look at neighbours valence
					QMap<int, Node*> valMap;
					foreach(Link* e, edges){
						Node * j = e->otherNode(n->id);
						valMap[actual->valence(j)] = j;
					}

					// Move my edges to my replacment
					Node* replacment = valMap[valMap.keys().back()];
					foreach(Link* e, edges)
					{
						QString otherID = e->otherNode(n->id)->id;
						if(otherID == replacment->id) continue;
						if(actual->shareEdge(otherID, replacment->id)) continue;

						Vec4d c(0,0,0,0);

						// More accurate but too slow..
						//Vec4d c = replacment->approxCoordinates( replacment->approxProjection(e->position(n->id)) );

						if(replacment->type() == Structure::CURVE)
						{
							double t = ((Structure::Curve*)replacment)->curve.fastTimeAt(e->position(n->id));
							c = Vec4d(t,0,0,0);
						}
						else
						{
							c = ((Structure::Sheet*)replacment)->surface.fastTimeAt(e->position(n->id));
						}

						Array1D_Vector4d coord(1, c);

						e->replace(n->id, replacment, coord);
					}
				}

				foreach(Link* e, actual->getEdges(n->id))
				{
					actual->removeEdge(e->n1, e->n2);
				}

				if( isSpliting )
				{
					QString siblingKeep;
					QVector<QString> siblings = actual->groupsOf(n->id).back();

					// Change one sibling in order to keep it
					foreach(QString nid, siblings)
					{
						if(nid == n->id) continue;
						
						// Override the chosen one
						actual->getNode(nid)->property["taskTypeReal"] = Task::MORPH;
						actual->removeGroup(siblings);
						break;
					}
				}

				stillCleaning = true;

				break;
			}
		}
	}

	// Remove edges kept after a merge
	QVector<Link*> toRemove;
	foreach(Link * link, actual->edges){
		if(link->property["mergedEdge"].toBool())
			toRemove.push_back(link);
	}
	foreach(Link * link, toRemove){
		actual->removeEdge( link->property["uid"].toInt() );
	}

	// Remove any disconnected nodes
	QVector<QString> nodesToRemove;
	foreach(Node * n, actual->nodes) if(actual->getEdges(n->id).size() == 0) nodesToRemove.push_back(n->id);
	foreach(QString nid, nodesToRemove) actual->removeNode(nid);

	return actual;
}
