#include "TaskSheet.h"
#include "AbsoluteOrientation.h"

#include "StructureCurve.h"

using namespace Structure;

void TaskSheet::prepareSheet()
{
    switch(type)
    {
    case GROW:
    case SHRINK:
        prepareGrowShrinkSheet();
        break;
    case SPLIT:
    case MERGE:
    case MORPH:
        prepareMorphSheet();
        break;
    }
}

void TaskSheet::executeSheet(double t)
{
    switch(type)
    {
    case GROW:
    case SHRINK:
        executeGrowShrinkSheet(t);
        break;
    case MORPH:
		{
			if( property["isCrossing"].toBool() )
				executeCrossingSheet(t);
			else
				executeMorphSheet(t);
		}
        break;
    }
}

void TaskSheet::prepareSheetOneEdge( Structure::Link * l )
{
    Structure::Node * n = node(), *tn = targetNode();
    Structure::Node * base = l->otherNode(n->id);
    Structure::Sheet* structure_sheet = ((Structure::Sheet*)n);
	Vector4d coord = l->getCoord(n->id).front();
	Vector4d otherCoord = l->getCoordOther(n->id).front();

    // Placement:
    structure_sheet->moveBy( base->position(otherCoord) - n->position(coord) );

    // Sheet folding:
	Structure::Sheet targetCopy (*((Structure::Sheet*)tn));
    Array2D_Vector3 deltas = targetCopy.foldTo( Array1D_Vector4d(1, coord), (this->type == GROW) );
    if (this->type != GROW) deltas = inverseVectors3(deltas);

    // Growing / shrinking instructions
    property["deltas"].setValue( deltas );
    property["orgCtrlPoints"].setValue( structure_sheet->surface.mCtrlPoint );
}

void TaskSheet::prepareSheetTwoEdges( Structure::Link * linkA, Structure::Link * linkB )
{
    Structure::Node * n = node();
	Structure::Node * tn = targetNode();

	Structure::Link * tlinkA = target->getEdge( linkA->property["correspond"].toInt() );
	Structure::Link * tlinkB = target->getEdge( linkB->property["correspond"].toInt() );

    Vector3d pointA = linkA->positionOther(n->id);
    Vector3d pointB = linkB->positionOther(n->id);

    // Geodesic distance between two link positions on the active graph excluding the running tasks
    QVector<QString> exclude = active->property["activeTasks"].value< QVector<QString> >();
    GraphDistance gd( active, exclude );
    gd.computeDistances( pointA, DIST_RESOLUTION );
    QVector< GraphDistance::PathPointPair > path;
    gd.smoothPathCoordTo(pointB, path);

    // Otherwise, self expand / contract
    if(path.size() < 3)
    {
        QVector<Link*> twoEdges(2);
        twoEdges[0] = linkA; twoEdges[1] = linkB;
		Structure::Link * pickedLink = preferredEnd(n, twoEdges, active);
        prepareSheetOneEdge( pickedLink );

		QVector<Link*> edges; edges.push_back(pickedLink);
		property["edges"].setValue( active->getEdgeIDs(edges) );
        return;
    }

    // Use the center of the path as the start point
    GraphDistance::PathPointPair startPointCoord = path[path.size() / 2];
    Vector3d startPoint = startPointCoord.position( active );

    // Separate the path into two for linkA and linkB
    int N = path.size(), hN = N / 2;
    if (N %2 == 0) path.insert(hN, path[hN]);

    QVector<GraphDistance::PathPointPair> pathA, pathB;
    for (int i = 0; i < hN; i++)
    {
        pathA.push_back(path[hN+1+i]);
        pathB.push_back(path[hN-1-i]);
    }

    // Record path
    property["pathA"].setValue( pathA );
    property["pathB"].setValue( pathB );

    // Encode sheet on a line segment
	RMF rmf( GraphDistance::positionalPath(active, pathA, 1) );
	property["rmf"].setValue( rmf );
	if(!rmf.count()) return;

	Vector3 X = rmf.U.back().r, Y = rmf.U.back().s, Z = rmf.U.back().t;

    SheetEncoding cpCoords = encodeSheetAsCurve((Structure::Sheet*)tn, tlinkA->position(tn->id), tlinkB->position(tn->id));
	property["cpCoords"].setValue( cpCoords );

    // DEBUG
    node()->property["rmf"].setValue( rmf );
    node()->property["rmf2"].setValue( RMF ( GraphDistance::positionalPath(active, pathB, 3) ) );

    if(this->type == GROW)
    {
        // Initial position and geometry
        n->setControlPoints( Array1D_Vector3(cpCoords.size(), startPoint) );
    }
}

void TaskSheet::prepareGrowShrinkSheet()
{
    Structure::Node * n = node();
    QVector<Structure::Link*> edges = filterEdges(n, active->getEdges(n->id));

	property["edges"].setValue( active->getEdgeIDs(edges) );

    if (edges.size() == 1)
    {
        prepareSheetOneEdge( edges.front() );
    }

    if (edges.size() == 2)
    {
        prepareSheetTwoEdges( edges.front(), edges.back() );
    }
}

void TaskSheet::executeGrowShrinkSheet(double t)
{
	Structure::Sheet* structure_sheet = ((Structure::Sheet*)node());
	QVector<Link*> edges = active->getEdges( property["edges"].value< QVector<int> >() );

	/// Single edge case
	if ( property.contains("deltas") )
	{
		Array2D_Vector3 cpts = property["orgCtrlPoints"].value<Array2D_Vector3>();
		Array2D_Vector3 deltas = property["deltas"].value<Array2D_Vector3>();

		// Grow sheet
		for(size_t u = 0; u < structure_sheet->surface.mNumUCtrlPoints; u++)
			for(size_t v = 0; v < structure_sheet->surface.mNumVCtrlPoints; v++)
				structure_sheet->surface.mCtrlPoint[u][v] = cpts[u][v] + (deltas[u][v] * t);
	}

	/// Two edges case
	if( property.contains("pathA") && property.contains("pathB") && property.contains("cpCoords") )
	{
		QVector< GraphDistance::PathPointPair > pathA = property["pathA"].value< QVector< GraphDistance::PathPointPair > >();
		QVector< GraphDistance::PathPointPair > pathB = property["pathB"].value< QVector< GraphDistance::PathPointPair > >();
		if(pathA.size() == 0 || pathB.size() == 0)	return;

		double dt = t;
		if(type == SHRINK) dt = 1 - t;

		double decodeT = qMin(1.0, dt * 2.0);

		int idxA = dt * (pathA.size() - 1);
		int idxB = dt * (pathB.size() - 1);

		// Move to next step
		Vector3 pointA = pathA[idxA].position(active);
		Vector3 pointB = pathB[idxB].position(active);

		Structure::Link *linkA = edges.front(), *linkB = edges.back();
		if (type == GROW){
			linkA = target->getEdge(linkA->property["correspond"].toInt());
			linkB = target->getEdge(linkB->property["correspond"].toInt());
		}

		Vector3d deltaA = linkA->property["delta"].value<Vector3d>() * dt;
		Vector3d deltaB = linkB->property["delta"].value<Vector3d>() * dt;

		Array1D_Vector3 decoded = Curve::decodeCurve(property["cpCoords"].value<CurveEncoding>(), pointA + deltaA, pointB + deltaB, decodeT);
		structure_sheet->setControlPoints( decoded );
	}
}


void TaskSheet::prepareMorphSheet()
{
	Structure::Node * n = node();
	Structure::Node * tn = targetNode();

	// Check if the sheet is crossing 
	{
		QVector<Structure::Link*> edges = filterEdges(n, active->getEdges(n->id));

		// check if neighbors on selected edges will change
		foreach(Link * l, edges){
			if(l->property["modified"].toBool()) continue;
			Structure::Link* tl = target->getEdge(l->property["correspond"].toInt());
			if (l->otherNode(n->id)->property["correspond"].toString() != tl->otherNode(tn->id)->id)
				property["isCrossing"] = true;
		}
	}

	// Morph parameters
	{
		Structure::Sheet * sheet = (Structure::Sheet *)n;
		Structure::Sheet * tsheet = (Structure::Sheet *)tn;

		// Get source and target frames
		RMF::Frame sframe = sheetFrame( sheet );
		RMF::Frame tframe = sheetFrame( tsheet );

		// Compute rotations between source and target sheet frames
		Eigen::Quaterniond rotation, eye = Eigen::Quaterniond::Identity();
		AbsoluteOrientation::compute(sframe.r,sframe.s,sframe.t,
			tframe.r,tframe.s,tframe.t, rotation);

		// Parameters needed for morphing
		QVector<double> rot; 
		rot.push_back(rotation.w());
		rot.push_back(rotation.x());
		rot.push_back(rotation.y());
		rot.push_back(rotation.z());

		property["rotation"].setValue( rot );
		property["sframe"].setValue( sframe );
		property["tframe"].setValue( tframe );

		property["cpCoords"].setValue( Structure::Sheet::encodeSheet(sheet, sframe.center, sframe.r,sframe.s,sframe.t) );
		property["cpCoordsT"].setValue( Structure::Sheet::encodeSheet(tsheet, tframe.center, tframe.r,tframe.s,tframe.t) );

		// Visualization
		std::vector<RMF::Frame> frames = smoothRotateFrame(sframe, rotation, 100);
		property["frames"].setValue( frames );
		n->property["frames"].setValue( frames );
		n->property["frame"].setValue( sframe );
		tn->property["frame"].setValue( tframe );
	}

	if ( property["isCrossing"].toBool() ) 
		prepareCrossingSheet();
	
	prepareMorphEdges();
}

void TaskSheet::executeMorphSheet( double t )
{
    Structure::Node * n = node();
    Structure::Sheet * sheet = (Structure::Sheet *)n;

    // Decode
    SheetEncoding cpCoords = property["cpCoords"].value<SheetEncoding>();
    SheetEncoding cpCoordsT = property["cpCoordsT"].value<SheetEncoding>();

    RMF::Frame sframe = property["sframe"].value<RMF::Frame>();
    RMF::Frame tframe = property["tframe"].value<RMF::Frame>();
	QVector<double> rot = property["rotation"].value< QVector<double> >();
    Eigen::Quaterniond rotation(rot[0],rot[1],rot[2],rot[3]),
        eye = Eigen::Quaterniond::Identity();

    // Source sheet
    Eigen::Vector3d R = (sframe.r), S = (sframe.s), T = (sframe.t);
    R = eye.slerp(t, rotation) * R;
    S = eye.slerp(t, rotation) * S;
    T = eye.slerp(t, rotation) * T;
    RMF::Frame curFrame ((R),(S),(T));

    curFrame.center = tframe.center = AlphaBlend(t, sframe.center, tframe.center);

    Array1D_Vector3 newPnts = Sheet::decodeSheet( cpCoords, curFrame.center, curFrame.r, curFrame.s, curFrame.t );
    Array1D_Vector3 newPntsT = Sheet::decodeSheet( cpCoordsT, tframe.center, tframe.r, tframe.s, tframe.t );

    Array1D_Vector3 blendedPnts;

    for(int i = 0; i < (int)newPnts.size(); i++)
    {
        blendedPnts.push_back( AlphaBlend(t, newPnts[i], newPntsT[i]) );
    }

    sheet->setControlPoints( blendedPnts );
}

void TaskSheet::encodeSheet( const Vector4d& coordinateA, const Vector4d& coordinateB )
{
	Node * n = node(), *tn = targetNode();

	//// Parameters needed for morphing
	Vector3 sourceA = n->position(coordinateA), sourceB = n->position(coordinateB);
	Vector3 targetA = tn->position(coordinateA), targetB = tn->position(coordinateB);

	// Two ends
	property["sourceA"].setValue( sourceA ); property["sourceB"].setValue( sourceB );
	property["targetA"].setValue( targetA ); property["targetB"].setValue( targetB );

	// Encoding
	property["cpCoords"].setValue( encodeSheetAsCurve((Structure::Sheet*)n, sourceA, sourceB) );
	property["cpCoordsT"].setValue( encodeSheetAsCurve((Structure::Sheet*)tn, targetA, targetB) );
}

void TaskSheet::prepareCrossingSheet()
{
	Structure::Node * n = node();

	//QVector<Structure::Link*> edges = filterEdges(n, active->getEdges(n->id));

	// For now, let us only use one edge:
	QVector<Structure::Link*> edges;
	edges.push_back( filterEdges(n, active->getEdges(n->id)).front() );

	if (edges.size() == 1)
	{
		// Start and end
		Link * link = edges.front();
		Vector3d start = link->positionOther(n->id);
		NodeCoord futureNodeCord = futureLinkCoord(link);
		Vector3d end = active->position(futureNodeCord.first,futureNodeCord.second);

		// Geodesic distances on the active graph excluding the running tasks
		QVector< GraphDistance::PathPointPair > path;
		QVector<QString> exclude = active->property["activeTasks"].value< QVector<QString> >();
		GraphDistance gd( active, exclude );
		gd.computeDistances( end, DIST_RESOLUTION );  
		gd.smoothPathCoordTo( start, path );

		// Check
		if(path.size() < 2) { path.clear(); path.push_back(GraphDistance::PathPointPair( PathPoint(futureNodeCord.first, futureNodeCord.second))); }

		property["path"].setValue( GraphDistance::positionalPath(active, path) );

		// Save links paths
		path.back() = GraphDistance::PathPointPair( PathPoint(futureNodeCord.first, futureNodeCord.second)  );
		link->setProperty("path", path);
	}

	property["edges"].setValue( active->getEdgeIDs(edges) );
}

void TaskSheet::executeCrossingSheet( double t )
{
	Node *n = node(), *tn = targetNode();

	QVector<Link*> edges = active->getEdges( property["edges"].value< QVector<int> >() );

	if (property.contains("path"))
	{
		Vector3 p0 = n->position(Vec4d(0,0,0,0));
		Vector3 delta(0,0,0);

		if(!edges.isEmpty())
		{
			p0 = edges.front()->position(n->id);

			// Blend Deltas
			Structure::Link *slink = edges.front();
			Structure::Link* tlink = target->getEdge(slink->property["correspond"].toInt());
			Vector3 d1 = slink->position(n->id) - slink->positionOther(n->id);
			Vector3 d2 = tlink->position(tn->id) - tlink->positionOther(tn->id);
			delta = AlphaBlend(t, d1, d2);
		}

		executeMorphSheet(t);
		
		// Cancel any absolute movement
		if(!edges.isEmpty())
			n->moveBy(p0 - edges.front()->position(n->id));

		// Move it to the correct position
		Array1D_Vector3 path = property["path"].value< Array1D_Vector3 >();
		int idx = t * (path.size() - 1);

		Vector3 oldPosOnMe = p0;
		Vector3 newPosOnMe = path[idx] + delta;

		n->moveBy(newPosOnMe - p0);
	}
}

SheetEncoding TaskSheet::encodeSheetAsCurve( Structure::Sheet * sheet, Vector3 start, Vector3 end)
{
    Array1D_Vector3 controlPoints = sheet->controlPoints();
    return Structure::Curve::encodeCurve(controlPoints,start,end);
}

Structure::Sheet * TaskSheet::targetSheet()
{
    Structure::Node* n = targetNode();
    if(!n || n->type() != Structure::SHEET) return NULL;
    return (Structure::Sheet *)n;
}
