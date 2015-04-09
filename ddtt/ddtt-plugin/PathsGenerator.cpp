#pragma warning(disable:4267)

#include <random>
#include <algorithm>
#include <iterator>

#include "SymmetryAnalysis.h"
Q_DECLARE_METATYPE( QVector<SymmetryAnalysis*> )

#include "PathsGenerator.h"

#include <QApplication>
#include <QGLWidget>
#include <QGLFramebufferObjectFormat>
#include <QTableWidget>
#include <QHeaderView>
#include <QDebug>
#include "GenericGraph.h"
#include "StructureGraph.h"
#include "GraphDistance.h"
using namespace Structure;

// Help generate combinations
#include "next_combination.h"
#include "cartesian.h"

#include "Hungarian.h"
#undef max
#undef min

static inline QString shortName(QString name){
    if(name.length() < 3) return name;
    return QString("%1%2%3").arg(name.at(0)).arg(name.at(1)).arg(name.at(name.length()-1));
}

QStringList toQStringList( const QVector<QString> & v ){
    QStringList l;
    for(auto & s : v) l << s;
    return l;
}

void PathsGenerator::measure_position(Structure::Graph * graph)
{
    Eigen::AlignedBox3d bbox = graph->bbox();
    Vector3 bbox_min = bbox.min();
    Vector3 bbox_max = bbox.max();

    foreach(Node * n, graph->nodes)
    {
        if(n->property.contains("taskType") && n->property["taskType"].toInt() == Task::GROW) continue;

        Vector3 p = n->bbox().center();
        std::vector<Vector3> curve = n->discretizedAsCurve( DIST_RESOLUTION );

        n->meta["geo_x"] = (p[0] - bbox_min[0]) / (bbox_max[0] - bbox_min[0]);
        n->meta["geo_y"] = (p[1] - bbox_min[1]) / (bbox_max[1] - bbox_min[1]);
        n->meta["geo_z"] = (p[2] - bbox_min[2]) / (bbox_max[2] - bbox_min[2]);
    }
}

void PathsGenerator::measure_ground_parts(Structure::Graph * graph)
{
    /// Find set of ground nodes
    std::vector<Structure::Node*> floorNodes;
    {
        Eigen::AlignedBox3d box = graph->bbox( true );
        double bottom = box.min().z();
        double height = box.max().z() - bottom;

        double threshold = height * 0.1;

        foreach(Node * n, graph->nodes)
        {
            double minz = n->bbox().min().z();
            double dist = abs(minz - bottom);

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

    /// Assign based on distance to ground
    {
        QMap<Node*, double> dists;

        foreach(Node * n, graph->nodes)
            dists[n] = g.NodeDistance(floorNode, n->property["index"].toUInt());

        QList<double> dists_sorted = dists.values();
        qSort(dists_sorted);

        foreach(Node * n, graph->nodes)
        {
            if(n->property.contains("taskType") && n->property["taskType"].toInt() == Task::GROW) continue;

            std::vector<Vector3> curve = n->discretizedAsCurve( DIST_RESOLUTION );

            double range = (dists_sorted.back() - dists_sorted.front());
            double val = (dists[n] - dists_sorted.front()) / range;

            if(range == 0.0) val = 0;

            n->meta["structure_ground"] = val;
        }
    }
}

void PathsGenerator::compute_part_measures(Structure::Graph * graph)
{
    //measure_ground_parts( graph );
    measure_position( graph );
}

double PathsGenerator::partDifference(QString sid, QString tid, Structure::Graph * source, Structure::Graph * target )
{
    std::vector<double> sim( target->nodes.size(), DBL_MAX );

    Structure::Node * sn = source->getNode(sid);
    Structure::Node * tn = target->getNode(tid);

    int k = sn->meta.size();
    QList<QString> keys = sn->meta.keys();
    QList<QString> tkeys = tn->meta.keys();
    keys = (keys.toSet().intersect(tkeys.toSet())).toList(); // make sure only test shared attributes

    // Source vector
    Eigen::VectorXd Vs( k ), Vt ( k );

    for(int i = 0; i < k; i++){
        Vs(i) = sn->meta.value(keys.at(i)).toDouble();
        Vt(i) = tn->meta.value(keys.at(i)).toDouble();
    }

    return (Vs - Vt).norm();
}

double PathsGenerator::partDifference2(QString sid, QString tid, Structure::Graph * source, Structure::Graph * target )
{
    Array1D_Vector3 sPoints = source->getNode(sid)->controlPoints();
    Array1D_Vector3 tPoints = target->getNode(tid)->controlPoints();
    return HausdorffDistance(sPoints, tPoints);
}

std::vector<std::vector<double> > PathsGenerator::buildDifferenceMatrix( Structure::Graph * source, Structure::Graph * target )
{
    int N = source->nodes.size();
    int M = target->nodes.size();

    int extra_N = (M > N) ? M-N : 0;
    int extra_M = (N > M) ? N-M : 0;

    std::vector<std::vector<double> > m( N + extra_N, std::vector<double>(M + extra_M, DBL_MAX) );

    double minVal = DBL_MAX, maxVal = -DBL_MAX;

    for(int i = 0; i < N; i++){
        for(int j = 0; j < M; j++){
            double val = partDifference2( source->nodes[i]->id, target->nodes[j]->id, source, target );
            m[i][j] = val;

            // Track limits
            minVal = std::min(minVal, val);
            maxVal = std::max(maxVal, val);
        }
    }

    // Normalize
    for(int i = 0; i < (int)m.size(); i++){
        for(int j = 0; j < (int)m.front().size(); j++){
            if(m[i][j] == DBL_MAX) m[i][j] = maxVal;
            m[i][j] = (m[i][j] - minVal) / (maxVal - minVal);
        }
    }

    return m;
}

void PathsGenerator::visualizeDifferenceMatrix(std::vector< std::vector<double> > m, Structure::Graph * sg, Structure::Graph * tg)
{
    int rows = sg->nodes.size(), cols = tg->nodes.size();
    QTableWidget * tw = new QTableWidget(rows, cols);

    // Cells
    int fixedSize = 40;
    tw->horizontalHeader()->setDefaultSectionSize(fixedSize);
    tw->verticalHeader()->setDefaultSectionSize(fixedSize);
    tw->horizontalHeader()->setSectionResizeMode(QHeaderView::Fixed);
    tw->verticalHeader()->setSectionResizeMode(QHeaderView::Fixed);

    // Labels
    QStringList hlabels, vlabels;
    for(auto n : sg->nodes) vlabels << shortName(n->id); tw->setVerticalHeaderLabels(vlabels);
    for(auto n : tg->nodes) hlabels << shortName(n->id); tw->setHorizontalHeaderLabels(hlabels);

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            double val = abs( m[i][j] );
            tw->setItem(i, j, new QTableWidgetItem( QString::number(val, 'g', 1) ));
            QColor c = starlab::qtJetColor( val );
            tw->item(i,j)->setBackground( starlab::qtJetColor( val ) );
            tw->item(i,j)->setTextColor( c.lighter() );
        }
    }

	tw->setMinimumSize(900,900);
    tw->show();
    qApp->processEvents();
}

VectorPairStrings PathsGenerator::pairsDebugging( QVector<Pairing> pairs )
{
    VectorPairStrings result;
    for(auto p : pairs) result.push_back( qMakePair(toQStringList(p.first).join(","), toQStringList(p.second).join(",")) );
    return result;
}

QVector<QString> PathsGenerator::nothingSet(){
    return QVector<QString>() << "NOTHING";
}

VectorPairings PathsGenerator::splitMany( QVector<QString> A, QVector<QString> B )
{
    VectorPairings result;

    // Pick first and last from smaller set
    QPair< QString, QString > smaller( A.front(), A.back() );

    // Divide larger into two
    QPair< QVector<QString>, QVector<QString> > larger;
    std::vector<QString> b1(B.begin(), B.begin() + B.size()/2);
    std::vector<QString> b2(B.begin() + B.size()/2, B.end());
    for(auto n : b1) larger.first.push_back(n);
    for(auto n : b2) larger.second.push_back(n);

    result.push_back( qMakePair(QVector<QString>() << smaller.first, larger.first) );
    result.push_back( qMakePair(QVector<QString>() << smaller.second, larger.second) );

    return result;
}

void PathsGenerator::prepareSymmetryAnalysis(Structure::Graph * sg, Structure::Graph * tg)
{
	QVector< QVector<SymmetryAnalysis*> > analysis;
	QVector<Structure::Graph*> graphs; graphs << sg << tg;

	for(auto g : graphs)
	{
		QVector<SymmetryAnalysis*> sanalysis;

		for(int i = 0; i < g->groups.size(); i++)
		{
			auto group = g->groups.at(i);
			std::vector<SurfaceMesh::SurfaceMeshModel*> group_meshes;
			for( auto element : group ) group_meshes.push_back( g->getMesh(element) );
			sanalysis.push_back( new SymmetryAnalysis( group_meshes ) );
		}

		analysis.push_back( sanalysis );
	}

	sg->property["symmetryAnalysis"].setValue( analysis.front() );
	tg->property["symmetryAnalysis"].setValue( analysis.back() );
}

QVector<VectorPairings> PathsGenerator::manyToMany2(Structure::Graph * sg, Structure::Graph * tg, QVector<QString> snodes, QVector<QString> tnodes)
{
	QVector<VectorPairings> allResolved;

	SymmetryAnalysis *sanalysis = NULL, *tanalysis = NULL;
	for(auto a : sg->property["symmetryAnalysis"].value< QVector<SymmetryAnalysis*> >()){
		if(a->hasMesh(snodes.front())){
			sanalysis = a;
			break;
		}
	}
	for(auto a : tg->property["symmetryAnalysis"].value< QVector<SymmetryAnalysis*> >()){
		if(a->hasMesh(tnodes.front())){
			tanalysis = a;
			break;
		}
	}
	
	assert(sanalysis && tanalysis);

	bool isSwapped = false;
	if( tnodes.size() > snodes.size() )
	{
		std::swap(sg, tg);
		std::swap(snodes, tnodes);
		std::swap(sanalysis, tanalysis);
		isSwapped = true;
	}

	VectorPairings resolved;

	for(size_t i = 0; i < snodes.size(); i++)
	{
		QString sid = snodes[i], tid;
		double minDist = DBL_MAX;
		Vector3 sCenter = sanalysis->centers.row(i);

		for(size_t j = 0; j < tnodes.size(); j++)
		{
			double dist = (sCenter - tanalysis->centers.row(j).transpose()).norm();
			if(dist < minDist){
				minDist = dist;
				tid = tnodes[j];
			}
		}

		if( isSwapped ) std::swap( sid, tid );
		
		resolved.push_back( Pairing( QVector<QString>()<< sid, QVector<QString>()<< tid ) );
	}

	allResolved.push_back( resolved );

	return allResolved;
}

QVector<VectorPairings> PathsGenerator::manyToMany(Structure::Graph * sg, Structure::Graph * tg, QVector<QString> snodes, QVector<QString> tnodes)
{
    QVector<VectorPairings> allResolved;

    // Sort based on 'x' axis
    QVector< QPair<double, Structure::Node*> > sortedS, sortedT;
    for(auto nid : snodes) {
        Structure::Node * n = sg->getNode(nid);
        sortedS.push_back( qMakePair(n->meta.contains("geo_x") ? n->meta["geo_x"].toDouble() : sortedS.size(), n) );
    }
    for(auto nid : tnodes) {
        Structure::Node * n = tg->getNode(nid);
        sortedT.push_back( qMakePair(n->meta.contains("geo_x") ? n->meta["geo_x"].toDouble() : sortedT.size(), n) );
    }
    qSort(sortedS);
    qSort(sortedT);

    bool isReversed = false;

    QVector<QString> sortedA, sortedB;
    for(auto nv : sortedS) sortedA.push_back(nv.second->id);
    for(auto nv : sortedT) sortedB.push_back(nv.second->id);

    if(sortedS.size() > sortedT.size()){
        isReversed = true;
        std::swap(sortedA, sortedB);
    }

    // Match using 'splitMany'
    {
        VectorPairings resolved;
        VectorPairings split = splitMany(sortedA, sortedB);
        if( isReversed ) for(auto & s : split) std::swap(s.first, s.second);
        for(auto p : split) resolved.push_back(p);
        allResolved.push_back(resolved);
    }

    // Try a one-to-one
    if(sortedA.size() == sortedB.size())
    {
        VectorPairings resolved;

        for(auto sid: sortedA){
            //Array1D_Vector3 sPoints = sg->getNode(sid)->controlPoints();

            Vector3 ps = sg->getNode(sid)->controlPoints().front();

            double minDist = DBL_MAX;
            QString bestTID;

            for(auto tid: sortedB){
                //Array1D_Vector3 tPoints = tg->getNode(tid)->controlPoints();
                //double dist = HausdorffDistance(sPoints, tPoints);
                Vector3 pt = tg->getNode(tid)->controlPoints().front();
                double dist = (ps - pt).norm();

                // Closest
                if(dist < minDist){
                    minDist = dist;
                    bestTID = tid;
                }
            }

            resolved.push_back( Pairing( QVector<QString>()<< sid, QVector<QString>()<< bestTID ) );
        }

        allResolved.push_back( resolved );
    }

    return allResolved;
}

Assignments PathsGenerator::allAssignments( QVector< QVector<QString> > sgroups, QVector< QVector<QString> > tgroups,
                           mat m, Structure::Graph * sg, Structure::Graph * tg, int K)
{
    Assignments assignments;

    // Map to index
    std::map< int, QVector<QString> > sourceItem, targetItem;
    for(auto items : sgroups) sourceItem[sourceItem.size()] = items;
    for(auto items : tgroups) targetItem[targetItem.size()] = items;
    sourceItem[sourceItem.size()] = nothingSet();
    targetItem[targetItem.size()] = nothingSet();

    // Bounds
    K = qMin(K, tgroups.size());

    // Record all initial costs
    std::map< int, QVector< QPair<double, int> > > candidates;
    std::map< QString, QString > candidates_debug;
    for(size_t i = 0; i < sgroups.size(); i++)
    {
        QVector< QPair<double, int> > diffs;

        // Score based on minimum matching of the source and target group
        for(size_t jj = 0; jj < tgroups.size(); jj++)
        {
            double minVal = DBL_MAX;
            int si = 0;
            int tj = 0;

            for(auto sid : sgroups[i])
            {
                int i = sg->indexOfNode( sg->getNode(sid) );

                for(auto tid : tgroups[jj])
                {
                    int j = tg->indexOfNode( tg->getNode(tid) );
                    double val = abs(m[i][j]);

                    if(val < minVal){
                        si = i;
                        tj = j;
                        minVal = val;
                    }
                }
            }

            diffs.push_back( qMakePair( m[si][tj], jj ) );
        }

        std::sort(diffs.begin(), diffs.end());
        //std::reverse(diffs.begin(), diffs.end());
        candidates[i] = diffs;

        // DEBUG:
        QString curGroup = toQStringList(sgroups[i]).join("");
        for(size_t j = 0; j < tgroups.size(); j++){
            QString targetGroup = tgroups[diffs[j].second].front();
            candidates_debug[ curGroup ] += targetGroup + "-";
        }
    }

    // Collect only 'k' nearest candidates + nothing
    std::map< int, std::vector<int> > allowed;
    for(size_t i = 0; i < sgroups.size(); i++){
        allowed[i] = std::vector<int>();
        for(int j = 0; j < K; j++) allowed[i].push_back(candidates[i][j].second); // Good candidates
        allowed[i].push_back( tgroups.size() ); // Nothing
    }

    // DEBUG:
    std::map<QString, QString> allowed_debug;
    for(size_t i = 0; i < sgroups.size(); i++){
        for(int j = 0; j < K; j++)
            allowed_debug[ toQStringList(sgroups[i]).join("") ] += tgroups[ allowed[i][j] ].front() + "-";
    }

    // Build full set of good assignments
    std::vector< std::vector<int> > goodCombs;
    {
        // Flatten sets
        std::vector< std::vector<int> > setAllowed;

        for(auto a : allowed){
            std::vector<int> s;
            for(auto e : a.second) s.push_back(e);
            setAllowed.push_back(s);
        }

        int NOTHING_ITEM = tgroups.size();

        blitz::CartesianProduct< std::vector<int>, std::vector<int> > product( setAllowed );

        // Go over all possible combinations, only keep good ones (i.e. single matching to targets)
        int totalCombs = 0;
        for(auto comb : product){
            totalCombs++;
            bool isGood = true;
            std::set<int> titems;
            for(auto i : comb) isGood = isGood && ((i == NOTHING_ITEM) || titems.insert(i).second);
            if( isGood ) goodCombs.push_back( comb );
        }
    }

    // Generate final assignments
    for(auto assignment : goodCombs)
    {
        VectorPairings paring;
        QVector<Pairing> manyManyCases;

        for(size_t i = 0; i < assignment.size(); i++)
        {
            QVector<QString> snodes = sourceItem[i];
            QVector<QString> tnodes = targetItem[assignment[i]];

            // Add all but many-to-many
            if(snodes.size() > 1 && tnodes.size() > 1)
                manyManyCases.push_back( qMakePair(snodes, tnodes) );
            else
                paring.push_back( qMakePair(snodes,tnodes) );
        }

        // Now resolve many-to-many using different possibilities
        if( manyManyCases.size() )
        {
            QVector< QVector<VectorPairings> > manyManyResolution;

            for(auto candidate : manyManyCases)
                manyManyResolution.push_back( manyToMany2(sg, tg, candidate.first, candidate.second) );

            blitz::CartesianProduct< QVector<VectorPairings>, QVector<VectorPairings> > product( manyManyResolution );

            for(auto comb : product){
                VectorPairings curParing = paring;
                for(auto pairs : comb)
                    for(auto pair : pairs)
                        curParing.push_back(pair);

                assignments.push_back( curParing );
            }
        }
        else
            assignments.push_back( paring );
    }

    return assignments;
}

PathsGenerator::PathsGenerator(Structure::Graph * source, Structure::Graph * target, int K)
{
    // Generate set of random colors
    QVector<QColor> colors;
    for(int i = 0; i < (int)source->nodes.size() * 2; i++) colors.push_back(starlab::qRandomColor2());
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(colors.begin(), colors.end(), g);

    QVector<Structure::Graph*> graphs; graphs << source << target;

    /// Graphs preprocessing
    for(auto g : graphs) compute_part_measures( g );
    mat m = buildDifferenceMatrix(source, target);

    // Visualization
    if(true) visualizeDifferenceMatrix(m, source, target);

	// Prepare symmetry analysis
	prepareSymmetryAnalysis( source, target );

    // Get all assignments
    Assignments assignments = allAssignments(source->nodesAsGroups(), target->nodesAsGroups(), m, source, target, K);

	// Check against Ground-truth, did we get the good assignment?
	VectorPairings groundTruth;
    {
		// Check for existing correspondence file
		QString path = QFileInfo(source->property["name"].toString()).absolutePath() + "/", ext = ".txt";
		QString g1n = source->name(), g2n = target->name();

		QStringList files; 
		files << (path+g1n+"_"+g2n+ext) << (path+g2n+"_"+g1n+ext);
		files << source->property["correspondenceFile"].toString();

		int fileidx = -1;
		for(int i = 0; i < files.size(); i++){
			QString f = files[i];
			if(!QFileInfo (f).exists()) continue;
			fileidx = i;
		}
		if( fileidx != -1 ){

			bool corrReversed = (fileidx == 1) ? true : false;
			GraphCorresponder * bestCorrespond = new GraphCorresponder(source, target);
			bestCorrespond->loadCorrespondences(files[fileidx], corrReversed);

			QVector< PART_LANDMARK > vp = bestCorrespond->correspondences;
			for(auto p : vp) groundTruth.push_back( Pairing(p.first, p.second) );
		}
    }

    // Prepare deformation paths
    for( auto a : assignments )
    {
        int pathIndex = paths.size();

        paths.push_back( DeformationPath() );
        DeformationPath & path = paths.back();

        path.pairs = a;
        path.pairsDebug = pairsDebugging( path.pairs );

		Structure::Graph * sourceCopy = new Structure::Graph(*source);
		Structure::Graph * targetCopy = new Structure::Graph(*target);

        path.gcorr = new GraphCorresponder(sourceCopy, targetCopy);

        int i = 0;

        for(auto p : path.pairs)
        {
            path.gcorr->addCorrespondences( p.first, p.second, -1 );

            // Nodes coloring
            if( !p.first.contains("NOTHING") && !p.second.contains("NOTHING") )
            {
                for(auto sid : p.first) path.scolors[sid] = colors[i];
                for(auto tid : p.second) path.tcolors[tid] = colors[i];

                i++;
            }
        }

        path.gcorr->isReady = true;
        path.gcorr->correspondAllNodes();

        path.i = pathIndex;

		// Check for ground truth
		VectorPairings v1 = path.pairs;
		VectorPairings v2 = groundTruth;
		if(v1.size() == v2.size() && std::is_permutation(v1.begin(), v1.end(), v2.begin())){
			path.property["isGroundTruth"].setValue( true );
		}
    }
}

