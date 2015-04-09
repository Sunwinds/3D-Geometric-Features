#include "Corresponder.h"
#include "StructureGraph.h"
#include "RenderObjectExt.h"
#include "GraphExplorer.h"
#include "GraphDistance.h"

inline void bounds(double val, double & maximum, double & minimum){
	maximum = qMax(maximum, val);
	minimum = qMin(minimum, val);
}

template<class T>
inline T general_mean(const std::vector<T> & data, T default_val = T()){
	T center = default_val;
	for(auto t : data) center += t;
	return center /= data.size();
}

inline Vector3 mean( const std::vector<Vector3> & data ){
	return general_mean(data, Vector3(0,0,0));
}

double length(const std::vector<Eigen::Vector3d> & data){
	double l = 0.0;
	for(size_t i = 0; i + 1 < data.size(); i++) l += (data[i+1] - data[i]).norm();
	return l;
}

Eigen::MatrixXd pointsToMatrix( const std::vector<Eigen::Vector3d> & data, bool isCols = true ){
	Eigen::MatrixXd m(data.front().size(), data.size());
	for(size_t i = 0; i < data.size(); i++){
		if(isCols) m.col(i) = data[i];
		else m.row(i) = data[i];
	}
	return m;
}

inline Eigen::Vector3d noise(double s = 1.0){
	return Eigen::Vector3d(starlab::uniformRand(), starlab::uniformRand(), starlab::uniformRand()) * s;
}

std::vector<Eigen::Vector3d> withNoise(const std::vector<Eigen::Vector3d> & data, double s = 1.0){
	std::vector<Eigen::Vector3d> result;
	for(auto p : data) result.push_back(p + noise(s));
	return result;
}

QString matToString(Eigen::MatrixXd m, QString name = ""){
	QStringList s;
	for(int i = 0; i < m.rows(); i++){
		QStringList r;
		for(int j = 0; j < m.cols(); j++)
			r << QString::number(m(i,j));
		s << r.join("\t");
	}
	return  "\n\n" + name + "\n\n" + s.join("\n") + "\n";
}

template<typename T>
std::vector< std::vector<T> > chunks(const std::vector<T> & data, int num_chunks){
	std::vector< std::vector<T> > chunks;
	if(num_chunks <= 1) { chunks.push_back(data); return chunks; }
	size_t k = data.size() / num_chunks;
	std::vector<T>::const_iterator chunk_begin = data.begin(), chunk_end = chunk_begin;
	std::vector<T>::const_iterator end = data.end();
	do{
		if(std::distance(chunk_end, end) < (int)k) chunk_end = end;
		else std::advance(chunk_end, k);
		chunks.push_back( std::vector<T>(chunk_begin,chunk_end) );
		chunk_begin = chunk_end;
	} while(std::distance(chunk_begin,end) > 0);
	return chunks;
}

Corresponder::Corresponder(Structure::Graph * g1, Structure::Graph * g2, int nidx, Metric metric, int splitTo)
	: source(g1), target(g2), selected_nidx(nidx), selected_metric(metric), num_splits(splitTo)
{
	n = m = 0;
}

Eigen::MatrixXd Corresponder::geometricDistance()
{
	Eigen::MatrixXd g = Eigen::MatrixXd::Zero(n, m);
	return g;
}

Eigen::MatrixXd Corresponder::structuralDistance()
{
	Eigen::MatrixXd s = Eigen::MatrixXd::Zero(n, m);

	return s;
}

void Corresponder::compute( PropertyMap prop )
{
	double RES = source->bbox().sizes().norm() * 0.025;

	QVector<Structure::Graph *> graphs; graphs << source << target;

	std::vector< std::vector< std::vector<Vector3> > > skeletons;
	std::vector< std::vector< double > > closeness;
	std::vector< std::vector< Vector3 > > centerCoord;
	std::vector< std::vector< double > > lengths;
	std::vector< std::vector< double > > agdVals;

	int numSkelSamples = 20;
	int num_chunks = qMax(1, num_splits & ~1);

	foreach(Structure::Graph * g, graphs)
	{
		/// Sample the skeletons
		{
			skeletons.push_back(std::vector< std::vector<Vector3> >());

			for(size_t i = 0; i < g->nodes.size(); i++)
			{
				std::vector<Vector3> samples;

				Structure::Node * n = g->nodes[int(i)];

				if( n->type() == Structure::CURVE ){
					Structure::Curve * c = (Structure::Curve *)n;
					c->curve.SubdivideByLength(numSkelSamples, samples);
				}
				if( n->type() == Structure::SHEET ){
					samples = refineByNumber(n->discretizedAsCurve(RES), numSkelSamples);
				}

				// Split skeletons, add the chunks
				std::vector< std::vector<Vector3> > subSkeletons = chunks(samples, num_chunks);
				for(auto chunk : subSkeletons) skeletons.back().push_back(chunk);
			}
		}

		/// Compute graph structure metrics
		{
			// Closeness
			/*
			{
				g->computeAllCloseness();

				// Normalize
				bool isNormalizeCloseness = true;
				if( isNormalizeCloseness )
				{
					double maxClose = -DBL_MAX, minClose = DBL_MAX;
					for(size_t i = 0; i < g->nodes.size(); i++)
					{
						bounds(g->nodes[i]->property["ClosenessCenter"].toDouble(), maxClose, minClose);
					}
					for(size_t i = 0; i < g->nodes.size(); i++)
					{
						// normalized
						Structure::Node * n = g->nodes[i];
						double d = (n->property["ClosenessCenter"].toDouble() - minClose) / (maxClose - minClose);
						n->property["ClosenessCenter"].setValue( d );
					}
				}

				// Store values
				closeness.push_back( std::vector< double >() );

				for(size_t i = 0; i < g->nodes.size(); i++)
				{
					closeness.back().push_back( g->nodes[i]->property["ClosenessCenter"].toDouble() );
				}
			}

			GraphDistance gd(g);
			std::vector<double> agd = gd.averageGeodesicDistance(RES);
			double agd_max = -DBL_MAX, agd_min = DBL_MAX;
			for(int i = 0; i < agd.size() - 1 ; i++) bounds(agd[i], agd_max, agd_min);
			
			agdVals.push_back( std::vector<double>() );

			for(size_t i = 0; i < skeletons.back().size(); i++)
			{			
				Vector3 V = mean( skeletons.back()[i] );

				double val = agd[ gd.closestPointIndex( V ) ];
				agdVals.back().push_back( (val - agd_min) / (agd_max-agd_min) );
			}*/

			// Debug AGD
			/*if( false ){
				starlab::PointSoup * ps = new starlab::PointSoup;
				for(int i = 0; i < agd.size(); i++){
					QColor c = starlab::qtJetColor((agd[i] - agd_min) / (agd_max - agd_min));
					ps->addPoint(gd.allPoints[i], c);
				}
				debug.push_back(ps);
				return Eigen::MatrixXd();
			}*/
		}

		/// Compute relative lengths w.r.t. bounding box
		{
			Eigen::AlignedBox3d bbox = g->bbox();
			lengths.push_back(std::vector< double >());

			for(size_t i = 0; i < skeletons.back().size(); i++)
			{
				lengths.back().push_back( length( skeletons.back()[i] ) );
			}
		}

		/// Compute relative position w.r.t. shape's bounding box
		{
			Eigen::AlignedBox3d bbox = g->bbox();

			Vector3 center = bbox.center();
			Vector3 diag = bbox.max() - center;

			double Height = bbox.sizes().z();
			double Length = bbox.sizes().y();
			double Width = bbox.sizes().x();

			Vec3d mP = bbox.min(); // this is the origin of our local frame
			Vec3d mS = Width * Vec3d(1,0,0);
			Vec3d mT = Length * Vec3d(0,1,0);
			Vec3d mU = Height * Vec3d(0,0,1);

			centerCoord.push_back(std::vector< Vector3 >());

			for(size_t i = 0; i < skeletons.back().size(); i++)
			{
				// Using center of the node
				Vector3 V = mean( skeletons.back()[i] );

				double s = 0, t = 0, u = 0;

				Vector3 TXU = cross(mT,mU);
				s = dot(TXU, Vector3(V - mP)) / dot(TXU, mS);

				Vector3 SXU = cross(mS,mU);
				t = dot(SXU, Vector3(V - mP)) / dot(SXU, mT);

				Vector3 SXT = cross(mS,mT);
				u = dot(SXT, Vector3(V - mP)) / dot(SXT, mU);

				centerCoord.back().push_back( Vector3(s,t,u) );
			}
		}
	}

	n = (int)lengths.front().size();
	m = (int)lengths.back().size();

	Eigen::MatrixXd g = Eigen::MatrixXd::Zero(n, m);

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m ; j++)
		{
			// Translation
			Vector3 meanSrc = centerCoord[0][i];
			Vector3 meanTgt = centerCoord[1][j];
			Vector3 t = meanSrc - meanTgt;
			double translation = t.norm();

			// Relative height
			double diffZ = abs(meanSrc.z() - meanTgt.z());

			// Rotation
			Vector3 srcDir = skeletons[0][i].back() - skeletons[0][i].front();
			Vector3 tgtDir = skeletons[1][j].back() - skeletons[1][j].front();
			double angle = 1.0 - abs( dot( srcDir.normalized(), tgtDir.normalized() ) );
			
			// Scaling
			double len = abs( lengths[0][i] - lengths[1][j] );

			// Structure
			//double structure = closeness[0][i] - closeness[1][j];
			double structure = abs(agdVals[0][i] - agdVals[1][j]);

			// Distortion
			double distortion = 0;

			// Energy:
			switch (selected_metric){
			case Corresponder::ALL:			g(i,j) = translation + diffZ + angle + len + structure + distortion; break;
			case Corresponder::POSITIONAL:	g(i,j) = translation;	break;
			case Corresponder::HEIGHT:		g(i,j) = diffZ;			break;
			case Corresponder::ROTATIONAL:	g(i,j) = angle;			break;
			case Corresponder::SCALING:		g(i,j) = len;			break;
			case Corresponder::STRUCTURE:	g(i,j) = structure;		break;
			case Corresponder::DISTORTION:	g(i,j) = distortion;	break;
			default: break;
			}
		}
	}

	// Colorize based on energy of source[idx]
	{
		int idx = qMin(selected_nidx, (int)g.rows() - 1);

		if(prop["forceBest"].toBool())
		{
			QMap<int,double> costs;
			for(auto i = 0; i < g.rows(); i++) costs[i] = g.row(i).minCoeff();
			QList< QPair<double, int> > sorted = sortQMapByValue(costs);

			idx = (*(sorted.begin() + idx)).second;
		}

		for(int i = 0; i < n; i++)
		{
			starlab::LineSegments * ls = new starlab::LineSegments(3);
			ls->addLines(skeletons[0][i], i == idx ? Qt::green : Qt::black );
			ls->translation += Vector3(source->property["posX"].toDouble(), 0, 0);
			debug.push_back(ls);
		}

		double minCoeff = g.row(idx).minCoeff();
		double maxCoeff = g.row(idx).maxCoeff();

		for(int j = 0; j < m ; j++)
		{
			QColor c = starlab::qtJetColor((g(idx,j) - minCoeff) / (maxCoeff - minCoeff));

			starlab::LineSegments * ls = new starlab::LineSegments(3);
			ls->addLines(skeletons[1][j], c);
			ls->translation += Vector3(target->property["posX"].toDouble(), 0, 0);
			debug.push_back(ls);
		}
	}

	if( false )
	{
		GraphExplorer * ge = new GraphExplorer();
		ge->update(target);
		ge->show();

		GraphExplorer * ge2 = new GraphExplorer();
		ge2->update(source);
		ge2->show();
	}
}
