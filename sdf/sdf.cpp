#include <QElapsedTimer>
#include <omp.h>

#include "sdf.h"
#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#include <Eigen/Geometry>

#define DEG2RAD 0.01745329

void sdf::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichInt("numCones", 4, "Num. cones"));
	pars->addParam(new RichInt("raysInCone", 8, "Rays in cone"));
	pars->addParam(new RichFloat("coneSeparation", 20.0f,"Cone separation"));
	pars->addParam(new RichBool("gaussianWeights", true, "Gaussian weights"));
	pars->addParam(new RichInt("Smooth", 1, "Smooth"));

	pars->addParam(new RichBool("Visualize", true, "Visualize"));
	pars->addParam(new RichBool("VisualizeRays", true, "Visualize rays"));
}

static inline Vector3 orthogonalVector(const Vector3& n) {
	if ((abs(n.y()) >= 0.9 * abs(n.x())) &&
		abs(n.z()) >= 0.9 * abs(n.x())) return Vector3(0.0, -n.z(), n.y());
	else if ( abs(n.x()) >= 0.9 * abs(n.y()) &&
		abs(n.z()) >= 0.9 * abs(n.y()) ) return Vector3(-n.z(), 0.0, n.x());
	else return Vector3(-n.y(), n.x(), 0.0);
}

Ray get_normal_opposite(const Vector3 & p, const Vector3 & n)
{
	Vector3 direction = n * -1;	
	return Ray(p + direction * 1e-6, direction);
}

std::list<Ray> get_cone_rays(const Vector3& p, const Vector3& n, const double degree, const int spread){
	std::list<Ray> list;

	Ray ray = get_normal_opposite(p, n);

	double theta = 2.0 * M_PI / spread;

	Eigen::Quaterniond q( Eigen::AngleAxisd(degree, orthogonalVector(n) ) );
	Eigen::Quaterniond qs( Eigen::AngleAxisd(theta, n ) );

	Vector3 spinMe = q * ray.direction.normalized();

	for (int i = 0; i < spread; i++){
		list.push_back(Ray(p + spinMe * 1e-6, spinMe));
		spinMe = qs * spinMe;
	}

	return list;
}

double sdf::compute_volume(const Vector3 & p, const Vector3 & n, int numCones, 
						   double coneSeparation, int raysInCone, bool gaussianWeights, 
						   Octree * octree, QList<Ray> * allRays)
{
	double result = 0.0;

	double distanceCounter = 0.0;
	double intersectionCounter = 0.0;
	double distance = DBL_MAX;

	// basically shoot rays opposite the normal and do something with the distances we receive
	Ray ray = get_normal_opposite(p, n);

	std::vector<double> values;
	std::vector<double> weights;
	values.reserve(numCones * raysInCone + 1);

	// first check straight opposite to normal
	if( octree->isIntersectsWithRay(ray, &distance) ) 
	{
		values.push_back(distance);
		weights.push_back(1.0);

		if( allRays ){
			ray.direction *= distance;
			allRays->push_back(ray);
		}
	}
	
	const double gaussianVar = 120.0 * DEG2RAD;	
	double gaussianDiv2 = 1.0 / (2.0 * pow(gaussianVar,2));	

	for (int i = 0; i < numCones; i++) 
	{
		// now check cone
		double degree = coneSeparation * DEG2RAD * (i+1);
		std::list<Ray> list = get_cone_rays(p, n, degree, raysInCone);

		for (std::list<Ray>::iterator it = list.begin(); it != list.end(); it++) 
		{
			distance = DBL_MAX;
			ray = *it;

			if ( octree->isIntersectsWithRay(ray, &distance) ) 
			{
				double weight = 1.0;

				if (gaussianWeights) weight = /*gaussianDiv1 * */expf(-pow(degree, 2) * gaussianDiv2);

				values.push_back(distance);
				weights.push_back(weight);

				if( allRays ){
					ray.direction *= distance;
					allRays->push_back(ray);
				}
			}			
		}
	}

	// calculate the average	
	for (int i = 0; i < values.size(); i++) {
		distanceCounter += values[i] * weights[i];
		intersectionCounter += weights[i];
	}
	
	result = distanceCounter / intersectionCounter;

	return result;
}

void sdf::applyFilter(RichParameterSet *pars)
{
	drawArea()->clear();

	int N = mesh()->n_vertices();
	Vector3VertexProperty points = mesh()->vertex_coordinates();
	Vector3VertexProperty normals = mesh()->vertex_normals();

	ScalarVertexProperty m_sdf = mesh()->vertex_property<Scalar>("v:sdf");

	// Parameters
	int m_numCones = pars->getInt("numCones");
	int m_raysInCone = pars->getInt("raysInCone");
	double m_coneSeperation = pars->getFloat("coneSeparation");
	bool m_gaussianWeights = pars->getBool("gaussianWeights");

	Octree octree( mesh() );

	// Visualization
	QList<Ray> allRays;
	starlab::VectorSoup * vs = new starlab::VectorSoup();

	#pragma omp parallel for
	for(int i = 0; i < N; i++)
	{
		//QList<Ray> curRays;

		Vertex v(i);
        qDebug() << i;
		m_sdf[v] = sdf::compute_volume( points[v], normals[v], 
			m_numCones, m_coneSeperation, m_raysInCone, m_gaussianWeights, &octree /*, curRays*/ );

		//allRays += curRays;
		//foreach(Ray r, allRays) vs->addVector(r.origin, r.direction);
	}

	// Smooth the function
	SurfaceMeshHelper(mesh()).smoothVertexProperty<Scalar>("v:sdf", pars->getInt("Smooth"));

	double maxVal = -DBL_MAX, minVal = DBL_MAX;
    foreach(Vertex v, mesh()->vertices()) { maxVal = qMax(m_sdf[v], maxVal); minVal = qMin(m_sdf[v], minVal);}
    foreach(Vertex v, mesh()->vertices()) m_sdf[v] = (m_sdf[v] - minVal) / (maxVal - minVal);

    foreach(Vertex v, mesh()->vertices()) drawArea()->drawPoint( points[v], 8, starlab::qtJetColor(m_sdf[v]));

	drawArea()->addRenderObject(vs);
}
