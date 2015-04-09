#include "SurfaceMeshModel.h"
#include "principal_curvature.h"

using namespace igl;


inline CurvatureCalculator::CurvatureCalculator()
{
	this->localMode=false;
	this->projectionPlaneCheck=false;
	this->sphereRadius=5;
	this->st=K_RING_SEARCH;
	this->nt=PROJ_PLANE;
	this->montecarlo=false;
	this->montecarloN=0;
	this->kRing=2;
	this->svd=true;
	this->zeroDetCheck=true;
	this->curvatureComputed=false;
	this->expStep=true;
}

inline void CurvatureCalculator::fitQuadric (Eigen::Vector3d v, std::vector<Eigen::Vector3d> ref, const std::vector<int>& vv, Quadric *q)
{
	std::vector<Eigen::Vector3d> points;
	points.reserve (vv.size());

	for (unsigned int i = 0; i < vv.size(); ++i) {

		Eigen::Vector3d  cp = vertices.row(vv[i]);

		// vtang non e` il v tangente!!!
		Eigen::Vector3d  vTang = cp - v;

		double x = vTang.dot(ref[0]);
		double y = vTang.dot(ref[1]);
		double z = vTang.dot(ref[2]);
		points.push_back(Eigen::Vector3d (x,y,z));
	}
	*q = Quadric::fit (points, zeroDetCheck, svd);
}

inline void CurvatureCalculator::finalEigenStuff (size_t i, std::vector<Eigen::Vector3d> ref, Quadric q)
{

	double a = q.a();
	double b = q.b();
	double c = q.c();
	double d = q.d();
	double e = q.e();

	//  if (fabs(a) < 10e-8 || fabs(b) < 10e-8)
	//  {
	//    std::cout << "Degenerate quadric: " << i << std::endl;
	//  }

	double E = 1.0 + d*d;
	double F = d*e;
	double G = 1.0 + e*e;

	Eigen::Vector3d n = Eigen::Vector3d(-d,-e,1.0).normalized();

	double L = 2.0 * a * n[2];
	double M = b * n[2];
	double N = 2 * c * n[2];


	// ----------------- Eigen stuff
	Eigen::Matrix2d m;
	m << L*G - M*F, M*E-L*F, M*E-L*F, N*E-M*F;
	m = m / (E*G-F*F);
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(m);

	Eigen::Vector2d c_val = eig.eigenvalues();
	Eigen::Matrix2d c_vec = eig.eigenvectors();

	// std::cerr << "c_val:" << c_val << std::endl;
	// std::cerr << "c_vec:" << c_vec << std::endl;

	// std::cerr << "c_vec:" << c_vec(0) << " "  << c_vec(1) << std::endl;

	c_val = -c_val;

	Eigen::Vector3d v1, v2;
	v1[0] = c_vec(0);
	v1[1] = c_vec(1);
	v1[2] = 0; //d * v1[0] + e * v1[1];

	v2[0] = c_vec(2);
	v2[1] = c_vec(3);
	v2[2] = 0; //d * v2[0] + e * v2[1];


	// v1 = v1.normalized();
	// v2 = v2.normalized();

	Eigen::Vector3d v1global = ref[0] * v1[0] + ref[1] * v1[1] + ref[2] * v1[2];
	Eigen::Vector3d v2global = ref[0] * v2[0] + ref[1] * v2[1] + ref[2] * v2[2];

	v1global.normalize();
	v2global.normalize();

	v1global *= c_val(0);
	v2global *= c_val(1);

	if (c_val[0] > c_val[1])
	{
		curv[i]=std::vector<double>(2);
		curv[i][0]=c_val(1);
		curv[i][1]=c_val(0);
		curvDir[i]=std::vector<Eigen::Vector3d>(2);
		curvDir[i][0]=v2global;
		curvDir[i][1]=v1global;
	}
	else
	{
		curv[i]=std::vector<double>(2);
		curv[i][0]=c_val(0);
		curv[i][1]=c_val(1);
		curvDir[i]=std::vector<Eigen::Vector3d>(2);
		curvDir[i][0]=v1global;
		curvDir[i][1]=v2global;
	}
	// ---- end Eigen stuff
}

inline void CurvatureCalculator::getKRing(const size_t start, const double r, std::vector<int>&vv)
{
	size_t bufsize = vertices.rows();
	vv.reserve(bufsize);
	std::list<std::pair<int,int> > queue;
	bool* visited = (bool*)calloc(bufsize,sizeof(bool));
	queue.push_back(std::pair<int,int>((int)start,0));
	visited[start]=true;
	while (!queue.empty())
	{
		int toVisit=queue.front().first;
		int distance=queue.front().second;
		queue.pop_front();
		vv.push_back(toVisit);
		if (distance<(int)r)
		{
			for (unsigned int i=0; i<vertex_to_vertices[toVisit].size(); i++)
			{
				int neighbor=vertex_to_vertices[toVisit][i];
				if (!visited[neighbor])
				{
					queue.push_back(std::pair<int,int> (neighbor,distance+1));
					visited[neighbor]=true;
				}
			}
		}
	}
	free(visited);
	return;
}

inline void CurvatureCalculator::getSphere(const size_t start, const double r, std::vector<int> &vv, int min)
{
	size_t bufsize=vertices.rows();
	vv.reserve(bufsize);
	std::list<size_t>* queue= new std::list<size_t>();
	bool* visited = (bool*)calloc(bufsize,sizeof(bool));
	queue->push_back(start);
	visited[start]=true;
	Eigen::Vector3d me=vertices.row(start);
	std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comparer >* extra_candidates= new  std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >, comparer >();
	while (!queue->empty())
	{
		size_t toVisit=queue->front();
		queue->pop_front();
		vv.push_back((int)toVisit);
		for (unsigned int i=0; i<vertex_to_vertices[toVisit].size(); i++)
		{
			int neighbor=vertex_to_vertices[toVisit][i];
			if (!visited[neighbor])
			{
				Eigen::Vector3d neigh=vertices.row(neighbor);
				double distance=(me-neigh).norm();
				if (distance<r)
					queue->push_back(neighbor);
				else if ((int)vv.size()<min)
					extra_candidates->push(std::pair<int,double>(neighbor,distance));
				visited[neighbor]=true;
			}
		}
	}
	while (!extra_candidates->empty() && (int)vv.size()<min)
	{
		std::pair<int, double> cand=extra_candidates->top();
		extra_candidates->pop();
		vv.push_back(cand.first);
		for (unsigned int i=0; i<vertex_to_vertices[cand.first].size(); i++)
		{
			int neighbor=vertex_to_vertices[cand.first][i];
			if (!visited[neighbor])
			{
				Eigen::Vector3d neigh=vertices.row(neighbor);
				double distance=(me-neigh).norm();
				extra_candidates->push(std::pair<int,double>(neighbor,distance));
				visited[neighbor]=true;
			}
		}
	}
	free(extra_candidates);
	free(queue);
	free(visited);
}

inline Eigen::Vector3d CurvatureCalculator::project(Eigen::Vector3d v, Eigen::Vector3d  vp, Eigen::Vector3d ppn)
{
	return (vp - (ppn * ((vp - v).dot(ppn))));
}

inline void CurvatureCalculator::computeReferenceFrame(size_t i, Eigen::Vector3d normal, std::vector<Eigen::Vector3d>& ref )
{

	Eigen::Vector3d longest_v=Eigen::Vector3d::Zero();
	longest_v=Eigen::Vector3d(vertices.row(vertex_to_vertices[i][0]));

	longest_v=(project(vertices.row(i),longest_v,normal)-Eigen::Vector3d(vertices.row(i))).normalized();

	/* L'ultimo asse si ottiene come prodotto vettoriale tra i due
	* calcolati */
	Eigen::Vector3d y_axis=(normal.cross(longest_v)).normalized();
	ref[0]=longest_v;
	ref[1]=y_axis;
	ref[2]=normal;
}

inline void CurvatureCalculator::getAverageNormal(size_t j, std::vector<int> vv, Eigen::Vector3d& normal)
{
	normal=(vertex_normals.row(j)).normalized();
	if (localMode)
		return;

	for (unsigned int i=0; i<vv.size(); i++)
	{
		normal+=vertex_normals.row(vv[i]).normalized();
	}
	normal.normalize();
}

inline void CurvatureCalculator::getProjPlane(size_t j, std::vector<int> vv, Eigen::Vector3d& ppn)
{
	int nr;
	double a, b, c;
	double nx, ny, nz;
	double abcq;

	a = b = c = 0;

	if (localMode)
	{
		for (unsigned int i=0; i<vertex_to_faces.at(j).size(); ++i)
		{
			Eigen::Vector3d faceNormal=face_normals.row(vertex_to_faces.at(j).at(i));
			a += faceNormal[0];
			b += faceNormal[1];
			c += faceNormal[2];
		}
	}
	else
	{
		for (unsigned int i=0; i<vv.size(); ++i)
		{
			a+= vertex_normals.row(vv[i])[0];
			b+= vertex_normals.row(vv[i])[1];
			c+= vertex_normals.row(vv[i])[2];
		}
	}
	nr = rotateForward (&a, &b, &c);
	abcq = a*a + b*b + c*c;
	nx = sqrt (a*a / abcq);
	ny = sqrt (b*b / abcq);
	nz = sqrt (1 - nx*nx - ny*ny);
	rotateBackward (nr, &a, &b, &c);
	rotateBackward (nr, &nx, &ny, &nz);

	ppn = chooseMax (Eigen::Vector3d(nx, ny, nz), Eigen::Vector3d (a, b, c), a * b);
	ppn.normalize();
}


inline double CurvatureCalculator::getAverageEdge()
{
	double sum = 0;
	int count = 0;

	for (int i = 0; i<faces.rows(); i++)
	{
		for (short unsigned j=0; j<3; j++)
		{
			Eigen::Vector3d p1=vertices.row(faces.row(i)[j]);
			Eigen::Vector3d p2=vertices.row(faces.row(i)[(j+1)%3]);

			double l = (p1-p2).norm();

			sum+=l;
			++count;
		}
	}

	return (sum/(double)count);
}


inline void CurvatureCalculator::applyProjOnPlane(Eigen::Vector3d ppn, std::vector<int> vin, std::vector<int> &vout)
{
	for (std::vector<int>::iterator vpi = vin.begin(); vpi != vin.end(); ++vpi)
		if (vertex_normals.row(*vpi) * ppn > 0.0f)
			vout.push_back (*vpi);
}

inline void CurvatureCalculator::applyMontecarlo(std::vector<int>& vin, std::vector<int> *vout)
{
	if (montecarloN >= vin.size ())
	{
		*vout = vin;
		return;
	}

	float p = ((float) montecarloN) / (float) vin.size();
	for (std::vector<int>::iterator vpi = vin.begin(); vpi != vin.end(); ++vpi)
	{
		float r;
		if ((r = ((float)rand () / RAND_MAX)) < p)
		{
			vout->push_back (*vpi);
		}
	}
}

inline void CurvatureCalculator::computeCurvature()
{
	using namespace std;

	//CHECK che esista la mesh
	size_t vertices_count=vertices.rows() ;

	if (vertices_count <=0)
		return;

	curvDir=std::vector< std::vector<Eigen::Vector3d> >(vertices_count);
	curv=std::vector<std::vector<double> >(vertices_count);


	scaledRadius = getAverageEdge() * sphereRadius;
	kRing = kRing * sphereRadius;

	std::vector<int> vv;
	std::vector<int> vvtmp;
	Eigen::Vector3d normal;

	//double time_spent;
	//double searchtime=0, ref_time=0, fit_time=0, final_time=0;

	for (size_t i=0; i<vertices_count; ++i)
	{
		vv.clear();
		vvtmp.clear();
		Eigen::Vector3d me=vertices.row(i);
		switch (st)
		{
		case SPHERE_SEARCH:
			getSphere(i, scaledRadius, vv, 6);
			break;
		case K_RING_SEARCH:
			getKRing(i,kRing,vv);
			break;
		default:
			fprintf(stderr,"Error: search type not recognized");
			return;
		}

		std::vector<Eigen::Vector3d> ref(3);
		if (vv.size()<6)
		{
			std::cerr << "Could not compute curvature of radius " << scaledRadius << endl;
			return;
		}
		switch (nt)
		{
		case AVERAGE:
			getAverageNormal(i,vv,normal);
			break;
		case PROJ_PLANE:
			getProjPlane(i,vv,normal);
			break;
		default:
			fprintf(stderr,"Error: normal type not recognized");
			return;
		}

		if (projectionPlaneCheck)
		{
			vvtmp.reserve (vv.size ());
			applyProjOnPlane (normal, vv, vvtmp);
			if (vvtmp.size() >= 6)
				vv = vvtmp;
		}

		if (vv.size()<6)
		{
			std::cerr << "Could not compute curvature of radius " << scaledRadius << endl;
			return;
		}
		if (montecarlo)
		{
			if(montecarloN<6)
				break;
			vvtmp.reserve(vv.size());
			applyMontecarlo(vv,&vvtmp);
			vv=vvtmp;
		}

		if (vv.size()<6)
			return;
		computeReferenceFrame(i,normal,ref);

		Quadric q;
		fitQuadric (me, ref, vv, &q);
		finalEigenStuff(i,ref,q);
	}

	lastRadius=sphereRadius;
	curvatureComputed=true;
}

inline void CurvatureCalculator::printCurvature(std::string outpath)
{
	using namespace std;
	if (!curvatureComputed)
		return;

	std::ofstream of;
	of.open(outpath.c_str());

	if (!of)
	{
		fprintf(stderr, "Error: could not open output file %s\n", outpath.c_str());
		return;
	}

	size_t vertices_count=vertices.rows();
	of << vertices_count << endl;
	for (int i=0; i<vertices_count; i++)
	{
		of << curv[i][0] << " " << curv[i][1] << " " << curvDir[i][0][0] << " " << curvDir[i][0][1] << " " << curvDir[i][0][2] << " " <<
			curvDir[i][1][0] << " " << curvDir[i][1][1] << " " << curvDir[i][1][2] << endl;
	}

	of.close();

}

void igl::principal_curvature(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	const Eigen::MatrixXd& NV,
	const Eigen::MatrixXd& NF,
	Eigen::MatrixXd& K1,
	Eigen::MatrixXd& K2,
	Eigen::MatrixXd& PD1,
	Eigen::MatrixXd& PD2,
	unsigned radius, double * maxValue)
{
	using namespace std;

	// Preallocate memory
	K1.resize(V.rows(),1);
	K2.resize(V.rows(),1);
	PD1.resize(V.rows(),3);
	PD2.resize(V.rows(),3);

	CurvatureCalculator cc;
	cc.vertices = V;
	cc.faces = F;
	cc.vertex_normals = NV;
	cc.face_normals = NF;
	cc.sphereRadius = radius;

	// Adjacency
	adjacency_list(F, cc.vertex_to_vertices);
	vf(V, F, cc.vertex_to_faces, cc.vertex_to_faces_index);

	// Compute
	cc.computeCurvature();

	if(maxValue) *maxValue = 0.0;

	// Copy it back
	for (unsigned i=0; i<V.rows(); i++)
	{
		Eigen::Vector3d d1;
		Eigen::Vector3d d2;
		d1 << cc.curvDir[i][0][0], cc.curvDir[i][0][1], cc.curvDir[i][0][2];
		d2 << cc.curvDir[i][1][0], cc.curvDir[i][1][1], cc.curvDir[i][1][2];
		d1.normalize();
		d2.normalize();
		d1 *= cc.curv[i][0];
		d2 *= cc.curv[i][1];

		K1(i,0) = cc.curv[i][0];
		K2(i,0) = cc.curv[i][1];
		PD1.row(i) = d1;
		PD2.row(i) = d2;

		// Keep track of maximum
		if(maxValue)
		{
			if(cc.curv[i][0] > *maxValue) *maxValue = cc.curv[i][0];
			if(cc.curv[i][1] > *maxValue) *maxValue = cc.curv[i][1];
		}

		if (PD1.row(i) * PD2.row(i).transpose() > 10e-6)
		{
			cerr << "Something is wrong, vertex: i" << endl;
			PD1.row(i) *= 0;
			PD2.row(i) *= 0;
		}
	}
}
