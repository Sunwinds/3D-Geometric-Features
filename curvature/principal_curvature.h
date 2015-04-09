// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include <Eigen/Geometry>
#include <Eigen/Dense>

#include <vector>
#include <stdio.h>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <list>
#include <cmath>
#include <limits>

#include <Eigen/SparseCholesky>

// Compute the principal curvature directions and magnitude of the given triangle mesh
//   DerivedV derived from vertex positions matrix type: i.e. MatrixXd
//   DerivedF derived from face indices matrix type: i.e. MatrixXi
// Inputs:
//   V       eigen matrix #V by 3
//   F       #F by 3 list of mesh faces (must be triangles)
//   radius  controls the size of the neighbourhood used, 1 = average edge lenght
//
// Outputs:
//   PD1 #V by 3 maximal curvature direction for each vertex.
//   PD2 #V by 3 minimal curvature direction for each vertex.  
//  Note that to maximal/minimal curvature value is the rowise norm of PD1/PD2
//
// See also: moveVF, moveFV
//
// This function has been developed by: Nikolas De Giorgis, Luigi Rocca and Enrico Puppo.
// The algorithm is based on:
// Efficient Multi-scale Curvature and Crease Estimation
// Daniele Panozzo, Enrico Puppo, Luigi Rocca
// GraVisMa, 2010
namespace igl{

typedef enum
{
	SPHERE_SEARCH,
	K_RING_SEARCH
} searchType;

typedef enum
{
	AVERAGE,
	PROJ_PLANE
} normalType;

class comparer
{
public:
	inline bool operator() (const std::pair<int, double>& lhs, const std::pair<int, double>&rhs) const
	{
		return lhs.second>rhs.second;
	}
};

class CurvatureCalculator
{
public:
	/* Row number i represents the i-th vertex, whose columns are:
	curv[i][0] : K2
	curv[i][1] : K1
	curvDir[i][0] : PD1
	curvDir[i][1] : PD2
	*/
	std::vector< std::vector<double> > curv;
	std::vector< std::vector<Eigen::Vector3d> > curvDir;
	bool curvatureComputed;
	class Quadric
	{
	public:

		inline Quadric ()
		{
			a() = b() = c() = d() = e() = 1.0;
		}

		inline Quadric(double av, double bv, double cv, double dv, double ev)
		{
			a() = av;
			b() = bv;
			c() = cv;
			d() = dv;
			e() = ev;
		}

		inline double& a() { return data[0];}
		inline double& b() { return data[1];}
		inline double& c() { return data[2];}
		inline double& d() { return data[3];}
		inline double& e() { return data[4];}

		double data[5];

		inline double evaluate(double u, double v)
		{
			return a()*u*u + b()*u*v + c()*v*v + d()*u + e()*v;
		}

		inline double du(double u, double v)
		{
			return 2.0*a()*u + b()*v + d();
		}

		inline double dv(double u, double v)
		{
			return 2.0*c()*v + b()*u + e();
		}

		inline double duv(double , double )
		{
			return b();
		}

		inline double duu(double , double )
		{
			return 2.0*a();
		}

		inline double dvv(double , double )
		{
			return 2.0*c();
		}


		inline static Quadric fit(std::vector<Eigen::Vector3d> &VV, bool , bool)
		{
			using namespace std;
			assert(VV.size() >= 5);
			if (VV.size() < 5)
			{
				cerr << "ASSERT FAILED!" << endl;
				exit(0);
			}

			Eigen::MatrixXd A(VV.size(),5);
			Eigen::MatrixXd b(VV.size(),1);
			Eigen::MatrixXd sol(5,1);

			for(unsigned int c=0; c < VV.size(); ++c)
			{
				double u = VV[c][0];
				double v = VV[c][1];
				double n = VV[c][2];

				A(c,0) = u*u;
				A(c,1) = u*v;
				A(c,2) = v*v;
				A(c,3) = u;
				A(c,4) = v;

				b(c) = n;
			}

			sol=A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

			return Quadric(sol(0),sol(1),sol(2),sol(3),sol(4));
		}
	};

public:

	Eigen::MatrixXd vertices;
	// Face list of current mesh    (#F x 3) or (#F x 4)
	// The i-th row contains the indices of the vertices that forms the i-th face in ccw order
	Eigen::MatrixXi faces;

	std::vector<std::vector<int> > vertex_to_vertices;
	std::vector<std::vector<int> > vertex_to_faces;
	std::vector<std::vector<int> > vertex_to_faces_index;
	Eigen::MatrixXd face_normals;
	Eigen::MatrixXd vertex_normals;

	/* Size of the neighborhood */
	double sphereRadius;
	int kRing;

	bool localMode; /* Use local mode */
	bool projectionPlaneCheck; /* Check collected vertices on tangent plane */
	bool montecarlo;
	bool svd; /* Use svd calculation instead of pseudoinverse */
	bool zeroDetCheck; /* Check if the determinant is close to zero */
	unsigned int montecarloN;

	searchType st; /* Use either a sphere search or a k-ring search */
	normalType nt;

	double lastRadius;
	double scaledRadius;
	std::string lastMeshName;

	/* Benchmark related variables */
	bool expStep; /* True if we want the radius to increase exponentially */
	int step;  /* If expStep==false, by how much rhe radius increases on every step */
	int maxSize; /* The maximum limit of the radius in the benchmark */

	inline CurvatureCalculator();

	inline void finalEigenStuff (size_t, std::vector<Eigen::Vector3d>, Quadric );
	inline void fitQuadric (Eigen::Vector3d, std::vector<Eigen::Vector3d> ref, const  std::vector<int>& , Quadric *);
	inline void applyProjOnPlane(Eigen::Vector3d, std::vector<int>, std::vector<int>&);
	inline void getSphere(const size_t, const double, std::vector<int>&, int min);
	static inline void getSphere(const Eigen::MatrixXd & vertices, const std::vector<std::vector<int> > & vertex_to_vertices, const size_t start, const double r, std::vector<int> &vv, int min)
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

	inline void getKRing(const size_t, const double,std::vector<int>&);
	inline Eigen::Vector3d project(Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d);
	inline void computeReferenceFrame(size_t, Eigen::Vector3d, std::vector<Eigen::Vector3d>&);
	inline void getAverageNormal(size_t, std::vector<int>, Eigen::Vector3d&);
	inline void getProjPlane(size_t, std::vector<int>, Eigen::Vector3d&);
	inline void applyMontecarlo(std::vector<int>&,std::vector<int>*);
	inline void computeCurvature();
	inline void printCurvature(std::string outpath);
	inline double getAverageEdge();

	inline static int rotateForward (double *v0, double *v1, double *v2)
	{
		double t;

		if (abs(*v2) >= abs(*v1) && abs(*v2) >= abs(*v0))
			return 0;

		t = *v0;
		*v0 = *v2;
		*v2 = *v1;
		*v1 = t;

		return 1 + rotateForward (v0, v1, v2);
	}

	inline static void rotateBackward (int nr, double *v0, double *v1, double *v2)
	{
		double t;

		if (nr == 0)
			return;

		t = *v2;
		*v2 = *v0;
		*v0 = *v1;
		*v1 = t;

		rotateBackward (nr - 1, v0, v1, v2);
	}

	inline static Eigen::Vector3d chooseMax (Eigen::Vector3d n, Eigen::Vector3d abc, double ab)
	{
		int i, max_i;
		double max_sp;
		Eigen::Vector3d nt[8];

		n.normalize ();
		abc.normalize ();

		max_sp = - std::numeric_limits<double>::max();

		for (i = 0; i < 4; i++)
		{
			nt[i] = n;
			if (ab > 0)
			{
				switch (i)
				{
				case 0:
					break;

				case 1:
					nt[i][2] = -n[2];
					break;

				case 2:
					nt[i][0] = -n[0];
					nt[i][1] = -n[1];
					break;

				case 3:
					nt[i][0] = -n[0];
					nt[i][1] = -n[1];
					nt[i][2] = -n[2];
					break;
				}
			}
			else
			{
				switch (i)
				{
				case 0:
					nt[i][0] = -n[0];
					break;

				case 1:
					nt[i][1] = -n[1];
					break;

				case 2:
					nt[i][0] = -n[0];
					nt[i][2] = -n[2];
					break;

				case 3:
					nt[i][1] = -n[1];
					nt[i][2] = -n[2];
					break;
				}
			}

			if (nt[i].dot(abc) > max_sp)
			{
				max_sp = nt[i].dot(abc);
				max_i = i;
			}
		}

		return nt[max_i];
	}
};


template <typename FaceArray, typename Index>
inline void adjacency_list(const FaceArray & F,std::vector<std::vector<Index> >& A)
{
	A.clear(); 
	A.resize(F.maxCoeff()+1);

	// Loop over faces
	for(int i = 0;i<F.rows();i++)
	{
		// Loop over this face (triangle)
		for(int j = 0;j<3;j++)
		{
			// Get indices of edge: s --> d
			int s = F(i,j);
			int d = F(i,(j+1) % 3);
			A.at(s).push_back(d);
			A.at(d).push_back(s);
		}
	}

	// Remove duplicates
	for(int i=0; i<(int)A.size();++i)
	{
		std::sort(A[i].begin(), A[i].end());
		A[i].erase(std::unique(A[i].begin(), A[i].end()), A[i].end());
	}
}

template <typename DerivedV, typename DerivedF, typename IndexType>
inline void vf( const Eigen::PlainObjectBase<DerivedV>& V,const Eigen::PlainObjectBase<DerivedF>& F,
			   std::vector<std::vector<IndexType> >& VF,std::vector<std::vector<IndexType> >& VFi)
{
	VF.clear(); VF.resize(V.rows());
	VFi.clear(); VFi.resize(V.rows());

	for(int fi=0; fi<F.rows(); ++fi){
		for(int i = 0; i < F.cols(); ++i){
			VF[F(fi,i)].push_back(fi);
			VFi[F(fi,i)].push_back(i);
		}
	}
}

void principal_curvature( const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& NV, 
						 const Eigen::MatrixXd& NF, Eigen::MatrixXd& K1, Eigen::MatrixXd& K2, Eigen::MatrixXd& PD1, 
						 Eigen::MatrixXd& PD2, unsigned radius, double * maxValue);
}
