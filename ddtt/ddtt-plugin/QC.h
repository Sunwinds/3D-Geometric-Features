#pragma once

#include <vector>
#include <Eigen/Core>

struct isEqualZero{
	std::vector<int>* GT;
	isEqualZero(std::vector<int> *g){GT = g;}
	void operator()(double i){static int it = 0; if(i == 0) GT->push_back(it); it++;}
};

/* The Quadratic-Chi Histogram Distance Family. Ofir Pele, Michael Werman. [ECCV 2010] */
struct QC
{
	static double distance( std::vector<double> histogramP, std::vector<double> histogramQ )
	{
		// Defaults
		int THRESHOLD = 3;
		double m = 0.9;
		int N = (int)histogramP.size();

		// Convert to Eigen vectors
		Eigen::Map<Eigen::VectorXd> P( &histogramP[0], N );
		Eigen::Map<Eigen::VectorXd> Q( &histogramQ[0], N );

		// The sparse bin-similarity matrix

		// PLACMENT / NOT A GOOD SIMILAIRTY MATRIX
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);
		for(int i = 0; i < N; i++)
			for(int j = std::max(0, i-THRESHOLD); j < std::min(N, i+THRESHOLD); j++)
				A(i,j) = 1.0 - (double(abs(i-j)) / THRESHOLD);
			
		Eigen::VectorXd Z = (P+Q).transpose() * A;
		
		// Z(Z==0) = 1
		// 1.0 can be any number as Z_i==0 iff D_i=0
		std::vector<int> idx;
		std::for_each(Z.data(), Z.data() + Z.rows()*Z.cols(), isEqualZero(&idx));
		for(int i=0;i<idx.size();++i) Z(idx[i]) = 1.0;

		// Z = Z.^m
		for(int i = 0; i < Z.size(); i++) Z(i) = pow(Z(i), m);

		// D = (P-Q)./Z
		Eigen::VectorXd D(N);
		for(int i = 0; i < D.size(); i++) D(i) = (P(i) - Q(i)) / Z(i);

		double DADt = (D.transpose() * A * D).maxCoeff(); // its a single number..
		return std::sqrt( std::max(0.0, DADt) );
	}

	static double L2distance( std::vector<double> histogramP, std::vector<double> histogramQ )
	{
		int N = (int)histogramP.size();
		if(N == 0 || histogramP.size() != histogramQ.size()) return 1e12;

		// Convert to Eigen vectors
		Eigen::Map<Eigen::VectorXd> P( &histogramP[0], N );
		Eigen::Map<Eigen::VectorXd> Q( &histogramQ[0], N );
		 
		return ( P - Q ).norm();
	}
};
