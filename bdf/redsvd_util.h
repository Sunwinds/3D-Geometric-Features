/* 
*  Copyright (c) 2010 Daisuke Okanohara
* 
*   Redistribution and use in source and binary forms, with or without
*   modification, are permitted provided that the following conditions
*   are met:
* 
*   1. Redistributions of source code must retain the above Copyright
*      notice, this list of conditions and the following disclaimer.
*
*   2. Redistributions in binary form must reproduce the above Copyright
*      notice, this list of conditions and the following disclaimer in the
*      documentation and/or other materials provided with the distribution.
*
*   3. Neither the name of the authors nor the names of its contributors
*      may be used to endorse or promote products derived from this
*      software without specific prior written permission.
*/
#pragma once

#include <iostream>

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

namespace REDSVD {

	const float SVD_EPS = 0.0001f;

	typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrixXf;
	typedef std::vector<std::pair<int, float> > fv_t;

	class Util{
	public:
		static void convertFV2Mat(const std::vector<fv_t>& fvs, SMatrixXf& A);
		static void sampleGaussianMat(Eigen::MatrixXf& x);
		static void processGramSchmidt(Eigen::MatrixXf& mat);

	private:
		static void sampleTwoGaussian(float& f1, float& f2);
	};

	void Util::sampleTwoGaussian(float& f1, float& f2){
		float v1 = (float)(rand() + 1.f) / ((float)RAND_MAX+2.f);
		float v2 = (float)(rand() + 1.f) / ((float)RAND_MAX+2.f);
		float len = sqrt(-2.f * log(v1));
		f1 = len * cos(2.f * M_PI * v2);
		f2 = len * sin(2.f * M_PI * v2);
	}

	void Util::sampleGaussianMat(MatrixXf& mat){
		for (int i = 0; i < mat.rows(); ++i){
			int j = 0;
			for ( ; j+1 < mat.cols(); j += 2){
				float f1, f2;
				sampleTwoGaussian(f1, f2);
				mat(i,j  ) = f1;
				mat(i,j+1) = f2;
			}
			for (; j < mat.cols(); j ++){
				float f1, f2;
				sampleTwoGaussian(f1, f2);
				mat(i, j)  = f1;
			}
		}
	} 


	void Util::processGramSchmidt(MatrixXf& mat){
		for (int i = 0; i < mat.cols(); ++i){
			for (int j = 0; j < i; ++j){
				float r = mat.col(i).dot(mat.col(j));
				mat.col(i) -= r * mat.col(j);
			}
			float norm = mat.col(i).norm();
			if (norm < SVD_EPS){
				for (int k = i; k < mat.cols(); ++k){
					mat.col(k).setZero();
				} 
				return;
			}
			mat.col(i) *= (1.f / norm);
		}
	}

	void Util::convertFV2Mat(const vector<fv_t>& fvs, REDSVD::SMatrixXf& A){
		size_t maxID = 0;
		size_t nonZeroNum = 0;
		for (size_t i = 0; i < fvs.size(); ++i){
			const fv_t& fv(fvs[i]);
			for (size_t j = 0; j < fv.size(); ++j){
				maxID = max(size_t(fv[j].first+1), maxID);
			}
			nonZeroNum += fv.size();
		}
		A.resize((int)fvs.size(), (int)maxID);
		A.reserve((int)nonZeroNum);
		for (size_t i = 0; i < fvs.size(); ++i){
			A.startVec((int)i);
			const fv_t& fv(fvs[i]);
			for (size_t j = 0; j < fv.size(); ++j){
				A.insertBack((int)i, fv[j].first) = fv[j].second;
			}
		}
		A.finalize();
	}
}
