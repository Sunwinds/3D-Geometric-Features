/****************************************************************************
Compute conformal factors for mesh

Reference:
Ben-Chen, M., & Gostman, C.(2008).Characterizing shape using
conformal factors. In Proceedings of Eurographics workshop on
shape retrieval 2008.

URL:
http://www.cs.technion.ac.il/~gotsman/AmendedPubl/Miri/Shape08.pdf


Note:
According to author's paper, conformal factor is a good feature of shape for
object retrieval. It has several advantage properties, such as pose invariance
and un-sensitivity to noise, but it is sensitive to topology.

Lubin Fan
Nov. 24th, 2011
/***************************************************************************/

#pragma once

#include "SurfaceMeshPlugins.h"
#include "SurfaceMeshHelper.h"

#include "Eigen/LU"
#include "Eigen/Geometry"
#include "Eigen/Eigenvalues"
#include <Eigen/Sparse>

#include <QVector>

#include <vector>
#include <map>
#include <cmath>
#include <cstdio>

#define Pi 3.14159265
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

template<typename _Matrix_Type_>
bool pinv(const _Matrix_Type_ &a, _Matrix_Type_ &result, double	epsilon = std::numeric_limits<typename _Matrix_Type_::Scalar>::epsilon())
{
    if(a.rows() < a.cols())	return false;
    Eigen::JacobiSVD< _Matrix_Type_ > svd = a.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    typename _Matrix_Type_::Scalar tolerance = epsilon * std::max(a.cols(),	a.rows()) * svd.singularValues().array().abs().maxCoeff();
    result = svd.matrixV() * ( (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0)).matrix().asDiagonal() * svd.matrixU().adjoint();
    return true;
}


class CConformalFactor
{
//protected:
//	Tri::CTriMesh*		m_pTriMesh;
//	double				m_fCurvSum;

//	std::vector<double>	m_vFactors;

//public:
//	int					calcConformalFactor();
//	std::vector<double>	getFactors();

//protected:
//	void				calcDiscreteGaussianCurv(std::vector<double>& vGauCurv);
//	int					calcTargetUniformGaussianCurv(std::vector<double>& vTargetGauCurv);
//	double				cotan(Math::Vector3f v1, Math::Vector3f v2);
//	double				arccos(Math::Vector3f v1, Math::Vector3f v2);
//	void				normalization(std::vector<double>& vData);

//	void				debug();

public:
    CConformalFactor(){}
    ~CConformalFactor(){}

    QVector<double> calConformalFactor(SurfaceMeshModel *mesh);

private:
    void calcDiscreteGaussianCurv(SurfaceMeshModel *mesh, QVector<double>& vGauCurv, double &fCurvSum);
    int	calcTargetUniformGaussianCurv(SurfaceMeshModel *mesh, QVector<double>& vTargetGauCurv, double fCurvSum);

    bool has_halfedge(SurfaceMeshModel *mesh, Vertex start, Vertex end)
    {
        Halfedge h = mesh->halfedge(start);
        const Halfedge hh = h;
        if (h.is_valid()){
            do{
                if (mesh->to_vertex(h) == end) return true;
                h = mesh->cw_rotated_halfedge(h);
            } while (h != hh);
        }
        return false;
    }

    Eigen::SparseMatrix<double> Laplacian(SurfaceMeshModel *mesh)
    {
        Vector3VertexProperty points = mesh->vertex_coordinates();

        // Efficient sparse matrix construction
        typedef Eigen::Triplet<double> T;
        std::vector< T > L_c;

        // Fill as cotan operator
        foreach(Edge e, mesh->edges()){
            Scalar cot_alpha = 0, cot_beta = 0;

            Vertex vi = mesh->vertex(e, 0);
            Vertex vj = mesh->vertex(e, 1);

            Vertex v_a = mesh->to_vertex(mesh->next_halfedge(mesh->halfedge(e, 0)));
            if(has_halfedge(mesh, v_a, vj)) cot_alpha = (points[vi]-points[v_a]).dot(points[vj]-points[v_a]) / (points[vi]-points[v_a]).cross(points[vj]-points[v_a]).norm();

            Vertex v_b = mesh->to_vertex(mesh->next_halfedge(mesh->halfedge(e, 1)));
            if(has_halfedge(mesh, v_b, vi)) cot_beta  = (points[vi]-points[v_b]).dot(points[vj]-points[v_b]) / (points[vi]-points[v_b]).cross(points[vj]-points[v_b]).norm();

            Scalar cots = (0.5 * (cot_alpha + cot_beta));

            if(abs(cots) == 0) continue;

            L_c.push_back(T(vi.idx(), vj.idx(), -cots));
            L_c.push_back(T(vj.idx(), vi.idx(), -cots));
            L_c.push_back(T(vi.idx(), vi.idx(), cots));
            L_c.push_back(T(vj.idx(), vj.idx(), cots));
        }

        // Initialize a sparse matrix
        Eigen::SparseMatrix<Scalar> Lc_mat(mesh->n_vertices(), mesh->n_vertices());
        Lc_mat.setFromTriplets(L_c.begin(), L_c.end());

        return 2 * Lc_mat;
    }

    QVector<double> SolveSparseSystem(Eigen::SparseMatrix<double> &mat, QVector<double> &vDiffCurv)
    {
        Eigen::MatrixXd Linv;

        pinv(mat.toDense(), Linv);

        int n = vDiffCurv.size();
        Eigen::VectorXd diffCurv(n);
        for(int i = 0; i < n; i++){
            diffCurv[i] = vDiffCurv[i];
        }

        Eigen::VectorXd phi = Linv * diffCurv;
        n = phi.size();
        QVector<double> res;
        for(int i = 0; i < n; i++){
            res.push_back(phi[i]);
        }

        return res;
    }

    double cotan(Vector3 v1, Vector3 v2)
    {
        double ds = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
        double dt = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
        double dd = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

        // Divide by zero check
        double lenProduct = ds*dt;
        if(lenProduct < 1e-6f)
            lenProduct = 1e-6f;

        double dcos = dd/lenProduct;
        dcos = MIN(1.0,MAX(dcos,-1.0));
        double an = std::acos(dcos);
        double dsin = std::sin(an);

        // Divide by zero check
        double dcot = (dsin!=0.0) ? (dcos/dsin) : 1.0e6;

        return dcot;
    }

    double arccos(Vector3 v1, Vector3 v2)
    {
        double dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
        double d1 = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
        double d2 = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);

        // Divide by zero check
        double lenProduct = d1*d2;
        if(lenProduct < 1e-6f)
            lenProduct = 1e-6f;

        double dcos = dot/lenProduct;
        dcos = MIN(1, MAX(-1, dcos));
        return std::acos(dcos);
    }

    void normalization(QVector<double>& vData, int nVrt)
    {
        double fMin = 1.0e10;
        double fMax = -1.0e10;

        for (int iv=0; iv<nVrt; iv++)
        {
            double data = vData[iv];
            fMin = (fMin > data) ? data : fMin;
            fMax = (fMax < data) ? data : fMax;
        }

        double fRange = fMax - fMin;
        for (int iv=0; iv<nVrt; iv++)
        {
            double data = vData[iv];
            data = (data - fMin) / fRange;
            vData[iv] = data;
        }
    }

};

