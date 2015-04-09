#include <QElapsedTimer>

#include "bdf.h"
#include "SurfaceMeshHelper.h"
#include "RenderObjectExt.h"

#include "redsvd.h"

void bdf::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichInt("StartVertex", 0, "Start vertex"));
	pars->addParam(new RichBool("UsePsudeoInverse", true, "UsePsudeoInverse"));
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

bool bdf::has_halfedge(Vertex start, Vertex end){
	Halfedge h = mesh()->halfedge(start);
	const Halfedge hh = h;
	if (h.is_valid()){
		do{
			if (mesh()->to_vertex(h) == end) return true;
			h = mesh()->cw_rotated_halfedge(h);
		} while (h != hh);
	}
	return false;
}

Eigen::SparseMatrix<double> bdf::Lc()
{
	Vector3VertexProperty points = mesh()->vertex_coordinates();

	// Efficient sparse matrix construction
	typedef Eigen::Triplet<double> T;
	std::vector< T > L_c;

	// Fill as cotan operator
	foreach(Edge e, mesh()->edges()){
		Scalar cot_alpha = 0, cot_beta = 0;

		Vertex vi = mesh()->vertex(e, 0);
		Vertex vj = mesh()->vertex(e, 1);

		Vertex v_a = mesh()->to_vertex(mesh()->next_halfedge(mesh()->halfedge(e, 0)));
		if(has_halfedge(v_a, vj)) cot_alpha = (points[vi]-points[v_a]).dot(points[vj]-points[v_a]) / (points[vi]-points[v_a]).cross(points[vj]-points[v_a]).norm();

		Vertex v_b = mesh()->to_vertex(mesh()->next_halfedge(mesh()->halfedge(e, 1)));
		if(has_halfedge(v_b, vi)) cot_beta  = (points[vi]-points[v_b]).dot(points[vj]-points[v_b]) / (points[vi]-points[v_b]).cross(points[vj]-points[v_b]).norm();

		Scalar cots = (0.5 * (cot_alpha + cot_beta));

		if(abs(cots) == 0) continue;

		L_c.push_back(T(vi.idx(), vj.idx(), -cots));
		L_c.push_back(T(vj.idx(), vi.idx(), -cots));
		L_c.push_back(T(vi.idx(), vi.idx(), cots));
		L_c.push_back(T(vj.idx(), vj.idx(), cots));
	}

	// Initialize a sparse matrix
	Eigen::SparseMatrix<Scalar> Lc_mat(mesh()->n_vertices(), mesh()->n_vertices());
	Lc_mat.setFromTriplets(L_c.begin(), L_c.end());

	return 2 * Lc_mat;
}

inline void printMatrix(Eigen::MatrixXd L, QString name = "")
{
	qDebug() << name << "\n";

	for(int i = 0; i < L.rows(); i++){
		QStringList line;
		for(int j = 0; j < L.cols(); j++) line << QString::number( L.coeffRef(i,j) );
		qDebug() << qPrintable(line.join(", "));
	}

	qDebug() << "\n";
}

template<typename _Matrix_Type_>
bool pinv(const _Matrix_Type_ &a, _Matrix_Type_ &result, double	epsilon = std::numeric_limits<typename _Matrix_Type_::Scalar>::epsilon())
{
	if(a.rows() < a.cols())	return false;
	Eigen::JacobiSVD< _Matrix_Type_ > svd = a.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	typename _Matrix_Type_::Scalar tolerance = epsilon * std::max(a.cols(),	a.rows()) * svd.singularValues().array().abs().maxCoeff();
	result = svd.matrixV() * ( (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0)).matrix().asDiagonal() * svd.matrixU().adjoint();
	return true;
}

void bdf::applyFilter(RichParameterSet *pars)
{
	drawArea()->clear();

	Vector3VertexProperty points = mesh()->vertex_coordinates();

	// Laplace matrix
	Eigen::SparseMatrix<double> L = Lc();

	// Vertex areas
	Eigen::VectorXd AV = Eigen::VectorXd::Zero(mesh()->n_vertices());
	SurfaceMeshHelper h(mesh());
	ScalarFaceProperty farea = h.computeFaceAreas();
    foreach(Face f, mesh()->faces()) foreach(Vertex v, mesh()->vertices(f)) AV[v.idx()] += (2.0 * farea[f]);

	// Diagonal entries
	for(unsigned int i = 0; i < AV.rows(); i++) AV[i] = 6.0 / AV[i];

	// System
    Eigen::SparseMatrix<double> X = L * AV.asDiagonal() * L;

	mainWindow()->setStatusBarMessage("System built, now solving...");
	qApp->processEvents();
	QElapsedTimer timer; timer.start();

	MatrixXd GD;

	bool isUsePsudeoInverse = pars->getBool("UsePsudeoInverse");

	if( isUsePsudeoInverse )
	{	
        pinv(X.toDense(), GD);
	}
	else
	{
		MatrixXf gd;
		REDSVD::RedSVD svd( X.cast<float>() );
		VectorXf ds = svd.singularValues();
		for(unsigned int i = 0; i < ds.rows(); i++) if(ds[i] > 1e-8) ds[i] = 1.0 / ds[i];
		gd = svd.matrixV() * ds.asDiagonal() * svd.matrixU().transpose();
		GD = gd.cast<double>();
	}

	// DEBUG:
	//printMatrix(GD, "gd");

    //mainWindow()->setStatusBarMessage(QString("pinv time (%1 ms)").arg( timer.elapsed() ));

	// Distance matrix
	MatrixXd dB( mesh()->n_vertices(), mesh()->n_vertices() );
	for(unsigned int i = 0; i < mesh()->n_vertices(); i++)
		for(unsigned int j = 0; j < mesh()->n_vertices(); j++)
			dB(i,j) = sqrt( GD(i,i) + GD(j,j) - 2.0 * GD(i,j) );

	dB /= dB.maxCoeff();

	// Visualize
	if(pars->getBool("Visualize"))
    {
		int startVertex = pars->getInt("StartVertex");

        foreach(Vertex v, mesh()->vertices())
		{
            QColor c = starlab::qtJetColor(dB(startVertex, v.idx()));
			drawArea()->drawPoint( points[v], 8, c );
		}
	}
}
