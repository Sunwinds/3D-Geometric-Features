#include "BatchProc.h"
#include "SegMeshLoader.h"

#include "../bdf/redsvd.h"

#include "emd.h"

#include "RenderObjectExt.h"
#include "GeoHeatHelper.h"

#include "principal_curvature.h"
#include "Monge_via_jet_fitting.h"

#include <cmath>
#include <QFileDialog>
#include <QTextStream>
#include <QStringList>

void BatchProc::initParameters(RichParameterSet *pars)
{
    pars->addParam(new RichString("FilePath", "C:\\Project\\EASC\\data\\Container_off\\Large_off", "FilePath"));
    pars->addParam(new RichFloat("Delta", 2.0, "Delta"));
    // AGD
    pars->addParam(new RichBool("IsSubset", false, "Subset of vertices"));
    pars->addParam(new RichInt("NumSubset", 30, "Count"));
    // BDF
    pars->addParam(new RichInt("StartVertex", 0, "Start vertex"));
    pars->addParam(new RichBool("UsePsudeoInverse", true, "UsePsudeoInverse"));
    // Gaussian Curvature
    pars->addParam(new RichInt("Radius", 3, "Radius"));
    pars->addParam(new RichInt("Smoothing", 0, "Smoothing"));
    pars->addParam(new RichBool("Fitting", false, "By fitting"));
    pars->addParam(new RichBool("Clamp", true, "Clamp"));
    pars->addParam(new RichInt("ActiveItem", 0, "Item (mean, ..etc)"));
    // SDF
    pars->addParam(new RichInt("numCones", 4, "Num. cones"));
    pars->addParam(new RichInt("raysInCone", 8, "Rays in cone"));
    pars->addParam(new RichFloat("coneSeparation", 20.0f,"Cone separation"));
    pars->addParam(new RichBool("gaussianWeights", true, "Gaussian weights"));
    pars->addParam(new RichInt("Smooth", 1, "Smooth"));
}

void BatchProc::applyFilter(RichParameterSet *pars)
{
    //FilterPlugin * gd_plugin = pluginManager()->getFilter("Batch Process");

    BatchloadMesh(pars->getString("FilePath"));

    calAffintyMatrix(pars);

    //qApp->processEvents();

    //gd_plugin->applyFilter( NULL );
}

void BatchProc::BatchloadMesh(QString filepath)
{
    QStringList filters, filelist;
	QString ext = filepath.right(3);
	QString fileName;

	SegMeshLoader sLoader;

	if(ext == "obj"){
		filters << "*.obj";
		filelist = QDir(filepath).entryList(filters);

		int mSize = filelist.size();

		QString id = filelist[0].left(3);
		SegMesh sMesh;
		for(int i = 0; i < mSize; i++){
			if(filelist[i].left(3) == id){
				fileName = filepath + "\\" + filelist[i];
				sMesh.push_back(sLoader.readOBJ(fileName));
			}
			else{
				id = filelist[i].left(3);
				meshList.push_back(sMesh);
				sMesh.clear();
				fileName = filepath + "\\" + filelist[i];
				sMesh.push_back(sLoader.readOBJ(fileName));
			}
		}
	}
	else if(ext == "off"){
		filters << "*.off";
		filelist = QDir(filepath).entryList(filters);

		int mSize = filelist.size();

		//qDebug() << "Start Loading" << filepath << mSize;

		for (int i = 0; i < mSize; i++)
		{
			fileName = filepath + "\\" + filelist[i];
			//qDebug() << fileName;
			meshList.push_back(loadMesh(fileName));
		}
	}
    //qDebug() << "#Mesh" << meshList.size();
}

SegMesh BatchProc::loadMesh(QString filename)
{
    SegMeshLoader sLoader;
    //SegMesh sMesh;
    QString segF = filename.mid(0, filename.size() - 3) + "seg";
    //qDebug() << segF;
    return sLoader.loadGroupMeshFromOFF(filename, segF);
}

void BatchProc::calAffintyMatrix(RichParameterSet *pars)
{
    // Calculate feature descriptor
	int nMesh = meshList.size();
    QVector<QVector<QVector<double>>> descs;
	QVector<int> mIdperPart;
	//descs.resize(nMesh);
	int nPart = 0;
    
    for(int i = 0; i < nMesh; i++){
        int nSeg = meshList[i].size();
		//descs[i].resize(nSeg);
		nPart += nSeg;
        for(int j = 0; j < nSeg; j++){
            qDebug() << "Shape" << i << "Seg" << j;
			mIdperPart.push_back(i);
            descs.push_back(calDescriptor(meshList[i][j], pars));
        }
    }

    double delta = pars->getFloat("Delta");

    // Write to file
    QFile file("AffinityMatrix.csv");
    file.open(QFile::WriteOnly|QFile::Truncate);
    QTextStream stream(&file);

    // Calculate affinity matrix
    QVector<QVector<double>> mat;
    //int nSize = descs.size();
    mat.resize(nPart);
    for(int i = 0; i < nPart; i++){
        mat[i].resize(nPart);
        for(int j = 0; j < nPart; j++){
			if(mIdperPart[i] == mIdperPart[j] && i != j){
				mat[i][j] = 0.0;
			}
			else if(i == j){ 
				mat[i][j] = 1.0;
			}
			else{
				mat[i][j] = std::pow(NATURALEXP, -calDist(descs[i], descs[j]) / (2.0*delta));
			}
            stream << mat[i][j] << ',' ; // this writes first line with two columns
        }
        stream << "\n";
    }
    file.close();
}

float dist(feature_t *F1, feature_t *F2)
{
   return ((*F1) - (*F2))*((*F1) - (*F2));
}

double BatchProc::calEMD(QVector<double> &d1, QVector<double> &d2)
{
    int size = d1.size();
    QVector<float> w(size, 0.01);
    signature_t s1 = {size, &(d1[0]), &(w[0])}, s2 = {size, &(d2[0]), &(w[0])};
    return emd(&s1, &s2, dist, 0, 0);
}

double BatchProc::calDist(QVector<QVector<double>> &d1, QVector<QVector<double>> &d2)
{
    // Add EMD
    double dist = 0.0;
    int nSize = d1.size();
    for(int i = 0; i < nSize - 1; i++){
        dist += calEMD(d1[i], d2[i]);
        //dist += (d1[i] - d2[i]) * (d1[i] - d2[i]);
    }
    // Add PCA descriptors
    for(int i = 0; i < 3; i++){
        dist += (d1[nSize-1][i] - d2[nSize-1][i]) * (d1[nSize-1][i] - d2[nSize-1][i]);
    }
    dist = std::sqrt(dist);
    return dist;
}

QVector<QVector<double>> BatchProc::calDescriptor(SurfaceMeshModel *mesh, RichParameterSet *pars)
{
    QVector<QVector<double>> desc;
//    RichParameterSet *pars1;
//    RichParameterSet *pars2;
//    RichParameterSet *pars3;
//    RichParameterSet *pars4;
    desc.push_back(calAGDDesc(mesh, pars));
    qDebug() << "AGD done!" ;
    desc.push_back(calBDFDesc(mesh, pars));
    qDebug() << "BDF done!" ;
    desc.push_back(calSDFDesc(mesh, pars));
    qDebug() << "SDF done!";
    desc.push_back(calCurvatureDesc(mesh, pars));
    qDebug() << "Curvature done!";
    desc.push_back(calConformalFactor(mesh));
    qDebug() << "Conformal Factor done!";
    // Special Case, not using EMD to measure
    desc.push_back(calPCADesc(mesh));
    qDebug() << "PCA done!";
    return desc;
}

// AGD
void BatchProc::geoheatHelper(SurfaceMeshModel *mesh, int srcVert)
{
	QSet<Vertex> src;

	src.insert(Vertex(srcVert));
	// Source points are tagged on mesh vertices
	//BoolVertexProperty src_points = mesh->get_vertex_property<bool>("v:geo_src_points");
	//if(!src_points.is_valid()) return;
	//foreach(Vertex v, mesh->vertices())
	//	if(src_points[v]) src.insert(v);

	GeoHeatHelper * h = new GeoHeatHelper(mesh);

	ScalarVertexProperty distance = h->getUniformDistance(src);

	delete h;
}

QVector<double> BatchProc::calAGDDesc(SurfaceMeshModel *mesh, RichParameterSet *pars)
{
    QVector<double> desc(100,0.0);

	/*if(application()->document()->selectedModel())
	application()->document()->deleteModel(application()->document()->selectedModel());
	application()->document()->addModel(mesh);
	application()->document()->setSelectedModel(mesh);*/

    //FilterPlugin * gd_plugin = pluginManager()->getFilter("Geodesic distance: heat kernel");

    BoolVertexProperty src_points = mesh->vertex_property<bool>("v:geo_src_points", false);
    ScalarVertexProperty geo_dist = mesh->vertex_property<Scalar>("v:uniformDistance", 0);

    // Consider a subset of vertices
    if(pars && pars->getBool("IsSubset"))
    {
        // Starting point
        src_points[Vertex(0)] = true;

        int numIter = pars->getInt("NumSubset");

        for(int itr = 0; itr < numIter; itr++)
        {
            //gd_plugin->applyFilter( NULL );
			geoheatHelper(mesh, 0);

            // Find furthest point, add it to source set
            double maxDist = -DBL_MAX;
            int maxIDX = -1;
            foreach(Vertex v, mesh->vertices())
            {
                if(geo_dist[v] > maxDist){
                    maxDist = geo_dist[v];
                    maxIDX = v.idx();
                }
            }

            src_points[Vertex(maxIDX)] = true;

            //mainWindow()->setStatusBarMessage( QString("iter (%1) out of (%2)").arg(itr).arg(numIter)  );

            //qApp->processEvents();
        }
    }
    else
    {
        foreach(Vertex v, mesh->vertices()) src_points[v] = true;
    }

    //qApp->processEvents();

    // Compute vertex-to-all distances
    ScalarVertexProperty all_dist = mesh->vertex_property<Scalar>("v:agd", 0);

    //RichParameterSet * gd_params = new RichParameterSet;
    //gd_plugin->initParameters( gd_params );

    foreach(Vertex v, mesh->vertices())
    {
        if(!src_points[v]) continue;

        //gd_params->setValue("sourceVertex", int( v.idx() ));
        //gd_plugin->applyFilter( gd_params );
		geoheatHelper(mesh, int( v.idx() ));
        ScalarVertexProperty geo_dist = mesh->vertex_property<Scalar>("v:uniformDistance", 0);

        foreach(Vertex v, mesh->vertices())
        {
            all_dist[v] += geo_dist[v];
        }
    }

    double minVal = DBL_MAX;
    double maxVal = -DBL_MAX;

    foreach(Vertex v, mesh->vertices())
    {
        minVal = qMin( minVal, all_dist[v] );
        maxVal = qMax( maxVal, all_dist[v] );
    }

    foreach(Vertex v, mesh->vertices())
    {
        all_dist[v] = (all_dist[v] - minVal) / (maxVal - minVal);
    }

    int N = mesh->n_vertices();

    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        Vertex v(i);
        for(int i = 0; i < 100; i++){
            double lb = float(i) / 100.0;
            if(all_dist[v] >= lb && all_dist[v] < lb+0.01){
                desc[i]++;
                break;
            }
        }
    }

    for(int i = 0; i < 100; i++){
        desc[i] /= float(N);
    }
    return desc;
}

// BDF
Eigen::SparseMatrix<double> BatchProc::Lc(SurfaceMeshModel *mesh)
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

bool BatchProc::has_halfedge(SurfaceMeshModel *mesh, Vertex start, Vertex end){
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

template<typename _Matrix_Type_>
bool pinv(const _Matrix_Type_ &a, _Matrix_Type_ &result, double	epsilon = std::numeric_limits<typename _Matrix_Type_::Scalar>::epsilon())
{
    if(a.rows() < a.cols())	return false;
    Eigen::JacobiSVD< _Matrix_Type_ > svd = a.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
    typename _Matrix_Type_::Scalar tolerance = epsilon * std::max(a.cols(),	a.rows()) * svd.singularValues().array().abs().maxCoeff();
    result = svd.matrixV() * ( (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0)).matrix().asDiagonal() * svd.matrixU().adjoint();
    return true;
}

QVector<double> BatchProc::calBDFDesc(SurfaceMeshModel *mesh, RichParameterSet *pars)
{
    QVector<double> desc(100,0.0);

    // Vector3VertexProperty points = mesh->vertex_coordinates();

    // Laplace matrix
    Eigen::SparseMatrix<double> L = Lc(mesh);

    // Vertex areas
    Eigen::VectorXd AV = Eigen::VectorXd::Zero(mesh->n_vertices());
    SurfaceMeshHelper h(mesh);
    ScalarFaceProperty farea = h.computeFaceAreas();
    foreach(Face f, mesh->faces()) foreach(Vertex v, mesh->vertices(f)) AV[v.idx()] += (2.0 * farea[f]);

    // Diagonal entries
    for(unsigned int i = 0; i < AV.rows(); i++) AV[i] = 6.0 / AV[i];

    // System
    Eigen::SparseMatrix<double> X = L * AV.asDiagonal() * L;

    qApp->processEvents();

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

    // Distance matrix
    MatrixXd dB( mesh->n_vertices(), mesh->n_vertices() );
    for(unsigned int i = 0; i < mesh->n_vertices(); i++)
        for(unsigned int j = 0; j < mesh->n_vertices(); j++)
            dB(i,j) = sqrt( GD(i,i) + GD(j,j) - 2.0 * GD(i,j) );

    dB /= dB.maxCoeff();

    int startVertex = pars->getInt("StartVertex");

    double maxVal = -DBL_MAX, minVal = DBL_MAX;
    foreach(Vertex v, mesh->vertices()) { maxVal = qMax(dB(startVertex, v.idx()), maxVal); minVal = qMin(dB(startVertex, v.idx()), minVal);}
    foreach(Vertex v, mesh->vertices()) dB(startVertex, v.idx()) = (dB(startVertex, v.idx()) - minVal) / (maxVal - minVal);

    int N = mesh->n_vertices();

    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        Vertex v(i);
        for(int i = 0; i < 100; i++){
            double lb = float(i) / 100.0;
            if(dB(startVertex, v.idx()) >= lb && dB(startVertex, v.idx()) < lb+0.01){
                desc[i]++;
                break;
            }
        }
    }

//    foreach(Vertex v, mesh->vertices()){
//        for(int i = 0; i < 100; i++){
//            double lb = float(i) / 100.0;
//            if( dB(startVertex, v.idx()) >= lb && dB(startVertex, v.idx()) < lb+0.01){
//                desc[i]++;
//                break;
//            }
//        }
//    }
    for(int i = 0; i < 100; i++){
        desc[i] /= float(N);
    }

    return desc;
}

// SDF
#define DEG2RAD 0.01745329

Vector3 orthogonalVector(const Vector3& n) {
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

double BatchProc::compute_volume(const Vector3 & p, const Vector3 & n, int numCones,
                           double coneSeparation, int raysInCone, bool gaussianWeights,
                           Octree * octree, QList<Ray> * allRays)
{
    double result = 0.0;

    double distanceCounter = 0.0;
    double intersectionCounter = 0.0;
    double distance = DBL_MAX;

    //qDebug() << "Vol computation";
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

QVector<double> BatchProc::calSDFDesc(SurfaceMeshModel *mesh, RichParameterSet *pars)
{
    QVector<double> desc(100,0.0);

    int N = mesh->n_vertices();

    Vector3VertexProperty points = mesh->vertex_coordinates();
	mesh->update_vertex_normals();
    Vector3VertexProperty normals = mesh->vertex_normals();

    ScalarVertexProperty m_sdf = mesh->vertex_property<Scalar>("v:sdf");

    // Parameters
    int m_numCones = pars->getInt("numCones");
    int m_raysInCone = pars->getInt("raysInCone");
    double m_coneSeperation = pars->getFloat("coneSeparation");
    bool m_gaussianWeights = pars->getBool("gaussianWeights");

    //qDebug() << "Start SDF" << N;

    Octree octree( mesh );
    //qDebug() << "Octree built";

    // Visualization
    //QList<Ray> allRays;
    //starlab::VectorSoup * vs = new starlab::VectorSoup();

    //#pragma omp parallel for
    for(int i = 0; i < N; i++)
    {
        //QList<Ray> curRays;
        Vertex v(i);
        //qDebug() << i;
        m_sdf[v] = compute_volume( points[v], normals[v],
            m_numCones, m_coneSeperation, m_raysInCone, m_gaussianWeights, &octree /*, curRays*/ );

        //allRays += curRays;
        //foreach(Ray r, allRays) vs->addVector(r.origin, r.direction);
    }

    // Smooth the function
    SurfaceMeshHelper(mesh).smoothVertexProperty<Scalar>("v:sdf", pars->getInt("Smooth"));

    double maxVal = -DBL_MAX, minVal = DBL_MAX;
    foreach(Vertex v, mesh->vertices()) { maxVal = qMax(m_sdf[v], maxVal); minVal = qMin(m_sdf[v], minVal);}
    foreach(Vertex v, mesh->vertices()) m_sdf[v] = (m_sdf[v] - minVal) / (maxVal - minVal);

    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        Vertex v(i);
        for(int i = 0; i < 100; i++){
            double lb = float(i) / 100.0;
            if( m_sdf[v] >= lb && m_sdf[v] < lb+0.01){
                desc[i]++;
                break;
            }
        }
    }
    for(int i = 0; i < 100; i++){
        desc[i] /= float(N);
    }

    return desc;
}

// Curvature
#ifdef WIN32
namespace std{  bool isnan(double x){ return _isnan(x); } }
#endif

void BatchProc::collectEnoughRings(SurfaceMeshModel *mesh, Vertex v, const size_t min_nb, std::vector<int> & all)
{
    std::vector<Vertex> current_ring, next_ring;
    SurfaceMeshModel::Vertex_property<int> visited_map = mesh->vertex_property<int>("v:visit_map", -1);

    foreach(Vertex v, mesh->vertices()) visited_map[v] = -1;

    //initialize
    visited_map[v] = 0;
    current_ring.push_back(v);
    all.push_back(v.idx());

    int i = 1;

    while ( (all.size() < min_nb) &&  (current_ring.size() != 0) )
    {
        // collect ith ring
        std::vector<Vertex>::iterator it = current_ring.begin(), ite = current_ring.end();

        for(;it != ite; it++){
            SurfaceMeshModel::Halfedge_around_vertex_circulator hedgeb = mesh->halfedges(*it), hedgee = hedgeb;

            do{
                Vertex vj = mesh->to_vertex(hedgeb);

                if (visited_map[vj] == -1){
                    visited_map[vj] = i;
                    next_ring.push_back(vj);
                    all.push_back(vj.idx());

                    if(all.size() > min_nb) return;
                }

                ++hedgeb;
            } while(hedgeb != hedgee);
        }

        //next round must be launched from p_next_ring...
        current_ring = next_ring;
        next_ring.clear();

        i++;
    }
}

QVector<double> BatchProc::calCurvatureDesc(SurfaceMeshModel *mesh, RichParameterSet *pars)
{
    QVector<double> desc(100,0.0);

    Eigen::MatrixXd V(mesh->n_vertices(), 3);
    Eigen::MatrixXi F(mesh->n_faces(), 3);

    Eigen::MatrixXd NV(mesh->n_vertices(), 3);
    Eigen::MatrixXd NF(mesh->n_faces(), 3);

    Vector3VertexProperty points = mesh->vertex_coordinates();
    Vector3VertexProperty normals = mesh->vertex_normals(true);
    Vector3FaceProperty fnormals = mesh->face_normals(true);

    // Fill in vertices
    int vi = 0;
    foreach(Vertex v, mesh->vertices()) { V.row(vi) = points[v]; NV.row(vi) = normals[v]; vi++; }

    // Fill in triangles
    int fi = 0;
    foreach(Face f, mesh->faces()){
        int i = 0;
        foreach(Vertex v, mesh->vertices(f)) F(fi,i++) = v.idx();
        NF.row(fi) = fnormals[f];
        fi++;
    };

    // Save results to mesh
    ScalarVertexProperty k1 = mesh->vertex_property("v:k1", 0.0);
    ScalarVertexProperty k2 = mesh->vertex_property("v:k2", 0.0);
    Vector3VertexProperty pd1 = mesh->vertex_property("v:pd1", Vector3(0,0,0));
    Vector3VertexProperty pd2 = mesh->vertex_property("v:pd2", Vector3(0,0,0));

    ScalarVertexProperty cmax = mesh->vertex_property("v:cmax", 0.0);
    ScalarVertexProperty cmin = mesh->vertex_property("v:cmin", 0.0);
    ScalarVertexProperty cmean = mesh->vertex_property("v:cmean", 0.0);
    ScalarVertexProperty cgauss = mesh->vertex_property("v:cgauss", 0.0);

    // Compute average edge length
    double avgEdge = 0;
    SurfaceMeshHelper h(mesh);
    ScalarEdgeProperty elengths = h.computeEdgeLengths();
    foreach(Edge e, mesh->edges()) avgEdge += elengths[e];
    avgEdge /= mesh->n_edges();

    if( !pars->getBool("Fitting") )
    {
        // Temporary results
        Eigen::MatrixXd K1(mesh->n_vertices(), 1);
        Eigen::MatrixXd K2(mesh->n_vertices(), 1);

        Eigen::MatrixXd PD1(mesh->n_vertices(), 3);
        Eigen::MatrixXd PD2(mesh->n_vertices(), 3);

        igl::principal_curvature( V, F, NV, NF, K1, K2, PD1, PD2, pars->getInt("Radius"), 0 );

        foreach(Vertex v, mesh->vertices()){
            int vidx = v.idx();

            k1[v] = K1(vidx, 0);
            k2[v] = K2(vidx, 0);

            pd1[v] = PD1.row(vidx);
            pd2[v] = PD2.row(vidx);
        }
    }
    else
    {
        std::vector< Vector3 > pnts;
        std::vector< std::vector<int> > neighbours(mesh->n_vertices());
        std::vector< std::vector<int> > vertex_to_vertices;
        //double scaledRadius = avgEdge * pars->getInt("Radius");
        //igl::adjacency_list(F, vertex_to_vertices);

        foreach(Vertex v, mesh->vertices())
        {
            // Do a sphere search to collect neighbours
            //igl::CurvatureCalculator::getSphere(V, vertex_to_vertices, v.idx(), scaledRadius, neighbours[v.idx()], 6);
            collectEnoughRings(mesh, v, qMax(7, 2 * pars->getInt("Radius")), neighbours[v.idx()]);

            pnts.push_back(points[v]);
        }

        // Fit around each point
        std::vector<Monge_via_jet_fitting::Monge_form> fit = Monge_via_jet_fitting::fit_points(pnts, neighbours);

        // Extract results
        foreach(Vertex v, mesh->vertices())
        {
            Monge_via_jet_fitting::Monge_form m = fit[v.idx()];
            m.comply_wrt_given_normal(normals[v]);

            k1[v] = m.principal_curvatures(0);
            k2[v] = m.principal_curvatures(1);

            pd1[v] = m.maximal_principal_direction() * k1[v];
            pd2[v] = m.minimal_principal_direction() * k2[v];
        }
    }

    std::vector<double> maxVals(4, -DBL_MAX);
    std::vector<double> minVals(4,  DBL_MAX);

    double minVal = DBL_MAX;
    double maxVal = -DBL_MAX;

    foreach(Vertex v, mesh->vertices()){
        if(std::isnan(pd1[v][0])) pd1[v] = Vector3(0,0,0);
        if(std::isnan(pd2[v][0])) pd2[v] = Vector3(0,0,0);

        cmax[v] = k1[v];
        cmin[v] = k2[v];
        cmean[v] = (cmax[v] + cmin[v]) / 2.0;
        cgauss[v] = cmin[v] * cmax[v];
        minVal = qMin( minVal, cgauss[v] );
        maxVal = qMax( maxVal, cgauss[v] );
    }

    foreach(Vertex v, mesh->vertices())
    {
        cgauss[v] = (cgauss[v] - minVal) / (maxVal - minVal);
    }

//    foreach(Vertex v, mesh->vertices()){
//        for(int i = 0; i < 100; i++){
//            double lb = float(i) / 100.0;
//            if( cgauss[v] >= lb && cgauss[v] < lb+0.01){
//                desc[i]++;
//                break;
//            }
//        }
//    }
    int N = mesh->n_vertices();

    #pragma omp parallel for
    for(int i = 0; i < N; i++){
        Vertex v(i);
        for(int i = 0; i < 100; i++){
            double lb = float(i) / 100.0;
            if( cgauss[v] >= lb && cgauss[v] < lb+0.01){
                desc[i]++;
                break;
            }
        }
    }

    for(int i = 0; i < 100; i++){
        desc[i] /= float(N);
    }

    return desc;
}

// PCA
QVector<double> BatchProc::calPCADesc(SurfaceMeshModel *mesh)
{
    QVector<double> desc;

    int n = mesh->vertices_size();
    VertexIterator begin = mesh->vertices_begin();
    VertexIterator end = mesh->vertices_end();
    Points pts = mesh->vertex_property<Point>(VPOINT);

    FT x, y, z,
        sumX = 0., sumY = 0., sumZ = 0.,
        sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
        sumXY = 0., sumXZ = 0., sumYZ = 0.,
        xx, yy, zz, xy, xz, yz;

    for (; begin != end; begin++){
        x = pts[begin][0];
        y = pts[begin][1];
        z = pts[begin][2];
        sumX += x / n;
        sumY += y / n;
        sumZ += z / n;
        sumX2 += x * x / n;
        sumY2 += y * y / n;
        sumZ2 += z * z / n;
        sumXY += x * y / n;
        sumXZ += x * z / n;
        sumYZ += y * z / n;
    }
    xx = sumX2 - sumX * sumX;
    yy = sumY2 - sumY * sumY;
    zz = sumZ2 - sumZ * sumZ;
    xy = sumXY - sumX * sumY;
    xz = sumXZ - sumX * sumZ;
    yz = sumYZ - sumY * sumZ;

    Eigen::Matrix3d covariance;
    covariance << xx, xy, xz,
        xy, yy, yz,
        xz, yz, zz;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(covariance);
    Eigen::Vector3d eigenVals = es.eigenvalues();

    Vector3 scalar;

    FT sum = eigenVals[0] + eigenVals[1] + eigenVals[2];
    scalar[0] = (eigenVals[2] - eigenVals[1]) / sum;
    scalar[1] = 2 * (eigenVals[1] - eigenVals[0]) / sum;
    scalar[2] = 3 * eigenVals[0]/ sum;

    desc.push_back(scalar[0]);
    desc.push_back(scalar[1]);
    desc.push_back(scalar[2]);

    return desc;
}

double BatchProc::cotan(Vector3 v1, Vector3 v2)
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

double BatchProc::arccos(Vector3 v1, Vector3 v2)
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

void BatchProc::normalization(QVector<double> &vData, int nVrt)
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

void BatchProc::calcDiscreteGaussianCurv(SurfaceMeshModel *mesh, QVector<double> &vGauCurv, double &fCurvSum)
{
    fCurvSum = 0.0;
    double _2Pi = 2.0 * Pi;
    Vector3VertexProperty points = mesh->vertex_coordinates();

    int j = 0;
    foreach(Vertex v, mesh->vertices()){
        double fAngles = 0.0;
        Vector3 vCenter = points[v];

        // Collect neighbors
        Surface_mesh::Vertex_around_vertex_circulator vvit, vvend;
        vvit = vvend = mesh->vertices(v);
        QVector<Vertex> vits;

        do{
            vits.push_back(vvit);
        } while(++vvit != vvend);

        int nNbrVrt = vits.size();
        //qDebug() << "Start computation" << j << ' ' << nNbrVrt;

        for(int i=0; i<nNbrVrt; i++)
        {
            Vector3 v1 = points[vits[i]] - vCenter;
            Vector3 v2 = points[vits[(i+1+nNbrVrt) % nNbrVrt]] - vCenter;
            fAngles += arccos(v1, v2);
            //qDebug() << "#Iter = " << i;
        }

        vGauCurv[j] = _2Pi - fAngles;
        fCurvSum += vGauCurv[j];
        j++;

    }
}

int BatchProc::calcTargetUniformGaussianCurv(SurfaceMeshModel *mesh, QVector<double> &vTargetGauCurv, double fCurvSum)
{
    SurfaceMeshHelper h(mesh);
    Scalar area = 0.0;
    ScalarFaceProperty  fArea = h.computeFaceAreas();

    foreach(Face f, mesh->faces()){
        area += fArea[f];
    }
    double fCurv = fCurvSum / area;

    int i = 0;
    foreach(Vertex v, mesh->vertices())
    {
        double fInfluArea = 0.0;
        Surface_mesh::Face_around_vertex_circulator fvit, fvend;
        fvit = fvend = mesh->faces(v);
        do{
           fInfluArea += fArea[fvit];
        } while(++fvit != fvend);

        fInfluArea /= 3.0;
        vTargetGauCurv[i] = fCurv * fInfluArea;
        i++;
    }
    return 1;
}

QVector<double> BatchProc::SolveSparseSystem(Eigen::SparseMatrix<double> &mat, QVector<double> &vDiffCurv)
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

QVector<double> BatchProc::calConformalFactor(SurfaceMeshModel *mesh)
{
    int nVrt = mesh->n_vertices();
    double fCurvSum;

    // discrete gaussian curvature
    QVector<double> vGauCurv(nVrt, 0.0);
    calcDiscreteGaussianCurv(mesh, vGauCurv, fCurvSum);

    // target gaussian curvature
    QVector<double> vTargetGauCurv(nVrt, 0.0);
    calcTargetUniformGaussianCurv(mesh, vTargetGauCurv, fCurvSum);

    // K_T - K_org
    QVector<double> vDiffCurv(nVrt, 0.0);
    for (int iv=0; iv<nVrt; iv++)
    {
        vDiffCurv[iv] = vTargetGauCurv[iv] - vGauCurv[iv];
    }

    // solve problem
    Eigen::SparseMatrix<double> matLaplacian = Lc(mesh);
    QVector<double> Phi = SolveSparseSystem(matLaplacian, vDiffCurv);

    // normalization
    normalization(Phi, nVrt);

    QVector<double> desc(100,0.0);

    #pragma omp parallel for
    for(int j = 0; j < nVrt; j++){
        for(int i = 0; i < 100; i++){
            double lb = float(i) / 100.0;
            if( Phi[j] >= lb && Phi[j] < lb+0.01){
                desc[i]++;
                break;
            }
        }
    }
    for(int i = 0; i < 100; i++){
        desc[i] /= float(nVrt);
    }

    return desc;
}






