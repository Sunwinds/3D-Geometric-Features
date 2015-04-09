#pragma once

//#include <QWidget>
//#include "ui_BatchProc.h"

//class BatchProcWidget : public QWidget
//{
//    Q_OBJECT

//public:
//    explicit BatchProcWidget(FdPlugin *fp, QWidget *parent = 0);
//    ~BatchProcWidget();

//private:
//    Ui::BatchProc *ui;
//    BatchProc *plugin;
//    QString filePath;

//public slots:
//    void loadData();
//    void
//};

#include "interfaces/ModePluginDockWidget.h"

#include "SurfaceMeshPlugins.h"
#include "SurfaceMeshHelper.h"
#include "RichParameterSet.h"
#include "SurfaceMeshModel.h"

#include "Octree.h"

#include <QVector>
#include <QString>

#include "Eigen/LU"
#include "Eigen/Geometry"
#include "Eigen/Eigenvalues"
#include <Eigen/Sparse>

#define NATURALEXP 2.71828
#define Pi 3.14159265
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

typedef QVector<SurfaceMeshModel *> SegMesh;
typedef double											FT;
typedef Surface_mesh::Vertex_iterator                   VertexIterator;
typedef Surface_mesh::Vertex_property<Surface_mesh::Point> Points;

class BatchProc: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "BatchProc.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Batch Process"; }
    QString description() { return "Process done"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);

private:
    void BatchloadMesh(QString filepath);
    SegMesh loadMesh(QString filename);

    // Calculate Affinity Matrix and
    void calAffintyMatrix(RichParameterSet *pars);
    double calEMD(QVector<double> &d1, QVector<double> &d2);
    double calDist(QVector<QVector<double>> &d1, QVector<QVector<double>> &d2);
    QVector<QVector<double>> calDescriptor(SurfaceMeshModel *mesh, RichParameterSet *pars);

    QVector<SegMesh> meshList;

    // AGD
	void geoheatHelper(SurfaceMeshModel *mesh, int srcVert);
    QVector<double> calAGDDesc(SurfaceMeshModel *mesh, RichParameterSet *pars);
    // BDF
    Eigen::SparseMatrix<double> Lc(SurfaceMeshModel *mesh);
    bool has_halfedge(SurfaceMeshModel *mesh, Vertex start, Vertex end);
    QVector<double> calBDFDesc(SurfaceMeshModel *mesh, RichParameterSet *pars);
    // SDF
    double compute_volume(const Vector3 & p, const Vector3 & n,
           int numCones, double coneSeparation, int raysInCone, bool gaussianWeights,
           Octree * octree, QList<Ray> * allRays = 0);
    QVector<double> calSDFDesc(SurfaceMeshModel *mesh, RichParameterSet *pars);
    // Curvature
    void collectEnoughRings(SurfaceMeshModel *mesh, Vertex v, const size_t min_nb, std::vector<int> & all);
    QVector<double> calCurvatureDesc(SurfaceMeshModel *mesh, RichParameterSet *pars);
    // PCA
    QVector<double> calPCADesc(SurfaceMeshModel *mesh);
    // CF
    double cotan(Vector3 v1, Vector3 v2);
    double arccos(Vector3 v1, Vector3 v2);

    void normalization(QVector<double>& vData, int nVrt);

    void calcDiscreteGaussianCurv(SurfaceMeshModel *mesh, QVector<double>& vGauCurv, double &fCurvSum);
    int	calcTargetUniformGaussianCurv(SurfaceMeshModel *mesh, QVector<double>& vTargetGauCurv, double fCurvSum);

    QVector<double> SolveSparseSystem(Eigen::SparseMatrix<double> &mat, QVector<double> &vDiffCurv);

    QVector<double> calConformalFactor(SurfaceMeshModel *mesh);

};
