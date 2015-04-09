#include "CalConformalFactor.h"


QVector<double> CConformalFactor::calConformalFactor(SurfaceMeshModel *mesh)
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
    Eigen::SparseMatrix<double> matLaplacian = Laplacian(mesh);
    QVector<double> Phi = SolveSparseSystem(matLaplacian, vDiffCurv);

    // normalization
    normalization(Phi, nVrt);

    qDebug() << "Done" ;

    return Phi;
}

void CConformalFactor::calcDiscreteGaussianCurv(SurfaceMeshModel *mesh, QVector<double>& vGauCurv, double &fCurvSum)
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
    qDebug() << "Done" ;
}

int CConformalFactor::calcTargetUniformGaussianCurv(SurfaceMeshModel *mesh, QVector<double>& vTargetGauCurv, double fCurvSum)
{
    SurfaceMeshHelper h(mesh);
    Scalar area = 0.0;
    ScalarFaceProperty  fArea = h.computeFaceAreas();

    foreach(Face f, mesh->faces()){
        area += fArea[f];
    }
    double fCurv = fCurvSum / area;
    qDebug() << "Area = " << area;

    //Vector3VertexProperty points = mesh->vertex_coordinates();

    //Tri::CVertex* pVrt = mesh->vertices();
    //Tri::CTriangle* pTri = mesh->m_pTriangles;

    //int nVrt = mesh->n_vertices();
    //for (int iv=0; iv<nVrt; iv++)
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
    qDebug() << "Done" ;
    return 1;
}
