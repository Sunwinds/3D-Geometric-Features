#include "PCAShapeEst.h"

#include <QDebug>

void PCAShapeEst::initParameters(RichParameterSet *pars)
{

}

void PCAShapeEst::applyFilter(RichParameterSet *pars)
{
    Vector_3 scalar = computePCAScalar(mesh()->vertices_begin(), mesh()->vertices_end(), mesh()->vertex_property<Point>(VPOINT));

    QString text = QString::number(scalar[0]) + ", " + QString::number(scalar[1]) + ", " + QString::number(scalar[2]);

   // mainWindow()->setStatusBarMessage("PCAScalar: " + text);
}

Vector_3 PCAShapeEst::computePCAScalar(VertexIterator begin, VertexIterator end, Points pts)
{
    int n = mesh()->vertices_size();

    qDebug() << "Vertices_Size = "<< n;

    FT x, y, z,
        sumX = 0., sumY = 0., sumZ = 0.,
        sumX2 = 0., sumY2 = 0., sumZ2 = 0.,
        sumXY = 0., sumXZ = 0., sumYZ = 0.,
        xx, yy, zz, xy, xz, yz;

    for (; begin != end; begin++)
    {
        //Point_3 lp = Point;
        x = pts[begin][0];//lp.x();
        y = pts[begin][1];//lp.y();
        z = pts[begin][2];//lp.z();
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

    // assemble covariance matrix as a
    // semi-definite matrix.
    // Matrix numbering:
    // 0
    // 1 2
    // 3 4 5
    Eigen::Matrix3d covariance;
    covariance << xx, xy, xz,
        xy, yy, yz,
        xz, yz, zz;

    // solve for eigenvalues and eigenvectors.
    // eigen values are sorted in ascending order,
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(covariance);

    Eigen::Vector3d eigenVals = es.eigenvalues();

     qDebug()<< "EigenValues: " <<eigenVals[0]<<','<<eigenVals[1]<<','<<eigenVals[2];

    Vector_3 scalar;

    FT sum = eigenVals[0] + eigenVals[1] + eigenVals[2];
    scalar[0] = (eigenVals[2] - eigenVals[1]) / sum;
    scalar[1] = 2 * (eigenVals[1] - eigenVals[0]) / sum;
    scalar[2] = 3 * eigenVals[0]/ sum;

    qDebug()<< "PCAScalar: "<<scalar[0]<<','<<scalar[1]<<','<<scalar[2];

    return scalar;
}
