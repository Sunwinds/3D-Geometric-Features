#include "ConformalFactor.h"
#include "CalConformalFactor.h"

#include "SurfaceMeshHelper.h"

#include <QDebug>


void ConformalFactor::initParameters(RichParameterSet *pars)
{

}

void ConformalFactor::applyFilter(RichParameterSet *pars)
{
    // Clear any past visualizations
    drawArea()->deleteAllRenderObjects();

    FilterPlugin * gd_plugin = pluginManager()->getFilter("Conformal Factor");

    drawArea()->updateGL();
    qApp->processEvents();

    Vector3VertexProperty points = mesh()->vertex_coordinates();
    int i = 0;
    CConformalFactor calCFactor;
    QVector<double> cfactor = calCFactor.calConformalFactor(mesh());
    qDebug() << "Conformal Factor: " << cfactor.size();

    double minVal = DBL_MAX;
    double maxVal = -DBL_MAX;

    foreach(double cf, cfactor)
    {
        minVal = qMin( minVal, cf );
        maxVal = qMax( maxVal, cf );
    }

    foreach(double cf, cfactor)
    {
        cf = (cf - minVal) / (maxVal - minVal);
    }

    foreach(Vertex v, mesh()->vertices())
    {

        drawArea()->drawPoint( points[v], 8.0, qtJetColorMap(1.0 - cfactor[i]) );
        i++;
        if(i == cfactor.size())
            break;
    }
}
