#include "symmetry.h"
#include "RenderObjectExt.h"

#include "SymmetryAnalysis.h"

void symmetry::initParameters(RichParameterSet *pars)
{
    pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

void symmetry::applyFilter(RichParameterSet *pars)
{
	SymmetryAnalysis sanalysis( mesh() );

	mainWindow()->setStatusBarMessage( sanalysis.description() );

    bool isVisualize = pars->getBool("Visualize");
    if( isVisualize )
    {
		if( sanalysis.isRotational )
		{
			starlab::SphereSoup * centerSphere = new starlab::SphereSoup;
			centerSphere->addSphere( starlab::QVector3(sanalysis.boundingVolume.center()), 0.1f );
			drawArea()->addRenderObject(centerSphere);
			drawArea()->drawSegment(sanalysis.boundingVolume.center(), sanalysis.boundingVolume.center() + sanalysis.axis);
		}

		if( sanalysis.isFitLine )
		{
			drawArea()->drawSegment( sanalysis.origin - sanalysis.axis, sanalysis.origin + sanalysis.axis, 10 );
		}

		if( sanalysis.isReflecitonal && !sanalysis.isFitLine )
		{
			drawArea()->drawSegment( sanalysis.origin - sanalysis.primaryAxis, sanalysis.origin + sanalysis.primaryAxis, 10, Qt::green );
			drawArea()->drawSegment( sanalysis.origin - sanalysis.secondaryAxis, sanalysis.origin + sanalysis.secondaryAxis, 10, Qt::blue );
		}
    }

    foreach(auto m, sanalysis.meshes)
		document()->addModel(m);
}
