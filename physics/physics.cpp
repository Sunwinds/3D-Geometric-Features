#include "physics.h"
#include "PhysicsHelper.h"

#include "RenderObjectExt.h"

void physics::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichBool("Visualize", true, "Visualize"));
}

void physics::applyFilter(RichParameterSet *pars)
{
	Eigen::Matrix3d J;

	Vector3 com = PhysicsHelper(mesh()).centerOfMass(J, 1.0);

	if(pars->getBool("Visualize"))
	{
		starlab::SphereSoup * ss = new starlab::SphereSoup;
		ss->addSphere(com, mesh()->bbox().sizes().minCoeff() * 0.05);
		drawArea()->addRenderObject(ss);
		drawArea()->drawSegment(com, com + Vector3(0,0,-1));
	}
}
