#include "hacd-plugin.h"
#include "hacdlib.h"
#include "RenderObjectExt.h"

void hacd::initParameters(RichParameterSet *pars)
{
	pars->addParam(new RichBool("Visualize", true, "Visualize"));

	pars->addParam(new RichInt("Depth", 0, "Depth"));
	pars->addParam(new RichFloat("Concavity", 0.2f, "Concavity"));
	pars->addParam(new RichInt("MaxHullCount", 256, "MaxHullCount"));
	pars->addParam(new RichInt("MaxMergeHullCount", 256, "MaxMergeHullCount"));
	pars->addParam(new RichInt("MaxHullVertices", 64, "MaxHullVertices"));
	pars->addParam(new RichFloat("SmallClusterThreshold", 0.0f, "SmallClusterThreshold"));
	pars->addParam(new RichFloat("BackFaceDistanceFactor", 0.2f, "BackFaceDistanceFactor"));
	pars->addParam(new RichBool("NormalizeInputMesh", false, "NormalizeInputMesh"));
	pars->addParam(new RichBool("RemoveDuplicateVertices", true, "RemoveDuplicateVertices"));
	pars->addParam(new RichBool("UseFastVersion", false, "UseFastVersion"));
}

void hacd::applyFilter(RichParameterSet *pars)
{
	// Clear previous hulls
	for(auto m : document()->models()){
		if(m == mesh()) continue;
		document()->deleteModel(m);
	}

    std::vector<SurfaceMesh::SurfaceMeshModel*> pieces = HACDlib::decompose( mesh(), 
		pars->getInt("Depth"), pars->getFloat("Concavity"), pars->getInt("MaxHullCount"), 
		pars->getInt("MaxMergeHullCount"), pars->getInt("MaxHullVertices"),
		pars->getFloat("SmallClusterThreshold"), pars->getFloat("BackFaceDistanceFactor"), 
		pars->getBool("NormalizeInputMesh"), pars->getBool("RemoveDuplicateVertices"), pars->getBool("UseFastVersion"));

    for(auto m : pieces)
	{
		m->color = starlab::qRandomColor2();
        document()->addModel(m);
	}
}
