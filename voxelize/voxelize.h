#pragma once

#include "SurfaceMeshPlugins.h"
#include "RichParameterSet.h"

class voxelize: public SurfaceMeshFilterPlugin{
    Q_OBJECT
    Q_PLUGIN_METADATA(IID "voxelize.plugin.starlab")
    Q_INTERFACES(FilterPlugin)

public:
    QString name() { return "Voxelize mesh"; }
    QString description() { return "Voxelize mesh"; }

    void initParameters(RichParameterSet* pars);
    void applyFilter(RichParameterSet* pars);
};

// Helper functions:
inline double ClosestPointTriangle(Vector3 p, Vector3 a, Vector3 b, Vector3 c, Vector3 & closest)
{
	// Check if P in vertex region outside A
	Vec3d ab = b - a;
	Vec3d ac = c - a;
	Vec3d ap = p - a;
	double d1 = dot(ab, ap);
	double d2 = dot(ac, ap);
	if (d1 <= 0 && d2 <= 0)
	{
		closest = a;
		return (p-closest).squaredNorm(); // barycentric coordinates (1,0,0)
	}
	// Check if P in vertex region outside B
	Vec3d bp = p - b;
	double d3 = dot(ab, bp);
	double d4 = dot(ac, bp);
	if (d3 >= 0 && d4 <= d3)
	{
		closest = b;
		return (p-closest).squaredNorm(); // barycentric coordinates (0,1,0)
	}
	// Check if P in edge region of AB, if so return projection of P onto AB
	double vc = d1*d4 - d3*d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) {
		double v = d1 / (d1 - d3);
		closest = a + v * ab;
		return (p - closest).squaredNorm(); // barycentric coordinates (1-v,v,0)
	}
	// Check if P in vertex region outside C
	Vec3d cp = p - c;
	double d5 = dot(ab, cp);
	double d6 = dot(ac, cp);
	if (d6 >= 0 && d5 <= d6)
	{
		closest = c;
		return (p-closest).squaredNorm(); // barycentric coordinates (0,0,1)
	}
	// Check if P in edge region of AC, if so return projection of P onto AC
	double vb = d5*d2 - d1*d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) {
		double w = d2 / (d2 - d6);
		closest = a + w * ac;
		return (p - closest).squaredNorm(); // barycentric coordinates (1-w,0,w)
	}
	// Check if P in edge region of BC, if so return projection of P onto BC
	double va = d3*d6 - d5*d4;
	if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
		double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		closest = b + w * (c - b);
		return (p - closest).squaredNorm(); // barycentric coordinates (0,1-w,w)
	}
	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	double denom = 1.0 / (va + vb + vc);
	double v = vb * denom;
	double w = vc * denom;
	closest = a + ab * v + ac * w;
	return (p - closest).squaredNorm(); // = u*a + v*b + w*c, u = va * denom = 1 - v - w
}
