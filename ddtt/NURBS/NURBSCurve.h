// Geometric Tools, LLC
// Copyright (c) 1998-2012
// Distributed under the Boost Software License, Version 1.0.

#pragma once

#include "NURBSGlobal.h"
#include "SingleCurve.h"
#include "BSplineBasis.h"

namespace NURBS
{

template <typename Real>
class NURBSCurve : public SingleCurve<Real>
{
public:
    // Construction and destruction.  The caller is responsible for deleting
    // the input arrays if they were dynamically allocated.  Internal copies
    // of the arrays are made, so to dynamically change control points,
    // control weights, or knots, you must use the 'SetControlPoint',
    // 'GetControlPoint', 'SetControlWeight', 'GetControlWeight', and 'Knot'
    // member functions.

    // The homogeneous input points are (x,y,z,w) where the (x,y,z) values are
    // stored in the ctrlPoint array and the w values are stored in the
    // ctrlWeight array.  The output points from curve evaluations are of
    // the form (x',y',z') = (x/w,y/w,z/w).

    // Uniform spline.  The number of control points is n+1 >= 2.  The degree
    // of the spline is d and must satisfy 1 <= d <= n.  The knots are
    // implicitly calculated in [0,1].  If open is 'true', the spline is
    // open and the knots are
    //   t[i] = 0,               0 <= i <= d
    //          (i-d)/(n+1-d),   d+1 <= i <= n
    //          1,               n+1 <= i <= n+d+1
    // If open is 'false', the spline is periodic and the knots are
    //   t[i] = (i-d)/(n+1-d),   0 <= i <= n+d+1
    // If loop is 'true', extra control points are added to generate a closed
    // curve.  For an open spline, the control point array is reallocated and
    // one extra control point is added, set to the first control point
    // C[n+1] = C[0].  For a periodic spline, the control point array is
    // reallocated and the first d points are replicated.  In either case the
    // knot array is calculated accordingly.

    NURBSCurve (){}
	virtual ~NURBSCurve() {} 

    NURBSCurve (const Array1D_Vector3 & ctrlPoint, const Array1D_Real & ctrlWeight, int degree = 3, bool loop = false, bool open = true);

    static Array1D_Real uniformWeights(size_t n){ return Array1D_Real(n,1.0); }
	static NURBSCurve<Real> createCurveFromPoints( const Array1D_Vector3 & ctrlPoint ) 
	{ return NURBSCurve<Real>( ctrlPoint, uniformWeights( ctrlPoint.size() ) ); }

    static NURBSCurve<Real> createCurve(Vector3 from, Vector3 to, int steps = 5);

    // Open, nonuniform spline.  The knot array must have n-d elements.  The
    // elements must be nondecreasing.  Each element must be in [0,1].
    NURBSCurve (int numCtrlPoints, const Array1D_Vector3 & ctrlPoint, const Array1D_Real & ctrlWeight, int degree, bool loop, const Real* knot);

    size_t GetNumCtrlPoints() const;
    int GetDegree () const;
    bool IsOpen () const;
    bool IsUniform () const;
    bool IsLoop () const;

    // Control points and weights may be changed at any time.  The input index
    // should be valid (0 <= i <= n).  If it is invalid, the return value of
    // GetControlPoint is a vector whose components are all MAX_REAL, and the
    // return value of GetControlWeight is MAX_REAL.
    // undefined.
    void SetControlPoint (int i, const Vector3& ctrl);
    Vector3 GetControlPoint (int i) const;
    void SetControlWeight (int i, Real weight);
    Real GetControlWeight (int i) const;

    // The knot values can be changed only if the basis function is nonuniform
    // and the input index is valid (0 <= i <= n-d-1).  If these conditions
    // are not satisfied, GetKnot returns MAX_REAL.
    void SetKnot (int i, Real knot);
    Real GetKnot (int i) const;

	Array1D_Real GetKnotVector(bool isInnerOnly = false);

	// Refinement:
	void refine( Array1D_Real & insKnts, Array1D_Vector3 & Qw, Array1D_Real & Ubar );
    int findSpan(int n, int p, Real u, Array1D_Real U );
	int findSpan( Real u );
	
	Array1D_Vector3 midPointRefined();
	void computeMidPointRefine( Array1D_Vector3 & Qw, Array1D_Real & Ubar );

	Array1D_Real insertKnot(Real u, int k, int s, int r, Array1D_Vector3 & Qw, int * uq = NULL);

	Array1D_Vector3 simpleRefine(int k); 

	Array1D_Vector3 removeKnots( int iterations = 1 );
	Real KnotRemovalError( int r, int s );

    // The spline is defined for 0 <= t <= 1.  If a t-value is outside [0,1],
    // an open spline clamps t to [0,1].  That is, if t > 1, t is set to 1;
    // if t < 0, t is set to 0.  A periodic spline wraps to to [0,1].  That
    // is, if t is outside [0,1], then t is set to t-floor(t).
    virtual Vector3 GetPosition (Real t) ;
    virtual Vector3 GetFirstDerivative (Real t) ;
    virtual Vector3 GetSecondDerivative (Real t) ;
    virtual Vector3 GetThirdDerivative (Real t) ;

    // If you need position and derivatives at the same time, it is more
    // efficient to call these functions.  Pass the addresses of those
    // quantities whose values you want.  You may pass 0 in any argument
    // whose value you do not want.
    void Get (Real t, Vector3* pos, Vector3* der1,
        Vector3* der2, Vector3* der3);

    // Access the basis function to compute it without control points.  This
    // is useful for least squares fitting of curves.
    BSplineBasis<Real>& GetBasis ();

    std::vector<Vector3> getControlPoints();

    Real timeAt(const Vector3 &pos);
	Real fastTimeAt(const Vector3 &pos);

    std::vector<std::vector<Vector3> > toSegments(Scalar resolution);

    // Modify
    void translate(const Vector3 &delta);
    void scale(Scalar scaleFactor);
    void scaleInPlace(Scalar scaleFactor, int placeCtrlPoint);
    void translateTo(const Vector3 &newPos, int cpIDX);

public:
    // Replicate the necessary number of control points when the Create
    // function has loop equal to true, in which case the spline curve must
    // be a closed curve.
    void CreateControl (const Array1D_Vector3 & ctrlPoint, const Array1D_Real & ctrlWeight);

    size_t mNumCtrlPoints;
    Array1D_Vector3 mCtrlPoint;  // ctrl[n+1]
    Array1D_Real mCtrlWeight;    // weight[n+1]
    bool mLoop;
    BSplineBasis<Real> mBasis;
    int mReplicate;  // the number of replicated control points

	// Misc
	Array1D_Vector3 misc_points;
};

typedef NURBSCurve<float> NURBSCurvef;
typedef NURBSCurve<double> NURBSCurved;

}
