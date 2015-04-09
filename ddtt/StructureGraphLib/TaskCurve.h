#pragma once
#include "Task.h"

class TaskCurve : public Task
{
    Q_OBJECT
public:
    TaskCurve( Structure::Graph * activeGraph, Structure::Graph * targetGraph, TaskType taskType, int ID ) :
        Task(activeGraph, targetGraph, taskType, ID){}

    void prepareCurve();
    void executeCurve(double t);

    /// Prepare Task
    void prepareGrowCurve();
    void prepareShrinkCurve();
    void prepareMorphCurve();

    /// Prepare sub-routines
    void prepareShrinkCurveOneEdge( Structure::Link* link );
	void prepareGrowCurveOneEdge( Structure::Link * tlink );
	void prepareCrossingMorphCurve();
	void encodeCurve(const Vector4d& coordinateA, const Vector4d& coordinateB);

    /// Execute
    void executeCrossingCurve( double t );
    void executeMorphCurve( double t );
    void foldCurve( double t );

    // Quick access
    Structure::Curve * targetCurve();

};
