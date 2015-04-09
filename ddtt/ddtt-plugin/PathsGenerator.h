#pragma once
#include "DeformationPath.h"

class PathsGenerator
{
public:
    PathsGenerator(Structure::Graph *source, Structure::Graph *target, int K);
    std::vector<DeformationPath> paths;

    Assignments allAssignments(QVector<QVector<QString> > sgroups,
                               QVector<QVector<QString> > tgroups,
                               std::vector< std::vector<double> > m, Structure::Graph *sg, Structure::Graph *tg, int K);

    QVector<VectorPairings> manyToMany(Structure::Graph *sg, Structure::Graph *tg,
										QVector<QString> snodes, QVector<QString> tnodes);

	void prepareSymmetryAnalysis(Structure::Graph * sg, Structure::Graph * tg);
	QVector<VectorPairings> manyToMany2(Structure::Graph *sg, Structure::Graph *tg,
										QVector<QString> snodes, QVector<QString> tnodes);

    VectorPairings splitMany(QVector<QString> A, QVector<QString> B);

    QVector<QString> nothingSet();
    VectorPairStrings pairsDebugging(QVector<Pairing> pairs);

    double partDifference2(QString sid, QString tid, Structure::Graph *source, Structure::Graph *target);
    double partDifference(QString sid, QString tid, Structure::Graph *source, Structure::Graph *target);

    void compute_part_measures(Structure::Graph *graph);
    void measure_position(Structure::Graph *graph);
    void measure_ground_parts(Structure::Graph *graph);

    std::vector< std::vector<double> > buildDifferenceMatrix(Structure::Graph *source, Structure::Graph *target);
    void visualizeDifferenceMatrix(std::vector< std::vector<double> > m, Structure::Graph *sg, Structure::Graph *tg);

};
