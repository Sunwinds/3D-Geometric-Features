#pragma once
#include "SurfaceMeshModel.h"

#include <stack>

enum SymmetryType { TRANSLATION, REFLECTION, ROTATION, INVALID };
struct SymmetryGroup {
    SymmetryType type;
    Eigen::Vector3d axis, delta;
    std::vector<SurfaceMesh::SurfaceMeshModel *> members;
    SymmetryGroup() : type(INVALID), axis(Eigen::Vector3d(0,0,0)), delta(Eigen::Vector3d(0,0,0)) {}
};

class SymmetryAnalysis
{
public:
    SymmetryAnalysis( SurfaceMesh::SurfaceMeshModel * mesh );
    SymmetryAnalysis( std::vector<SurfaceMesh::SurfaceMeshModel *> meshes );

    void analyze();
	bool hasMesh(QString mesh_name);

	Eigen::MatrixXd centers;
	Eigen::AlignedBox3d boundingVolume;
	Eigen::Vector3d axis, primaryAxis, secondaryAxis, origin;
	bool isFitLine, isRotational, isDiagonalTranslation, isReflecitonal;
	
	static std::vector< SurfaceMesh::SurfaceMeshModel* > connectedPieces(SurfaceMesh::SurfaceMeshModel * mesh);
	std::vector<SurfaceMesh::SurfaceMeshModel *> meshes;

	QString description();
};

std::vector< SurfaceMesh::SurfaceMeshModel* > SymmetryAnalysis::connectedPieces(SurfaceMesh::SurfaceMeshModel * mesh) {
    std::vector< SurfaceMesh::SurfaceMeshModel* > pieces;
    std::vector< std::vector< SurfaceMesh::Face > > pieceFaces;

    SurfaceMesh::Vector3VertexProperty points = mesh->vertex_coordinates();
    SurfaceMesh::BoolFaceProperty fvisisted = mesh->face_property<bool>("f:visisted", false);
    SurfaceMesh::ScalarFaceProperty fsegment = mesh->face_property<SurfaceMesh::Scalar>("f:segment", -1);
    for(auto f : mesh->faces()){
        fvisisted[f] = false;
        fsegment[f] = -1;
    }

    // Track pieces vertex count
    std::vector<int> sizes;

    // Visit all faces
    for( SurfaceMesh::Face f : mesh->faces() )
    {
        if( fvisisted[f] ) continue;

        pieceFaces.push_back( std::vector< SurfaceMesh::Face >() );

        std::stack<SurfaceMesh::Face> unprocessed;
        unprocessed.push(f);

        while( !unprocessed.empty() )
        {
            SurfaceMesh::Face curFace = unprocessed.top();
            unprocessed.pop();
            if(fvisisted[curFace]) continue;

            fvisisted[curFace] = true;

            for(auto v : mesh->vertices(curFace))
            {
                for(auto fj : mesh->faces(v))
                {
                    if(!mesh->is_valid(fj) || fvisisted[fj]) continue;
                    unprocessed.push(fj);
                }
            }

            fsegment[curFace] = (pieceFaces.size() - 1);

            pieceFaces.back().push_back( curFace );
        }

        SurfaceMesh::SurfaceMeshModel * piece = new SurfaceMesh::SurfaceMeshModel;
        piece->reserve( uint(pieceFaces.back().size() * 3), uint(pieceFaces.back().size() * 6), uint(pieceFaces.back().size()) );

        for(auto v : mesh->vertices())
        {
            piece->add_vertex( points[v] );
        }

        for(auto f : pieceFaces.back()){
            std::vector<SurfaceMesh::Vertex> face;
            for(auto v : mesh->vertices(f)) face.push_back(v);
            piece->add_face(face);
        }

        for(auto v : piece->vertices()) if(piece->is_isolated(v)) piece->remove_vertex(v);
        piece->garbage_collection();

        sizes.push_back( piece->n_vertices() );
        pieces.push_back( piece );
    }

    return pieces;
}

SymmetryAnalysis::SymmetryAnalysis(SurfaceMesh::SurfaceMeshModel *mesh)
{
    meshes = connectedPieces( mesh );
    analyze();
}

SymmetryAnalysis::SymmetryAnalysis(std::vector<SurfaceMesh::SurfaceMeshModel *> meshes) : meshes( meshes )
{
    analyze();
}

void SymmetryAnalysis::analyze()
{
	if( meshes.empty() ) return;

	bool dbg = false;

	int N = int( meshes.size() );

	// Compute bounding volume and keep track of the boxes
	centers = Eigen::MatrixXd ( N, 3 );
	std::vector<Eigen::AlignedBox3d> mboxes( N );

	for(int i = 0; i < N; i++) 
	{
		auto m = meshes[i];
		m->updateBoundingBox();
		Eigen::AlignedBox3d mbbox = m->bbox();
		boundingVolume = boundingVolume.extend( mbbox );
		centers.row(i) = mbbox.center();
		mboxes[i] = mbbox;
	}

	// Check for concentric parts
	{
		bool isHasConcentric = false;

		// Tolerance to fitting centers to a line
		double tolerance = boundingVolume.sizes().minCoeff() * 0.2;

		for(int i = 0; i < N; i++){
			for(int j = i+1; j < N; j++){
				double dist = (centers.row(i) - centers.row(j)).norm();
				if( dist < tolerance ){
					isHasConcentric = true;
					i = j = N; // break all loops
				}
			}
		}

		// Change centers to top most point
		if( isHasConcentric ){
			for(int i = 0; i < N; i++){
				SurfaceMesh::Vector3VertexProperty points = meshes[i]->vertex_coordinates();
				for(auto v : meshes[i]->vertices())
					if(points[v].z() > centers.row(i).z()) 
						centers.row(i) = points[v];
			}
		}
	}

	origin = centers.colwise().mean();

	// Check if all centers fall within a straight line
	isFitLine = true;
	{
		// Tolerance to fitting centers to a line
		double tolerance = boundingVolume.sizes().minCoeff() * 0.1;
	
		// Compute PCA
		Eigen::MatrixXd centered = centers.rowwise() - origin.transpose();
		Eigen::MatrixXd cov = centered.adjoint() * centered;
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

		// Compute distances to fitting line
		axis = eig.eigenvectors().col(2).normalized();
		primaryAxis = axis;
		secondaryAxis = eig.eigenvectors().col(1).normalized();

		for(int i = 0; i < N; i++) 
		{
			double dist = axis.cross( Eigen::Vector3d(centers.row(i) - origin.transpose()) ).norm();
			if(dist > tolerance){
				isFitLine = false;
				if(dbg)qDebug() << "Line fit test: large deviation " << dist << " > " << tolerance;
				axis = eig.eigenvectors().col(0).normalized();
				break;
			}
		}
	}

	// Check for rotational symmetries
	isRotational = false;
	bool isNonFourGraterTwo = N != 4 && N > 2;
	if( !isFitLine && isNonFourGraterTwo )
	{
		isRotational = true;

		Eigen::Vector3d v0 = centers.row(0) - origin.transpose();
		double v0_dist = v0.norm();
		double tolerance = v0_dist * 0.05;

		for(int i = 1; i < N; i++) 
		{
			double dist = abs( (centers.row(i) - origin.transpose()).norm() - v0_dist );
			if(dist > tolerance){
				isRotational = false;
				if(dbg)qDebug() << "Rotational test: large deviation " << dist << " > " << tolerance;
				break;
			}
		}
	}

	// Check if axis is not globally aligned
	isDiagonalTranslation = false;
	if( isFitLine && !isRotational )
	{
		double tol = Eigen::Vector3d(1,1,1).normalized().z() * 0.5;
		isDiagonalTranslation = abs(axis.x()) > tol && abs(axis.y()) > tol && abs(axis.z()) > tol;
		
		if( isDiagonalTranslation )
			if(dbg)qDebug() << "Diagonal translation test: axis is not globally aligned > " << tol;
	}

	isReflecitonal = false;
	if( !isRotational && !isDiagonalTranslation )
	{
		isReflecitonal = true;
	}
}

QString SymmetryAnalysis::description()
{
	QStringList desc;

	desc << QString("LINE [%1]").arg( isFitLine );
	desc << QString("ROTA [%1]").arg( isRotational );
	desc << QString("REFL [%1]").arg( isReflecitonal );
	desc << QString("DIAG [%1]").arg( isDiagonalTranslation );

	if( isFitLine ) desc << QString("LINE AXIS (%1, %2, %3)").arg(axis.x()).arg(axis.y()).arg(axis.z());
	if( isRotational ) desc << QString("ROT AXIS (%1, %2, %3)").arg(axis.x()).arg(axis.y()).arg(axis.z());
	if( isReflecitonal ) {
		desc << QString("PRIME AXIS (%1, %2, %3)").arg(primaryAxis.x()).arg(primaryAxis.y()).arg(primaryAxis.z());
		desc << QString("SECON AXIS (%1, %2, %3)").arg(secondaryAxis.x()).arg(secondaryAxis.y()).arg(secondaryAxis.z());
	}
	
	return desc.join("\n");
}

bool SymmetryAnalysis::hasMesh(QString mesh_name)
{
	for(auto m : meshes) if(m->name == mesh_name) return true;
	return false;
}
