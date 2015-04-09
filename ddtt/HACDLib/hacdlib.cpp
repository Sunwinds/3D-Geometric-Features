#include "hacdlib.h"

#pragma warning(disable:4996 4100)
#include "PlatformConfigHACD.h"
#include "HACD.h"
#include "JobSwarm.h"

std::vector<SurfaceMesh::SurfaceMeshModel *> HACDlib::decompose(SurfaceMesh::SurfaceMeshModel *fromMesh,
                                                                int mDecompositionDepth, float mConcavity,
                                                                float mMaxHullCount, float mMaxMergeHullCount,
                                                                float mMaxHullVertices, float mSmallClusterThreshold,
                                                                float mBackFaceDistanceFactor, bool mNormalizeInputMesh,
                                                                bool mRemoveDuplicateVertices, bool mUseFastVersion)
{
    HACD::HACD_API * gHACD = HACD::createHACD_API();
    HACD::HACD_API::Desc desc;

    /// Parameters
    {
        desc.mMaxHullVertices = mMaxHullVertices;
        desc.mMaxHullCount = mMaxHullCount;
        desc.mMaxMergeHullCount = mMaxMergeHullCount;
        desc.mSmallClusterThreshold = mSmallClusterThreshold;
        desc.mConcavity = mConcavity;
        desc.mBackFaceDistanceFactor = mBackFaceDistanceFactor;
        desc.mDecompositionDepth = mDecompositionDepth;
        desc.mNormalizeInputMesh = mNormalizeInputMesh;
        desc.mRemoveDuplicateVertices = mRemoveDuplicateVertices;
        desc.mUseFastVersion = mUseFastVersion;
    }

    desc.mTriangleCount = fromMesh->n_faces();
    desc.mVertexCount = fromMesh->n_vertices();

    std::vector<hacd::HaF32> vertices;
    std::vector<hacd::HaU32> faces;
    SurfaceMesh::Vector3VertexProperty points = fromMesh->vertex_coordinates();

    // Fill vertices
    for(auto v : fromMesh->vertices())
        for(int i = 0; i < 3; i++)
            vertices.push_back(points[v][i]);

    // Fill face info
    for(auto f : fromMesh->faces())
        for(auto v : fromMesh->vertices(f))
            faces.push_back(v.idx());

    desc.mVertices = &(vertices[0]);
    desc.mIndices = &(faces[0]);

    printf("Performing HACD on %d input triangles.\r\n", desc.mTriangleCount );
    printf("\r\n");
    printf("TriangleCount            : %d\r\n", desc.mTriangleCount );
    printf("VertexCount              : %d\r\n", desc.mVertexCount );
    printf("Concavity                : %0.2f\r\n", desc.mConcavity );
    printf("Max Hull Vertex Count    : %3d\r\n", desc.mMaxHullVertices );
    printf("Max Convex Hulls         : %3d\r\n", desc.mMaxHullCount );
    printf("Max Merged Convex Hulls  : %3d\r\n", desc.mMaxMergeHullCount );
    printf("Merge Threshold          : %0.2f\r\n", desc.mSmallClusterThreshold);
    printf("Back Face Distance Factor: %0.2f\r\n", desc.mBackFaceDistanceFactor );
    printf("Small Cluster Threshold  : %0.2f\r\n", desc.mSmallClusterThreshold );
    printf("DecompositionDepth       : %d\r\n", desc.mDecompositionDepth );

    std::vector<SurfaceMesh::SurfaceMeshModel *> hulls;

    hacd::HaU32 hullCount = gHACD->performHACD(desc);

    for (hacd::HaU32 i=0; i<hullCount; i++)
    {
        const HACD::HACD_API::Hull *hull = gHACD->getHull(i);
        if ( !hull ) continue;

        QString hullName = QString("hull_%1").arg(i);
        SurfaceMesh::SurfaceMeshModel * mesh = new SurfaceMesh::SurfaceMeshModel(hullName + ".obj", hullName);

        for (hacd::HaU32 i=0; i<hull->mVertexCount; i++)
        {
            const hacd::HaF32 *p = &hull->mVertices[i*3];
            mesh->add_vertex(SurfaceMesh::Vector3(p[0], p[1], p[2]));
        }

        for (hacd::HaU32 j=0; j<hull->mTriangleCount; j++)
        {
            hacd::HaU32 i1 = hull->mIndices[j*3+0];
            hacd::HaU32 i2 = hull->mIndices[j*3+1];
            hacd::HaU32 i3 = hull->mIndices[j*3+2];

            mesh->add_triangle(SurfaceMesh::Vertex(i1), SurfaceMesh::Vertex(i2), SurfaceMesh::Vertex(i3));
        }

        mesh->updateBoundingBox();
        mesh->update_face_normals();
        mesh->update_vertex_normals();

		hulls.push_back( mesh );
    }

    return hulls;
}
