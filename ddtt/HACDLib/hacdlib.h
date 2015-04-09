#pragma once

#include "SurfaceMeshModel.h"

namespace HACDlib{
    std::vector<SurfaceMesh::SurfaceMeshModel*> decompose(SurfaceMesh::SurfaceMeshModel* fromMesh,
                                                          int mDecompositionDepth = 0,
                                                          float mConcavity = 0.2f,
                                                          float mMaxHullCount = 256,
                                                          float mMaxMergeHullCount = 256,
                                                          float mMaxHullVertices = 64,
                                                          float mSmallClusterThreshold = 0.0f,
                                                          float mBackFaceDistanceFactor = 0.2f,
                                                          bool mNormalizeInputMesh = false,
                                                          bool mRemoveDuplicateVertices = true,
                                                          bool mUseFastVersion = false);
}
