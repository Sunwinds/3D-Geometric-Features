#pragma once

#include "SurfaceMeshPlugins.h"
#include "SurfaceMeshModel.h"

#include <QVector>

#include <cstdio>

template <typename T>
void read(FILE *in, T &t)
{
    int err = 0;
    err = fread(&t, 1, sizeof(t), in);
}

class SegMeshLoader
{
public:
    SegMeshLoader(){
        this->entireMesh = NULL;
    }
    ~SegMeshLoader(){
        if(entireMesh)
            delete entireMesh;
    }

    // Entire mesh
    SurfaceMeshModel* entireMesh;
    Vector3VertexProperty entirePoints;

    // Groups load from .seg file
    QMap< QString, QVector<int> > groupFaces;

    QVector<SurfaceMeshModel*> loadGroupMeshFromOFF(QString meshF, QString segF);
    SurfaceMeshModel* extractSegMesh( QString gid );

	SurfaceMeshModel* readOBJ(QString filename);

private:
    // help function
    SurfaceMeshModel* readOFF(QString filename);
    bool readSeg(QString filename);

	bool saveOBJ(QString filename, SurfaceMeshModel *mesh);
    bool saveOFF(QString filename, SurfaceMeshModel *mesh);
    void normalize(SurfaceMeshModel *mesh);
};
