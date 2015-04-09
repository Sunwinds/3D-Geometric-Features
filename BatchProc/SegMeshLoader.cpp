#include "SegMeshLoader.h"

#include <QFile>
#include <QTextStream>

#include "Eigen/LU"
#include "Eigen/Geometry"
#include "Eigen/Eigenvalues"

#include "Decimater.h"

QVector<SurfaceMeshModel*> SegMeshLoader::loadGroupMeshFromOFF(QString meshF, QString segF)
{
    //qDebug() << "Start read OFF";
    entireMesh = readOFF(meshF);
	//qDebug() << "Start normalization";
	//normalize(entireMesh);
    entirePoints = entireMesh->vertex_coordinates();

    //qDebug() << "Start read Seg";
    readSeg(segF);

    QVector<SurfaceMeshModel*> meshes;
    foreach (QString gid, groupFaces.keys()){
		SurfaceMeshModel * submesh = extractSegMesh(gid);
        meshes.push_back(submesh);
		//saveOFF(meshF, submesh);
		saveOBJ(meshF, submesh);
	}

    //qDebug() << "#Seg" << meshes.size();
    return meshes;
}

SurfaceMeshModel *SegMeshLoader::extractSegMesh(QString gid)
{
    SurfaceMeshModel* subMesh = new SurfaceMeshModel(gid + ".obj", gid);

    QVector<int> part = groupFaces[gid];
    //qDebug() << "#Faces" << part.size();

    // vertex set from face
    QSet<int> vertSet;
    foreach(int fidx, part)
    {
        Surface_mesh::Vertex_around_face_circulator vit = entireMesh->vertices(Face(fidx)),vend=vit;
		/*for(auto v: entireMesh->vertices(Face(fidx))){
			vertSet.insert(Vertex(vit).idx());
		}*/
        do{ 
			int vidx = Vertex(vit).idx();
			vertSet.insert(vidx); 
		} while(++vit != vend);
    }

	//int i = 0;

    // entire vid => segment vid
    //QMap<Vertex,Vertex> vmap;
	QMap<int,Vertex> vmap;
    foreach(int vidx, vertSet)
    {
        //vmap[Vertex(vidx)] = Vertex(vmap.size());
		vmap[vidx] = Vertex(vmap.size());
        subMesh->add_vertex( entirePoints[Vertex(vidx)] );
		//qDebug() << i++ << vidx << entirePoints[Vertex(vidx)][0] << entirePoints[Vertex(vidx)][1] << entirePoints[Vertex(vidx)][2];
        //qDebug() << entirePoints[Vertex(vidx)][0] << entirePoints[Vertex(vidx)][1] << entirePoints[Vertex(vidx)][2];
    }

    // copy vertices to subMesh
    foreach(int fidx, part){
        std::vector<Vertex> pnts;
        Surface_mesh::Vertex_around_face_circulator vit = entireMesh->vertices(Face(fidx)),vend=vit;
        do{ 
			int vidx = Vertex(vit).idx();
			//int cid = vmap[vidx].idx();
			pnts.push_back(vmap[vidx]); 
		} while(++vit != vend);
        subMesh->add_face(pnts);
    }

    //qDebug() << subMesh->n_vertices() << subMesh->n_faces();

    subMesh->updateBoundingBox();
    subMesh->isVisible = false;

	// Simplification to save space
	int nVert = subMesh->n_vertices();
	//if( nVert > 6000){
	//	qDebug() << "Before simplification" << nVert;
	//	float percent = 6000.0f / float(nVert);
	//	Decimater::simplify(subMesh, percent);
	//	//subMesh->update_vertex_normals();
	//	qDebug() << "After simplification" << subMesh->n_vertices();
	//}
	
    return subMesh;
}

SurfaceMeshModel* SegMeshLoader::readOBJ(QString filename)
{
	SurfaceMeshModel * mesh = new SurfaceMeshModel(filename);
	bool has_loaded = mesh->read(qPrintable(filename));
	//qDebug() << "#Vert" << mesh->n_vertices() << "#Face" << mesh->n_faces();
	if(!has_loaded)
		delete mesh;
	//throw StarlabException("surfacemesh_io_off::open failed");
	normalize(mesh);
	return mesh;
}

SurfaceMeshModel* SegMeshLoader::readOFF(QString filename)
{
    SurfaceMeshModel * mesh = new SurfaceMeshModel(filename);
    bool has_loaded = mesh->read(qPrintable(filename));
    //qDebug() << "#Vert" << mesh->n_vertices() << "#Face" << mesh->n_faces();
    if(!has_loaded)
        delete mesh;
    //throw StarlabException("surfacemesh_io_off::open failed");
    return mesh;
}

bool SegMeshLoader::readSeg(QString filename)
{
    QFile file(filename);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return false;
    QTextStream in(&file);

    int fidx = 0;
    while(!in.atEnd()){
        QString gid = in.readLine();
        groupFaces[gid].push_back(fidx);
        fidx++;
    }

    file.close();
    return true;
}

void SegMeshLoader::normalize(SurfaceMeshModel *mesh)
{
    //qDebug() << "#Vert" << entireMesh->n_vertices() << "#Face" << entireMesh->n_faces();
    mesh->updateBoundingBox();
    Eigen::AlignedBox3d bbox = mesh->bbox();
    Vector3 offset = bbox.center();

    /// Normalize to have longest side size = 1
    Vector3 s = bbox.diagonal();
    Scalar scale = qMax(s.x(),qMax(s.y(),s.z()));

    Vector3VertexProperty points = mesh->vertex_coordinates();
    foreach(Vertex v, mesh->vertices()){
         Point& p = points[v];
         p.x() -= offset.x();
         p.y() -= offset.y();
         p.z() -= offset.z();
         p.x() /= scale;
         p.y() /= scale;
         p.z() /= scale;
    }

    /// And update it after, so we can reset the viewpoint
    mesh->updateBoundingBox();
}

bool SegMeshLoader::saveOBJ( QString filename,  SurfaceMeshModel *mesh)
{
	if(filename.isEmpty()) return false;
	filename.chop(4);
	QString filepath = filename + '_' + mesh->path;
	FILE* fid = fopen( filepath.toStdString().c_str(), "w" );
	//if( fid == NULL ) throw StarlabException("the file cannot be opened");

	Vector3VertexProperty p = mesh->vertex_coordinates();
	//Vector3VertexProperty n = mesh->vertex_normals();

	foreach(Vertex v, mesh->vertices()){
		//if(n) fprintf( fid, "vn %.10f %.10f %.10f\n", n[v].x(), n[v].y(), n[v].z() );        
		fprintf( fid, "v %.10f %.10f %.10f\n", p[v].x(), p[v].y(), p[v].z() );
	}

	foreach(Face f, mesh->faces()) {
		// int nV = mesh->valence(f);
		fprintf(fid, "f");
		foreach(Vertex v, mesh->vertices(f)){
			int ni = ((Surface_mesh::Vertex)v).idx()+1;
			//if(n) fprintf(fid, " %d//%d", ni, ni);
			fprintf(fid, " %d", ni);
		}
		fprintf(fid,"\n");
	}

	fclose(fid);
	return true;
}

bool SegMeshLoader::saveOFF( QString filename, SurfaceMeshModel *mesh )
{
	filename.chop(4);
	QString filepath = filename + '_' + mesh->name + ".off";
	return mesh->write(filepath.toStdString());
}
