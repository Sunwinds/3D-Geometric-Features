#include "spherelib.h"
#include <random>

#include <QPainter>
#include <QGLWidget>

int Spherelib::SphereMaker::addVertex(const Spherelib::SphereMaker::Vector3 &p)
{
    geometry.add_vertex( p.normalized() );
    return index++;
}

int Spherelib::SphereMaker::getMiddlePoint(int p1, int p2)
{
    // first check if we have it already
    bool firstIsSmaller = p1 < p2;
    size_t smallerIndex = firstIsSmaller ? p1 : p2;
    size_t greaterIndex = firstIsSmaller ? p2 : p1;
    size_t key = (smallerIndex << 32) + greaterIndex;

    if(middlePointIndexCache.find(key) != middlePointIndexCache.end())
        return middlePointIndexCache[key];

    // not in cache, calculate it
	Vector3 v1 = geometry.vertex_property<Vector3>("v:point")[Surface_mesh::Vertex(p1)];
	Vector3 v2 = geometry.vertex_property<Vector3>("v:point")[Surface_mesh::Vertex(p2)];
    Vector3 middle = ( v1 + v2 ) * 0.5;

    // add vertex makes sure point is on unit sphere
    int i = addVertex( middle );

    // store it, return index
    middlePointIndexCache[key] = i;
    return i;
}

Spherelib::SphereMaker::SphereMaker(int recursionLevel)
{
    index = 0;

    // create 12 vertices of a icosahedron
    auto t = (1.0 + std::sqrt(5.0)) / 2.0;

    addVertex(Vector3(-1,  t,  0));
    addVertex(Vector3( 1,  t,  0));
    addVertex(Vector3(-1, -t,  0));
    addVertex(Vector3( 1, -t,  0));

    addVertex(Vector3( 0, -1,  t));
    addVertex(Vector3( 0,  1,  t));
    addVertex(Vector3( 0, -1, -t));
    addVertex(Vector3( 0,  1, -t));

    addVertex(Vector3( t,  0, -1));
    addVertex(Vector3( t,  0,  1));
    addVertex(Vector3(-t,  0, -1));
    addVertex(Vector3(-t,  0,  1));

    // create 20 triangles of the icosahedron
    std::list<TriangleIndices> faces;

    // 5 faces around point 0
    faces.push_back(TriangleIndices(0, 11, 5));
    faces.push_back(TriangleIndices(0, 5, 1));
    faces.push_back(TriangleIndices(0, 1, 7));
    faces.push_back(TriangleIndices(0, 7, 10));
    faces.push_back(TriangleIndices(0, 10, 11));

    // 5 adjacent faces
    faces.push_back(TriangleIndices(1, 5, 9));
    faces.push_back(TriangleIndices(5, 11, 4));
    faces.push_back(TriangleIndices(11, 10, 2));
    faces.push_back(TriangleIndices(10, 7, 6));
    faces.push_back(TriangleIndices(7, 1, 8));

    // 5 faces around point 3
    faces.push_back(TriangleIndices(3, 9, 4));
    faces.push_back(TriangleIndices(3, 4, 2));
    faces.push_back(TriangleIndices(3, 2, 6));
    faces.push_back(TriangleIndices(3, 6, 8));
    faces.push_back(TriangleIndices(3, 8, 9));

    // 5 adjacent faces
    faces.push_back(TriangleIndices(4, 9, 5));
    faces.push_back(TriangleIndices(2, 4, 11));
    faces.push_back(TriangleIndices(6, 2, 10));
    faces.push_back(TriangleIndices(8, 6, 7));
    faces.push_back(TriangleIndices(9, 8, 1));

    // refine triangles
    for (int i = 0; i < recursionLevel; i++)
    {
        std::list<TriangleIndices> faces2;
        for (auto tri : faces)
        {
            // replace triangle by 4 triangles
            int a = getMiddlePoint(tri.v1, tri.v2);
            int b = getMiddlePoint(tri.v2, tri.v3);
            int c = getMiddlePoint(tri.v3, tri.v1);

            faces2.push_back(TriangleIndices(tri.v1, a, c));
            faces2.push_back(TriangleIndices(tri.v2, b, a));
            faces2.push_back(TriangleIndices(tri.v3, c, b));
            faces2.push_back(TriangleIndices(a, b, c));
        }
        faces = faces2;
    }

    // done, now add triangles to mesh
    for (auto tri : faces)
    {
        geometry.add_triangle( Surface_mesh::Vertex(tri.v1),
                               Surface_mesh::Vertex(tri.v2),
                               Surface_mesh::Vertex(tri.v3) );
    }
}

Spherelib::Sphere::Sphere( int resolution, Vector3 center, double radius ) : center(center), radius(radius)
{
	// Copy sphere geometry to our SurfaceMeshModel
	geometry = new SurfaceMesh::SurfaceMeshModel;
	Surface_mesh tempGeometry = SphereMaker::makeSphere( resolution );
    for(int i = 0; i < (int)tempGeometry.n_vertices(); i++) geometry->add_vertex( 
		tempGeometry.vertex_property<Vector3>("v:point")[Surface_mesh::Vertex(i)] );
	Surface_mesh::Face_iterator fit, fend = tempGeometry.faces_end();
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;
	for (fit=tempGeometry.faces_begin(); fit!=fend; ++fit){
		fvit = fvend = tempGeometry.vertices(fit);
		std::vector<Surface_mesh::Vertex> verts;
		do{ verts.push_back(fvit); } while (++fvit != fvend);
		geometry->add_face(verts);
	}

	// Save number of points
	numPoints = geometry->n_vertices();
}

std::vector<Spherelib::Sphere::Vector3> Spherelib::Sphere::rays() const
{
	std::vector<Spherelib::Sphere::Vector3> result;
	for(int v = 0; v < numPoints; v++) result.push_back(geometry->vertex_property<Vector3>("v:point")[Surface_mesh::Vertex(v)]);
	return result;
}

void Spherelib::Sphere::fillRandom()
{
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	std::random_device rand_dev;          // Use random_device to get a random seed.
	std::mt19937 rand_engine(rand_dev()); 

	Surface_mesh::Vertex_property<float> values = geometry->vertex_property<float>("v:values", 0.0);
	for(int i = 0; i < 10; i++)
	{
		Surface_mesh::Vertex rv( unif(rand_engine) * (numPoints-1) );
		values[rv] = unif(rand_engine);
	}
}

void Spherelib::Sphere::fillPattern()
{
	static int pattern = 0;

	Surface_mesh::Vertex_property<float> values = geometry->vertex_property<float>("v:values", 0.0);
	Surface_mesh::Vertex_property<Vector3> points = geometry->vertex_property<Vector3>("v:point");

	for( auto v : geometry->vertices() )
	{
		double x = points[v].x(), y = points[v].y(), z = points[v].z();
		double val = 0;

		switch (pattern){
			case 0: val = pow(x,2) + pow(y,2); break;
			case 1: val = pow(x,2) + pow(z,2); break;
			case 2: val = pow(y,2) + pow(z,2); break;
			case 3: val = cos(pow(x,2) + pow(z,2)); break;
			case 4: val = sin(pow(y,2) + pow(z,2)); break;
			default: break;
		}

		values[v] = val;
	}

	pattern++;
}

void Spherelib::Sphere::smoothValues( int iterations )
{
	Surface_mesh::Vertex_property<float> vprop = geometry->vertex_property<float>("v:values");

	for(int i = 0; i < iterations; i++)
	{
		std::vector<float> newValues(geometry->n_vertices(), 0.0);

		// average the values of neighbors
		for( auto v : geometry->vertices() ){
			for( auto vj : geometry->onering_hedges(v) )
				newValues[ v.idx() ] += vprop[ geometry->to_vertex(vj) ];
			newValues[ v.idx() ] /= geometry->valence(v);
		}

		// copy results back to property
		for(auto v : geometry->vertices()) vprop[v] = newValues[v.idx()];
	}
}

void Spherelib::Sphere::normalizeValues()
{
	Surface_mesh::Vertex_property<float> vprop = geometry->vertex_property<float>("v:values");
	std::pair<float,float> bounds( FLT_MAX, -FLT_MAX );
	for( auto v : geometry->vertices() ){
		bounds.first = std::min(bounds.first, vprop[v]);
		bounds.second = std::max(bounds.second, vprop[v]);
	}

	float range = bounds.second - bounds.first;
	for( auto v : geometry->vertices() ) vprop[v] = (vprop[v]-bounds.first) / range;
}

std::vector<float> Spherelib::Sphere::values() const
{
	std::vector<float> values;
	const Surface_mesh::Vertex_property<float> vprop = geometry->vertex_property<float>("v:values");
	for( auto v : geometry->vertices() ) values.push_back(vprop[v]);
	return values;
}

void Spherelib::Sphere::setValues(const std::vector<float>& newValues)
{
	Surface_mesh::Vertex_property<float> vprop = geometry->vertex_property<float>("v:values");
	for( auto v : geometry->vertices() ) vprop[v] = newValues[v.idx()];
}

void Spherelib::Sphere::assignLocalFrame( int /*tracks*/, int /*sectors*/ )
{
	//Surface_mesh::Vertex_property<Vector3> points = geometry->vertex_property<Vector3>("v:point");
	//Surface_mesh::Vertex_property<Vector3> projected = geometry->vertex_property<Vector3>("v:projected", Vector3(0,0,0));

	//std::vector<float> values = this->values();
	//SurfaceMesh::Vertex major( std::max_element(values.begin(),values.end()) - values.begin() );

	//Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(points[major], Vector3(0,0,1));

	// Discretized
	//grid = Spherelib::Sphere::createGrid(values, tracks, sectors);
	//grid.majorAxis = points[major];

	// Align to computed local frame
	//grid.align();
}

Spherelib::RadialGrid Spherelib::Sphere::createGrid( const std::vector<double> & fromValues, int tracks, int sectors )
{
	const Surface_mesh::Vertex_property<Vector3> points = geometry->vertex_property<Vector3>("v:point");

	SurfaceMesh::Vertex majorVert( std::max_element(fromValues.begin(),fromValues.end()) - fromValues.begin() );
	Vector3 majorAxis = points[majorVert];

	Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(majorAxis, Vector3(0,0,1));

	std::vector<Vector3> projected( numPoints, Vector3(0,0,0) );

	// Project rays to global 2D plane
	for( auto v : geometry->vertices() )
	{
		// Rotated
		Vector3 pr = q * points[v]; 

		// Replace 'z' coordinate
		Vector3 pj = pr;
		pj.z() = fromValues[v.idx()];

		// Encode side
		if(points[v].dot(majorAxis) < 0) pj.z() *= -1;

		projected[v.idx()] = pj;
	}

	return Spherelib::RadialGrid::createGrid (projected, tracks, sectors, majorAxis);
}

void Spherelib::Sphere::draw()
{
	geometry->updateBoundingBox();
	geometry->update_face_normals();
	geometry->update_vertex_normals();

	Surface_mesh * mesh = geometry;
	Surface_mesh::Vertex_property<Vector3> points = mesh->vertex_property<Vector3>("v:point");
	Surface_mesh::Face_property<Vector3> fnormals = mesh->face_property<Vector3>("f:normal");
	Surface_mesh::Face_iterator fit, fend = mesh->faces_end();
	Surface_mesh::Vertex_around_face_circulator fvit, fvend;

	std::vector<float> values = this->values();

	glBegin(GL_TRIANGLES);
	for (fit=mesh->faces_begin(); fit!=fend; ++fit){
		glNormal3dv( fnormals[fit].data() );
		fvit = fvend = mesh->vertices(fit);
		do{ 
			int vidx = SurfaceMesh::Vertex(fvit).idx();
			float d = values[vidx];
			glColor3d( d, d, d );
			Vector3 p = points[fvit];
			//p += (p * (d - 0.5) * 0.2);
			glVertex3dv( p.data() );
		} while (++fvit != fvend);
	}
	glEnd();
}

Spherelib::Sphere::~Sphere()
{
	delete geometry;
}

Spherelib::RadialGrid Spherelib::RadialGrid::createGrid(const std::vector<Vector3> & projected, int numTracks, int numSectors, Vector3 majorAxis)
{
	RadialGrid grid;
	grid.tracks = numTracks;
	grid.sectors = numSectors;
	grid.majorAxis = majorAxis;

	grid.front_values.resize( numTracks, std::vector<double>(numSectors, 0) );
	grid.back_values = grid.front_values;

	int count = 0;
	double sum = 0;

	for(auto p : projected)
	{
		double r = Eigen::Vector2d(p.x(),p.y()).norm();

		int track = (r == 0) ? 0 : std::ceil(r*numTracks) - 1;
		track = std::min(numTracks-1, track);

		double angle = std::atan2(p.x(),p.y());
		if(angle < 0) angle += (M_PI * 2);

		double s = angle / (M_PI * 2);
		int sector = std::min(numSectors-1, int(numSectors * s));

		// Clock-wise
		//sector = (numSectors-1) - sector;

		double val = p.z();

		// Front or back
		if( val > 0 )
			grid.front_values[track][sector] = std::max(grid.front_values[track][sector], val);
		else
			grid.back_values[track][sector] = std::max(grid.back_values[track][sector], std::abs(val));

		sum += std::abs(val);
		count++;
	}

	grid.average = sum / count;

	return grid;
}

template<typename VectorType>
inline void rotateVector( VectorType & v, int steps ){
	if(steps < 0) std::rotate(v.begin(),v.begin()+(v.size()+steps),v.end());
	else std::rotate(v.begin(),v.begin()+steps,v.end());
}

void Spherelib::RadialGrid::rotate(int steps){
	for(auto & track : front_values) rotateVector( track, steps );
	for(auto & track : back_values) rotateVector( track, steps );
}

void Spherelib::RadialGrid::align()
{
	std::vector< std::vector< std::vector<double> >* > values_set;
	values_set.push_back( &front_values );
	values_set.push_back( &back_values );

	int final_step = 1;

	while(abs(final_step) > 0)
	{
		std::vector<double> data, weights;

		for(auto values : values_set)
		{
			for(int track = 0; track < tracks; track++)
			{
				int sector = std::max_element(values->at(track).begin(), values->at(track).end()) - values->at(track).begin();
				double val = values->at(track).at(sector);

				// Weight based on distance to equator
				double weight = double(track) / (tracks-1);
				weight = pow(weight, 2);
				val *= weight;

				// Consider above average features
				if(val < average) 
					continue;

				// Shortest step size
				if(sector > sectors * 0.5) 
					sector = -((sectors-1) - sector);

				data.push_back( sector );
				weights.push_back( val );
			}
		}

		double numerator = 0, denominator = 1;

		for(size_t i = 0; i < data.size(); i++)
		{
			numerator += (data[i] * weights[i]);
			denominator += weights[i];
		}

		if(denominator == 0)
		{
			final_step = 0;
		}
		else
		{
			double weighted_step = (numerator / denominator);
			final_step = int(weighted_step);
			//qDebug() << "Weighted step:" << QString::number(weighted_step) << ", int Step: " << QString::number(final_step);

			rotate( final_step );
		}
	}
}

void Spherelib::RadialGrid::draw(QPainter & painter)
{
	int width = 300, height = 300;
	int thickness = 10;

	int padding = thickness;
	int deltaY = 0;

	std::vector< std::vector< std::vector<double> > > values_set;
	values_set.push_back( front_values );
	values_set.push_back( back_values );

	for(auto values : values_set)
	{
		for(int track = 0; track < tracks; track++)
		{
			QRect trackRect(0,0, padding + ((thickness-1) * 2 * track), padding + ((thickness-1) * 2 * track));
			trackRect.moveCenter( QPoint(width*0.5,height*0.5) + QPoint(0,deltaY) );

			for(int sector = 0; sector < sectors; sector++)
			{
				double val = values[track][sector];

				QColor c = QColor::fromRgbF(val,val,val);

				QPen pen(c, thickness);
				pen.setCapStyle(Qt::FlatCap);
				painter.setPen( pen );

				int spanAngle = 360 / sectors;
				int startAngle = sector * spanAngle;

				painter.drawArc( trackRect, startAngle * 16, spanAngle * 16 );
			}

			// Highlight largest
			int largestSector = std::max_element(values[track].begin(), values[track].end()) - values[track].begin();
			double val = values[track][largestSector];

			if( val > average )
			{
				QColor c = QColor::fromRgbF(0,val,0);
				QPen pen(c, thickness);
				pen.setCapStyle(Qt::FlatCap);
				painter.setPen( pen );
				int spanAngle = 360 / sectors;
				int startAngle = largestSector * spanAngle;
				painter.drawArc( trackRect, startAngle * 16, spanAngle * 16 );
			}
		}

		deltaY += height;
	}
}

std::vector<double> Spherelib::RadialGrid::alignedValues()
{
	this->align();

	std::vector<double> result;
    for(auto & track : front_values)result.insert(result.end(), track.begin(), track.end());
    for(auto & track : back_values)result.insert(result.end(), track.begin(), track.end());
	return result;
}
