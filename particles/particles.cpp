#include "particles.h"
#include "particles-widget.h"
#include "ui_particles-widget.h"
#include <omp.h>
#include <numeric>

#include "SurfaceMeshModel.h"
#include "SurfaceMeshHelper.h"
using namespace SurfaceMesh;

#include <QOpenGLShaderProgram>
#include <QWidget>
#include <QFileDialog>
#include <QElapsedTimer>

#define AlphaBlend(alpha, start, end) ( ((1-alpha) * start) + (alpha * end) )
QTimer * timer = NULL;

#include "ParticleMesh.h"
#include "myglobals.h"

#include "Raytracing.h"

#include "kmeans.h"

#include "spherelib.h"
#include "SphericalHarmonic.h"

std::vector<Spherelib::Sphere> spheres;

void particles::create()
{
	if( widget ) return;

	//drawArea()->setAxisIsDrawn(true);

    ModePluginDockWidget * dockwidget = new ModePluginDockWidget("Particles", mainWindow());
	ParticlesWidget * pw = new ParticlesWidget();
    widget = pw;

    dockwidget->setWidget( widget );
    mainWindow()->addDockWidget(Qt::RightDockWidgetArea, dockwidget);

	// General Tests
	connect(pw->ui->testButton, &QPushButton::released, [=]{
		//for(auto p : sphere_fibonacci_points( 100 ))drawArea()->drawPoint(p, 5);
		//drawArea()->update();
	});

	// Load and process shapes:
	connect(pw->ui->loadShapes, &QPushButton::released, [=]{
		QStringList files = QFileDialog::getOpenFileNames(nullptr, "Open Shapes", "", "All Supported (*.obj *.off)");

		for(auto filename : files)
		{
			SurfaceMeshModel fromMesh( filename, QFileInfo(filename).baseName() );
			fromMesh.read( filename.toStdString() );

			pw->pmeshes.push_back( new ParticleMesh( &fromMesh, pw->ui->gridsize->value() ) );
		}

		emit( pw->shapesLoaded() );
	});

	// Post-processing
	connect(pw, SIGNAL(shapesLoaded()), SLOT(processShapes()));
	connect(pw->ui->processShapesButton, SIGNAL(clicked()), SLOT(processShapes()));
}

void particles::processShapes()
{
	if(!widget) return;

	ParticlesWidget * pw = (ParticlesWidget *) widget;
	pw->isReady = false;

	// Show voxelized meshes
	//for(auto & s : pw->pmeshes) document()->addModel(s->surface_mesh->clone());

	QElapsedTimer allTimer; allTimer.start();

	mainWindow()->setStatusBarMessage("Shapes loaded, now processing..");
	qApp->processEvents();

	/// Fixed set of ray directions:
	//std::vector< Eigen::Vector3d > sampledRayDirections = sphere_fibonacci_points( perSampleRaysCount );

	// Uniform sampling on the sphere
	Spherelib::Sphere sphere( pw->ui->sphereResolution->value() );
	std::vector< Eigen::Vector3d > sampledRayDirections  = sphere.rays();
	size_t perSampleRaysCount = sampledRayDirections.size();

	// DEBUG:
	spheres.clear();
	int selectedParticle = pw->ui->selectedParticle->value();
	if(selectedParticle >= 0) spheres.push_back( sphere );

	// Rotation invariant descriptor
	SphericalHarmonic<Vector3,float> sh( std::max(1,pw->ui->bands->value()) );
	std::vector< SHSample<Vector3,float> > sh_samples;
	sh.SH_setup_spherical( sampledRayDirections, sh_samples );

	if( true )
	{
		for(auto & s : pw->pmeshes)
		{
			QElapsedTimer timer; timer.start();

			// Hit results are saved as vector of distances
			std::vector< std::vector<float> > & descriptor = (s->desc = std::vector< std::vector<float> >( 
				s->particles.size(), std::vector<float>( perSampleRaysCount ) ));
				
			int rayCount = 0;

			// Accelerated raytracing
			{
				raytracing::Raytracing<Eigen::Vector3d> rt( s->surface_mesh );

				// Shoot rays around all particles
				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++)
				{
					auto & p = s->particles[pi];

					int r = 0;
					for(auto d : sampledRayDirections)
					{
						raytracing::RayHit hit = rt.hit( p.pos, d );

						descriptor[pi][r++] = hit.distance;
						rayCount++;
					}
				}
			}

			// Average of neighbors
			int avgNeighIters = pw->ui->avgNeighIters->value();
			if( avgNeighIters )
			{
				size_t gridsize = s->grid.gridsize;

				for(int it = 0; it < avgNeighIters; it++)
				{
					auto smoothDesc = descriptor;
					{
						#pragma omp parallel for
						for(int pi = 0; pi < (int)s->particles.size(); pi++)
						{
							const auto & p = s->particles[pi];

							std::vector<float> sum( descriptor[p.id].size(), 0 );
							int count = 0;

							unsigned int x,y,z;
							mortonDecode(p.morton, x, y, z);

							for(int u = -1; u <= 1; u++){
								for(int v = -1; v <= 1; v++){
									for(int w = -1; w <= 1; w++)
									{
										Eigen::Vector3i c(x + u, y + v, z + w);

										// Skip outside grid
										if(c.x() < 0 || c.y() < 0 || c.z() < 0) continue;
										if(c.x() > gridsize-1 || c.y() > gridsize-1 || c.z() > gridsize-1) continue;

										uint64_t mcode = mortonEncode_LUT(c.x(),c.y(),c.z());
										if(!s->grid.occupied[ mcode ]) continue;

										size_t nj = s->mortonToParticleID[ mcode ];

										std::vector<float> nei = descriptor[nj];

										// Add this neighbor to sum
										for(size_t j = 0; j < sum.size(); j++) sum[j] += nei[j];
										count++;
									}
								}
							}

							// divide over count
							for(size_t j = 0; j < sum.size(); j++) sum[j] /= count;

							smoothDesc[p.id] = sum;
						}
					}
					descriptor = smoothDesc;
				}
			}

			// Report
			mainWindow()->setStatusBarMessage( QString("Ray tracing: rays (%1) / time (%2 ms)").arg( rayCount ).arg( timer.elapsed() ) );
			timer.restart();

			// Inner voxels more visible
			double global_min = DBL_MAX;
			for(auto & p : s->particles)
			{
				std::vector<float> & desc = descriptor[p.id];
				double min_desc = *std::min_element(desc.begin(),desc.end());
				global_min = std::min(global_min, min_desc);
			}

			// Extract one of the descriptors
			if(selectedParticle >= 0)
			{
				int idx = std::min(selectedParticle, (int)s->particles.size()- 1);
				spheres.back().setValues( descriptor[idx] );
				spheres.back().normalizeValues();
			}
            else
                spheres.clear();

			// Rotation invariant description
			if( pw->ui->bands->value() > 0 )
			{
				#pragma omp parallel for
				for(int pi = 0; pi < (int)s->particles.size(); pi++)
				{
					auto & p = s->particles[pi];

					std::vector<float> & desc = descriptor[p.id];

					p.direction = sampledRayDirections[ std::max_element(desc.begin(),desc.end()) - desc.begin() ];

					// Normalize response
					if( pw->ui->normalizeSphereFn->isChecked() )
					{
						double max_desc = *std::max_element(desc.begin(),desc.end());
						double min_desc = *std::min_element(desc.begin(),desc.end());
						for(auto & d : desc) d /= max_desc;

						p.alpha = min_desc / global_min;
					}

					std::vector<float> coeff;
					sh.SH_project_function(desc, sh_samples, coeff);
						
					desc = sh.SH_signature(coeff);

					//auto grid = sphere.createGrid(desc,pw->ui->tracks->value(),pw->ui->sectors->value());
					//desc = grid.alignedValues();
				}
			}

			// Report
			mainWindow()->setStatusBarMessage( QString("Alignment took (%1 ms)").arg( timer.elapsed() ) );
			timer.restart();

			// Debug
			{
				/*starlab::LineSegments * vs = new starlab::LineSegments;
				for(auto particle : s->particles){		
					Eigen::Vector4d color(abs(particle.direction[0]), abs(particle.direction[1]), abs(particle.direction[2]), particle.alpha);
					color[0] *= color[0];color[1] *= color[1];color[2] *= color[2];
					vs->addLine(particle.pos,  Vector3(particle.pos + particle.direction * 0.005), QColor::fromRgbF(color[0],color[1],color[2],1));
				}
				s->debug.push_back(vs);*/
			}
		}
	}

	// k-means clustering
	if( true )
	{
		// Collect descriptors
		std::vector< std::vector<float> > allDesc;
		for(auto & s : pw->pmeshes)
			allDesc.insert(allDesc.end(), s->desc.begin(), s->desc.end());

		// Cluster
		typedef clustering::l2norm_squared< std::vector<float> > dist_fn;
		clustering::kmeans< std::vector< std::vector<float> >, dist_fn > km(allDesc, pw->ui->kclusters->value());
		km.run(100, 0.01);

		// Assign classes
		int pi = 0;
		for(auto & s : pw->pmeshes)
			for(auto & p : s->particles)
				p.flag = (int) km.clusters()[pi++];
	}

	/// [ Correspondence ] match particles
	if( false )
	{
		for(size_t i = 0; i < pw->pmeshes.size(); i++){
			for(size_t j = i+1; j < pw->pmeshes.size(); j++){
				NanoKdTree * itree = pw->pmeshes[i]->kdtree;
				NanoKdTree * jtree = pw->pmeshes[j]->kdtree;

				for(auto & iparticle : pw->pmeshes[i]->particles)
				{
					if( true )
					{
						iparticle.correspondence = jtree->closest( iparticle.relativePos );
					}
					else
					{
						// Experiment
						KDResults matches;
						jtree->ball_search( iparticle.relativePos, 0.2, matches );

						QMap<double, size_t> measures;
						for(auto p : matches) 
						{
							double weight = p.second;
							//double dist = dist_fn()( pw->pmeshes[i]->desc[iparticle.id], pw->pmeshes[j]->desc[p.first] );
							double dist = 1;
							measures[ weight * dist ] = p.first;
						}
						iparticle.correspondence = measures[ measures.keys().front() ];
					}
				}

				for(auto & jparticle : pw->pmeshes[j]->particles)
				{
					if( true )
					{
						jparticle.correspondence = itree->closest( jparticle.relativePos );
					}
					else
					{
						// Experiment
						KDResults matches;
						itree->ball_search( jparticle.relativePos, 0.2, matches );

						QMap<double, size_t> measures;
						for(auto p : matches) 
						{
							double weight = p.second;
							//double dist = dist_fn()( pw->pmeshes[j]->desc[jparticle.id], pw->pmeshes[i]->desc[p.first] );
							double dist = 1;
							measures[ weight * dist ] = p.first;
						}
						jparticle.correspondence = measures[ measures.keys().front() ];
					}
				}
			}
		}
	}

	// Report
	mainWindow()->setStatusBarMessage( QString("All time (%1 ms)").arg( allTimer.elapsed() ) );

	pw->isReady = true;

	drawArea()->update();
}

void particles::decorate()
{
	ParticlesWidget * pwidget = (ParticlesWidget*) widget;
	if(!pwidget || !pwidget->isReady || pwidget->pmeshes.size() < 1) return;

	for(auto & sphere : spheres)
	{
		sphere.draw();
	}

	//if(pwidget->pmeshes.size() < 2) 
	{
		glPushMatrix();

		// Evaluation
		for( auto s : pwidget->pmeshes )
		{
			s->drawParticles();
			s->drawDebug( *drawArea() );

			glTranslated(s->bbox.sizes().x() * 1.1, 0, 0);
		}

		glPopMatrix();

		return;
	}

	// Experimental
	static bool isForward = true;
	static double alpha = 0;
	if(isForward) alpha += 0.005; else alpha -= 0.005;
	if(alpha > 1.0){alpha = 1.0;isForward = false;}
	if(alpha < 0.0){alpha = 0.0;isForward = true;}

	// Prepare scene once
	Eigen::AlignedBox3d largeBox;
	for(auto pmesh : pwidget->pmeshes) largeBox.extend(pmesh->bbox);
	if(timer == NULL){
		drawArea()->setSceneRadius( largeBox.sizes().norm() * 2 );
		drawArea()->showEntireScene();

		timer = new QTimer;
		connect(timer, &QTimer::timeout, [=]() { drawArea()->update(); });
		timer->start(20);
	}

	//for(auto pmesh : pwidget->pmeshes) pmesh->drawParticles();
	for(auto pmesh : pwidget->pmeshes) pmesh->drawDebug( *drawArea() );

	// Collect points
	std::vector<Eigen::Vector3f> mixedPoints;
	for(size_t i = 0; i < pwidget->pmeshes.size(); i++)
	{
		for(size_t j = i+1; j < pwidget->pmeshes.size(); j++)
		{
			ParticleMesh * imesh = pwidget->pmeshes[i];
			ParticleMesh * jmesh = pwidget->pmeshes[j];

			for(auto particle : imesh->particles){
				Vector3 p = AlphaBlend(alpha, particle.pos, jmesh->particles[particle.correspondence].pos);
				mixedPoints.push_back( p.cast<float>() );
			}

			for(auto particle : jmesh->particles){
				Vector3 p = AlphaBlend((1.0-alpha), particle.pos, imesh->particles[particle.correspondence].pos);
				mixedPoints.push_back( p.cast<float>() );
			}
		}
	}

	// Setup shader
	static QMap<QString, int> location;
	static QOpenGLShaderProgram * program = NULL;
	if( !program )
	{
		program = new QOpenGLShaderProgram;
		program->addShaderFromSourceCode(QOpenGLShader::Vertex, 
			"attribute highp vec4 vertex;\n"
			"uniform highp mat4 matrix;\n"
			"void main(void)\n"
			"{\n"
			"   gl_Position = vertex * matrix;\n"
			"}");
		program->addShaderFromSourceCode(QOpenGLShader::Fragment, 
			"#version 330\n"
			"out vec4 vFragColor;\n"
			"uniform mediump vec4 color;\n"
			"uniform mediump vec3 lightDir;\n"
			"void main(void)\n"
			"{\n" 
			"vec3 N;\n"
			"N.xy = gl_PointCoord* 2.0 - vec2(1.0);    \n"
			"float mag = dot(N.xy, N.xy);\n"
			"if (mag > 1.0) discard;   // kill pixels outside circle\n"
			"N.z = sqrt(1.0-mag);\n"
			"float diffuse = max(0.0, dot(lightDir, N));\n"
			"vFragColor = color * diffuse;\n"
			"}"
			);
		program->link();
		program->bind();

		location["vertex"] = program->attributeLocation("vertex");
		location["matrix"] = program->uniformLocation("matrix");
		location["color"] = program->uniformLocation("color");
		location["lightDir"] = program->uniformLocation("lightDir");
	}
	else
	{
		GLdouble m[16];
		drawArea()->camera()->getModelViewProjectionMatrix(m);
		QMatrix4x4 pmvMatrix(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7],m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);

		QColor color(255,0,0,255);
		QVector3D lightdir(0,0,1);

		program->bind();
		program->setAttributeArray(location["vertex"], mixedPoints.front().data(), 3);
		program->enableAttributeArray(location["vertex"]);
		program->setUniformValue(location["matrix"], pmvMatrix);
		program->setUniformValue(location["color"], color);
		program->setUniformValue(location["lightDir"], lightdir);
	}

	glPointSize(10);
	glEnable(GL_POINT_SPRITE);

	glDrawArrays(GL_POINTS, 0, (int)mixedPoints.size());

	if( program )
	{
		program->disableAttributeArray( location["vertex"] );
		program->release();
	}
}

bool particles::keyPressEvent(QKeyEvent*e)
{
	if(e->key() == Qt::Key_Space)
		spheres.clear();

	drawArea()->update();
	return true;
}
