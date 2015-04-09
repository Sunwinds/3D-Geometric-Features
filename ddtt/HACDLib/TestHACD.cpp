#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <iostream>
/*!
**
** Copyright (c) 20011 by John W. Ratcliff mailto:jratcliffscarab@gmail.com
**
**
** The MIT license:
**
** Permission is hereby granted, free of charge, to any person obtaining a copy
** of this software and associated documentation files (the "Software"), to deal
** in the Software without restriction, including without limitation the rights
** to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
** copies of the Software, and to permit persons to whom the Software is furnished
** to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in all
** copies or substantial portions of the Software.

** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
** WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
** CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

#pragma warning(disable:4996 4100)

#include "PlatformConfigHACD.h"
#include "wavefront.h"
#include "HACD.h"
#include "JobSwarm.h"

namespace HACD
{
	HACD_API *gHACD = NULL;
};


float getFloatArg(int arg,int argc,const char **argv)
{
	float ret = 0;
	if ( arg < argc )
	{
		ret = (float)atof(argv[arg]);
	}
	else
	{
		printf("Error: Missing input argument value at argument location %d.\r\n",arg+1);
	}
	return ret;
}

int getIntArg(int arg,int argc,const char **argv)
{
	int ret = 0;
	if ( arg < argc )
	{
		ret = atoi(argv[arg]);
	}
	else
	{
		printf("Error: Missing input argument value at argument location %d.\r\n",arg+1);
	}
	return ret;
}

class MyCallback : public hacd::ICallback
{
public:
	MyCallback(void)
	{
		mLastPercent = 0xFFFFFFFF;
		mLastMessage = NULL;
	}

	virtual bool Cancelled()
	{
		// Don't have a cancel button in the test console app.
		return false;
	}

	virtual void ReportProgress(const char* message, hacd::HaF32 progress)
	{
		hacd::HaU32 percent = (hacd::HaU32)(progress*100.0f);
		if ( percent != mLastPercent || message != mLastMessage )
		{
			if ( message == mLastMessage )
			{
				printf("%c%c%c%c",8,8,8,8);
				printf("%3d%c", percent,'%' );
			}
			else
			{
				printf("\r\n%s : %3d", message, percent);
			}
			mLastMessage = message;
			mLastPercent = percent;
		}
	}
	const char *	mLastMessage;
	hacd::HaU32		mLastPercent;
};

void TestHACD_main(int argc,const char ** argv)
{
	if ( argc == 1 )
	{
		printf("Usage: TestHACD <wavefront.obj> (options)\r\n");
		printf("\r\n");
		printf("Options:\r\n");
		printf("-v	<count>	: Max Hull Vertices (default 64)\r\n");
		printf("-m	<count>	: Maximum number of hulls output from HACD (default 256)\r\n");
		printf("-merge <count> : Maximum number of hulls after merging the HACD result.\r\n");
		printf("-mergethreshold <volume> : Threshold below which hulls are merged if they are smaller than the given volume.\r\n");
		printf("-c <concavity> : Between 0 and 1 are good ranges to try; default is 0.2.  The smaller the number, the more convex hulls are produced.\r\n");
		printf("-b <backFaceDistanceFactor : The back face distance factor, default is 0.2\r\n");
		printf("-t <threadCount> : Specifies the number of threads to use for a multi-threaded implementation\r\n");
		printf("-d <decompositionDepth> : Specifies the decomposition depth; uses legacy ACD instead of HACD\r\n");
		printf("-n : Normalize the input mesh.\r\n");
		printf("-f : Use legacy 'fast' version of HACD\r\n");
		printf("\r\n");
		printf("Example: TestHACD hornbug.obj -m 40 -v 64\r\n");
		printf("\r\n");
	}
	else
	{
		HACD::HACD_API::Desc desc;
		const char *wavefront = argv[1];
		int scan = 2;
		int threadCount = 0;


		while ( scan < argc )
		{
			const char *option = argv[scan];
			if ( strcmp(option,"-v") == 0 )
			{
				desc.mMaxHullVertices = getIntArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-t") == 0 )
			{
				threadCount  = getIntArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-d") == 0 )
			{
				desc.mDecompositionDepth  = getIntArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-b") == 0 )
			{
				desc.mBackFaceDistanceFactor = getFloatArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-m") == 0 )
			{
				desc.mMaxHullCount = getIntArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-merge") == 0 )
			{
				desc.mMaxMergeHullCount = getIntArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-mergethreshold") == 0 )
			{
				desc.mSmallClusterThreshold = getFloatArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-c") == 0 )
			{
				desc.mConcavity = getFloatArg(scan+1,argc,argv);
				scan+=2;
			}
			else if ( strcmp(option,"-n") == 0 )
			{
				desc.mNormalizeInputMesh = true;
				scan++;
			}
			else if ( strcmp(option,"-f") == 0 )
			{
				desc.mUseFastVersion = true;
				scan++;
			}
			else
			{
				scan++;
			}
		}
		JOB_SWARM::JobSwarmContext *jobSwarmContext = NULL;
		if ( threadCount )
		{
			jobSwarmContext = JOB_SWARM::createJobSwarmContext(threadCount);
//			jobSwarmContext->setUseThreads(false);
		}
		HACD::gHACD = HACD::createHACD_API();
		if  ( HACD::gHACD )
		{
			MyCallback callback;
			WavefrontObj obj;
			unsigned int tcount = obj.loadObj(wavefront);
			if ( tcount )
			{
				desc.mTriangleCount = obj.mTriCount;
				desc.mVertexCount = obj.mVertexCount;
				desc.mIndices = (hacd::HaU32 *)obj.mIndices;
				desc.mVertices = obj.mVertices;
			}
			desc.mCallback = static_cast<hacd::ICallback*>(&callback);
			if ( desc.mTriangleCount )
			{
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
				if ( desc.mDecompositionDepth )
				{
					printf("DecompositionDepth       : %d\r\n", desc.mDecompositionDepth );
				}

				if ( jobSwarmContext )
				{
					printf("Parallel Solver TheadCount: %d\r\n", threadCount );
				}

				if ( desc.mNormalizeInputMesh )
				{
					printf("Normalizing input mesh.\r\n");
				}

				if ( desc.mUseFastVersion )
				{
					printf("Using 'fast' version of HACD.\r\n");
				}

				printf("\r\n");

				desc.mJobSwarmContext = jobSwarmContext;

				hacd::HaU32 hullCount = HACD::gHACD->performHACD(desc);

				if ( hullCount != 0 )
				{
					printf("\r\n");
					printf("Produced %d output convex hulls. Saving output to 'ConvexDecomposition.obj' for review.\r\n", hullCount );
					FILE *fph = fopen("ConvexDecomposition.obj", "wb");
					if ( fph )
					{
						fprintf(fph,"# Input mesh '%s' produced %d convex hulls.\r\n", wavefront, hullCount );
						hacd::HaU32 *baseVertex = new hacd::HaU32[hullCount];
						hacd::HaU32 vertexCount = 0;
						for (hacd::HaU32 i=0; i<hullCount; i++)
						{
							const HACD::HACD_API::Hull *hull = HACD::gHACD->getHull(i);
							if ( hull )
							{
								baseVertex[i] = vertexCount;
								fprintf(fph,"## Hull %d has %d vertices.\r\n", i+1, hull->mVertexCount );
								for (hacd::HaU32 i=0; i<hull->mVertexCount; i++)
								{
									const hacd::HaF32 *p = &hull->mVertices[i*3];
									fprintf(fph,"v %0.9f %0.9f %0.9f\r\n", p[0], p[1], p[2] );
								}
								vertexCount+=hull->mVertexCount;
							}
						}
						for (hacd::HaU32 i=0; i<hullCount; i++)
						{
							const HACD::HACD_API::Hull *hull = HACD::gHACD->getHull(i);
							if ( hull )
							{
								hacd::HaU32 startVertex = baseVertex[i];
								fprintf(fph,"# Convex Hull %d contains %d triangles and %d vertices.  Starting vertex index is: %d It has a volume of: %0.9f\r\n", i+1, hull->mTriangleCount, hull->mVertexCount, startVertex);
								for (hacd::HaU32 j=0; j<hull->mTriangleCount; j++)
								{
									hacd::HaU32 i1 = hull->mIndices[j*3+0]+startVertex+1;
									hacd::HaU32 i2 = hull->mIndices[j*3+1]+startVertex+1;
									hacd::HaU32 i3 = hull->mIndices[j*3+2]+startVertex+1;
									fprintf(fph,"f %d %d %d\r\n", i1, i2, i3 );
								}
							}
						}
					}
					else
					{
						printf("Failed to open output file.\r\n");
					}
				}
			}
			else
			{
				printf("Failed to load Wavefront OBJ file '%s'\r\n",wavefront);
			}
		}
		else
		{
			printf("Failed to load the HACD DLL\r\n");
		}
		if ( jobSwarmContext )
		{
			JOB_SWARM::releaseJobSwarmContext(jobSwarmContext);
		}
	}
}
