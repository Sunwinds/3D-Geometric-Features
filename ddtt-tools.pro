TEMPLATE = subdirs
CONFIG += ordered

# Libraries
#SUBDIRS += NURBS

# Plugins
#SUBDIRS += voxel_resampler

SUBDIRS  = agd          # Average geodesic distance
SUBDIRS += bdf   	# Biharmonic distance function
SUBDIRS += sdf   	# Shape Diameter Function
SUBDIRS += curvature    # Curvature
SUBDIRS += repair	# Basic mesh repair
#SUBDIRS += segmentation # Segmentation via skeletons
SUBDIRS += symmetry	# Basic symmetry analysis
SUBDIRS += PCAShapeEst  # Segment level PCA based shape estimation
SUBDIRS += ConformalFactor # Comformal factor

#SUBDIRS += test        # Performance test
SUBDIRS += BatchProc    # Batch process (Import parts, Normalization, Calculate Features, Output Affinity Matrix W(i,j) = exp(-D(S(i),S(j))/2*sigma))

# Dependency map
