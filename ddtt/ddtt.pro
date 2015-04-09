TEMPLATE = subdirs

# Libraries
SUBDIRS += NURBS
SUBDIRS += GlSplatRendererLib
SUBDIRS += Reconstruction
SUBDIRS += StructureGraphLib
SUBDIRS += AuctionLIB
SUBDIRS += bowlib

# Main plugin
SUBDIRS += ddtt-plugin functional

# Aux. plugins
SUBDIRS += empty-mesh

ddtt-plugin.depends = AuctionLIB StructureGraphLib bowlib
