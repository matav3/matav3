#############################################################################
# Makefile for building: matav3
# Generated by qmake (2.01a) (Qt 4.8.1) on: Sun Nov 18 19:36:08 2012
# Project:  matav3.pro
# Template: app
# Command: /usr/bin/qmake -o Makefile matav3.pro
#############################################################################

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_WEBKIT -DQT_NO_DEBUG -DQT_GUI_LIB -DQT_CORE_LIB -DQT_SHARED
CFLAGS        = -m64 -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)
CXXFLAGS      = -fpermissive -w -g -O2 -Wall -W -D_REENTRANT $(DEFINES)
INCPATH       = -I/usr/share/qt4/mkspecs/linux-g++-64 -I. -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I/usr/include/qt4 -I. -Isrc -Isrc/treeViewer -I. -I.
LINK          = g++
LFLAGS        = -m64 -Wl,-O1
LIBS          = $(SUBLIBS)  -L/usr/lib/x86_64-linux-gnu -lDIGFile -lxerces-c -lQtGui -lQtCore -lpthread 
AR            = ar cqs
RANLIB        = 
QMAKE         = /usr/bin/qmake
TAR           = tar -cf
COMPRESS      = gzip -9f
COPY          = cp -f
SED           = sed
COPY_FILE     = $(COPY)
COPY_DIR      = $(COPY) -r
STRIP         = strip
INSTALL_FILE  = install -m 644 -p
INSTALL_DIR   = $(COPY_DIR)
INSTALL_PROGRAM = install -m 755 -p
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cpp \
		mainwindow.cpp \
		src/treeViewer/Edge.cpp \
		src/treeViewer/Node.cpp \
		src/treeViewer/TreeManager.cpp \
		src/treeViewer/TreeScene.cpp \
		src/treeViewer/TreeView.cpp moc_mainwindow.cpp \
		moc_TreeManager.cpp \
		moc_TreeScene.cpp \
		moc_TreeView.cpp \
		qrc_application.cpp
OBJECTS       = main.o \
		mainwindow.o \
		Edge.o \
		Node.o \
		TreeManager.o \
		TreeScene.o \
		TreeView.o \
		moc_mainwindow.o \
		moc_TreeManager.o \
		moc_TreeScene.o \
		moc_TreeView.o \
		qrc_application.o
DIST          = /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/modules/qt_webkit_version.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		matav3.pro
QMAKE_TARGET  = matav3
DESTDIR       = 
TARGET        = matav3

first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

all: Makefile $(TARGET)

$(TARGET): ui_mainwindow.h $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: matav3.pro  /usr/share/qt4/mkspecs/linux-g++-64/qmake.conf /usr/share/qt4/mkspecs/common/unix.conf \
		/usr/share/qt4/mkspecs/common/linux.conf \
		/usr/share/qt4/mkspecs/common/gcc-base.conf \
		/usr/share/qt4/mkspecs/common/gcc-base-unix.conf \
		/usr/share/qt4/mkspecs/common/g++-base.conf \
		/usr/share/qt4/mkspecs/common/g++-unix.conf \
		/usr/share/qt4/mkspecs/qconfig.pri \
		/usr/share/qt4/mkspecs/modules/qt_webkit_version.pri \
		/usr/share/qt4/mkspecs/features/qt_functions.prf \
		/usr/share/qt4/mkspecs/features/qt_config.prf \
		/usr/share/qt4/mkspecs/features/exclusive_builds.prf \
		/usr/share/qt4/mkspecs/features/default_pre.prf \
		/usr/share/qt4/mkspecs/features/release.prf \
		/usr/share/qt4/mkspecs/features/default_post.prf \
		/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf \
		/usr/share/qt4/mkspecs/features/warn_on.prf \
		/usr/share/qt4/mkspecs/features/qt.prf \
		/usr/share/qt4/mkspecs/features/unix/thread.prf \
		/usr/share/qt4/mkspecs/features/moc.prf \
		/usr/share/qt4/mkspecs/features/resources.prf \
		/usr/share/qt4/mkspecs/features/uic.prf \
		/usr/share/qt4/mkspecs/features/yacc.prf \
		/usr/share/qt4/mkspecs/features/lex.prf \
		/usr/share/qt4/mkspecs/features/include_source_dir.prf \
		/usr/lib/x86_64-linux-gnu/libQtGui.prl \
		/usr/lib/x86_64-linux-gnu/libQtCore.prl
	$(QMAKE) -o Makefile matav3.pro
/usr/share/qt4/mkspecs/common/unix.conf:
/usr/share/qt4/mkspecs/common/linux.conf:
/usr/share/qt4/mkspecs/common/gcc-base.conf:
/usr/share/qt4/mkspecs/common/gcc-base-unix.conf:
/usr/share/qt4/mkspecs/common/g++-base.conf:
/usr/share/qt4/mkspecs/common/g++-unix.conf:
/usr/share/qt4/mkspecs/qconfig.pri:
/usr/share/qt4/mkspecs/modules/qt_webkit_version.pri:
/usr/share/qt4/mkspecs/features/qt_functions.prf:
/usr/share/qt4/mkspecs/features/qt_config.prf:
/usr/share/qt4/mkspecs/features/exclusive_builds.prf:
/usr/share/qt4/mkspecs/features/default_pre.prf:
/usr/share/qt4/mkspecs/features/release.prf:
/usr/share/qt4/mkspecs/features/default_post.prf:
/usr/share/qt4/mkspecs/features/unix/gdb_dwarf_index.prf:
/usr/share/qt4/mkspecs/features/warn_on.prf:
/usr/share/qt4/mkspecs/features/qt.prf:
/usr/share/qt4/mkspecs/features/unix/thread.prf:
/usr/share/qt4/mkspecs/features/moc.prf:
/usr/share/qt4/mkspecs/features/resources.prf:
/usr/share/qt4/mkspecs/features/uic.prf:
/usr/share/qt4/mkspecs/features/yacc.prf:
/usr/share/qt4/mkspecs/features/lex.prf:
/usr/share/qt4/mkspecs/features/include_source_dir.prf:
/usr/lib/x86_64-linux-gnu/libQtGui.prl:
/usr/lib/x86_64-linux-gnu/libQtCore.prl:
qmake:  FORCE
	@$(QMAKE) -o Makefile matav3.pro

dist: 
	@$(CHK_DIR_EXISTS) .tmp/matav31.0.0 || $(MKDIR) .tmp/matav31.0.0 
	$(COPY_FILE) --parents $(SOURCES) $(DIST) .tmp/matav31.0.0/ && $(COPY_FILE) --parents mainwindow.h src/constants.hpp src/CPPHeaders.hpp src/deftypes.hpp src/DIG.h src/DIGEndian.h src/DIGFile.h src/dsf_for_tree_bgnd_gen.hpp src/dsf_for_tree_gen.hpp src/edgytree.hpp src/historytree.hpp src/historytree_bgnd.hpp src/image3d.hpp src/matrix3d.hpp src/perbin.hpp src/treeViewer/Edge.h src/treeViewer/Node.h src/treeViewer/TreeManager.h src/treeViewer/TreeScene.h src/treeViewer/TreeView.h .tmp/matav31.0.0/ && $(COPY_FILE) --parents application.qrc .tmp/matav31.0.0/ && $(COPY_FILE) --parents main.cpp mainwindow.cpp src/treeViewer/Edge.cpp src/treeViewer/Node.cpp src/treeViewer/TreeManager.cpp src/treeViewer/TreeScene.cpp src/treeViewer/TreeView.cpp .tmp/matav31.0.0/ && $(COPY_FILE) --parents mainwindow.ui .tmp/matav31.0.0/ && (cd `dirname .tmp/matav31.0.0` && $(TAR) matav31.0.0.tar matav31.0.0 && $(COMPRESS) matav31.0.0.tar) && $(MOVE) `dirname .tmp/matav31.0.0`/matav31.0.0.tar.gz . && $(DEL_FILE) -r .tmp/matav31.0.0


clean:compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


####### Sub-libraries

distclean: clean
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


check: first

mocclean: compiler_moc_header_clean compiler_moc_source_clean

mocables: compiler_moc_header_make_all compiler_moc_source_make_all

compiler_moc_header_make_all: moc_mainwindow.cpp moc_TreeManager.cpp moc_TreeScene.cpp moc_TreeView.cpp
compiler_moc_header_clean:
	-$(DEL_FILE) moc_mainwindow.cpp moc_TreeManager.cpp moc_TreeScene.cpp moc_TreeView.cpp
moc_mainwindow.cpp: ui_mainwindow.h \
		src/CPPHeaders.hpp \
		src/constants.hpp \
		src/deftypes.hpp \
		src/DIGFile.h \
		src/matrix3d.hpp \
		src/edgytree.hpp \
		src/historytree.hpp \
		src/historytree_bgnd.hpp \
		src/perbin.hpp \
		src/dsf_for_tree_gen.hpp \
		src/dsf_for_tree_bgnd_gen.hpp \
		src/image3d.hpp \
		src/treeViewer/TreeManager.h \
		src/treeViewer/TreeView.h \
		src/treeViewer/Node.h \
		src/treeViewer/TreeScene.h \
		mainwindow.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) mainwindow.h -o moc_mainwindow.cpp

moc_TreeManager.cpp: src/CPPHeaders.hpp \
		src/constants.hpp \
		src/deftypes.hpp \
		src/edgytree.hpp \
		src/historytree.hpp \
		src/historytree_bgnd.hpp \
		src/perbin.hpp \
		src/treeViewer/TreeManager.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) src/treeViewer/TreeManager.h -o moc_TreeManager.cpp

moc_TreeScene.cpp: src/treeViewer/TreeScene.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) src/treeViewer/TreeScene.h -o moc_TreeScene.cpp

moc_TreeView.cpp: src/treeViewer/TreeView.h
	/usr/bin/moc-qt4 $(DEFINES) $(INCPATH) src/treeViewer/TreeView.h -o moc_TreeView.cpp

compiler_rcc_make_all: qrc_application.cpp
compiler_rcc_clean:
	-$(DEL_FILE) qrc_application.cpp
qrc_application.cpp: application.qrc \
		img/open.ico \
		img/zoomin.png \
		img/deleteselected.png \
		img/arrow.png \
		img/zoomout.png \
		img/addvalue.png \
		img/hand.png
	/usr/bin/rcc -name application application.qrc -o qrc_application.cpp

compiler_image_collection_make_all: qmake_image_collection.cpp
compiler_image_collection_clean:
	-$(DEL_FILE) qmake_image_collection.cpp
compiler_moc_source_make_all:
compiler_moc_source_clean:
compiler_uic_make_all: ui_mainwindow.h
compiler_uic_clean:
	-$(DEL_FILE) ui_mainwindow.h
ui_mainwindow.h: mainwindow.ui
	/usr/bin/uic-qt4 mainwindow.ui -o ui_mainwindow.h

compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: compiler_moc_header_clean compiler_rcc_clean compiler_uic_clean 

####### Compile

main.o: main.cpp mainwindow.h \
		ui_mainwindow.h \
		src/CPPHeaders.hpp \
		src/constants.hpp \
		src/deftypes.hpp \
		src/DIGFile.h \
		src/matrix3d.hpp \
		src/edgytree.hpp \
		src/historytree.hpp \
		src/historytree_bgnd.hpp \
		src/perbin.hpp \
		src/dsf_for_tree_gen.hpp \
		src/dsf_for_tree_bgnd_gen.hpp \
		src/image3d.hpp \
		src/treeViewer/TreeManager.h \
		src/treeViewer/TreeView.h \
		src/treeViewer/Node.h \
		src/treeViewer/TreeScene.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

mainwindow.o: mainwindow.cpp mainwindow.h \
		ui_mainwindow.h \
		src/CPPHeaders.hpp \
		src/constants.hpp \
		src/deftypes.hpp \
		src/DIGFile.h \
		src/matrix3d.hpp \
		src/edgytree.hpp \
		src/historytree.hpp \
		src/historytree_bgnd.hpp \
		src/perbin.hpp \
		src/dsf_for_tree_gen.hpp \
		src/dsf_for_tree_bgnd_gen.hpp \
		src/image3d.hpp \
		src/treeViewer/TreeManager.h \
		src/treeViewer/TreeView.h \
		src/treeViewer/Node.h \
		src/treeViewer/TreeScene.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o mainwindow.o mainwindow.cpp

Edge.o: src/treeViewer/Edge.cpp src/treeViewer/Edge.h \
		src/treeViewer/Node.h \
		src/CPPHeaders.hpp \
		src/constants.hpp \
		src/deftypes.hpp \
		src/edgytree.hpp \
		src/historytree.hpp \
		src/historytree_bgnd.hpp \
		src/perbin.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Edge.o src/treeViewer/Edge.cpp

Node.o: src/treeViewer/Node.cpp src/treeViewer/Node.h \
		src/CPPHeaders.hpp \
		src/constants.hpp \
		src/deftypes.hpp \
		src/edgytree.hpp \
		src/historytree.hpp \
		src/historytree_bgnd.hpp \
		src/perbin.hpp \
		src/treeViewer/Edge.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o Node.o src/treeViewer/Node.cpp

TreeManager.o: src/treeViewer/TreeManager.cpp src/treeViewer/TreeManager.h \
		src/CPPHeaders.hpp \
		src/constants.hpp \
		src/deftypes.hpp \
		src/edgytree.hpp \
		src/historytree.hpp \
		src/historytree_bgnd.hpp \
		src/perbin.hpp \
		src/treeViewer/TreeScene.h \
		src/treeViewer/Node.h \
		src/treeViewer/Edge.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o TreeManager.o src/treeViewer/TreeManager.cpp

TreeScene.o: src/treeViewer/TreeScene.cpp src/treeViewer/TreeScene.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o TreeScene.o src/treeViewer/TreeScene.cpp

TreeView.o: src/treeViewer/TreeView.cpp src/treeViewer/TreeView.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o TreeView.o src/treeViewer/TreeView.cpp

moc_mainwindow.o: moc_mainwindow.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_mainwindow.o moc_mainwindow.cpp

moc_TreeManager.o: moc_TreeManager.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_TreeManager.o moc_TreeManager.cpp

moc_TreeScene.o: moc_TreeScene.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_TreeScene.o moc_TreeScene.cpp

moc_TreeView.o: moc_TreeView.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o moc_TreeView.o moc_TreeView.cpp

qrc_application.o: qrc_application.cpp 
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o qrc_application.o qrc_application.cpp

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:    GL