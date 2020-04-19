FROM ubuntu:14.04

### Install prerequistes
RUN apt update --fix-missing
RUN apt update && apt upgrade -y && apt install -y build-essential cmake gfortran git flex bison \
	wget curl python python-sip-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev hdf5-tools linux-headers-generic \
	libvtk5-dev libvtk5-qt4-dev python-vtk tcl-vtk lvtk-dev lvtk-tools vtkdata libhdf5-dev


### Install PETSC
RUN wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.5.4.tar.gz
RUN tar -xvf petsc-3.5.4.tar.gz
WORKDIR petsc-3.5.4
RUN ./configure --with-mpi=0 \
  		--with-debugging=0 --with-shared-libraries=0 --with-x=0 --with-pic=1 \
  		--download-f2cblaslapack=1 \
  		--download-metis=1 \
  		--download-superlu=1 \
  		--COPTFLAGS="-O2" --CXXOPTFLAGS="-O2" --FOPTFLAGS="-O2"
RUN make


### Install CGNS
WORKDIR /
RUN curl -O -J -L https://sourceforge.net/projects/cgns/files/cgnslib_2.5/Release%205/cgnslib_2.5-5.tar.gz/download
RUN tar -xvf cgnslib_2.5-5.tar.gz
WORKDIR cgnslib_2.5
RUN chmod +x configure
RUN ./configure --enable-shared
RUN make && make install


### Install Genius
WORKDIR /
RUN git clone https://github.com/cogenda/Genius-TCAD-Open.git
WORKDIR Genius-TCAD-Open

RUN echo "--- wscript.orig	2020-04-18 19:02:15.583526950 +0000" >> lib.patch
RUN echo "+++ wscript	2020-04-18 19:06:15.895709943 +0000" > lib.patch
RUN echo "@@ -805,7 +805,7 @@" >> lib.patch
RUN echo "   def config_vtk():" >> lib.patch
RUN echo "     libs = '''vtkRendering vtkGraphics vtkImaging vtkIO vtkFiltering" >> lib.patch
RUN echo "-              vtkCommon vtksys vtkDICOMParser vtkpng vtktiff vtkzlib" >> lib.patch
RUN echo "+              vtkCommon vtksys vtkDICOMParser" >> lib.patch
RUN echo "-              vtkjpeg vtkexpat vtkftgl vtkfreetype'''.split()" >> lib.patch
RUN echo "+              vtkmetaio'''.split()" >> lib.patch
RUN echo "" >> lib.patch
RUN echo "     if platform=='Linux':" >> lib.patch
RUN echo "       libs.append('pthread')" >> lib.patch
RUN patch wscript < lib.patch

RUN ./waf --prefix=$PWD \
	--with-petsc-dir=/petsc-3.5.4 \
	--with-petsc-arch=arch-linux2-c-opt \
	--with-vtk-ver=vtk-5.8 \
	--with-vtk-dir=/usr \
	configure build install

# set the enviroment enviroments
ENV GENIUS_DIR /Genius-TCAD-Open
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH/:$GENIUS_DIR/lib:/usr/local/lib/
