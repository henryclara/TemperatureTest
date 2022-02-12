module load intel/18.0.4
module load intelmpi/2018.4.274

export ELMER_HOME="/home/mpim/m300792/ElmerIce/elmerice/Elmer_devel_04-29-20"
export ELMER_SOLVER_HOME="$ELMER_HOME/bin"
export PATH=/home/mpim/m300792/ElmerIce/elmerice/Elmer_devel_04-29-20/bin:$PATH
export \
PATH=/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/include:/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-intel14/lib/libnetcdf.so:/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/lib/libnetcdff.so:/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/lib:$PATH
export LD_LIBRARY_PATH=/home/mpim/m300792/ElmerIce/elmerice/Elmer_devel_04-29-20/share/elmersolver/lib/:$LD_LIBRARY_PATH

export \
LD_LIBRARY_PATH=/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/include:/sw/rhel6-x64/netcdf/netcdf_c-4.3.2-intel14/lib/libnetcdf.so:/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/lib/libnetcdff.so:/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14/lib:$LD_LIBRARY_PATH

export PATH="/home/mpim/m300792/ElmerIce/gmsh-2.12.0-source/build/:$PATH"
