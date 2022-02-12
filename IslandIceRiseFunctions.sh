#!/bin/bash
###############################################################################
### This function creates a mesh file for elmer for a rectangular domain ######
###############################################################################
CreateMeshGeoRectangleInit() {
				rm Mesh.geo
				HalfWidth=$(($1/2))
				DomLength=$2
				Res=$3
				echo $HalfWidth
cat > Mesh.geo << EOF
Mesh.Algorithm=5;
Point(1)={0,$HalfWidth,0.0,$Res};
Point(2)={0,-$HalfWidth,0.0,$Res};
Point(3)={$DomLength,-$HalfWidth,0.0,$Res};
Point(4)={$DomLength,$HalfWidth,0.0,$Res};
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
Line Loop(5)={1,2,3,4};
Physical Line(6)={1};
Physical Line(7)={2};
Physical Line(8)={3};
Physical Line(9)={4};
Plane Surface(10)={5};
Physical Surface(11)={10};
EOF
}
###############################################################################
### This function creates a makefile for the intial remeshing #################
###############################################################################
CreateMakefileRemeshInit() {
rm Makefile
SifFileName=$1
echo $SifFileName
NumProcs=1
NumProcx=1
NumProcy=1
cat > Makefile << EOF
# Makefile for Elmer
# ----------------------------------------
# Use external Geometry to create mesh
# Calculate Depth and Height for Paraview

EXECUTABLES = src/DistanceSolverRD 


NumProcs=$NumProcs
NumProcx=$NumProcx
NumProcy=$NumProcy

InputSif=$SifFileName



.SUFFIXES: .f90

all: clean ini grid submit

grid:
	gmsh Mesh.geo -1 -2
	ElmerGrid 14 2 Mesh.msh  -autoclean

submit: ini

	mpirun -n $NumProcs ElmerSolver_mpi

compile:
	elmerf90 src/AgeSolverRD.f90 -o src/AgeSolverRD  
	elmerf90 src/DistanceSolveRD.f90 -o src/DistanceSolveRD
	elmerf90 src/BedrockBump.f90 -o src/BedrockBump
	elmerf90 src/GroundedMaskSolver.f90 -o src/GroundedMaskSolver
	elmerf90 src/SeaLevel.f90 -o src/SeaLevel
	elmerf90 src/USF_BMB.f90 -o src/USF_BMB;  elmerf90 -o src/USF_Contact src/USF_Contact.f90 src/USF_Sliding.f90

clean:
	rm -fr results/*

ini:
	echo $SifFileName > ELMERSOLVER_STARTINFO

.f90:
	elmerf90  -o $@ $<
.c:
	gcc  -o $@ $< -lm
EOF
}
###############################################################################
### This function creates an Elmer sif file for the intial remeshing ##########
###############################################################################
CreateElmerSifRemeshInit() {
SifFileName=$1
BumpApl=$2
Sigmax=$3
Sigmay=$4
xCenter=$5
yCenter=$6
ZsInit=$7
ZbInit=$8
BedInit=$9
RiseIntoShelf=${10}
BackgroundRes=${11}
MaxRes=${12}
RefineRadius=${13}
Incline=${14}
InclineShelf=${15}
Alpha=${16}
G=${17}
A=${18}
rho=${19}

rm $SifFileName
cat > "$SifFileName" << EOF
!------------------------------------------------------------------
! Isotropic mesh adaptation:
!------------------------------------------------------------------
\$Bedrock=$BedInit

Header
Mesh DB "." "Mesh"
End

Constants
Bump Amplitude = Real $BumpApl
Sigmax = Real $Sigmax
Sigmay = Real $Sigmay
x0 = Real $xCenter
y0 = Real $yCenter
ZbInit = Real $ZbInit
ZsInit = Real $ZsInit
RiseIntoShelf = Real $RiseIntoShelf
Incline = Real $Incline
InclineShelf = Real $InclineShelf
Alpha = Real $Alpha
G = Real $G
A = Real $A
rho = Real $rho
End


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Simulation
Coordinate System  = Cartesian 2D
Simulation Type = Steady

Steady State Min Iterations = 6
Steady State Max Iterations = 6

max output level = 30
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Body 1
Equation = 1
Body Force = 1
Initial Condition = 1
End
Initial Condition 1
BedInit = Real \$Bedrock
Bedrock = Variable BedInit
Real Procedure "src/BedrockBump" "BedrockBump"
Zb = Variable Bedrock
Real Procedure "src/BedrockBump" "ZbAdj"
Zs = Variable Zb
Real Procedure "src/BedrockBump" "ZsAdj"
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Body Force 1
Distance = Real 0.0
Distance Condition = Variable GroundedMask
Real MATC "tx"
ElementSize = Variable Distance
REAL MATC "if (tx(0)<=$RefineRadius) {$MaxRes} else {$BackgroundRes}"
END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Solver 1
!Exec Solver = Never
!Exec Solver = "Before All"
!Exec Solver = "Before Simulation"
Equation = GroundedMask
Procedure = "src/GroundedMaskSolver" "GroundedMaskSolver"
Variable = GroundedMask
Variable DOFs = 1
End

Solver 2
!Exec Solver = Never
!Exec Solver = Before All
Equation = "SolveDistance"

Procedure = "src/DistanceSolveRD" "DistanceSolver1"
!Procedure = "Executables/DistanceSolveRD" "DistanceSolverInit"
Variable = Distance

H scale = real 2
Distance Pseudo DT = Real 100
! Nonlinear System Relaxation Factor = 0.25

Nonlinear System Max Iterations = 50
Nonlinear System Convergence Tolerance = 1.0e-5

! Linear System Solver = Direct
! Linear System Direct Method = UMFPack
Linear System Solver = "Iterative"
Linear System Iterative Method = "BiCGStab"
Linear System Max Iterations = 300
Linear System Convergence Tolerance = 1.0E-09
Linear System Abort Not Converged = False
Linear System Preconditioning = "ILU1"
Linear System Residual Output = 1
Steady State Convergence Tolerance = 1.0e-4

Dummy Distance Computation = Logical False

End
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Solver 3
Equation = "Initialise fn"
Procedure = "ElmerIceSolvers" "UpdateExport"

Exported Variable 1 = Zs
Exported Variable 2 = Zb
Exported Variable 3 = Bedrock
Exported Variable 4 = Distance
Exported Variable 5 = ElementSize
Exported Variable 6 = BedInit
Exported Variable 7 = RiseIntoShelf
!Exported Variable 8 = Sea level
End

Solver 4
!! mandatory else Model % Mesh % Changed reset to .FALSE. in coupled simulations
!Exec Solver = Never
!Exec Solver = after timestep

Equation = "MMG"
Procedure = "ElmerIce_MeshAdapt2D" "MMG2DSolver"

Output file name = "square_iso"
Metric Variable Name = String "ElementSize"
hausd = Real 5000.0 !Hausdorff parameter (controls the refinement near boundaries)
hgrad = Real 1.3  !gradation value (controls the ratio between two adjacent edges)
End

Solver 5
Exec Solver = After Saving
Equation = "result output"
Procedure = "ResultOutputSolve" "ResultOutputSolver"
Save Geometry Ids = Logical True ! add this line if you want to access boundaries in Paraview
Output File Name = File "EkstroemLGM"
Output Format = String vtu
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Equation 1
Active Solvers(5) = 1 2 3 4 5 
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Boundary Condition 1
Target Boundaries = 1
End
Boundary Condition 2
Target Boundaries = 2
End
Boundary Condition 3
Target Boundaries = 3
End
Boundary Condition 4
Target Boundaries = 4
End
EOF
}

###############################################################################
### This function creates a submit script for the initial remeshing ###########
###############################################################################
CreateSLURMSubmitScriptRemeshInit(){
				rm Submit.sh
 				Path2Dir=$1
				InitDir=$2
				ForwardDir=$3
				SuperComputer=$4
				Email=$5
cat > Submit.sh << EOF
#!/bin/bash
#SBATCH -o $Path2Dir/SLURM_job.%j.%N.out
#SBATCH -e $Path2Dir/SLURM_job.%j.%N.err
#SBATCH -D $Path2Dir
#SBATCH -J RemeshInit
#SBATCH --get-user-env
EOF
if [ ! -z "$Email" ]; then
cat >> Submit.sh << EOF
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$Email
EOF
fi
cat >> Submit.sh << EOF
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --time=00:15:00
EOF
if [ "$SuperComputer" = "Mistral" ]; then
cat >> Submit.sh << EOF
#SBATCH --partition=compute,compute2
#SBATCH --account=bm1164
EOF
else
cat >> Submit.sh << EOF
#SBATCH --partition=test
#SBATCH --account=pn56pe
#source /etc/profile.d/modules.sh
module load slurm_setup
EOF
fi

cat >> Submit.sh << EOF
#=================================================================================================================
set -e
source ModulesPlusPaths${SuperComputer}.sh
echo Here comes the Nodelist:
echo \$SLURM_JOB_NODELIST

echo Here comes the partition the job runs in:
echo \$SLURM_JOB_PARTITION
cd \$SLURM_SUBMIT_DIR

make compile
make ini
make grid
make submit

cp -r square_iso_N6/mesh*  $InitDir/Mesh/
#cp -r square_iso_N6/mesh*  $ForwardDir/Mesh/

cd $InitDir
sbatch Submit.sh 0

EOF
}
###############################################################################
### This function copies all necessary files to the source directory ###########
###############################################################################
CopyFiles2Source(){
				Path2SourceDir=$1
  			cp -r $Path2SourceDir/ $PWD
}
###############################################################################
### This function creates a makefile for the intial forward run ###############
###############################################################################
CreateMakefileForwardInit() {
				rm Makefile
				SifFileName=$1
				echo $SifFileName
				NumProcx=$2
				NumProcy=$3
				NumProcz=$4
				NumProcs=$5
cat > Makefile << EOF
# Makefile for Elmer
# ----------------------------------------
# Use external Geometry to create mesh
# Calculate Depth and Height for Paraview

EXECUTABLES = src/DistanceSolverRD 


NumProcs=$NumProcs
NumProcx=$NumProcx
NumProcy=$NumProcy
NumProcz=$NumProcz

InputSif=$SifFileName



.SUFFIXES: .f90

all: clean ini grid submit

grid:
	#ElmerGrid 2 2 Mesh -partition $NumProcx $NumProcy $NumProcz -autoclean
	ElmerGrid 2 2 Mesh -metis $NumProcs 4  -autoclean

submit: ini

	mpirun -n $NumProcs ElmerSolver_mpi

compile:  
	elmerf90 src/DistanceSolveRD.f90 -o src/DistanceSolveRD
	elmerf90 src/BedrockBump.f90 -o src/BedrockBump
	elmerf90 src/GroundedMaskSolver.f90 -o src/GroundedMaskSolver
	elmerf90 src/SeaLevel.f90 -o src/SeaLevel
	elmerf90 src/USF_BMB.f90 -o src/USF_BMB; elmerf90 -o src/USF_Contact src/USF_Contact.f90 src/USF_Sliding.f90

clean:
	rm -fr results/*

ini:
	echo $SifFileName > ELMERSOLVER_STARTINFO

.f90:
	elmerf90  -o $@ $<
.c:
	gcc  -o $@ $< -lm
EOF
}
###############################################################################
### This function creates a submit script for the initial forward sim #########
###############################################################################
CreateSLURMSubmitScriptForwardInit(){
				rm Submit.sh
 				Path2Dir=$1
				Nodes=$2
				PrcoNo=$3
				Queue=$4
				RunTime=$5
				JobNameInit=$6
				SimLength=$7
				TotalSimLength=$8
				CopyYes=$9
				ForwardDir=${10}
				SuperComputer=${11}
				Email=${12}
				echo $RunTime

cat > Submit.sh << EOF
#!/bin/bash
#SBATCH -o $Path2Dir/SLURM_job.%j.%N.out
#SBATCH -e $Path2Dir/SLURM_job.%j.%N.err
#SBATCH -D $Path2Dir
#SBATCH -J $JobNameInit
#SBATCH --get-user-env
EOF
if [ ! -z "$Email" ]; then
cat >> Submit.sh << EOF
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$Email
EOF
fi
if [ "$SuperComputer" = "Mistral" ]; then
cat >> Submit.sh << EOF
#SBATCH --account=bm1164
#SBATCH --ntasks=$ProcNo
#SBATCH --time=$RunTime
EOF
else
cat >> Submit.sh << EOF
#SBATCH --account=pn56pe
#SBATCH --nodes=$Nodes
#SBATCH --ntasks=$ProcNo
#SBATCH --ntasks-per-node=48
#SBATCH --export=NONE
#SBATCH --time=$RunTime
#source /etc/profile.d/modules.sh
module load slurm_setup
EOF
fi
cat >> Submit.sh << EOF
#SBATCH --partition=$Queue
#=================================================================================================================
set -e
echo Here comes the Nodelist:
echo \$SLURM_JOB_NODELIST

echo Here comes the partition the job runs in:
echo \$SLURM_JOB_PARTITION
cd \$SLURM_SUBMIT_DIR

source ModulesPlusPaths${SuperComputer}.sh

cp \$ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so \
src/MyFreeSurfaceSolver.so
echo \$ELMER_HOME
echo \$ELMER_SOLVER_HOME

YearCounter=\$1
YearCounterFormatted=\$(printf %06d \$YearCounter)
sed -i "s/FORMAT/\${YearCounterFormatted}/g" Init.sif
echo YearCounter is: $YearCounter
make compile
make ini
make grid
EOF
if [ "$SuperComputer" = "Mistral" ]; then
cat >> Submit.sh << EOF
srun -l --export=ALL --cpu_bind=cores --distribution=block:cyclic -n $ProcNo ElmerSolver_mpi
EOF
else
cat >> Submit.sh << EOF
make submit
EOF
fi
cat >> Submit.sh << EOF
if [ "\${YearCounter}" -lt "${TotalSimLength}" ]; then
	if [ "${CopyYes}" -eq 1 ]; then
					cp Mesh/*result* $ForwardDir/Mesh/
					cp Mesh/mesh* $ForwardDir/Mesh/
					cd $ForwardDir
					sbatch Submit.sh \$YearCounter
	fi
fi
EOF
}
###############################################################################
### This function creates an Elmer sif file for the intial remeshing ##########
###############################################################################
CreateElmerSifForwardInit() {
				SifFileName=$1
                BumpApl=$2
				Sigmax=$3
				Sigmay=$4
				xCenter=$5
				yCenter=$6
				ZsInit=$7
				ZbInit=$8
				BedInit=$9
				RiseIntoShelf=${10}
				MeshLayers=${11}
				Rhoi=${12}
				Rhow=${13}
				RateFactor=${14}
				GlenExponent=${15}
				BasalFrictionCoeff=${16}
				SMB=${17}
				FluxInit=${18}
				OutputName=${19}
				IceTemp=${20}
				Incline=${21}
				InclineShelf=${22}
				Alpha=${23}
				G=${24}
				A=${25}
				rho=${26}
                SlidStr=${27}
                SlidExp=${28}
				rm ${SifFileName}
cat > "${SifFileName}" << EOF
!!--------------------------------------------------------!!
!  Island ice rise setup for initial step
!!--------------------------------------------------------!!

check keywords warn
!
! working units are MPa, a, m
!
\$yearinsec = 365.25*24*60*60
\$rhoi = ${Rhoi}/(1.0e6*yearinsec^2)
\$rhow = ${Rhow}/(1.0e6*yearinsec^2)
\$A = ${RateFactor}*yearinsec*1.0e18
\$n = ${GlenExponent}
\$eta = 1.0/(2.0*A)^(1.0/n)
\$gravity = -9.8*yearinsec^2
\$C = ${BasalFrictionCoeff}/(1.0e6*yearinsec^(1.0/${SlidExp}))

\$Bedrock=$BedInit
\$Incline=$Incline
\$InclineShelf=$InclineShelf

$ function BedTopo(xcoord) import Bedrock,InclineShelf {\\
				_BedTopo = Bedrock + 1/100.0 * xcoord * tan(InclineShelf*pi/180)\\
}

Header
  Mesh DB "." "Mesh"
End

Constants
  Water Density = Real \$rhow
  Gas Constant = Real 8.314 !Joule/mol x  K
	Bump Amplitude = Real $BumpApl
	Sigmax = Real $Sigmax
	Sigmay = Real $Sigmay
  x0 = Real $xCenter
	y0 = Real $yCenter
	ZbInit = Real $ZbInit
	ZsInit = Real $ZsInit
	RiseIntoShelf = Real $RiseIntoShelf
	Incline = Real \$Incline
	Alpha = Real $Alpha
	G = Real $G
	A = Real $A
	rho = Real $rho
  ! For SeaSpring/SeaPressure
End

!---------------------------------------------------
!---------------- SIMULATION -----------------------
!---------------------------------------------------

Simulation
  Coordinate System  = Cartesian 3D
  Simulation Type = transient
  Extruded Mesh Levels = Integer $MeshLayers

  Timestepping Method = "bdf"
  BDF Order = 1
  Timestep Intervals = 1
  Output Intervals = 1
  Timestep Sizes = 0.01

  Initialize Dirichlet Conditions = Logical False
  Steady State Max Iterations = 6
  Steady State Min Iterations = 6

	!Specify name of result file. Used for restarts!!
  Output File = "${OutputName}FORMAT.result"
  max output level = 30
End

!---------------------------------------------------
!---------------- BODIES ---------------------------
!---------------------------------------------------

! the ice
Body 1
  Name = "ice"
  Equation = 1
  Body Force = 1
  Material = 1
  Initial Condition = 1
End

! The upper surface
Body 2
  Name= "top free surface"
  Equation = 2
  Material = 1
  Body Force = 2
  Initial Condition = 2
End

! the lower surface
Body 3
  Name= "free surface sea/ice-shelf"
  Equation = 3
  Material = 1
  Body Force = 3
  Initial Condition = 3
End

!---------------------------------------------------
!---------------- INITIAL CONDITIONS ---------------
!---------------------------------------------------

!! for ice
Initial Condition 1
  Pressure = Real 0.0
  Velocity 1 = Real 0.0
  Velocity 2 = Real 0.0
  Velocity 3 = Real 0.0
	BedInit = Variable Coordinate 1
		Real MATC "BedTopo(tx)"
	FluxInit = Real $FluxInit
	BedBump = Variable BedInit
			Real Procedure "src/BedrockBump" "BedrockBump"
End

!! for top free surface
Initial Condition 2
	Zs = Variable BedBump
			Real Procedure "src/BedrockBump" "ZsAdj"
End

!! for free surface sea/ice-shelf
Initial Condition 3
	Zb = Variable BedBump
			Real Procedure "src/BedrockBump" "ZbAdj"
End

!---------------------------------------------------
!---------------- BODY FORCES ----------------------
!---------------------------------------------------

Body Force 1
  Flow BodyForce 1 = Real 0.0
  Flow BodyForce 2 = Real 0.0
  Flow BodyForce 3 = Real \$gravity
End

!! accumulation flux in m/year
Body Force 2
   Zs Accumulation Flux 1 = Real 0.0e0
   Zs Accumulation Flux 2 = Real 0.0e0 !m/a
   Zs Accumulation Flux 3 = Real $SMB
End

!! no melting/accretion under ice/shelf
Body Force 3
  !Zb Accumulation = Real $BMB
	Zb Accumulation = Variable Time
			Real Procedure "src/USF_BMB" "GetBMB"
End

!---------------------------------------------------
!---------------- MATERIALS ------------------------
!---------------------------------------------------

!! ice material properties in MPa - m - a system
Material 1
  Viscosity Model = String "power law"
  Density = Real \$rhoi
  Viscosity = Real \$eta
  Viscosity Exponent = Real \$1.0/n
  Critical Shear Rate = Real 1.0e-15

  Sea level = Real 0.0

  Glen Enhancement Factor = Real 1.0
! the temperature to switch between the
! two regimes in the flow law
  Limit Temperature = Real -10.0
! In case there is no temperature variable
  Constant Temperature = Real $IceTemp

  Min Zs = Variable "Bottom Zb"
    Real MATC "tx + 10.0"
  Max Zs = Real 1.0e6

  !! Bed condition
  Min Zb = Equals BedBump
  Max Zb = Real 1.0e6
  !Cauchy = Logical True
End

!---------------------------------------------------
!---------------- SOLVERS --------------------------
!---------------------------------------------------
!! Initialisation of the Grounded Mask
Solver 1
  !Exec Solver = Never
  Exec Solver = Before All
  Equation = "MapCoordinate"
  Procedure = "StructuredMeshMapper" "StructuredMeshMapper"

  Active Coordinate = Integer 3
  Mesh Velocity Variable = String "dSdt"
  Mesh Update Variable = String "dS"
  Mesh Velocity First Zero = Logical True

  Top Surface Variable Name = String "Zs"
  Bottom Surface Variable Name = String "Zb"

  Displacement Mode = Logical False
  Correct Surface = Logical True
  Minimum Height = Real 1.0
End

Solver 2
  !Exec Solver = Never
  Equation = GroundedMaskIni
  Procedure = "ElmerIceSolvers" "GroundedSolver"
  Variable = GroundedMask
  Variable DOFs = 1

  Toler = Real 1.0e-3
  Bedrock Variable = String "BedBump"
End


Solver 3
  !Exec Solver = Never
  Equation = "NormalVector"
  Procedure = "ElmerIceSolvers" "ComputeNormalSolver"
  Variable = String "Normal Vector"
  Variable DOFs = 3

  ComputeAll = Logical False
  Optimize Bandwidth = Logical False
End

Solver 4
  !Exec Solver = Never
  Equation = Fw
  Procedure = "ElmerIceSolvers" "GetHydrostaticLoads"
  Variable = Fw[Fwater:3]
  Variable DOFs = 3
End

Solver 5
  !Exec Solver = Never
  Equation = "Navier-Stokes"
	Optimize Bandwidth = Logical True
  Linear System Solver = Direct
  Linear System Direct Method = "Mumps"
	Mumps percentage increase working space = Integer 1600

  Nonlinear System Max Iterations = 50
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Newton After Iterations = 50
  Nonlinear System Newton After Tolerance = 1.0e-05
  Nonlinear System Relaxation Factor = 1.00
  Nonlinear System Reset Newton = Logical True

  Steady State Convergence Tolerance = Real 5.0e-5

  Stabilization Method = String Stabilized!Bubbles

  Exported Variable 1 = Flow Solution Loads[Stress Vector:3 CEQ Residual:1]
  Calculate Loads = Logical True

  Exported Variable 2 = -dofs 1 "dSdt"
  Exported Variable 3 = -dofs 1 "dS"
  Exported Variable 4 = -dofs 1 "BedInit"
  Exported Variable 5 = -dofs 1 "RiseIntoShelf"
  Exported Variable 6 = -dofs 1 "BedBump"
  Exported Variable 7 = -dofs 1 "FluxInit"
  Flow Model = String "Stokes"
End

Solver 6
  Equation = String "StressSolver"
  Procedure =  File "ElmerIceSolvers" "ComputeDevStress"
  ! this is just a dummy, hence no output is needed
  !-----------------------------------------------------------------------
  Variable = -nooutput "Sij"
  Variable DOFs = 1
  ! the name of the variable containing the flow solution (U,V,W,Pressure)
  !-----------------------------------------------------------------------
  Flow Solver Name = String "Flow Solution"
  ! no default value anymore for "Stress Variable Name"
  Stress Variable Name = String "Stress"
  !-----------------------------------------------------------------------
  Exported Variable 1 = "Stress" ! [Sxx, Syy, Szz, Sxy] in 2D
                                 ! [Sxx, Syy, Szz, Sxy, Syz, Szx] in 3D
  Exported Variable 1 DOFs = 6   ! 4 in 2D, 6 in 3D
  Linear System Solver = "Iterative"
  Linear System Iterative Method = "BiCGStab"
  Linear System Max Iterations = 300
  Linear System Convergence Tolerance = 1.0E-09
  Linear System Abort Not Converged = True
  Linear System Preconditioning = "ILU0"
  Linear System Residual Output = 1
End

Solver 7
  !Exec Solver = Never
  Exec Solver = Before All
  Equation = "HeightDepth"
  Procedure = "StructuredProjectToPlane" "StructuredProjectToPlane"
  Active Coordinate = Integer 3
  Dot Product Tolerance = Real 1.0e-3

  Operator 1 = Depth
  Operator 2 = Height
  Variable 3 = Zb
  Operator 3 = Bottom
End

Solver 8
  !Exec Solver = Never
   Equation = "SolveDistance"

   Procedure = "src/DistanceSolveRD" "DistanceSolver1"
   Variable = Distance

   H scale = real 2
   Distance Pseudo DT = Real 100
! Nonlinear System Relaxation Factor = 0.25

   Nonlinear System Max Iterations = 50
   Nonlinear System Convergence Tolerance = 1.0e-5

 ! Linear System Solver = Direct
 ! Linear System Direct Method = UMFPack
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 300
   Linear System Convergence Tolerance = 1.0E-09
   Linear System Abort Not Converged = False
   Linear System Preconditioning = "ILU1"
   Linear System Residual Output = 1
   Steady State Convergence Tolerance = 1.0e-4

   Dummy Distance Computation = Logical False

End

Solver 9
  !Exec Solver = Never
  Equation = "Free Surface Top"
  Procedure =  "./src/MyFreeSurfaceSolver" "FreeSurfaceSolver"
  !Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Variable = "Zs"
  Variable DOFs =  1
  Exported Variable 1 = "Zs Residual"
  Exported Variable 1 DOFs = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  !Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 10
  !Exec Solver = Never
  Equation = "Free Surface Sea/Shelf"
  Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Variable = "Zb"
  Variable DOFS =  1
  Exported Variable 1 = "Zb Residual"
  Exported Variable 1 DOFs = 1

  Nonlinear Update Exported Variables = Logical True

  Exported Variable 2 = "Zb Accumulation "
  Exported Variable 2 DOFS = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 11
  Exec Solver = After Saving
  Equation = "result output"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Save Geometry Ids = Logical True ! add this line if you want to access boundaries in Paraview
  Output File Name = File "${OutputName}FORMAT"
  Output Format = String vtu
End


!---------------------------------------------------
!---------------- EQUATIONS ------------------------
!---------------------------------------------------

Equation 1
  Active Solvers (6) = 1 3 5 6 7 11 
End

Equation 2
  Active Solvers(1) = 9
  Flow Solution Name = String "Flow Solution"
  Convection = String Computed
End

Equation 3
  Active Solvers(4) = 2 4 8 10 
  Flow Solution Name = String "Flow Solution"
  Convection = String Computed
End

!---------------------------------------------------
!---------------- BOUNDARY CONDITIONS --------------
!---------------------------------------------------

!! Back
Boundary Condition 1
  Name = "back"
  Target Boundaries = 1
	Velocity 1 = Variable FluxInit
		Real Procedure "src/BedrockBump" "IceFluxAtBack"
	Velocity 2 = Real 0.0
End

Boundary Condition 2
  Name = "Looking Downhill Left"
  Target Boundaries = 2

  Velocity 2 = Real 0.0
End

!! BC Lateral Ice-Shelf (air or sea contact)
Boundary Condition 3
  Name = "front"
  Target Boundaries = 3


  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Compute Sea Pressure = Logical True
  ComputeNormal = Logical False

End

Boundary Condition 4
  Name = "Looking Downhill Right"
  Target Boundaries = 4

  Velocity 2 = Real 0.0
End

Boundary Condition 5
  Name = "bottom"
  Target Boundaries = 5
  Body Id = 3

  Normal-Tangential Velocity = Logical True
  Flow Force BC = Logical True

!
! Condition where the bed is stuck
!
  Zb = Equals BedBump
  Zb Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Bedrock conditions
!
  Slip Coefficient 2 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"
  Slip Coefficient 3 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"

  Sliding Law = String "${SlidStr}"
EOF
if [ "$SlidStr" = "Weertman" ]; then
cat >> "${SifFileName}" << EOF
  Weertman Friction Coefficient = Real \$C
  Weertman Exponent = Real \$(1.0/${SlidExp})
  Weertman Linear Velocity = Real 0.001
EOF
elif [ "$SlidStr" = "Coulomb" ]; then
cat >> "${SifFileName}" << EOF
  Friction Law Sliding Coefficient = Real 4.1613e5 ! was C
  Friction Law Post-Peak Exponent = Real 1.0
  Friction Law Maximum Value = Real 0.5
  Friction Law PowerLaw Exponent = Real ${SlidExp}
  Friction Law Linear Velocity = Real 0.001
EOF
fi
cat >> "${SifFileName}" << EOF
  ! Options are 'Last Grounded' (default), 'First Floating' or 'Discontinuous'
   Grounding Line Definition = String "First Floating"
   Test Contact Tolerance = real 1.0e-3
   Non Detachment Inland Distance = Real 5000.0 ! distance from the GL where nodes

  Velocity 1 = Real 0.0
  Velocity 1 Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Shelf conditions
!
  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Slip Coefficient 1 = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaSpring"

  ComputeNormal Condition = Variable GroundedMask
    Real MATC "tx + 0.5"

  Compute Sea Pressure = Logical True
  Compute Sea Spring = Logical True

  Distance = Real 0.0
  Distance Condition = Variable GroundedMask
    Real MATC "tx"
End

!! BC Lateral Ice-Shelf (air or sea contact)
!! BC  Free surface Top
Boundary Condition 6
  Name = "top"
  Target Boundaries = 6
  Body Id = 2
  ComputeNormal = Logical False
End
EOF
}
###############################################################################
### This function creates a submit script for the initial remeshing ###########
###############################################################################
CreateSLURMSubmitScriptForward(){
				rm Submit.sh
 				Path2Dir=$1
				Nodes=$2
				ProcNo=$3
				Queue=$4
				RunTime=$5
				JobName=$6
				SimLength=$7
				TotalSimLength=$8
				ForwardSif=${9}
				SuperComputer=${10}
				Email=${11}
				echo $RunTime

cat > Submit.sh << EOF
#!/bin/bash
#SBATCH -o $Path2Dir/SLURM_job.%j.%N.out
#SBATCH -e $Path2Dir/SLURM_job.%j.%N.err
#SBATCH -D $Path2Dir
#SBATCH -J $JobName
#SBATCH --get-user-env
EOF
if [ ! -z "$Email" ]; then
cat >> Submit.sh << EOF
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$Email
EOF
fi
if [ "$SuperComputer" = "Mistral" ]; then
cat >> Submit.sh << EOF
#SBATCH --account=bm1164
#SBATCH --ntasks=$ProcNo
#SBATCH --time=$RunTime
EOF
else
cat >> Submit.sh << EOF
#SBATCH --account=pn56pe
#SBATCH --nodes=$Nodes
#SBATCH --ntasks=$ProcNo
#SBATCH --ntasks-per-node=48
#SBATCH --export=NONE
#SBATCH --time=$RunTime
#source /etc/profile.d/modules.sh
module load slurm_setup
EOF
fi
cat >> Submit.sh << EOF
#SBATCH --partition=$Queue
#=================================================================================================================
set -e
echo Here comes the Nodelist:
echo \$SLURM_JOB_NODELIST

echo Here comes the partition the job runs in:
echo \$SLURM_JOB_PARTITION
cd \$SLURM_SUBMIT_DIR

source ModulesPlusPaths${SuperComputer}.sh

cp \$ELMER_HOME/share/elmersolver/lib/FreeSurfaceSolver.so \
src/MyFreeSurfaceSolver.so
YearCounter=\$1
echo YearCounter is: \$YearCounter
if [ "\${YearCounter}" -lt "${TotalSimLength}" ]; then
	YearCounterFormatted=\$(printf %06d \$YearCounter)
	YearCounter=\$((\$YearCounter+$SimLength))
	YearCounterFormattedNew=\$(printf %06d \$YearCounter)
	cp ${ForwardSif}.bak ${ForwardSif}
	sed -i "s/START/\${YearCounterFormatted}/g" $ForwardSif
	sed -i "s/END/\${YearCounterFormattedNew}/g" $ForwardSif
	echo \$YearCounter
	make compile
	make ini
	make grid
EOF
if [ "$SuperComputer" = "Mistral" ]; then
cat >> Submit.sh << EOF
srun -l --export=ALL --cpu_bind=cores --distribution=block:cyclic -n $ProcNo ElmerSolver_mpi
EOF
else
cat >> Submit.sh << EOF
make submit
EOF
fi
cat >> Submit.sh << EOF
	mv Mesh Output\${YearCounterFormattedNew}
	mkdir -p Mesh
	cp Output\${YearCounterFormattedNew}/mesh* Mesh/
	cp Output\${YearCounterFormattedNew}/Forward\${YearCounterFormattedNew}*.result* \
		 Mesh
	sbatch Submit.sh \$YearCounter
fi
EOF
}

###############################################################################
### This function creates an Elmer sif file for the forward simulation ########
###############################################################################
CreateElmerSifForward() {
				SifFileName=$1
				MeshLayers=${2}
				Rhoi=${3}
				Rhow=${4}
				RateFactor=${5}
				GlenExponent=${6}
				BasalFrictionCoeff=${7}
				IceTemp=${8}
				SMB=${9}
				FluxInit=${10}
				OutputName=${11}
				RestartName=${12}
				NoOfTimeSteps=${13}
				OutputInterval=${14}
				TimeStepSize=${15}
				Alpha=${16}
				G=${17}
				A=${18}
				rho=${19}
                SlidStr=${20}
                SlidExp=${21}
				rm ${SifFileName}.bak
cat > "${SifFileName}.bak" << EOF
!!--------------------------------------------------------!!
!  Island ice rise setup for forwards simulation
!!--------------------------------------------------------!!

check keywords warn
!
! working units are MPa, a, m
!
\$yearinsec = 365.25*24*60*60
\$rhoi = ${Rhoi}/(1.0e6*yearinsec^2)
\$rhow = ${Rhow}/(1.0e6*yearinsec^2)
\$A = ${RateFactor}*yearinsec*1.0e18
\$n = ${GlenExponent}
\$eta = 1.0/(2.0*A)^(1.0/n)
\$gravity = -9.8*yearinsec^2
\$C = ${BasalFrictionCoeff}/(1.0e6*yearinsec^(1.0/${SlidExp}))


Header
  Mesh DB "." "Mesh"
End

Constants
  Water Density = Real \$rhow
  Gas Constant = Real 8.314 !Joule/mol x  K
	Alpha = Real $Alpha
	G = Real $G
	A = Real $A
	rho = Real $rho
  ! For SeaSpring/SeaPressure
End

!---------------------------------------------------
!---------------- SIMULATION -----------------------
!---------------------------------------------------

Simulation
  Coordinate System  = Cartesian 3D
  Simulation Type = transient
  Extruded Mesh Levels = Integer $MeshLayers

  Timestepping Method = "bdf"
  BDF Order = 1
  Timestep Intervals = $NoOfTimeSteps
  Output Intervals = $OutputInterval
  Timestep Sizes = $TimeStepSize

  Initialize Dirichlet Conditions = Logical False
  Steady State Max Iterations = 1
  Steady State Min Iterations = 1

	Restart File="${RestartName}START.result"
	Restart Before Initial Conditions = Logical True
	!Specify name of result file. Used for restarts!!
  Output File = "${OutputName}END.result"
  max output level = 30
End

!---------------------------------------------------
!---------------- BODIES ---------------------------
!---------------------------------------------------

! the ice
Body 1
  Name = "ice"
  Equation = 1
  Body Force = 1
  Material = 1
  Initial Condition = 1
End

! The upper surface
Body 2
  Name= "top free surface"
  Equation = 2
  Material = 1
  Body Force = 2
  Initial Condition = 2
End

! the lower surface
Body 3
  Name= "free surface sea/ice-shelf"
  Equation = 3
  Material = 1
  Body Force = 3
  Initial Condition = 3
End

!---------------------------------------------------
!---------------- INITIAL CONDITIONS ---------------
!---------------------------------------------------

!! for ice
Initial Condition 1
End

!! for top free surface
Initial Condition 2
End

!! for free surface sea/ice-shelf
Initial Condition 3
End

!---------------------------------------------------
!---------------- BODY FORCES ----------------------
!---------------------------------------------------

Body Force 1
  Flow BodyForce 1 = Real 0.0
  Flow BodyForce 2 = Real 0.0
  Flow BodyForce 3 = Real \$gravity
End

!! accumulation flux in m/year
Body Force 2
   Zs Accumulation Flux 1 = Real 0.0e0
   Zs Accumulation Flux 2 = Real 0.0e0 !m/a
   Zs Accumulation Flux 3 = Real $SMB
End

!! no melting/accretion under ice/shelf
Body Force 3
  !Zb Accumulation = Real $BMB
	Zb Accumulation = Variable Time
			Real Procedure "src/USF_BMB" "GetBMB"
End

!---------------------------------------------------
!---------------- MATERIALS ------------------------
!---------------------------------------------------

!! ice material properties in MPa - m - a system
Material 1
  Viscosity Model = String "power law"
  Density = Real \$rhoi
  Viscosity = Real \$eta
  Viscosity Exponent = Real \$1.0/n
  Critical Shear Rate = Real 1.0e-15
  !Sea level = Real 0.0
  Sea level = Variable Time
    Real Procedure "src/SeaLevel" "getSeaLevel"
  Glen Enhancement Factor = Real 1.0
! the temperature to switch between the
! two regimes in the flow law
  Limit Temperature = Real -10.0
! In case there is no temperature variable
  Constant Temperature = Real $IceTemp

  Min Zs = Variable "Bottom Zb"
    Real MATC "tx + 10.0"
  Max Zs = Real 1.0e6

  !! Bed condition
  Min Zb = Equals BedBump
  Max Zb = Real 1.0e6
  Cauchy = Logical True
End

!---------------------------------------------------
!---------------- SOLVERS --------------------------
!---------------------------------------------------
Solver 1
  Exec Solver = Before Simulation
  Equation = "MapCoordinateInit"
  Procedure = "StructuredMeshMapper" "StructuredMeshMapper"

  Active Coordinate = Integer 3
  Mesh Velocity Variable = String "dSdt"
  Mesh Update Variable = String "dS"
  Mesh Velocity First Zero = Logical True

  Top Surface Variable Name = String "Zs"
  Bottom Surface Variable Name = String "Zb"

  Displacement Mode = Logical False
  Correct Surface = Logical True
  Minimum Height = Real 1.0
End

Solver 2
  !Exec Solver = Never
  Equation = "NormalVector"
  Procedure = "ElmerIceSolvers" "ComputeNormalSolver"
  Variable = String "Normal Vector"
  Variable DOFs = 3

  ComputeAll = Logical False
  Optimize Bandwidth = Logical False
End

Solver 3
  !Exec Solver = Never
  Equation = Fw
  Procedure = "ElmerIceSolvers" "GetHydrostaticLoads"
  Variable = Fw[Fwater:3]
  Variable DOFs = 3
End

Solver 4
  !Exec Solver = Never
  Equation = "Navier-Stokes"
	Optimize Bandwidth = Logical True
  Linear System Solver = Direct
  Linear System Direct Method = "Mumps"
	Mumps percentage increase working space = Integer 1600

  Nonlinear System Max Iterations = 50
  Nonlinear System Convergence Tolerance  = 1.0e-6
  Nonlinear System Newton After Iterations = 50
  Nonlinear System Newton After Tolerance = 1.0e-05
  Nonlinear System Relaxation Factor = 1.00
  Nonlinear System Reset Newton = Logical True

  Steady State Convergence Tolerance = Real 5.0e-5

  Stabilization Method = String Stabilized!Bubbles

  Exported Variable 1 = Flow Solution Loads[Stress Vector:3 CEQ Residual:1]
  Calculate Loads = Logical True

  Exported Variable 2 = -dofs 1 "dSdt"
  Exported Variable 3 = -dofs 1 "dS"
  Exported Variable 4 = -dofs 1 "BedInit"
  Exported Variable 5 = -dofs 1 "RiseIntoShelf"
  Exported Variable 6 = -dofs 1 "BedBump"
  Exported Variable 7 = -dofs 1 "FluxInit"
  Flow Model = String "Stokes"
End

Solver 5
  Equation = String "StressSolver"
  Procedure =  File "ElmerIceSolvers" "ComputeDevStress"
  ! this is just a dummy, hence no output is needed
  Variable = -nooutput "Sij"
  Variable DOFs = 1
  ! the name of the variable containing the flow solution
  !(U,V,W,Pressure)
  Flow Solver Name = String "Flow Solution"
  ! no default value anymore for "Stress Variable Name"
  Stress Variable Name = String "Stress"
  !-----------------------------------------------------------------------
   Exported Variable 1 = "Stress"
   Exported Variable 1 DOFs = 6   ! 4 in 2D, 6 in 3D
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 300
   Linear System Convergence Tolerance = 1.0E-09
   Linear System Abort Not Converged = True
   Linear System Preconditioning = "ILU3"
   Linear System Residual Output = 1
End

Solver 6
  !Exec Solver = Never
  Equation = "HeightDepth"
  Procedure = "StructuredProjectToPlane" "StructuredProjectToPlane"
  Active Coordinate = Integer 3
  Dot Product Tolerance = Real 1.0e-3
  Operator 1 = Depth
  Operator 2 = Height
  Variable 3 = Zb
  Operator 3 = Bottom
End

Solver 7
  !Exec Solver = Never
   Equation = "SolveDistance"

   Procedure = "src/DistanceSolveRD" "DistanceSolver1"
   Variable = Distance

   H scale = real 2
   Distance Pseudo DT = Real 100
! Nonlinear System Relaxation Factor = 0.25

   Nonlinear System Max Iterations = 50
   Nonlinear System Convergence Tolerance = 1.0e-5

 ! Linear System Solver = Direct
 ! Linear System Direct Method = UMFPack
   Linear System Solver = "Iterative"
   Linear System Iterative Method = "BiCGStab"
   Linear System Max Iterations = 300
   Linear System Convergence Tolerance = 1.0E-09
   Linear System Abort Not Converged = False
   Linear System Preconditioning = "ILU1"
   Linear System Residual Output = 1
   Steady State Convergence Tolerance = 1.0e-4

   Dummy Distance Computation = Logical False

End

Solver 8
  !Exec Solver = Never
  Equation = "Free Surface Top"
  Procedure =  "./src/MyFreeSurfaceSolver" "FreeSurfaceSolver"
  !Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Variable = "Zs"
  Variable DOFs =  1
  Exported Variable 1 = "Zs Residual"
  Exported Variable 1 DOFs = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  !Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 9
  !Exec Solver = Never
  Equation = "Free Surface Sea/Shelf"
  Procedure =  "FreeSurfaceSolver" "FreeSurfaceSolver"
  Variable = "Zb"
  Variable DOFS =  1
  Exported Variable 1 = "Zb Residual"
  Exported Variable 1 DOFs = 1

  Nonlinear Update Exported Variables = Logical True

  Exported Variable 2 = "Zb Accumulation "
  Exported Variable 2 DOFS = 1

  !Before Linsolve = "EliminateDirichlet" "EliminateDirichlet"

  Linear System Solver = Iterative
  Linear System Direct Method = UMFPACK
  Linear System Max Iterations = 1500
  Linear System Iterative Method = BiCGStab
  Linear System Preconditioning = ILU0
  Linear System Convergence Tolerance = Real 1.0e-6
  Linear System Abort Not Converged = False
  Linear System Residual Output = 1

  Nonlinear System Max Iterations = 100
  Nonlinear System Convergence Tolerance  = 1.0e-5
  Nonlinear System Relaxation Factor = 1.00

  Steady State Convergence Tolerance = 1.0e-03

  Stabilization Method = Stabilized
  Apply Dirichlet = Logical True

  Relaxation Factor = Real 1.0
End

Solver 10
  Equation = "MapCoordinate"
  Procedure = "StructuredMeshMapper" "StructuredMeshMapper"

  Active Coordinate = Integer 3
  Mesh Velocity Variable = String "dSdt"
  Mesh Update Variable = String "dS"
  !Mesh Velocity First Zero = Logical True

  Top Surface Variable Name = String "Zs"
  Bottom Surface Variable Name = String "Zb"

  Displacement Mode = Logical False
  Correct Surface = Logical True
  Minimum Height = Real 1.0
End

Solver 11
  !Exec Solver = Never
  Equation = GroundedMask
  Procedure = "ElmerIceSolvers" "GroundedSolver"
  Variable = GroundedMask
  Variable DOFs = 1

  Toler = Real 1.0e-3
  Bedrock Variable = String "BedBump"
End

Solver 12
  Procedure = "SaveData" "SaveMaterials"
  Parameter 1 = String "Sea Level"
End

Solver 13
  Exec Solver = After Saving
  Equation = "result output"
  Procedure = "ResultOutputSolve" "ResultOutputSolver"
  Save Geometry Ids = Logical True ! add this line if you want to access boundaries in Paraview
  Output File Name = File "${OutputName}END"
  Output Format = String vtu
End

Solver 14
  Equation = "Flowdepth"
   Exec Solver = "Never"
   Procedure = File "ElmerIceSolvers" "FlowDepthSolver"
   Variable = String "Depth"
   Variable DOFs = 1
   Linear System Solver = "Direct"
   Linear System Direct Method = "UMFPACK"
   ! this sets the direction
   ! -1 is negative z-direction (upside down)
   ! +1 is positive (downside up)
   Gradient = Real -1.0E00
  ! switch that to True, if you want to have
  ! free surface gradients to be computed
  !------------------------------------
  Calc Free Surface = Logical True
  ! the name for the exported (if not existing) added variable
  ! the gradients will be stored in variables with the base
  ! name given and "Grad1" and (in 3 dimensions) "Grad2" added,
  ! so in our case "FreeSurfGrad1" and "FreeSurfGrad2"
  ! again, if those variables did not exist, they will be
  ! automatically created
  !-----------------------------------------------------------
  Freesurf Name = String "FreeSurf"
End

!---------------------------------------------------
!---------------- EQUATIONS ------------------------
!---------------------------------------------------

Equation 1
  Active Solvers (9) = 1 2 4 5 6 10 12 13 14
End

Equation 2
  Active Solvers(1) = 8
  Flow Solution Name = String "Flow Solution"
  Convection = String Computed
End

Equation 3
  Active Solvers(4) = 3 7 9 11
  Flow Solution Name = String "Flow Solution"
  Convection = String Computed
End

!---------------------------------------------------
!---------------- BOUNDARY CONDITIONS --------------
!---------------------------------------------------

!! Back
Boundary Condition 1
  Name = "back"
  Target Boundaries = 1
	Velocity 1 = Variable FluxInit
		Real Procedure "src/BedrockBump" "IceFluxAtBack"
	Velocity 2 = Real 0.0
End

Boundary Condition 2
  Name = "Looking Downhill Left"
  Target Boundaries = 2

  Velocity 2 = Real 0.0
End

!! BC Lateral Ice-Shelf (air or sea contact)
Boundary Condition 3
  Name = "front"
  Target Boundaries = 3


  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Compute Sea Pressure = Logical True
  ComputeNormal = Logical False

End

Boundary Condition 4
  Name = "Looking Downhill Right"
  Target Boundaries = 4

  Velocity 2 = Real 0.0
End

Boundary Condition 5
  Name = "bottom"
  Target Boundaries = 5
  Body Id = 3

  Normal-Tangential Velocity = Logical True
  Flow Force BC = Logical True

!
! Condition where the bed is stuck
!
  Zb = Equals BedBump
  Zb Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Bedrock conditions
!
  Slip Coefficient 2 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"
  Slip Coefficient 3 = Variable Coordinate 1
    Real Procedure "src/USF_Contact" "SlidCoef_Contact"
  Sliding Law = String "${SlidStr}"
EOF
if [ "$SlidStr" = "Weertman" ]; then
cat >> "${SifFileName}.bak" << EOF
  Weertman Friction Coefficient = Real \$C
  Weertman Exponent = Real \$(1.0/${SlidExp})
  Weertman Linear Velocity = Real 0.001
EOF
elif [ "$SlidStr" = "Coulomb" ]; then
cat >> "${SifFileName}.bak" << EOF
  Friction Law Sliding Coefficient = Real 4.1613e5 ! Was C
  Friction Law Post-Peak Exponent = Real 1.0
  Friction Law Maximum Value = Real 0.5
  Friction Law PowerLaw Exponent = Real ${SlidExp}
  Friction Law Linear Velocity = Real 0.001
EOF
fi
cat >> "${SifFileName}.bak" << EOF
  ! Options are 'Last Grounded' (default), 'First Floating' or 'Discontinuous'
    Grounding Line Definition = String "First Floating"
  Test Contact Tolerance = real 1.0e-3
  Non Detachment Inland Distance = Real 5000.0 ! distance from the GL where nodes

  Velocity 1 = Real 0.0
  Velocity 1 Condition = Variable GroundedMask
    Real MATC "tx + 0.5"
!
! Shelf conditions
!
  External Pressure = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaPressure"

  Slip Coefficient 1 = Variable Coordinate 3
     Real Procedure "ElmerIceUSF" "SeaSpring"

  ComputeNormal Condition = Variable GroundedMask
    Real MATC "tx + 0.5"

  Compute Sea Pressure = Logical True
  Compute Sea Spring = Logical True

  Distance = Real 0.0
  Distance Condition = Variable GroundedMask
    Real MATC "tx"
End

!! BC Lateral Ice-Shelf (air or sea contact)
!! BC  Free surface Top
Boundary Condition 6
  Name = "top"
  Target Boundaries = 6
  Body Id = 2
  ComputeNormal = Logical False
End
EOF
}
###############################################################################
### This function checks if the total CPU number specified is consistent with
### the specified CPU numbers in x,y, and z-direction ########################
###############################################################################
CheckCPUNumberConsistency() {
	ProcNo=$1
	Procx=$2
	Procy=$3
	Procz=$4
	ProcProd=$(( $Procx * $Procy * $Procz ))
	if [ "$ProcProd" -eq "$ProcNo" ]; then
		echo Good. Processor number match. Continuing...
	else
		echo Error. Processor numbers do not match. 
		echo Procx \* Procy \* Procz = $ProcProd 
		echo But ProcNo = $ProcNo
		echo Needs to be the same number. Aborting script ...
		exit 1
	fi
}

###############################################################################
### This function checks if the total CPU number specified is a multiple of the
### available CPUs at each node. It makes sure that no resources are wasted.
###############################################################################

CheckNodeNumberConsistency(){
	ProcNo=$1
	CPUsOnSingleNode=$2
	Result=$(( $ProcNo % $CPUsOnSingleNode ))
	if [ "$Result" -eq "0" ]; then
		echo Good. ProcNo is multiple of CPUsOnSinlgeNode. Continuing ...
	else
		echo Error. ProcNo is not a multiple of CPUsOnSingleNode.
		echo ProcNo = $ProcNo
		echo CPUsOnSingleNode = $CPUsOnSingleNode
		echo Please change this. Script aborting ... 
		exit 1
	fi
}
