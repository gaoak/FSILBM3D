**FSILBM3D**: Fluid-structure interaction solver based on Lattice Boltzmann method

## Table of contents

  * [Installation](#installation)
  * [Run solver](#run-solver)
  * [Input file description](#input-file-description)
  * [Simulation control parameters](#simulation-control-parameters)
  * [Output file description](#Output-file-description)
  * [License](#license)

## Installation
### Linux
Clone FSIBLM3D to a local machine,
```
git clone git@github.com:gaoak/FSILBM3D.git
cd FSILBM3D
```
To compile FSILBM3D in linux system, run
```
make
```
The executable file is **FSILBM3D**

The default compiler is intel fortran. To use a GNU gfortran compiler **gcc**, you need to change the first line in the Makefile
```
CMP = intel# intel,gcc
```
to
```
CMP = gcc# intel,gcc
```

## Run solver
### prepare work direcotry and input files
Copy fluid mesh (FluidMesh.dat in examples), body mesh (Plate.dat in examples) files to the work directory.

Create data folders in the work directory by running bash script (cleanfiles.sh)
```
for folder in DatTemp DatBodyIB DatBodySpan DatInfo DatOthe DatBody DatFlow
do
  echo $folder 
  if [[ -d $folder ]]; then
    rm $folder/*
  else
    mkdir $folder
  fi
done
```

### On local linux machine
Run the solve using bash script
```
nohup /home/gaoak/LBM/FSILBM3D/FSILBM3D &>log
```
### On compute node
submit the PBS script (2 threads for example)
```
#!/bin/sh 
#PBS -N LBM2D
#PBS -e stderr.txt 
#PBS -o stdout.txt
#PBS -l nodes=1:ppn=2
#PBS -l walltime=1:59:00
#PBS -q CPU_Small
PBS_WDIR=`pwd`/code/FSILBM3D
PBS_ENAME=FSILBM3D
#########################################
echo Working directory is $PBS_O_WORKDIR
echo test is $PBS_O_HOST
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`
NPROCS=`wc -l < $PBS_NODEFILE`
echo This job has allocated $NPROCS nodes
cd $PBS_O_WORKDIR
$PBS_WDIR/$PBS_ENAME
```

## Input file description

- **inFlow.dat**     Simulation control parameters
- **FluidMesh.dat**  Fluid mesh (use CRLF)
- **BMeshNum.dat**   Body mesh info (Removed)
- **Plate.dat**      Body mesh
- **cleanfiles.sh**  Bash script to create/clean output directories

## Simulation control parameters
- **Calculate control parameters**
1. *npsize* : The core number used in the simulaion
2. *isRelease* : Determining the content of the output log
    - 0 : Output more detail (body states) 
    - 1 : Brief output to save memory
3. *isConCmpt* : Determining new simulation or continue simulaion
    - 0 : Calculate from the beginning
    - 1 : Calculate from last output (DatTemp/conwr.dat)
4. *iCollidModel* : Determines the LBM model used in simualtion
    - 1 : *LBGK*
    - 2 : *MRT-LBGK*
5. *isKB* : Determining the kind of parameters of the flexible bodies
    - 0 : *KB, KS*
    - 1 : *EmR, TcR*

- **Time control**
1. *timeSimTotl*  : Dimensionless total simulation time (time/Tref <= timeSimTotl)
2. *timeOutTemp*  : Dimensionless time interval for restart file output (for continue simulation)
3. *timeOutFlow*  : Dimensionless time interval for flow field output
4. *timeOutBody*  : Dimensionless time interval for body output
5. *timeOutInfo*  : Dimensionless time interval for DatInfo output (force, velocity *et. al.*)
6. *timeOutBegin* : Dimensionless time to start flow field and body output
7. *timeOutEnd*   : Dimensionless time to end flow field and body output (output only in timeOutBegin<= time/Tref <= timeOutEnd)
8. *dt* : For uniform grid dt=dx=dy=dz. For non-uniform grid $(dt<min(dx, dy,dz))$

- **Reference values**
1. *Lref* : The reference length, whihc is the mimimum chord length of all bodies
3. *Uref* : The reference velocity, which is determined by *RefVelocity*
4. *RefVelocity* : Determining the definition of reference velocity
    - 0  : X-incoming flow velocity $U_\infty$;
    - 1  : Y-incoming flow velocity $V_\infty$;
    - 2  : Incoming flow velocity magnitude $\sqrt{U_\infty^2 + V_\infty^2}$
    - 3  : Maximum moving wall velocity
    - 4  : Velocity amplitude of Concussion velocity
    - 10 : Flapping frequency velocity $L f$
    - 11 : Maximum plunging velocity $2\pi f a$
    - 12 : Twice maximum plunging velocity used by Park et al. (2017) PoF $2\pi f a * 2$
    - *else* : The input value $(Uref)$ in parameter file inflow.dat
5. *RefTime* : Determining the definition of reference time
    - 0  : Caculated by referece length and reference velocity $(Tref = Lref / Uref)$
    - 1  : Caculated by the maximum frequency of the bodies $(Tref = 1 / max(frequency))$
    - *else*  : The input value $(Tref)$ in parameter file inflow.dat
5. *Tref=Lref/Uref* : The reference time, when Tref < 0 is defined as $(Lref/Uref)$, otherwise, it equals its input value
6. *Freq, St* : The flapping frequence for flexible bodies;
7. *Frod*: Inverse square of the Froude number $Frod = g L / U_{ref}^2$. Determining the gravity force exerted on the body

- **Initial conditions and boundary conditions**
1. *uIn* : The incoming Velocity ($U_\infty$, $V_\infty$, $W_\infty$), determined by the boundary kind
2. *boundaryConditions* : Boundary conditions on four boundaries: $(xmin, xmax, ymin, ymax, zmin, zmax)$ 
    - 101 : Symmetric boundary
    - 200 : Fixed wall
    - 201 : Moving wall, only for the top and bottom boundaries $(ymin, ymax)$
    - 300 : Dirichlet boundary (DirecletUP)
    - 301 : Dirichlet boundary (DirecletUU)
    - 302 : Neumann boundary (Advection1)
    - 303 : Neumann boundary (Advection2)
    - 304 : Periodic boundary (When calculating infinite bodies, the fluid grid should be one grid less than the solid grid)
3. *VelocityKind* : Especially for Direclet velocity conditions on the left and right boundaries
    - 0, Uniform or uniform shear flow, the boundary velocity is 
    $[uIn(1) = uuuIn(1) + 0 * shearRateIn(1) + y * shearRateIn(2) + z * shearRateIn(3)]$;
    - 2, Vibration flow, the boundary velocity is 
    $[uIn(1) + VelocityAmp * dcos(2 * pi * VelocityFreq * time), uIn(2)]$
4. *shearRateIn* : Parameters for uniform shear flow
5. *VelocityAmp,VelocityFreq* : Parameters for vibration flow
6. *MovingKind* : Especially for moving wall conditions $(201)$ on the bottom wall (1) and upper wall (2)
    - 0, Passive moving wall, the velocity is $uIn(1)$; 
    - 1, Couette moving wall, the velocity is $movingVel$; 
    - 2, Stokes moving wall, the velocity is $MovingVel * dcos(2 * pi * MovingFreq * time)$;
8. *movingVel,movingFreq* : Parameters for moving wall boundary
9. *AmplInitDist, FreqInitDist*
    - Parameters for initial velocity disturbance, the velocity is 
    $[Uin(1) + AmplInitDist(1) * dsin(2.0 * pi * FreqInitDist(1)),$
     $Uin(2) + AmplInitDist(2) * dsin(2.0 * pi * FreqInitDist(2)),$
     $Uin(3) + AmplInitDist(3) * dsin(2.0 * pi * FreqInitDist(3))]$;
10. *VolumeForceIn,VolumeForceAmp,VolumeForceFreq,VolumeForcePhi* : Parameters for Volume Force
    According to the NS equation : $[rho * Du/Dt = rho * f_x - dp/dx]$
    $[VolumeForceIn(1) + VolumeForceAmp * dsin(2 * pi * VolumeForceFreq * time + VolumeForcePhi), VolumeForceIn(2), VolumeForceIn(3)]$;

- **Moving grid**
1. *isMoveGrid* : The computation domain moves with first body if this equals 1
2. *iMoveDimX,Y,Z* : The grid moves space dimension, if this equals 1
3. *isMoveOutputX,Y,Z* : Outputting the relative flow grid and body, i.e. in the moving frame of reference, if this equals 1

- **Data probe set**
1. *numSampFlow* : The number of data detection points in flow field
2. *isFluidOutput* : Only when it equals 1 the detected fluid data will be output
3. *SampFlowPint* : The coordinates of detection points
4. *numSampBody* : The number of data detection points in bodies
5. *isBodyOutput * : Only when it equals 1 the detected body data will be output
6. *SampBodyNode* : The serial number of each detection points in the specific body

- **Body parameters input**
1. *dspan* : The virtual grid length of the bodies along the span
2. *theta* : The angle between the body span direction and the x-axis (have not yet implemented)
3. *Nspan* : The number of virtual grids along the span of the bodies
4. *FishKind* : Numbers of the solid bodies type
5. *nFish* : Total number of the solid bodies
6. *Fishnum* : The number of the solid bodies in the specific type
7. *FEmeshName* : The name of the body mesh file  in the specific type
8. *iBodyModel* : Choose for the body models
    - 1, rigid body;
    - 2, flexible body;
9. *isMotionGiven* : Degrees of freedom in six directions
10. *denR* : Density ratio, $rho_b * h / rho_f * L$, $h$ is plate thickness
11. *psR* : Poisson ratio
12. *KB, KS* : Dimensionless tension rigidity and bending rigidity
13. *waittingTime* : The dimensionaless time stayed at the peak and trough in flapping period $(t/T)$
14. *XYZointial* : The initial position of the first point of the bodies
15. *dXYZo* : If there are more than one bodies in a type, *dXYZo* determines the interval between front and rear solids
16. *XYZAmpl, XYZPhi* : Parameters for body flapping
    - $XYZ = XYZAmpl * dcos(2.0 * pi * Freq * time + XYZPhi)$
17. *AoAo, AoAAmpl, AoAPhi* : Parameters for body rotation
    - $Theta = AoAo + AoAAmpl * dcos(2.0 * pi * Freq * time + AoAPhi)$

## Output file description
- **Output files**
1. *DatBody*     : Folder of body results
2. *DatFlow*     : Floder of flow field results
3. *DatInfo*     : Floder of force and power *et. al.* results
4. *DatTemp*     : Floder of continuation document
5. *DatBodyIB*   ：Folder of body results for immersed boundary method when $Palpha \gt 0$
6. *DatBodySpan* ：Folder of spanwise-extension body results
7. *Check.dat*   : Parameter record file

- **Files description**
1. *SampBodyAngular.plt* 
    - *Hy*: *y* Coordinate of the first point $(y1)$
    - *Ty*: *y* Coordinate of the last point $(y2)$
    - *Ty-TH* : Height difference between the first and last points $(y2 - y1)$
    - *AoA* :  The deflection angle of the body $(y2-y1)/(x2-x1)$;
2. *SampBodyMean.plt* : The average information $(x, y, u, v, ax, ay)$ of the bodies 
3. *SampBodyBegin.plt* : The same information of the first point of the bodies 
4. *SampBodyEnd.plt* : The same information of the last point of the bodies 
5. *SampFlowPint.plt* : The information of the detected points in the flow field
6. *SampBodyNode.plt* : The information of the detected points in the bodies
7. *MaxValBody.plt* : The max velocity in all bodies points
8. *ForceDirect.plt* :The forces exerted on the bodies
9. *Energy.plt* : The energy of the bodies
10. *Power.plt* : The input power of the bodies
11. *Converg.plt* : The convergence of the simulation
12. *MaMax.plt* : The max Mach number in flow field
13. *Area.plt* : The surface area of the bodies

## License
Licese is owned by X-Y Lu research group. No distribution is allowed without a writing permission from Prof. Lu.
