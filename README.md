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
- **Plate.dat**      Body mesh
- **cleanfiles.sh**  Bash script to create/clean output directories

## Simulation control parameters

- **Parallel**

  1. *npsize* : The core number used in the simulaion

- **Communication**

  1. *npairs*：Number of parent-child relationships (0 represents only one fluid block)
  2.  *fatherId*: The ID of the parent block
  3. *sonId*: The ID of the son block

- **FlowCondition**

  1. *isConCmpt* : Determining new simulation or continue simulaion

     + 0 : Calculate from the beginning
   + 1: Calculate from last output
  
  2. *numsubstep*: Number of sub-steps for solid time-stepping solution

  3. *timeSimTotal*: Dimensionless total simulation time (time/Tref <= timeSimTotl)

  4. *timeContiDelta*: Dimensionless time interval for restart file output (for continue simulation)

  5. *timeWriteBegin*: Dimensionless time to start flow field and body output

  6. *timeWriteEnd*: Dimensionless time to end flow field and body output (output only in timeOutBegin<= time/Tref <= timeOutEnd)

  7. *timeWriteFlow*: Dimensionless time interval for flow field output

  8. *timeWriteBody*  : Dimensionless time interval for body output

  9. *timeWriteInfo*  : Dimensionless time interval for DatInfo output (force, velocity *et. al.*)

  10. *Re*: Dimensionless Reynolds Number 

  11. *denIn*: Fluid density           

  12. *uvwIn*: The incoming Velocity 
    $$
      U_\infty, V_\infty, W_\infty
      $$
       (determined by the boundary kind)
  
  13. *shearRateIn*: Parameters for uniform shear flow 

  14. *VolumeForceIn,VolumeForceAmp,VolumeForceFreq,VolumeForcePhi* : Parameters for Volume Force

       According to the NS equation : 
    $$
      rho * Du/Dt = rho * f_x - dp/dx
      $$
  
      $$
    $[VolumeForceIn(1) + VolumeForceAmp * dsin(2 * pi * VolumeForceFreq * time + VolumeForcePhi), VolumeForceIn(2), VolumeForceIn(3)]$
      $$
  
  15. *TrefType*: Determining the definition of reference time

      * 0  : Caculated by referece length and reference velocity 
      $$
        (Tref = Lref / Uref)
        $$
  
      * 1  : Caculated by the maximum frequency of the bodies 
      $$
        Tref = 1 / max(frequency)
        $$
  
      * *else*  : The input value Tref in parameter file inflow.dat

  16. *UrefType*:  Determining the definition of reference velocity

      - 0  : X-incoming flow velocity
      $$
        U_\infty
        $$
  
      - 1  : Y-incoming flow velocity 
      $$
        V_\infty
        $$
  
      - 2  : Incoming flow velocity magnitude
      $$
        \sqrt{U_\infty^2 + V_\infty^2}
        $$
  
      - 3  : Maximum moving wall velocity

      - 4  : Velocity amplitude of Concussion velocity

      - 10 : Flapping frequency velocity 
      $$
        L f
        $$
  
      - 11 : Maximum plunging velocity 
      $$
        2\pi f a
        $$
  
      - 12 : Twice maximum plunging velocity used by Park et al. (2017) PoF 
      $$
        2\pi f a * 2
        $$
  
      - *else* : The input value Uref in parameter file inflow.dat

  17. *ntolLBM*: Maximum number of iterations for the LBM method

  18. *dtolLBM*: Error tolerance for the LBM method

- **FluidBlocks**

  1. *nblock*: Number of fluid grid partitions
  2. *ID*: The ID of  fluid grid
  3. *iCollideModel*:  Determines the LBM model used in simualtion  
     + 1 : *LBGK* ：Single Relaxation Time
     + 2 : *MRT-LBGK*：Multiple Relaxation Time Lattice Boltzmann Method
  4. *offsetOutput*: The computation domain moves with first body if this equals 1
  5. *isoutput*: Outputting the relative flow grid and body, i.e. in the moving frame of reference, if this equals 1
  6.  *xDim,yDim,zDim*: Number of nodes in the x, y, and z directions of the fluid block
  7. *dh*: For uniform grid dh=dx=dy=dz
  8. *xmin,ymin,zmin*: Starting position of fluid block
  9. *boundaryConditions(1:6) (xmin,xmax,ymin,ymax,zmin,zmax)*: Boundary conditions on six boundaries
     + 101 : Symmetric boundary
     + 103：Periodic Boundary
     + 200 : Fixed wall
     + 201 : Moving wall, only for the top and bottom boundaries $(ymin, ymax)$
     + 300 : Dirichlet boundary (DirecletUP)
     + 301 : Dirichlet boundary (DirecletUU)
     + 302 : Neumann boundary (Advection1)
     + 303 : Neumann boundary (Advection2)
     + 304 : Periodic boundary (When calculating infinite bodies, the fluid grid should be one grid less than the solid grid, i.e. fluid grid containes one end whereas solid grid contains both ends)
  10. *params*: Pending Parameters

- **SolidBody**

  1. *IBPenaltyalpha*:Velocity correction parameter of the penalty function in the IBM method

  2. *alphaf*:Parameter for correcting torsion in the mass matrix

  3. *NewmarkGamma,NewmarkBeta*:Parameters of the Newmark method

  4. *dampK,dampM*: Stiffness damping K, mass damping M

  5. *dtolFEM,ntolFEM*:Maximum number of iterations and error tolerance for the FEM method

  6. *nFish* : Total number of the solid bodies

  7. *nfishGroup*: Numbers of the solid bodies type

  8. *carrierFluidId*: Usually it's 2

  9. *isKB*:  Determining the kind of parameters of the flexible bodies 

     * 0 : *KB, KS*(Bending Stiffness, Stretching Stiffness)
     * 1 : *EmR, TcR* (Elastic Modulus Ratio, Characteristic Time Ratio)

  10. *fishnum  (fishGroup）*: The number of the solid bodies in the specific type

  11. *numX,numY,numZ*: Number of arrangements of multiple bodies in the x, y, and z directions   !!!!!!!!!!!

  12. *FEmeshName(iFish)*: The name of the body mesh file  in the specific type

  13. *iBodyModel(iFish)*: Choose for the body models

      + 1: rigid body
      + 2: flexible body

  14. *iBodyType(iFish)*: Type of virtual object

  15. *isMotionGiven*: Degrees of freedom in six directions

  16. *denR（iFish）* : Density ratio, 
      $$
      rho_b * h / rho_f * L
      $$
      L is plate thickness

  17. *psR(iFish)*: Poisson ratio

  18. *Freq（iFish)* : The flapping frequence for flexible bodies

  19. *firstXYZ*: The initial position of the first point of the bodies

  20. *deltaXYZ*: If there are more than one bodies in a type, *deltaXYZ* determines the interval between front and rear solids

  21. *XYZAmpl, XYZPhi* : Parameters for body flapping

       - $$
         XYZ = XYZAmpl * dcos(2.0 * pi * Freq * time + XYZPhi)
         $$

  22. *AoAo, AoAAmpl, AoAPhi* : Parameters for body rotation

       - $$
         Theta = AoAo + AoAAmpl * dcos(2.0 * pi * Freq * time + AoAPhi)
         $$

- **ProbingFluid**
  1. *fluidProbingNum*: The number of data detection points in flow field
  2. *fluidProbingCoords* : The coordinates of detection points
- **ProbingSolid**
  1. *solidProbingNum*:  The number of data detection points in bodies
  2. *solidProbingNode*: The serial number of each detection points in the specific body

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
