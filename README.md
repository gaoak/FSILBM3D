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

The default compiler is **intel fortran**. To use a GNU gfortran compiler **gcc**, you need to change the first line in the Makefile

```
CMP = intel# intel,gcc
```

to

```
CMP = gcc# intel,gcc
```

## Run solver

### prepare work direcotry and input files

Copy input parameter file (**inFlow.dat** in examples), body mesh file (**plate.dat** in examples) and create or clean directories bash script (**cleanfiles.sh** in examples) to the work directory.

Create data folders in the work directory by running bash script (**cleanfiles.sh**)

```
for folder in DatBodySpan DatInfo DatBody DatFlow 
do
  if [[ -d $folder ]]; then
    rm -r $folder
    mkdir $folder
    echo Clear $folder
  else
    mkdir $folder
    echo Create $folder
  fi
done

for folder in DatContinue 
do
  if [[ ! -d $folder ]]; then
    mkdir $folder
    echo Create $folder
  fi
done

for file in check.dat log.* *.txt
do
  if [[ -e $file ]]; then
    rm $file
    echo Delete $file
  fi
done
```

### Running on local linux machine

Run the solve using bash script

```
nohup /your/local/path/FSILBM3D/FSILBM3D &>log
```

### Runing on the compute node

Submit the PBS script (2 threads for example)

```
#!/bin/sh 
#PBS -N LBM3D-Test
#PBS -e stderr.txt 
#PBS -o stdout.txt
#PBS -l nodes=1:ppn=2
#PBS -l walltime=1:59:00
#PBS -q CPU_Small
PBS_WDIR=`pwd`/FSILBM3D/FSILBM3D
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
- **plate.dat**      Body mesh file
- **cleanfiles.sh**  Bash script to create/clean output directories

## Simulation parameters description

(A line starting with # indicates a comment)

- **Parallel**

  1. *npsize* : The core number used in the simulaion

- **FlowCondition**

  1. *isConCmpt* : Determining new simulation or continue simulaion
     + 0 : Calculate from the beginning
     + 1 : Calculate from the continue file (in **DatContinue**)

  2. *numsubstep*: Number of sub-steps for solid time-stepping solution

  3. *timeSimTotal*: Dimensionless total simulation time (time/Tref <= timeSimTotl)

  4. *timeContiDelta*: Dimensionless time interval for writing continue files (continue** in **DatContinue**) which used for continue simulation

  5. *timeWriteBegin*: Dimensionless time to start fluid mesh (in DatFlow) and body mesh (in DatBody) writing, also is the start time for fluid averaging

  6. *timeWriteEnd*: Dimensionless time to end fluid mesh (in DatFlow) and body mesh (in DatBody) writing (timeOutBegin should less than timeOutEnd)

  7. *timeWriteFlow*: Dimensionless time interval for fluid mesh (in DatFlow) writing

  8. *timeWriteBody*  : Dimensionless time interval for body mesh (in DatBody) writing

  9. *timeWriteInfo*  : Dimensionless time interval for post procesing information (such as force, velocity *et. al.* in DatInfo)

  10. *Re*: Dimensionless Reynolds Number 

  11. *denIn*: Fluid density  (Usually it's 1 )  

  12. *uvwIn*: The given incoming Velocity $(U_\infty, V_\infty, W_\infty)$

      (determined by the velocity control parameters **velocityKind**)

  13. *shearRateIn*, *velocityKind* : the velocity control parameters
   * 0 : if *velocityKind* equals 0 , *shearRateIn* are the shear rate for incoming flow

     $U_\infty = U_{in} + 0 * shearRate(1) + y * shearRate(2) + z * shearRate(3)$

     $V_\infty = V_{in} + x * shearRate(1) + 0 * shearRate(2) + z * shearRate(3)$

     $W_\infty = W_{in} + x * shearRate(1) + y * shearRate(2) + 0 * shearRate(3)$

   * 1 : if *velocityKind* equals 2 , *shearRateIn* are the oscillatory parameter for incoming flow

     $U_\infty = U_{in} + shearRate(1) * cos(2 * pi * shearRate(2) * time + shearRate(3) * pi / 180)$

     $V_\infty = V_{in}$

     $W_\infty = W_{in}$

  14. *VolumeForceIn,VolumeForceAmp,VolumeForceFreq,VolumeForcePhi* : Parameters for Volume Force

      According to the NS equation : 

      $rho * Du/Dt = rho * f_x - dp/dx$

      The volume force is calculated as following:

      $volumeForce(1) = VolumeForceIn(1) + VolumeForceAmp * sin(2 * pi * VolumeForceFreq * time + VolumeForcePhi)$

      $volumeForce(2) = VolumeForceIn(2)$

      $volumeForce(3) = VolumeForceIn(3)$

  15. *TrefType*: Determining the definition of reference time

      * 0  : Caculated by referece length and reference velocity 

        $Tref = Lref / Uref$

      * 1  : Caculated by the maximum frequency of the bodies 

        $Tref = 1 / max(frequency)$

      * *else*  : The input value Tref in parameter file inflow.dat

  16. *UrefType*:  Determining the definition of reference velocity

      * 0  : X-incoming flow velocity : $U_\infty$

      * 1  : Y-incoming flow velocity : $V_\infty$

      * 2  : Z-incoming flow velocity : $W_\infty$

      * 3  : Incoming flow velocity magnitude : $\sqrt{U_\infty^2 + V_\infty^2}$

      * 4  : The velocity amplitude only for **velocityKind**  2 : $|shearRateIn(1)|$
      * 5  : Flapping frequency velocity : $L f$

      - 6  : Maximum plunging velocity : $2\pi f a$

      - 7  : Twice maximum plunging velocity used by Park et al. (2017) PoF : $2\pi f a * 2$

      - *else* : The input value Uref in parameter file inflow.dat

  17. *ntolLBM*: Maximum number of iterations for the LBM method

  18. *dtolLBM*: Error tolerance for the LBM method

  19. *interpolateScheme*: Interpolate Scheme for multi-block communication.

      - 1  : Linear interpolation.

      - 2  : Third-(boundary) and fourth-order(inner) mix interpolation.
  
- **FluidBlocks**

  1. *nblock*: Number of fluid grid partitions
  2. *ID*: The ID of  fluid grid
  3. *iCollideModel*:  Determines the LBM model used in simualtion  
     + 1 : *SRT* : Single Relaxation Time
     + 2 : *TRT* : Double Relaxation Time (only for single fluid block)
     + 3 : *MRT* : Multiple Relaxation Time (only for single fluid block)
     + 11 : *SMAG-SRT*: Single Relaxation Time With Smagorinsky Model
     + 12 : *Regular-SRT*: Regularised Single Relaxation Time
     + 13 : *ELBM-SRT*: Single Relaxation Time In ELBM 
     + 14 : *WALE-SRT*: WALE Single Relaxation Time
     + 15 : *Vremann-SRT*: Vremann Single Relaxation Time
  4. *offsetOutput*: The computation domain moves with first body if this equals 1
  5. *outputtype*: Determines output type of the relative flow grid and body
     + 0 : no fluid output
     + 1 : output fluid instantaneous velocity
     + 2 : output fluid average velocity
     + 3 : output fluid instantaneous and average velocity
  6. *xDim,yDim,zDim*: Number of nodes in the x, y, and z directions of the fluid block
  7. *dh*: For uniform grid dh=dx=dy=dz
  8. *xmin,ymin,zmin*: Starting position of the fluid block
  9. *boundaryConditions(1:6) (xmin,xmax,ymin,ymax,zmin,zmax)*: Boundary conditions parameters on six directions
     + 0:  Fluid boundary conditions
     + 101: Dirichlet boundary condition (velocity equal to a specified value)
     + 102: Dirichlet boundary condition (The value of velocity at the boundary is not a fixed constant, but a non-uniform distribution of function values)
     + 103: First order extrapolation boundary conditions
     + 104: Second-order extrapolation boundary conditions
     + 201: Fixed wall(full way bounce back)
     + 202: Moving wall, only for the top and bottom boundaries $(ymin, ymax)$
     + 203: Fixed wall(half way bounce back, good choice for turbulent in LES compute)
     + 301: Periodic Boundary(For the son block, the xDim/yDim/zDim should have have one more layer than the father block)
     + 302: Symmetric boundary
  10. *params*: Pending Parameters

- **SolidBody**

  1. *IBPenaltyalpha* : Velocity correction parameter of the penalty function in the IBM method

  2. *alphaf* : Parameter for correcting torsion in the mass matrix

  3. *NewmarkGamma, NewmarkBeta*: Parameters of the Newmark method

  4. *dampK, dampM* : Stiffness damping K, mass damping M

  5. *dtolFEM, ntolFEM* : Maximum number of iterations and error tolerance for the FEM method

  6. *nFish* : Total number of the solid bodies

  7. *nfishGroup* : Numbers of the solid bodies type

  8. *isKB* :  Determining the kind of parameters of the flexible bodies 

     * 0 : *KB, KS*(Bending Stiffness, Stretching Stiffness)
     * 1 : *EmR, TcR* (Elastic Modulus Ratio, Characteristic Time Ratio)

  9. *fishnum  (fishGroup)* : The number of the solid bodies in the specific type

  10. *numX,numY,numZ* : Number of arrangements of multiple bodies in the x, y, and z directions  

  11. *FEmeshName(iFish)* : The name of the body mesh file  in the specific type

  12. *iBodyModel(iFish)* : Choose for the body models

      + 1: rigid body
      + 2: flexible body

  13. *iBodyType(iFish)* : Type of virtual object

  14. *isMotionGiven* : Degrees of freedom in six directions

  15. *denR(iFish)* : Density ratio $(rho_b * h / rho_f * L)$, where $L$ is the virtual plate thickness

  16. *psR(iFish)* : Poisson ratio

  17. *Freq(iFish)* : The flapping frequence for flexible bodies

  18. *firstXYZ* : The initial position of the first point of the bodies

  19. *deltaXYZ* : If there are more than one bodies in a type, *deltaXYZ* determines the interval between front and rear solids

  19. *initXYZVel* : The given initial translatory velocity for the flexible body and the given translational velocity for the rigid body.

  21. *XYZAmpl, XYZPhi* : Parameters for body flapping

       - $XYZ = XYZAmpl * dcos(2.0 * pi * Freq * time + XYZPhi)$

  22. *AoAo, AoAAmpl, AoAPhi* : Parameters for body rotation

       - $Theta = AoAo + AoAAmpl * dcos(2.0 * pi * Freq * time + AoAPhi)$

- **ProbingFluid**

  1. *fluidProbingNum*: The number of data detection points in flow field
  2. *fluidProbingCoords* : The coordinates of detection points

- **ProbingSolid**

  1. *solidProbingNum*:  The number of data detection points in bodies
  2. *solidProbingNode*: The serial number of each detection points in the specific body

## Output file description

- **Output files and directories**

1. *DatBody*      : Folder of body results
2. *DatBodySpan*  : Folder of spanwise-extension body results
3. *DatFlow*      : Floder of flow field results
4. *DatInfo*      : Floder of force and power *et. al.* results
5. *Check.dat*    : Parameters record file for checking
6. *continue.dat* : Results record file for continue calculating

- **Files in DatInfo description**

1. *FishAngular.plt* 
   - *AoA* :  The deflection angle of the body $(y2-y1)/(x2-x1)$;
   - *Ty-TH* : Height difference between the first and last points $(y2 - y1)$
   - *Hy* : *y* Coordinate of the first point $(y1)$
   - *Ty* : *y* Coordinate of the last point $(y2)$
2. *FishPower.plt* : The power information of the bodies
   - *Ptot*  : *Px* + *Py* + *Pz*
   - *Px*,*Py*,*Pz* : Output power calculated by ($force * velcity$)
3. *FishEnergy.plt* : The energy of the bodies
   - *E_s* : Streching strain energy
   - *E_b* : Bending strain energy
   - *E_p* : Total strain energy calculated by (*E_s* + *E_b*)
   - *E_w* : Kinetic energy
   - *E_t* : total energy calculated by (*E_p* + *E_w*)
   - *E_k* : Output energy calculated by ($E_k = E_k + t * ptot$)
4. *FishForce.plt* :The forces exerted on the bodies
5. *FishNodeBegin.plt* : The information of the first point of the bodies 
6. *FishNodeCenter.plt* : The information of the center point of the bodies 
7. *FishNoEnd.plt* : The information of the last point of the bodies 
8. *FishNodeMean.plt* : The average information of the all points of the bodies 

