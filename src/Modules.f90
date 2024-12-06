!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Lattice Boltzman model module
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULE LBModel
!    D3Q19model*************************************************************************************
    INTEGER, PARAMETER:: SpcDim = 3, LBMDim = 18
!   Directions*************************************************************************************
    integer, parameter:: ee(0:lbmDim,1:3)=reshape([&
                                        !0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
                                         0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, &
                                         0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, &
                                         0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1  &
                                                 ],[lbmDim+1,SpcDim])
!    Opposite directions****************************************************************************
    integer, parameter:: oppo(0:lbmDim)=[0, 2, 1, 4, 3, 6, 5,10, 9, 8, 7,14,13,12,11,18,17,16,15]
!    Weights****************************************************************************************
    real(8), parameter:: wt(0:lbmDim) = [&
       1.0d0/3.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0, &
                   1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0, &
                   1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0  ]
!   MRT model Matrix*******************************************************************************
    real(8)::M_COLLID(0:lbmDim,0:lbmDim),M_FORCE(0:lbmDim,0:lbmDim)
    !real(8),parameter::s0=1.0d0,s1=1.0d0,s2=1.0d0,s4=1.0d0,s10=1.0d0,s16=1.0d0
    real(8), parameter::s0=1.0d0,s1=1.19d0,s2=1.4d0,s4=1.2d0,s10=1.4d0,s16=1.98d0
!   ***********************************************************************************************
    END MODULE LBModel

    MODULE BoundCondParams
    integer, parameter:: fluid = 0, wall = 200, movingWall=201
    integer, parameter:: DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304, Symmetric=101
                        !given balance function,unbalanced extrapolation,1st order extrapolate,2nd order extrapolate,periodic
    END MODULE BoundCondParams

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Fluid structure interaction module
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULE simParam
    USE LBModel
    USE BoundCondParams
    integer:: npsize
    real(8), parameter:: Pi=3.141592653589793d0,eps=1.0d-5
!   ***********************************************************************************************
    integer:: step
    real(8):: time

    integer:: isConCmpt,iCollidModel,iStreamModel,iForce2Body,iKB,isRelease,RefVelocity,RefTime,isFluidOutput,isBodyOutput
    integer:: move(1:SpcDim),offsetOutput
    integer:: isMoveGrid,isMoveDimX,isMoveOutputX,isMoveDimY,isMoveOutputY,isMoveDimZ,isMoveOutputZ,MoveOutputIref(1:3)
    logical:: isUniformGrid(1:SpcDim)
    integer:: IXref,IYref,IZref,ntolLBM,ntolFEM,ntolFSI,numsubstep,numSampFlow,numSampBody
    integer:: boundaryConditions(1:6),MovingKind1,MovingKind2,VelocityKind
    real(8):: VelocityAmp,VelocityFreq,VelocityPhi,MovingVel1,MovingVel2,MovingFreq1,MovingFreq2
    real(8):: VolumeForce(1:SpcDim),VolumeForceAmp,VolumeForceFreq,VolumeForcePhi,VolumeForceIn(1:SpcDim)
    integer, allocatable:: SampBodyNode(:,:), iBodyModel(:), iBodyType(:)
    real(8), allocatable:: SampFlowPint(:,:)
    integer:: nFish,NDref
    real(8):: Xref,Yref,Zref
    real(8):: timeSimTotl,timeOutTemp,timeOutBody,timeOutFlow,timeOutInfo,timeOutBegin,timeOutEnd
    real(8):: dtolLBM,Palpha,Pbeta,Pramp,uMax,dtolFEM,dtolFSI,deltat,subdeltat
    real(8):: uuuIn(1:SpcDim),shearRateIn(1:SpcDim),denIn,g(1:SpcDim)
    real(8):: AmplInitDist(1:SpcDim),waveInitDist,AmplforcDist(1:SpcDim),FreqforcDist
    real(8):: posiForcDist(1:SpcDim),begForcDist,endForcDist
    real(8):: Re,AR,Frod(1:SpcDim)
    real(8):: dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,alpham,alphap
    real(8):: Uref,Lref,Tref,Aref,Fref,Eref,Pref,Lthck,Lchod,Lspan,Asfac

    real(8):: upxc0, upxcm, upxcmm
    real(8):: upyc0, upycm, upycmm
    real(8):: upzc0, upzcm, upzcmm
    real(8):: cnxc0, cnxcm, cnxcp
    real(8):: cnyc0, cnycm, cnycp
    real(8):: cnzc0, cnzcm, cnzcp

!   ***********************************************************************************************
    character (LEN=40):: LBmeshName
    integer:: xDim,yDim,zDim
    integer:: xMinBC,xMaxBC,yMinBC,yMaxBC,zMinBC,zMaxBC,iBC
    real(8):: dh,dt,ratio,dxmin,dymin,dzmin,dxmax,dymax,dzmax
    real(8):: cptxMin,cptxMax,cptyMin,cptyMax,cptzMin,cptzMax
    real(8):: Omega,tau,Cs2,nu,Mu
    real(8), allocatable:: dx(:), dy(:), dz(:)
    real(8), allocatable:: xGrid0(:), yGrid0(:), zGrid0(:), xGrid(:), yGrid(:), zGrid(:)
    real(8), allocatable:: fIn(:,:,:,:), fInTemp(:,:,:)
    real(8), allocatable:: uuu(:,:,:,:), force(:,:,:,:), den(:,:,:), prs(:,:,:)

!   ***********************************************************************************************
!   ***********************************************************************************************
    integer, parameter:: DOFDim=6
    END MODULE simParam

    MODULE PartitionXDim
        integer:: npsize_copy, xDim_copy
        integer, allocatable:: partition(:), parindex(:),eid(:)
        real(8), allocatable:: edge(:,:,:)
    END MODULE

    MODULE OutFlowWorkspace
        real(4), allocatable:: oututmp(:,:,:),outvtmp(:,:,:),outwtmp(:,:,:)
        real(4):: offsetMoveGrid(1:3)
    END MODULE
