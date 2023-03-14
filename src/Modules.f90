!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Lattice Boltzman model module
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULE LBModel
!    D3Q19model*************************************************************************************
    INTEGER, PARAMETER:: SpcDim = 3, LBMDim = 18
!   Directions*************************************************************************************
    integer, parameter:: ee(0:lbmDim,1:3)=reshape([&
                                        !0  1  2     3    4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
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

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Fluid structure interaction module
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULE simParam
    USE LBModel
    integer:: npsize
    real(8), parameter:: Pi=3.141592653589793d0,eps=1.0d-5
!   ***********************************************************************************************
    integer, parameter:: fluid = 0, wall = 200, movingWall=201
    integer, parameter:: DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304
                        !given balance function,unbalanced extrapolation,1st order extrapolate,2nd order extrapolate,periodic
!   ***********************************************************************************************
!   ***********************************************************************************************
    integer:: step
    real(8):: time

    integer:: isConCmpt,iCollidModel,iStreamModel,iForce2Body,iKB,isRelease,RefVelocity
    integer:: iChordDirection,move(1:SpcDim),numOutput
    integer:: isMoveGrid,isMoveDimX,isMoveOutputX,isMoveDimY,isMoveOutputY,isMoveDimZ,isMoveOutputZ
    integer:: IXref,IYref,IZref,ntolLBM,ntolFEM,ntolFSI,numsubstep,numSampFlow,numSampBody
    integer:: boundaryConditions(1:6)
    integer, allocatable:: SampBodyNode(:,:), iBodyModel(:)
    real(8), allocatable:: SampFlowPint(:,:)
    real(8):: Xref,Yref,Zref
    real(8):: timeSimTotl,timeOutTemp,timeOutBody,timeOutFlow,timeOutInfo,timeOutFlBg,timeOutFlEd
    real(8):: dtolLBM,Palpha,Pbeta,Pramp,uMax,dtolFEM,dtolFSI,subdeltat
    real(8):: uuuIn(1:SpcDim),shearRateIn(1:SpcDim),denIn,g(1:SpcDim)
    real(8):: AmplInitDist(1:SpcDim),waveInitDist,AmplforcDist(1:SpcDim),FreqforcDist
    real(8):: posiForcDist(1:SpcDim),begForcDist,endForcDist
    real(8):: Re,AR,Frod(1:SpcDim)
    real(8), allocatable:: denR(:),KB(:),KS(:),EmR(:),psR(:),tcR(:),St(:)
    real(8):: dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,alpham,alphap
    real(8):: Uref,Lref,Tref,Aref,Fref,Eref,Pref,Lthck,Lchod,Lspan,Asfac
    real(8), allocatable:: nAsfac(:),nLchod(:),nLspan(:)
    real(8):: UPre,UNow,Et,Ek,Ep,Es,Eb,Ew

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
    real(8), allocatable:: elmax(:),elmin(:)
    real(8):: dh,dt,ratio,dxmin,dymin,dzmin,dxmax,dymax,dzmax
    real(8):: cptxMin,cptxMax,cptyMin,cptyMax,cptzMin,cptzMax
    real(8):: Omega,tau,Cs2,nu,Mu
    integer, allocatable:: image(:,:,:)
    real(8), allocatable:: dx(:), dy(:), dz(:) 
    real(8), allocatable:: xGrid0(:), yGrid0(:), zGrid0(:), xGrid(:), yGrid(:), zGrid(:)
    real(8), allocatable:: fIn(:,:,:,:), fInTemp(:,:,:)
    real(8), allocatable:: uuu(:,:,:,:), force(:,:,:,:), den(:,:,:), prs(:,:,:)
    
!   *********************************************************************************************** 
!   *********************************************************************************************** 
    integer, parameter:: idat=12, DOFDim=6
    character (LEN=40), allocatable:: FEmeshName(:)
!   ***********************************************************************************************     
    real(8):: deltaT
    real(8), allocatable:: Freq(:)
    real(8), allocatable:: XYZ(:,:),XYZo(:,:),XYZAmpl(:,:),XYZPhi(:,:),XYZd(:,:),UVW(:,:)
    real(8), allocatable:: AoA(:,:),AoAo(:,:),AoAAmpl(:,:),AoAPhi(:,:),AoAd(:,:),WWW1(:,:),WWW2(:,:),WWW3(:,:)
    real(8), allocatable:: TTT00(:,:,:),TTT0(:,:,:),TTTnow(:,:,:),TTTnxt(:,:,:)
!   ***********************************************************************************************
    integer:: nFish,nEL_all,nND_all,nND_max,nEL_max,nMT_max,nEQ_max,NDref
    integer,allocatable:: ele_all(:,:)
    integer, allocatable:: nND(:),nEL(:),nEQ(:),nMT(:),nBD(:),nSTF(:) 
    integer, allocatable:: NDtl(:,:),NDhd(:,:),NDct(:)
    real(8), allocatable:: xyzful_all(:,:),velful_all(:,:),xyzfulIB_all(:,:),extful1_all(:,:),extful2_all(:,:),isMotionGiven(:,:)
!   ===============================================================================================
    integer, allocatable:: ele(:,:,:),jBC(:,:,:),nloc(:,:),nprof(:,:),nprof2(:,:)
    real(8), allocatable:: xyzful00(:,:,:),mssful(:,:,:),vBC(:,:,:),prop(:,:,:),mss(:,:),areaElem00(:,:),areaElem(:,:)
    real(8), allocatable:: lodful(:,:,:),repful(:,:,:),extful(:,:,:),extful1(:,:,:),extful2(:,:,:),grav(:,:,:),streI(:,:),bendO(:,:)

    real(8), allocatable:: xyzful0(:,:,:),xyzfulnxt(:,:,:),dspful(:,:,:),accful(:,:,:)
    real(8), allocatable:: xyzful(:,:,:),xyzfulIB(:,:,:),velful(:,:,:)
    real(8), allocatable:: triad_nn(:,:,:,:),triad_ee(:,:,:,:),triad_e0(:,:,:,:)
    real(8), allocatable:: triad_n1(:,:,:,:),triad_n2(:,:,:,:),triad_n3(:,:,:,:)
!   ***********************************************************************************************   
    END MODULE simParam
