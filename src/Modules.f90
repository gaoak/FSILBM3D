!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULE LBModel
!    D2Q9model
    integer, parameter:: SpcDim = 2, LBMDim=8 
!   Directions*************************************************************************************
    integer, parameter:: ee(0:LBMDim,1:SpcDim)  = reshape(  [&
                                          !0  1  2  3  4  5  6  7  8
                                           0, 1, 0,-1, 0, 1,-1,-1, 1,&
                                           0, 0, 1, 0,-1, 1, 1,-1,-1 ],[LBMDim+1,SpcDim])
!    Opposite directions****************************************************************************
    integer, parameter:: oppo(0:LBMDim) = [0, 3, 4, 1, 2, 7, 8, 5, 6]
!    Weights****************************************************************************************
    real(8), parameter:: wt(0:LBMDim) = [&
                                4.0d0/ 9.0d0,1.0d0/ 9.0d0,1.0d0/ 9.0d0,1.0d0/ 9.0d0,1.0d0/ 9.0d0, &
                                             1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0  ]
!   MRT model Matrix*******************************************************************************
    real(8)::M_COLLID(0:lbmDim,0:lbmDim),M_FORCE(0:lbmDim,0:lbmDim)
    !real(8),parameter::s0=1.00d0,s1=1.00d0,s2=1.40d0,s4=1.00d0
    !real(8),parameter::s0=1.00d0,s1=1.63d0,s2=1.14d0,s4=1.92d0
    real(8), parameter::s0=1.00d0,s1=1.40d0,s2=1.40d0,s4=1.20d0
    END MODULE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MODULE simParam
    USE LBModel
    integer:: npsize 
!    ***********************************************************************************************
    real(8), parameter:: pi=3.141592653589793d0, eps=1.0d-5
!    ***********************************************************************************************
    integer, parameter:: fluid = 0, wall = 200, MovingWall=201
    integer, parameter:: DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304
                        !given equilibrium function, non-equilibrium extrapolation, first order extrapolation, second order extrapolation, periodic
!   ***********************************************************************************************
!   ***********************************************************************************************
    integer:: isConCmpt, isMoveGrid, isMoveOutput, RefVelocity, RefLength, isKB, isRelease, isVelInit
    integer:: iMoveDim, iCollidModel, iStreamModel, iForce2Flow, iForce2Body, move(1:SpcDim)
    integer:: IXref, IYref, IZref, ntolLBM, ntolFEM,  numsubstep, numSampFlow,numSampBody
    real(8):: Xref,  Yref,  Zref,  dtolLBM, dtolFEM, subdeltat
    integer, allocatable:: SampBodyNode(:,:)
    real(8), allocatable:: SampFlowPint(:,:)
    integer:: Step 
    real(8):: time, Tref, Uref, Aref, Lref, Fref, Eref, Pref, Lthck, Lleng, Palpha, Pbeta,Pramp, uMax
    real(8), allocatable:: Lchod(:)   !The distance between the two farthest points on the body
    real(8):: timeSimTotl, timeOutTemp, timeOutBody, timeOutFlow, timeOutInfo, timeOutFlBg, timeOutFlEd
    real(8):: uIn(1:SpcDim), denIn, g(1:3),AmplInitDist(1:SpcDim),AmplTranDist(1:SpcDim),FreqInitDist(1:SpcDim)
    real(8):: Re, St, Frod(1:3)
    real(8), allocatable::denR(:), tcR(:), EmR(:), psR(:), KB(:), KS(:)
    real(8):: dampK,dampM,Newmarkdelta,Newmarkalpha,alphaf,alpham,alphap

    real(8):: UPre,UNow,Et,Ek,Ep,Es,Eb,Ew
!    ***********************************************************************************************
    character (LEN=40):: LBmeshName
    integer:: boundaryConditions(1:4),MovingKind1,MovingKind2,VelocityKind
    real(8):: VelocityAmp,VelocityFreq,MovingVel1,MovingVel2,MovingFreq1,MovingFreq2,shearRateIn(1:SpcDim),ParabolicParameter(1:2)
    integer:: xDim, yDim
    integer:: xMinBC,xMaxBC,yMinBC,yMaxBC,iBC
    real(8):: dh, dt, ratio
    real(8):: dxmin,dymin,dxmax,dymax,cptxmin,cptxmax,cptymin,cptymax
    real(8):: Cs2, Nu, Mu, Tau, Omega
    real(8):: upxc0, upxcm, upxcmm, upyc0, upycm, upycmm
    real(8):: cnxcm, cnxc0, cnxcp,  cnycm, cnyc0, cnycp
    integer, allocatable:: image(:,:)
    real(8), allocatable:: xGrid(:), yGrid(:), xGrid0(:), yGrid0(:), dx(:), dy(:)
    real(8), allocatable:: fIn(:,:,:), fInTemp(:,:)
    real(8), allocatable:: uuu(:,:,:), force(:,:,:), den(:,:), prs(:,:), vor(:,:)    
     
       
!    ***********************************************************************************************
!   *********************************************************************************************** 
!   *********************************************************************************************** 
    integer, parameter:: idat=12, DOFDim=6  ! Data reading format; Number of body freedom degrees 
    character(LEN=40), allocatable:: FEmeshName(:)
!   ***********************************************************************************************     
    real(8):: deltaT,Freq
    real(8), allocatable:: waittingTime(:,:)      !Dimensionless number for Platehead's waitting time at wave crest and wave trough (t/T)
    real(8), allocatable:: XYZ(:,:),XYZo(:,:),XYZAmpl(:,:),XYZPhi(:,:),XYZd(:,:),UVW(:,:)
    real(8), allocatable:: AoA(:,:),AoAo(:,:),AoAAmpl(:,:),AoAPhi(:,:),AoAd(:,:),WWW1(:,:),WWW2(:,:),WWW3(:,:)
    real(8), allocatable:: TTT00(:,:,:),TTT0(:,:,:),TTTnow(:,:,:),TTTnxt(:,:,:)
!   ***********************************************************************************************
    integer:: nFish
    integer:: nND_max,nEL_max,nMT_max,nEQ_max
    ! nND(1:nFish), number of nodes in each fish
    ! nEL(1:nFish), number of elements in each fish
    ! nMT(1:nFish), 1 for 2D
    integer, allocatable:: nND(:),nEL(:),nEQ(:),nMT(:),nBD(:),maxstiff(:)
    integer, allocatable:: isMotionGiven(:,:),iBodyModel(:)
    integer:: NDref 
!   ===============================================================================================
    integer:: nEL_all, nND_all
    integer,allocatable:: ele_all(:,:)
    real(8),allocatable:: xyzful_all(:,:),velful_all(:,:),xyzfulIB_all(:,:),extful_all(:,:)
!   ===============================================================================================    
    integer, allocatable:: ele(:,:,:),nloc(:,:),nprof(:,:),nprof2(:,:),jBC(:,:,:)
    ! xyzful00(1:nFish, nND(fish), 1:6) coordinates of nodes in each fish from the input file, only (1:2) are used
    real(8), allocatable:: xyzful00(:,:,:),mssful(:,:,:),vBC(:,:,:),mss(:,:),prop(:,:,:)
    real(8), allocatable:: lodful(:,:,:),extful(:,:,:),extful1(:,:,:),extful2(:,:,:),grav(:,:,:),repful(:,:,:)

    ! xyzful0, the initial coordinate of body nodes. Translation and rotation are applied
    real(8), allocatable:: xyzful0(:,:,:),xyzfulnxt(:,:,:),dspful(:,:,:),accful(:,:,:)
    ! xyzful the current coordinates of nodes in each fish
    ! velful the current velocity of nodes in each fish
    ! xyzfulIB the virtual IB nodes
    real(8), allocatable:: xyzful(:,:,:),xyzfulIB(:,:,:),velful(:,:,:)
    
    real(8), allocatable:: triad_nn(:,:,:,:),triad_ee(:,:,:,:),triad_e0(:,:,:,:)
    real(8), allocatable:: triad_n1(:,:,:,:),triad_n2(:,:,:,:),triad_n3(:,:,:,:)
!    ***********************************************************************************************

    END MODULE
