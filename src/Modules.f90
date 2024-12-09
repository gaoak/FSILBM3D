MODULE ConstParams
    ! LBM module
    !    D3Q19model
    INTEGER, PARAMETER:: SpaceDim = 3, LBMDim = 18
    !   Directions
    integer, parameter:: ee(0:lbmDim,1:3)=reshape([&
                                        !0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
                                         0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, &
                                         0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, &
                                         0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1  &
                                                 ],[lbmDim+1,SpaceDim])
!    Opposite directions
    integer, parameter:: oppo(0:lbmDim)=[0, 2, 1, 4, 3, 6, 5,10, 9, 8, 7,14,13,12,11,18,17,16,15]
    integer, parameter:: positivedirs(1:lbmDim/2)=[1, 3, 5, 7, 8, 11, 12, 15, 16]
    integer, parameter:: negativedirs(1:lbmDim/2)=[2, 4, 6,10, 9, 14, 13, 18, 17]
!    Weights
    real(8), parameter:: wt(0:lbmDim) = [&
       1.0d0/3.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0, &
                   1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0, &
                   1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0  ]
!   MRT model Matrix
    real(8)::M_COLLID(0:lbmDim,0:lbmDim),M_FORCE(0:lbmDim,0:lbmDim)
    !real(8),parameter::s0=1.0d0,s1=1.0d0,s2=1.0d0,s4=1.0d0,s10=1.0d0,s16=1.0d0
    real(8), parameter::s0=0.0d0,s1=1.19d0,s2=1.4d0,s4=1.2d0,s10=1.4d0,s16=1.98d0
!   

    !MODULE BoundCondParams
    integer, parameter:: BCfluid = 0, BCwall = 200, BCmovingWall=201
    integer, parameter:: BCDirecletUP=300,BCDirecletUU=301,BCAdvection1=302,BCAdvection2=303,BCPeriodic=304, BCSymmetric=101
                        !given balance function,unbalanced extrapolation,1st order extrapolate,2nd order extrapolate,periodic
    
    real(8), parameter:: Pi=3.141592653589793d0,eps=1.0d-5

    integer, parameter:: DOFDim=6

END MODULE ConstParams
