module ConstParams
    !   LBM module
    !   D3Q19model
    integer, parameter:: SpaceDim = 3, lbmDim = 18
    !   Directions
    integer, parameter:: ee(0:lbmDim,1:3) = reshape([&
                                        !0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
                                         0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0, &
                                         0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1, &
                                         0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1  &
                                                 ],[lbmDim+1,SpaceDim])
    !   Opposite directions
    integer, parameter:: oppo(0:lbmDim) = [0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15]
    integer, parameter:: positivedirs(1:lbmDim/2) = [1, 3, 5, 7, 8, 11, 12, 15, 16]
    integer, parameter:: negativedirs(1:lbmDim/2) = [2, 4, 6,10, 9, 14, 13, 18, 17]
    !   Weights
    real(8), parameter:: wt(0:lbmDim) = [&
       1.0d0/3.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0,1.0d0/18.0d0, &
                   1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0, &
                   1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0  ]
    !   MRT model Matrix
    !real(8),parameter::s0=1.0d0,s1=1.0d0,s2=1.0d0,s4=1.0d0,s10=1.0d0,s16=1.0d0
    real(8), parameter::s0=0.0d0,s1=1.19d0,s2=1.4d0,s4=1.2d0,s10=1.4d0,s16=1.98d0

    !   MODULE BoundCondParams
    integer, parameter:: Eq_DirecletU = 101,nEq_DirecletU = 102,order1_Extrapolate = 103,order2_Extrapolate = 104
                        !given balance function,unbalanced extrapolation,1st order extrapolate,2nd order extrapolate
    integer, parameter:: stationary_Wall = 201, moving_Wall = 202
    integer, parameter:: Periodic = 301,Symmetric = 302

    real(8), parameter:: Pi = 3.141592653589793d0,eps = 1.0d-5
    integer, parameter:: DOFDim = 6

    real(8), parameter:: Cs2 = 1.d0/3.0d0

end module ConstParams
