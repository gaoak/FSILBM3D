! SPDX-License-Identifier: GPL-3.0-or-later
!
! FSILBM3D
! Copyright (C) 2025-2026 Ankang Gao and contributors

module ConstParams
    integer, parameter :: sp = selected_real_kind(6, 37)
#ifdef FSILBM_SINGLE_PRECISION
    integer, parameter :: rp = sp
#else
    integer, parameter :: rp = selected_real_kind(15, 307)
#endif
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
    real(rp), parameter:: wt(0:lbmDim) = [&
       1.0e0_rp/3.0e0_rp,1.0e0_rp/18.0e0_rp,1.0e0_rp/18.0e0_rp,1.0e0_rp/18.0e0_rp,1.0e0_rp/18.0e0_rp,1.0e0_rp/18.0e0_rp,1.0e0_rp/18.0e0_rp, &
                   1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp, &
                   1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp,1.0e0_rp/36.0e0_rp  ]
    !   MRT model Matrix
    !real(rp),parameter::s0=1.0e0_rp,s1=1.0e0_rp,s2=1.0e0_rp,s4=1.0e0_rp,s10=1.0e0_rp,s16=1.0e0_rp
    real(rp), parameter::s0=0.0e0_rp,s1=1.19e0_rp,s2=1.4e0_rp,s4=1.2e0_rp,s10=1.4e0_rp,s16=1.98e0_rp

    !   MODULE BoundCondParams
    integer, parameter:: BCEq_DirecletU = 101,BCnEq_DirecletU = 102,BCorder1_Extrapolate = 103,BCorder2_Extrapolate = 104
                        !given balance function,unbalanced extrapolation,1st order extrapolate,2nd order extrapolate
    integer, parameter:: BCstationary_Wall = 201, BCmoving_Wall = 202, BCstationary_Wall_halfway = 203, BCmoving_Wall_halfway = 204
    integer, parameter:: BCPeriodic = 301,BCSymmetric = 302,BCfluid = 0,BCfluid_father = 1

    real(rp), parameter:: Pi = 3.141592653589793e0_rp,eps = 1.0e-5_rp,MachineTolerace = 1.0e-12_rp
    integer, parameter:: DOFDim = 6

    real(rp), parameter:: Cs2 = 1.e0_rp/3.0e0_rp
    real(rp), parameter:: Csmag = 0.17e0_rp
    ! real(rp), parameter:: CsmagConst = 16.e0_rp * sqrt(2.e0_rp) / (3.e0_rp * Pi * Pi)
    real(rp), parameter:: CsmagConst = 2.0e0_rp * Csmag * Csmag * sqrt(2.e0_rp) * 9.e0_rp
    real(rp), parameter:: CWALE = 0.50e0_rp
    real(rp), parameter:: CWALEConst = CWALE*CWALE
    real(rp), parameter:: CvremConst = 2.5e0_rp*Csmag*Csmag
end module ConstParams
