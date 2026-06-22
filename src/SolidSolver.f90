! SPDX-License-Identifier: GPL-3.0-or-later
!
! FSILBM3D
! Copyright (C) 2025-2026 Ankang Gao and contributors

! program algorithm: ISBN 9781441929105 James F. Doyle. P354
module SegmentStructure
    !use mkl
    use ConstParams, only: rp
    implicit none
    private
    integer, parameter :: nElmtDofs = 12
    public :: Segment, Segment_get_angle_triad
    type :: Segment
        integer :: node0, node1, m_localToGlobal(1:nElmtDofs), itype, bc(1:nElmtDofs), Nspan
        real(rp) :: x00(1:nElmtDofs),x0(1:nElmtDofs),x1(1:nElmtDofs),xnxt(1:nElmtDofs)
        real(rp) :: dx0, dy0, dz0, dx1, dy1, dz1, xll0, xmm0, xnn0, xll1, xmm1, xnn1, len0, len1
        real(rp) :: Lspan, spanlen, dirc00(3), dirc0(3), dirc1(3), dircnxt(3)
        real(rp) :: geoFRM
        real(rp) :: areaElem00
        real(rp) :: strainEnergy(2)
        real(rp) :: triad_ee(3,3),triad_n1(3,3),triad_n2(3,3)
        real(rp) :: m_property(1:8)
        real(rp) :: m_coefMat(1:nElmtDofs, 1:nElmtDofs)
        real(rp) :: m_tanMat(1:nElmtDofs, 1:nElmtDofs)
        real(rp) :: m_stfMat(1:nElmtDofs, 1:nElmtDofs)
        real(rp) :: m_masMat(1:nElmtDofs, 1:nElmtDofs)
        real(rp) :: m_geoMat(1:nElmtDofs, 1:nElmtDofs)
        real(rp) :: m_rotMat(1:3, 1:3)
    contains
        procedure :: Build => Segment_Build
        procedure :: Init => Segment_Init
        procedure :: cptdxyz1 => Segment_cptdxyz1
        procedure :: Multiply => Segment_Multiply
        procedure :: UpdateMatrix => Segment_UpdateMatrix
        procedure :: ScalarPreconditioned => Segment_ScalarPreconditioned
        procedure :: BlockPreconditioned => Segment_BlockPreconditioned
        procedure :: UpdateLoad => Segment_UpdateLoad
        procedure :: LocToGlobal => Segment_LocToGlobal
        procedure :: GlobalToLoc => Segment_GlobalToLoc
        procedure :: FormMassMatrix => Segment_FormMassMatrix
        procedure :: FormStiffMatrix => Segment_FormStiffMatrix
        procedure :: FormGeomMatrix => Segment_FormGeomMatrix
        procedure :: RotateMatrix => Segment_RotateMatrix
        procedure :: MassMultiply => Segment_MassMultiply
        procedure :: RKR => Segment_RKR
        procedure :: BoundaryCond => Segment_BoundaryCond
        procedure :: BodyStress_D => Segment_BodyStress_D
        procedure :: StrainEnergy_D => Segment_StrainEnergy_D
        procedure :: InitTriad_D => Segment_InitTriad_D
        procedure :: RigidUpdateTriad_D => Segment_RigidUpdateTriad_D
        procedure :: BuildAxisDirTriad => Segment_BuildAxisDirTriad
        procedure :: UpdateTriad_D => Segment_UpdateTriad_D
        procedure :: MakeTriad_ee => Segment_MakeTriad_ee
        procedure :: MapReferencePosToCurrent => Segment_MapReferencePosToCurrent
        procedure :: MapReferenceDirToCurrent => Segment_MapReferenceDirToCurrent
    end type Segment

  contains

    subroutine Segment_Build(this, p0Id, p1Id, nND, itype_, Nspan_, xyz, material, boundary)
        implicit none
        class(Segment), intent(inout) :: this
        integer, intent(in) :: p0Id, p1Id, nND
        real(rp), intent(in) :: xyz(1:8, 1:nND), material(1:8)
        integer, intent(in) :: boundary(1:6, 1:nND), itype_, Nspan_
        integer :: i, offset0, offset1
        real(rp) :: dirc_norm
        this%node0 = p0Id
        this%node1 = p1Id
        ! material property
        this%m_property(1:8) = material(1:8)
        !local dof to global dof mapping
        offset0 = (p0Id - 1) * 6
        offset1 = (p1Id - 1) * 6
        do i=1,6
            this%m_localToGlobal(i) = offset0 + i
            this%m_localToGlobal(i+6) = offset1 + i
        enddo
        this%x00(1:3) = xyz(1:3, p0Id)
        this%x00(7:9) = xyz(1:3, p1Id)
        this%x00(4:6) = 0.0e0_rp
        this%x00(10:12) = 0.0e0_rp
        this%areaElem00 = sqrt((this%x00(7)-this%x00(1))**2+&
                                (this%x00(8)-this%x00(2))**2+&
                                (this%x00(9)-this%x00(3))**2)
        this%bc(1:6) = boundary(1:6, p0Id)
        this%bc(7:12) = boundary(1:6, p1Id)
        this%itype = itype_
        this%Nspan = Nspan_
        this%Lspan = 0.5e0_rp * (xyz(4,p0Id) + xyz(4,p1Id))
        this%spanlen = 0.5e0_rp * (xyz(5,p0Id) + xyz(5,p1Id)) + this%Lspan
        this%dirc00(1:3) = 0.5e0_rp * (xyz(6:8,p0Id) + xyz(6:8,p1Id))
        dirc_norm = sqrt(sum(this%dirc00**2))
        if (dirc_norm .gt. 1e-10_rp) then
            this%dirc00 = this%dirc00 / dirc_norm
        else
            write(*,*) xyz(6:8,p0Id), "and", xyz(6:8,p1Id), "are opposite directions; no unique bisector exists."
            stop
        endif
    end subroutine Segment_Build

    subroutine Segment_Init(this)
        implicit none
        class(Segment), intent(inout) :: this
        this%x1(1:12) = this%x0(1:12)
        this%xnxt(1:12) = this%x1(1:12)
        this%dirc1(1:3) = this%dirc0(1:3)
        this%dircnxt(1:3) = this%dirc1(1:3)
        this%dx0  = this%x0(7) - this%x0(1)
        this%dy0  = this%x0(8) - this%x0(2)
        this%dz0  = this%x0(9) - this%x0(3)
        this%len0 = sqrt(this%dx0*this%dx0+this%dy0*this%dy0+this%dz0*this%dz0)
        this%xll0 = this%dx0/this%len0
        this%xmm0 = this%dy0/this%len0
        this%xnn0 = this%dz0/this%len0
        call this%cptdxyz1()
    end subroutine Segment_Init

    subroutine Segment_cptdxyz1(this)
        implicit none
        class(Segment), intent(inout) :: this
        this%dx1  = this%x1(7) - this%x1(1)
        this%dy1  = this%x1(8) - this%x1(2)
        this%dz1  = this%x1(9) - this%x1(3)
        this%len1 = sqrt(this%dx1*this%dx1+this%dy1*this%dy1+this%dz1*this%dz1)
        this%xll1 = this%dx1/this%len1
        this%xmm1 = this%dy1/this%len1
        this%xnn1 = this%dz1/this%len1
    end subroutine Segment_cptdxyz1

    subroutine Segment_UpdateMatrix(this, coeffs, gamma, dampM, dampK)
        ! Form the tangent stiffness matrix
        ! ISBN 9781441929105 James F. Doyle. P268
        ! [C]=dampM*[M]+dampK*[K]
        ! ISBN 9781441929105 James F. Doyle. P193,196
        ! K[T] = K[E] + gamma*K[G]
        ! ISBN 9781441929105 James F. Doyle. P355
        ! K[T] = K[E] + K[G]
        ! K = K[T] + coeffs(1)*[C] + coeffs(0)*[M]
        class(Segment), intent(inout) :: this
        real(rp), intent(in) :: coeffs(0:7), gamma, dampM, dampK
        !update m_coefMat
        call this%FormGeomMatrix

        call this%RotateMatrix

        call this%RKR(this%m_stfMat)

        call this%RKR(this%m_geoMat)
        ! The RKR of this%m_stfMat and this%m_geoMat cannot be merged, calculate dampK use this%m_stfMat after RKR !
        this%m_tanMat = this%m_stfMat + gamma * this%m_geoMat
        this%m_coefMat = this%m_tanMat + coeffs(0) * this%m_masMat + coeffs(1) * dampM * this%m_masMat
        if(dampK .gt. 0.0e0_rp) then
            this%m_coefMat = this%m_coefMat + coeffs(1) * dampK * this%m_stfMat
        endif
        return
    end subroutine Segment_UpdateMatrix

    subroutine Segment_UpdateLoad(this, coeffs, dampM, dampK, nND, gEQ, dspO, dsp, vel, acc, lodEffe)
        ! Form the effective load vector {F}_(t+delta t)
        ! ISBN 9787302388333 Xiong Zhang. P113-114
        ! [C]=dampM*[M]+dampK*[K]
        ! ISBN 9781441929105 James F. Doyle. P268
        ! Newton-Raphson method
        ! ISBN 9781441929105 James F. Doyle. P353 (5.11): The book shows the full Newton-Raphson method with alpha=0.25 and delta=0.5
        ! ISBN 9787576318555 Dong Chunying. P142 (8.84)-(8.89)
        implicit none
        class(Segment), intent(inout) :: this
        integer, intent(in) :: nND, gEQ
        real(rp), intent(in) :: coeffs(0:7), dampM, dampK
        real(rp), intent(in) :: dspO(1:6, 1:nND)
        real(rp), intent(in) :: dsp(1:6, 1:nND), vel(1:6, 1:nND), acc(1:6, 1:nND)
        real(rp), intent(inout) :: lodEffe(1:gEQ)
        real(rp) :: qM(1:nElmtDofs), qC(1:nElmtDofs), qMC(1:nElmtDofs)
        real(rp) :: massLoad(1:nElmtDofs),dampKLoad(1:nElmtDofs)
        integer :: i,node0,node1
        node0 = this%node0
        node1 = this%node1
        do i = 1, 6
            qM(i)   = coeffs(0)*(dspO(i,node0)-dsp(i,node0)) + coeffs(2)*vel(i,node0) + coeffs(3)*acc(i,node0)
            qM(i+6) = coeffs(0)*(dspO(i,node1)-dsp(i,node1)) + coeffs(2)*vel(i,node1) + coeffs(3)*acc(i,node1)
            qC(i)   = coeffs(1)*(dspO(i,node0)-dsp(i,node0)) + coeffs(4)*vel(i,node0) + coeffs(5)*acc(i,node0)
            qC(i+6) = coeffs(1)*(dspO(i,node1)-dsp(i,node1)) + coeffs(4)*vel(i,node1) + coeffs(5)*acc(i,node1)
        enddo
        qMC = qM + dampM*qC
        call this%MassMultiply(qMC, massLoad)
        call this%LocToGlobal(massLoad, lodEffe, gEQ)
        if (dampK .gt. 0.0e0_rp) then
            dampKLoad = dampK * matmul(this%m_stfMat, qC)
            call this%LocToGlobal(dampKLoad, lodEffe, gEQ)
        endif
        return
    end subroutine Segment_UpdateLoad

    subroutine Segment_MassMultiply(this, q, mq)
        ! Lumped beam mass matrix in 3-by-3 block form:
        !
        ! DOFs: [u1 v1 w1 | rx1 ry1 rz1 | u2 v2 w2 | rx2 ry2 rz2]
        !
        !               1:3      4:6      7:9     10:12
        !            +--------+--------+--------+--------+
        !      1:3   |  Mt1   |   0    |   0    |   0    |
        !            +--------+--------+--------+--------+
        !      4:6   |   0    |  Mr1   |   0    |   0    |
        ! M =        +--------+--------+--------+--------+
        !      7:9   |   0    |   0    |  Mt2   |   0    |
        !            +--------+--------+--------+--------+
        !     10:12  |   0    |   0    |   0    |  Mr2   |
        !            +--------+--------+--------+--------+
        !
        ! Mt1 and Mt2 are isotropic translational blocks, Mt1=m1*I_3, Mt2=m2*I_3.
        ! Mr1 and Mr2 are rotational inertia blocks and may be full 3-by-3 matrices after rotation to the global frame.
        !
        ! Hence:
        !
        !   M*q = [ Mt1*q(1:3),
        !           Mr1*q(4:6),
        !           Mt2*q(7:9),
        !           Mr2*q(10:12) ].
        implicit none
        class(Segment), intent(in) :: this
        real(rp), intent(in)  :: q(1:nElmtDofs)
        real(rp), intent(out) :: mq(1:nElmtDofs)
        integer :: i
        mq(:) = 0.0e0_rp
        ! translational lumped massLoad blocks: m I_3
        do i = 1, 3
            mq(i)   = this%m_masMat(i,i)     * q(i)
            mq(i+6) = this%m_masMat(i+6,i+6) * q(i+6)
        enddo
        ! rotational inertia blocks: full 3x3 after rotation
        mq(4:6)   = matmul(this%m_masMat(4:6,4:6), q(4:6))
        mq(10:12) = matmul(this%m_masMat(10:12,10:12), q(10:12))
    end subroutine Segment_MassMultiply

    subroutine Segment_BoundaryCond(this, iter, gEQ, x, fixed, vBC)
        ! Build fixed-DOF mask and prescribed displacement increment.
        ! For Newton iteration 1, fixed DOFs take the non-homogeneous prescribed increment vBC. 
        ! For later Newton iterations, the correction on fixed DOFs is homogeneous, i.e. zero.
        implicit none
        class(Segment), intent(in) :: this
        integer, intent(in) :: iter, gEQ
        real(rp), intent(in) :: vBC(1:gEQ)
        real(rp), intent(inout) :: x(1:gEQ)
        logical, intent(inout) :: fixed(1:gEQ)
        integer :: i, gid
        do i=1,nElmtDofs
            if (this%bc(i) .gt. 0) then
                gid = this%m_localToGlobal(i)
                fixed(gid) = .true.
                if (iter .eq. 1) then
                    x(gid) = vBC(gid)
                else
                    x(gid) = 0.0e0_rp
                endif
            endif
        enddo
    end subroutine Segment_BoundaryCond

    subroutine Segment_ScalarPreconditioned(this, M, gEQ)
        class(Segment), intent(in) :: this
        integer, intent(in) :: gEQ
        real(rp) :: Melmts(1:nElmtDofs),M(1:gEQ)
        integer :: i
        do i=1,nElmtDofs
            Melmts(i)=this%m_coefMat(i,i)
        enddo
        call this%LocToGlobal(Melmts, M, gEQ)
    end subroutine Segment_ScalarPreconditioned

    subroutine Segment_BlockPreconditioned(this, blockM, nND)
        implicit none
        class(Segment), intent(in) :: this
        integer, intent(in) :: nND
        real(rp), intent(inout) :: blockM(6,6,nND)
        blockM(1:6,1:6,this%node0) = blockM(1:6,1:6,this%node0) + this%m_coefMat(1:6,1:6)
        blockM(1:6,1:6,this%node1) = blockM(1:6,1:6,this%node1) + this%m_coefMat(7:12,7:12)
    end subroutine Segment_BlockPreconditioned

    subroutine Segment_Multiply(this, x, b, gEQ)
        class(Segment), intent(in) :: this
        integer, intent(in) :: gEQ
        real(rp), intent(in) :: x(1:gEQ)
        real(rp), intent(inout) :: b(1:gEQ)
        real(rp) :: lx(1:nElmtDofs), lb(1:nElmtDofs)
        lx = 0.0e0_rp
        lb = 0.0e0_rp
        call this%GlobalToLoc(x, lx, gEQ)
        lb = matmul(this%m_coefMat, lx)
        call this%LocToGlobal(lb, b, gEQ)
        return
    end subroutine Segment_Multiply

    subroutine Segment_GlobalToLoc(this, x, lx, gEQ)
        class(Segment), intent(in) :: this
        integer, intent(in) :: gEQ
        real(rp), intent(in) :: x(1:gEQ)
        real(rp), intent(out) :: lx(1:nElmtDofs)
        integer :: i
        do i=1,nElmtDofs
            lx(i) = x(this%m_localToGlobal(i))
        enddo
        return
    end subroutine Segment_GlobalToLoc

    subroutine Segment_LocToGlobal(this, lx, x, gEQ)
        class(Segment), intent(in) :: this
        integer, intent(in) :: gEQ
        real(rp), intent(in) :: lx(1:nElmtDofs)
        real(rp), intent(inout) :: x(1:gEQ)
        integer :: i
        do i=1,nElmtDofs
            x(this%m_localToGlobal(i)) = x(this%m_localToGlobal(i)) + lx(i)
        enddo
        return
    end subroutine Segment_LocToGlobal

    subroutine Segment_FormMassMatrix(this)
        ! ELeMent MASs matrix for 3D Timoshenko FRaMe
        ! Same as Abaqus B31 Timoshenko frame
        ! Lumped massLoad matrix
        ! ISBN 9781441929105 James F. Doyle. P273
        ! ISBN 9780792312086 James F. Doyle. P423
        !
        ! Rotary inertia about the local x-axis uses polar second moment Ip = Iy + Iz.
        ! This is not generally equal to the St. Venant torsion constant Jt, except for circular sections.
        implicit none
        class(Segment), intent(inout) :: this
        real(rp):: area,rho,zix,ziy,ziz,length
        real(rp):: roal

        area = this%m_property(3) ! A    : cross-sectional area 
        rho  = this%m_property(4) ! rho  : density
        zix  = this%m_property(6) ! Jt   : St. Venant torsion constant, used in torsional stiffness G*Jt
        ziy  = this%m_property(7) ! Iy   : second moment of area about local y-axis
        ziz  = this%m_property(8) ! Iz   : second moment of area about local z-axis
        length  = this%len0

        this%m_masMat(1:12,1:12) = 0.0e0_rp
        roal = rho*area*length/2.0e0_rp
        ! Translational lumped mass at node 0.
        this%m_masMat(1,1)     = roal
        this%m_masMat(2,2)     = roal
        this%m_masMat(3,3)     = roal

        ! Rotary inertia per node.
        ! The rotational inertia about the local beam axis uses the polar second
        ! moment of area Ip = Iy + Iz, not the St. Venant torsion constant Jt.
        ! For circular sections, Jt = Ip.
        ! For rectangular sections, Jt generally differs from Ip, so (Iy + Iz) should be used here.
        this%m_masMat(4,4)     = roal*(ziy+ziz)/area
        this%m_masMat(5,5)     = roal*ziy/area
        this%m_masMat(6,6)     = roal*ziz/area

        ! Lumped mass at node 1.
        this%m_masMat(7,7)     = this%m_masMat(1,1)
        this%m_masMat(8,8)     = this%m_masMat(2,2)
        this%m_masMat(9,9)     = this%m_masMat(3,3)
        this%m_masMat(10,10)   = this%m_masMat(4,4)
        this%m_masMat(11,11)   = this%m_masMat(5,5)
        this%m_masMat(12,12)   = this%m_masMat(6,6)
        return
    end subroutine Segment_FormMassMatrix

    subroutine Segment_FormStiffMatrix(this)
        ! ELeMent STiFfness for 3D Timoshenko FRaMe
        ! calculates the element stiffness matrices.
        ! https://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
        ! Henri Gavin, Department of Civil and Environmental Engineering, Duke University
        ! For Euler-Bernoulli :
        ! ISBN 9787040258417 Zeng Pan. P70
        ! ISBN 9780792312086 James F. Doyle. P81
        !
        ! The default shear correction factors ksy=ksz=5/6 are for rectangular sections.
        ! Set phiy = phiz = 0 for Euler-Bernoulli beam behavior.
    
        implicit none
        class(Segment), intent(inout) :: this
    
        real(rp):: emod,gmod,area,zix,ziy,ziz,length
        real(rp):: Invlength
        real(rp):: ksy,ksz,phiy,phiz
        real(rp):: ky1,ky2,ky3,ky4
        real(rp):: kz1,kz2,kz3,kz4
    
        emod = this%m_property(1)   ! E   : Young's modulus
        gmod = this%m_property(2)   ! G   : Shear modulus
        area = this%m_property(3)   ! A   : cross-sectional area
        zix  = this%m_property(6)   ! Jt  : St. Venant torsion constant, used in torsional stiffness G*Jt
        ziy  = this%m_property(7)   ! Iy  : second moment of area about local y-axis
        ziz  = this%m_property(8)   ! Iz  : second moment of area about local z-axis
        ! Note that Jt is generally not equal to Iy+Iz except for circular sections.

        length = this%len0
    
        this%m_stfMat(1:12,1:12)=0.0e0_rp
    
        Invlength = 1.0e0_rp/length
    
        ! Shear correction factors; 5/6 is the rectangular-section default.
        ! Replace with section-specific values, e.g. 6/7 for circular sections.
        ksy = 5.0e0_rp/6.0e0_rp
        ksz = 5.0e0_rp/6.0e0_rp
    
        ! Timoshenko shear parameters
        ! v-theta_z plane: bending about local z, uses Iz and ksy.
        ! w-theta_y plane: bending about local y, uses Iy and ksz.
        ! For an Euler-Bernoulli beam, set phiy = phiz = 0.0e0_rp.
        phiy = 12.0e0_rp*emod*ziz/(ksy*gmod*area*length*length)
        phiz = 12.0e0_rp*emod*ziy/(ksz*gmod*area*length*length)
    
        ! v-theta_z plane, bending about local z, use Iz
        ky1 = 12.0e0_rp*emod*ziz/(length**3*(1.0e0_rp+phiy))
        ky2 =  6.0e0_rp*emod*ziz/(length**2*(1.0e0_rp+phiy))
        ky3 = (4.0e0_rp+phiy)*emod*ziz/(length*(1.0e0_rp+phiy))
        ky4 = (2.0e0_rp-phiy)*emod*ziz/(length*(1.0e0_rp+phiy))
    
        ! w-theta_y plane, bending about local y, use Iy
        kz1 = 12.0e0_rp*emod*ziy/(length**3*(1.0e0_rp+phiz))
        kz2 =  6.0e0_rp*emod*ziy/(length**2*(1.0e0_rp+phiz))
        kz3 = (4.0e0_rp+phiz)*emod*ziy/(length*(1.0e0_rp+phiz))
        kz4 = (2.0e0_rp-phiz)*emod*ziy/(length*(1.0e0_rp+phiz))
    
        ! diagonal terms
        this%m_stfMat(1,1)   = area*emod*Invlength
        this%m_stfMat(2,2)   = ky1
        this%m_stfMat(3,3)   = kz1
        this%m_stfMat(4,4)   = gmod*zix*Invlength
        this%m_stfMat(5,5)   = kz3
        this%m_stfMat(6,6)   = ky3
    
        this%m_stfMat(7,7)   = this%m_stfMat(1,1)
        this%m_stfMat(8,8)   = this%m_stfMat(2,2)
        this%m_stfMat(9,9)   = this%m_stfMat(3,3)
        this%m_stfMat(10,10) = this%m_stfMat(4,4)
        this%m_stfMat(11,11) = this%m_stfMat(5,5)
        this%m_stfMat(12,12) = this%m_stfMat(6,6)
    
        ! upper triangular terms
        this%m_stfMat(1,7)   = -this%m_stfMat(1,1)
    
        this%m_stfMat(2,6)   =  ky2
        this%m_stfMat(2,8)   = -ky1
        this%m_stfMat(2,12)  =  ky2
        this%m_stfMat(6,8)   = -ky2
        this%m_stfMat(6,12)  =  ky4
        this%m_stfMat(8,12)  = -ky2
    
        this%m_stfMat(3,5)   = -kz2
        this%m_stfMat(3,9)   = -kz1
        this%m_stfMat(3,11)  = -kz2
        this%m_stfMat(5,9)   =  kz2
        this%m_stfMat(5,11)  =  kz4
        this%m_stfMat(9,11)  =  kz2
    
        this%m_stfMat(4,10)  = -this%m_stfMat(4,4)
    
        ! symmetric terms
        this%m_stfMat(7,1)   = this%m_stfMat(1,7)
    
        this%m_stfMat(6,2)   = this%m_stfMat(2,6)
        this%m_stfMat(8,2)   = this%m_stfMat(2,8)
        this%m_stfMat(12,2)  = this%m_stfMat(2,12)
        this%m_stfMat(8,6)   = this%m_stfMat(6,8)
        this%m_stfMat(12,6)  = this%m_stfMat(6,12)
        this%m_stfMat(12,8)  = this%m_stfMat(8,12)
    
        this%m_stfMat(5,3)   = this%m_stfMat(3,5)
        this%m_stfMat(9,3)   = this%m_stfMat(3,9)
        this%m_stfMat(11,3)  = this%m_stfMat(3,11)
        this%m_stfMat(9,5)   = this%m_stfMat(5,9)
        this%m_stfMat(11,5)  = this%m_stfMat(5,11)
        this%m_stfMat(11,9)  = this%m_stfMat(9,11)
    
        this%m_stfMat(10,4)  = this%m_stfMat(4,10)
    
        return
    end subroutine Segment_FormStiffMatrix

    subroutine Segment_FormGeomMatrix(this)
        ! ELeMent GEOMetric stiffness matrix for 3D Timoshenko FRaMe
        ! https://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
        ! Henri Gavin, Department of Civil and Environmental Engineering, Duke University
        ! For Euler-Bernoulli :
        ! ISBN 9781441929105 James F. Doyle. P217,228,229,405
        ! ISBN 9780792312086 James F. Doyle. P129,424
        !
        ! Uses the same section convention and shear parameters as Segment_FormStiffMatrix.
        !
        ! DOF order:
        ! [u1,v1,w1,tx1,ty1,tz1,u2,v2,w2,tx2,ty2,tz2]
        implicit none
        class(Segment), intent(inout) :: this
    
        real(rp):: emod,gmod,area,zix,ziy,ziz,length
        real(rp):: s
        real(rp):: ksy,ksz,phiy,phiz
        real(rp):: gy1,gy2,gy3,gy4
        real(rp):: gz1,gz2,gz3,gz4
        real(rp):: gt
    
        s = this%geoFRM

        emod = this%m_property(1)   ! E   : Young's modulus
        gmod = this%m_property(2)   ! G   : Shear modulus
        area = this%m_property(3)   ! A   : cross-sectional area
        zix  = this%m_property(6)   ! Jt  : St. Venant torsion constant, used in torsional stiffness G*Jt
        ziy  = this%m_property(7)   ! Iy  : second moment of area about local y-axis
        ziz  = this%m_property(8)   ! Iz  : second moment of area about local z-axis
        ! Note that Jt is generally not equal to Iy+Iz except for circular sections.

        length = this%len0
    
        ! initialize all geometric stiffness terms to zero
        this%m_geoMat(1:12,1:12) = 0.0e0_rp
    
        ! Shear correction factors; 5/6 is the rectangular-section default.
        ! Replace with section-specific values, e.g. 6/7 for circular sections.
        ksy = 5.0e0_rp/6.0e0_rp
        ksz = 5.0e0_rp/6.0e0_rp
    
        ! Timoshenko shear parameters
        ! v-theta_z plane: bending about local z, uses Iz and ksy.
        ! w-theta_y plane: bending about local y, uses Iy and ksz.
        ! For an Euler-Bernoulli beam, set phiy = phiz = 0.0e0_rp.
        phiy = 12.0e0_rp*emod*ziz/(ksy*gmod*area*length*length)
        phiz = 12.0e0_rp*emod*ziy/(ksz*gmod*area*length*length)
    
        ! ------------------------------------------------------------
        ! v-theta_z plane, DOFs 2,6,8,12
        ! ------------------------------------------------------------
        gy1 = s/length * (6.0e0_rp/5.0e0_rp + 2.0e0_rp*phiy + phiy*phiy) / (1.0e0_rp + phiy)**2
        gy2 = s/length * (length/10.0e0_rp) / (1.0e0_rp + phiy)**2
        gy3 = s/length * (2.0e0_rp*length*length/15.0e0_rp + phiy*length*length/6.0e0_rp + phiy*phiy*length*length/12.0e0_rp) / (1.0e0_rp + phiy)**2
        gy4 = s/length * (-length*length/30.0e0_rp - phiy*length*length/6.0e0_rp - phiy*phiy*length*length/12.0e0_rp) / (1.0e0_rp + phiy)**2
    
        this%m_geoMat(2,2)   =  gy1
        this%m_geoMat(6,6)   =  gy3
        this%m_geoMat(8,8)   =  gy1
        this%m_geoMat(12,12) =  gy3
    
        this%m_geoMat(2,6)   =  gy2
        this%m_geoMat(2,8)   = -gy1
        this%m_geoMat(2,12)  =  gy2
        this%m_geoMat(6,8)   = -gy2
        this%m_geoMat(6,12)  =  gy4
        this%m_geoMat(8,12)  = -gy2
    
        ! ------------------------------------------------------------
        ! w-theta_y plane, DOFs 3,5,9,11
        ! sign convention follows your elastic stiffness matrix
        ! ------------------------------------------------------------
        gz1 = s/length * (6.0e0_rp/5.0e0_rp + 2.0e0_rp*phiz + phiz*phiz) / (1.0e0_rp + phiz)**2
        gz2 = s/length * (length/10.0e0_rp) / (1.0e0_rp + phiz)**2
        gz3 = s/length * (2.0e0_rp*length*length/15.0e0_rp + phiz*length*length/6.0e0_rp + phiz*phiz*length*length/12.0e0_rp) / (1.0e0_rp + phiz)**2
        gz4 = s/length * (-length*length/30.0e0_rp - phiz*length*length/6.0e0_rp - phiz*phiz*length*length/12.0e0_rp) / (1.0e0_rp + phiz)**2
    
        this%m_geoMat(3,3)   =  gz1
        this%m_geoMat(5,5)   =  gz3
        this%m_geoMat(9,9)   =  gz1
        this%m_geoMat(11,11) =  gz3
    
        this%m_geoMat(3,5)   = -gz2
        this%m_geoMat(3,9)   = -gz1
        this%m_geoMat(3,11)  = -gz2
        this%m_geoMat(5,9)   =  gz2
        this%m_geoMat(5,11)  =  gz4
        this%m_geoMat(9,11)  =  gz2
    
        ! ------------------------------------------------------------
        ! torsional geometric stiffness
        ! standard Timoshenko/frame form
        ! ------------------------------------------------------------
        gt = s*zix/(area*length)
    
        this%m_geoMat(4,4)   =  gt
        this%m_geoMat(4,10)  = -gt
        this%m_geoMat(10,10) =  gt
    
        ! ------------------------------------------------------------
        ! symmetric terms
        ! ------------------------------------------------------------
        this%m_geoMat(6,2)   = this%m_geoMat(2,6)
        this%m_geoMat(8,2)   = this%m_geoMat(2,8)
        this%m_geoMat(12,2)  = this%m_geoMat(2,12)
        this%m_geoMat(8,6)   = this%m_geoMat(6,8)
        this%m_geoMat(12,6)  = this%m_geoMat(6,12)
        this%m_geoMat(12,8)  = this%m_geoMat(8,12)
    
        this%m_geoMat(5,3)   = this%m_geoMat(3,5)
        this%m_geoMat(9,3)   = this%m_geoMat(3,9)
        this%m_geoMat(11,3)  = this%m_geoMat(3,11)
        this%m_geoMat(9,5)   = this%m_geoMat(5,9)
        this%m_geoMat(11,5)  = this%m_geoMat(5,11)
        this%m_geoMat(11,9)  = this%m_geoMat(9,11)
    
        this%m_geoMat(10,4)  = this%m_geoMat(4,10)
    
        return
    end subroutine Segment_FormGeomMatrix

    subroutine Segment_InitTriad_D(this)
        implicit none
        class(Segment), intent(inout) :: this
        call this%BuildAxisDirTriad(this%xll0,this%xmm0,this%xnn0,this%dirc0,this%triad_n1)
        ! all element triads have same initial orientation
        this%triad_n2(1:3,1:3)=this%triad_n1(1:3,1:3)
        this%triad_ee(1:3,1:3)=this%triad_n1(1:3,1:3)
        return
    end subroutine Segment_InitTriad_D

    subroutine Segment_RigidUpdateTriad_D(this)
        ! Used for rigid-body triad updates.
        ! Rebuild the full triad from the current beam axis and the current material direction.
        implicit none
        class(Segment), intent(inout) :: this
        call this%BuildAxisDirTriad(this%xll1,this%xmm1,this%xnn1,this%dirc1,this%triad_n1)
        ! In rigid-body motion, the two nodal triads and the element triad share the same current orientation.
        this%triad_n2(1:3,1:3)=this%triad_n1(1:3,1:3)
        this%triad_ee(1:3,1:3)=this%triad_n1(1:3,1:3)
        return
    end subroutine Segment_RigidUpdateTriad_D

    subroutine Segment_BuildAxisDirTriad(this,l,m,n,d,triad)
        ! If span direction vector is not 0
        ! ex: beam axis direction
        ! ey: span direction
        ! ez: ex cross product ey
        ! If span direction vector is 0, initial triad method follows Doyle
        ! Based on the rotation matrix [R]
        ! [triad] = [R]^T
        ! ISBN 9780792312086 James F. Doyle. P157,158
        implicit none
        class(Segment), intent(in) :: this
        real(rp), intent(in) :: l,m,n,d(3)
        real(rp), intent(out) :: triad(3,3)
        real(rp) :: ex(3), ey(3), ez(3), dir(3), dd, proj

        ! ex: beam axis direction
        ex(1) = l
        ex(2) = m
        ex(3) = n
        ! Normalization
        dd = sqrt(dot_product(ex, ex))
        if (dd .gt. 1.0e-14_rp) then
            ex = ex / dd
        endif

        dir(1:3) = d(1:3)
        ! Normalization
        dd = sqrt(dot_product(dir, dir))
        if (dd .gt. 1.0e-14_rp) then
            dir = dir / dd
        endif
        proj = dot_product(dir, ex)
        ! ey: span direction
        ! Gram-Schmidt
        ey = dir - proj * ex
        dd = sqrt(dot_product(ey, ey))

        ! if no span direction, use 
        if (dd .le. 1.0e-10_rp) then
            if (abs(ex(3)) .gt. 0.995e0_rp) then
                ey(1) = 0.0e0_rp
                ey(2) = 1.0e0_rp
                ey(3) =  0.0e0_rp
                ez(1) = -ex(3)
                ez(2) = 0.0e0_rp
                ez(3) =  0.0e0_rp
            else
                dd = sqrt(ex(1)*ex(1)+ex(2)*ex(2))
                ey(1) = -ex(2)/dd
                ey(2) =  ex(1)/dd
                ey(3) =  0.0e0_rp
                ez(1) = -ex(1)*ex(3)/dd
                ez(2) = -ex(2)*ex(3)/dd
                ez(3) =  dd
            endif
        else
            ! Normalization
            ey = ey / dd
            ! ex cross product ey
            ez(1) = ex(2)*ey(3) - ex(3)*ey(2)
            ez(2) = ex(3)*ey(1) - ex(1)*ey(3)
            ez(3) = ex(1)*ey(2) - ex(2)*ey(1)
            ! Normalization
            dd = sqrt(dot_product(ez, ez))
            if (dd .gt. 1.0e-14_rp) then
                ez = ez / dd
            endif
        endif

        triad(1:3,1) = ex(1:3)
        triad(1:3,2) = ey(1:3)
        triad(1:3,3) = ez(1:3)
        return
    end subroutine Segment_BuildAxisDirTriad

    subroutine Segment_RotateMatrix(this)
        ! ISBN 9780792312086 James F. Doyle. P157,158
        ! Based on the element triad [triad_ee]
        ! [R] = [triad_ee]^T
        implicit none
        class(Segment), intent(inout) :: this
        this%m_rotMat(1:3,1:3) = transpose(this%triad_ee(1:3,1:3))
        return
    end subroutine Segment_RotateMatrix

    subroutine Segment_RKR(this,ek)
        implicit none
        class(Segment), intent(inout) :: this
        real(rp), intent(inout) :: ek(12,12)
        real(rp):: r(3,3),rt(3,3),ktemp(12,12)
        integer:: i,j,k,j1,j2,ii,jj,in,jn
        r = 0.0e0_rp
        r = this%m_rotMat
        do  in=1,3
        do  jn=1,3
            rt(jn,in)=this%m_rotMat(in,jn)
        enddo
        enddo
        ! take [Rtrans][K][R] using the nature of [R] to speed computation.
        ! k is sectioned off into 3x3s then multiplied [rtrans][k][r]
        !
        do  i=0,3
        do  j=0,3
            j1=i*3
            j2=j*3
            do    k=1,3
            do    ii=1,3
                ktemp(j1+k,j2+ii)=0.0e0_rp
                do     jj=1,3
                ktemp(j1+k,j2+ii)=ktemp(j1+k,j2+ii)+ek(j1+k,j2+jj)*r(jj,ii)
                enddo
            enddo
            enddo
            do  k=1,3
            do  ii=1,3
                ek(j1+k,j2+ii)=0.0e0_rp
                do  jj=1,3
                    ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rt(k,jj)*ktemp(j1+jj,j2+ii)
                enddo
            enddo
            enddo
        enddo
        enddo
        return
    end subroutine Segment_RKR

    subroutine Segment_BodyStress_D(this, lodInte, gEQ)
        ! Element nodal force
        ! Calculate the lateral buckling instability caused by axial force
        ! https://people.duke.edu/~hpgavin/cee421/frame-finite-def.pdf
        ! Henri Gavin, Department of Civil and Environmental Engineering, Duke University
        ! ISBN 9781441929105 James F. Doyle. P214,P353
        implicit none
        class(Segment), intent(inout) :: this
        integer, intent(in) :: gEQ
        real(rp), intent(inout) :: lodInte(1:gEQ)
        real(rp) :: triad_00(3,3),triad_11(3,3),triad_22(3,3)
        real(rp) :: ub(12),dl
        real(rp) :: rr(3,3)
        real(rp) :: force(12),forceb(12)
        real(rp) :: du,dv,dw
        real(rp) :: tx1,tx2,ty1,ty2,tz1,tz2,tx,ty,tz
        real(rp) :: fxx,emod,area
        integer :: i,j

        du = this%dx1 - this%dx0
        dv = this%dy1 - this%dy0
        dw = this%dz1 - this%dz0
        ! Define the the xyz directional distance of the two nodes of the beam element after deformation
        ! as well as the length, xl, of the element after deformation
        dl = ( (this%dx0+this%dx1)*du +(this%dy0+this%dy1)*dv +(this%dz0+this%dz1)*dw )/ (this%len0+this%len1)
        ! get twisting angles
        do    i=1,3
        do    j=1,3
            triad_00(i,j)=this%triad_ee(i,j)
            triad_11(i,j)=this%triad_ee(i,j)
            triad_22(i,j)=this%triad_n1(i,j)
        enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        call Segment_global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)
        !
        do    i=1,3
        do    j=1,3
            triad_11(i,j)=this%triad_ee(i,j)
            triad_22(i,j)=this%triad_n2(i,j)
        enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        call Segment_global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

        ! non-zero ty1 tz1 u2 tx2 ty2 tz2
        ub(1)=0.0e0_rp
        ub(2)=0.0e0_rp
        ub(3)=0.0e0_rp
        ub(4)=tx1
        ub(5)=ty1
        ub(6)=tz1
        !
        ub(7)=dl
        ub(8)=0.0e0_rp
        ub(9)=0.0e0_rp
        ub(10)=tx2
        ub(11)=ty2
        ub(12)=tz2
        !
        ! compute axial force
        emod = this%m_property(1)
        area = this%m_property(3)
        fxx=dl*emod*area/this%len0
        ! save local force for geo stiff
        ! write(igeoFRM,'(D25.15)') fxx
        this%geoFRM=fxx
        ! nodal forces in local coords
        ! {F}=[k]{u}
        forceb(1:12) = matmul(this%m_stfMat,ub)
        ! transform to global
        do  i=1,3
        do  j=1,3
            rr(i,j)=this%triad_ee(i,j)
        enddo
        enddo
        do  i=1,nElmtDofs
            force(i)=0.0e0_rp
        enddo
        do    i=1,3
        do    j=1,3
            force(0+i) = force(0+i) + rr(i,j)*forceb(0+j)
            force(3+i) = force(3+i) + rr(i,j)*forceb(3+j)
            force(6+i) = force(6+i) + rr(i,j)*forceb(6+j)
            force(9+i) = force(9+i) + rr(i,j)*forceb(9+j)
        enddo
        enddo
        call this%LocToGlobal(force, lodInte,  gEQ)
        return
    end subroutine Segment_BodyStress_D

    subroutine Segment_StrainEnergy_D(this)
        implicit none
        class(Segment), intent(inout) :: this
        real(rp) :: Strech(nElmtDofs,nElmtDofs),BendTor(nElmtDofs,nElmtDofs)
        real(rp) :: triad_00(3,3),triad_11(3,3),triad_22(3,3)
        real(rp) :: ub(12),dl
        real(rp) :: du,dv,dw
        real(rp) :: tx1,tx2,ty1,ty2,tz1,tz2,tx,ty,tz
        integer :: i,j

        du = this%dx1 - this%dx0
        dv = this%dy1 - this%dy0
        dw = this%dz1 - this%dz0
        ! Define the the xyz directional distance of the two nodes of the beam element after deformation
        ! as well as the length, xl, of the element after deformation
        dl = ( (this%dx0+this%dx1)*du +(this%dy0+this%dy1)*dv +(this%dz0+this%dz1)*dw )/ (this%len0+this%len1)
        ! get twisting angles
        do    i=1,3
        do    j=1,3
            triad_00(i,j)=this%triad_ee(i,j)
            triad_11(i,j)=this%triad_ee(i,j)
            triad_22(i,j)=this%triad_n1(i,j)
        enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        call Segment_global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)
        !
        do    i=1,3
        do    j=1,3
            triad_11(i,j)=this%triad_ee(i,j)
            triad_22(i,j)=this%triad_n2(i,j)
        enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        call Segment_global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

        ! non-zero ty1 tz1 u2 tx2 ty2 tz2
        ub(1)=0.0e0_rp
        ub(2)=0.0e0_rp
        ub(3)=0.0e0_rp
        ub(4)=tx1
        ub(5)=ty1
        ub(6)=tz1
        !
        ub(7)=dl
        ub(8)=0.0e0_rp
        ub(9)=0.0e0_rp
        ub(10)=tx2
        ub(11)=ty2
        ub(12)=tz2

        Strech(:,:)=0.0e0_rp
        Strech(1,1)=this%m_stfMat(1,1)
        Strech(1,7)=this%m_stfMat(1,7)
        Strech(7,7)=this%m_stfMat(7,7)
        Strech(7,1)=this%m_stfMat(7,1)

        BendTor(:,:)=this%m_stfMat(:,:)
        BendTor(1,1)=0.0e0_rp
        BendTor(1,7)=0.0e0_rp
        BendTor(7,7)=0.0e0_rp
        BendTor(7,1)=0.0e0_rp

        ! nodal strain in local coords
        this%strainEnergy(1)=0.5e0_rp*sum(matmul(Strech,ub)*ub)
        this%strainEnergy(2)=0.5e0_rp*sum(matmul(BendTor,ub)*ub)
    end subroutine Segment_StrainEnergy_D

    subroutine Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        ! GET ANGLE of between TRIADs
        ! ISBN 9781441929105 James F. Doyle. P186,187 Equ.(3.7) (3.8)
        implicit none
        real(rp):: triad_11(3,3),triad_22(3,3)
        real(rp):: rr(3,3)
        real(rp):: tx,ty,tz, dtx,dty,dtz,theta,sint,trace_rr,factor
        integer:: i,j,k
        !
        ! get angle between two triads
        do    i=1,3
            do    j=1,3
                rr(i,j)=0.0e0_rp
                do    k=1,3
                    rr(i,j)=rr(i,j) + triad_22(i,k)*triad_11(j,k)
                enddo
            enddo
        enddo

        dtx = (rr(3,2)-rr(2,3))/2.0e0_rp
        dty = (rr(1,3)-rr(3,1))/2.0e0_rp
        dtz = (rr(2,1)-rr(1,2))/2.0e0_rp

        trace_rr = rr(1,1) + rr(2,2) + rr(3,3)
        trace_rr = (trace_rr - 1.0e0_rp)/2.0e0_rp
        if (trace_rr .gt. 1.0e0_rp) trace_rr = 1.0e0_rp
        if (trace_rr .lt. -1.0e0_rp) trace_rr = -1.0e0_rp

        sint = sqrt(dtx*dtx+dty*dty+dtz*dtz)
        theta = acos(trace_rr)

        if (sint .lt. 1.0e-10_rp .or. theta .lt. 1.0e-10_rp) then
            tx = dtx
            ty = dty
            tz = dtz
        else
            factor = theta/sint
            tx = factor*dtx
            ty = factor*dty
            tz = factor*dtz
        endif

        return
    end subroutine Segment_get_angle_triad

    subroutine Segment_global_to_local(triad,tx,ty,tz,tx2,ty2,tz2)
        implicit none
        real(rp):: triad(3,3)
        real(rp):: tx,ty,tz,tx2,ty2,tz2
        ! [R] = [triad]^T
        ! [g_angle] = [R]*[l_angle]
        tx2 = triad(1,1)*tx+triad(2,1)*ty+triad(3,1)*tz
        ty2 = triad(1,2)*tx+triad(2,2)*ty+triad(3,2)*tz
        tz2 = triad(1,3)*tx+triad(2,3)*ty+triad(3,3)*tz
        return
    end subroutine Segment_global_to_local

    subroutine Segment_UpdateTriad_D(this,dspnn,nND)
        ! update angle of triads
        ! Finite Rotation(Rodrique's formula) 
        ! ISBN 9781441929105 James F. Doyle. P184 Equ.(3.4)
        ! triad_n1: the triad of node 1 in beam element
        ! triad_n2: the triad of node 2 in beam element
        ! For each beam element, calculate the current orientation triad: [triad] = [R][triad]
        implicit none
        class(Segment), intent(inout) :: this
        integer, intent(in) :: nND
        real(rp):: dspnn(1:6,1:nND)
        real(rp):: rr(3,3)
        real(rp):: dtx1,dty1,dtz1
        integer:: node0,node1

        node0 = this%node0
        node1 = this%node1
        ! n1 node
        dtx1=dspnn(4,node0)
        dty1=dspnn(5,node0)
        dtz1=dspnn(6,node0)
        call Segment_FiniteRot(dtx1,dty1,dtz1,rr)
        this%triad_n1=matmul(rr,this%triad_n1)
        ! n2 node
        dtx1=dspnn(4,node1)
        dty1=dspnn(5,node1)
        dtz1=dspnn(6,node1)
        call Segment_FiniteRot(dtx1,dty1,dtz1,rr)
        this%triad_n2=matmul(rr,this%triad_n2)

        return
    end subroutine Segment_UpdateTriad_D

    subroutine Segment_MakeTriad_ee(this)
        ! Get orientation of element
        ! ISBN 9781441929105 James F. Doyle. P189-191
        ! triad_n1: the triad of node 1 in beam element
        ! triad_n2: the triad of node 2 in beam element
        ! triad_ee: the triad of beam element(beam center)
        implicit none
        class(Segment), intent(inout) :: this
        real(rp):: triad_aa(3,3)
        real(rp):: rr(3,3),tx,ty,tz
        real(rp):: triad_11(3,3),triad_22(3,3)
        real(rp):: xll,xmm,xnn,dd,r2e1,r3e1
        integer:: i,j,k

        xll=this%xll1
        xmm=this%xmm1
        xnn=this%xnn1
        dd =sqrt(xll*xll+xmm*xmm)
        do    i=1,3
            do    j=1,3
                this%triad_ee(i,j)=0.0e0_rp
            enddo
        enddo
        this%triad_ee(1,1)=xll
        this%triad_ee(2,1)=xmm
        this%triad_ee(3,1)=xnn
        !
        ! get angle between two triads
        do    i=1,3
            do    j=1,3
                triad_11(i,j)=this%triad_n1(i,j)
                triad_22(i,j)=this%triad_n2(i,j)
            enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        !
        ! rotate n1 to intermediate
        tx=tx/2.0e0_rp
        ty=ty/2.0e0_rp
        tz=tz/2.0e0_rp
        call Segment_FiniteRot(tx,ty,tz,rr)
        triad_aa(1:3,1:3)=matmul(rr(1:3,1:3),this%triad_n1(1:3,1:3))
        !
        ! vectors e2 e3
        r2e1 = 0.0e0_rp
        r3e1 = 0.0e0_rp
        do    k=1,3
            r2e1 = r2e1 + triad_aa(k,2)*this%triad_ee(k,1)
            r3e1 = r3e1 + triad_aa(k,3)*this%triad_ee(k,1)
        enddo
        do    j=1,3
            this%triad_ee(j,2)=triad_aa(j,2) - r2e1*(triad_aa(j,1)+this%triad_ee(j,1))/2.0e0_rp
            this%triad_ee(j,3)=triad_aa(j,3) - r3e1*(triad_aa(j,1)+this%triad_ee(j,1))/2.0e0_rp
        enddo
        !
        ! Gram-Schmidt
        this%triad_ee(:,2) = this%triad_ee(:,2) - dot_product(this%triad_ee(:,2), this%triad_ee(:,1)) * this%triad_ee(:,1)
        dd = sqrt(dot_product(this%triad_ee(:,2), this%triad_ee(:,2)))
        if (dd > 1.0e-14_rp) this%triad_ee(:,2) = this%triad_ee(:,2) / dd
        !
        ! e3 = e1 × e2
        this%triad_ee(:,3) = (/ this%triad_ee(2,1)*this%triad_ee(3,2)-this%triad_ee(3,1)*this%triad_ee(2,2), &
                                this%triad_ee(3,1)*this%triad_ee(1,2)-this%triad_ee(1,1)*this%triad_ee(3,2), &
                                this%triad_ee(1,1)*this%triad_ee(2,2)-this%triad_ee(2,1)*this%triad_ee(1,2) /)
        return
    end subroutine Segment_MakeTriad_ee

    subroutine Segment_FiniteRot(t1,t2,t3,rr)
        ! Finite rotation(Rodrique's formula)
        ! ISBN 9781441929105 James F. Doyle. P184 Equ.(3.4)
        implicit none
        real(rp):: rr(3,3),rr1(3,3),rr2(3,3),rr3(3,3)
        real(rp):: t1,t2,t3,tt,ss,cc,c1,c2
        !
        integer:: i,j
        !
        tt=sqrt( t1**2 + t2**2 + t3**2 )
        ss=sin(tt)
        cc=cos(tt)
        !
        rr1(1,1)=1.0e0_rp
        rr1(1,2)=0.0e0_rp
        rr1(1,3)=0.0e0_rp
        rr1(2,1)=0.0e0_rp
        rr1(2,2)=1.0e0_rp
        rr1(2,3)=0.0e0_rp
        rr1(3,1)=0.0e0_rp
        rr1(3,2)=0.0e0_rp
        rr1(3,3)=1.0e0_rp
        !
        rr2(1,1)=0.0e0_rp
        rr2(1,2)=-t3
        rr2(1,3)= t2
        rr2(2,1)= t3
        rr2(2,2)=0.0e0_rp
        rr2(2,3)=-t1
        rr2(3,1)=-t2
        rr2(3,2)= t1
        rr2(3,3)=0e0_rp
        !
        rr3(1,1)=-t3*t3-t2*t2
        rr3(1,2)= t2*t1
        rr3(1,3)= t3*t1
        rr3(2,1)= t1*t2
        rr3(2,2)=-t3*t3-t1*t1
        rr3(2,3)= t3*t2
        rr3(3,1)= t1*t3
        rr3(3,2)= t2*t3
        rr3(3,3)=-t2*t2-t1*t1
        !
        do    i=1,3
        do    j=1,3
            if    (tt .lt. 1.0e-10_rp) then
                c1=1.0e0_rp
                ! c2=1.0
                c2=0.5e0_rp
            else
                c1 = ss/tt
                c2 = (1.0e0_rp-cc)/tt**2
            endif
            rr(i,j) = rr1(i,j) + rr2(i,j)*c1 + rr3(i,j)*c2
        enddo
        enddo
        !
        return
    end subroutine Segment_FiniteRot

    subroutine Segment_MapReferencePosToCurrent(this, coordsOut, TTT, XYZ, AoA)
        implicit none
        class(Segment), intent(in) :: this
        real(rp), intent(out) :: coordsOut(12)
        real(rp), intent(in) :: TTT(3,3)
        real(rp), intent(in) :: XYZ(3)
        real(rp), intent(in) :: AoA(3)
        coordsOut(1:3) = matmul(TTT, this%x00(1:3)) + XYZ
        coordsOut(4:6) = AoA
        coordsOut(7:9) = matmul(TTT, this%x00(7:9)) + XYZ
        coordsOut(10:12) = AoA
    end subroutine Segment_MapReferencePosToCurrent

    subroutine Segment_MapReferenceDirToCurrent(this, dirc, TTT)
        implicit none
        class(Segment), intent(in) :: this
        real(rp), intent(out) :: dirc(3)
        real(rp), intent(in) :: TTT(3,3)
        dirc = matmul(TTT, this%dirc00)
    end subroutine Segment_MapReferenceDirToCurrent

end module SegmentStructure


module SolidSolver
    use SegmentStructure
    use ConstParams, only: rp
    implicit none
    private
    integer, parameter:: m_idat=12
    real(rp):: m_dampK,m_dampM,m_GeoGamma,m_NewmarkGamma,m_NewmarkBeta
    real(rp):: m_dtolFEM,m_pi
    integer:: m_ntolFEM,m_isKB
    real(rp):: m_g(3)
    public :: BeamSolver,Set_SolidSolver_Params
    type :: BeamSolver
        type(Segment), allocatable :: m_elements(:)

        integer :: nND,nEL,nMT,gEQ
        real(rp) :: FishInfo(3)
        real(rp), allocatable :: pos(:,:), dsp(:,:), vel(:,:), acc(:,:), mss(:,:)
        real(rp), allocatable :: lodInte(:), lodExte(:), lodEffe(:), lodFlow(:), lodRepl(:), lodGrav(:), vBC(:)
        real(rp) :: coeffs(0:7)

        character (LEN=40):: FEmeshName
        integer:: iBodyModel
        real(rp):: Freq,denR,KB,KS,EmR,psR,tcR,St
        real(rp):: elmax,elmin
        real(rp):: XYZ(3),XYZo(3),initXYZVel(3),XYZAmpl(3),XYZPhi(3),XYZd(3),UVW(3)
        real(rp):: AoA(3),AoAo(3),AoAAmpl(3),AoAPhi(3),AoAd(3),WWW1(3),WWW2(3),WWW3(3)
        real(rp):: TTT00(3,3),TTT0(3,3),TTTnxt(3,3)
        integer:: isMotionGiven(6)
    contains
        procedure :: SetSolver => Beam_SetSolver
        procedure :: ReadBuild => Beam_ReadBuild
        procedure :: Initialise => Beam_Initialise
        procedure :: InitLoad => Beam_InitLoad
        procedure :: InitPosDspVelAcc => Beam_InitPosDspVelAcc
        procedure :: UpdateVelFromPosAngular => Beam_UpdateVelFromPosAngular
        procedure :: InitTriadANDFormMass => Beam_InitTriadANDFormMass
        procedure :: InitGrav => Beam_InitGrav
        procedure :: UpdateStrainEnergy => Beam_UpdateStrainEnergy
        procedure :: UpdateNewmarkCoeffs => Beam_UpdateNewmarkCoeffs
        procedure :: Solver => Beam_Solver
        procedure :: InitDspVelAccATTimeT => Beam_InitDspVelAccATTimeT
        procedure :: UpdateMatrixANDLoad => Beam_UpdateMatrixANDLoad
        procedure :: CG_Solve => Beam_CG_Solve
        procedure :: MatrixMultipy => Beam_MatrixMultipy
        procedure :: BoundaryCond => Beam_BoundaryCond
        procedure :: ApplyFixedZero => Beam_ApplyFixedZero
        procedure :: ApplyFixedValue => Beam_ApplyFixedValue
        procedure :: InitPreconditioner => Beam_InitPreconditioner
        procedure :: ApplyPreconditioner => Beam_ApplyPreconditioner
        procedure :: ScalarPreconditioned => Beam_ScalarPreconditioned
        procedure :: BlockPreconditioned => Beam_BlockPreconditioned
        procedure :: InitScalarJacobiPreconditioner => Beam_InitScalarJacobiPreconditioner
        procedure :: ApplyScalarJacobiPreconditioner => Beam_ApplyScalarJacobiPreconditioner
        procedure :: InitBlockJacobiPreconditioner => Beam_InitBlockJacobiPreconditioner
        procedure :: ApplyBlockJacobiPreconditioner => Beam_ApplyBlockJacobiPreconditioner
        procedure :: CheckCGBreakdown => Beam_CheckCGBreakdown
        procedure :: UpdateDspANDTride => Beam_UpdateDspANDTride
        procedure :: UpdateVelAcc => Beam_UpdateVelAcc
        procedure :: UpdateIterInfo => Beam_UpdateIterInfo
        procedure :: calculate_angle_material => Beam_calculate_angle_material
        procedure :: write_solid => Beam_write_solid
        procedure :: write_solid_temp => Beam_write_solid_temp
        procedure :: read_solid_temp => Beam_read_solid_temp
        procedure :: write_solid_params => Beam_write_solid_params
        procedure :: write_solid_materials => Beam_write_solid_materials
        procedure :: write_solid_info => Beam_write_solid_info
        procedure :: write_solid_probes => Beam_write_solid_probes
        procedure :: structure => Beam_structure
    end type BeamSolver
  contains
    subroutine Beam_SetSolver(this,FEmeshName_,iBodyModel_,isMotionGiven_,denR_,KB_,KS_,EmR_,psR_,tcR_,St_, &
                            Freq_,initXYZVel_,XYZo_,XYZAmpl_,XYZPhi_,AoAo_,AoAAmpl_,AoAPhi_)
        implicit none
        class(BeamSolver), intent(inout) :: this
        character (LEN=40),intent(in):: FEmeshName_
        integer,intent(in):: iBodyModel_,isMotionGiven_(6)
        real(rp),intent(in):: denR_,KB_,KS_,EmR_,psR_,tcR_,St_
        real(rp),intent(in):: Freq_
        real(rp),intent(in):: initXYZVel_(3),XYZo_(3),XYZAmpl_(3),XYZPhi_(3)
        real(rp),intent(in):: AoAo_(3),AoAAmpl_(3),AoAPhi_(3)
        
        this%FEmeshName = FEmeshName_
        this%iBodyModel = iBodyModel_
        this%isMotionGiven(1:6)=isMotionGiven_(1:6)
        this%denR= denR_
        this%psR = psR_
        this%EmR = EmR_
        this%tcR = tcR_
        this%KB  = KB_
        this%KS  = KS_
        this%Freq=Freq_
        this%St  =St_
        this%XYZo(1:3)     = XYZo_(1:3)
        this%initXYZVel(1:3) = initXYZVel_(1:3)
        this%XYZAmpl(1:3)  = XYZAmpl_(1:3)
        this%XYZPhi(1:3)   = XYZPhi_(1:3)
        this%AoAo(1:3)     = AoAo_(1:3)
        this%AoAAmpl(1:3)  = AoAAmpl_(1:3)
        this%AoAPhi(1:3)   = AoAPhi_(1:3)
    end subroutine Beam_SetSolver
    subroutine Set_SolidSolver_Params(dampK,dampM,GeoGamma,NewmarkGamma,NewmarkBeta,dtolFEM,ntolFEM,isKB,g)
        implicit none
        real(rp),intent(in):: dampK,dampM,GeoGamma,NewmarkGamma,NewmarkBeta
        real(rp),intent(in):: dtolFEM,g(3)
        integer,intent(in):: ntolFEM,isKB
        m_dampK = dampK
        m_dampM = dampM
        m_GeoGamma = GeoGamma
        m_NewmarkGamma = NewmarkGamma
        m_NewmarkBeta = NewmarkBeta
        m_dtolFEM = dtolFEM
        m_ntolFEM = ntolFEM
        m_pi = 3.141592653589793e0_rp
        m_isKB = isKB
        m_g = g
    end subroutine Set_SolidSolver_Params

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ structural Datafile
!   Allocate memory for solid simulation
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Beam_ReadBuild(this,nAsfac,nLchod)
        implicit none
        class(BeamSolver), intent(inout) :: this
        character (LEN=40):: filename
        integer :: fileiD = 996
        character (LEN=1000):: buffer
        real(rp), allocatable :: xyz(:, :), material(:, :)
        integer, allocatable :: boundary(:, :)
        real(rp), intent(out):: nAsfac,nLchod
        ! load mesh information
        filename = this%FEmeshName
        open(unit=fileiD, file = filename )
            read(fileiD,*) buffer
            read(fileiD,*) this%nND, this%nEL, this%nMT
            this%gEQ = this%nND * 6
            allocate(this%vBC(1:this%gEQ))
            ! load points data
            allocate(xyz(1:8, 1:this%nND))
            call Beam_ReadPoints()
            ! load matieral data
            ! property data will be overwritten if isKB = 0 or 1
            allocate(material(1:8, 1:this%nMT))
            call Beam_ReadMaterials()
            ! load boundary condition
            allocate(boundary(1:6, 1:this%nND))
            call Beam_ReadBoundary()
            ! load and build elements
            allocate(this%m_elements(1:this%nEL))
            call Beam_ReadBuildElements()
        close(fileiD)
        call Beam_cptAsfac()
        call Beam_adjustBC()

    contains

        subroutine Beam_ReadPoints()
            implicit none
            integer :: tmpid, i, nCheck
            call FindSection('POINT')
            read(fileiD,*) nCheck
            if (nCheck /= this%nND) stop 'ERROR: POINT number inconsistent'
            do i = 1,this%nND
                read(fileiD,*) tmpid, xyz(1:8,i)
            enddo
        end subroutine Beam_ReadPoints

        subroutine Beam_ReadMaterials()
            implicit none
            integer :: tmpid, i, nCheck
            call FindSection('MATERIAL')
            read(fileiD,*) nCheck
            if (nCheck /= this%nMT) stop 'ERROR: MATERIAL number inconsistent'
            do i = 1,this%nMT
                read(fileiD,*) tmpid, material(1:8,i)
            enddo
        end subroutine Beam_ReadMaterials

        subroutine Beam_ReadBoundary()
            implicit none
            integer :: tmpid, i, nCheck
            call FindSection('CONSTRAINT')
            read(fileiD,*) nCheck
            if (nCheck /= this%nND) stop 'ERROR: CONSTRAINT number inconsistent'
            do i = 1,this%nND
                read(fileiD,*) tmpid, boundary(1:6,i)
            enddo
        end subroutine Beam_ReadBoundary

        subroutine Beam_ReadBuildElements()
            implicit none
            integer :: tmpid, n, i, j, k, imat, itype, Nspan, nCheck
            call FindSection('ELEMENT')
            read(fileiD,*) nCheck
            if (nCheck /= this%nEL) stop 'ERROR: ELEMENT number inconsistent'
            do n=1,this%nEL
                read(fileiD,*) tmpid, i, j, k, itype, imat, Nspan
                if(1.le.tmpid .and. tmpid.le.this%nEL) then
                    call this%m_elements(tmpid)%Build(i, j, this%nND, itype, Nspan, xyz, material(1:8, imat), boundary)
                endif
            enddo
        end subroutine Beam_ReadBuildElements

        subroutine FindSection(sectionName)
            implicit none
            character(len=*), intent(in) :: sectionName
            character(len=1000) :: line
            integer :: ios
        
            rewind(fileiD)
        
            do
                read(fileiD, *, iostat=ios) line
                if (ios /= 0) exit
        
                if (trim(adjustl(line)) == trim(sectionName)) return
            enddo
        
            write(*,*) 'ERROR: cannot find section ', trim(sectionName), ' in ', trim(filename)
            stop
        end subroutine FindSection

        subroutine Beam_cptAsfac()
            implicit none
            real(rp) :: lentemp
            integer :: i
            nAsfac = sum([(this%m_elements(i)%areaElem00, i=1,this%nEL)])
            this%elmax = maxval([(this%m_elements(i)%areaElem00, i=1,this%nEL)])
            this%elmin = minval([(this%m_elements(i)%areaElem00, i=1,this%nEL)])
            !calculate spanwise length, chord length, aspect ratio
            nLchod  = maxval([(this%m_elements(i)%x00(1), this%m_elements(i)%x00(7), i=1,this%nEL)]) - &
                      minval([(this%m_elements(i)%x00(1), this%m_elements(i)%x00(7), i=1,this%nEL)])
            lentemp = maxval([(this%m_elements(i)%x00(2), this%m_elements(i)%x00(8), i=1,this%nEL)]) - &
                      minval([(this%m_elements(i)%x00(2), this%m_elements(i)%x00(8), i=1,this%nEL)])
            if(lentemp .gt. nLchod) nLchod = lentemp
        end subroutine

        subroutine Beam_adjustBC()
            implicit none
            integer :: n
            do n = 1, this%nEL
                if (this%m_elements(n)%bc(1) == 1) then
                    this%m_elements(n)%bc(1:6) = this%isMotionGiven(1:6)
                end if
                if (this%m_elements(n)%bc(7) == 1) then
                    this%m_elements(n)%bc(7:12) = this%isMotionGiven(1:6)
                end if
            enddo
        end subroutine

    end subroutine Beam_ReadBuild

    subroutine Beam_Initialise(this,time)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(in):: time
        integer:: n
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   initialize solid field
    !   compute initial values
    !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        this%TTT00(:,:)=0.0e0_rp
        this%TTT00(1,1)=1.0e0_rp
        this%TTT00(2,2)=1.0e0_rp
        this%TTT00(3,3)=1.0e0_rp

        this%XYZ(1:3)=this%XYZo(1:3)+this%XYZAmpl(1:3)*cos(2.0e0_rp*m_pi*this%Freq*time+this%XYZPhi(1:3)) + this%initXYZVel(1:3)*time
        this%AoA(1:3)=this%AoAo(1:3)+this%AoAAmpl(1:3)*cos(2.0e0_rp*m_pi*this%Freq*time+this%AoAPhi(1:3))

        call AoAtoTTT(this%AoA(1:3),this%TTT0(1:3,1:3))
        call AoAtoTTT(this%AoA(1:3),this%TTTnxt(1:3,1:3))
        call Segment_get_angle_triad(this%TTT0(1:3,1:3),this%TTTnxt(1:3,1:3),this%AoAd(1),this%AoAd(2),this%AoAd(3))

        do n=1,this%nEL
            call this%m_elements(n)%MapReferencePosToCurrent(this%m_elements(n)%x0, this%TTT0, this%XYZ, this%AoAd)
            call this%m_elements(n)%MapReferenceDirToCurrent(this%m_elements(n)%dirc0, this%TTT0)
            call this%m_elements(n)%Init()
        enddo

        call this%InitLoad()
        call this%InitPosDspVelAcc()
        
        if(this%iBodyModel.eq.1)then
            this%UVW(1:3) =-2.0e0_rp*m_pi*this%Freq*this%XYZAmpl(1:3)*sin(2.0e0_rp*m_pi*this%Freq*time+this%XYZPhi(1:3)) + this%initXYZVel(1:3) !time=0
            !rotational velocity
            this%WWW1(1:3)=-2.0e0_rp*m_pi*this%Freq*this%AoAAmpl(1:3)*sin(2.0e0_rp*m_pi*this%Freq*time+this%AoAPhi(1:3))
            this%WWW2(1:3)=[this%WWW1(1)*cos(this%AoA(2))+this%WWW1(3),    &
                            this%WWW1(1)*sin(this%AoA(2))*sin(this%AoA(3))+this%WWW1(2)*cos(this%AoA(3)),   &
                            this%WWW1(1)*sin(this%AoA(2))*cos(this%AoA(3))-this%WWW1(2)*sin(this%AoA(3))    ]
            this%WWW3(1:3)=matmul(this%TTT0(1:3,1:3),this%WWW2(1:3))

            call this%UpdateVelFromPosAngular(this%WWW3, this%UVW)
        endif

        call this%InitTriadANDFormMass()
        call this%InitGrav()

    end subroutine Beam_Initialise

    subroutine Beam_InitLoad(this)
        implicit none
        class(BeamSolver), intent(inout) :: this
        allocate(this%lodInte(1:this%gEQ),this%lodExte(1:this%gEQ))
        allocate(this%lodEffe(1:this%gEQ),this%lodFlow(1:this%gEQ))
        allocate(this%lodRepl(1:this%gEQ),this%lodGrav(1:this%gEQ))
        this%lodInte(:) = 0.0e0_rp
        this%lodExte(:) = 0.0e0_rp
        this%lodEffe(:) = 0.0e0_rp
        this%lodFlow(:) = 0.0e0_rp
        this%lodRepl(:) = 0.0e0_rp
        this%lodGrav(:) = 0.0e0_rp
    end subroutine Beam_InitLoad

    subroutine Beam_InitPosDspVelAcc(this)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer :: n
        allocate(this%pos(1:6, 1:this%nND),this%dsp(1:6, 1:this%nND),this%vel(1:6, 1:this%nND),this%acc(1:6, 1:this%nND),this%mss(1:3, 1:this%nND))
        this%pos(:, :) = 0.0e0_rp
        do n = 1,this%nEL
            this%pos(1:6,this%m_elements(n)%node0)=this%m_elements(n)%x1(1:6)
            this%pos(1:6,this%m_elements(n)%node1)=this%m_elements(n)%x1(7:12)
        enddo
        this%dsp(:, :) = 0.0e0_rp
        this%vel(:, :) = 0.0e0_rp
        this%acc(:, :) = 0.0e0_rp
        this%mss(:, :) = 0.0e0_rp
    end subroutine Beam_InitPosDspVelAcc

    subroutine Beam_UpdateVelFromPosAngular(this, WWW, UVW)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(in) :: WWW(3), UVW(3)
        real(rp) :: rel(3)
        integer :: n
        do n = 1, this%nND
            ! Position relative to the reference point.
            rel(1:3) = this%pos(1:3,n) - this%XYZ(1:3)
            ! v = U + omega x r
            this%vel(1:3, n) = [ WWW(2)*rel(3) - WWW(3)*rel(2), &
                                 WWW(3)*rel(1) - WWW(1)*rel(3), &
                                 WWW(1)*rel(2) - WWW(2)*rel(1) ] + UVW(1:3)
            this%vel(4:6, n) = WWW(1:3)
        end do
    
    end subroutine Beam_UpdateVelFromPosAngular

    subroutine Beam_InitTriadANDFormMass(this)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer:: i

        do i = 1, this%nEL
            ! InitTriad
            call this%m_elements(i)%InitTriad_D
            ! FormMass
            call this%m_elements(i)%FormMassMatrix
            call this%m_elements(i)%RotateMatrix
            call this%m_elements(i)%RKR(this%m_elements(i)%m_masMat)
        enddo

        return
    end subroutine Beam_InitTriadANDFormMass

    subroutine Beam_InitGrav(this)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer :: i, j, node0, node1
        real(rp) :: gvec(1:12), grav(1:12)
        gvec(:) = 0.0e0_rp
        gvec(1:3) = m_g(1:3)
        gvec(7:9) = m_g(1:3)
        do i = 1, this%nEL
            call this%m_elements(i)%MassMultiply(gvec, grav)
            call this%m_elements(i)%LocToGlobal(grav, this%lodGrav, this%gEQ)
            node0 = this%m_elements(i)%node0
            node1 = this%m_elements(i)%node1
            do j = 1, 3
                this%mss(j,node0) = this%mss(j,node0) + this%m_elements(i)%m_masMat(j,j)
                this%mss(j,node1) = this%mss(j,node1) + this%m_elements(i)%m_masMat(j+6,j+6)
            enddo
        enddo
    end subroutine Beam_InitGrav

    subroutine Beam_calculate_angle_material(this, Lref, Uref, denIn, uMax, uuuIn, nLthck)
    implicit none
    class(BeamSolver), intent(inout) :: this
    real(rp),intent(in) :: Lref, Uref, denIn, uuuIn(3)
    real(rp),intent(inout) :: uMax
    real(rp),intent(out) :: nLthck
    real(rp) :: xmax, ymax, zmax
    real(rp) :: rRot(3)
    integer :: i
    real(rp) :: len,ratio

    this%St = Lref * this%Freq / Uref
    ! angle to radian
    this%AoAo(1:3)=this%AoAo(1:3)/180.0e0_rp*m_pi
    this%AoAAmpl(1:3)=this%AoAAmpl(1:3)/180.0e0_rp*m_pi
    this%AoAPhi(1:3)=this%AoAPhi(1:3)/180.0e0_rp*m_pi
    this%XYZPhi(1:3)=this%XYZPhi(1:3)/180.0e0_rp*m_pi
    xmax = maxval([(abs(this%m_elements(i)%x00(1)), abs(this%m_elements(i)%x00(7)), i=1,this%nEL)])
    ymax = maxval([(abs(this%m_elements(i)%x00(2)), abs(this%m_elements(i)%x00(8)), i=1,this%nEL)])
    zmax = maxval([(abs(this%m_elements(i)%x00(3)), abs(this%m_elements(i)%x00(9)), i=1,this%nEL)])
    ! effective rotation radius for rotations about x, y, z
    rRot(1) = max(ymax, zmax)   ! rotation about x
    rRot(2) = max(xmax, zmax)   ! rotation about y
    rRot(3) = max(xmax, ymax)   ! rotation about z
    uMax = maxval([ uMax, &
                maxval(abs(uuuIn(1:3))), &
                2.0e0_rp*m_pi*maxval(abs(this%XYZAmpl(1:3)))*this%Freq, &
                2.0e0_rp*m_pi*maxval(abs(this%AoAAmpl(1:3))*rRot(1:3))*this%Freq ])

    ! Standard section-property convention used by the beam element:
    ! len = b     : local y-axis spanwise width
    ! nLthck = h  : local z-axis thickness, with b > h
    ! 1 E         : Young's modulus
    ! 2 G         : shear modulus
    ! 3 A         : cross-sectional area
    ! 4 rho       : density
    ! 5 gamma     : self-rotation angle in degrees, unused
    ! 6 Jt        : St. Venant torsion constant about local x-axis
    ! 7 Iy        : second moment of area about local y-axis
    ! 8 Iz        : second moment of area about local z-axis
    !
    ! For the automatically generated rectangular plate section:
    ! A = b*h, Iy = b*h^3/12, Iz = h*b^3/12.

    ! property data will use the parameters read from the file if isKB != 0 or 1.
    ! calculate material parameters
    if(m_isKB==0)then
        do i = 1, this%nEL
            len = this%m_elements(i)%spanlen
            this%m_elements(i)%m_property(1) = this%EmR * denIn * Uref**2
            this%m_elements(i)%m_property(2) = this%m_elements(i)%m_property(1) / (2.0e0_rp*(1.0e0_rp+this%psR))
            nLthck = this%tcR*Lref
            this%m_elements(i)%m_property(3) = len * nLthck
            this%m_elements(i)%m_property(4) = this%denR*len*Lref*denIn / this%m_elements(i)%m_property(3)
            ratio = nLthck/len
            this%m_elements(i)%m_property(6) = len * nLthck**3/3.0e0_rp * (1e0_rp-0.63e0_rp*ratio+0.052e0_rp*(ratio)**5)
            this%m_elements(i)%m_property(7) = len * nLthck**3/12.0e0_rp
            this%m_elements(i)%m_property(8) = nLthck * len**3/12.0e0_rp
        enddo
        len = this%m_elements(1)%spanlen
        this%KB = this%m_elements(1)%m_property(1)*this%m_elements(1)%m_property(7) / (denIn*Uref**2*Lref**3*len)
        this%KS = this%m_elements(1)%m_property(1)*this%m_elements(1)%m_property(3) / (denIn*Uref**2*Lref*len)
    endif

    if(m_isKB==1)then
        do i = 1, this%nEL
            len = this%m_elements(i)%spanlen
            nLthck =  sqrt(this%KB/this%KS*12.0e0_rp) * Lref
            this%m_elements(i)%m_property(3) = len * nLthck
            this%m_elements(i)%m_property(1) = this%KS * denIn * Uref**2 * Lref * len / this%m_elements(i)%m_property(3)
            this%m_elements(i)%m_property(2) = this%m_elements(i)%m_property(1) / (2.0e0_rp*(1.0e0_rp+this%psR))
            this%m_elements(i)%m_property(4) = this%denR*len*Lref*denIn / this%m_elements(i)%m_property(3)
            ratio = nLthck / len
            this%m_elements(i)%m_property(6) = len * nLthck**3/3.0e0_rp * (1.0e0_rp - 0.63e0_rp*ratio + 0.052e0_rp*ratio**5)
            this%m_elements(i)%m_property(7) = len * nLthck**3/12.0e0_rp
            this%m_elements(i)%m_property(8) = nLthck * len**3/12.0e0_rp
        enddo
        len = this%m_elements(1)%spanlen
        nLthck = this%m_elements(1)%m_property(3) / len
        this%EmR = this%m_elements(1)%m_property(1)/(denIn*Uref**2)
        this%tcR = nLthck/Lref
    endif

    end subroutine Beam_calculate_angle_material

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Beam_write_solid(this,Lref,Uref,Aref,Fref,iFish,idfile)
    implicit none
    class(BeamSolver), intent(inout) :: this
    real(rp):: Lref,Uref,Aref,Fref
    integer:: i,iFish,ElmType,idfile
    integer,parameter::nameLen=10,numVar=15
    character(len=nameLen):: idstr

    !Write zone information
    write(idstr,  '(I4.4)') iFish ! assume iFish < 10000
    write(idfile, '(A,A,A)') ' ZONE T = "fish',trim(idstr),'"'
    write(idfile, '(A)') ' STRANDID=0, SOLUTIONTIME=0'
    write(idfile, '(A,I8,A,I8,A)', advance='no') ' Nodes=',this%nND,', Elements=',this%nEL,', ZONETYPE='
    ElmType = this%m_elements(1)%itype
    if(ElmType.eq.2) then
        write(idfile, '(A)') 'FELINESEG'
    elseif (ElmType.eq.3) then
        ! write(idfile, '(A)') 'FETRIANGLE'
    endif
    write(idfile, '(A)') ' DATAPACKING=POINT'
    write(idfile, '(A)', advance='no') ' DT=('
    do i=1,numVar-1
        write(idfile, '(A)', advance='no') 'SINGLE '
    enddo
    write(idfile, '(A)') 'SINGLE )'
    !Write node data
    do i = 1, this%nND
        write(idfile,'(15E28.18)') this%pos(1:3,i)/Lref, &
                                   this%vel(1:3,i)/Uref, &
                                   this%acc(1:3,i)/Aref, &
                                   this%lodFlow((i-1)*6+1:(i-1)*6+3)/Fref, &
                                   this%lodRepl((i-1)*6+1:(i-1)*6+3)/Fref
    enddo
    !Write element data
    if(ElmType.eq.2) then
        do i = 1, this%nEL
            write(idfile, *) this%m_elements(i)%node0,this%m_elements(i)%node1
        enddo
    elseif (ElmType.eq.3) then
    !     do i = 1, this%nEL
    !         write(idfile, *) this%ele(i,1),this%ele(i,2),this%ele(i,3)
    !     enddo
    endif
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write check point file for restarting simulation
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Beam_write_solid_temp(this,fid)
    implicit none
    class(BeamSolver), intent(inout) :: this
    integer,intent(in) :: fid
        write(fid) this%m_elements
        write(fid) this%nND,this%nEL,this%nMT,this%gEQ
        write(fid) this%pos,this%dsp,this%vel,this%acc,this%mss
        write(fid) this%lodInte,this%lodExte,this%lodEffe
        write(fid) this%lodFlow,this%lodRepl,this%lodGrav,this%vBC
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read check point file for restarting simulation
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Beam_read_solid_temp(this,fid)
    implicit none
    class(BeamSolver), intent(inout) :: this
    integer,intent(in) :: fid
        read(fid) this%m_elements
        read(fid) this%nND,this%nEL,this%nMT,this%gEQ
        read(fid) this%pos,this%dsp,this%vel,this%acc,this%mss
        read(fid) this%lodInte,this%lodExte,this%lodEffe
        read(fid) this%lodFlow,this%lodRepl,this%lodGrav,this%vBC
    end subroutine

    subroutine Beam_write_solid_params(this,fid)
    implicit none
    class(BeamSolver), intent(inout) :: this
    integer,intent(in) :: fid
        !write(fid,'(A,2F20.10)')'Freq,St     =',this%Freq,this%St
        !write(fid,'(A,2F20.10)')'denR,psR    =',this%denR,this%psR
        !write(fid,'(A,2F20.10)')'KB,  KS     =',this%KB,this%KS
        write(fid,'(A,2F20.10)')'EmR, tcR    =',this%EmR,this%tcR
        !write(fid,'(A,1x,3F20.10,2x)')'XYZo(1:3)   =',this%XYZo(1:3)
        !write(fid,'(A,1x,3F20.10,2x)')'XYZAmpl(1:3)=',this%XYZAmpl(1:3)
        !write(fid,'(A,1x,3F20.10,2x)')'XYZPhi(1:3) =',this%XYZPhi(1:3)
        !write(fid,'(A,1x,3F20.10,2x)')'AoAo(1:3)   =',this%AoAo(1:3)
        !write(fid,'(A,1x,3F20.10,2x)')'AoAAmpl(1:3)=',this%AoAAmpl(1:3)
        !write(fid,'(A,1x,3F20.10,2x)')'AoAPhi(1:3) =',this%AoAPhi(1:3)
        write(fid,'(2(A,1x,I8,2x))')'nND =',this%nND,'nEL =',this%nEL
        write(fid,'(2(A,1x,I8,2x))')'nMT =',this%nMT,'gEQ =',this%gEQ
    end subroutine

    subroutine Beam_write_solid_materials(this,fid)
    implicit none
    class(BeamSolver), intent(inout) :: this
    character(len=4):: IDstr
    integer,intent(in) :: fid
        write(IDstr,'(I4.4)') 1
        write(fid,'(A,A,A  )')'------------------------------- m_element( ',IDstr,' ) ------------------------------'
        write(fid,'(A,E20.10 )')'E    =',this%m_elements(1)%m_property(1)
        write(fid,'(A,E20.10 )')'G    =',this%m_elements(1)%m_property(2)
        write(fid,'(A,E20.10 )')'A    =',this%m_elements(1)%m_property(3)
        write(fid,'(A,E20.10 )')'rho  =',this%m_elements(1)%m_property(4)
        write(fid,'(A,E20.10 )')'gamma=',this%m_elements(1)%m_property(5)
        write(fid,'(A,E20.10 )')'Jt   =',this%m_elements(1)%m_property(6)
        write(fid,'(A,E20.10 )')'Iy   =',this%m_elements(1)%m_property(7)
        write(fid,'(A,E20.10 )')'Iz   =',this%m_elements(1)%m_property(8)
    end subroutine

    subroutine Beam_write_solid_info(this,groupNum,XYZo,Lref,Uref,Aref,Fref,Pref,Eref)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp),intent(in):: XYZo(1:3),Lref,Uref,Aref,Fref,Pref,Eref
        character(LEN=3),intent(in):: groupNum
        integer:: idfile=100,i,node0,node1
        real(rp):: msum(3),xcm(3),vcm(3),acm(3)
        real(rp):: Ptot,Pax,Pay,Paz
        real(rp):: Etot,Ev,Ep,Es,Eb
        real(rp)::ve(12)
        
        ! write begin information
        open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_firstNode.dat',position='append')
        write(idfile,'(15E20.10)')XYZo(1:3)/Lref,(this%pos(1:3,1)-XYZo(1:3))/Lref,this%pos(4:6,1),this%vel(1:3,1)/Uref,this%acc(1:3,1)/Aref
        close(idfile)
        ! write end information titles
        open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_lastNode.dat',position='append')
        write(idfile,'(15E20.10)')XYZo(1:3)/Lref,(this%pos(1:3,this%nND)-XYZo(1:3))/Lref,this%pos(4:6,this%nND),this%vel(1:3,this%nND)/Uref,this%acc(1:3,this%nND)/Aref
        close(idfile)
        ! write center information titles
        open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_centerNode.dat',position='append')
        write(idfile,'(15E20.10)')XYZo(1:3)/Lref,(this%pos(1:3,(this%nND+1)/2)-XYZo(1:3))/Lref,this%pos(4:6,(this%nND+1)/2),this%vel(1:3,(this%nND+1)/2)/Uref,this%acc(1:3,(this%nND+1)/2)/Aref
        close(idfile)
        ! write mean information titles
        open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_nodeAverage.dat',position='append')
        msum=sum(this%mss(1:3,1:this%nND),2)
        xcm=sum(this%pos(1:3,1:this%nND)*this%mss(1:3,1:this%nND),2)/msum
        vcm=sum(this%vel(1:3,1:this%nND)*this%mss(1:3,1:this%nND),2)/msum
        acm=sum(this%acc(1:3,1:this%nND)*this%mss(1:3,1:this%nND),2)/msum
        write(idfile,'(12E20.10)') XYZo(1:3)/Lref,(xcm-XYZo(1:3))/Lref,vcm/Uref,acm/Aref
        close(idfile)
        ! write forces
        open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_forces.dat',position='append')
        write(idfile,'(6E20.10)')XYZo(1:3)/Lref,sum(this%lodFlow(1:this%gEQ:6))/Fref, &
                                                sum(this%lodFlow(2:this%gEQ:6))/Fref, &
                                                sum(this%lodFlow(3:this%gEQ:6))/Fref
        close(idfile)
        ! write power
        Pax=sum(this%lodFlow(1:this%gEQ:6)*this%vel(1,1:this%nND))/Pref
        Pay=sum(this%lodFlow(2:this%gEQ:6)*this%vel(2,1:this%nND))/Pref
        Paz=sum(this%lodFlow(3:this%gEQ:6)*this%vel(3,1:this%nND))/Pref
        Ptot=Pax+Pay+Paz
        open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_power.dat',position='append')
        write(idfile,'(7E20.10)')XYZo(1:3)/Lref,Ptot,Pax,Pay,Paz
        close(idfile)

        call this%UpdateStrainEnergy()
        Es=0.0e0_rp
        Eb=0.0e0_rp
        do i=1,this%nEL
            Es=Es+this%m_elements(i)%strainEnergy(1)
            Eb=Eb+this%m_elements(i)%strainEnergy(2)
        enddo
        Es=Es/Eref
        Eb=Eb/Eref
        Ep=Es+Eb
        Ev=0.0e0_rp
        do i=1, this%nEL
            node0=this%m_elements(i)%node0
            node1=this%m_elements(i)%node1
            ve(1:6) =this%vel(1:6,node0)
            ve(7:12)=this%vel(1:6,node1)
            Ev=Ev+0.5e0_rp*dot_product(ve, matmul(this%m_elements(i)%m_masMat, ve))
        enddo
        Ev=Ev/Eref
        Etot=Ev+Ep
        ! write energy title
        open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_energy.dat',position='append')
        write(idfile,'(8E20.10)')XYZo(1:3)/Lref,Etot,Ev,Ep,Es,Eb
        close(idfile)
    end subroutine

    subroutine Beam_write_solid_probes(this,groupNum,XYZo,solidProbingNum,solidProbingNode,Lref,Uref,Aref)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer,intent(in):: solidProbingNum,solidProbingNode(solidProbingNum)
        real(rp),intent(in):: XYZo(1:3),Lref,Uref,Aref
        character(LEN=3),intent(in):: groupNum
        character (LEN=4):: probeNum
        integer:: i,idfile=100
        do  i=1,solidProbingNum
            write(probeNum,'(I4.4)') i
            open(idfile,file='./DatInfo/Group'//trim(groupNum)//'_solidProbes_'//trim(probeNum)//'.dat',position='append')
            write(idfile,'(15E20.10)')XYZo(1:3)/Lref,(this%pos(1:3,solidProbingNode(i))-XYZo(1:3))/Lref, &
                                                                              this%pos(4:6,solidProbingNode(i)), &
                                                                              this%vel(1:3,solidProbingNode(i))/Uref, &
                                                                              this%acc(1:3,solidProbingNode(i))/Aref
            close(idfile)
        enddo
    end subroutine

    subroutine Beam_structure(this,iFish,time,isubstep,deltat,subdeltat)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer:: iFish,isubstep,i
        real(rp):: deltat,subdeltat,time
        !!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iND)
            if(this%iBodyModel.eq.1)then     ! rigid body
                !======================================================
                !prescribed motion
                !------------------------------------------------------
                !translational displacement
                this%XYZ(1:3)=this%XYZo(1:3)+this%XYZAmpl(1:3)*cos(2.0e0_rp*m_pi*this%Freq*(time-deltat+real(isubstep, rp)*subdeltat)+this%XYZPhi(1:3)) + this%initXYZVel(1:3) * (time-deltat+real(isubstep, rp)*subdeltat)
                !rotational displacement
                this%AoA(1:3)=this%AoAo(1:3)+this%AoAAmpl(1:3)*cos(2.0e0_rp*m_pi*this%Freq*(time-deltat+real(isubstep, rp)*subdeltat)+this%AoAPhi(1:3))
                call AoAtoTTT(this%AoA(1:3),this%TTTnxt(1:3,1:3))
                call Segment_get_angle_triad(this%TTT0(1:3,1:3),this%TTTnxt(1:3,1:3),this%AoAd(1),this%AoAd(2),this%AoAd(3))
                !given displacement
                do i=1,this%nEL
                    call this%m_elements(i)%MapReferencePosToCurrent(this%m_elements(i)%xnxt, this%TTTnxt, this%XYZ, this%AoAd)
                    call this%m_elements(i)%MapReferenceDirToCurrent(this%m_elements(i)%dircnxt, this%TTTnxt)
                    this%m_elements(i)%x1(:)=this%m_elements(i)%xnxt(:)
                    this%m_elements(i)%dirc1(:)=this%m_elements(i)%dircnxt(:)
                    this%pos(1:6,this%m_elements(i)%node0)=this%m_elements(i)%x1(1:6)
                    this%pos(1:6,this%m_elements(i)%node1)=this%m_elements(i)%x1(7:12)
                    call this%m_elements(i)%cptdxyz1()
                    call this%m_elements(i)%RigidUpdateTriad_D()
                enddo
                !------------------------------------------------------
                !translational velocity
                this%UVW(1:3) =-2.0e0_rp*m_pi*this%Freq*this%XYZAmpl(1:3)*sin(2.0e0_rp*m_pi*this%Freq*(time-deltat+real(isubstep, rp)*subdeltat)+this%XYZPhi(1:3)) + this%initXYZVel(1:3)
                !rotational velocity
                this%WWW1(1:3)=-2.0e0_rp*m_pi*this%Freq*this%AoAAmpl(1:3)*sin(2.0e0_rp*m_pi*this%Freq*(time-deltat+real(isubstep, rp)*subdeltat)+this%AoAPhi(1:3))
                this%WWW2(1:3)=[this%WWW1(1)*cos(this%AoA(2))+this%WWW1(3),    &
                                this%WWW1(1)*sin(this%AoA(2))*sin(this%AoA(3))+this%WWW1(2)*cos(this%AoA(3)), &
                                this%WWW1(1)*sin(this%AoA(2))*cos(this%AoA(3))-this%WWW1(2)*sin(this%AoA(3))]
                this%WWW3(1:3)=matmul(this%TTTnxt(1:3,1:3),this%WWW2(1:3))
                !given velocity
                call this%UpdateVelFromPosAngular(this%WWW3, this%UVW)
                !-------------------------------------------------------
            elseif(this%iBodyModel.eq.2)then !elastic model
                !translational displacement
                this%XYZ(1:3)=this%XYZo(1:3)+this%XYZAmpl(1:3)*cos(2.0e0_rp*m_pi*this%Freq*(time-deltat+real(isubstep, rp)*subdeltat)+this%XYZPhi(1:3)) + this%initXYZVel(1:3) * (time-deltat+real(isubstep, rp)*subdeltat)
                !rotational displacement
                this%AoA(1:3)=this%AoAo(1:3)+this%AoAAmpl(1:3)*cos(2.0e0_rp*m_pi*this%Freq*(time-deltat+real(isubstep, rp)*subdeltat)+this%AoAPhi(1:3))
                call AoAtoTTT(this%AoA(1:3),this%TTTnxt(1:3,1:3))
                call Segment_get_angle_triad(this%TTT0(1:3,1:3),this%TTTnxt(1:3,1:3),this%AoAd(1),this%AoAd(2),this%AoAd(3))
                !given displacement
                do i=1,this%nEL
                    call this%m_elements(i)%MapReferencePosToCurrent(this%m_elements(i)%xnxt, this%TTTnxt, this%XYZ, this%AoAd)
                enddo
                !-----------------------------------------
                call this%UpdateNewmarkCoeffs(subdeltat)
                CALL this%Solver(iFish)
            else
                stop 'no define body model'
            endif
            !!$OMP END PARALLEL DO
    end subroutine

    subroutine Beam_UpdateNewmarkCoeffs(this,dt)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(in) :: dt
        ! coefficients in the Newmark-beta method
        this%coeffs(2) = 1.0e0_rp / (m_NewmarkBeta * dt)
        this%coeffs(1) = m_NewmarkGamma * this%coeffs(2)
        this%coeffs(0) = this%coeffs(2) / dt
        this%coeffs(3) = 0.5e0_rp / m_NewmarkBeta - 1.0e0_rp
        this%coeffs(4) = m_NewmarkGamma / m_NewmarkBeta - 1.0e0_rp
        this%coeffs(5) = dt * (0.5e0_rp * m_NewmarkGamma / m_NewmarkBeta - 1.0e0_rp)
        this%coeffs(6) = dt * (1.0e0_rp - m_NewmarkGamma)
        this%coeffs(7) = dt * m_NewmarkGamma
    end subroutine Beam_UpdateNewmarkCoeffs

    subroutine Beam_Solver(this,iFish)
        ! Nonlinear Algorithm: full Newton-Raphson method
        ! ISBN 9781441929105 James F. Doyle. P354
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp) :: dspO(1:6, 1:this%nND), velO(1:6, 1:this%nND), accO(1:6, 1:this%nND)
        real(rp) :: dspn(1:this%gEQ)
        integer :: iter, iFish
        real(rp) :: dnorm
        dnorm = 1.0e0_rp
        ! solve the next dispalce, velocity and acceleration using CG method
        call this%InitDspVelAccATTimeT(dspO, velO, accO)
        do iter = 1, m_ntolFEM ! m_ntolFEM is maxNewtonRaphson
            call this%UpdateMatrixANDLoad(dspO)
            call this%CG_Solve(dspn, this%lodEffe, iter)
            call this%UpdateDspANDTride(iter, dspn, dnorm)
            if (dnorm .le. m_dtolFEM) exit
        enddo
        call this%UpdateVelAcc(dspO, velO, accO)
        call this%UpdateIterInfo(iFish, iter, dnorm)
    end subroutine Beam_Solver

    subroutine Beam_InitDspVelAccATTimeT(this, dspO, velO, accO)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(out) :: dspO(1:6, 1:this%nND), velO(1:6, 1:this%nND), accO(1:6, 1:this%nND)
        integer :: i
        this%vBC(:) = 0.0e0_rp
        do i=1,this%nEL
            ! displacement condition
            this%vBC(this%m_elements(i)%m_localToGlobal(1:12))=this%m_elements(i)%xnxt(1:12)-this%m_elements(i)%x1(1:12)
        enddo
        this%lodExte = this%lodFlow + this%lodGrav + this%lodRepl
        dspO(1:6, 1:this%nND) = this%dsp(1:6, 1:this%nND)
        velO(1:6, 1:this%nND) = this%vel(1:6, 1:this%nND)
        accO(1:6, 1:this%nND) = this%acc(1:6, 1:this%nND)
    end subroutine Beam_InitDspVelAccATTimeT

    subroutine Beam_UpdateMatrixANDLoad(this, dspO)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(inout) :: dspO(1:6, 1:this%nND)
        integer :: i
        this%lodInte = 0.0e0_rp
        do i = 1, this%nEL
            call this%m_elements(i)%FormStiffMatrix()
            call this%m_elements(i)%BodyStress_D(this%lodInte, this%gEQ)
            call this%m_elements(i)%UpdateMatrix(this%coeffs, m_GeoGamma, m_dampM, m_dampK)
        enddo
        this%lodEffe = this%lodExte - this%lodInte
        do i = 1, this%nEL
            call this%m_elements(i)%UpdateLoad(this%coeffs, m_dampM, m_dampK, this%nND, this%gEQ, dspO, this%dsp, this%vel, this%acc, this%lodEffe)
        enddo
    end subroutine Beam_UpdateMatrixANDLoad

    subroutine Beam_CG_Solve(this, x, b, iterNR)
        ! Matrix-free Jacobi-preconditioned CG solver.
        ! Iterative solution for displacement vector.
        ! https://www.detailedpedia.com/wiki-Conjugate_gradient_method
        !
        ! Non-homogeneous Dirichlet boundary conditions are handled by lifting:
        !       x = xFixed + xFree,
        ! where xFixed satisfies the prescribed displacement increments and xFree is zero on fixed DOFs. Therefore,
        !       A*xFree = b - A*xFixed.
        ! The CG iteration is projected onto the free-DOF subspace by enforcing
        ! r(fixed)=0, z(fixed)=0, p(fixed)=0, and Ap(fixed)=0.
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer, intent(in) :: iterNR
        real(rp), intent(out) :: x(1:this%gEQ)
        real(rp), intent(in)  :: b(1:this%gEQ)
        real(rp) :: r(1:this%gEQ), p(1:this%gEQ), Ap(1:this%gEQ)
        real(rp) :: z(1:this%gEQ), M(1:this%gEQ)
        real(rp) :: blockM(6,6,this%nND)
        real(rp) :: xFixed(1:this%gEQ)
        logical :: fixed(1:this%gEQ)
        integer :: iter, max_iter, precondType
        real(rp) :: alpha, beta, rsold, rsnew, err, pAp, residual_norm
        ! select preconditioner
        precondType = 2
        ! maximum iteration count and tolerance
        max_iter = 10000
        err = 1.0e-6_rp
        ! lifting of non-homogeneous Dirichlet boundary conditions
        call this%BoundaryCond(iterNR, xFixed, fixed)
        ! initial guess
        x(1:this%gEQ) = xFixed(1:this%gEQ)
        ! initial residual: r = b - A*x
        call this%MatrixMultipy(x, Ap)
        r = b - Ap
        call this%ApplyFixedZero(r, fixed)
        call this%InitPreconditioner(precondType, M, blockM, fixed)
        call this%ApplyPreconditioner(precondType, r, z, M, blockM, fixed)
        residual_norm = sqrt(dot_product(r, r))
        if (residual_norm .le. err) return
        p = z
        rsold = dot_product(r, z)
        call this%CheckCGBreakdown(rsold, residual_norm, 1.0e-30_rp, 'rsold', 'initial dot_product(r,z)')
        do iter = 1, max_iter
            call this%MatrixMultipy(p, Ap)
            call this%ApplyFixedZero(Ap, fixed)
            pAp = dot_product(p, Ap)
            call this%CheckCGBreakdown(pAp, residual_norm, 1.0e-30_rp, 'pAp', 'dot_product(p,Ap) before alpha', iter)
            alpha = rsold / pAp
            x = x + alpha * p
            ! In exact arithmetic, p(fixed)=0 keeps x(fixed)=xFixed(fixed).
            ! This projection is kept as a safeguard against numerical drift or future code changes.
            call this%ApplyFixedValue(x, xFixed, fixed)
            r = r - alpha * Ap
            call this%ApplyFixedZero(r, fixed)
            residual_norm = sqrt(dot_product(r, r))
            if (residual_norm .le. err) exit
            call this%ApplyPreconditioner(precondType, r, z, M, blockM, fixed)
            rsnew = dot_product(r, z)
            call this%CheckCGBreakdown(rsold, residual_norm, 1.0e-30_rp, 'rsold', 'rsold before beta=rsnew/rsold', iter)
            beta = rsnew / rsold
            p = z + beta * p
            call this%ApplyFixedZero(p, fixed)
            rsold = rsnew
        enddo
        if (residual_norm .gt. err) then
            write(*,*) 'WARNING: CG did not fully converge.'
            write(*,*) 'residual_norm = ', residual_norm, ' err = ', err
            write(*,*) 'max_iter = ', max_iter
        endif
        return
    end subroutine Beam_CG_Solve

    subroutine Beam_BoundaryCond(this, iter, x, fixed)
        class(BeamSolver), intent(in) :: this
        integer, intent(in) :: iter
        real(rp), intent(out) :: x(1:this%gEQ)
        logical, intent(out) :: fixed(1:this%gEQ)
        integer :: i
        fixed(:)  = .false.
        x(:) = 0.0e0_rp
        do i=1,this%nEL
            call this%m_elements(i)%BoundaryCond(iter, this%gEQ, x, fixed, this%vBC)
        enddo
    end subroutine Beam_BoundaryCond

    subroutine Beam_MatrixMultipy(this, x, b)
        ! Matrix Free Method
        ! Implement matrix-vector multiplication for individual elements, without assembling the global matrix.
        ! Element‐by‐Element(Hughes, Thomas J. R.(1983).doi:10.1061/(ASCE)0733-9399(1983)109:2(576))
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(in)  :: x(1:this%gEQ)
        real(rp), intent(out) :: b(1:this%gEQ)
        integer :: i
        b(1:this%gEQ) = 0.0e0_rp
        do i=1,this%nEL
            call this%m_elements(i)%Multiply(x, b, this%gEQ)
        enddo
        return
    end subroutine Beam_MatrixMultipy

    subroutine Beam_ApplyFixedZero(this, x, fixed)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(inout) :: x(1:this%gEQ)
        logical, intent(in) :: fixed(1:this%gEQ)
        integer :: i
        do i = 1, this%gEQ
            if (fixed(i)) x(i) = 0.0e0_rp
        enddo
    end subroutine Beam_ApplyFixedZero

    subroutine Beam_ApplyFixedValue(this, x, xFixed, fixed)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(inout) :: x(1:this%gEQ)
        real(rp), intent(in) :: xFixed(1:this%gEQ)
        logical, intent(in) :: fixed(1:this%gEQ)
        integer :: i
        do i = 1, this%gEQ
            if (fixed(i)) x(i) = xFixed(i)
        enddo
    end subroutine Beam_ApplyFixedValue

    subroutine Beam_InitPreconditioner(this, precondType, M, blockM, fixed)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer, intent(in) :: precondType
        real(rp), intent(inout) :: M(1:this%gEQ)
        real(rp), intent(inout) :: blockM(6,6,this%nND)
        logical, intent(in) :: fixed(1:this%gEQ)
        select case (precondType)
        case (1)
            ! Scalar Jacobi preconditioner.
            call this%ScalarPreconditioned(M)
            call this%InitScalarJacobiPreconditioner(M, fixed)
        case (2)
            ! Node-wise 6x6 block Jacobi preconditioner.
            call this%BlockPreconditioned(blockM)
            call this%InitBlockJacobiPreconditioner(blockM, fixed)
        case default
            write(*,*) 'ERROR: unknown preconditioner type: ', precondType
            stop
        end select
    end subroutine Beam_InitPreconditioner

    subroutine Beam_ApplyPreconditioner(this, precondType, r, z, M, blockM, fixed)
        implicit none
        class(BeamSolver), intent(in) :: this
        integer, intent(in) :: precondType
        real(rp), intent(in) :: r(1:this%gEQ)
        real(rp), intent(out) :: z(1:this%gEQ)
        real(rp), intent(in) :: M(1:this%gEQ)
        real(rp), intent(in) :: blockM(6,6,this%nND)
        logical, intent(in) :: fixed(1:this%gEQ)
        select case (precondType)
        case (1)
            call this%ApplyScalarJacobiPreconditioner(r, z, M, fixed)
        case (2)
            call this%ApplyBlockJacobiPreconditioner(r, z, blockM, fixed)
        case default
            write(*,*) 'ERROR: unknown preconditioner type: ', precondType
            stop
        end select
    end subroutine Beam_ApplyPreconditioner

    subroutine Beam_ScalarPreconditioned(this, M)
        ! Scalar Jacobi Preconditioned
        class(BeamSolver), intent(in) :: this
        real(rp) :: M(1:this%gEQ)
        integer :: i
        M=0.0e0_rp
        do i=1,this%nEL
            call this%m_elements(i)%ScalarPreconditioned(M, this%gEQ)
        enddo
    end subroutine Beam_ScalarPreconditioned

    subroutine Beam_BlockPreconditioned(this, blockM)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(out) :: blockM(6,6,this%nND)
        integer :: i
        blockM(:,:,:) = 0.0e0_rp
        do i = 1, this%nEL
            call this%m_elements(i)%BlockPreconditioned(blockM, this%nND)
        enddo
    end subroutine Beam_BlockPreconditioned

    subroutine Beam_InitScalarJacobiPreconditioner(this, M, fixed)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(inout) :: M(1:this%gEQ)
        logical, intent(in) :: fixed(1:this%gEQ)
        integer :: i
        do i = 1, this%gEQ
            if (fixed(i)) then
                M(i) = 0.0e0_rp
            else
                if (abs(M(i)) .gt. 1.0e-30_rp) then
                    M(i) = 1.0e0_rp / M(i)
                else
                    write(*,*) 'WARNING: near-zero Jacobi diagonal at DOF ', i
                    M(i) = 1.0e0_rp
                endif
            endif
        enddo
    end subroutine Beam_InitScalarJacobiPreconditioner

    subroutine Beam_ApplyScalarJacobiPreconditioner(this, r, z, M, fixed)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(in) :: r(1:this%gEQ), M(1:this%gEQ)
        real(rp), intent(out) :: z(1:this%gEQ)
        logical, intent(in) :: fixed(1:this%gEQ)
        integer :: i
        do i = 1, this%gEQ
            if (fixed(i)) then
                z(i) = 0.0e0_rp
            else
                z(i) = M(i) * r(i)
            endif
        enddo
    end subroutine Beam_ApplyScalarJacobiPreconditioner

    subroutine Beam_InitBlockJacobiPreconditioner(this, blockM, fixed)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(inout) :: blockM(6,6,this%nND)
        logical, intent(in) :: fixed(1:this%gEQ)
        integer :: node, i, gid
        real(rp) :: invB(6,6)
        integer :: info
        do node = 1, this%nND
            ! For fixed DOFs, remove their coupling inside the local 6x6 block
            ! and set unit diagonal. Since r(fixed)=0, this keeps z(fixed)=0.
            do i = 1, 6
                gid = (node-1)*6 + i
                if (fixed(gid)) then
                    blockM(i,1:6,node) = 0.0e0_rp
                    blockM(1:6,i,node) = 0.0e0_rp
                    blockM(i,i,node)   = 1.0e0_rp
                endif
            enddo
            call Invert6x6(blockM(1:6,1:6,node), invB, info)
            if (info /= 0) then
                write(*,*) 'WARNING: block Jacobi inverse failed at node ', node
                write(*,*) 'Use identity block for this node.'
                blockM(1:6,1:6,node) = 0.0e0_rp
                do i = 1, 6
                    blockM(i,i,node) = 1.0e0_rp
                enddo
            else
                blockM(1:6,1:6,node) = invB(1:6,1:6)
            endif
        enddo
    end subroutine Beam_InitBlockJacobiPreconditioner

    subroutine Beam_ApplyBlockJacobiPreconditioner(this, r, z, blockM, fixed)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(in) :: r(1:this%gEQ)
        real(rp), intent(out) :: z(1:this%gEQ)
        real(rp), intent(in) :: blockM(6,6,this%nND)
        logical, intent(in) :: fixed(1:this%gEQ)
        integer :: node, i, gid
        real(rp) :: rnode(6), znode(6)
        z(:) = 0.0e0_rp
        do node = 1, this%nND
            do i = 1, 6
                gid = (node-1)*6 + i
                rnode(i) = r(gid)
            enddo
            znode = matmul(blockM(1:6,1:6,node), rnode)
            do i = 1, 6
                gid = (node-1)*6 + i
                if (fixed(gid)) then
                    z(gid) = 0.0e0_rp
                else
                    z(gid) = znode(i)
                endif
            enddo
        enddo
    end subroutine Beam_ApplyBlockJacobiPreconditioner

    subroutine Invert6x6(A, Ainv, info)
        implicit none
        real(rp), intent(in) :: A(6,6)
        real(rp), intent(out) :: Ainv(6,6)
        integer, intent(out) :: info
        real(rp) :: aug(6,12)
        real(rp) :: pivot, factor, tmp
        integer :: i, j, k, p
        info = 0
        aug(:,1:6) = A(:,:)
        aug(:,7:12) = 0.0e0_rp
        do i = 1, 6
            aug(i,6+i) = 1.0e0_rp
        enddo
        do i = 1, 6
            p = i
            pivot = abs(aug(i,i))
            do k = i+1, 6
                if (abs(aug(k,i)) .gt. pivot) then
                    pivot = abs(aug(k,i))
                    p = k
                endif
            enddo
            if (pivot .le. 1.0e-30_rp) then
                info = i
                Ainv(:,:) = 0.0e0_rp
                return
            endif
            if (p /= i) then
                do j = 1, 12
                    tmp = aug(i,j)
                    aug(i,j) = aug(p,j)
                    aug(p,j) = tmp
                enddo
            endif
            pivot = aug(i,i)
            aug(i,1:12) = aug(i,1:12) / pivot
            do k = 1, 6
                if (k /= i) then
                    factor = aug(k,i)
                    aug(k,1:12) = aug(k,1:12) - factor * aug(i,1:12)
                endif
            enddo
        enddo
        Ainv(:,:) = aug(:,7:12)
    end subroutine Invert6x6

    subroutine Beam_CheckCGBreakdown(this, value, residual_norm, threshold, valueName, locationName, iter)
        implicit none
        class(BeamSolver), intent(in) :: this
        real(rp), intent(in) :: value, residual_norm, threshold
        character(len=*), intent(in) :: valueName, locationName
        integer, intent(in), optional :: iter
    
        if (abs(value) .le. threshold) then
            write(*,*) 'ERROR: CG solver breakdown.'
            write(*,*) 'XYZo x y z    : ', this%XYZo(1:3)
            write(*,*) 'Location      : ', trim(locationName)
            write(*,*) 'Quantity      : ', trim(valueName)
            write(*,*) 'Value         : ', value
            write(*,*) 'Threshold     : ', threshold
            write(*,*) 'residual_norm : ', residual_norm
            if (present(iter)) then
                write(*,*) 'CG iter       : ', iter
            endif
            stop
        endif
    
    end subroutine Beam_CheckCGBreakdown

    subroutine Beam_UpdateDspANDTride(this, iter, dspn, dnorm)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(out) :: dnorm
        real(rp) :: dspn(1:this%gEQ), dspnn(1:6, 1:this%nND)
        real(rp) :: beta0,beta,zi,z0
        integer :: iter,i,maxramp
        beta0 = 1.0e0_rp
        maxramp = 1
        if    (iter <= maxramp) then
            zi=2.0e0_rp**(iter)
            z0=2.0e0_rp**(maxramp)
            beta=zi/z0*beta0
        else
            beta=1.0e0_rp*beta0
        endif
        do i = 1, this%nND
            dspnn(1:6,i)= beta*dspn((i-1)*6+1:(i-1)*6+6)
        enddo
        this%dsp(1:6, 1:this%nND) = this%dsp(1:6, 1:this%nND) + dspnn(1:6, 1:this%nND)
        do i = 1, this%nEL
            ! UpdateElementPosition
            this%m_elements(i)%x1(1:6)  = this%m_elements(i)%x0(1:6)  + this%dsp(1:6, this%m_elements(i)%node0)
            this%m_elements(i)%x1(7:12) = this%m_elements(i)%x0(7:12) + this%dsp(1:6, this%m_elements(i)%node1)
            this%pos(1:6,this%m_elements(i)%node0)=this%m_elements(i)%x1(1:6)
            this%pos(1:6,this%m_elements(i)%node1)=this%m_elements(i)%x1(7:12)
            call this%m_elements(i)%cptdxyz1
            ! UpdateNodeTriad
            call this%m_elements(i)%UpdateTriad_D(dspnn,this%nND)
            ! MakeElementTriad
            call this%m_elements(i)%MakeTriad_ee
        enddo
        dnorm=abs(maxval((beta*dspn(1:this%gEQ))**2))
        return
    end subroutine Beam_UpdateDspANDTride

    subroutine Beam_UpdateStrainEnergy(this)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer :: i
        do i = 1, this%nEL
            call this%m_elements(i)%FormStiffMatrix
            call this%m_elements(i)%StrainEnergy_D
        enddo
    end subroutine Beam_UpdateStrainEnergy

    subroutine Beam_UpdateVelAcc(this, dspO, velO, accO)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(rp), intent(inout) :: dspO(1:6, 1:this%nND), velO(1:6, 1:this%nND), accO(1:6, 1:this%nND)
        integer :: nND
        nND=this%nND
        this%acc(1:6, 1:nND)  = this%coeffs(0)*(this%dsp(1:6, 1:nND) - dspO(1:6, 1:nND)) -this%coeffs(2)*velO(1:6, 1:nND) - this%coeffs(3)*accO(1:6, 1:nND)
        this%vel(1:6, 1:nND)  = velO(1:6, 1:nND) + this%coeffs(6)*accO(1:6, 1:nND) + this%coeffs(7)*this%acc(1:6, 1:nND)
        return
    end subroutine Beam_UpdateVelAcc

    subroutine Beam_UpdateIterInfo(this, iFish, iter, dnorm)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer :: iFish, iter
        real(rp) :: dnorm
        this%FishInfo(1)=real(iFish, rp)
        this%FishInfo(2)=real(iter, rp)
        this%FishInfo(3)=dnorm
    end subroutine Beam_UpdateIterInfo

    subroutine Beam_ReportDispFieldStat(nND, field, fileName)
        implicit none
        integer, intent(in) :: nND
        real(rp), intent(in) :: field(1:6, 1:nND)
        character(LEN=*), intent(in) :: fileName
        integer, parameter :: statUnit = 112

        open(unit=statUnit, file=fileName, status='replace')
        call Beam_ReportDispGroup(nND, statUnit, field(1:3,1:nND), 'DISP_TRANS', (/ 'Ux', 'Uy', 'Uz' /))
        call Beam_ReportDispGroup(nND, statUnit, field(4:6,1:nND), 'DISP_ROT',   (/ 'Rx', 'Ry', 'Rz' /))
        close(statUnit)
    end subroutine Beam_ReportDispFieldStat

    subroutine Beam_ReportDispGroup(nND, fileUnit, fieldGroup, groupName, dofName)
        implicit none
        integer, intent(in) :: nND
        integer, intent(in) :: fileUnit
        real(rp), intent(in) :: fieldGroup(1:3, 1:nND)
        character(LEN=*), intent(in) :: groupName
        character(LEN=2), intent(in) :: dofName(3)
        real(rp) :: l2(3), linfty(3)
        integer :: i

        do i = 1, 3
            l2(i) = sqrt(sum(fieldGroup(i,1:nND) * fieldGroup(i,1:nND)) / real(nND, rp))
            linfty(i) = maxval(abs(fieldGroup(i,1:nND)))
            write(fileUnit,'(A,1X,A,1X,A,1X,ES24.16)') 'FIELDSTAT', trim(groupName), 'L2 ' // trim(dofName(i)), l2(i)
        enddo
        do i = 1, 3
            write(fileUnit,'(A,1X,A,1X,A,1X,ES24.16)') 'FIELDSTAT', trim(groupName), 'Linfinity ' // trim(dofName(i)), linfty(i)
        enddo
    end subroutine Beam_ReportDispGroup

    subroutine AoAtoTTT(AoA,TTT)
        implicit none
        real(rp), parameter:: pi=3.141562653589793e0_rp,eps=1.0e-5_rp
        real(rp)::AoA(3),TTT(3,3),rrx(3,3),rry(3,3),rrz(3,3),vcos,vsin
        TTT(:,:)=0.0e0_rp
        TTT(1,1)=1.0e0_rp
        TTT(2,2)=1.0e0_rp
        TTT(3,3)=1.0e0_rp
    
        vcos=cos(AoA(1))
        vsin=sin(AoA(1))
        if(abs(AoA(1))<eps)then
        vcos=1.0e0_rp
        vsin=0.0e0_rp
        endif
        if(abs(AoA(1)-0.5e0_rp*pi)<eps)then
        vcos=0.0e0_rp
        vsin=1.0e0_rp
        endif
        if(abs(AoA(1)+0.5e0_rp*pi)<eps)then
        vcos=0.0e0_rp
        vsin=-1.0e0_rp
        endif
    
        rrx(1:3,1:3)=reshape([  1.0e0_rp,0.0e0_rp,0.0e0_rp,  &
                                0.0e0_rp,vcos,vsin, &
                                0.0e0_rp,-vsin,vcos],[3,3])
        vcos=cos(AoA(2))
        vsin=sin(AoA(2))
        if(abs(AoA(2))<eps)then
        vcos=1.0e0_rp
        vsin=0.0e0_rp
        endif
        if(abs(AoA(2)-0.5e0_rp*pi)<eps)then
        vcos=0.0e0_rp
        vsin=1.0e0_rp
        endif
        if(abs(AoA(2)+0.5e0_rp*pi)<eps)then
        vcos=0.0e0_rp
        vsin=-1.0e0_rp
        endif
    
        rry(1:3,1:3)=reshape([  vcos,0.0e0_rp,-vsin,  &
                                0.0e0_rp,1.0e0_rp,0.0e0_rp, &
                               vsin,0.0e0_rp,vcos],[3,3])
    
        vcos=cos(AoA(3))
        vsin=sin(AoA(3))
        if(abs(AoA(3))<eps)then
        vcos=1.0e0_rp
        vsin=0.0e0_rp
        endif
        if(abs(AoA(3)-0.5e0_rp*pi)<eps)then
        vcos=0.0e0_rp
        vsin=1.0e0_rp
        endif
        if(abs(AoA(3)+0.5e0_rp*pi)<eps)then
        vcos=0.0e0_rp
        vsin=-1.0e0_rp
        endif
    
        rrz(1:3,1:3)=reshape([ vcos,vsin,0.0e0_rp, &
                              -vsin,vcos,0.0e0_rp, &
                               0.0e0_rp,0.0e0_rp,1.0e0_rp],[3,3])
        TTT=matmul(rrz,TTT)
        TTT=matmul(rry,TTT)
        TTT=matmul(rrx,TTT)
    
    end subroutine

end module SolidSolver
