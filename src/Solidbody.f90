module SolidBody
    use SolidSolver
    implicit none
    private
    ! Immersed boundary method parameters
    integer:: m_nFish, m_ntolLBM, m_zDim, m_yDim, m_xDim, m_Kb, m_Jb, m_Ib, m_Kmask, m_Jmask, m_Imask
    real(8):: m_dtolLBM, m_Pbeta, m_dt, m_dh, m_denIn, m_uuuIn(3), m_Aref, m_Eref, m_Fref, m_Lref, m_Pref, m_Tref, m_Uref
    integer:: m_boundaryConditions(1:6)
    ! nFish     number of bodies
    ! ntolLBM maximum number of iterations for IB force calculation
    ! dtolIBM   tolerance for IB force calculation
    ! Pbeta     coefficient in penalty force calculation
    public :: VBodies,read_solid_file,Initialise_bodies,allocate_solid_memory,Write_solid_v_bodies,FSInteraction_force, &
              Initialise_Calculate_Solid_params,Solver,Write_cont,Read_cont,write_solid_field,Write_solid_Check,Write_solid_Data,Write_SampBodyNode
    type :: VirtualBody
        type(BeamSolver):: rbm
        !!!virtual infomation
        !!!virtual body surface
        integer :: v_nelmts,v_type ! v_type 0 (solid body), 1 (plate), 2 (rod)
        real(8), allocatable :: v_Exyz(:, :) ! element center (x, y, z)
        real(8), allocatable :: v_Ea(:) ! element area
        !area center with equal weight on both sides
        real(8), allocatable :: v_Evel(:, :)
        real(8), allocatable :: v_Eforce(:, :) ! element center (x, y, z)
        integer(4), allocatable :: v_Ei(:, :) ! element stencial integer index [ix-1,ix,ix1,ix2, iy-1,iy,iy1,iy2, iz-1,iz,iz1,iz2]
        real(8), allocatable :: v_Ew(:, :) ! element stential weight [wx-1, wx, wx1, wx2, wy-1, wy, wy1, wy2, wz-1, wz, wz1, wz2]
        !calculated using central linear and angular velocities
        integer,allocatable :: vtor(:)!of size fake_npts
        integer,allocatable :: rtov(:)! of size real_npts+1
        integer,allocatable :: array_i(:)
        real(8),allocatable :: array_u(:,:)
    contains
        procedure :: Initialise => Initialise_
        procedure :: PlateBuild => PlateBuild_
        procedure :: UpdatePosVelArea => UpdatePosVelArea_
        procedure :: PlateUpdatePosVelArea => PlateUpdatePosVelArea_
        procedure :: Write_body => Write_body_
        procedure :: PlateWrite_body => PlateWrite_body_
        procedure :: UpdateElmtInterp => UpdateElmtInterp_
        procedure :: getuuu => getuuu_
        procedure :: setuuu => setuuu_
        procedure :: PenaltyForce => PenaltyForce_
        procedure :: FluidVolumeForce => FluidVolumeForce_
    end type VirtualBody
    type(VirtualBody), allocatable :: VBodies(:)
  contains
    subroutine read_solid_file(nFish,FEmeshName,iBodyModel,isMotionGiven,denR,KB,KS,EmR,psR,tcR,St, &
                               Freq,XYZo,XYZAmpl,XYZPhi,AoAo,AoAAmpl,AoAPhi, &
                               ntolLBM,dtolLBM,Pbeta,dt,denIn,uuuIn,boundaryConditions, &
                               dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,dtolFEM,ntolFEM,iForce2Body,iKB)
        integer,intent(in):: nFish
        character (LEN=40),intent(in):: FEmeshName(nFish)
        integer,intent(in):: iBodyModel(nFish),isMotionGiven(6,nFish)
        real(8),intent(in):: denR(nFish),KB(nFish),KS(nFish),EmR(nFish),psR(nFish),tcR(nFish),St(nFish)
        real(8),intent(in):: Freq(nFish)
        real(8),intent(in):: XYZo(3,nFish),XYZAmpl(3,nFish),XYZPhi(3,nFish)
        real(8),intent(in):: AoAo(3,nFish),AoAAmpl(3,nFish),AoAPhi(3,nFish)
        integer,intent(in):: ntolLBM,boundaryConditions(6)
        real(8),intent(in):: dtolLBM,Pbeta,dt,denIn,uuuIn(3)
        real(8),intent(in):: dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,dtolFEM
        integer,intent(in):: ntolFEM,iForce2Body,iKB
        integer:: iFish
        
        m_nFish = nFish

        m_ntolLBM = ntolLBM
        m_dtolLBM = dtolLBM
        m_Pbeta = Pbeta
        m_dt = dt
        m_denIn = denIn
        m_uuuIn = uuuIn
        m_boundaryConditions(1:6) = boundaryConditions(1:6)

        allocate(VBodies(m_nFish))

        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%Read_inFlow(FEmeshName(iFish),iBodyModel(iFish),isMotionGiven(1:6,iFish), &
                                                denR(iFish),KB(iFish),KS(iFish),EmR(iFish),psR(iFish),tcR(iFish),St(iFish), &
                                                Freq(iFish),XYZo(1:3,iFish),XYZAmpl(1:3,iFish),XYZPhi(1:3,iFish), &
                                                AoAo(1:3,iFish),AoAAmpl(1:3,iFish),AoAPhi(1:3,iFish))
        enddo

        CALL Read_SolidSolver_Params(dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,dtolFEM,ntolFEM,iForce2Body,iKB)
    end subroutine read_solid_file

    subroutine allocate_solid_memory(Asfac,Lchod,Lspan,AR)
        integer :: iFish,maxN
        real(8),intent(out):: Asfac,Lchod,Lspan,AR
        real(8) :: nAsfac(m_nFish),nLchod(m_nFish)
        write(*,'(A)') '=============================================================================='
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%Allocate_solid(nAsfac(iFish),nLchod(iFish))
        enddo
        !Use the object with the largest area as the reference object
        maxN  = maxloc(nAsfac, dim=1)
        Asfac = nAsfac(maxN)
        Lchod = nLchod(maxN)
        if (VBodies(maxN)%v_type.eq.1) then
            Lspan = sum(VBodies(maxN)%rbm%r_Lspan(:)+VBodies(maxN)%rbm%r_Rspan(:))/dble(VBodies(maxN)%rbm%nND)
        else
            Lspan = maxval(VBodies(maxN)%rbm%xyzful00(:,3))-minval(VBodies(maxN)%rbm%xyzful00(:,3))
        endif
        if((Lchod-1.0d0)<=1.0d-2)Lchod=1.0d0
        if((Lspan-1.0d0)<=1.0d-2)Lspan=1.0d0
        AR    = Lspan**2/Asfac
    end subroutine allocate_solid_memory

    subroutine Initialise_bodies(time,zDim,yDim,xDim,dh,g)
        implicit none
        real(8),intent(in):: time,dh,g(3)
        integer,intent(in):: zDim,yDim,xDim
        integer :: iFish
        m_zDim = zDim
        m_yDim = yDim
        m_xDim = xDim
        m_dh = dh
        m_Kb = ceiling(log(dble(m_zDim))/log(2.0))
        m_Jb = ceiling(log(dble(m_yDim))/log(2.0))
        m_Ib = ceiling(log(dble(m_xDim))/log(2.0))
        m_Kmask = ishft(1,m_Kb)-1
        m_Jmask = ishft(1,m_Jb+m_Kb)-1
        m_Imask = ishft(1,m_Ib+m_Jb+m_Kb)-1
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%Initialise(time,g)
            call VBodies(iFish)%Initialise(VBodies(iFish)%rbm%iBodyType)
        enddo
    end subroutine Initialise_bodies

    subroutine Initialise_(this,iBodyType)
        ! read beam central line file and allocate memory
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer, intent(in) :: iBodyType
        this%v_type = iBodyType
        if (this%v_type .eq. 1) then
            call this%PlateBuild()
        else
            write(*,*) 'not implemented body type', this%v_type
        endif
    end subroutine Initialise_

    subroutine Initialise_Calculate_Solid_params(Aref,Eref,Fref,Lref,Pref,Tref,Uref,uMax,Lthck)
        implicit none
        real(8),intent(in):: Aref,Eref,Fref,Lref,Pref,Tref,Uref
        real(8),intent(out):: Lthck,uMax
        integer:: iFish
        real(8):: nLthck(m_nFish)
        m_Aref = Aref
        m_Eref = Eref
        m_Fref = Fref
        m_Lref = Lref
        m_Pref = Pref
        m_Tref = Tref
        m_Uref = Uref
        uMax = 0.d0
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%calculate_angle_material(m_Lref, m_Uref, m_denIn, uMax, m_uuuIn, nLthck(iFish))
        enddo
        Lthck = maxval(nLthck)
    end subroutine Initialise_Calculate_Solid_params

    SUBROUTINE Solver(time,isubstep,deltat,subdeltat)
        implicit none
        integer,intent(in):: isubstep
        real(8),intent(in):: time,deltat,subdeltat
        integer:: iFish
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iFish)
        do iFish=1,m_nFish
            call VBodies(iFish)%rbm%structure(iFish,time,isubstep,deltat,subdeltat)
        enddo !do iFish=1,nFish
        !$OMP END PARALLEL DO
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Write_solid_v_bodies(time)
        implicit none
        real(8),intent(in):: time
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%Write_body(iFish,time)
        enddo
    end subroutine Write_solid_v_bodies

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write check point file for restarting simulation
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Write_cont(fid)
        implicit none
        integer,intent(in):: fid
        integer:: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%write_solid_temp(fid)
        enddo
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read check point file for restarting simulation
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Read_cont(fid)
        implicit none
        integer,intent(in):: fid
        integer:: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%read_solid_temp(fid)
        enddo
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_solid_field(time)
        implicit none
        real(8),intent(in):: time
        real(8):: timeTref
    !   -------------------------------------------------------
        integer:: i,iFish
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName
        !==================================================================================================
        timeTref = time/m_Tref
        write(fileName,'(I10)') nint(timeTref*1d5)
        fileName = adjustr(fileName)
        DO  i=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        END DO
    
        do iFish=1,m_nFish
            call VBodies(iFish)%rbm%write_solid(m_Lref,m_Uref,m_Aref,m_Fref,iFish,fileName)
        enddo
    !   =============================================
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write solid parameters for checking, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Write_solid_Check(fid)
        implicit none
        integer,intent(in):: fid
        integer:: iFish
        do iFish=1,m_nFish
            write(fid,'(A      )')'===================================='
            write(fid,'(A,I20.10)')'Fish number is',iFish
            write(fid,'(A      )')'===================================='
            call VBodies(iFish)%rbm%write_solid_params(fid)
        enddo
        do iFish=1,m_nFish
            call VBodies(iFish)%rbm%write_solid_materials(fid,iFish)
        enddo
    end subroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write solid data
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Write_solid_Data(fid,time,timeOutInfo,Asfac)
        implicit none
        integer,intent(in):: fid
        real(8),intent(in):: time,timeOutInfo,Asfac
        integer:: iFish
        do iFish=1,m_nFish
            call VBodies(iFish)%rbm%write_solid_info(fid,iFish,time,timeOutInfo,m_Tref,m_Lref,m_Uref,m_Aref,m_Fref,m_Pref,m_Eref,Asfac)
        enddo
    end subroutine

    subroutine Write_SampBodyNode(fid,time,numSampBody,SampBodyNode)
        implicit none
        integer,intent(in):: fid,numSampBody,SampBodyNode(numSampBody,m_nFish)
        real(8),intent(in):: time
        integer:: iFish
        do iFish=1,m_nFish
            call VBodies(iFish)%rbm%write_solid_SampBodyNode(fid,iFish,time,numSampBody,SampBodyNode(1:numSampBody,iFish),m_Tref,m_Lref,m_Uref,m_Aref)
        enddo !nFish
    end subroutine
    
    subroutine FSInteraction_force(xGrid,yGrid,zGrid,uuu,force)
        implicit none
        real(8),intent(in):: xGrid(m_xDim),yGrid(m_yDim),zGrid(m_zDim)
        real(8),intent(inout)::uuu(m_zDim,m_yDim,m_xDim,1:3)
        real(8),intent(out)::force(m_zDim,m_yDim,m_xDim,1:3)
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%UpdatePosVelArea()
        enddo
        call calculate_interaction_force(xGrid,yGrid,zGrid,uuu,force)
    end subroutine

    subroutine PlateUpdatePosVelArea_(this)
        !   compute displacement, velocity, area at surface element center
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        integer :: i,s,cnt,i1,i2
        real(8) :: tmpxyz(3), tmpvel(3), tmpdx(3)
        real(8) :: dh, left, len, dl, ls, area

        do i = 1,this%rbm%nEL
            i1 = this%rbm%ele(i,1)
            i2 = this%rbm%ele(i,2)
            tmpxyz = 0.5d0 * (this%rbm%xyzful(i1,1:3) + this%rbm%xyzful(i2,1:3))
            tmpvel = 0.5d0 * (this%rbm%velful(i1,1:3) + this%rbm%velful(i2,1:3))
            left = 0.5d0 * (this%rbm%r_Lspan(i1) + this%rbm%r_Lspan(i2))
            len = 0.5d0 * (this%rbm%r_Rspan(i1) + this%rbm%r_Rspan(i2)) + left
            dl = len / dble(this%rbm%r_Nspan(i))
            tmpdx = this%rbm%xyzful(i2,1:3) - this%rbm%xyzful(i1,1:3)
            dh = dsqrt(tmpdx(1)*tmpdx(1)+tmpdx(2)*tmpdx(2)+tmpdx(3)*tmpdx(3))
            area = dl * dh
            cnt = this%rtov(i) - 1
            do s=1,this%rbm%r_Nspan(i)
                ls = dl * (0.5d0 + dble(s-1)) - left
                this%v_Exyz(cnt+s, 1:3) = tmpxyz + this%rbm%r_dirc*ls
                this%v_Evel(cnt+s,1:3) = tmpvel
                this%v_Ea(cnt+s) = area
            enddo
        enddo
    endsubroutine PlateUpdatePosVelArea_

    subroutine UpdatePosVelArea_(this)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        if (this%v_type .eq. 1) then
            call this%PlateUpdatePosVelArea()
        else
            write(*,*) 'body type not implemented', this%v_type
        endif
    end subroutine UpdatePosVelArea_

    FUNCTION vtor_(this,x)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer, intent(in) :: x
        integer :: vtor_
        vtor_ = this%vtor(x)
    ENDFUNCTION vtor_
    FUNCTION rtov_(this,x)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer, intent(in) :: x
        integer :: rtov_
        rtov_ = this%rtov(x)
    ENDFUNCTION rtov_

    subroutine getuuu_(this,uuu)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in)::uuu(m_zDim,m_yDim,m_xDim,1:3)
        integer:: iEL,ix(-1:2),jy(-1:2),kz(-1:2),x,y,z,num
        integer:: index
        allocate(this%array_i(this%v_nelmts*64),this%array_u(this%v_nelmts*64,3))
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(iEL,1:4)
            jy = this%v_Ei(iEL,5:8)
            kz = this%v_Ei(iEL,9:12)
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        call findindex(ix(x),jy(y),kz(z),index)
                        num = (iEl-1)*64 + (x+1)*16 + (y+1)*4 + (z+1) + 1
                        this%array_i(num) = index
                        this%array_u(num,1:3) = uuu(kz(z),jy(y),ix(x),1:3)
                    enddo
                enddo
            enddo
        enddo
        call sort_indices(this%array_i, this%array_u, this%v_nelmts*64)
    end subroutine

    subroutine sort_indices(A, B, n)
        integer:: n, A(n)
        real(8):: B(n,1:3)
        integer:: indices(n)
        integer:: i, j, temp

        do i = 1, size(A)
            indices(i) = i
        end do

        do i = 1, size(A) - 1
            do j = i + 1, size(A)
                if (A(indices(i)) > A(indices(j))) then
                    temp = indices(i)
                    indices(i) = indices(j)
                    indices(j) = temp
                end if
            end do
        end do

        call reorder_array()
        contains
        subroutine reorder_array()
            integer:: tempA(n)
            real(8):: tempB(n,1:3)
    
            do i = 1, size(A)
                tempA(i) = A(indices(i))
                tempB(i,1:3) = B(indices(i),1:3)
            end do
    
            A = tempA
            B = tempB
        end subroutine reorder_array
    end subroutine sort_indices

    subroutine setuuu_(this,uuu)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(inout)::uuu(m_zDim,m_yDim,m_xDim,1:3)
        integer:: index,i,j,k,ijk
        allocate(this%array_i(this%v_nelmts*64),this%array_u(this%v_nelmts*64,3))
        do  index = 1,this%v_nelmts*64
            ijk = this%array_i(index)
            call findijk(ijk,i,j,k)
            uuu(k,j,i,1:3) = this%array_u(index,1:3)
        enddo
        call sort_indices(this%array_i, this%array_u, this%v_nelmts*64)
    end subroutine

    subroutine findindex(i,j,k,index)
        IMPLICIT NONE
        integer,intent(out):: index
        integer:: i,j,k
        i = ishft((i-1),M_Jb+M_Kb)
        j = ishft((j-1),M_Kb)
        index = ior((k-1),ior(j,i))
    end subroutine
    subroutine findijk(index,i,j,k)
        IMPLICIT NONE
        integer,intent(in):: index
        integer,intent(out):: i,j,k
        i = iand(index,m_Imask)+1
        j = iand(index,m_Jmask)+1
        k = iand(index,m_Kmask)+1
    end subroutine

    subroutine UpdateElmtInterp_(this,xGrid,yGrid,zGrid)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: xGrid(m_xDim),yGrid(m_yDim),zGrid(m_zDim)
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2)
        real(8)::x0,y0,z0,detx,dety,detz,invdh
        integer::i0,j0,k0,iEL,x,y,z,i,j,k
        !==================================================================================================
        invdh = 1.D0/m_dh
        call my_minloc(this%v_Exyz(1,1), xGrid, m_xDim, .false., i0)
        call my_minloc(this%v_Exyz(1,2), yGrid, m_yDim, .false., j0)
        call my_minloc(this%v_Exyz(1,3), zGrid, m_zDim, .false., k0)
        x0 = xGrid(i0)
        y0 = yGrid(j0)
        z0 = zGrid(k0)
        ! compute the velocity of IB nodes at element center
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,k,x,y,z,rx,ry,rz,detx,dety,detz,ix,jy,kz)
        do  iEL=1,this%v_nelmts
            call minloc_fast(this%v_Exyz(iEL,1), x0, i0, invdh, i, detx)
            call minloc_fast(this%v_Exyz(iEL,2), y0, j0, invdh, j, dety)
            call minloc_fast(this%v_Exyz(iEL,3), z0, k0, invdh, k, detz)
            call trimedindex(i, m_xDim, ix, m_boundaryConditions(1:2))
            call trimedindex(j, m_yDim, jy, m_boundaryConditions(3:4))
            call trimedindex(k, m_zDim, kz, m_boundaryConditions(5:6))
            this%v_Ei(iEL,1:4) = ix
            this%v_Ei(iEL,5:8) = jy
            this%v_Ei(iEL,9:12) = kz
            do x=-1,2
                rx(x)=Phi(dble(x)-detx)
            enddo
            do y=-1,2
                ry(y)=Phi(dble(y)-dety)
            enddo
            do z=-1,2
                rz(z)=Phi(dble(z)-detz)
            enddo
            this%v_Ew(iEL,1:4) = rx
            this%v_Ew(iEL,5:8) = ry
            this%v_Ew(iEL,9:12) = rz
        enddo
        !$OMP END PARALLEL DO

        contains
        ! return the array(index) <= x < array(index+1)
        ! assum uniform grid around x(x=i0) = x0
        SUBROUTINE minloc_fast(x_, x0_, i0_, invdh_, index_, offset_)
            implicit none
            integer, intent(in):: i0_
            real(8), intent(in):: x_, x0_, invdh_
            integer, intent(out):: index_
            real(8), intent(out):: offset_
            offset_ = (x_ - x0_)*invdh_
            index_ = floor(offset_)
            offset_ = offset_ - dble(index_)
            index_ = index_ + i0_
        END SUBROUTINE
        FUNCTION Phi(x_)
            IMPLICIT NONE
            real(8)::Phi,x_,r
            r=dabs(x_)
            if(r<1.0d0)then
                Phi=(3.d0-2.d0*r+dsqrt( 1.d0+4.d0*r*(1.d0-r)))*0.125d0
            elseif(r<2.0d0)then
                Phi=(5.d0-2.d0*r-dsqrt(-7.d0+4.d0*r*(3.d0-r)))*0.125d0
            else
                Phi=0.0d0
            endif
        ENDFUNCTION Phi
        SUBROUTINE trimedindex(i_, xDim_, ix_, boundaryConditions_)
            USE BoundCondParams
            implicit none
            integer, intent(in):: boundaryConditions_(1:2)
            integer, intent(in):: i_, xDim_
            integer, intent(out):: ix_(-1:2)
            integer:: k_
            do k_=-1,2
                ix_(k_) = i_ + k_
                if (ix_(k_)<1) then
                    if(boundaryConditions_(1).eq.Periodic) then
                        ix_(k_) = ix_(k_) + xDim_
                    else if((boundaryConditions_(1).eq.SYMMETRIC .or. boundaryConditions_(1).eq.wall) .and. ix_(k_).eq.0) then
                        ix_(k_) = 2
                    else
                        write(*,*) 'index out of xmin bound', ix_(k_)
                        stop
                    endif
                else if(ix_(k_)>xDim_) then
                    if(boundaryConditions_(2).eq.Periodic) then
                        ix_(k_) = ix_(k_) - xDim_
                    else if((boundaryConditions_(2).eq.SYMMETRIC .or. boundaryConditions_(2).eq.wall) .and. ix_(k_).eq.xDim_+1) then
                        ix_(k_) = xDim_ - 1
                    else
                        write(*,*) 'index out of xmax bound', ix_(k_)
                        stop
                    endif
                endif
            enddo
        END SUBROUTINE trimedindex
        SUBROUTINE my_minloc(x_, array_, len_, uniform_, index_) ! return the array(index) <= x < array(index+1)
            implicit none
            integer:: len_, index_, count, step, it
            real(8):: x_, array_(len_)
            logical:: uniform_
            if (.not.uniform_) then
                if (x_<array_(1) .or. x_>array_(len_)) then
                    write(*, *) 'index out of bounds when searching my_minloc', x_, '[', array_(1), array_(len_), ']'
                    stop
                endif
                index_ = 1
                count = len_
                do while(count > 0)
                    step = count / 2
                    it = index_ + step
                    if (array_(it) < x_) then
                        index_ = it + 1
                        count = count - (step + 1)
                    else
                        count = step
                    endif
                enddo
                if (array_(index_)>x_) then
                    index_ = index_ - 1
                endif
            else
                index_ = 1 + int((x_ - array_(1))/(array_(len_)-array_(1))*dble(len_-1))
                !int -1.1 -> -1; 1.1->1; 1.9->1
                if (index_<1 .or. index_>len_) then
                    write(*, *) 'index out of bounds when searching my_minloc', x_, '[', array_(1), array_(len_), ']'
                    stop
                endif
            endif
        END SUBROUTINE
    end subroutine UpdateElmtInterp_

    SUBROUTINE calculate_interaction_force(xGrid,yGrid,zGrid,uuu,force)
        ! calculate elements interaction force using IB method
        IMPLICIT NONE
        real(8),intent(in):: xGrid(m_xDim),yGrid(m_yDim),zGrid(m_zDim)
        real(8),intent(inout)::uuu(m_zDim,m_yDim,m_xDim,1:3)
        real(8),intent(out)::force(m_zDim,m_yDim,m_xDim,1:3)
        !================================
        integer:: iFish
        integer:: iterLBM
        real(8):: dmaxLBM,dsum
        real(8)::tol,ntol
        ! update virtual body shape and velocity
        do iFish = 1, m_nFish
            call VBodies(iFish)%UpdateElmtInterp(xGrid,yGrid,zGrid)
            call VBodies(iFish)%getuuu(uuu)
            VBodies(iFish)%v_Eforce = 0.0d0
        enddo
        ! calculate interaction force using immersed-boundary method
        iterLBM=0
        dmaxLBM=1d10
        do  while( iterLBM<m_ntolLBM .and. dmaxLBM>m_dtolLBM)
            dmaxLBM = 0.d0
            dsum=0.0d0
            do iFish=1,m_nFish
                call VBodies(iFish)%PenaltyForce(tol,ntol)
                dmaxLBM = dmaxLBM + tol
                dsum = dsum + ntol
            enddo
            dmaxLBM=dmaxLBM/(dsum * m_Uref)
            iterLBM=iterLBM+1
        enddo
        ! update body load and fluid force
        do iFish=1,m_nFish
            VBodies(iFish)%rbm%extful = 0.0d0
            ! to do, consider gravity
        enddo
        do iFish=1,m_nFish
            call VBodies(iFish)%setuuu(uuu)
            call VBodies(iFish)%FluidVolumeForce(force)
        enddo
    END SUBROUTINE

    SUBROUTINE FluidVolumeForce_(this,force)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(out)::force(m_zDim,m_yDim,m_xDim,1:3)
        !==================================================================================================
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2),forcetemp(1:3)
        real(8):: forceElemTemp(3),invh3
        !==================================================================================================
        integer::x,y,z,iEL,i1,i2
        !==================================================================================================
        invh3 = (1.d0/m_dh)**3
        ! compute the velocity of IB nodes at element center
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(iEL,1:4)
            jy = this%v_Ei(iEL,5:8)
            kz = this%v_Ei(iEL,9:12)
            rx = this%v_Ew(iEL,1:4)
            ry = this%v_Ew(iEL,5:8)
            rz = this%v_Ew(iEL,9:12)
            forceElemTemp = this%v_Eforce(iEL,1:3)
            ! update beam load, momentum is not included
            i1=this%rbm%ele(this%vtor(iEL),1)
            i2=this%rbm%ele(this%vtor(iEL),2)
            this%rbm%extful(i1,1:3) = this%rbm%extful(i1,1:3) + 0.5d0 * forceElemTemp
            this%rbm%extful(i2,1:3) = this%rbm%extful(i2,1:3) + 0.5d0 * forceElemTemp
            forceElemTemp(1:3) = forceElemTemp(1:3) * invh3
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        forceTemp(1:3) = -forceElemTemp(1:3)*rx(x)*ry(y)*rz(z)
                        ! add flow body force
                        force(kz(z),jy(y),ix(x),1:3) = force(kz(z),jy(y),ix(x),1:3) + forceTemp(1:3)
                    enddo
                enddo
            enddo
        enddo
    END SUBROUTINE FluidVolumeForce_

    SUBROUTINE PenaltyForce_(this,tolerance,ntolsum)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(out)::tolerance, ntolsum
        !==================================================================================================
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2),forcetemp(1:3)
        real(8):: velElem(3),velElemIB(this%v_nelmts,3),forceElemTemp(this%v_nelmts,3),invh3
        real(8):: uuu(3)
        !==================================================================================================
        integer::x,y,z,iEL,index,index0
        !==================================================================================================
        invh3 = 0.5d0*m_dt*(1.d0/m_dh)**3/m_denIn
        tolerance = 0.d0
        ntolsum = dble(this%v_nelmts)
        ! compute the velocity of IB nodes at element center
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(iEL,1:4)
            jy = this%v_Ei(iEL,5:8)
            kz = this%v_Ei(iEL,9:12)
            rx = this%v_Ew(iEL,1:4)
            ry = this%v_Ew(iEL,5:8)
            rz = this%v_Ew(iEL,9:12)
            velElemIB(iEL,1:3)=0.0d0
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        call findindex(ix(x),jy(y),kz(z),index0)
                        call my_minloc(index0,this%array_i,this%v_nelmts*64,index)
                        uuu(1:3) = this%array_u(index,1:3)
                        velElemIB(iEL,1:3)=velElemIB(iEL,1:3)+uuu(1:3)*rx(x)*ry(y)*rz(z)
                    enddo
                enddo
            enddo
        enddo
        do  iEL=1,this%v_nelmts
            velElem = this%v_Evel(iEL,1:3)
            forceElemTemp(iEL,1:3) = -m_Pbeta* 2.0d0*m_denIn*(velElem(1:3)-velElemIB(iEL,1:3))*this%v_Ea(iEL)
            if ((.not. IEEE_IS_FINITE(forceElemTemp(iEL,1))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(iEL,2))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(iEL,3)))) then
                write(*, *) 'Nan found in forceElemTemp', forceElemTemp
                ix(0) = this%v_Ei(iEL,2)
                jy(0) = this%v_Ei(iEL,6)
                kz(0) = this%v_Ei(iEL,10)
                write(*, *) 'Nan found at (ix, jy, kz)', ix(0), jy(0), kz(0)
                stop
            endif
            tolerance = tolerance + dsqrt(sum((velElem(1:3)-velElemIB(iEL,1:3))**2))
            this%v_Eforce(iEL,1:3) = this%v_Eforce(iEL,1:3) + forceElemTemp(iEL,1:3)
            forceElemTemp(iEL,1:3) = forceElemTemp(iEL,1:3) * invh3
        enddo
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(iEL,1:4)
            jy = this%v_Ei(iEL,5:8)
            kz = this%v_Ei(iEL,9:12)
            rx = this%v_Ew(iEL,1:4)
            ry = this%v_Ew(iEL,5:8)
            rz = this%v_Ew(iEL,9:12)
            ! correct velocity
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        forceTemp(1:3) = -forceElemTemp(iEL,1:3)*rx(x)*ry(y)*rz(z)
                        call findindex(ix(x),jy(y),kz(z),index0)
                        call my_minloc(index0,this%array_i,this%v_nelmts*64,index)
                        uuu(1:3) = this%array_u(index,1:3)
                        uuu(1:3)  = uuu(1:3)+forceTemp(1:3)
                    enddo
                enddo
            enddo
        enddo
        contains
        SUBROUTINE my_minloc(x_, array_, len_, index_) ! return the array(index) <= x < array(index+1)
            implicit none
            integer:: x_, array_(len_),len_, index_, count, step, it
                if (x_<array_(1) .or. x_>array_(len_)) then
                    write(*, *) 'index out of bounds when searching my_minloc', x_, '[', array_(1), array_(len_), ']'
                    stop
                endif
                index_ = 1
                count = len_
                do while(count > 0)
                    step = count / 2
                    it = index_ + step
                    if (array_(it) < x_) then
                        index_ = it + 1
                        count = count - (step + 1)
                    else
                        count = step
                    endif
                enddo
                if (array_(index_)>x_) then
                    index_ = index_ - 1
                endif
        END SUBROUTINE
    END SUBROUTINE PenaltyForce_

    subroutine PlateBuild_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer:: i, j
        allocate(this%rtov(this%rbm%nEL+1))
        this%v_nelmts = 0
        do i=1,this%rbm%nEL
            this%rtov(i) = this%v_nelmts + 1
            this%v_nelmts = this%v_nelmts + this%rbm%r_Nspan(i)
        enddo
        this%rtov(this%rbm%nEL+1) = this%v_nelmts + 1
        allocate(this%vtor(this%v_nelmts))
        do i=1,this%rbm%nEL
            do j= this%rtov(i), this%rtov(i+1)-1
                this%vtor(j) = i
            enddo
        enddo
        allocate(this%v_Exyz(this%v_nelmts,3), this%v_Ea(this%v_nelmts), this%v_Eforce(this%v_nelmts,3))
        allocate(this%v_Evel(this%v_nelmts,3), this%v_Ei(this%v_nelmts,12), this%v_Ew(this%v_nelmts,12))
    end subroutine PlateBuild_

    subroutine Write_body_(this,iFish,time)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer,intent(in) :: iFish
        real(8),intent(in) :: time
        if(this%v_type.eq.1) then
            call this%PlateWrite_body(iFish,time)
        else
            write(*, *) "body type not implemented", this%v_type
        endif
    endsubroutine Write_body_

    subroutine PlateWrite_body_(this,iFish,time)
        ! to do: generate a temporary mesh
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer,intent(in) :: iFish
        real(8),intent(in) :: time
        !   -------------------------------------------------------
        real(8):: timeTref
        integer:: i, i1, i2
        real(8) :: tmpxyz(3)
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName,idstr
        integer,parameter:: idfile=100
        !==========================================================================
        timeTref = time/m_Tref
        write(fileName,'(I10)') nint(timeTref*1d5)
        fileName = adjustr(fileName)
        do  I=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo
        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        !==========================================================================
        open(idfile, FILE='./DatBodySpan/BodyFake'//trim(idstr)//'_'//trim(filename)//'.dat')
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        write(idfile, '(A,I7,A,I7,A)') 'ZONE N=',2*this%rbm%nND,', E=',this%rbm%nEL,', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'

        do i = 1,this%rbm%nNd
            tmpxyz = this%rbm%xyzful(i,1:3)
            write(idfile, *) (tmpxyz - this%rbm%r_Lspan(i) * this%rbm%r_dirc)/m_Lref
            write(idfile, *) (tmpxyz + this%rbm%r_Rspan(i) * this%rbm%r_dirc)/m_Lref
        enddo
        do  i=1,this%rbm%nEL
            i1 = this%rbm%ele(i,1)
            i2 = this%rbm%ele(i,2)
            write(idfile, *) 2*i1 -1, 2*i1,2*i2,2*i2-1
        enddo
        close(idfile)
    end subroutine PlateWrite_body_

end module SolidBody