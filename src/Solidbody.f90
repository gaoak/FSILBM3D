module SolidBody
    use SolidSolver
    implicit none
    private
    ! Immersed boundary method parameters
    integer:: m_nFish, m_maxIterIB, m_zDim, m_yDim, m_xDim
    real(8):: m_dtolLBM, m_Pbeta, m_dt, m_h, m_denIn, m_Eref, m_Fref, m_Lref, m_Pref, m_Tref, m_Uref
    integer:: m_boundaryConditions(1:6)
    ! nFish     number of bodies
    ! maxIterIB maximum number of iterations for IB force calculation
    ! dtolIBM   tolerance for IB force calculation
    ! Pbeta     coefficient in penalty force calculation
    public :: VirtualBody,Initialise_bodies,Write_solid_bodies,FSInteraction_force, &
              Write_cont,Read_cont,Write_BeamChk,Write_params
    type :: VirtualBody
        type(BeamSolver):: rbm
        !!!virtual infomation
        !!!virtual body surface
        integer :: v_nelmts,v_type ! v_type 0 (solid body), 1 (plate), 2 (rod)
        real(8), allocatable :: v_Exyz(:, :) ! element center (x, y, z)
        real(8), allocatable :: v_Ea(:) ! element area
        !area center with equal weight on both sides
        real(8), allocatable :: v_Evel(:, :)
        integer(2), allocatable :: v_Ei(:, :) ! element stencial integer index [ix-1,ix,ix1,ix2, iy-1,iy,iy1,iy2, iz-1,iz,iz1,iz2]
        real(4), allocatable :: v_Ew(:, :) ! element stential weight [wx-1, wx, wx1, wx2, wy-1, wy, wy1, wy2, wz-1, wz, wz1, wz2]
        !calculated using central linear and angular velocities
        integer,allocatable :: vtor(:)!of size fake_npts
        integer,allocatable :: rtov(:)! of size real_npts+1
    contains
        procedure :: Initialise => Initialise_
        procedure :: PlateBuild => PlateBuild_
        procedure :: UpdatePosVelArea => UpdatePosVelArea_
        procedure :: PlateUpdatePosVelArea => PlateUpdatePosVelArea_
        procedure :: Write_body => Write_body_
        procedure :: PlateWrite_body => PlateWrite_body_
        procedure :: UpdateElmtInterp => UpdateElmtInterp_
        procedure :: PenaltyForce => PenaltyForce_
    end type VirtualBody
    type(VirtualBody), allocatable :: VBodies(:)
  contains
    
    subroutine Initialise_bodies(zDim,yDim,xDim,nFish,filenames,iBodyModel,maxIterIB,dtolLBM,Pbeta,dt,h,denIn,Eref,Fref,Lref,Pref,Tref,Uref,BCs)
        USE simParam
        implicit none
        integer,intent(in)::zDim,yDim,xDim,nFish,maxIterIB,iBodyModel(nFish),BCs(6)
        real(8),intent(in):: dtolLBM,Pbeta,dt,h,denIn,Eref,Fref,Lref,Pref,Tref,Uref
        character(LEN=40), intent(in):: filenames(nFish)
        integer :: iFish
        m_zDim = zDim
        m_yDim = yDim
        m_xDim = xDim
        m_nFish = nFish
        m_maxIterIB = maxIterIB
        m_dtolLBM = dtolLBM
        m_Pbeta = Pbeta
        m_dt = dt
        m_h = h
        m_denIn = denIn
        m_Uref = Uref
        m_Lref = Lref
        m_Tref = Tref
        m_Eref = Eref
        m_Fref = Fref
        m_Lref = Lref
        m_Pref = Pref
        m_Tref = Tref
        m_Uref = Uref
        m_boundaryConditions(1:6) = BCs(1:6)
        allocate(VBodies(nFish))
        allocate(nAsfac(1:m_nFish),nLchod(1:m_nFish),lentemp(1:m_nFish))
        do iFish = 1,nFish
            call VBodies(iFish)%Initialise(filenames(iFish),iBodyModel(iFish))
            call VBody(iFish)%rbm%Allocate_solid(this,filenames(iFish),nAsfac(iFish),nLchod(iFish),lentemp(iFish))
        enddo
        maxN  = maxloc(nAsfac, dim=1)
        Asfac = nAsfac(maxN)
        Lchod = nLchod(maxN)
        if((Lchod-1.0d0)<=1.0d-2)Lchod=1.0d0
        if((maxval(Lspan)-1.0d0)<=1.0d-2)Lspan(maxloc(Lspan))=1.0d0
        AR    = maxval(Lspan)**2/Asfac
    end subroutine Initialise_bodies

    subroutine Initialise_(this,filename,iBodyModel)
        ! read beam central line file and allocate memory
        implicit none
        class(VirtualBody), intent(inout) :: this
        character(LEN=40), intent(in):: filename
        integer, intent(in) :: iBodyModel
        !call this%rbm%Allocate_solid(filename)
        this%v_type = iBodyModel
        if (this%v_type .eq. 1) then
            call this%PlateBuild()
        else
            write(*,*) 'not implemented body type', this%v_type
        endif
    end subroutine Initialise_

    subroutine FSInteraction_force(xGrid,yGrid,zGrid,uuu,den,force)
        implicit none
        real(8),intent(in):: xGrid(m_xDim),yGrid(m_yDim),zGrid(m_zDim),den(m_zDim,m_yDim,m_xDim)
        real(8),intent(inout)::uuu(m_zDim,m_yDim,m_xDim,1:3)
        real(8),intent(out)::force(m_zDim,m_yDim,m_xDim,1:3)
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%UpdatePosVelArea()
        enddo
        call calculate_interaction_force(xGrid,yGrid,zGrid,uuu,den,force)
    end subroutine

    subroutine PlateUpdatePosVelArea_(this)
        !   compute displacement, velocity, area at surface element center
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        integer :: i,s,cnt
        real(8) :: tmpxyz(3), tmpvel(3), tmpdx(3)
        real(8) :: dh, left, len, dl, ls, area

        do i = 1,this%rbm%nEL
            tmpxyz = 0.5d0 * (this%rbm%xyzful(i,1:3) + this%rbm%xyzful(i+1,1:3))
            tmpvel = 0.5d0 * (this%rbm%velful(i,1:3) + this%rbm%velful(i+1,1:3))
            left = 0.5d0 * (this%rbm%r_Lspan(i) - this%rbm%r_Lspan(i+1))
            len = 0.5d0 * (this%rbm%r_Rspan(i) + this%rbm%r_Rspan(i+1)) - left
            dl = len / dble(this%rbm%r_Nspan(i))
            tmpdx = this%rbm%xyzful(i+1,1:3) - this%rbm%xyzful(i,1:3)
            dh = dsqrt(tmpdx(1)*tmpdx(1)+tmpdx(2)*tmpdx(2)+tmpdx(3)*tmpdx(3))
            area = dl * dh
            cnt = this%rtov(i) - 1
            do s=1,this%rbm%r_Nspan(i)
                ls = left + dl * (0.5d0 + dble(s-1))
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

    subroutine UpdateElmtInterp_(this,xGrid,yGrid,zGrid)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: xGrid(m_xDim),yGrid(m_yDim),zGrid(m_zDim)
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2),Phi
        real(8)::x0,y0,z0,detx,dety,detz,invdh
        integer::i0,j0,k0,iEL,x,y,z,i,j,k
        !==================================================================================================
        invdh = 1.D0/m_h
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
                    endif
                else if(ix_(k_)>xDim_) then
                    if(boundaryConditions_(2).eq.Periodic) then
                        ix_(k_) = ix_(k_) - xDim_
                    else if((boundaryConditions_(2).eq.SYMMETRIC .or. boundaryConditions_(2).eq.wall) .and. ix_(k_).eq.xDim_+1) then
                        ix_(k_) = xDim_ - 1
                    else
                        write(*,*) 'index out of xmax bound', ix_(k_)
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

    SUBROUTINE calculate_interaction_force(xGrid,yGrid,zGrid,uuu,den,force)
        ! calculate elements interaction force using IB method
        IMPLICIT NONE
        real(8),intent(in):: xGrid(m_xDim),yGrid(m_yDim),zGrid(m_zDim),den(m_zDim,m_yDim,m_xDim)
        real(8),intent(inout)::uuu(m_zDim,m_yDim,m_xDim,1:3)
        real(8),intent(out)::force(m_zDim,m_yDim,m_xDim,1:3)
        !================================
        integer:: iFish
        integer:: i,j,k,iEL,nt,iterLBM
        real(8):: dmaxLBM,dsum
        real(8)::tol,tolsum, ntol
        !================================
        do iFish = 1, m_nFish
            call VBodies(iFish)%UpdateElmtInterp(xGrid,yGrid,zGrid)
        enddo
        tolsum=0.d0
        iterLBM=0
        do iFish=1,m_nFish
            VBodies(iFish)%rbm%lodful = 0.0d0
        enddo
        do  while( iterLBM<m_maxIterIB .and. dmaxLBM>m_dtolLBM)
            dmaxLBM=0.0d0
            dsum=0.0d0
            do iFish=1,m_nFish
                call VBodies(iFish)%PenaltyForce(uuu,den,force,tol,ntol)
                tolsum = tolsum + dsqrt(tol)
                dsum = dsum + m_Uref*ntol
            enddo
            dmaxLBM=tolsum/dsum
            iterLBM=iterLBM+1
        enddo
    END SUBROUTINE

    SUBROUTINE PenaltyForce_(this,uuu,den,force,tolerance,ntolsum)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: den(m_zDim,m_yDim,m_xDim)
        real(8),intent(inout)::uuu(m_zDim,m_yDim,m_xDim,1:3)
        real(8),intent(out)::force(m_zDim,m_yDim,m_xDim,1:3)
        real(8),intent(out)::tolerance, ntolsum
        !==================================================================================================
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2),Phi,invdh,forcetemp(1:3)
        real(8):: velElem(3),velElemIB(3),forceElemTemp(3)
        !==================================================================================================
        real(8)::x0,y0,z0,detx,dety,detz
        integer::x,y,z,iEL
        !==================================================================================================
        tolerance = 0.d0
        ntolsum = dble(this%v_nelmts)
        ! compute the velocity of IB nodes at element center
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,x,y,z,rx,ry,rz,ix,jy,kz,velElem,velElemIB,forceElemTemp,forceTemp)  REDUCTION(+:tolerance)
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(iEL,1:4)
            jy = this%v_Ei(iEL,5:8)
            kz = this%v_Ei(iEL,9:12)
            rx = this%v_Ew(iEL,1:4)
            ry = this%v_Ew(iEL,5:8)
            rz = this%v_Ew(iEL,9:12)
            velElemIB(1:3)=0.0d0
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        velElemIB(1:3)=velElemIB(1:3)+uuu(kz(z),jy(y),ix(x),1:3)*rx(x)*ry(y)*rz(z)
                    enddo
                enddo
            enddo
            velElem = this%v_Evel(iEL,1:3)
            forceElemTemp(1:3) = -m_Pbeta* 2.0d0*m_denIn*(velElem-velElemIB)/m_dt*this%v_Ea(iEL)*m_h
            if ((.not. IEEE_IS_FINITE(forceElemTemp(1))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(2))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(3)))) then
                write(*, *) 'Nan found in forceElemTemp', forceElemTemp
                write(*, *) 'Nan found at (ix, jy, kz)', ix(0), jy(0), kz(0)
                stop
            endif
            tolerance = tolerance + sum((velElem(1:3)-velElemIB(1:3))**2)
            !$OMP critical
            ! update beam load, momentum is not included
            this%rbm%lodful(this%vtor(iEL)  ,1:3) = this%rbm%lodful(this%vtor(iEL)  ,1:3) + 0.5d0 * forceElemTemp
            this%rbm%lodful(this%vtor(iEL)+1,1:3) = this%rbm%lodful(this%vtor(iEL)+1,1:3) + 0.5d0 * forceElemTemp
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        forceTemp(1:3) = -forceElemTemp(1:3)*rx(x)*ry(y)*rz(z)
                        ! correct velocity
                        uuu(kz(z),jy(y),ix(x),1:3)  = uuu(kz(z),jy(y),ix(x),1:3)+0.5d0*m_dt*forceTemp(1:3)/den(kz(z),jy(y),ix(x))
                        ! add flow body force
                        force(kz(z),jy(y),ix(x),1:3) = force(kz(z),jy(y),ix(x),1:3) + forceTemp(1:3)
                    enddo
                enddo
            enddo
            !$OMP end critical
        enddo
        !$OMP END PARALLEL DO
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
        allocate(this%v_Exyz(this%v_nelmts,3), this%v_Ea(this%v_nelmts))
        allocate(this%v_Evel(this%v_nelmts,3))
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
        integer:: i,j,r,Nspanpts,ElmType,n_theta,Ea_A,Ea_D
        real(8) :: tmpxyz(3)
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName,idstr
        integer,parameter:: idfile=100
        real(8) :: invl0,cos_xyz(3),dspan,xyz1(3),xyz2(3),low_L(3),low_R(3),upp_R(3),upp_L(3)
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
            write(idfile, *) (tmpxyz - this%rbm%r_Rspan(i) * this%rbm%r_dirc)/m_Lref
            write(idfile, *) (tmpxyz + this%rbm%r_Lspan(i) * this%rbm%r_dirc)/m_Lref
        enddo
        do  i=1,this%rbm%nEL
            write(idfile, *) 2*i -1, 2*i,2*i+2,2*i+1
        enddo
        close(idfile)
    end subroutine PlateWrite_body_

    subroutine Write_solid_bodies(time,Lref,Tref)
        implicit none
        real(8) :: Lref,time,Tref
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%Write_body(iFish,time)
        enddo
    end subroutine Write_solid_bodies

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE Write_BeamChk(time)
        implicit none
        real(8):: time
    !   -------------------------------------------------------
        integer:: i,iFish
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName
        !==================================================================================================
        write(fileName,'(I10)') nint(time*1d5)
        fileName = adjustr(fileName)
        DO  i=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        END DO
    
        do iFish=1,m_nFish
            call VBody(iFish)%rbm%write_solid(m_Lref,m_Uref,m_Aref,m_Fref,iFish,fileName)
        enddo
    !   =============================================
    END SUBROUTINE
end module SolidBody