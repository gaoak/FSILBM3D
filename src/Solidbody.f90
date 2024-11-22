module SolidBody
    use SolidSolver
    implicit none
    private
    ! Immersed boundary method parameters
    integer:: m_nFish, m_maxIterIB
    real(8):: m_dtolLBM, m_Pbeta, m_dt, m_h, m_denIn, m_Uref
    integer:: m_boundaryConditions(1:6)
    ! nFish     number of bodies
    ! maxIterIB maximum number of iterations for IB force calculation
    ! dtolIBM   tolerance for IB force calculation
    ! Pbeta     coefficient in penalty force calculation
    public :: VirtualBody,Initialise_bodies,Write_solid_bodies,FSInteraction_force
    type :: VirtualBody
        !!!virtual infomation
        integer :: r_npts,r_nelmts
        !!!virtual body surface
        integer :: v_npts,v_nelmts,v_type ! v_type 0 (solid body), 1 (plate), 2 (rod)
        real(8) :: v_dirc(3)
        real(8), allocatable :: v_Exyz(:, :) ! element center (x, y, z)
        integer, allocatable :: v_elmt(:, :) ! element ID
        integer(2), allocatable :: v_Ei(:, :) ! element stencial integer index [ix-1,ix,ix1,ix2, iy-1,iy,iy1,iy2, iz-1,iz,iz1,iz2]
        real(4), allocatable :: v_Ew(:, :) ! element stential weight [wx-1, wx, wx1, wx2, wy-1, wy, wy1, wy2, wz-1, wz, wz1, wz2]
        real(8), allocatable :: v_Ea(:) ! element area
        !area center with equal weight on both sides
        real(8), allocatable :: v_Evel(:, :)
        !calculated using central linear and angular velocities
        integer(2),allocatable :: vtor(:)!of size fake_npts
        integer(2),allocatable :: rtov(:,:)! of size real_npts, 1:2
    contains
        procedure :: Initialise => Initialise_
        procedure :: BuildPlate => BuildPlate_
        procedure :: Beam_ReadUnstructured => Beam_ReadUnstructured_
        procedure :: Readreal => Readreal_
        procedure :: Unstruc_Virtual_Elmts => Unstruc_Virtual_Elmts_
        procedure :: Struc_Virtual_Elmts => Struc_Virtual_Elmts_
        procedure :: Struc_Rib_Elmts => Struc_Rib_Elmts_
        procedure :: Struc_Cylinder_Elmts => Struc_Cylinder_Elmts_

        procedure :: UpdateElmtInfo => UpdateElmtInfo_
        procedure :: UpdateElmtPosVelArea => UpdateElmtPosVelArea_
        procedure :: UpdateElmtInterp => UpdateElmtInterp_
        procedure :: PenaltyForce => PenaltyForce_

        procedure :: Write_body => Write_body_
    end type VirtualBody
    type(VirtualBody), allocatable :: VBodies(:)
  contains
    
    subroutine Initialise_bodies(nFish,iBodyModel,maxIterIB,dtolLBM,Pbeta,dt,h,denIn,Uref,BCs,beams)
        implicit none
        integer,intent(in)::nFish,maxIterIB,iBodyModel(nFish),BCs(6)
        real(8),intent(in):: dtolLBM,Pbeta,dt,h,denIn,Uref
        type(BeamSolver), intent(in):: beams(nFish)
        integer :: iFish
        m_nFish = nFish
        m_maxIterIB = maxIterIB
        m_dtolLBM = dtolLBM
        m_Pbeta = Pbeta
        m_dt = dt
        m_h = h
        m_denIn = denIn
        m_Uref = Uref
        m_boundaryConditions(1:6) = BCs(1:6)
        allocate(VBodies(nFish))
        do iFish = 1,nFish
            call VBodies(iFish)%Initialise(iBodyModel(iFish), beams(iFish))
        enddo
    end subroutine Initialise_bodies

    subroutine Initialise_(this,iBodyModel,beam)
        ! read beam central line file and allocate memory
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer, intent(in) :: iBodyModel
        type(BeamSolver), intent(in):: beam
        this%v_type = iBodyModel
        if (this%v_type .eq. 1) then
            call this%BuildPlate(beam)
        else
            write(*,*) 'not implemented body type', this%v_type
        endif
    end subroutine Initialise_

    subroutine FSInteraction_force(dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
        implicit none
        real(8),intent(in):: dt,dh,denIn,Uref
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: xGrid(xDim),yGrid(yDim),zGrid(zDim),den(zDim,yDim,xDim)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%UpdateElmtInfo(dh,zDim,yDim,xDim,xGrid,yGrid,zGrid)
        enddo
        call calculate_interaction_force(dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
    end subroutine

    subroutine UpdateElmtInfo_(this,dh,zDim,yDim,xDim,xGrid,yGrid,zGrid) ! update element position, velocity, weight
        implicit none
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: dh
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: xGrid(xDim),yGrid(yDim),zGrid(zDim)
        call this%UpdateElmtPosVelArea()
        call this%UpdateElmtInterp(dh,zDim,yDim,xDim,xGrid,yGrid,zGrid)
    end subroutine UpdateElmtInfo_

    subroutine UpdateElmtPosVelArea_(this)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        integer :: i,iEL,n_theta
        integer :: Ea_A,Ea_B,Ea_C,Ea_D
        real(8) :: E_left,E_right
        real(8) :: pi
        real(8) :: dxyz0(3),dxyz(3),xlmn0(3),xlmn(3),dl0,dl,temp_xyz(3),temp_xyzoffset(3)
        real(8) :: r_rotMat(3,3),r_self_rotMat(3,3)
        real(8) :: self_rot_omega(this%r_nelmts,3),temp_vel(3)
        integer :: s,begin,end
        real(8) :: xyz1(3),xyz2(3)
        real(8) :: dspan,invl0,cos_xyz(3)
        pi=4.0d0*datan(1.0d0)
        !   compute displacement, velocity, area at surface element center
        if (this%v_type .eq. 0) then
            this%v_Evel = 0.0d0
        elseif ((this%v_type .eq. 1)) then
            invl0 = 1 / dsqrt(sum(this%v_dirc(1:3)**2))
            cos_xyz = this%v_dirc(1:3) * invl0
            do i = 1, this%r_nelmts
                dspan = this%r_Lspan(i)/this%r_Nspan(i)
                xyz1(1:3) = this%r_xyz(this%r_elmt(i,1),1:3)
                xyz2(1:3) = this%r_xyz(this%r_elmt(i,2),1:3)
                do s = 1,this%r_Nspan(i)
                    this%v_Exyz(i,1:3) = (xyz1(1:3)+xyz2(1:3))*0.5d0 + dspan*(s-0.5)*cos_xyz(1:3)
                    call UpdateElmtArea_()
                enddo
            enddo
        elseif ((this%v_type .eq. 2)) then
            do i = 1,this%r_nelmts
                ! update r_rotMat and self_rotMar
                dxyz0(1:3)= this%r_xyz0(this%r_elmt(i,1),1:3)-this%r_xyz0(this%r_elmt(i,2),1:3)
                dxyz(1:3) = this%r_xyz0(this%r_elmt(i,1),1:3)-this%r_xyz0(this%r_elmt(i,2),1:3)
                dl0       = dsqrt(sum(dxyz0(1:3)**2))
                dl        = dsqrt(sum(dxyz(1:3)**2))
                xlmn0(1:3)= dxyz0(1:3)/dl0
                xlmn(1:3) = dxyz(1:3)/dl
                call RotateMatrix(xlmn0,xlmn,r_rotMat)
                call Self_RotateMatrix(this%r_xyz(i,4),xlmn,r_self_rotMat)
                n_theta  = maxval(this%r_Nspan(:))
                self_rot_omega(i,1:3) = xlmn(1:3)*(this%r_vel(this%r_elmt(i,1),4)+this%r_vel(this%r_elmt(i,2),4))*0.5d0

                if (i .eq. this%r_nelmts) then
                    begin = this%rtov(i)
                    end = this%v_nelmts
                else
                    begin = this%rtov(i)
                    end = this%rtov(i+1)-1
                endif

                do iEL = begin,end
                    ! update xyz
                    temp_xyz(1:3) = matmul(r_rotMat(:,:), (/this%v_Exyz0(iEL,1), 0.0d0, this%v_Exyz0(iEL,3)/))
                    temp_xyz(1:3) = matmul(r_self_rotMat(:,:), temp_xyz)
                    ! Default The fake point on the plane (this%v_Exyz0(iEL,1:3) = (x,y,z)) is in the same plane as the point on the centre axis(this%r_xyz0(i,1:3) = (0,y,0)).
                    temp_xyzoffset(1:3) = (this%r_xyz(this%r_elmt(i,1),1:3)+this%r_xyz(this%r_elmt(i,2),1:3))*0.5d0 - &
                                        (this%r_xyz0(this%r_elmt(i,1),1:3)+this%r_xyz0(this%r_elmt(i,2),1:3))*0.5d0
                    this%v_Exyz(iEL,1:3) = temp_xyz(1:3)+temp_xyzoffset(1:3)
                    this%v_Exyz(iEL,2) = this%v_Exyz(iEL,2) + this%v_Exyz0(iEL,2)
                    ! update vel
                    call cross_product(self_rot_omega(i,:),temp_xyz,temp_vel)
                    this%v_Evel(iEL,1:3) = temp_vel(1:3)+(this%r_vel(this%r_elmt(i,1),4)+this%r_vel(this%r_elmt(i,2),4))*0.5d0
                    enddo
            enddo
        endif
        contains
        subroutine UpdateElmtArea_()
            IMPLICIT NONE
            integer :: r
            if (this%v_type .eq. 0) then
            elseif ((this%v_type .eq. 1)) then
                ! only rectangle
                this%v_Ea(i) = dspan * dsqrt(((xyz1(1)-xyz2(1))**2)+((xyz1(2)-xyz2(2))**2)+((xyz1(3)-xyz2(3))**2))
            elseif ((this%v_type .eq. 2)) then
                n_theta  = maxval(this%r_Nspan(:))
                do iEL = 2,this%v_nelmts-1
                    r = mod((iEL-1),n_theta)
                    !update area
                    Ea_A = iEL-1
                    if (r.eq.1) Ea_A = iEL+(n_theta-1)
                    Ea_B = iEL-n_theta
                    if (Ea_B .le. 1) Ea_B = 1
                    Ea_C = iEL+1
                    if (r.eq.0) Ea_C = iEL-(n_theta-1)
                    Ea_D = iEL+n_theta
                    if (Ea_D .ge. this%v_nelmts) Ea_D = this%v_nelmts

                    call cpt_area(this%v_Exyz(Ea_A,1:3),this%v_Exyz(Ea_B,1:3),this%v_Exyz(Ea_D,1:3),E_left)
                    call cpt_area(this%v_Exyz(Ea_B,1:3),this%v_Exyz(Ea_C,1:3),this%v_Exyz(Ea_D,1:3),E_right)
                    this%v_Ea(iEL) = (E_left+E_right)*0.5d0
                enddo
                this%v_Ea(1) = 0.d0
                this%v_Ea(this%v_nelmts) = 0.d0
            endif
        end subroutine UpdateElmtArea_
    end subroutine UpdateElmtPosVelArea_
    
    FUNCTION vtor_(this,x)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer(2) :: vtor_(2),x
        integer :: i
        real(8) :: y1,y2
        if (this%v_type .eq. 1) then
            do i =1,this%r_nelmts
                if ((x-sum(this%r_Nspan(1:i))) .le. 0) then
                    vtor_(1) = i
                    vtor_(2) = i+1
                endif
            enddo
        elseif (this%v_type .eq. 2) then
            if (dabs(this%v_Exyz0(x,2) - this%r_xyz0(1,2)) .lt. 1e-5) then
                vtor_(1:2) = 1
            endif
            if (dabs(this%v_Exyz0(x,2) - this%r_xyz0(this%r_npts,2)) .lt. 1e-5) then
                vtor_(1:2) = this%r_npts
            endif
            do i =1,this%r_nelmts
                y1 = this%r_xyz0(this%r_elmt(i,1),2)
                y2 = this%r_xyz0(this%r_elmt(i,2),2)
                if ((this%v_Exyz0(x,2) .gt. y1) .and. (this%v_Exyz0(x,2) .lt. y2)) then
                    vtor_(1) = i
                    vtor_(2) = i+1
                endif
            enddo
        endif
    ENDFUNCTION vtor_
    FUNCTION rtov_(this,x)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer(2) :: rtov_(2),x
        real(8) :: pi,dspan
        integer :: n_theta,n_radius
        pi=4.0d0*datan(1.0d0)
        if (this%v_type .eq. 1) then
            rtov_(1) = sum(this%r_Nspan(1:(x-1)))+1
            rtov_(2) = sum(this%r_Nspan(1:(x)))
        elseif (this%v_type .eq. 2) then
            n_theta  = maxval(this%r_Nspan(:))
            if (n_theta  .le. 3) n_theta = 3
            dspan = (2*pi*this%r_Lspan(maxloc(this%r_Nspan(:),1))/n_theta)
            n_radius = floor(this%r_Lspan(maxloc(this%r_Nspan(:),1))/dspan)
            if (n_radius .eq. 1) n_radius = 2
            if (x .eq. 1) then
                rtov_(1) = 1
                rtov_(2) = n_radius*n_theta+1+n_theta
            elseif (x .eq. this%r_nelmts) then
                rtov_(1) = this%v_nelmts-n_radius*n_theta-n_theta
                rtov_(2) = this%v_nelmts
            else
                rtov_(1) = n_radius*n_theta+1+n_theta+n_theta*(x-2)+1
                rtov_(2) = n_radius*n_theta+1+n_theta+n_theta*(x-1)+1
            endif
        endif
    ENDFUNCTION rtov_

    subroutine UpdateElmtInterp_(this,dh,zDim,yDim,xDim,xGrid,yGrid,zGrid)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: dh
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: xGrid(xDim),yGrid(yDim),zGrid(zDim)
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2),Phi
        real(8)::x0,y0,z0,detx,dety,detz,invdh
        integer::i0,j0,k0,iEL,x,y,z,i,j,k
        !==================================================================================================
        invdh = 1.D0/dh
        call my_minloc(this%v_Exyz(1,1), xGrid, xDim, .false., i0)
        call my_minloc(this%v_Exyz(1,2), yGrid, yDim, .false., j0)
        call my_minloc(this%v_Exyz(1,3), zGrid, zDim, .false., k0)
        x0 = xGrid(i0)
        y0 = yGrid(j0)
        z0 = zGrid(k0)
        ! compute the velocity of IB nodes at element center
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,k,x,y,z,rx,ry,rz,detx,dety,detz,ix,jy,kz)
        do  iEL=1,this%v_nelmts
            call minloc_fast(this%v_Exyz(iEL,1), x0, i0, invdh, i, detx)
            call minloc_fast(this%v_Exyz(iEL,2), y0, j0, invdh, j, dety)
            call minloc_fast(this%v_Exyz(iEL,3), z0, k0, invdh, k, detz)
            call trimedindex(i, xDim, ix, m_boundaryConditions(1:2))
            call trimedindex(j, yDim, jy, m_boundaryConditions(3:4))
            call trimedindex(k, zDim, kz, m_boundaryConditions(5:6))
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

    SUBROUTINE calculate_interaction_force(dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
        ! calculate elements interaction force using IB method
        IMPLICIT NONE
        real(8),intent(in):: dt,dh,denIn,Uref
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: xGrid(xDim),yGrid(yDim),zGrid(zDim),den(zDim,yDim,xDim)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        !================================
        integer:: iFish
        integer:: i,j,k,iEL,nt,iterLBM
        real(8):: dmaxLBM,dsum
        real(8)::tol,tolsum, ntol
        !================================
        tolsum=0.d0
        iterLBM=0
        do iFish=1,m_nFish
            allocate(VBodies(iFish)%r_force(VBodies(iFish)%r_npts,6))
            VBodies(iFish)%r_force(:,:) = 0.0d0
        enddo
        do  while( iterLBM<m_maxIterIB .and. dmaxLBM>m_dtolLBM)
            dmaxLBM=0.0d0
            dsum=0.0d0
            do iFish=1,m_nFish
                call VBodies(iFish)%PenaltyForce(dh,dt,denIn,zDim,yDim,xDim,uuu,den,force,tol,ntol)
                tolsum = tolsum + dsqrt(tol)
                dsum = dsum + Uref*ntol
            enddo
            dmaxLBM=tolsum/dsum
            iterLBM=iterLBM+1
        enddo
    END SUBROUTINE

    SUBROUTINE PenaltyForce_(this,dh,dt,denIn,zDim,yDim,xDim,uuu,den,force,tolerance,ntolsum)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: dh,dt,denIn
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: den(zDim,yDim,xDim)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
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
            forceElemTemp(1:3) = -m_Pbeta* 2.0d0*denIn*(velElem-velElemIB)/dt*this%v_Ea(iEL)*dh
            if ((.not. IEEE_IS_FINITE(forceElemTemp(1))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(2))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(3)))) then
                write(*, *) 'Nan found in forceElemTemp', forceElemTemp
                write(*, *) 'Nan found at (ix, jy, kz)', ix(0), jy(0), kz(0)
                stop
            endif
            tolerance = tolerance + sum((velElem(1:3)-velElemIB(1:3))**2)
            !$OMP critical
            ! update beam load, momentum is not included
            this%r_force(this%vtor(iEL,1),1:3) = this%r_force(this%vtor(iEL,1),1:3) + 0.5d0 * forceElemTemp
            this%r_force(this%vtor(iEL,2),1:3) = this%r_force(this%vtor(iEL,2),1:3) + 0.5d0 * forceElemTemp
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        forceTemp(1:3) = -forceElemTemp(1:3)*rx(x)*ry(y)*rz(z)
                        ! correct velocity
                        uuu(kz(z),jy(y),ix(x),1:3)  = uuu(kz(z),jy(y),ix(x),1:3)+0.5d0*dt*forceTemp(1:3)/den(kz(z),jy(y),ix(x))
                        ! add flow body force
                        force(kz(z),jy(y),ix(x),1:3) = force(kz(z),jy(y),ix(x),1:3) + forceTemp(1:3)
                    enddo
                enddo
            enddo
            !$OMP end critical
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE PenaltyForce_

    subroutine BuildPlate_(this, beam)
        implicit none
        class(VirtualBody), intent(inout) :: this
        type(BeamSolver), intent(in):: beam
        
        
        integer(2),allocatable :: vtor(:)!of size fake_npts
        integer :: r_npts,r_nelmts
        integer(2),allocatable :: rtov(:,:)! of size real_npts, 1:2


        call this%Readreal()
        call this%Struc_Virtual_Elmts()
    end subroutine BuildPlate_

    subroutine Beam_ReadUnstructured_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        call this%Unstruc_Virtual_Elmts()
        this%v_type = 0
        this%r_npts = 1
        allocate(this%r_xyz0(1,3))
        this%r_xyz0(1,1:3) = this%v_Pxyz0(1,1:3)
        this%r_Lspan = 0.0d0
        this%r_Nspan = 0



        allocate(this%v_Exyz(1:this%v_nelmts,1:6))
        this%v_Exyz = 0.0d0
        allocate(this%v_Evel(1:this%v_nelmts,1:6))
        this%v_Evel = 0.0d0
        allocate(this%v_Ei(1:this%v_nelmts,1:12))
        this%v_Ei = 0.0d0
        allocate(this%v_Ew(1:this%v_nelmts,1:12))
        this%v_Ew = 0.0d0
        allocate(this%v_Ea(1:this%v_nelmts))
        this%v_Ea = 0.0d0
    end subroutine Beam_ReadUnstructured_

    subroutine Section_RotateMatrix(lmn0,lmn,r_rotMat)
        implicit none
        ! real(8):: lmn0(3),lmn(3),angle,pi
        ! real(8):: v_0(3),e1(3),e2(3),e3(3),r1(3,3)
        ! real(8):: r_rotMat(3,3)
        ! pi=4.0*datan(1.0d0)
        ! v_0= lmn0
        ! e1 = lmn
        ! e2 = 0.0d0
        ! e3 = 0.0d0
        ! call angle_between_vectors(v_0,e1)
        ! if (angle .gt. 1e-3) then
        !     call normal_vector(v_0,e1,e2)
        !     call normal_vector(e1,e2,e3)
        ! else
        !     e1 = (/1.0d0,0.0d0,0.0d0/)
        !     e2 = (/0.0d0,1.0d0,0.0d0/)
        !     e3 = (/0.0d0,0.0d0,1.0d0/)
        ! endif
        ! r1(1,1)  =   e1(1)
        ! r1(1,2)  =   e1(2)
        ! r1(1,3)  =   e1(3)
        ! r1(2,1)  =   e2(1)
        ! r1(2,2)  =   e2(2)
        ! r1(2,3)  =   e2(3)
        ! r1(3,1)  =   e3(1)
        ! r1(3,2)  =   e3(2)
        ! r1(3,3)  =   e3(3)
        ! r_rotMat(:,:) = r1
        ! contains
        real(8), intent(in):: lmn0(3),lmn(3)
        real(8), intent(out):: r_rotMat(3,3)
        real(8):: angle,pi
        real(8):: v_1(3),v_2(3),nn(3),r1(3,3)
        pi=4.0*datan(1.0d0)
        v_1=lmn0
        v_2=lmn
        nn=0.0d0
        angle=0.0d0
        call angle_between_vectors(v_1,v_2)
        call normal_vector(v_1,v_2,nn)
        call quaternion_rotate(angle,nn,r1)
        r_rotMat(:,:) = r1(:,:)
        contains
        subroutine normal_vector(A,B,C)
            implicit none
            real(8):: A(3),B(3),C(3)
            real(8) :: norm_C
            call cross_product(A,B,C)
            norm_C = sqrt(C(1)**2 + C(2)**2 + C(3)**2)
            if (abs(norm_C - 0.0d0).gt.1e-10) then
                C = C / norm_C
            endif
            return
        end subroutine normal_vector
        subroutine angle_between_vectors(A,B)
            implicit none
            real(8) :: A(3),B(3)
            real(8) :: dot_product, magnitude1, magnitude2
            dot_product = A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
            magnitude1 = sqrt(A(1)**2+A(2)**2+A(3)**2)
            magnitude2 = sqrt(B(1)**2+B(2)**2+B(3)**2)
            angle = acos(dot_product / (magnitude1 * magnitude2))
            return
        end subroutine angle_between_vectors
    end subroutine Section_RotateMatrix
    subroutine Section_Self_RotateMatrix(angle,lmn,r_self_rotMat)
        implicit none
        real(8), intent(in):: angle,lmn(3)
        real(8), intent(out):: r_self_rotMat(3,3)
        real(8):: r1(3,3)
        call quaternion_rotate(angle,lmn,r1)
        r_self_rotMat(:,:) = r1(:,:)
    end subroutine Section_Self_RotateMatrix
    subroutine quaternion_rotate(angle,nn,r1)
        implicit none
        real(8), intent(in) :: angle,nn(3)
        real(8), intent(out) :: r1(3,3)
        real(8) :: lamda(0:3)
        lamda(0)=cos(angle*0.5)
        lamda(1)=nn(1)*sin(angle*0.5)
        lamda(2)=nn(2)*sin(angle*0.5)
        lamda(3)=nn(3)*sin(angle*0.5)
        r1(1,1) = 2*(lamda(0)**2+lamda(1)**2)-1
        r1(1,2) = 2*(lamda(1)*lamda(2)-lamda(0)*lamda(3))
        r1(1,3) = 2*(lamda(1)*lamda(3)+lamda(0)*lamda(2))
        r1(2,1) = 2*(lamda(2)*lamda(1)+lamda(0)*lamda(3))
        r1(2,2) = 2*(lamda(0)**2+lamda(2)**2)-1
        r1(2,3) = 2*(lamda(2)*lamda(3)-lamda(0)*lamda(1))
        r1(3,1) = 2*(lamda(3)*lamda(1)-lamda(0)*lamda(2))
        r1(3,2) = 2*(lamda(3)*lamda(2)+lamda(0)*lamda(1))
        r1(3,3) = 2*(lamda(0)**2+lamda(3)**2)-1
    end subroutine quaternion_rotate
    subroutine cross_product(w,r,v)
        implicit none
        real(8), intent(in) :: w(3), r(3)
        real(8), intent(out) :: v(3)
        v(1) = w(2)*r(3) - w(3)*r(2)
        v(2) = w(3)*r(1) - w(1)*r(3)
        v(3) = w(1)*r(2) - w(2)*r(1)
    end subroutine cross_product

    subroutine Unstruc_Virtual_Elmts_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer :: fileiD = 111, num, i, temp_nelmts, temp_prop(4)
        open(unit=fileiD, file = this%filename )! read *.msh file
            ! read nodes
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*) num
            this%v_npts = num
            allocate(this%v_Pxyz0(1:this%v_npts, 1:6))
            this%v_Pxyz0 = 0.0d0
            do i = 1,this%v_npts
                read(fileiD,*)num,this%v_Pxyz0(i,1),this%v_Pxyz0(i,2),this%v_Pxyz0(i,3)
            enddo
            ! read element
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*) num
            temp_nelmts = num
            do i = 1,temp_nelmts
                read(fileiD,*)num,temp_prop(1:4)
                if((temp_prop(1)==2) .and. (temp_prop(2)==2) .and. (temp_prop(3)==0)) then
                    this%v_nelmts = temp_nelmts-num+1
                    exit
                endif
            enddo

            backspace(fileiD)

            allocate(this%v_elmt(1:this%v_nelmts, 1:3))
            this%v_elmt = 0
            allocate(this%v_Ea(1:this%v_nelmts))
            this%v_Ea = 0.0d0
            allocate(this%v_Exyz0(1:this%v_nelmts, 1:3))
            this%v_Exyz0 = 0.0d0
            do i = 1,this%v_nelmts
                read(fileiD,*)num,temp_prop(1:4),this%v_elmt(i,1),this%v_elmt(i,2),this%v_elmt(i,3)
                call cpt_area(this%v_Pxyz0(this%v_elmt(i,1),1:3),this%v_Pxyz0(this%v_elmt(i,2),1:3),this%v_Pxyz0(this%v_elmt(i,3),1:3),this%v_Ea(i))
                call cpt_incenter(this%v_Pxyz0(this%v_elmt(i,1),1:3),this%v_Pxyz0(this%v_elmt(i,2),1:3),this%v_Pxyz0(this%v_elmt(i,3),1:3),this%v_Exyz0(i,1:3))
            enddo

        close(fileiD)
    end subroutine Unstruc_Virtual_Elmts_
    subroutine cpt_incenter(NodeA,NodeB,NodeC,Exyz0)
        implicit none
        real(8), intent(in) :: NodeA(3),NodeB(3),NodeC(3)
        real(8), intent(out) :: Exyz0(3)
        real(8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3
        real(8) :: a, b, c, invC
            x1=NodeA(1)
            x2=NodeB(1)
            x3=NodeC(1)
            y1=NodeA(2)
            y2=NodeB(2)
            y3=NodeC(2)
            z1=NodeA(3)
            z2=NodeB(3)
            z3=NodeC(3)
            a = sqrt((x2 - x3)**2 + (y2 - y3)**2 + (z2 - z3)**2)
            b = sqrt((x3 - x1)**2 + (y3 - y1)**2 + (z3 - z1)**2)
            c = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
            invC = 1 / (a + b + c)
            Exyz0(1) = (a*x1 + b*x2 + c*x3) * invC
            Exyz0(2) = (a*y1 + b*y2 + c*y3) * invC
            Exyz0(3) = (a*z1 + b*z2 + c*z3) * invC
    end subroutine cpt_incenter
    subroutine cpt_area(A,B,C,area)
        implicit none
        real(8), intent(in) :: A(3),B(3),C(3)
        real(8), intent(out) :: area
        real(8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3, ax, ay, az
            x1=A(1)
            x2=B(1)
            x3=C(1)
            y1=A(2)
            y2=B(2)
            y3=C(2)
            z1=A(3)
            z2=B(3)
            z3=B(3)
            ax =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.0d0
            ay =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.0d0
            az =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.0d0
            area=dsqrt( ax*ax + ay*ay + az*az)
    end subroutine cpt_area
    subroutine Readreal_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer:: i, node, fileiD = 111, temp
        open(unit=fileiD, file = this%filename)
        rewind(fileiD)
        read(fileiD,*)
        read(fileiD,*) this%r_npts,this%r_nelmts,this%v_type,this%v_dirc(1:3)
        allocate(this%r_xyz0(this%r_npts,3))
        allocate(this%r_elmt(this%r_nelmts,5))
        allocate(this%r_Lspan(this%r_nelmts),this%r_Nspan(this%r_nelmts))
        read(fileiD,*)
        read(fileiD,*)
        do    i= 1, this%r_npts
            read(fileiD,*) node,this%r_xyz0(node,1),this%r_xyz0(node,2),this%r_xyz0(node,3)
        enddo
        read(fileiD,*)
        read(fileiD,*)
        do  i= 1, this%r_nelmts
            read(fileiD,*) node,this%r_elmt(node,1:5),this%r_Lspan(node),this%r_Nspan(node)
        enddo
        close(fileiD)
    end subroutine Readreal_
    subroutine Struc_Virtual_Elmts_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        if (this%v_type .eq. 1) call this%Struc_Rib_Elmts()
        if (this%v_type .eq. 2) call this%Struc_Cylinder_Elmts()
    end subroutine Struc_Virtual_Elmts_
    subroutine Struc_Rib_Elmts_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer(2) :: i
        integer :: s,num
        real(8) :: xyz1(3),xyz2(3)
        real(8) :: dspan,invl0,cos_xyz(3)
        this%v_nelmts = sum(this%r_Nspan)
        allocate(this%v_Exyz0(1:this%v_nelmts, 1:3))
        this%v_Exyz0 = 0.0d0
        allocate(this%vtor(1:this%v_nelmts,1:2))
        this%vtor = 0
        invl0 = 1 / dsqrt(sum(this%v_dirc(1:3)**2))
        cos_xyz = this%v_dirc(1:3) * invl0
        num = 1
        do i = 1, this%r_nelmts
            dspan = this%r_Lspan(i)/this%r_Nspan(i)
            xyz1(1:3) = this%r_xyz0(this%r_elmt(i,1),1:3)
            xyz2(1:3) = this%r_xyz0(this%r_elmt(i,2),1:3)
            do s = 1,this%r_Nspan(i)
                this%v_Exyz0(i,1:3) = (xyz1(1:3)+xyz2(1:3))*0.5d0 + dspan*(s-0.5)*cos_xyz(1:3)
                this%vtor(num,1:2) = this%r_elmt(i,1:2)
                num = num + 1
            enddo
        enddo
    end subroutine Struc_Rib_Elmts_
    subroutine Struc_Cylinder_Elmts_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        real(8) :: pi
        real(8) :: dspan, theta, rho, height, dradius,dtheta,y1,y2
        integer :: n_height,n_theta,n_radius
        integer(2) :: i
        integer :: j,k,num
        integer :: n_circle, n_rectangle
        pi=4.0d0*datan(1.0d0)

        n_height = this%r_nelmts
        n_theta  = maxval(this%r_Nspan(:))
        if (n_theta  .le. 3) n_theta = 3
        dspan = (2*pi*this%r_Lspan(maxloc(this%r_Nspan(:),1))/n_theta)
        n_radius = floor(this%r_Lspan(maxloc(this%r_Nspan(:),1))/dspan)
        if (n_radius .eq. 1) n_radius = 2
        dtheta = 2.0 * pi / n_theta

        n_circle = (n_radius * n_theta) + 1
        n_rectangle = n_height * n_theta
        this%v_nelmts = n_circle * 2 + n_rectangle

        allocate(this%v_Exyz0(1:this%v_nelmts, 1:3))
        this%v_Exyz0 = 0.0d0
        allocate(this%vtor(1:this%v_nelmts,1:2))
        this%vtor = 0
        allocate(this%rtov(1:this%r_nelmts))
        this%rtov = 0

        ! Given node v_Exyz0 coordinate value
        ! lower circle
        num = 1
        this%v_Exyz0(num,1:3) = this%r_xyz0(1,1:3)
        this%vtor(num,1:2) = 1
        num = num + 1
        dradius = this%r_Lspan(1)/n_radius
        height = this%r_xyz0(1,2)
        do j = 1, n_radius
            rho = (j - 0.5) * dradius
            do k = 1, n_theta
                theta = (k - 0.5) * dtheta
                this%v_Exyz0(num,1) = rho * cos(theta)
                this%v_Exyz0(num,2) = height
                this%v_Exyz0(num,3) = rho * sin(theta)
                this%vtor(num,1:2) = 1
                num = num + 1
            enddo
        enddo
        ! side
        do i = 1, this%r_nelmts
            if (i.eq.1) then
                this%rtov(1) = 1
            else
                this%rtov(i) = num
            endif
            rho = this%r_Lspan(i)
            y1 = this%r_xyz0(this%r_elmt(i,1),2)
            y2 = this%r_xyz0(this%r_elmt(i,2),2)
            height = (y1+y2)*0.5d0
            do k = 1, n_theta
                theta = (k - 0.5) * dtheta
                this%v_Exyz0(num,1) = rho * cos(theta)
                this%v_Exyz0(num,2) = height
                this%v_Exyz0(num,3) = rho * sin(theta)
                this%vtor(num,1) = this%r_elmt(i,1)
                this%vtor(num,2) = this%r_elmt(i,2)
                num = num+1
            enddo
        enddo
        ! upper circle
        dradius = this%r_Lspan(this%r_nelmts)/n_radius
        height = this%r_xyz0(this%r_npts,2)
        do j = 1, n_radius
            rho = (n_radius - j + 0.5) * dradius
            do k = 1, n_theta
                theta = (k - 0.5) * dtheta
                this%v_Exyz0(num,1) = rho * cos(theta)
                this%v_Exyz0(num,2) = height
                this%v_Exyz0(num,3) = rho * sin(theta)
                this%vtor(num,1:2) = this%r_npts
                num = num + 1
            enddo
        enddo
        this%v_Exyz0(num,1:3) = this%r_xyz0(i,1:3)
        this%vtor(num,1:2) = this%r_npts
    end subroutine Struc_Cylinder_Elmts_

    subroutine Write_body_(this,iFish,time,Lref,Tref)
        ! to do: generate a temporary mesh
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer,intent(in) :: iFish
        real(8),intent(in) :: time,Lref,Tref
        !   -------------------------------------------------------
        real(8):: pi,timeTref
        integer:: i,j,r,Nspanpts,ElmType,n_theta,Ea_A,Ea_D
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName,idstr
        integer,parameter:: idfile=100
        real(8) :: invl0,cos_xyz(3),dspan,xyz1(3),xyz2(3),low_L(3),low_R(3),upp_R(3),upp_L(3)
        pi=4.0d0*datan(1.0d0)
        !==========================================================================
        timeTref = time/Tref
        !==========================================================================
        write(fileName,'(I10)') nint(timeTref*1d5)
        fileName = adjustr(fileName)
        do  I=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo
    
        ! ElmType = 3

        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        open(idfile, FILE='./DatBodySpan/BodyFake'//trim(idstr)//'_'//trim(filename)//'.dat')
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        if (this%v_type .eq. 0) then
            write(idfile, '(A,I7,A,I7,A)') 'ZONE N=',this%v_npts,', E=',this%v_nelmts,', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
        elseif (this%v_type .eq. 1) then
            write(idfile, '(A,I7,A,I7,A)') 'ZONE N=',4*this%v_nelmts,', E=',this%v_nelmts,', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
        elseif (this%v_type .eq. 2) then
            write(idfile, '(A,I7,A,I7,A)') 'ZONE N=',this%v_nelmts,', E=',this%v_nelmts-2,', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
        endif
        ! if(ElmType.eq.2) then
        !     write(idfile, '(A)') 'FELINESEG'
        ! elseif (ElmType.eq.3) then
        !    write(idfile, '(A)') 'FETRIANGLE'
        ! elseif(ElmType.eq.4) then
        !     write(idfile, '(A)') 'FEQUADRILATERAL'
        ! endif
        if (this%v_type .eq. 0) then
            do  i=1,this%v_npts
                write(idfile, *)  this%v_Pxyz0(i,1:3)!/Lref
            enddo
            do  i=1,this%v_nelmts
                ! if(ElmType.eq.2) then
                !     write(idfile, *) this%v_elmt(i,1),this%v_elmt(i,2)
                ! elseif(ElmType.eq.3) then
                    write(idfile, *) this%v_elmt(i,1),this%v_elmt(i,2),this%v_elmt(i,3)
                ! elseif(ElmType.eq.4) then
                !     write(idfile, *) this%v_elmt(i,1),this%v_elmt(i,2),this%v_elmt(i,3),this%v_elmt(i,4)
                ! else
                ! endif
            enddo
        elseif (this%v_type .eq. 1)then
            invl0 = 1 / dsqrt(sum(this%v_dirc(1:3)**2))
            cos_xyz(1:3) = this%v_dirc(1:3) * invl0
            do i = 1, this%r_nelmts
                dspan = this%r_Lspan(i)/this%r_Nspan(i)
                Nspanpts = this%r_Nspan(i)+1
                xyz1(1:3) = this%r_xyz(this%r_elmt(i,1),1:3)
                xyz2(1:3) = this%r_xyz(this%r_elmt(i,2),1:3)
                do j = 1,Nspanpts
                    low_L(1:3) = xyz1(1:3) + dspan*(j-1)*cos_xyz(1:3)
                    low_R(1:3) = xyz1(1:3) + dspan*(j)*cos_xyz(1:3)
                    upp_R(1:3) = xyz2(1:3) + dspan*(j)*cos_xyz(1:3)
                    upp_L(1:3) = xyz2(1:3) + dspan*(j-1)*cos_xyz(1:3)
                    write(idfile, *) low_L(1:3)!/Lref
                    write(idfile, *) low_R(1:3)!/Lref
                    write(idfile, *) upp_R(1:3)!/Lref
                    write(idfile, *) upp_L(1:3)!/Lref
                enddo
            enddo
            do  i=1,this%v_nelmts
                j = 4*(i-1)
                write(idfile, *) j+1,j+2,j+3,j+4
            enddo
        elseif (this%v_type .eq. 2)then
            do  i=1,this%v_nelmts
                write(idfile, *)  this%v_Exyz(i,1:3)!/Lref
            enddo
            n_theta  = maxval(this%r_Nspan(:))
            do  i=2,this%v_nelmts-1
                r = mod((i-1),n_theta)
                Ea_A = i-1
                if (r.eq.1) Ea_A = i+(n_theta-1)
                Ea_D = i+n_theta
                if (Ea_D .ge. this%v_nelmts) Ea_D = this%v_nelmts

                write(idfile, *) Ea_A,i,Ea_D
            enddo
        endif
        close(idfile)
        !   =============================================
    end subroutine

    subroutine Write_solid_bodies(time,Lref,Tref)
        implicit none
        real(8) :: Lref,time,Tref
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%Write_body(iFish,time,Lref,Tref)
        enddo
    end subroutine Write_solid_bodies

end module SolidBody