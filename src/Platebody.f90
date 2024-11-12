module SolidBody
    implicit none
    private
    ! Immersed boundary method parameters
    integer:: m_nFish, m_maxIterIB
    real(8):: m_dtolLBM, m_Pbeta
    integer:: boundaryConditions(1:6)
    ! nFish     number of bodies
    ! maxIterIB maximum number of iterations for IB force calculation
    ! dtolIBM   tolerance for IB force calculation
    ! Pbeta     coefficient in penalty force calculation
    public :: BeamBody,Initialise_bodies,Write_solid_bodies,FSInteraction_force
    type :: BeamBody
        character(LEN=100) :: filename
        !!!centeral line beam
        integer :: r_npts
        real(8), allocatable :: r_xyz0(:, :)
        real(8), allocatable :: r_xyz(:, :)
        real(8), allocatable :: r_vel(:, :)
        real(8), allocatable :: r_force(:, :)
        !!!virtual body surface
        integer :: v_npts,v_nelmts
        real(8), allocatable :: v_Pxyz0(:, :) ! initial Node (x0, y0, z0)
        real(8), allocatable :: v_Pxyz(:, :) ! Node (xt, yt, zt)
        real(8), allocatable :: v_Exyz0(:, :) ! initial element center (x, y, z)
        real(8), allocatable :: v_Exyz(:, :) ! element center (x, y, z)
        integer, allocatable :: v_elmt(:, :) ! element ID
        integer, allocatable :: v_Ei(:, :) ! element stencial integer index [ix-1,ix,ix1,ix2, iy-1,iy,iy1,iy2, iz-1,iz,iz1,iz2]
        real(8), allocatable :: v_Ew(:, :) ! element stential weight [wx-1, wx, wx1, wx2, wy-1, wy, wy1, wy2, wz-1, wz, wz1, wz2]
        real(8), allocatable :: v_Ea(:) ! element area
        !area center with equal weight on both sides
        real(8), allocatable :: v_Evel(:, :)
        !calculated using central linear and angular velocities
        integer(2),allocatable :: vtor(:)
        real(8), allocatable :: rotMat(:, :, :)
        real(8), allocatable :: self_rotMat(:, :, :)
    contains
        procedure :: Initialise => Initialise_
        procedure :: UpdateElmtInfo => UpdateElmtInfo_
        procedure :: PenaltyForce => PenaltyForce_
        procedure :: BuildStructured => Beam_BuildStructured
        procedure :: ReadUnstructured => Beam_ReadUnstructured

        procedure :: Readreal => Readreal_
        procedure :: Struc_Virtual_Elmts => Struc_Virtual_Elmts_
        procedure :: Struc_Cylinder_Elmts => Struc_Cylinder_Elmts_

        procedure :: Unstruc_Virtual_Elmts => Unstruc_Virtual_Elmts_


        procedure :: RotateMatrix => Section_RotateMatrix
        procedure :: Self_RotateMatrix => Section_Self_RotateMatrix

        procedure :: Write_body => Write_body_
    end type BeamBody
    type(BeamBody), allocatable :: Beam(:)
  contains
    subroutine Initialise_(this,filename,iBodyModel)
        ! read beam central line file and allocate memory
        implicit none
        class(BeamBody), intent(inout) :: this
        character (LEN=100), intent(in):: filename
        integer :: iBodyModel

        this%filename = filename

        if (iBodyModel .eq. 2) then
            call this%BuildStructured()
        elseif (iBodyModel .eq. 1) then
            call this%ReadUnstructured()
        endif

        allocate(this%rotMat(1:this%r_npts,1:3,1:3))
        this%rotMat=0.0d0
        allocate(this%self_rotMat(1:this%r_npts,1:3,1:3))
        this%self_rotMat=0.0d0
        allocate(this%v_Exyz(1:this%v_nelmts,1:6))
        allocate(this%v_Evel(1:this%v_nelmts,1:6))
        call this%UpdateElmtInfo()
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
            call Beam(iFish)%UpdateElmtInfo()
        enddo
        call calculate_interaction_force(dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
    end subroutine

    subroutine UpdateElmtInfo_(this) ! update element position, velocity, weight
        implicit none
        class(BeamBody), intent(inout) :: this
        integer :: i,j,k,iEL,n_theta
        integer :: Ea_A,Ea_B,Ea_C,Ea_D
        real(8) :: E_left,E_right
        real(8) :: pi
        real(8) :: dxyz0(3),dxyz(3),xlmn0(3),xlmn(3),dl0,dl,temp_xyz(3),temp_xyzoffset(3)
        real(8) :: self_rot_omega(3),temp_vel(3)
        pi=4.0d0*datan(1.0d0)
        !   compute displacement, velocity, area at surface element center
        if (this%r_npts .eq. 1) then
            this%v_Evel = 0.0d0
        else
            do i = 1,this%r_npts

                ! update rotMat and self_rotMar
                if ( i .eq. 1) then
                    dxyz0(1:3)= this%r_xyz0(i+1,1:3)-this%r_xyz0(i,1:3)
                    dxyz(1:3) = this%r_xyz(i+1,1:3)-this%r_xyz(i,1:3)
                elseif ( i .eq. this%r_npts) then
                    dxyz0(1:3)= this%r_xyz0(i,1:3)-this%r_xyz0(i-1,1:3)
                    dxyz(1:3) = this%r_xyz(i,1:3)-this%r_xyz(i-1,1:3)
                else
                    dxyz0(1:3)= this%r_xyz0(i+1,1:3)-this%r_xyz0(i-1,1:3)
                    dxyz(1:3) = this%r_xyz(i+1,1:3)-this%r_xyz(i-1,1:3)
                endif
                dl0       = dsqrt(dxyz0(1)**2+dxyz0(2)**2+dxyz0(3)**2)
                dl        = dsqrt(dxyz(1)**2+dxyz(2)**2+dxyz(3)**2)
                xlmn0(1:3)= dxyz0(1:3)/dl0
                xlmn(1:3) = dxyz(1:3)/dl
                call this%RotateMatrix(i,xlmn0,xlmn)
                call this%Self_RotateMatrix(i,this%r_xyz(i,4),xlmn)
                n_theta  = floor(2*pi*maxval(this%r_xyz0(:,4))/maxval(this%r_xyz0(:,5)))
                if (n_theta .le. 3) n_theta = 3
            enddo

            do iEL = 1,this%v_nelmts
                i = this%vtor(iEL)
                    ! update xyz
                    temp_xyz(1:3) = matmul(this%rotMat(i,:,:), (/this%v_Exyz0(iEL,1), 0.0d0, this%v_Exyz0(iEL,3)/))
                    temp_xyz(1:3) = matmul(this%self_rotMat(i,:,:), temp_xyz)
                    ! Default The fake point on the plane (this%v_Exyz0(iEL,1:3) = (x,y,z)) is in the same plane as the point on the centre axis(this%r_xyz0(i,1:3) = (0,y,0)).
                    temp_xyzoffset(1:3) = this%r_xyz(i,1:3) - this%r_xyz0(i,1:3)
                    this%v_Exyz(iEL,1:3) = temp_xyz(1:3)+temp_xyzoffset(1:3)
                    this%v_Exyz(iEL,2) = this%v_Exyz(iEL,2) + this%v_Exyz0(iEL,2)

                    ! update vel
                    self_rot_omega(1:3) = this%r_vel(i,4)*xlmn(1:3)
                    call cross_product(self_rot_omega,temp_xyz,temp_vel)
                    this%v_Evel(iEL,1:3) = temp_vel(1:3)+this%r_vel(i,1:3)
                !update area
                if ((iEL .ne. 1) .or. (iEL .ne. this%v_nelmts)) then
                    Ea_A = iEL-1
                    if ((i .ne. this%vtor(iEL-1)).or.(Ea_A .eq. 1)) Ea_A = iEL+(n_theta-1)
                    Ea_B = iEL-n_theta
                    if (Ea_B .le. 1) Ea_B = 1
                    Ea_C = iEL+1
                    if ((i .ne. this%vtor(iEL+1)).or.(Ea_C .eq. this%v_nelmts)) Ea_C = iEL-(n_theta-1)
                    Ea_D = iEL+n_theta
                    if (Ea_D .ge. this%v_nelmts) Ea_D = this%v_nelmts
                    call cpt_area(this%v_Exyz(Ea_A,1:3),this%v_Exyz(Ea_B,1:3),this%v_Exyz(Ea_D,1:3),E_left)
                    call cpt_area(this%v_Exyz(Ea_B,1:3),this%v_Exyz(Ea_C,1:3),this%v_Exyz(Ea_D,1:3),E_right)
                    this%v_Ea(iEL) = (E_left+E_right)*0.5d0
                endif
            enddo
            this%v_Ea(1) = 0.d0
            this%v_Ea(this%v_nelmts) = 0.d0


        endif
    end subroutine UpdateElmtInfo_

    subroutine UpdateElmtPosVel_(this,dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
        IMPLICIT NONE
        class(BeamBody), intent(inout) :: this
        real(8),intent(in):: dt,dh,denIn,Uref
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: xGrid(xDim),yGrid(yDim),zGrid(zDim),den(zDim,yDim,xDim)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        real(8)::x0,y0,z0,detx,dety,detz
        integer::i0,j0,k0
    end subroutine UpdateElmtPosVel_

    subroutine UpdateElmtInterp_(this,dh,zDim,yDim,xDim,xGrid,yGrid,zGrid)
        IMPLICIT NONE
        class(BeamBody), intent(inout) :: this
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
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,k,x,y,z,s,rx,ry,rz,detx,dety,detz,ix,jy,kz)
        do  iEL=1,this%v_nelmts
            call minloc_fast(this%v_Exyz(iEL,1), x0, i0, invdh, i, detx)
            call minloc_fast(this%v_Exyz(iEL,2), y0, j0, invdh, j, dety)
            call minloc_fast(this%v_Exyz(iEL,3), z0, k0, invdh, k, detz)
            call trimedindex(i, xDim, ix, boundaryConditions(1:2))
            call trimedindex(j, yDim, jy, boundaryConditions(3:4))
            call trimedindex(k, zDim, kz, boundaryConditions(5:6))
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
        SUBROUTINE minloc_fast(x, x0, i0, invdh, index, offset)
            implicit none
            integer, intent(in):: i0
            real(8), intent(in):: x, x0, invdh
            integer, intent(out):: index
            real(8), intent(out):: offset
            offset = (x - x0)*invdh
            index = floor(offset)
            offset = offset - dble(index)
            index = index + i0
        END SUBROUTINE
    
        SUBROUTINE trimedindex(i, xDim, ix, boundaryConditions)
            USE BoundCondParams
            implicit none
            integer, intent(in):: boundaryConditions(1:2)
            integer, intent(in):: i, xDim
            integer, intent(out):: ix(-1:2)
            integer:: k
            do k=-1,2
                ix(k) = i + k
                if (ix(k)<1) then
                    if(boundaryConditions(1).eq.Periodic) then
                        ix(k) = ix(k) + xDim
                    else if((boundaryConditions(1).eq.SYMMETRIC .or. boundaryConditions(1).eq.wall) .and. ix(k).eq.0) then
                        ix(k) = 2
                    else
                        write(*,*) 'index out of xmin bound', ix(k)
                    endif
                else if(ix(k)>xDim) then
                    if(boundaryConditions(2).eq.Periodic) then
                        ix(k) = ix(k) - xDim
                    else if((boundaryConditions(2).eq.SYMMETRIC .or. boundaryConditions(2).eq.wall) .and. ix(k).eq.xDim+1) then
                        ix(k) = xDim - 1
                    else
                        write(*,*) 'index out of xmax bound', ix(k)
                    endif
                endif
            enddo
        END SUBROUTINE trimedindex
        SUBROUTINE my_minloc(x, array, len, uniform, index) ! return the array(index) <= x < array(index+1)
            implicit none
            integer:: len, index, count, step, it
            real(8):: x, array(len)
            logical:: uniform
            if (.not.uniform) then
                if (x<array(1) .or. x>array(len)) then
                    write(*, *) 'index out of bounds when searching my_minloc', x, '[', array(1), array(len), ']'
                    stop
                endif
                index = 1
                count = len
                do while(count > 0)
                    step = count / 2
                    it = index + step
                    if (array(it) < x) then
                        index = it + 1
                        count = count - (step + 1)
                    else
                        count = step
                    endif
                enddo
                if (array(index)>x) then
                    index = index - 1
                endif
            else
                index = 1 + int((x - array(1))/(array(len)-array(1))*dble(len-1))
                !int -1.1 -> -1; 1.1->1; 1.9->1
                if (index<1 .or. index>len) then
                    write(*, *) 'index out of bounds when searching my_minloc', x, '[', array(1), array(len), ']'
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
            allocate(Beam(iFish)%r_force(Beam(iFish)%r_npts,6))
            Beam(iFish)%r_force(:,:) = 0.0d0
        enddo
        do  while( iterLBM<m_maxIterIB .and. dmaxLBM>m_dtolLBM)
            dmaxLBM=0.0d0
            dsum=0.0d0
            do iFish=1,m_nFish
                call Beam(iFish)%PenaltyForce(dh,dt,denIn,zDim,yDim,xDim,uuu,den,force,tol,ntol)
                tolsum = tolsum + tol
                dsum = dsum + Uref*ntol
            enddo
            dmaxLBM=tolsum/dsum
            iterLBM=iterLBM+1
        enddo
    END SUBROUTINE

    SUBROUTINE PenaltyForce_(this,dh,dt,denIn,zDim,yDim,xDim,uuu,den,force,tolerance,ntolsum)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        class(BeamBody), intent(inout) :: this
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
        !$OMP PARALLEL DO SCHEDULE(STATIC)
        !$OMP PRIVATE(iEL,x,y,z,rx,ry,rz,ix,jy,kz,velElem,velElemIB,forceElemTemp,forceTemp)
        !$OMP REDUCTION(+:tolerance)
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
            tolerance = tolerance + dsqrt(sum((velElem(1:3)-velElemIB(1:3))**2))
            !$OMP critical
            ! update beam load, momentum is not included
            this%r_force(this%vtor(iEL),1:3) = this%r_force(this%vtor(iEL),1:3) + forceElemTemp
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

    subroutine Beam_BuildStructured(this)
        implicit none
        class(BeamBody), intent(inout) :: this
        call this%Readreal()
        call this%Struc_Virtual_Elmts()
    end subroutine Beam_BuildStructured

    subroutine Beam_ReadUnstructured(this)
        implicit none
        class(BeamBody), intent(inout) :: this
        call this%Unstruc_Virtual_Elmts()
        this%r_npts = 1
        allocate(this%r_xyz0(1,5))
        this%r_xyz0(1,1:3) = this%v_Pxyz0(1,1:3)
        this%r_xyz0(1,4:5) = 0.0d0
    end subroutine Beam_ReadUnstructured

    subroutine Section_RotateMatrix(this,i,lmn0,lmn)
        implicit none
        class(BeamBody), intent(inout) :: this
        integer, intent(in) :: i
        ! real(8):: lmn0(3),lmn(3),angle,pi
        ! real(8):: v_0(3),e1(3),e2(3),e3(3),r1(3,3)
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
        ! this%rotMat(i,:,:) = r1
        ! contains
        real(8), intent(in):: lmn0(3),lmn(3)
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
        this%rotMat(i,:,:) = r1
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
    subroutine Section_Self_RotateMatrix(this,i,angle,lmn)
        implicit none
        class(BeamBody), intent(inout) :: this
        real(8), intent(in):: angle,lmn(3)
        integer, intent(in):: i
        real(8):: r1(3,3)
        call quaternion_rotate(angle,lmn,r1)
        this%self_rotMat(i,:,:) = r1
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
        class(BeamBody), intent(inout) :: this
        integer :: fileiD = 111, num, i, temp_nelmts, temp_prop(4)
        character(LEN=1000) :: buffer
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
        integer :: iEL, i, j, k
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
        integer :: iEL, i, j, k
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
        class(BeamBody), intent(inout) :: this
        integer:: i, node, fileiD = 111
        open(unit=fileiD, file = this%filename )
        rewind(fileiD)
        read(fileiD,*)
        read(fileiD,*)
        read(fileiD,*)
        read(fileiD,*)  this%r_npts
        allocate(this%r_xyz0(this%r_npts,5))
        do    i= 1, this%r_npts
            read(fileiD,*) node,this%r_xyz0(node,1),this%r_xyz0(node,2),this%r_xyz0(node,3),this%r_xyz0(node,4),this%r_xyz0(node,5)
        enddo
        close(fileiD)
    end subroutine Readreal_
    subroutine Struc_Virtual_Elmts_(this)
        implicit none
        class(BeamBody), intent(inout) :: this
        call this%Struc_Cylinder_Elmts()
    end subroutine Struc_Virtual_Elmts_
    subroutine Struc_Cylinder_Elmts_(this)
        implicit none
        class(BeamBody), intent(inout) :: this
        real(8) :: pi
        real(8) :: theta, rho, height
        integer :: n_height,n_theta,n_radius
        integer :: i,j,k,num
        integer :: n_circle, n_rectangle
        pi=4.0d0*datan(1.0d0)

        this%v_nelmts = 0
        n_height = this%r_npts
        n_radius = 1+floor(maxval(this%r_xyz0(:,4))/maxval(this%r_xyz0(:,5)))
        n_theta  = floor(2*pi*maxval(this%r_xyz0(:,4))/maxval(this%r_xyz0(:,5)))
        if (n_radius .eq. 1) n_radius = 2
        if (n_theta  .le. 3) n_theta = 3
        
        n_circle = ((n_radius-1) * n_theta) + 1
        n_rectangle = (n_height - 2) * n_theta
        this%v_nelmts = n_circle * 2 + n_rectangle

        allocate(this%v_Exyz0(1:this%v_nelmts, 1:3))
        this%v_Exyz0 = 0.0d0
        allocate(this%vtor(1:this%v_nelmts))
        this%vtor = 0

        ! Given node v_Exyz0 coordinate value
        num = 1
        do i = 1, this%r_npts
            if (i .eq. 1) then
                this%v_Exyz0(num,2) = this%r_xyz0(i,2)
                this%vtor(num) = i
                num = num + 1
                do j = 1, n_radius - 1
                    rho = j * this%r_xyz0(i,4) / (n_radius - 1)
                    height = this%r_xyz0(i,2)
                    do k = 1, n_theta
                        theta = (k - 1) * 2.0 * pi / n_theta
                        this%v_Exyz0(num,1) = rho * cos(theta)
                        this%v_Exyz0(num,2) = height
                        this%v_Exyz0(num,3) = rho * sin(theta)
                        this%vtor(num) = i
                        num = num+1
                    enddo
                enddo
            elseif (i .eq. this%r_npts) then
                do j = 1, n_radius - 1
                    rho = (n_radius - j) * this%r_xyz0(i,4) / (n_radius - 1)
                    height = this%r_xyz0(i,2)
                    do k = 1, n_theta
                        theta = (k - 1) * 2.0 * pi / n_theta
                        this%v_Exyz0(num,1) = rho * cos(theta)
                        this%v_Exyz0(num,2) = height
                        this%v_Exyz0(num,3) = rho * sin(theta)
                        this%vtor(num) = i
                        num = num+1
                    enddo
                enddo
                this%v_Exyz0(num,2) = this%r_xyz0(i,2)
                this%vtor(num) = i
            else
                rho = this%r_xyz0(i,4)
                height = this%r_xyz0(i,2)
                do k = 1, n_theta
                    theta = (k - 1) * 2.0 * pi / n_theta
                    this%v_Exyz0(num,1) = rho * cos(theta)
                    this%v_Exyz0(num,2) = height
                    this%v_Exyz0(num,3) = rho * sin(theta)
                    this%vtor(num) = i
                    num = num+1
                enddo
            endif
        enddo

    end subroutine Struc_Cylinder_Elmts_

    subroutine Write_body_(this,iFish,time,Lref,Tref)
        ! to do: generate a temporary mesh
        implicit none
        class(BeamBody), intent(inout) :: this
        integer,intent(in) :: iFish
        real(8),intent(in) :: time,Lref,Tref
        !   -------------------------------------------------------
        real(8):: timeTref
        integer:: i,ElmType
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName,idstr
        integer,parameter:: idfile=100
        !==========================================================================
        timeTref = time/Tref
        !==========================================================================
        write(fileName,'(I10)') nint(timeTref*1d5)
        fileName = adjustr(fileName)
        do  I=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo
    
        ! ElmType = this%fake_ele(1,4)

        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        open(idfile, FILE='./DatBodySpan/BodyFake'//trim(idstr)//'_'//trim(filename)//'.dat')
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        write(idfile, '(A,I7,A,I7,A)', advance='no') 'ZONE N=',this%v_npts,', E=',this%v_nelmts,', DATAPACKING=POINT, ZONETYPE='
        if(ElmType.eq.2) then
            write(idfile, '(A)') 'FELINESEG'
        elseif (ElmType.eq.3) then
            write(idfile, '(A)') 'FETRIANGLE'
        elseif(ElmType.eq.4) then
            write(idfile, '(A)') 'FEQUADRILATERAL'
        endif
        do  i=1,this%v_nelmts
            write(idfile, *)  this%v_Exyz(i,1:3)!/Lref
        enddo
        ! do  i=1,this%v_nelmts
        !     if(ElmType.eq.2) then
        !         write(idfile, *) this%fake_ele(i,1),this%fake_ele(i,2)
        !     elseif(ElmType.eq.3) then
        !         write(idfile, *) this%fake_ele(i,1),this%fake_ele(i,2),this%fake_ele(i,3)
        !     ! elseif(ElmType.eq.4) then
        !     !     write(idfile, *) this%fake_ele(i,1),this%fake_ele(i,2),this%fake_ele(i,3),this%fake_ele(i,4)
        !     else
        !     endif
        ! enddo
        close(idfile)
        !   =============================================
    end subroutine

    subroutine Initialise_bodies(nFish,maxIterIB,dtolLBM,Pbeta,filenames,iBodyModel)
        implicit none
        integer,intent(in)::nFish, maxIterIB,iBodyModel(nFish)
        real(8),intent(in):: dtolLBM, Pbeta
        character(LEN=100),intent(in):: filenames(nFish)
        integer :: iFish
        m_nFish = nFish
        m_maxIterIB = maxIterIB
        m_dtolLBM = dtolLBM
        m_Pbeta = Pbeta
        allocate(Beam(nFish))
        do iFish = 1,nFish
            call Beam(iFish)%Initialise(filenames(iFish),iBodyModel(iFish))
        enddo
    end subroutine Initialise_bodies
    subroutine Write_solid_bodies(time,Lref,Tref)
        implicit none
        integer :: nFish
        real(8) :: Lref,time,Tref
        integer :: iFish
        do iFish = 1,nFish
            call Beam(iFish)%Write_body(iFish,time,Lref,Tref)
        enddo
    end subroutine Write_solid_bodies

end module SolidBody