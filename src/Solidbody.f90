module SolidBody
    implicit none
    private
    ! Immersed boundary method parameters
    integer:: nFish, maxIterIB
    real(8):: dtolLBM, Pbeta
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
        real(8), allocatable :: r_force(:, :)
        !!!virtual body surface
        integer :: v_nelmts
        real(8), allocatable :: v_Exyz(:, :) ! element center (x, y, z)
        integer, allocatable :: v_Ei(:, :) ! element stencial integer index [ix-1,ix,ix1,ix2, iy-1,iy,iy1,iy2, iz-1,iz,iz1,iz2]
        real(8), allocatable :: v_Ew(:, :) ! element stential weight [wx-1, wx, wx1, wx2, wy-1, wy, wy1, wy2, wz-1, wz, wz1, wz2]
        real(8), allocatable :: v_Ea(:) ! element area
        !area center with equal weight on both sides
        real(8), allocatable :: v_Evel(:, :)
        !calculated using central linear and angular velocities
        integer(2),allocatable :: vtor(:,:)
    contains
        procedure :: Initialise => Initialise_
        procedure :: UpdateElmtInfo => UpdateElmtInfo_
        procedure :: PenaltyForce => PenaltyForce_
        procedure :: BuildStructured => Beam_BuildStructured
        procedure :: ReadUnstructured => Beam_ReadUnstructured
        procedure :: InitialSection => Beam_InitialSection

        procedure :: UpdateLoad => Beam_UpdateLoad

        procedure :: StrucFakNod => Structured_Fake_Nodes
        procedure :: StrucFakEle => Structured_Fake_Elements
        procedure :: Cylinder_Nodes => Structured_Cylinder_Nodes
        procedure :: Cylinder_Elements => Structured_Cylinder_Elements

        procedure :: UnstrucFakNod => Unstructured_Fake_Nodes
        procedure :: UnStrucFakEle => Unstructured_Fake_Elements
        procedure :: RotateMatrix => Section_RotateMatrix
        procedure :: Self_RotateMatrix => Section_Self_RotateMatrix

        procedure :: Write_body => Write_body_
    end type BeamBody
    type(BeamBody), allocatable :: Beam(:)
  contains
    subroutine Initialise_(this,filename)
        ! read beam central line file and allocate memory
        implicit none
        class(BeamBody), intent(inout) :: this
        character (LEN=100), intent(in):: filename
        integer :: iFish

        this%fake_meshfile = FakeBeamMeshName
        this%real_npts     = realnodes
        allocate(this%real_xyz0(1:this%real_npts,1:6))
        this%real_xyz0   = realxyzful

        if (ifUnstructured .eq. 0) then
            this%fake_type = fake_tp
            call this%BuildStructured(fake_r,fake_dh)
        elseif (ifUnstructured .eq. 1) then
            call this%ReadUnstructured()
        endif

        allocate(this%faketoreal(this%fake_npts),this%realtofake(this%real_npts),this%fakepids(this%fake_npts))
        this%faketoreal = 0
        this%realtofake = 0
        this%fakepids = 0
        call this%InitialSection()

        allocate(this%fake_xyz(1:this%fake_npts,1:6))
        do i = 1,this%fake_npts
            this%fake_xyz(i,1:3)=this%fake_xyz0(i,1:3)+realXYZ(1:3)
        enddo
        this%fake_xyz(:,4:6)=this%fake_xyz0(:,4:6)
        allocate(this%fake_vel(1:this%fake_npts,1:6))
        this%fake_vel=0.0d0
        allocate(this%fake_ext(1:this%fake_npts,1:6))
        this%fake_ext=0.0d0
        allocate(this%rotMat(1:this%fake_npts,1:3,1:3))
        this%rotMat=0.0d0
        allocate(this%self_rotMat(1:this%fake_npts,1:3,1:3))
        this%self_rotMat=0.0d0
    end subroutine Initialise_

    subroutine FSInteraction_force(dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
        implicit none
        real(8),intent(in):: dt,dh,denIn,Uref
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: xGrid(xDim),yGrid(yDim),zGrid(zDim),den(zDim,yDim,xDim)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        integer :: iFish
        do iFish = 1,nFish
            call Beam(iFish)%UpdateElmtInfo()
        enddo
        call calculate_interaction_force(dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
    end subroutine

    subroutine UpdateElmtInfo_(this) ! update element position, velocity, weight
        implicit none
        class(BeamBody), intent(inout) :: this
        integer :: i,j,k,iEL
        !   compute displacement, velocity, area at surface element center
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,s,nt,x1,x2,y1,y2,ax,ay)
        do  iEL=1,this%v_nelmts
            i=this%v_elmt(iEL,1)
            j=ele(iEL,2)
            nt=ele(iEL,4)

            x1=xyzful(i,1)
            x2=xyzful(j,1)
            y1=xyzful(i,2)
            y2=xyzful(j,2)
            if(nt/=2) write(*,*) 'only support line segments'
            do s=1,Nspan
                posElem(s,iEL,1)=(x1+x2)*0.5d0
                posElem(s,iEL,2)=(y1+y2)*0.5d0
                posElem(s,iEL,3)=xyzful(i,3) + dspan*(s-0.5)
                velElem(s,iEL,1:2)=(velful(i,1:2)+velful(j,1:2))*0.5d0
                velElem(s,iEL,3)=0.d0
            enddo
            ax =(x1-x2)
            ay =(y1-y2)
            areaElem(iEL)=dsqrt( ax*ax + ay*ay) * dspan
        enddo
        !$OMP END PARALLEL DO
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
        do  while( iterLBM<maxIterIB .and. dmaxLBM>dtolLBM)
            dmaxLBM=0.0d0
            dsum=0.0d0
            do iFish=1,nFish
                call Beam(iFish)%PenaltyForce(dh,dt,denIn,zDim,yDim,xDim,uuu,den,force,tol,ntol)
                tolsum=tolsum+Uref*tol
                dsum = dsum + ntol
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
            forceElemTemp(1:3) = -Pbeta* 2.0d0*denIn*(velElem-velElemIB)/dt*this%v_Ea(iEL)*dh
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

    subroutine Beam_InitialSection(this)
        implicit none
        class(BeamBody), intent(inout) :: this
        integer :: i,j,num
        num = 1
        do i = 1,this%real_npts
            this%realtofake(i) = num
            do j = 1,this%fake_npts
                if ( i .eq. 1) then
                    if ( this%fake_xyz0(j,2) .le. ((this%real_xyz0(i+1,2)+this%real_xyz0(i,2))/2)) then
                        this%faketoreal(j) = i
                        this%fakepids(num) = j
                        num = num + 1
                    endif
                elseif ( i .eq. this%real_npts) then
                    if ( this%fake_xyz0(j,2) .gt. ((this%real_xyz0(i,2)+this%real_xyz0(i-1,2))/2)) then
                        this%faketoreal(j) = i
                        this%fakepids(num) = j
                        num = num + 1
                    endif
                else
                    if ( this%fake_xyz0(j,2) .le. ((this%real_xyz0(i+1,2)+this%real_xyz0(i,2))/2) .and. this%fake_xyz0(j,2) .gt. ((this%real_xyz0(i,2)+this%real_xyz0(i-1,2))/2)) then
                        this%faketoreal(j) = i
                        this%fakepids(num) = j
                        num = num + 1
                    endif
                endif
            enddo
        enddo
    end subroutine Beam_InitialSection

    subroutine Beam_BuildStructured(this,r,dh)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(in) :: r,dh
        call this%StrucFakNod(r,dh)
        call this%StrucFakEle()
    end subroutine Beam_BuildStructured

    subroutine Beam_ReadUnstructured(this)
        implicit none
        class(Body), intent(inout) :: this
        call this%UnstrucFakNod()
        call this%UnstrucFakEle()
    end subroutine Beam_ReadUnstructured

    subroutine Beam_UpdateLoad(this,updaterealextful)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(inout):: updaterealextful(this%real_npts,6)
        integer :: i,j
        updaterealextful = 0.0d0
        do j = 1,this%fake_npts
            i = this%faketoreal(j)
                    updaterealextful(i,1:6) = updaterealextful(i,1:6) + this%fake_ext(j,1:6)
        enddo
    end subroutine Beam_UpdateLoad

    subroutine Section_RotateMatrix(this,i,lmn0,lmn)
        implicit none
        class(Body), intent(inout) :: this
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
        class(Body), intent(inout) :: this
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

    subroutine Unstructured_Fake_Nodes(this)
        implicit none
        class(Body), intent(inout) :: this
        integer :: fileiD = 111, num, i
        character(LEN=1000) :: buffer
        real(8) :: x,y,z
        open(unit=fileiD, file = this%fake_meshfile )
        do i=1,1000000
            read(fileiD,*) buffer
            buffer = trim(buffer)
            if(buffer(1:6) .eq. '$Nodes') exit
        enddo
        read(fileiD,*) num
        this%fake_npts = num
        allocate(this%fake_xyz0(1:this%fake_npts, 1:6))
        this%fake_xyz0 = 0.0d0
        do i = 1,this%fake_npts
            read(fileiD,*)num,x,y,z
            this%fake_xyz0(i,1) = x
            this%fake_xyz0(i,2) = y
            this%fake_xyz0(i,3) = z
        enddo
        close(fileiD)
    end subroutine Unstructured_Fake_Nodes
    subroutine Unstructured_Fake_Elements(this)
        implicit none
        class(Body), intent(inout) :: this
        integer :: fileiD = 111, num, i, temp_nelmts, temp_prop(4)
        character(LEN=1000) :: buffer
        integer :: l,m,n
        open(unit=fileiD, file = this%fake_meshfile )
            do i=1,1000000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:9) .eq. '$Elements') exit
            enddo
            read(fileiD,*) num
            temp_nelmts = num
            do i = 1,temp_nelmts
                read(fileiD,*)num,temp_prop(1:4)
                if((temp_prop(1)==2) .and. (temp_prop(2)==2) .and. (temp_prop(3)==0)) then
                    this%fake_nelmts = temp_nelmts-num+1
                    exit
                endif
            enddo

            backspace(fileiD)

            allocate(this%fake_ele(1:this%fake_nelmts, 1:5))
            this%fake_ele = 0
            do i = 1,this%fake_nelmts
                read(fileiD,*)num,temp_prop(1:4),l,m,n
                this%fake_ele(i,1) = l
                this%fake_ele(i,2) = m
                this%fake_ele(i,3) = n
            enddo
            this%fake_ele(:,4) = 3
            this%fake_ele(:,5) = 1
        close(fileiD)
    end subroutine Unstructured_Fake_Elements
    subroutine Structured_Fake_Nodes(this,r,dh)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(in) :: r,dh
        select case(this%fake_type)
            case(1)
                call this%Cylinder_Nodes(r,dh)
            case(2)
                ! call triangular_Nodes(this,r,dh)
            case(3)
                ! call quadrangular_Nodes(this,r,dh)
            case default
                write(*, *) 'undefined body shape ', this%fake_type
                stop
        end select
    end subroutine Structured_Fake_Nodes
    subroutine Structured_Fake_Elements(this)
        implicit none
        class(Body), intent(inout) :: this
        select case(this%fake_type)
            case(1)
                call this%Cylinder_Elements()
            case(2)
                ! call triangular_Elements(this)
            case(3)
                ! call quadrangular_Elements(this)
            case default
                write(*, *) 'undefined body shape ', this%fake_type
                stop
        end select
    end subroutine
    subroutine Structured_Cylinder_Nodes(this,r,dh)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(in) :: r,dh
        real(8) :: pi
        real(8) :: theta, ztemp, rho
        integer :: n_height,n_theta,n_radius
        integer :: i,j,k,num
        integer :: n_circle, n_rectangle
        pi=4.0d0*datan(1.0d0)
        n_height = this%real_npts
        n_radius = floor(r/dh)+1
        n_theta  = floor(2*pi*r/dh)
        if (n_radius .eq. 1) n_radius = 2
        if (n_theta  .le. 3) n_radius = 3

        this%n_h = n_height
        this%n_r = n_radius
        this%n_t = n_theta

        n_circle = ((n_radius-1) * n_theta) + 1 ! Node number on top(bottom) surface (circle) of cylinder
        n_rectangle = (n_height - 2) * n_theta ! Node number on side surface (rectangle) of cylinder
        this%fake_npts = n_circle * 2 + n_rectangle ! Total node number

        allocate(this%fake_xyz0(1:this%fake_npts, 1:6))
        this%fake_xyz0 = 0.0d0

        ! Given node fake_xyz0 coordinate value
        num = 1
        this%fake_xyz0(num,:)=0.0d0
        num = num+1
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = j * r / (n_radius-1)
                this%fake_xyz0(num,1) = rho * cos(theta)
                this%fake_xyz0(num,2) = 0.0d0
                this%fake_xyz0(num,3) = rho * sin(theta)
                num = num+1
            end do
        enddo
        do k = 1, n_height-2
            ztemp = this%real_xyz0(k+1,2)
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                this%fake_xyz0(num,1) = r * cos(theta)
                this%fake_xyz0(num,2) = ztemp
                this%fake_xyz0(num,3) = r * sin(theta)
                num = num+1
            end do
        end do
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = (n_radius-j)* r / (n_radius-1)
                this%fake_xyz0(num,1) = rho * cos(theta)
                this%fake_xyz0(num,2) = this%real_xyz0(this%real_npts,2)
                this%fake_xyz0(num,3) = rho * sin(theta)
                num = num+1
            end do
        enddo
        this%fake_xyz0(num,1) = 0.0d0
        this%fake_xyz0(num,2) = this%real_xyz0(this%real_npts,2)
        this%fake_xyz0(num,3) = 0.0d0
        
    end subroutine Structured_Cylinder_Nodes
    subroutine Structured_Cylinder_Elements(this)
        implicit none
        class(Body), intent(inout) :: this
        integer :: n_height,n_theta,n_radius
        integer :: i,j,num
        integer :: n_longitude

        n_height = this%n_h
        n_radius = this%n_r
        n_theta  = this%n_t

        n_longitude = n_height + n_radius*2 - 2 ! Node number on one longitude line (include top and bottom node)
        this%fake_nelmts = n_theta + n_theta * (n_longitude - 2 - 1) * 2 + n_theta ! Total element number
        allocate(this%fake_ele(1:this%fake_nelmts, 1:5))
        this%fake_ele = 0
        
        ! Given and output triangle element number (element face normal direction outward)
        
        num = 1
        do i = 1, n_theta-1
            this%fake_ele(num,1) = 1
            this%fake_ele(num,2) = i+2
            this%fake_ele(num,3) = i+1
            num = num+1 ! Element around bottom circle center
        enddo
        this%fake_ele(num,1) = 1
        this%fake_ele(num,2) = 2
        this%fake_ele(num,3) = n_theta+1
        num = num+1 ! Element around bottom circle center
        do j = 1, n_longitude-3
            do i = 1, n_theta-1
                this%fake_ele(num,1) = 1+((j-1)*n_theta+i)
                this%fake_ele(num,2) = 1+((j-1)*n_theta+i+1)
                this%fake_ele(num,3) = 1+( j   *n_theta+i)
                num = num+1 ! Lower triangula
                this%fake_ele(num,1) = 1+((j-1)*n_theta+i+1)
                this%fake_ele(num,2) = 1+( j   *n_theta+i+1)
                this%fake_ele(num,3) = 1+( j   *n_theta+i)
                num = num+1 ! Upper triangular
            enddo
            this%fake_ele(num,1) = 1+((j-1)*n_theta+n_theta)
            this%fake_ele(num,2) = 1+((j-1)*n_theta+1)
            this%fake_ele(num,3) = 1+( j   *n_theta+n_theta)
            num = num+1 ! Lower triangula
            this%fake_ele(num,1) = 1+((j-1)*n_theta+1)
            this%fake_ele(num,2) = 1+( j   *n_theta+1)
            this%fake_ele(num,3) = 1+( j   *n_theta+n_theta)
            num = num+1 ! Upper triangular
        enddo
        do i = 1, n_theta - 1
            this%fake_ele(num,1) = this%fake_npts
            this%fake_ele(num,2) = this%fake_npts-(n_theta-1-i)-2
            this%fake_ele(num,3) = this%fake_npts-(n_theta-1-i)-1
            num = num+1 ! Element around top circle center
        enddo
            this%fake_ele(num,1) = this%fake_npts
            this%fake_ele(num,2) = this%fake_npts-1
            this%fake_ele(num,3) = this%fake_npts-n_theta
        num = num+1 ! Element around top circle center

        this%fake_ele(:,4) = 3
        this%fake_ele(:,5) = 1

    end subroutine Structured_Cylinder_Elements

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
    
        ElmType = this%fake_ele(1,4)

        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        open(idfile, FILE='./DatBodySpan/BodyFake'//trim(idstr)//'_'//trim(filename)//'.dat')
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        write(idfile, '(A,I7,A,I7,A)', advance='no') 'ZONE N=',this%fake_npts,', E=',this%fake_nelmts,', DATAPACKING=POINT, ZONETYPE='
        if(ElmType.eq.2) then
            write(idfile, '(A)') 'FELINESEG'
        elseif (ElmType.eq.3) then
            write(idfile, '(A)') 'FETRIANGLE'
        elseif(ElmType.eq.4) then
            write(idfile, '(A)') 'FEQUADRILATERAL'
        endif
        do  i=1,this%fake_npts
            write(idfile, *)  this%fake_xyz(i,1:3)/Lref
        enddo
        do  i=1,this%fake_nelmts
            if(ElmType.eq.2) then
                write(idfile, *) this%fake_ele(i,1),this%fake_ele(i,2)
            elseif(ElmType.eq.3) then
                write(idfile, *) this%fake_ele(i,1),this%fake_ele(i,2),this%fake_ele(i,3)
            ! elseif(ElmType.eq.4) then
            !     write(idfile, *) this%fake_ele(i,1),this%fake_ele(i,2),this%fake_ele(i,3),this%fake_ele(i,4)
            else
            endif
        enddo
        close(idfile)
        !   =============================================
    end subroutine

    subroutine Initialise_bodies(filenames)
        implicit none
        character(LEN=100),intent(in):: filenames(nFish)
        integer :: iFish
        allocate(Beam(nFish))
        do iFish = 1,nFish
            call Beam(iFish)%Initialise(filenames(iFish))
        enddo
    end subroutine Initialise_bodies
    subroutine Write_solid_bodies(time,Lref,Tref)
        use BodyWorkSpace
        implicit none
        integer :: nFish
        real(8) :: Lref,time,Tref
        integer :: iFish
        do iFish = 1,nFish
            call Beam(iFish)%Write_body(iFish,time,Lref,Tref)
        enddo
    end subroutine Write_solid_bodies

end module SolidBody