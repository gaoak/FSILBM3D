module SolidBody
    implicit none
    private
    ! Immersed boundary method parameters
    integer:: nFish, maxIterIB
    real(8):: dtolLBM, Pbeta
    ! nFish     number of bodies
    ! maxIterIB maximum number of iterations for IB force calculation
    ! dtolIBM   tolerance for IB force calculation
    ! Pbeta     coefficient in penalty force calculation
    public :: Body,Initialise_bodies,Write_solid_bodies,FSInteraction_force
    type :: Body
        character(LEN=100) :: filename
        integer :: real_npts
        real(8), allocatable :: real_xyz0(:, :)

        integer :: fake_type
        integer :: fake_npts, fake_nelmts
        integer :: n_h, n_r, n_t
        real(8), allocatable :: fake_xyz0(:, :)
        real(8), allocatable :: fake_xyz(:, :)
        real(8), allocatable :: fake_vel(:, :)
        real(8), allocatable :: fake_ext(:, :)
        integer, allocatable :: fake_ele(:, :)
        integer, allocatable :: faketoreal(:) !of size fake_npts
        integer, allocatable :: realtofake(:) ! of size real_npts
        integer, allocatable :: fakepids(:) ! of size fake_npts
        real(8), allocatable :: rotMat(:, :, :)
        real(8), allocatable :: self_rotMat(:, :, :)
    contains
        procedure :: Initialise => Body_Initialise
        procedure :: BuildStructured => Beam_BuildStructured
        procedure :: ReadUnstructured => Beam_ReadUnstructured
        procedure :: InitialSection => Beam_InitialSection

        procedure :: UpdateInfo => Beam_UpdateInfo
        procedure :: Update_xyz_vel => Beam_Update_xyz_vel
        procedure :: UpdateLoad => Beam_UpdateLoad

        procedure :: StrucFakNod => Structured_Fake_Nodes
        procedure :: StrucFakEle => Structured_Fake_Elements
        procedure :: Cylinder_Nodes => Structured_Cylinder_Nodes
        procedure :: Cylinder_Elements => Structured_Cylinder_Elements

        procedure :: UnstrucFakNod => Unstructured_Fake_Nodes
        procedure :: UnStrucFakEle => Unstructured_Fake_Elements
        procedure :: RotateMatrix => Section_RotateMatrix
        procedure :: Self_RotateMatrix => Section_Self_RotateMatrix

        procedure :: write_solid_body => Body_Output
    end type Body
    type(Body), allocatable :: Beam(:)
  contains
    subroutine Body_Initialise(this,filename)
        implicit none
        class(Body), intent(inout) :: this
        character (LEN=100), intent(in):: filename
        integer, intent(in) :: realnodes
        real(8), intent(in) :: realxyzful(realnodes,6),realXYZ(3)
        real(8), intent(in) :: fake_r, fake_dh
        integer, intent(in) :: fake_tp, ifUnstructured
        integer :: i

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
    end subroutine Body_Initialise

    subroutine Beam_InitialSection(this)
        implicit none
        class(Body), intent(inout) :: this
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

    subroutine Beam_UpdateInfo(this,updaterealxyzful,updaterealvelful)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(in) :: updaterealxyzful(this%real_npts,6),updaterealvelful(this%real_npts,6)
        call this%Update_xyz_vel(updaterealxyzful,updaterealvelful)
    end subroutine Beam_UpdateInfo

    subroutine Beam_Update_xyz_vel(this,updaterealxyzful,updaterealvelful)
        implicit none
        class(Body), intent(inout) :: this
        real(8) :: updaterealxyzful(this%real_npts,6),updaterealvelful(this%real_npts,6)
        integer :: i,j,index,begin,end
        real(8) :: dxyz0(3),dxyz(3),xlmn0(3),xlmn(3),dl0,dl,temp_xyz(3),temp_xyzoffset(3)
        real(8) :: self_rot_omega(3),temp_vel(3)
        do i = 1,this%real_npts
            ! update rotMat and self_rotMar
            if ( i .eq. 1) then
                dxyz0(1:3)= this%real_xyz0(i+1,1:3)-this%real_xyz0(i,1:3)
                dxyz(1:3) = updaterealxyzful(i+1,1:3)-updaterealxyzful(i,1:3)
            elseif ( i .eq. this%real_npts) then
                dxyz0(1:3)= this%real_xyz0(i,1:3)-this%real_xyz0(i-1,1:3)
                dxyz(1:3) = updaterealxyzful(i,1:3)-updaterealxyzful(i-1,1:3)
            else
                dxyz0(1:3)= this%real_xyz0(i+1,1:3)-this%real_xyz0(i-1,1:3)
                dxyz(1:3) = updaterealxyzful(i+1,1:3)-updaterealxyzful(i-1,1:3)
            endif
            dl0       = dsqrt(dxyz0(1)**2+dxyz0(2)**2+dxyz0(3)**2)
            dl        = dsqrt(dxyz(1)**2+dxyz(2)**2+dxyz(3)**2)
            xlmn0(1:3)= dxyz0(1:3)/dl0
            xlmn(1:3) = dxyz(1:3)/dl
            call this%RotateMatrix(i,xlmn0,xlmn)
            call this%Self_RotateMatrix(i,updaterealxyzful(i,4),xlmn)

            if (i .eq. this%real_npts) then
                begin = this%realtofake(i)
                end = this%fake_npts
            else
                begin = this%realtofake(i)
                end = this%realtofake(i+1)-1
            endif

            do index = begin,end
                j = this%fakepids(index)
                    ! update xyz
                    temp_xyz(1:3) = matmul(this%rotMat(i,:,:), (/this%fake_xyz0(j,1), 0.0d0, this%fake_xyz0(j,3)/))
                    temp_xyz(1:3) = matmul(this%self_rotMat(i,:,:), temp_xyz)
                    ! Default The fake point on the plane (this%fake_xyz0(j,1:3) = (x,y,z)) is in the same plane as the point on the centre axis(this%real_xyz0(i,1:3) = (0,y,0)).
                    temp_xyzoffset(1:3) = updaterealxyzful(i,1:3) - this%real_xyz0(i,1:3)
                    this%fake_xyz(j,1:3) = temp_xyz(1:3)+temp_xyzoffset(1:3)
                    this%fake_xyz(j,2) = this%fake_xyz(j,2) + this%fake_xyz0(j,2)

                    ! update vel
                    self_rot_omega(1:3) = updaterealvelful(i,4)*xlmn(1:3)
                    call cross_product(self_rot_omega,temp_xyz,temp_vel)
                    this%fake_vel(j,1:3) = temp_vel(1:3)+updaterealvelful(i,1:3)
            enddo
        enddo
    end subroutine Beam_Update_xyz_vel

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

    subroutine Body_Output(this,iFish,time,Lref,Tref)
        implicit none
        class(Body), intent(inout) :: this
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
        if (maxval(isFake_ful) .eq. 1) then
            do iFish = 1,nFish
                call Beam(iFish)%Initialise(FakeBeamMeshName_ful(iFish), nND(iFish), xyzful00(1:nND(iFish),1:6,iFish), XYZ(1:3,iFish), &
                                            fake_r_ful(iFish),fake_dh_ful(iFish),fake_tp_ful(iFish),ifUnstructured_ful(iFish))
            enddo
            Beam_nND_max = maxval(Beam%fake_npts)
            Beam_nEL_max = maxval(Beam%fake_nelmts)
        endif
    end subroutine Initialise_bodies
    subroutine Write_solid_bodies(time,Lref,Tref)
        use BodyWorkSpace
        implicit none
        integer :: nFish
        real(8) :: Lref,time,Tref
        integer :: iFish
        do iFish = 1,nFish
            call Beam(iFish)%write_solid_body(iFish,time,Lref,Tref)
        enddo
    end subroutine Write_solid_bodies

    subroutine FSInteraction_force(dt,dh,denIn,Uref,zDim,yDim,xDim,xGrid,yGrid,zGrid,uuu,den,force)
        implicit none
        real(8),intent(in):: dt,dh,denIn,Uref
        integer,intent(in):: zDim,yDim,xDim
        real(8),intent(in):: xGrid(xDim),yGrid(yDim),zGrid(zDim),den(zDim,yDim,xDim)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        integer :: iFish
        do iFish = 1,nFish
            call Beam(iFish)%UpdateInfo(xyzful(1:nND(iFish),1:6,iFish),velful(1:nND(iFish),1:6,iFish))
        enddo
        call calculate_interaction_force(zDim,yDim,xDim,nFish,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                                         Pbeta,ntolLBM,dtolLBM,force,isUniformGrid)
        do iFish = 1,nFish
            call Beam(iFish)%UpdateLoad(extful(1:nND(iFish),1:6,iFish))
        enddo
    end subroutine

    SUBROUTINE calculate_interaction_force(zDim,yDim,xDim,nFish,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                                           Pbeta,ntolLBM,dtolLBM,force,isUniformGrid)
        use BodyWorkSpace
        IMPLICIT NONE
        integer,intent(in):: zDim,yDim,xDim,nFish,ntolLBM
        real(8),intent(in):: dh,Uref,denIn,dtolLBM,dt,Pbeta
        real(8),intent(in):: den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
        logical,intent(in):: isUniformGrid(1:3)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        !==================================================================================================
        integer:: iFish
        integer:: i,j,k,iEL,nt,iterLBM
        real(8):: dmaxLBM,dsum
        real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az
        real(8):: forceElem(Beam_nEL_max,3,nFish),forceNode(Beam_nND_max,3,nFish),areaElem(Beam_nEL_max,nFish)
        real(8):: posElem(Beam_nEL_max,3,nFish),velElem(Beam_nEL_max,3,nFish),velElemIB(Beam_nEL_max,3,nFish)
        !==================================================================================================
        !   compute displacement, velocity, area at surface element center
        do iFish=1,nFish
            do  iEL=1,Beam(iFish)%fake_nelmts
                i=Beam(iFish)%fake_ele(iEL,1)
                j=Beam(iFish)%fake_ele(iEL,2)
                k=Beam(iFish)%fake_ele(iEL,3)
                nt=Beam(iFish)%fake_ele(iEL,4)

                x1=Beam(iFish)%fake_xyz(i,1)
                x2=Beam(iFish)%fake_xyz(j,1)
                x3=Beam(iFish)%fake_xyz(k,1)
                y1=Beam(iFish)%fake_xyz(i,2)
                y2=Beam(iFish)%fake_xyz(j,2)
                y3=Beam(iFish)%fake_xyz(k,2)
                z1=Beam(iFish)%fake_xyz(i,3)
                z2=Beam(iFish)%fake_xyz(j,3)
                z3=Beam(iFish)%fake_xyz(k,3)

                if(nt==2)then
                    posElem(iEL,1:3,iFish)=(Beam(iFish)%fake_xyz(i,1:3)+Beam(iFish)%fake_xyz(j,1:3))/2.0d0
                    velElem(iEL,1:3,iFish)=(Beam(iFish)%fake_vel(i,1:3)+Beam(iFish)%fake_vel(j,1:3))/2.0d0
                    ax =(x1-x2)
                    ay =(y1-y2)
                    az =(z1-z2)
                    areaElem(iEL,iFish)=dsqrt( ax*ax + ay*ay + az*az)

                elseif(nt==3)then
                    posElem(iEL,1:3,iFish)=(Beam(iFish)%fake_xyz(i,1:3)+Beam(iFish)%fake_xyz(j,1:3)+Beam(iFish)%fake_xyz(k,1:3))/3.0d0
                    velElem(iEL,1:3,iFish)=(Beam(iFish)%fake_vel(i,1:3)+Beam(iFish)%fake_vel(j,1:3)+Beam(iFish)%fake_vel(k,1:3))/3.0d0
                    ax =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.0d0
                    ay =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.0d0
                    az =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.0d0
                    areaElem(iEL,iFish)=dsqrt( ax*ax + ay*ay + az*az)
                else
                        write(*,*)'cell type is not defined'
                endif
            enddo
        enddo

        !**************************************************************************************************
        !**************************************************************************************************
        forceElem(1:Beam_nEL_max,1:3,1:nFish)=0.0d0
        dmaxLBM=1.0d0
        iterLBM=0
        !   ***********************************************************************************************
        do  while( iterLBM<ntolLBM .and. dmaxLBM>dtolLBM)

            dmaxLBM=0.0d0
            dsum=0.0d0

            do iFish=1,nFish

                call calculate_interaction_force_core(zDim,yDim,xDim,Beam(iFish)%fake_nelmts,Beam(iFish)%fake_ele,dh,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                Pbeta,force,isUniformGrid,posElem(1:Beam(iFish)%fake_nelmts,1:3,iFish),velElem(1:Beam(iFish)%fake_nelmts,1:3,iFish), &
                areaElem(1:Beam(iFish)%fake_nelmts,iFish),forceElem(1:Beam(iFish)%fake_nelmts,1:3,iFish),velElemIB(1:Beam(iFish)%fake_nelmts,1:3,iFish))

                dsum=dsum+Uref*Beam(iFish)%fake_nelmts

                do iEL=1,Beam(iFish)%fake_nelmts
                    dmaxLBM=dmaxLBM+dsqrt(sum((velElem(iEL,1:3,iFish)-velElemIB(iEL,1:3,iFish))**2))
                enddo
        !   ***********************************************************************************************
            enddo
            dmaxLBM=dmaxLBM/dsum
            iterLBM=iterLBM+1
        enddo
        !write(*,'(A,I5,A,D20.10)')' iterLBM =',iterLBM,'    dmaxLBM =',dmaxLBM
        !**************************************************************************************************
        !**************************************************************************************************
        !   element force to nodal force
        forceNode(1:Beam_nND_max,1:3,1:nFish)=0.0d0
        do iFish=1,nFish
            Beam(iFish)%fake_ext(1:Beam(iFish)%fake_npts,1:6)=0.0d0
            do    iEL=1,Beam(iFish)%fake_nelmts
                i=Beam(iFish)%fake_ele(iEL,1)
                j=Beam(iFish)%fake_ele(iEL,2)
                k=Beam(iFish)%fake_ele(iEL,3)
                nt=Beam(iFish)%fake_ele(iEL,4)
                forceNode(i,1:3,iFish)=forceNode(i,1:3,iFish)+forceElem(iEL,1:3,iFish)/3.0d0
                forceNode(j,1:3,iFish)=forceNode(j,1:3,iFish)+forceElem(iEL,1:3,iFish)/3.0d0
                forceNode(k,1:3,iFish)=forceNode(k,1:3,iFish)+forceElem(iEL,1:3,iFish)/3.0d0
            enddo
            Beam(iFish)%fake_ext(1:Beam(iFish)%fake_npts,1:3) = forceNode(1:Beam(iFish)%fake_npts,1:3,iFish)
            Beam(iFish)%fake_ext(1:Beam(iFish)%fake_npts,4:6) = 0.0d0
        enddo
    END SUBROUTINE

    SUBROUTINE calculate_interaction_force_core(zDim,yDim,xDim,nEL,ele,dh,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                                                Pbeta,force,isUniformGrid,posElem,velElem,areaElem,forceElem,velElemIB)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        integer,intent(in):: zDim,yDim,xDim,nEL,ele(nEL,5)
        real(8),intent(in):: dh,denIn,dt,Pbeta
        real(8),intent(in):: den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
        logical,intent(in):: isUniformGrid(1:3)
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        !==================================================================================================
        real(8),intent(in):: posElem(nEL,3),velElem(nEL,3)
        real(8),intent(in):: areaElem(nEL)
        real(8),intent(inout)::forceElem(nEL,3)
        real(8),intent(out)::velElemIB(nEL,3)
        !==================================================================================================
        integer:: i,j,k,x,y,z,iEL
        real(8):: rx,ry,rz,Phi,invdh,forcetemp(1:3)
        real(8):: forceElemTemp(nEL,3)
        !==================================================================================================
        !   compute velocity and displacement at IB nodes
        invdh = 1.D0/dh
        !   ***********************************************************************************************
        !   compute the velocity of IB nodes at element center
        do  iEL=1,nEL
            call my_minloc(posElem(iEL,1), xGrid, xDim, isUniformGrid(1), i)
            call my_minloc(posElem(iEL,2), yGrid, yDim, isUniformGrid(2), j)
            call my_minloc(posElem(iEL,3), zGrid, zDim, isUniformGrid(3), k)
            velElemIB(iEL,1:3)=0.0d0
            do x=-1+i,2+i
                rx=Phi((posElem(iEL,1)-xGrid(x))*invdh)
                do y=-1+j,2+j
                    ry=Phi((posElem(iEL,2)-yGrid(y))*invdh)
                    do z=-1+k,2+k
                        rz=Phi((posElem(iEL,3)-zGrid(z))*invdh)
                        velElemIB(iEL,1:3)=velElemIB(iEL,1:3)+uuu(z,y,x,1:3)*rx*ry*rz
                    enddo
                enddo
            enddo
        enddo
        !   ***********************************************************************************************
        !   calculate interaction force
        do  iEL=1,nEL
            if(ele(iEL,4)==3) then
                forceElemTemp(iEL,1:3) = -Pbeta* 2.0d0*denIn*(velElem(iEL,1:3)-velElemIB(iEL,1:3))/dt*areaElem(iEL)*dh
            else
                forceElemTemp(iEL,1:3) = 0.0d0
            endif
            if ((.not. IEEE_IS_FINITE(forceElemTemp(iEL,1))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(iEL,2))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(iEL,3)))) then
                write(*, *) 'Nan found in forceElemTemp', forceElemTemp
                stop
            endif
        enddo
        !   ***********************************************************************************************
        !   calculate Eulerian body force
        do    iEL=1,nEL
            call my_minloc(posElem(iEL,1), xGrid, xDim, isUniformGrid(1), i)
            call my_minloc(posElem(iEL,2), yGrid, yDim, isUniformGrid(2), j)
            call my_minloc(posElem(iEL,3), zGrid, zDim, isUniformGrid(3), k)
            do x=-1+i,2+i
                rx=Phi((posElem(iEL,1)-xGrid(x))*invdh)
                do y=-1+j,2+j
                    ry=Phi((posElem(iEL,2)-yGrid(y))*invdh)
                    do z=-1+k,2+k
                        rz=Phi((posElem(iEL,3)-zGrid(z))*invdh)
                        forcetemp(1:3) = -forceElemTemp(iEL,1:3)*rx*ry*rz*invdh*invdh*invdh
                        ! update velocity
                        uuu(z,y,x,1:3)  = uuu(z,y,x,1:3)+0.5d0*dt*forceTemp(1:3)/den(z,y,x)
                        force(z,y,x,1:3)=force(z,y,x,1:3) + forcetemp(1:3)
                    enddo
                enddo
            enddo
        enddo
        forceElem(1:nEL,1:3) = forceElem(1:nEL,1:3)+forceElemTemp(1:nEL,1:3)
    END SUBROUTINE
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

end module SolidBody