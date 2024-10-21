module SolidBody
    implicit none
    private
    integer :: real_npts,nBeam
    character(LEN=100) :: fake_meshfile
    real(8), allocatable :: real_xyzful0(:, :)
    public :: Body,Beam,Beam_initialise,Beam_UpdateInfo,Beam_UpdateLoad
    type :: Body
        integer :: fake_type
        integer :: fake_npts, fake_nelmts
        integer :: n_h, n_r, n_t
        real(8), allocatable :: fake_xyz0(:, :)
        real(8), allocatable :: fake_xyz(:, :)
        real(8), allocatable :: fake_vel(:, :)
        real(8), allocatable :: fake_extful(:, :)
        integer, allocatable :: fake_ele(:, :)
        integer, allocatable :: fake_sec(:, :)
        real(8), allocatable :: rotMat(:, :, :)
    contains
        private
        procedure :: StrucFakNod => Structured_Fake_Nodes
        procedure :: StrucFakEle => Structured_Fake_Elements
        procedure :: UnstrucFakNod => Unstructured_Fake_Nodes
        procedure :: UnStrucFakEle => Unstructured_Fake_Elements
        procedure :: RotateMatrix => Section_RotateMatrix
    end type Body
    type(Body) :: Beam
  contains
    subroutine Beam_initialise(filename,realnodes,realxyzful,r,dh,tp,ifunstructured)
        implicit none
        character (LEN=100):: filename
        real(8) :: r, dh
        integer :: realnodes
        real(8) :: realxyzful(realnodes,6)
        integer :: tp, ifunstructured

        fake_meshfile = filename
        real_npts     = realnodes
        allocate(real_xyzful0(1:real_npts,1:6))
        real_xyzful0   = realxyzful

        if (ifunstructured.eq.0) then
            Beam%fake_type = tp
            call Beam_BuildStructured(r,dh)
        elseif (ifunstructured .eq. 1) then
            call Beam_ReadUnstructured()
        endif
        allocate(Beam%fake_sec(1:real_npts,Beam%fake_npts))
        Beam%fake_sec = 0
        call Beam_InitialSection()

        allocate(Beam%fake_xyz(1:Beam%fake_npts,1:3))
        Beam%fake_xyz=Beam%fake_xyz0
        allocate(Beam%fake_extful(1:Beam%fake_npts,1:6))
        Beam%fake_extful=0.0d0
        allocate(Beam%rotMat(1:Beam%fake_npts,1:3,1:3))
        Beam%rotMat=0.0d0
        allocate(Beam%fake_vel(1:Beam%fake_npts, 1:6))
        Beam%fake_vel=0.0d0
    end subroutine Beam_initialise

    subroutine Beam_InitialSection()
        implicit none
        integer :: i,j
        do j = 1,Beam%fake_npts
            do i = 1,real_npts
                if ( i .eq. 1) then
                    if ( Beam%fake_xyz0(j,2) .le. ((real_xyzful0(i+1,2)+real_xyzful0(i,2))/2)) then
                        Beam%fake_sec(i,j) = j
                    endif
                elseif ( i .eq. real_npts) then
                    if ( Beam%fake_xyz0(j,2) .gt. ((real_xyzful0(i,2)+real_xyzful0(i-1,2))/2)) then
                        Beam%fake_sec(i,j) = j
                    endif
                else
                    if ( Beam%fake_xyz0(j,2) .le. ((real_xyzful0(i+1,2)+real_xyzful0(i,2))/2) .and. Beam%fake_xyz0(j,2) .gt. ((real_xyzful0(i,2)+real_xyzful0(i-1,2))/2)) then
                        Beam%fake_sec(i,j) = j
                    endif
                endif
            enddo
        enddo
    end subroutine Beam_InitialSection

    subroutine Beam_BuildStructured(r,dh)
        implicit none
        real(8), intent(in) :: r,dh
        call Beam%StrucFakNod(r,dh)
        call Beam%StrucFakEle()
    end subroutine Beam_BuildStructured

    subroutine Beam_ReadUnstructured()
        implicit none
        call Beam%UnstrucFakNod()
        call Beam%UnstrucFakEle()
    end subroutine Beam_ReadUnstructured

    subroutine Beam_UpdateInfo(updaterealxyzful,updaterealvelful)
        implicit none
        real(8), intent(in) :: updaterealxyzful(real_npts,6),updaterealvelful(real_npts,6)
        call Beam_Updatexyz(updaterealxyzful)
        call Beam_Updatevel(updaterealvelful)
    end subroutine Beam_UpdateInfo

    subroutine Beam_Updatexyz(updaterealxyzful)
        implicit none
        real(8), intent(in):: updaterealxyzful(real_npts,6)
        integer :: i,j
        real(8) :: dxyz0(3),dxyz(3),xlmn0(3),xlmn(3),dl0,dl,temp(3),tempdxyz(3)
        do j = 1,Beam%fake_npts
            do i = 1,real_npts
                if (Beam%fake_sec(i,j) .ne. 0) then
                    if ( i .eq. 1) then
                        dxyz0(1:3)= real_xyzful0(i+1,1:3)-real_xyzful0(i,1:3)
                        dxyz(1:3) = updaterealxyzful(i+1,1:3)-updaterealxyzful(i,1:3)
                    elseif ( i .eq. real_npts) then
                        dxyz0(1:3)= real_xyzful0(i,1:3)-real_xyzful0(i-1,1:3)
                        dxyz(1:3) = updaterealxyzful(i,1:3)-updaterealxyzful(i-1,1:3)
                    else
                        dxyz0(1:3)= real_xyzful0(i+1,1:3)-real_xyzful0(i-1,1:3)
                        dxyz(1:3) = updaterealxyzful(i+1,1:3)-updaterealxyzful(i-1,1:3)
                    endif
                    dl0       = dsqrt(dxyz0(1)**2+dxyz0(2)**2+dxyz0(3)**2)
                    dl        = dsqrt(dxyz(1)**2+dxyz(2)**2+dxyz(3)**2)
                    xlmn0(1:3)= dxyz0(1:3)/dl0
                    xlmn(1:3) = dxyz(1:3)/dl
                    call Beam%RotateMatrix(i,xlmn0(1),xlmn0(2),xlmn0(3),xlmn(1),xlmn(2),xlmn(3))
                    temp(1:3) = matmul(Beam%rotMat(i,:,:), (/Beam%fake_xyz0(j,1), 0.0d0, Beam%fake_xyz0(j,3)/))
                    ! Default The fake point on the plane (Beam%fake_xyz0(j,1:3) = (x,y,z)) is in the same plane as the point on the centre axis(real_xyzful0(i,1:3) = (0,y,0)).
                    tempdxyz(1:3) = updaterealxyzful(i,1:3) - real_xyzful0(i,1:3)
                    Beam%fake_xyz(j,1:3) = temp(1:3)+tempdxyz(1:3)+real_xyzful0(i,1:3)
                endif
            enddo
        enddo
    end subroutine Beam_Updatexyz

    subroutine Beam_Updatevel(updaterealvelful)
        real(8), intent(in) :: updaterealvelful(real_npts,6)
        integer :: i,j
        do j = 1,Beam%fake_npts
            do i = 1,real_npts
                if (Beam%fake_sec(i,j) .ne. 0) then
                    Beam%fake_vel(j,1:6) = updaterealvelful(i,1:6)
                endif
            enddo
        enddo
    end subroutine Beam_Updatevel

    subroutine Beam_UpdateLoad(TwoDextful)
        real(8), intent(inout):: TwoDextful(real_npts,6)
        integer :: i,j
        TwoDextful = 0.0d0
        do j = 1,Beam%fake_npts
            do i = 1,real_npts
                if (Beam%fake_sec(i,j) .ne. 0) then
                    TwoDextful(i,1:6) = TwoDextful(i,1:6) + Beam%fake_extful(j,1:6)
                endif
            enddo
        enddo
    end subroutine Beam_UpdateLoad

    subroutine Section_RotateMatrix(this,i,l0,m0,n0,l,m,n)
        ! The three rotations of the standard quaternion rotation matrix are all clockwise or counterclockwise.
        ! In order to correspond to doyle's angular rotation matrix, one rotation is changed to the opposite 
        ! direction of the other two rotations, and the resulting quaternion rotation matrix is added with 
        ! a negative sign in the non-diagonal position.
        
        ! First rotation in the positive (counterclockwise) direction of the Z-axis
        ! Second rotation in the negative (clockwise) direction of the Y-axis
        ! Third rotation in the positive (counterclockwise) direction of the X-axis
        implicit none
        class(Body), intent(inout) :: this
        integer:: i
        ! real(8):: l0,m0,n0,l,m,n,angle,pi
        ! real(8):: v_0(3),e1(3),e2(3),e3(3),r1(3,3)
        ! pi=4.0*datan(1.0d0)
        ! v_0= (/l0,m0,n0/)
        ! e1 = (/l,m,n/)
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
        real(8):: l0,m0,n0,l,m,n,angle,pi
        real(8):: v_1(3),v_2(3),lamda(0:3),nn(3),r1(3,3)
        pi=4.0*datan(1.0d0)
        v_1=(/l0,m0,n0/)
        v_2=(/l,m,n/)
        nn=0.0d0
        lamda=0.0d0
        angle=0.0d0
        call angle_between_vectors(v_1,v_2)
        if (angle .gt. 1e-3) then
            call normal_vector(v_1,v_2,nn)
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
        else
            r1(1,:) = (/1.0d0,0.0d0,0.0d0/)
            r1(2,:) = (/0.0d0,1.0d0,0.0d0/)
            r1(3,:) = (/0.0d0,0.0d0,1.0d0/)
        endif
        this%rotMat(i,:,:) = r1
        contains
        subroutine normal_vector(A,B,C)
            real(8):: A(3),B(3),C(3)
            real(8) :: norm_C
            C(1)=A(2)*B(3) - A(3)*B(2)
            C(2)=A(3)*B(1) - A(1)*B(3)
            C(3)=A(1)*B(2) - A(2)*B(1)
            norm_C = sqrt(C(1)**2 + C(2)**2 + C(3)**2)
            if (abs(norm_C - 0.0d0).gt.1e-10) then
                C = C / norm_C
            endif
            return
        end subroutine normal_vector
        subroutine angle_between_vectors(A,B)
            real(8), intent(in) :: A(3),B(3)
            real(8) :: dot_product, magnitude1, magnitude2
            dot_product = A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
            magnitude1 = sqrt(A(1)**2+A(2)**2+A(3)**2)
            magnitude2 = sqrt(B(1)**2+B(2)**2+B(3)**2)
            angle = acos(dot_product / (magnitude1 * magnitude2))
            return
        end subroutine angle_between_vectors
    end subroutine Section_RotateMatrix

    subroutine Unstructured_Fake_Nodes(this)
        implicit none
        class(Body), intent(inout) :: this
        integer :: fileiD = 111, num, i
        character(LEN=1000) :: buffer
        real(8) :: x,y,z
        open(unit=fileiD, file = fake_meshfile )
        do i=1,10000
            read(fileiD,*) buffer
            buffer = trim(buffer)
            if(buffer(1:6) .eq. '$Nodes') exit
        enddo
        read(fileiD,*) num
        this%fake_npts = num
        allocate(this%fake_xyz0(1:this%fake_npts, 1:3))
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
        integer :: fileiD = 111, num, i, temp_nelmts, prop(4)
        character(LEN=1000) :: buffer
        integer :: l,m,n
        open(unit=fileiD, file = fake_meshfile )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:9) .eq. '$Elements') exit
            enddo
            read(fileiD,*) num
            temp_nelmts = num
            do i = 1,temp_nelmts
                read(fileiD,*)num,prop(1:4)
                if((prop(1)==2) .and. (prop(2)==2) .and. (prop(3)==0)) then
                    this%fake_nelmts = temp_nelmts-num+1
                    exit
                endif
            enddo
        close(fileiD)

        open(unit=fileiD, file = fake_meshfile )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:9) .eq. '$Elements') exit
            enddo
            read(fileiD,*) buffer
            do i=1,temp_nelmts-this%fake_nelmts
                read(fileiD,*) buffer
            enddo
            allocate(this%fake_ele(1:this%fake_nelmts, 1:5))
            do i = 1,this%fake_nelmts
                read(fileiD,*)num,prop(1:4),l,m,n
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
                call Cylinder_Nodes(this,r,dh)
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
                call Cylinder_Elements(this)
            case(2)
                ! call triangular_Elements(this)
            case(3)
                ! call quadrangular_Elements(this)
            case default
                write(*, *) 'undefined body shape ', this%fake_type
                stop
        end select
    end subroutine
    subroutine Cylinder_Nodes(this,r,dh)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(in) :: r,dh
        real(8) :: pi
        real(8) :: theta, ztemp, rho
        integer :: n_height,n_theta,n_radius
        integer :: i,j,k,num
        integer :: n_circle, n_rectangle
        pi=4.0d0*datan(1.0d0)
        n_height = real_npts
        n_radius = floor(r/dh)+1
        n_theta  = floor(2*pi*r/dh)
        if (n_radius == 1) n_radius = 2

        this%n_h = n_height
        this%n_r = n_radius
        this%n_t = n_theta

        n_circle = ((n_radius-1) * n_theta) + 1 ! Node number on top(bottom) surface (circle) of cylinder
        n_rectangle = (n_height - 2) * n_theta ! Node number on side surface (rectangle) of cylinder
        this%fake_npts = n_circle * 2 + n_rectangle ! Total node number

        allocate(this%fake_xyz0(1:this%fake_npts, 1:3))
    
        ! Given node fake_xyz0 coordinate value
        num = 1
        this%fake_xyz0(num,:)=0.0d0
        num = num+1
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = j * r / (n_radius-1)
                this%fake_xyz0(num,1) = rho * cos(theta)
                this%fake_xyz0(num,2) = rho * sin(theta)
                this%fake_xyz0(num,3) = 0.0d0
                num = num+1
            end do
        enddo
        do k = 1, n_height-2
            ztemp = real_xyzful0(k+1,3)
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                this%fake_xyz0(num,1) = r * cos(theta)
                this%fake_xyz0(num,2) = r * sin(theta)
                this%fake_xyz0(num,3) = ztemp
                num = num+1
            end do
        end do
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = (n_radius-j)* r / (n_radius-1)
                this%fake_xyz0(num,1) = rho * cos(theta)
                this%fake_xyz0(num,2) = rho * sin(theta)
                this%fake_xyz0(num,3) = real_xyzful0(real_npts,3)
                num = num+1
            end do
        enddo
        this%fake_xyz0(num,1) = 0.0d0
        this%fake_xyz0(num,2) = 0.0d0
        this%fake_xyz0(num,3) = real_xyzful0(real_npts,3)
        
    end subroutine Cylinder_Nodes
    subroutine Cylinder_Elements(this)
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

    end subroutine Cylinder_Elements

end module SolidBody
