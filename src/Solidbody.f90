module SolidBody
    implicit none
    private
    integer :: real_npts,nBeam
    character(LEN=100) :: fake_meshfile
    public :: real_npts,Body,Beam,Beam_initialise,Beam_Load
    type :: Body
        integer :: fake_type
        integer, allocatable :: fake_npts, fake_nelmts
        integer :: n_h, n_r, n_t
        real(8), allocatable :: facke_xyz(:, :)
        integer, allocatable :: facke_ele(:, :)
        integer, allocatable :: facke_sec(:, :)
    contains
        private
        procedure :: StrucFakNod => Structured_Fake_Nodes
        procedure :: StrucFakEle => Structured_Fake_Elements
        procedure :: UnstrucFakNod => Unstructured_Fake_Nodes
        procedure :: UnStrucFakEle => Unstructured_Fake_Elements
    end type Body
    type(Body) :: Beam
    real(8), allocatable :: real_xyzful(:, :)
  contains
    subroutine Beam_initialise(filename,realnodes,realxyzful,r,dh,tp,ifunstructured)
        implicit none
        character (LEN=100):: filename
        real(8) :: r, dh
        integer :: realnodes
        real(8) :: realxyzful(realnodes,3)
        integer :: tp, ifunstructured

        fake_meshfile = filename
        real_npts     = realnodes
        allocate(real_xyzful(1:real_npts,1:3))
        real_xyzful   = realxyzful

        if (ifunstructured.eq.0) then
            Beam%fake_type = tp
            call Beam_BuildStructured(r,dh)
        elseif (ifunstructured .eq. 1) then
            call Beam_ReadUnstructured()
        endif
        call Beam_Section()

    end subroutine Beam_initialise

    subroutine Beam_Section()
        implicit none
        integer :: i,j
        allocate(Beam%facke_sec(1:real_npts,Beam%fake_npts))
        Beam%facke_sec = 0
        do j = 1,Beam%fake_npts
            do i = 1,real_npts
                if ( i .eq. 1) then
                    if ( Beam%facke_xyz(j,2) .le. ((real_xyzful(i+1,2)+real_xyzful(i,2))/2)) then
                        Beam%facke_sec(i,j) = j
                    endif
                endif
                if ( Beam%facke_xyz(j,2) .le. ((real_xyzful(i+1,2)+real_xyzful(i,2))/2) .and. Beam%facke_xyz(j,2) .gt. ((real_xyzful(i,2)+real_xyzful(i-1,2))/2)) then
                    Beam%facke_sec(i,j) = j
                endif
                if ( i .eq. real_npts) then
                    if ( Beam%facke_xyz(j,2) .gt. ((real_xyzful(i,2)+real_xyzful(i-1,2))/2)) then
                        Beam%facke_sec(i,j) = j
                    endif
                endif
            enddo
        enddo
    end subroutine

    subroutine Beam_ReadPoints()
        implicit none
        integer :: fileiD = 111, num, i
        character(LEN=1000) :: buffer
        character (LEN=100):: filename
        
        open(unit=fileiD, file = filename )
            read(fileiD,*) buffer
            read(fileiD,*) real_npts
        close(fileiD)
        ! load points data
        allocate(real_xyzful(1:real_npts,1:3))
        open(unit=fileiD, file = filename )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:5) .eq. 'POINT') exit
            enddo
            read(fileiD,*) buffer
            do i = 1,real_npts
                read(fileiD,*)num,real_xyzful(i,1),real_xyzful(i,2),real_xyzful(i,3)
            enddo
        close(fileiD)
    end subroutine Beam_ReadPoints

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
    
    subroutine Beam_Load(ThreeDextful,TwoDextful)
        real(8), intent(in) :: ThreeDextful(Beam%fake_npts,6)
        real(8), intent(out):: TwoDextful(real_npts,6)
        integer :: i,j
        TwoDextful = 0.0d0
        do j = 1,Beam%fake_npts
            do i = 1,real_npts
                if (Beam%facke_sec(i,j) .ne. 0) then
                    TwoDextful(i,1:6) = TwoDextful(i,1:6) + ThreeDextful(j,6)
                endif
            enddo
        enddo
    end subroutine

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
        allocate(this%facke_xyz(1:this%fake_npts, 1:3))
        do i = 1,this%fake_npts
            read(fileiD,*)num,x,y,z
            this%facke_xyz(i,1) = x
            this%facke_xyz(i,2) = y
            this%facke_xyz(i,3) = z
        enddo
        close(fileiD)
    end subroutine
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
            allocate(this%facke_ele(1:this%fake_nelmts, 1:4))
            do i = 1,this%fake_nelmts
                read(fileiD,*)num,prop(1:4),l,m,n
                this%facke_ele(i,1) = l
                this%facke_ele(i,2) = m
                this%facke_ele(i,3) = n
            enddo
            this%facke_ele(i,4) = 3
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
    end subroutine
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

        allocate(this%facke_xyz(1:this%fake_npts, 1:3))
    
        ! Given node facke_xyz coordinate value
        num = 1
        this%facke_xyz(num,:)=0.0d0
        num = num+1
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = j * r / (n_radius-1)
                this%facke_xyz(num,1) = rho * cos(theta)
                this%facke_xyz(num,2) = rho * sin(theta)
                this%facke_xyz(num,3) = 0.0d0
                num = num+1
            end do
        enddo
        do k = 1, n_height-2
            ztemp = real_xyzful(k+1,3)
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                this%facke_xyz(num,1) = r * cos(theta)
                this%facke_xyz(num,2) = r * sin(theta)
                this%facke_xyz(num,3) = ztemp
                num = num+1
            end do
        end do
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = (n_radius-j)* r / (n_radius-1)
                this%facke_xyz(num,1) = rho * cos(theta)
                this%facke_xyz(num,2) = rho * sin(theta)
                this%facke_xyz(num,3) = real_xyzful(real_npts,3)
                num = num+1
            end do
        enddo
        this%facke_xyz(num,1) = 0.0d0
        this%facke_xyz(num,2) = 0.0d0
        this%facke_xyz(num,3) = real_xyzful(real_npts,3)
        
    end subroutine
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
        allocate(this%facke_ele(1:this%fake_nelmts, 1:4))
    
        
        ! Given and output triangle element number (element face normal direction outward)
        
        num = 1
        do i = 1, n_theta-1
            this%facke_ele(num,1) = 1
            this%facke_ele(num,2) = i+2
            this%facke_ele(num,3) = i+1
            num = num+1 ! Element around bottom circle center
        enddo
        this%facke_ele(num,1) = 1
        this%facke_ele(num,2) = 2
        this%facke_ele(num,3) = n_theta+1
        num = num+1 ! Element around bottom circle center
        do j = 1, n_longitude-3
            do i = 1, n_theta-1
                this%facke_ele(num,1) = 1+((j-1)*n_theta+i)
                this%facke_ele(num,2) = 1+((j-1)*n_theta+i+1)
                this%facke_ele(num,3) = 1+( j   *n_theta+i)
                num = num+1 ! Lower triangula
                this%facke_ele(num,1) = 1+((j-1)*n_theta+i+1)
                this%facke_ele(num,2) = 1+( j   *n_theta+i+1)
                this%facke_ele(num,3) = 1+( j   *n_theta+i)
                num = num+1 ! Upper triangular
            enddo
            this%facke_ele(num,1) = 1+((j-1)*n_theta+n_theta)
            this%facke_ele(num,2) = 1+((j-1)*n_theta+1)
            this%facke_ele(num,3) = 1+( j   *n_theta+n_theta)
            num = num+1 ! Lower triangula
            this%facke_ele(num,1) = 1+((j-1)*n_theta+1)
            this%facke_ele(num,2) = 1+( j   *n_theta+1)
            this%facke_ele(num,3) = 1+( j   *n_theta+n_theta)
            num = num+1 ! Upper triangular
        enddo
        do i = 1, n_theta - 1
            this%facke_ele(num,1) = this%fake_npts
            this%facke_ele(num,2) = this%fake_npts-(n_theta-1-i)-2
            this%facke_ele(num,3) = this%fake_npts-(n_theta-1-i)-1
            num = num+1 ! Element around top circle center
        enddo
            this%facke_ele(num,1) = this%fake_npts
            this%facke_ele(num,2) = this%fake_npts-1
            this%facke_ele(num,3) = this%fake_npts-n_theta
        num = num+1 ! Element around top circle center

        this%facke_ele(:,4) = 3

    end subroutine

end module SolidBody

! program main
!    use SolidBody
!    implicit none
!    character (LEN=100):: filename1
!    character (LEN=100):: filename2
!    real(8) :: r, dh
!    integer :: tp, ifstructured , i
!    real(8), allocatable :: triDextful(:,:),biDextful(:,:)
! !    filename = 'Beam.dat'
! !    r  = 0.5
! !    dh = 1
! !    tp = 1
! !    ifstructured = 1
! !    call Beam_initialise(filename, r, dh, tp, ifstructured)
!    filename1 = 'Beam.dat'
!    filename2 = 'sphere0025(x-axis-is-polar)(Version2).msh'
!    r  = 0.5
!    dh = 1
!    tp = 1
!    ifstructured = 0
!    call Beam_initialise(filename1,filename2, r, dh, tp, ifstructured)
!    allocate(triDextful(1:Beam%fake_npts,1:6),biDextful(1:real_npts,1:6))
!    triDextful(:,1:3) = 1.0d0
!    triDextful(:,4:6) = 0.0d0
!    call Beam_Load(triDextful,biDextful)
!    open(unit=444, file = '3.dat' )
!     do i = 1,real_npts
!         write(444,*)biDextful(i,1:6)
!     enddo
! close(444)
! end program
