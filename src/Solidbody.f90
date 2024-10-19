module SolidBody
    implicit none
    private
    integer :: nND,nBeam
    character(LEN=100) :: m_meshfile1
    character(LEN=100) :: m_meshfile2
    public :: nND,Beam,Beam_initialise,Beam_Load,Body
    type :: Body
        integer :: m_type
        integer, allocatable :: m_npts, m_nelmts
        integer :: n_h, n_r, n_t
        real(8), allocatable :: xyz(:, :)
        integer, allocatable :: ele(:, :)
        integer, allocatable :: sec(:, :)
    contains
        private
        procedure :: Nodes => Virtual_Nodes
        procedure :: Elements => Virtual_Elements
    end type Body
    type(Body) :: Beam
    real(8), allocatable :: xyzful(:, :)
  contains
    subroutine Beam_initialise(filename1,filename2,r,dh,tp,ifstructured)
        implicit none
        character (LEN=100):: filename1
        character (LEN=100):: filename2
        integer :: fileiD = 111
        character (LEN=1000):: buffer
        real(8) :: r, dh
        integer :: tp, ifstructured
        ! load mesh information
        m_meshfile1 = filename1
        m_meshfile2 = filename2
        open(unit=fileiD, file = m_meshfile1 )
            read(fileiD,*) buffer
            read(fileiD,*) nND
        close(fileiD)
        ! load points data
        allocate(xyzful(1:nND,1:3))
        call Beam_ReadPoints()

        if (ifstructured.eq.1) then
            Beam%m_type = tp
            ! build beam
            call Beam_Build(r,dh)
        elseif (ifstructured .eq. 0) then
            call Beam_ReadUnstructured()
        endif
        call Beam_Section()

    end subroutine Beam_initialise
    subroutine Beam_ReadPoints()
        implicit none
        integer :: fileiD = 111, tmpid, i
        character(LEN=1000) :: buffer
        open(unit=fileiD, file = m_meshfile1 )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:5) .eq. 'POINT') exit
            enddo
            read(fileiD,*) buffer
            do i = 1,nND
                read(fileiD,*)tmpid,xyzful(i,1),xyzful(i,2),xyzful(i,3)
            enddo
        close(fileiD)
    end subroutine Beam_ReadPoints
    subroutine Beam_Build(r,dh)
        implicit none
        real(8), intent(in) :: r,dh
        call Beam%Nodes(r,dh)
        call Beam%Elements()
    end subroutine Beam_Build
    subroutine Beam_ReadUnstructured()
        implicit none
        integer :: fileiD = 111, tmpid, i, temp_nelmts, prop(4)
        character(LEN=1000) :: buffer
        real(8) :: x,y,z
        integer :: l,m,n
        open(unit=fileiD, file = m_meshfile2 )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:6) .eq. '$Nodes') exit
            enddo
            read(fileiD,*) tmpid
            Beam%m_npts = tmpid
            allocate(Beam%xyz(1:Beam%m_npts, 1:3))
            do i = 1,Beam%m_npts
                read(fileiD,*)tmpid,x,y,z
                Beam%xyz(i,1) = x
                Beam%xyz(i,2) = y
                Beam%xyz(i,3) = z
            enddo

            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:9) .eq. '$Elements') exit
            enddo
            read(fileiD,*) tmpid
            temp_nelmts = tmpid
            do i = 1,temp_nelmts
                read(fileiD,*)tmpid,prop(1:4)
                if((prop(1)==2) .and. (prop(2)==2) .and. (prop(3)==0)) then
                    Beam%m_nelmts = temp_nelmts-tmpid+1
                    exit
                endif
            enddo
        close(fileiD)

        open(unit=fileiD, file = m_meshfile2 )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:9) .eq. '$Elements') exit
            enddo
            read(fileiD,*) buffer
            do i=1,temp_nelmts-Beam%m_nelmts
                read(fileiD,*) buffer
            enddo
            allocate(Beam%ele(1:Beam%m_nelmts, 1:4))
            do i = 1,Beam%m_nelmts
                read(fileiD,*)tmpid,prop(1:4),l,m,n
                Beam%ele(i,1) = l
                Beam%ele(i,2) = m
                Beam%ele(i,3) = n
            enddo
            Beam%ele(i,4) = 3
        close(fileiD)
        open(unit=333, file = '2.dat' )
            do i = 1,Beam%m_npts
                write(333,*)Beam%xyz(i,1),Beam%xyz(i,2),Beam%xyz(i,3)
            enddo
            do i = 1,Beam%m_nelmts
                write(333,*)Beam%ele(i,1),Beam%ele(i,2),Beam%ele(i,3)
            enddo
        close(333)
    end subroutine Beam_ReadUnstructured
    subroutine Beam_Section()
        implicit none
        integer :: i,j
        allocate(Beam%sec(1:nND,Beam%m_npts))
        Beam%sec = 0
        do j = 1,Beam%m_npts
            do i = 1,nND
                if ( i .eq. 1) then
                    if ( Beam%xyz(j,2) .le. ((xyzful(i+1,2)+xyzful(i,2))/2)) then
                        Beam%sec(i,j) = j
                    endif
                endif
                if ( Beam%xyz(j,2) .le. ((xyzful(i+1,2)+xyzful(i,2))/2) .and. Beam%xyz(j,2) .gt. ((xyzful(i,2)+xyzful(i-1,2))/2)) then
                    Beam%sec(i,j) = j
                endif
                if ( i .eq. nND) then
                    if ( Beam%xyz(j,2) .gt. ((xyzful(i,2)+xyzful(i-1,2))/2)) then
                        Beam%sec(i,j) = j
                    endif
                endif
            enddo
        enddo
    end subroutine
    subroutine Beam_Load(ThreeDextful,TwoDextful)
        real(8), intent(in) :: ThreeDextful(Beam%m_npts,6)
        real(8), intent(out):: TwoDextful(nND,6)
        integer :: i,j
        TwoDextful = 0.0d0
        do j = 1,Beam%m_npts
            do i = 1,nND
                if (Beam%sec(i,j) .ne. 0) then
                    TwoDextful(i,1:6) = TwoDextful(i,1:6) + ThreeDextful(j,6)
                endif
            enddo
        enddo
    end subroutine
    subroutine Virtual_Nodes(this,r,dh)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(in) :: r,dh
        select case(this%m_type)
            case(1)
                call Cylinder_Nodes(this,r,dh)
            case(2)
                ! call triangular_Nodes(this,xyzful,r,dh)
            case(3)
                ! call quadrangular_Nodes(this, xyzful, r, dh)
            case default
                write(*, *) 'undefined body shape ', this%m_type
                stop
        end select
    end subroutine
    subroutine Virtual_Elements(this)
        implicit none
        class(Body), intent(inout) :: this
        select case(this%m_type)
            case(1)
                call Cylinder_Elements(this)
            case(2)
                ! call triangular_Elements(this)
            case(3)
                ! call quadrangular_Elements(this, element)
            case default
                write(*, *) 'undefined body shape ', this%m_type
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
        n_height = nND
        n_radius = floor(r/dh)+1
        n_theta  = floor(2*pi*r/dh)
        if (n_radius == 1) n_radius = 2

        this%n_h = n_height
        this%n_r = n_radius
        this%n_t = n_theta

        n_circle = ((n_radius-1) * n_theta) + 1 ! Node number on top(bottom) surface (circle) of cylinder
        n_rectangle = (n_height - 2) * n_theta ! Node number on side surface (rectangle) of cylinder
        this%m_npts = n_circle * 2 + n_rectangle ! Total node number

        allocate(this%xyz(1:this%m_npts, 1:3))
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Given node xyz coordinate value
        num = 1
        this%xyz(num,:)=0.0d0
        num = num+1
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = j * r / (n_radius-1)
                this%xyz(num,1) = rho * cos(theta)
                this%xyz(num,2) = rho * sin(theta)
                this%xyz(num,3) = 0.0d0
                num = num+1
            end do
        enddo
        do k = 1, n_height-2
            ztemp = xyzful(k+1,3)
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                this%xyz(num,1) = r * cos(theta)
                this%xyz(num,2) = r * sin(theta)
                this%xyz(num,3) = ztemp
                num = num+1
            end do
        end do
        do j = 1, n_radius-1
            do i = 1, n_theta
                theta = (i - 1) * 2.0 * pi / n_theta
                rho = (n_radius-j)* r / (n_radius-1)
                this%xyz(num,1) = rho * cos(theta)
                this%xyz(num,2) = rho * sin(theta)
                this%xyz(num,3) = xyzful(nND,3)
                num = num+1
            end do
        enddo
        this%xyz(num,1) = 0.0d0
        this%xyz(num,2) = 0.0d0
        this%xyz(num,3) = xyzful(nND,3)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! open(unit=222,file='1.dat')
        ! do i = 1, this%m_npts
        !     write(222,"(F23.18,F23.18,F23.18)") this%xyz(i,1),this%xyz(i,2),this%xyz(i,3)
        ! end do
        ! close(222)
        
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
        this%m_nelmts = n_theta + n_theta * (n_longitude - 2 - 1) * 2 + n_theta ! Total element number
        allocate(this%ele(1:this%m_nelmts, 1:4))
    
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Given and output triangle element number (element face normal direction outward)
        
        num = 1
        do i = 1, n_theta-1
            this%ele(num,1) = 1
            this%ele(num,2) = i+2
            this%ele(num,3) = i+1
            num = num+1 ! Element around bottom circle center
        enddo
        this%ele(num,1) = 1
        this%ele(num,2) = 2
        this%ele(num,3) = n_theta+1
        num = num+1 ! Element around bottom circle center
        do j = 1, n_longitude-3
            do i = 1, n_theta-1
                this%ele(num,1) = 1+((j-1)*n_theta+i)
                this%ele(num,2) = 1+((j-1)*n_theta+i+1)
                this%ele(num,3) = 1+( j   *n_theta+i)
                num = num+1 ! Lower triangula
                this%ele(num,1) = 1+((j-1)*n_theta+i+1)
                this%ele(num,2) = 1+( j   *n_theta+i+1)
                this%ele(num,3) = 1+( j   *n_theta+i)
                num = num+1 ! Upper triangular
            enddo
            this%ele(num,1) = 1+((j-1)*n_theta+n_theta)
            this%ele(num,2) = 1+((j-1)*n_theta+1)
            this%ele(num,3) = 1+( j   *n_theta+n_theta)
            num = num+1 ! Lower triangula
            this%ele(num,1) = 1+((j-1)*n_theta+1)
            this%ele(num,2) = 1+( j   *n_theta+1)
            this%ele(num,3) = 1+( j   *n_theta+n_theta)
            num = num+1 ! Upper triangular
        enddo
        do i = 1, n_theta - 1
            this%ele(num,1) = this%m_npts
            this%ele(num,2) = this%m_npts-(n_theta-1-i)-2
            this%ele(num,3) = this%m_npts-(n_theta-1-i)-1
            num = num+1 ! Element around top circle center
        enddo
            this%ele(num,1) = this%m_npts
            this%ele(num,2) = this%m_npts-1
            this%ele(num,3) = this%m_npts-n_theta
        num = num+1 ! Element around top circle center

        this%ele(:,4) = 3
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! open(unit=222,file='1.dat',position='append')
        ! do i = 1, this%m_nelmts
        !     write(222,"(I5,I5,I5,I5)") this%ele(i,1),this%ele(i,2),this%ele(i,3),this%ele(i,4)
        ! end do
        ! close(222)

    end subroutine

    ! subroutine triangular_Nodes(this,xyzful,L_side,dh)
    !     implicit none
    !     class(Body), intent(inout) :: this
    !     real(8), intent(in) :: L_side,dh,xyzful(nND,3)
    !     real(8) :: pi
    !     real(8) :: slope,invslope,trislope,x,d_side,xc,yc,zc
    !     integer :: n_height,n_side,n_perimeter,n_triangle
    !     integer :: i,j,num
    !     pi=4.0d0*datan(1.0d0)
    !     slope       = sqrt(3.0d0)
    !     invslope    = 1/slope
    !     trislope    = slope/2

    !     xc=xyzful(1,1)
    !     yc=xyzful(1,2)
    !     zc=xyzful(1,3)

    !     n_height = nND
    !     n_side   = floor(L_side/dh)

    !     this%n_h = n_height
    !     this%n_r = n_side

    !     n_perimeter = n_side*3-3
    !     n_triangle = (n_side*(n_side+1))/2
    !     d_side = L_side/n_side

    !     this%m_npts = n_height*n_perimeter + n_triangle*2 - n_perimeter*2
    !     allocate(this%xyz(1:this%m_npts, 1:3))
    !     num = 1
    !     do j = 0,n_side-1
    !         if (j .eq. 0) then
    !             do i = 1,n_height
    !                 this%xyz(num,1) = xc
    !                 this%xyz(num,2) = xyzful(i,3)
    !                 this%xyz(num,3) = zc
    !                 num = num+1
    !             enddo
    !         else
    !             do i = 0,j-1
    !                 x=j*d_side*trislope
    !                 this%xyz(num,1) = xc+x
    !                 this%xyz(num,2) = yc
    !                 this%xyz(num,3) = zc-x*invslope+i*d_side
    !                 num = num+1
    !             enddo
    !             do i = 1,n_height
    !                 x=j*d_side*trislope
    !                 this%xyz(num,1) = xc+x
    !                 this%xyz(num,2) = xyzful(i,3)
    !                 this%xyz(num,3) = zc+x*invslope
    !                 num = num+1
    !             enddo
    !             do i = j-1,0
    !                 x=j*d_side*trislope
    !                 this%xyz(num,1) = xc+x
    !                 this%xyz(num,2) = xyzful(nND,3)
    !                 this%xyz(num,3) = zc-x*invslope+i*d_side
    !                 num = num+1
    !             enddo
    !             do i = n_height-1,2
    !                 x=j*d_side*trislope
    !                 this%xyz(num,1) = xc+x
    !                 this%xyz(num,2) = xyzful(i,3)
    !                 this%xyz(num,3) = zc+x*invslope
    !                 num = num+1
    !             enddo
    !         endif
    !     enddo
    !     do j = 1,n_side-2
    !         do i = 2,n_height-1
    !             x=L_side*trislope
    !             this%xyz(num,1) = xc+x
    !             this%xyz(num,2) = xyzful(i,3)
    !             this%xyz(num,3) = zc-x*invslope+j*d_side
    !             num = num+1
    !         enddo
    !     enddo
    ! end subroutine 
    ! subroutine triangular_Elements(this)
    !     implicit none
    !     class(Body), intent(inout) :: this
    !     integer :: n_tri_layer,n_height,n_side
    !     integer :: i,j,num
        
    !     n_height = this%n_h
    !     n_side   = this%n_r
    !     this%m_nelmts = (n_side-1)**2 + (n_height*n_side)
    !     allocate(this%xyz(1:this%m_nelmts, 1:5))
        
    ! end subroutine
end module SolidBody

program main
   use SolidBody
   implicit none
   character (LEN=100):: filename1
   character (LEN=100):: filename2
   real(8) :: r, dh
   integer :: tp, ifstructured , i
   real(8), allocatable :: triDextful(:,:),biDextful(:,:)
!    filename = 'Beam.dat'
!    r  = 0.5
!    dh = 1
!    tp = 1
!    ifstructured = 1
!    call Beam_initialise(filename, r, dh, tp, ifstructured)
   filename1 = 'Beam.dat'
   filename2 = 'sphere0025(x-axis-is-polar)(Version2).msh'
   r  = 0.5
   dh = 1
   tp = 1
   ifstructured = 0
   call Beam_initialise(filename1,filename2, r, dh, tp, ifstructured)
   allocate(triDextful(1:Beam%m_npts,1:6),biDextful(1:nND,1:6))
   triDextful(:,1:3) = 1.0d0
   triDextful(:,4:6) = 0.0d0
   call Beam_Load(triDextful,biDextful)
   open(unit=444, file = '3.dat' )
    do i = 1,nND
        write(444,*)biDextful(i,1:6)
    enddo
close(444)
end program
