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
        real(8), allocatable :: fake_xyzIB(:, :)
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

    end subroutine Beam_initialise

    subroutine Beam_InitialSection()
        implicit none
        integer :: i,j
        
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

    subroutine Beam_UpdateInfo(updaterealxyzful,updaterealvelful,updaterealxyzIBful)
        implicit none
        real(8), intent(in) :: updaterealxyzful(real_npts,6),updaterealvelful(real_npts,6),updaterealxyzIBful(real_npts,6)
        call Beam_Updatexyz(updaterealxyzful)
        call Beam_Updatevel(updaterealvelful)
        call Beam_UpdatexyzIB(updaterealxyzIBful)
    end subroutine Beam_UpdateInfo

    subroutine Beam_Updatexyz(updaterealxyzful)
        implicit none
        real(8), intent(in):: updaterealxyzful(real_npts,6)
        integer :: i,j
        real(8) :: dxyz0(3),dxyz(3),xlmn0(3),xlmn(3),dl0,dl,temp(3),tempdxyz(3)
        
    end subroutine Beam_Updatexyz

    subroutine Beam_Updatevel(updaterealvelful)
        real(8), intent(in) :: updaterealvelful(real_npts,6)
        integer :: i,j
        
    end subroutine Beam_Updatevel

    subroutine Beam_UpdatexyzIB(updaterealxyzIBful)
        real(8), intent(in) :: updaterealxyzIBful(real_npts,6)
        integer :: i,j
        
    end subroutine Beam_UpdatexyzIB

    subroutine Beam_UpdateLoad(TwoDextful)
        real(8), intent(inout):: TwoDextful(real_npts,6)
        integer :: i,j
        
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
        
        real(8):: l0,m0,n0,l,m,n,angle,pi
        real(8):: v_1(3),v_2(3),lamda(0:3),nn(3),r1(3,3)
        
    end subroutine Section_RotateMatrix

    subroutine Unstructured_Fake_Nodes(this)
        implicit none
        class(Body), intent(inout) :: this
        integer :: fileiD = 111, num, i
        character(LEN=1000) :: buffer
        real(8) :: x,y,z
        
    end subroutine Unstructured_Fake_Nodes
    subroutine Unstructured_Fake_Elements(this)
        implicit none
        class(Body), intent(inout) :: this
        integer :: fileiD = 111, num, i, temp_nelmts, prop(4)
        character(LEN=1000) :: buffer
        integer :: l,m,n
        
    end subroutine Unstructured_Fake_Elements
    subroutine Structured_Fake_Nodes(this,r,dh)
        implicit none
        class(Body), intent(inout) :: this
        real(8), intent(in) :: r,dh
        
    end subroutine Structured_Fake_Nodes
    subroutine Structured_Fake_Elements(this)
        implicit none
        class(Body), intent(inout) :: this
        
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
           
    end subroutine Cylinder_Nodes
    subroutine Cylinder_Elements(this)
        implicit none
        class(Body), intent(inout) :: this
        integer :: n_height,n_theta,n_radius
        integer :: i,j,num
        integer :: n_longitude

    end subroutine Cylinder_Elements

end module SolidBody
