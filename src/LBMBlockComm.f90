module LBMBlockComm
    use ConstParams
    implicit none
    !include 'mpif.h'
    public:: ReadCommPair,ExchangeFluidInterface
    type :: CommPair
        integer:: fatherId
        integer:: sonId
        integer:: sondirs(1:6) ! 0 no need communication; 1 need communication
        integer:: type ! local or mpi
    end type CommPair
    integer:: m_npairs
    type(CommPair),allocatable:: commpairs(:)

    contains

    SUBROUTINE ReadCommPair(filename)
        implicit none
        character(LEN=40),intent(in):: filename
        ! allocate and read commpairs from file
        allocate(commpairs(m_npairs))
    end subroutine ReadCommPair

    SUBROUTINE ExchangeFluidInterface()
        implicit none
        integer:: ip
        ! allocate and read commpairs from file
        do ip = 1,m_npairs
            call ExchangeDataSerial(commpairs(ip))
        enddo
    end subroutine ExchangeFluidInterface

    SUBROUTINE ExchangeDataSerial(pair)
        use FluidDomain
        implicit none
        type(CommPair),intent(in):: pair
        integer:: fxS,fyS,fzS,fxE,fyE,fzE,sxD,syD,szD
        sxD = LBMblks(pair%sonId)%xDim
        syD = LBMblks(pair%sonId)%yDim
        szD = LBMblks(pair%sonId)%zDim
        fxS = floor((LBMblks(pair%sonId)%xmin - LBMblks(pair%fatherId)%xmin) / LBMblks(pair%fatherId)%dh + 1.5d0)
        fyS = floor((LBMblks(pair%sonId)%ymin - LBMblks(pair%fatherId)%ymin) / LBMblks(pair%fatherId)%dh + 1.5d0)
        fzS = floor((LBMblks(pair%sonId)%zmin - LBMblks(pair%fatherId)%zmin) / LBMblks(pair%fatherId)%dh + 1.5d0)
        fxE = fxS + sxD - 1
        fyE = fyS + syD - 1
        fzE = fzS + szD - 1
        ! x direction
        if(pair%sondirs(1).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1:syD,1,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS:fyE,fxS,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+1):(fzE-1),(fyS+1):(fyE-1),(fxS+1),0:lbmDim) = LBMblks(pair%sonId)%fIn(2:(szD-1),2:(syD-1),2,0:lbmDim)
        endif
        if(pair%sondirs(2).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1:syD,sxD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS:fyE,fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+1):(fzE-1),(fyS+1):(fyE-1),(fxE-1),0:lbmDim) = LBMblks(pair%sonId)%fIn(2:(szD-1),2:(syD-1),sxD-1,0:lbmDim)
        endif
        ! y direction
        if(pair%sondirs(3).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1,1:szD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS,fxS:fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+1):(fzE-1),(fyS+1),(fxS+1):(fxE-1),0:lbmDim) = LBMblks(pair%sonId)%fIn(2:(szD-1),2,2:(sxD-1),0:lbmDim)
        endif
        if(pair%sondirs(4).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,syD,1:szD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyE,fxS:fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+1):(fzE-1),(fyE-1),(fxS+1):(fxE-1),0:lbmDim) = LBMblks(pair%sonId)%fIn(2:(szD-1),(syD-1),2:(sxD-1),0:lbmDim)
        endif
        ! z direction
        if(pair%sondirs(5).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1:syD,1,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS:fyE,fxS,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+1):(fzE-1),(fyS+1):(fyE-1),(fxS+1),0:lbmDim) = LBMblks(pair%sonId)%fIn(2:(szD-1),2:(syD-1),2,0:lbmDim)
        endif
        if(pair%sondirs(6).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1:syD,szD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS:fyE,fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+1):(fzE-1),(fyE-1),(fxS+1):(fxE-1),0:lbmDim) = LBMblks(pair%sonId).fIn(2:(szD-1),(syD-1),2:(sxD-1),0:lbmDim)
        endif
    end subroutine ExchangeDataSerial
end module LBMBlockComm

! subroutine testmpiomp
!     use omp_lib
!     use LBMBlockComm
!     integer:: ipvd,ie,rank,ie2,size,np,ie3,x
!     integer:: ompn, ompid
!     call MPI_Init_thread(MPI_THREAD_FUNNELED,ipvd,ie)
!     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ie2)
!     call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ie3)
!     ompn = omp_get_num_threads()
!     write(*,*) ompn
!     !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,ompid,ompn)
!     do x=1,4
!         ompn = omp_get_num_threads()
!         ompid = omp_get_thread_num()
!         write(*,*)rank,np,ompn,ompid,x
!     enddo
!     !$OMP END PARALLEL DO
!     call MPI_Finalize(ie)
! end subroutine testmpiomp
