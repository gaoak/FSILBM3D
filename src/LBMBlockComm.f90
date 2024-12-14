module LBMBlockComm
    use ConstParams
    use FluidDomain
    implicit none
    !include 'mpif.h'
    public:: Read_Comm_Pair,ExchangeFluidInterface
    type :: CommPair
        integer:: fatherId
        integer:: sonId
        integer:: sondirs(1:6) ! 0 no need communication; 1 need communication
        integer:: type ! local (0) or mpi (1)
    end type CommPair
    integer:: m_npairs
    type(CommPair),allocatable:: commpairs(:)

    contains

    SUBROUTINE Read_Comm_Pair(filename)
        implicit none
        character(LEN=40),intent(in):: filename
        character(LEN=40):: keywordstr
        character(LEN=256):: buffer
        integer:: i,j
        ! allocate and read commpairs from file
        open(unit=111, file=filename, status='old', action='read')
        keywordstr = 'Communication'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*) m_npairs
        allocate(commpairs(m_npairs))
        do i=1,m_npairs
            call readNextData(111, buffer)
            read(buffer,*) commpairs(i)%fatherId,commpairs(i)%sonId,commpairs(i)%type
            do j=1,6
                if(LBMblks(LBMblksIndex(commpairs(i)%sonId))%BndConds(j).eq.BCfluid) then
                    commpairs(i)%sondirs(j) = 1
                else
                    commpairs(i)%sondirs(j) = 0
                endif
            enddo
        enddo
        close(111)
    end subroutine Read_Comm_Pair

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
        integer:: fxS,fyS,fzS,fxE,fyE,fzE,sxD,syD,szD,f(6)
        sxD = LBMblks(pair%sonId)%xDim
        syD = LBMblks(pair%sonId)%yDim
        szD = LBMblks(pair%sonId)%zDim
        fxS = floor((LBMblks(pair%sonId)%xmin - LBMblks(pair%fatherId)%xmin) / LBMblks(pair%fatherId)%dh + 1.5d0)
        fyS = floor((LBMblks(pair%sonId)%ymin - LBMblks(pair%fatherId)%ymin) / LBMblks(pair%fatherId)%dh + 1.5d0)
        fzS = floor((LBMblks(pair%sonId)%zmin - LBMblks(pair%fatherId)%zmin) / LBMblks(pair%fatherId)%dh + 1.5d0)
        fxE = fxS + sxD - 1
        fyE = fyS + syD - 1
        fzE = fzS + szD - 1
        f = pair%sondirs
        ! x direction
        if(f(1).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1:syD,1,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS:fyE,fxS,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+f(6)):(fzE-f(5)),(fyS+f(4)):(fyE-f(3)),(fxS+1),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn((1+f(6)):(szD-f(5)),(1+f(4)):(syD-f(3)),2,0:lbmDim)
        endif
        if(f(2).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1:syD,sxD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS:fyE,fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+f(6)):(fzE-f(5)),(fyS+f(4)):(fyE-f(3)),(fxE-1),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn((1+f(6)):(szD-f(5)),(1+f(4)):(syD-f(3)),sxD-1,0:lbmDim)
        endif
        ! y direction
        if(f(3).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,1,1:sxD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyS,fxS:fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+f(6)):(fzE-f(5)),(fyS+1),(fxS+f(2)):(fxE-f(1)),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn((1+f(6)):(szD-f(5)),2,(1+f(2)):(sxD-f(1)),0:lbmDim)
        endif
        if(f(4).eq.1) then
            LBMblks(pair%sonId)%fIn(1:szD,syD,1:sxD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS:fzE,fyE,fxS:fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+f(6)):(fzE-f(5)),(fyE-1),(fxS+f(2)):(fxE-f(1)),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn((1+f(6)):(szD-f(5)),(syD-1),(1+f(2)):(sxD-f(1)),0:lbmDim)
        endif
        ! z direction
        if(f(5).eq.1) then
            LBMblks(pair%sonId)%fIn(1,1:syD,1:sxD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzS,fyS:fyE,fxS:fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzS+1),(fyS+f(4)):(fyE-f(3)),(fxS+f(2)):(fxE-f(1)),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn(2,(1+f(4)):(syD-f(3)),(1+f(2)):(sxD-f(1)),0:lbmDim)
        endif
        if(f(6).eq.1) then
            LBMblks(pair%sonId)%fIn(szD,1:syD,1:sxD,0:lbmDim) = LBMblks(pair%fatherId)%fIn(fzE,fyS:fyE,fxS:fxE,0:lbmDim)
            LBMblks(pair%fatherId)%fIn((fzE-1),(fyS+f(4)):(fyE-f(3)),(fxS+f(2)):(fxE-f(1)),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn((szD-1),(1+f(4)):(syD-f(3)),(1+f(2)):(sxD-f(1)),0:lbmDim)
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
