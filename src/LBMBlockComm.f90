module LBMBlockComm
    use ConstParams
    use FluidDomain
    implicit none
    !include 'mpif.h'
    public:: m_containSolidId
    public:: Read_Comm_Pair,ExchangeFluidInterface
    type :: CommPair
        integer:: fatherId
        integer:: sonId
        integer:: sds(1:6) ! 0 no need communication; 1(min),-1(max) need communication
        integer:: s(1:6),f(1:6),si(1:6),fi(1:6) ! s son's boundary layer; si son's first inner layer
        integer:: islocal ! local (0) or mpi (1)
    end type CommPair
    integer:: m_npairs,m_containSolidId
    type(CommPair),allocatable:: commpairs(:)

    contains

    SUBROUTINE Read_Comm_Pair(filename)
        implicit none
        character(LEN=40),intent(in):: filename
        character(LEN=40):: keywordstr
        character(LEN=256):: buffer
        integer:: i,j,sxD,syD,szD,ratio,fId,sId,islocal
        ! allocate and read commpairs from file
        open(unit=111, file=filename, status='old', action='read')
        keywordstr = 'Communication'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*) m_npairs,m_containSolidId
        allocate(commpairs(m_npairs))
        do i=1,m_npairs
            call readNextData(111, buffer)
            read(buffer,*) fId,sId,islocal
            commpairs(i)%fatherId = fId
            commpairs(i)%sonId = sId
            commpairs(i)%islocal = islocal
            do j=1,6
                if(LBMblks(LBMblksIndex(sId))%BndConds(j).eq.BCfluid) then
                    if(mod(j,2).eq.1) then
                        commpairs(i)%sds(j) = 1
                    else
                        commpairs(i)%sds(j) = -1
                    endif
                else
                    commpairs(i)%sds(j) = 0
                endif
            enddo
            sxD = LBMblks(sId)%xDim
            syD = LBMblks(sId)%yDim
            szD = LBMblks(sId)%zDim
            commpairs(i)%s([1,3,5]) = 1
            commpairs(i)%s(2) = sxD
            commpairs(i)%s(4) = syD
            commpairs(i)%s(6) = szD
            commpairs(i)%f(1) = floor((LBMblks(sId)%xmin - LBMblks(fId)%xmin) / LBMblks(fId)%dh + 1.5d0)
            commpairs(i)%f(3) = floor((LBMblks(sId)%ymin - LBMblks(fId)%ymin) / LBMblks(fId)%dh + 1.5d0)
            commpairs(i)%f(5) = floor((LBMblks(sId)%zmin - LBMblks(fId)%zmin) / LBMblks(fId)%dh + 1.5d0)
            ratio = floor(LBMblks(fId)%dh / LBMblks(sId)%dh + 0.5d0)
            sxD = (sxD - 1) / ratio
            syD = (syD - 1) / ratio
            szD = (szD - 1) / ratio
            commpairs(i)%f(2) = commpairs(i)%f(1) + sxD
            commpairs(i)%f(4) = commpairs(i)%f(3) + syD
            commpairs(i)%f(6) = commpairs(i)%f(5) + szD
            commpairs(i)%si(1:6) = commpairs(i)%s(1:6) + commpairs(i)%sds(1:6)
            commpairs(i)%fi(1:6) = commpairs(i)%f(1:6) + commpairs(i)%sds(1:6)
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
        integer:: s(1:6),f(1:6),si(1:6),fi(1:6)
        ! x direction
        s = pair%s
        f = pair%f
        si = pair%si
        fi = pair%fi
        if(pair%sds(1).eq.1) then
            LBMblks(pair%sonId)%fIn(s(5):s(6),s(3):s(4),s(1),0:lbmDim) =&
                LBMblks(pair%fatherId)%fIn(f(5):f(6),f(3):f(4),f(1),0:lbmDim)
            LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(3):fi(4),fi(1),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn(si(5):si(6),si(3):si(4),si(1),0:lbmDim)
        endif
        if(pair%sds(2).eq.-1) then
            LBMblks(pair%sonId)%fIn(s(5):s(6),s(3):s(4),s(2),0:lbmDim) =&
                LBMblks(pair%fatherId)%fIn(f(5):f(6),f(3):f(4),f(2),0:lbmDim)
            LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(3):fi(4),fi(2),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn(si(5):si(6),si(3):si(4),si(2),0:lbmDim)
        endif
        ! y direction
        if(pair%sds(3).eq.1) then
            LBMblks(pair%sonId)%fIn(s(5):s(6),s(3),s(1):s(2),0:lbmDim) =&
                LBMblks(pair%fatherId)%fIn(f(5):f(6),f(3),f(1):f(2),0:lbmDim)
            LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(3),fi(1):fi(2),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn(si(5):si(6),si(3),si(1):si(2),0:lbmDim)
        endif
        if(pair%sds(4).eq.-1) then
            LBMblks(pair%sonId)%fIn(s(5):s(6),s(4),s(1):s(2),0:lbmDim) =&
                LBMblks(pair%fatherId)%fIn(f(5):f(6),f(4),f(1):f(2),0:lbmDim)
            LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(4),fi(1):fi(2),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn(si(5):si(6),si(4),si(1):si(2),0:lbmDim)
        endif
        ! z direction
        if(pair%sds(5).eq.1) then
            LBMblks(pair%sonId)%fIn(s(5),s(3):s(4),s(1):s(2),0:lbmDim) =&
                LBMblks(pair%fatherId)%fIn(f(5),f(3):f(4),f(1):f(2),0:lbmDim)
            LBMblks(pair%fatherId)%fIn(fi(5),fi(3):fi(4),fi(1):fi(2),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn(si(5),si(3):si(4),si(1):si(2),0:lbmDim)
        endif
        if(pair%sds(6).eq.-1) then
            LBMblks(pair%sonId)%fIn(s(6),s(3):s(4),s(1):s(2),0:lbmDim) =&
                LBMblks(pair%fatherId)%fIn(f(6),f(3):f(4),f(1):f(2),0:lbmDim)
            LBMblks(pair%fatherId)%fIn(fi(6),fi(3):fi(4),fi(1):fi(2),0:lbmDim) =&
                LBMblks(pair%sonId)%fIn(si(6),si(3):si(4),si(1):si(2),0:lbmDim)
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
