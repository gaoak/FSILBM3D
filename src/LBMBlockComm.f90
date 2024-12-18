module LBMBlockComm
    use ConstParams
    use FluidDomain
    implicit none
    !include 'mpif.h'
    type :: blockTreeNode
        integer:: fatherId
        integer:: nsons
        integer:: nt
        integer,allocatable:: sons(:)
        integer,allocatable:: pairId(:)
    end type blockTreeNode
    type :: CommPair
        integer:: fatherId
        integer:: sonId
        integer:: sds(1:6) ! 0 no need communication; 1(min),-1(max) need communication
        integer:: s(1:6),f(1:6),si(1:6),fi(1:6) ! s son's boundary layer; si son's first inner layer
        integer:: islocal ! local (0) or mpi (1)
    end type CommPair
    integer:: m_npairs,m_gridDelta
    integer:: blockTreeRoot ! assume there is only one root
    type(blockTreeNode),allocatable:: blockTree(:)
    type(CommPair),allocatable:: commpairs(:)

    contains

    SUBROUTINE read_blocks_comunication(filename)
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
        read(buffer,*) m_npairs
        if(m_npairs .ne. 0) then
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
        endif
        close(111)
    end subroutine read_blocks_comunication

    subroutine bluid_block_tree()
        implicit none
        integer:: i,f,s,ns(1:m_nblock)
        ! initialise blocktree
        allocate(blockTree(1:m_nblock))
        do i=1,m_nblock
            blockTree(i)%fatherId = 0
            blockTree(i)%nsons = 0
        enddo
        ! get number of sons and father Id
        do i=1,m_npairs
            f = commpairs(i)%fatherId
            s = commpairs(i)%sonId
            blockTree(f)%nsons = blockTree(f)%nsons + 1
            blockTree(s)%fatherId = f
        enddo
        ! allocate sons and determine tree root
        blockTreeRoot = -1111
        do i=1,m_nblock
            allocate(blockTree(i)%sons(blockTree(i)%nsons),blockTree(i)%pairId(blockTree(i)%nsons))
            if(blockTree(i)%fatherId .eq. 0) then
                if(blockTreeRoot.ne.-1111) then
                    write(*,*) 'Error: there are more than one root block', i, blockTreeRoot
                    stop
                endif
                blockTreeRoot = i
            endif
        enddo
        ! set sons
        ns = 0
        do i=1,m_npairs
            f = commpairs(i)%fatherId
            s = commpairs(i)%sonId
            ns(f) = ns(f) + 1
            blockTree(f)%sons(ns(f)) = s
            blockTree(f)%pairId(ns(f)) = i
        enddo
        ! set tree n times
        call set_block_tree_nt(blockTreeRoot, 1)
    endsubroutine bluid_block_tree

    recursive subroutine set_block_tree_nt(treenode, nt)
        implicit none
        integer:: i, n, s, treenode, nt
        if(blockTree(treenode)%nsons.eq.0) then
            blockTree(treenode)%nt = nt
        else
            do i=1,blockTree(treenode)%nsons
                s = blockTree(treenode)%sons(i)
                if(dabs(LBMblks(treenode)%dh / LBMblks(s)%dh - 1.) .lt. 1d-6) then
                    call set_block_tree_nt(s, nt)
                    write(*,*) 'Warning: imbeded grid with the same grid size is not suggested'
                else
                    call set_block_tree_nt(s, 2*nt)
                endif
            enddo
        endif
    endsubroutine set_block_tree_nt

    recursive subroutine tree_collision_streaming(treenode)
        implicit none
        integer:: i, n, s, treenode, nt, n_gridDelta
        real(8):: fIn_F1(LBMblks(treenode)%zDim,LBMblks(treenode)%yDim,LBMblks(treenode)%xDim,0:lbmDim) !   n   time step father fIn
        real(8):: fIn_F2(LBMblks(treenode)%zDim,LBMblks(treenode)%yDim,LBMblks(treenode)%xDim,0:lbmDim) ! n + 1 time step father fIn
            if(blockTree(treenode)%nsons.eq.0) then
                call collision_block(treenode)
                call streaming_block(treenode)
            else
                do i=1,blockTree(treenode)%nsons
                    s = blockTree(treenode)%sons(i)
                    call deliver_son_to_father(commpairs(i)) ! to be check
                    fIn_F1 = LBMblks(treenode)%fIn
                    call collision_block(treenode)
                    fIn_F2 = LBMblks(treenode)%fIn
                    call streaming_block(treenode)
                    do n_gridDelta=0,m_gridDelta-1
                        call interpolation_father_to_son(commpairs(i),fIn_F1,fIn_F2,n_gridDelta)
                        call tree_collision_streaming(s)
                    enddo
                enddo
            endif
    endsubroutine tree_collision_streaming

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
        
    SUBROUTINE calculate_blocks_tau(flow_nu)
        ! ensure the Reynolds numbers of each block are the same
        implicit none
        integer:: i,nblock
        real(8):: flow_nu
        if(m_npairs .eq. 0) then
            ! no need to set
        elseif(m_npairs .ge. 1) then
            LBMblks(commpairs(1)%sonId)%tau = 0.50d0 + m_gridDelta*(LBMblks(commpairs(1)%fatherId)%tau - 0.50d0)
            if(m_npairs .ge. 2) then
                do i=2,m_npairs
                    LBMblks(commpairs(i)%sonId)%tau = 0.50d0 + m_gridDelta*(LBMblks(commpairs(i)%fatherId)%tau - 0.50d0)
                enddo
            endif
        else
            stop 'wrong block pairs (npairs) input'
        endif
    END SUBROUTINE

    ! verify the parameters of fluid blocks
    SUBROUTINE check_blocks_params(nblock)
        implicit none
        integer:: i,nblock
        logical:: flag
        flag = (1 .eq. 2) ! default flase
        if(nblock .ne. (m_npairs + 1)) stop 'the number of fluid blocks is wrong (should equal to n_pairs + 1) '
        if(m_npairs .ge. 1 .and. m_gridDelta .ne. 1) then
            do i=1,m_npairs
                flag = (LBMblks(commpairs(i)%fatherId)%dh .ne. LBMblks(commpairs(i)%sonId)%dh*m_gridDelta .or. mod(LBMblks(commpairs(i)%sonId)%xDim,m_gridDelta) .ne. 1 .or. &
                        mod(LBMblks(commpairs(i)%sonId)%yDim,m_gridDelta) .ne. 1 .or. mod(LBMblks(commpairs(i)%sonId)%zDim,m_gridDelta) .ne. 1)
                if(flag) stop 'grid points do not match between fluid blocks'
            enddo
        endif
    END SUBROUTINE

    SUBROUTINE deliver_son_to_father(pair)
        use FluidDomain
        implicit none
        type(CommPair),intent(in):: pair
        integer:: xS,yS,zS,zF,yF,xF
        ! x direction slices
        if(pair%sds(1).eq.1) then
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                xS=m_gridDelta+1
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,zF,yF,xF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
        endif
        if(pair%sds(2).eq.-1) then
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                xS=LBMblks(pair%sonId)%xDim - m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,zF,yF,xF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
        endif
        ! y direction slices
        if(pair%sds(3).eq.1) then
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                yS=m_gridDelta+1
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,zF,yF,xF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
        endif
        if(pair%sds(4).eq.-1) then
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                yS=LBMblks(pair%sonId)%yDim - m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,zF,yF,xF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
        endif
        ! z direction slices
        if(pair%sds(5).eq.1) then
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
                zS=m_gridDelta+1
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,zF,yF,xF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
        endif
        if(pair%sds(6).eq.-1) then
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
                zS=LBMblks(pair%sonId)%zDim - m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,zF,yF,xF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
        endif
    end subroutine 

    SUBROUTINE deliver_grid_distribution(father,son,xS,yS,zS,zF,yF,xF)
        implicit none
        integer:: father,son
        integer:: xS,yS,zS,xF,yF,zF
        real(8):: xCoordSon,yCoordSon,zCoordSon
        ! calcualte the pubilc grid coordinates in son block
        xCoordSon = LBMblks(son)%xmin + LBMblks(son)%dh*(xS - 1)
        yCoordSon = LBMblks(son)%ymin + LBMblks(son)%dh*(yS - 1)
        zCoordSon = LBMblks(son)%zmin + LBMblks(son)%dh*(zS - 1)
        ! calcualte the pubilc grid number in father block
        xF = NINT((xCoordSon - LBMblks(father)%xmin)/LBMblks(father)%dh + 1)
        yF = NINT((yCoordSon - LBMblks(father)%ymin)/LBMblks(father)%dh + 1)
        zF = NINT((zCoordSon - LBMblks(father)%zmin)/LBMblks(father)%dh + 1)
        ! deliver distribution function from son block to father block
        LBMblks(father)%fIn(zF,yF,xF,0:lbmDim) = LBMblks(son)%fIn(zS,yS,xS,0:lbmDim)
    end subroutine

    SUBROUTINE interpolation_father_to_son(pair,fIn_F1,fIn_F2,n_timeStep)
        implicit none
        type(CommPair),intent(in):: pair
        integer:: n_timeStep,xS,yS,zS
        real(8):: fIn_F1(LBMblks(pair%fatherId)%zDim,LBMblks(pair%fatherId)%yDim,LBMblks(pair%fatherId)%xDim,0:lbmDim) !   n   time step father fIn
        real(8):: fIn_F2(LBMblks(pair%fatherId)%zDim,LBMblks(pair%fatherId)%yDim,LBMblks(pair%fatherId)%xDim,0:lbmDim) ! n + 1 time step father fIn
        ! x direction
        if(pair%sds(1).eq.1) then
            do  yS=1,LBMblks(pair%sonId)%yDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                xS=1
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,fIn_F1,fIn_F2,n_timeStep)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        if(pair%sds(2).eq.-1) then
            do  yS=1,LBMblks(pair%sonId)%yDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                xS=  LBMblks(pair%sonId)%xDim
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,fIn_F1,fIn_F2,n_timeStep)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        ! y direction
        if(pair%sds(3).eq.1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                yS=1
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,fIn_F1,fIn_F2,n_timeStep)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        if(pair%sds(4).eq.-1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                yS=  LBMblks(pair%sonId)%yDim
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,fIn_F1,fIn_F2,n_timeStep)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        ! z direction
        if(pair%sds(5).eq.1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  yS=1,LBMblks(pair%sonId)%yDim
                zS=1
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,fIn_F1,fIn_F2,n_timeStep)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        if(pair%sds(6).eq.-1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  yS=1,LBMblks(pair%sonId)%yDim
                zS=  LBMblks(pair%sonId)%zDim
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,fIn_F1,fIn_F2,n_timeStep)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
    end subroutine

    SUBROUTINE interpolation_grid_distribution(father,son,xS,yS,zS,fIn_F1,fIn_F2,n_timeStep)
        ! interpolating distribution function from father block to son block
        implicit none
        integer:: father,son,n_timeStep,xS,yS,zS
        integer:: xF1,yF1,zF1,xF2,yF2,zF2
        real(8):: xCoordSon,yCoordSon,zCoordSon
        real(8):: dx1,dy1,dz1,dx2,dy2,dz2
        real(8):: coffe,c1,c2,c3,c4,c5,c6,c7,c8
        real(8):: fIn_S1(0:lbmDim),fIn_S2(0:lbmDim)
        real(8):: fIn_F1(LBMblks(father)%zDim,LBMblks(father)%yDim,LBMblks(father)%xDim,0:lbmDim) !   n   time step father fIn
        real(8):: fIn_F2(LBMblks(father)%zDim,LBMblks(father)%yDim,LBMblks(father)%xDim,0:lbmDim) ! n + 1 time step father fIn
        ! calcualte the pubilc grid coordinates in son block
        xCoordSon = LBMblks(son)%xmin + LBMblks(son)%dh*(xS - 1)
        yCoordSon = LBMblks(son)%ymin + LBMblks(son)%dh*(yS - 1)
        zCoordSon = LBMblks(son)%zmin + LBMblks(son)%dh*(zS - 1)
        ! calcualte the father grids around the son gird
        xF1 = FLOOR((xCoordSon - LBMblks(father)%xmin)/LBMblks(father)%dh + 1)
        yF1 = FLOOR((yCoordSon - LBMblks(father)%ymin)/LBMblks(father)%dh + 1)
        zF1 = FLOOR((zCoordSon - LBMblks(father)%zmin)/LBMblks(father)%dh + 1)
        xF2 = xF1 + 1
        yF2 = yF1 + 1
        zF2 = zF1 + 1
        ! coordinate difference between the son gird and around father grids
        dx1 = xCoordSon - (LBMblks(father)%xmin + LBMblks(father)%dh*(xF1 - 1))
        dy1 = yCoordSon - (LBMblks(father)%ymin + LBMblks(father)%dh*(yF1 - 1))
        dz1 = zCoordSon - (LBMblks(father)%zmin + LBMblks(father)%dh*(zF1 - 1))
        dx2 = LBMblks(father)%dh - dx1
        dy2 = LBMblks(father)%dh - dy1
        dz2 = LBMblks(father)%dh - dz1
        ! space interpolation coefficient
        coffe = 1.0d0/(LBMblks(father)%dh**3)
        c1 = dx2*dy2*dz2*coffe
        c2 = dx2*dy2*dz1*coffe
        c3 = dx2*dy1*dz2*coffe
        c4 = dx1*dy2*dz2*coffe
        c5 = dx2*dy1*dz1*coffe
        c6 = dx1*dy2*dz1*coffe
        c7 = dx1*dy1*dz2*coffe
        c8 = dx1*dy1*dz1*coffe
        ! space interpolation distribution function from father block to son block
        fIn_S1(0:lbmDim) = c1*fIn_F1(zF1,yF1,xF1,0:lbmDim) + c2*fIn_F1(zF2,yF1,xF1,0:lbmDim) + &
                           c3*fIn_F1(zF1,yF2,xF1,0:lbmDim) + c4*fIn_F1(zF1,yF1,xF2,0:lbmDim) + &
                           c5*fIn_F1(zF2,yF2,xF1,0:lbmDim) + c6*fIn_F1(zF2,yF1,xF2,0:lbmDim) + &
                           c7*fIn_F1(zF1,yF2,xF2,0:lbmDim) + c8*fIn_F1(zF2,yF2,xF2,0:lbmDim)
        fIn_S2(0:lbmDim) = c1*fIn_F2(zF1,yF1,xF1,0:lbmDim) + c2*fIn_F2(zF2,yF1,xF1,0:lbmDim) + &
                           c3*fIn_F2(zF1,yF2,xF1,0:lbmDim) + c4*fIn_F2(zF1,yF1,xF2,0:lbmDim) + &
                           c5*fIn_F2(zF2,yF2,xF1,0:lbmDim) + c6*fIn_F2(zF2,yF1,xF2,0:lbmDim) + &
                           c7*fIn_F2(zF1,yF2,xF2,0:lbmDim) + c8*fIn_F2(zF2,yF2,xF2,0:lbmDim)
        ! time interpolation distribution function from father block to son block
        LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) = (n_timeStep - 0)/m_gridDelta * fIn_S2(0:lbmDim) + (m_gridDelta - n_timeStep)/m_gridDelta * fIn_S1(0:lbmDim)
    end subroutine

    subroutine fIn_father_to_son(father,son,xS,yS,zS)! Dupius-Chopard method
        implicit none
        integer:: father,son,xS,yS,zS
        real(8):: fEq(0:lbmDim),coffe
        ! Guo 2008 P97 6.1.8
        call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
        coffe = (LBMblks(son)%tau / LBMblks(father)%tau) / m_gridDelta
        LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) = fEq(0:lbmDim) + coffe * (LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) - fEq(0:lbmDim))
    end subroutine

    subroutine fIn_son_to_father(father,son,xF,yF,zF)! Dupius-Chopard method
        implicit none
        integer:: father,son,xF,yF,zF
        real(8):: fEq(0:lbmDim),coffe
        ! Guo 2008 P97 6.1.8
        call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
        coffe = (LBMblks(father)%tau / LBMblks(son)%tau) * m_gridDelta
        LBMblks(father)%fIn(zF,yF,xF,0:lbmDim) = fEq(0:lbmDim) + coffe * (LBMblks(father)%fIn(zF,yF,xF,0:lbmDim) - fEq(0:lbmDim))
    end subroutine
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
