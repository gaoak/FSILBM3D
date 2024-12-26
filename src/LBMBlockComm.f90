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
    integer:: m_npairs,m_gridDelta=2! one divides intwo
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
        do i=1,m_npairs
            f = commpairs(i)%fatherId
            s = commpairs(i)%sonId
            if(LBMblks(s)%BndConds(1).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fx1t1(LBMblks(f)%zDim,LBMblks(f)%yDim,0:lbmDim))
                allocate(LBMblks(s)%fIn_Fx1t2(LBMblks(f)%zDim,LBMblks(f)%yDim,0:lbmDim))
                allocate(LBMblks(s)%uuu_Fx1t1(LBMblks(f)%zDim,LBMblks(f)%yDim,1:3     ))
            endif
            if(LBMblks(s)%BndConds(2).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fx2t1(LBMblks(f)%zDim,LBMblks(f)%yDim,0:lbmDim))
                allocate(LBMblks(s)%fIn_Fx2t2(LBMblks(f)%zDim,LBMblks(f)%yDim,0:lbmDim))
                allocate(LBMblks(s)%uuu_Fx2t1(LBMblks(f)%zDim,LBMblks(f)%yDim,1:3     ))
            endif
            if(LBMblks(s)%BndConds(3).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fy1t1(LBMblks(f)%zDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%fIn_Fy1t2(LBMblks(f)%zDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%uuu_Fy1t1(LBMblks(f)%zDim,LBMblks(f)%xDim,1:3     ))
            endif
            if(LBMblks(s)%BndConds(4).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fy2t1(LBMblks(f)%zDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%fIn_Fy2t2(LBMblks(f)%zDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%uuu_Fy2t1(LBMblks(f)%zDim,LBMblks(f)%xDim,1:3     ))
            endif
            if(LBMblks(s)%BndConds(5).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fz1t1(LBMblks(f)%yDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%fIn_Fz1t2(LBMblks(f)%yDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%uuu_Fz1t1(LBMblks(f)%yDim,LBMblks(f)%xDim,1:3     ))
            endif
            if(LBMblks(s)%BndConds(6).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fz2t1(LBMblks(f)%yDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%fIn_Fz2t2(LBMblks(f)%yDim,LBMblks(f)%xDim,0:lbmDim))
                allocate(LBMblks(s)%uuu_Fz2t1(LBMblks(f)%yDim,LBMblks(f)%xDim,1:3     ))
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
    endsubroutine bluid_block_tree

    recursive subroutine tree_collision_streaming_IBM_FEM(treenode,time_collision,time_streaming,time_IBM,time_FEM)
        use SolidBody, only: m_carrierFluidId
        implicit none
        integer:: i, s, treenode, n_timeStep
        real(8):: time_collision,time_streaming,time_IBM,time_FEM,time_begine2,time_end2
        ! calculate macro quantities for each blocks,must be ahead of collision(Huang Haibo 2024 P162)
        call LBMblks(treenode)%calculate_macro_quantities()
        call LBMblks(treenode)%update_volume_force(LBMblks(treenode)%blktime)
        call LBMblks(treenode)%ResetVolumeForce()
        if (blockTree(treenode)%nsons .eq. 0) then
            if (treenode .ne. m_carrierFluidId) then
                write(*,*) 'Error: solid is not inside in smallest block. smallest bolckID = ', treenode, 'solid in blockID = ', m_carrierFluidId
                stop
            else
                call IBM_FEM(time_IBM,time_FEM,LBMblks(m_carrierFluidId)%blktime)
            endif
        endif
        call LBMblks(treenode)%add_volume_force()
        ! extract interpolation layer for old time
        call extract_inner_layer(treenode,1)
        ! collision
        call get_now_time(time_begine2)
        call collision_block(treenode)
        call get_now_time(time_end2)
        time_collision = time_collision + (time_end2 - time_begine2)
        ! streaming
        call get_now_time(time_begine2)
        call streaming_block(treenode)
        call get_now_time(time_end2)
        time_streaming = time_streaming + (time_end2 - time_begine2)
        ! extract interpolation layer for new time
        call extract_inner_layer(treenode,2)
        !set boundary
        call set_boundary_conditions_block(treenode)
        ! tree cycle
        if(blockTree(treenode)%nsons.gt.0) then
            do i=1,blockTree(treenode)%nsons
                s = blockTree(treenode)%sons(i)
                do n_timeStep=0,m_gridDelta-1! one divides intwo
                    LBMblks(s)%blktime = LBMblks(s)%blktime + dble(n_timeStep) * LBMblks(s)%dh
                    call tree_collision_streaming_IBM_FEM(s,time_collision,time_streaming,time_IBM,time_FEM)
                    call interpolation_father_to_son(commpairs(blockTree(treenode)%pairId(i)),n_timeStep)
                enddo
                call deliver_son_to_father(commpairs(blockTree(treenode)%pairId(i)))
            enddo
        endif
    endsubroutine tree_collision_streaming_IBM_FEM

    subroutine IBM_FEM(time_IBM,time_FEM,time)
        use SolidBody, only: m_carrierFluidId,m_nFish,Solver,FSInteraction_force
        use FlowCondition, only: flow
        implicit none
        real(8):: time,dt_solid,time_IBM,time_FEM,time_begine2,time_end2
        integer:: isubstep=0
        dt_solid = LBMblks(m_carrierFluidId)%dh/dble(flow%numsubstep)       !time step of the solid
        call get_now_time(time_begine2)
        if (m_nFish .gt. 0) then
            call FSInteraction_force(LBMblks(m_carrierFluidId)%dh,LBMblks(m_carrierFluidId)%dh,LBMblks(m_carrierFluidId)%xmin,LBMblks(m_carrierFluidId)%ymin,LBMblks(m_carrierFluidId)%zmin, &
                                    LBMblks(m_carrierFluidId)%xDim,LBMblks(m_carrierFluidId)%yDim,LBMblks(m_carrierFluidId)%zDim,LBMblks(m_carrierFluidId)%uuu,LBMblks(m_carrierFluidId)%force)
        endif
        call get_now_time(time_end2)
        time_IBM = time_IBM + (time_end2 - time_begine2)
        call get_now_time(time_begine2)
        do isubstep=1,flow%numsubstep
            call Solver(time,isubstep,LBMblks(m_carrierFluidId)%dh,dt_solid)
        enddo !do isubstep=1,numsubstep
        call get_now_time(time_end2)
        time_FEM = time_FEM + (time_end2 - time_begine2)
    endsubroutine

    subroutine extract_inner_layer(treenode,time)
        integer:: treenode, time, f, s, i, ns
        integer:: xS, yS, zS, xF, yF, zF
        real(8):: xCoordSon, yCoordSon, zCoordSon
        f = treenode
        ns = blockTree(treenode)%nsons
        if(ns .ne. 0) then
            do i = 1,ns
                s = blockTree(treenode)%sons(i)
                if(LBMblks(s)%BndConds(1).eq.BCfluid) then
                    xS = 1
                    xCoordSon = LBMblks(s)%xmin + LBMblks(s)%dh*(xS - 1)
                    xF = NINT((xCoordSon - LBMblks(f)%xmin)/LBMblks(f)%dh + 1)
                    if (time .eq. 1) LBMblks(s)%fIn_Fx1t1(:,:,:) = LBMblks(f)%fIn(:,:,xF,:)
                    if (time .eq. 2) LBMblks(s)%fIn_Fx1t2(:,:,:) = LBMblks(f)%fIn(:,:,xF,:)
                    if (time .eq. 1) LBMblks(s)%uuu_Fx1t1(:,:,:) = LBMblks(f)%uuu(:,:,xF,:)
                endif
                if(LBMblks(s)%BndConds(2).eq.BCfluid) then
                    xS = LBMblks(s)%xDim
                    xCoordSon = LBMblks(s)%xmin + LBMblks(s)%dh*(xS - 1)
                    xF = NINT((xCoordSon - LBMblks(f)%xmin)/LBMblks(f)%dh + 1)
                    if (time .eq. 1) LBMblks(s)%fIn_Fx2t1(:,:,:) = LBMblks(f)%fIn(:,:,xF,:)
                    if (time .eq. 2) LBMblks(s)%fIn_Fx2t2(:,:,:) = LBMblks(f)%fIn(:,:,xF,:)
                    if (time .eq. 1) LBMblks(s)%uuu_Fx2t1(:,:,:) = LBMblks(f)%uuu(:,:,xF,:)
                endif
                if(LBMblks(s)%BndConds(3).eq.BCfluid) then
                    yS = 1
                    yCoordSon = LBMblks(s)%ymin + LBMblks(s)%dh*(yS - 1)
                    yF = NINT((yCoordSon - LBMblks(f)%ymin)/LBMblks(f)%dh + 1)
                    if (time .eq. 1) LBMblks(s)%fIn_Fy1t1(:,:,:) = LBMblks(f)%fIn(:,yF,:,:)
                    if (time .eq. 2) LBMblks(s)%fIn_Fy1t2(:,:,:) = LBMblks(f)%fIn(:,yF,:,:)
                    if (time .eq. 1) LBMblks(s)%uuu_Fy1t1(:,:,:) = LBMblks(f)%uuu(:,yF,:,:)
                endif
                if(LBMblks(s)%BndConds(4).eq.BCfluid) then
                    yS = LBMblks(s)%yDim
                    yCoordSon = LBMblks(s)%ymin + LBMblks(s)%dh*(yS - 1)
                    yF = NINT((yCoordSon - LBMblks(f)%ymin)/LBMblks(f)%dh + 1)
                    if (time .eq. 1) LBMblks(s)%fIn_Fy2t1(:,:,:) = LBMblks(f)%fIn(:,yF,:,:)
                    if (time .eq. 2) LBMblks(s)%fIn_Fy2t2(:,:,:) = LBMblks(f)%fIn(:,yF,:,:)
                    if (time .eq. 1) LBMblks(s)%uuu_Fy2t1(:,:,:) = LBMblks(f)%uuu(:,yF,:,:)
                endif
                if(LBMblks(s)%BndConds(5).eq.BCfluid) then
                    zS = 1
                    zCoordSon = LBMblks(s)%zmin + LBMblks(s)%dh*(zS - 1)
                    zF = NINT((zCoordSon - LBMblks(f)%zmin)/LBMblks(f)%dh + 1)
                    if (time .eq. 1) LBMblks(s)%fIn_Fz1t1(:,:,:) = LBMblks(f)%fIn(zF,:,:,:)
                    if (time .eq. 2) LBMblks(s)%fIn_Fz1t2(:,:,:) = LBMblks(f)%fIn(zF,:,:,:)
                    if (time .eq. 1) LBMblks(s)%uuu_Fz1t1(:,:,:) = LBMblks(f)%uuu(zF,:,:,:)
                endif
                if(LBMblks(s)%BndConds(6).eq.BCfluid) then
                    zS = LBMblks(s)%zDim
                    zCoordSon = LBMblks(s)%zmin + LBMblks(s)%dh*(zS - 1)
                    zF = NINT((zCoordSon - LBMblks(f)%zmin)/LBMblks(f)%dh + 1)
                    if (time .eq. 1) LBMblks(s)%fIn_Fz2t1(:,:,:) = LBMblks(f)%fIn(zF,:,:,:)
                    if (time .eq. 2) LBMblks(s)%fIn_Fz2t2(:,:,:) = LBMblks(f)%fIn(zF,:,:,:)
                    if (time .eq. 1) LBMblks(s)%uuu_Fz2t1(:,:,:) = LBMblks(f)%uuu(zF,:,:,:)
                endif
            enddo
        endif
    endsubroutine

    ! SUBROUTINE ExchangeFluidInterface()
    !     implicit none
    !     integer:: ip
    !     ! allocate and read commpairs from file
    !     do ip = 1,m_npairs
    !         call ExchangeDataSerial(commpairs(IP))
    !     enddo
    ! end subroutine ExchangeFluidInterface

    ! SUBROUTINE ExchangeDataSerial(pair)
    !     use FluidDomain
    !     implicit none
    !     type(CommPair),intent(in):: pair
    !     integer:: s(1:6),f(1:6),si(1:6),fi(1:6)
    !     ! x direction
    !     s = pair%s
    !     f = pair%f
    !     si = pair%si
    !     fi = pair%fi
    !     if(pair%sds(1).eq.1) then
    !         LBMblks(pair%sonId)%fIn(s(5):s(6),s(3):s(4),s(1),0:lbmDim) =&
    !             LBMblks(pair%fatherId)%fIn(f(5):f(6),f(3):f(4),f(1),0:lbmDim)
    !         LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(3):fi(4),fi(1),0:lbmDim) =&
    !             LBMblks(pair%sonId)%fIn(si(5):si(6),si(3):si(4),si(1),0:lbmDim)
    !     endif
    !     if(pair%sds(2).eq.-1) then
    !         LBMblks(pair%sonId)%fIn(s(5):s(6),s(3):s(4),s(2),0:lbmDim) =&
    !             LBMblks(pair%fatherId)%fIn(f(5):f(6),f(3):f(4),f(2),0:lbmDim)
    !         LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(3):fi(4),fi(2),0:lbmDim) =&
    !             LBMblks(pair%sonId)%fIn(si(5):si(6),si(3):si(4),si(2),0:lbmDim)
    !     endif
    !     ! y direction
    !     if(pair%sds(3).eq.1) then
    !         LBMblks(pair%sonId)%fIn(s(5):s(6),s(3),s(1):s(2),0:lbmDim) =&
    !             LBMblks(pair%fatherId)%fIn(f(5):f(6),f(3),f(1):f(2),0:lbmDim)
    !         LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(3),fi(1):fi(2),0:lbmDim) =&
    !             LBMblks(pair%sonId)%fIn(si(5):si(6),si(3),si(1):si(2),0:lbmDim)
    !     endif
    !     if(pair%sds(4).eq.-1) then
    !         LBMblks(pair%sonId)%fIn(s(5):s(6),s(4),s(1):s(2),0:lbmDim) =&
    !             LBMblks(pair%fatherId)%fIn(f(5):f(6),f(4),f(1):f(2),0:lbmDim)
    !         LBMblks(pair%fatherId)%fIn(fi(5):fi(6),fi(4),fi(1):fi(2),0:lbmDim) =&
    !             LBMblks(pair%sonId)%fIn(si(5):si(6),si(4),si(1):si(2),0:lbmDim)
    !     endif
    !     ! z direction
    !     if(pair%sds(5).eq.1) then
    !         LBMblks(pair%sonId)%fIn(s(5),s(3):s(4),s(1):s(2),0:lbmDim) =&
    !             LBMblks(pair%fatherId)%fIn(f(5),f(3):f(4),f(1):f(2),0:lbmDim)
    !         LBMblks(pair%fatherId)%fIn(fi(5),fi(3):fi(4),fi(1):fi(2),0:lbmDim) =&
    !             LBMblks(pair%sonId)%fIn(si(5),si(3):si(4),si(1):si(2),0:lbmDim)
    !     endif
    !     if(pair%sds(6).eq.-1) then
    !         LBMblks(pair%sonId)%fIn(s(6),s(3):s(4),s(1):s(2),0:lbmDim) =&
    !             LBMblks(pair%fatherId)%fIn(f(6),f(3):f(4),f(1):f(2),0:lbmDim)
    !         LBMblks(pair%fatherId)%fIn(fi(6),fi(3):fi(4),fi(1):fi(2),0:lbmDim) =&
    !             LBMblks(pair%sonId)%fIn(si(6),si(3):si(4),si(1):si(2),0:lbmDim)
    !     endif
    ! end subroutine ExchangeDataSerial
        
    ! SUBROUTINE calculate_blocks_tau(flow_nu)
    !     ! ensure the Reynolds numbers of each block are the same
    !     implicit none
    !     integer:: i,block
    !     real(8):: flow_nu
    !     if(m_npairs .eq. 0) then
    !         ! no need to set
    !     elseif(m_npairs .ge. 1) then
    !         do i=1,m_npairs
    !             LBMblks(commpairs(i)%sonId)%tau = 0.50d0 + dble(m_gridDelta) *(LBMblks(commpairs(i)%fatherId)%tau - 0.50d0)
    !         enddo
    !     else
    !         stop 'wrong block pairs (npairs) input'
    !     endif
    ! END SUBROUTINE

    ! verify the parameters of fluid blocks
    SUBROUTINE check_blocks_params(nblock)
        implicit none
        integer:: i,nblock
        logical:: flag
        flag = (1 .eq. 2) ! default flase
        if(nblock .ne. (m_npairs + 1)) stop 'the number of fluid blocks is wrong (should equal to n_pairs + 1) '
        if(m_npairs .ge. 1) then
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
        LBMblks(father)%uuu(zF,yF,xF,1:3     ) = LBMblks(son)%uuu(zS,yS,xS,1:3     )
    end subroutine

    SUBROUTINE interpolation_father_to_son(pair,n_timeStep)
        implicit none
        type(CommPair),intent(in):: pair
        integer:: n_timeStep,xS,yS,zS
        ! x direction
        if(pair%sds(1).eq.1) then
            do  yS=1,LBMblks(pair%sonId)%yDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                xS=1
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,n_timeStep, &
                                                    LBMblks(pair%fatherId)%yDim, LBMblks(pair%fatherId)%zDim, &
                                                    LBMblks(pair%sonId)%fIn_Fx1t1, LBMblks(pair%sonId)%fIn_Fx1t2, &
                                                    LBMblks(pair%sonId)%uuu_Fx1t1, 1)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        if(pair%sds(2).eq.-1) then
            do  yS=1,LBMblks(pair%sonId)%yDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                xS=  LBMblks(pair%sonId)%xDim
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,n_timeStep, &
                                                    LBMblks(pair%fatherId)%yDim, LBMblks(pair%fatherId)%zDim, &
                                                    LBMblks(pair%sonId)%fIn_Fx2t1, LBMblks(pair%sonId)%fIn_Fx2t2, &
                                                    LBMblks(pair%sonId)%uuu_Fx2t1, 2)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        ! y direction
        if(pair%sds(3).eq.1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                yS=1
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,n_timeStep, &
                                                    LBMblks(pair%fatherId)%xDim, LBMblks(pair%fatherId)%zDim, &
                                                    LBMblks(pair%sonId)%fIn_Fy1t1, LBMblks(pair%sonId)%fIn_Fy1t2, &
                                                    LBMblks(pair%sonId)%uuu_Fy1t1, 3)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        if(pair%sds(4).eq.-1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  zS=1,LBMblks(pair%sonId)%zDim
                yS=  LBMblks(pair%sonId)%yDim
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,n_timeStep, &
                                                    LBMblks(pair%fatherId)%xDim, LBMblks(pair%fatherId)%zDim, &
                                                    LBMblks(pair%sonId)%fIn_Fy2t1, LBMblks(pair%sonId)%fIn_Fy2t2, &
                                                    LBMblks(pair%sonId)%uuu_Fy2t1, 4)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        ! z direction
        if(pair%sds(5).eq.1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  yS=1,LBMblks(pair%sonId)%yDim
                zS=1
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,n_timeStep, &
                                                    LBMblks(pair%fatherId)%xDim, LBMblks(pair%fatherId)%yDim, &
                                                    LBMblks(pair%sonId)%fIn_Fz1t1, LBMblks(pair%sonId)%fIn_Fz1t2, &
                                                    LBMblks(pair%sonId)%uuu_Fz1t1, 5)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
        if(pair%sds(6).eq.-1) then
            do  xS=1,LBMblks(pair%sonId)%xDim
            do  yS=1,LBMblks(pair%sonId)%yDim
                zS=  LBMblks(pair%sonId)%zDim
                call interpolation_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,n_timeStep, &
                                                    LBMblks(pair%fatherId)%xDim, LBMblks(pair%fatherId)%yDim, &
                                                    LBMblks(pair%sonId)%fIn_Fz2t1, LBMblks(pair%sonId)%fIn_Fz2t2, &
                                                    LBMblks(pair%sonId)%uuu_Fz2t1, 6)
                call fIn_father_to_son(pair%fatherId,pair%sonId,xS,yS,zS)
            enddo
            enddo
        endif
    end subroutine

    SUBROUTINE interpolation_grid_distribution(father,son,xS,yS,zS,n_timeStep,aDim,bDim,fIn_t1,fIn_t2,uuu_t1,xyz)
        ! interpolating distribution function from father block to son block
        implicit none
        integer:: father,son,n_timeStep,xS,yS,zS,xyz
        integer:: xF1,yF1,zF1,xF2,yF2,zF2
        integer:: a1,a2,b1,b2,aDim,bDim
        real(8):: xCoordSon,yCoordSon,zCoordSon
        real(8):: dx1,dy1,dz1,dx2,dy2,dz2,da1,da2,db1,db2
        real(8):: coffe,c1,c2,c3,c4
        real(8):: fIn_S1(0:lbmDim),fIn_S2(0:lbmDim)
        real(8):: uuu_S_(1:3     )
        real(8):: fIn_t1(bDim,aDim,0:lbmDim),fIn_t2(bDim,aDim,0:lbmDim)
        real(8):: uuu_t1(bDim,aDim,1:3     )
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
        if (xF2 .gt. LBMblks(father)%xDim) xF2 = xF1
        if (yF2 .gt. LBMblks(father)%yDim) yF2 = yF1
        if (zF2 .gt. LBMblks(father)%zDim) zF2 = zF1
        if (xyz .eq. 1 .or. xyz .eq. 2) then
            a1 = yF1
            a2 = yF2
            b1 = zF1
            b2 = zF2
        endif
        if (xyz .eq. 3 .or. xyz .eq. 4) then
            a1 = xF1
            a2 = xF2
            b1 = zF1
            b2 = zF2
        endif
        if (xyz .eq. 5 .or. xyz .eq. 6) then
            a1 = xF1
            a2 = xF2
            b1 = yF1
            b2 = yF2
        endif
        ! coordinate difference between the son gird and around father grids
        dx1 = xCoordSon - (LBMblks(father)%xmin + LBMblks(father)%dh*(xF1 - 1))
        dy1 = yCoordSon - (LBMblks(father)%ymin + LBMblks(father)%dh*(yF1 - 1))
        dz1 = zCoordSon - (LBMblks(father)%zmin + LBMblks(father)%dh*(zF1 - 1))
        dx2 = LBMblks(father)%dh - dx1
        dy2 = LBMblks(father)%dh - dy1
        dz2 = LBMblks(father)%dh - dz1
        if (xyz .eq. 1 .or. xyz .eq. 2) then
            da1 = dy1
            da2 = dy2
            db1 = dz1
            db2 = dz2
        endif
        if (xyz .eq. 3 .or. xyz .eq. 4) then
            da1 = dx1
            da2 = dx2
            db1 = dz1
            db2 = dz2
        endif
        if (xyz .eq. 5 .or. xyz .eq. 6) then
            da1 = dx1
            da2 = dx2
            db1 = dy1
            db2 = dy2
        endif
        ! space interpolation coefficient
        coffe = 1.0d0/(LBMblks(father)%dh**2)
        c1 = da2*db2*coffe
        c2 = da2*db1*coffe
        c3 = da1*db2*coffe
        c4 = da1*db1*coffe
        ! space interpolation distribution function from father block to son block
        fIn_S1(0:lbmDim) = c1*fIn_t1(b1,a1,0:lbmDim) + c2*fIn_t1(b2,a1,0:lbmDim) + &
                           c3*fIn_t1(b1,a2,0:lbmDim) + c4*fIn_t1(b2,a2,0:lbmDim)
        fIn_S2(0:lbmDim) = c1*fIn_t2(b1,a1,0:lbmDim) + c2*fIn_t2(b2,a1,0:lbmDim) + &
                           c3*fIn_t2(b1,a2,0:lbmDim) + c4*fIn_t2(b2,a2,0:lbmDim)
        uuu_S_(1:3     ) = c1*uuu_t1(b1,a1,1:3     ) + c2*uuu_t1(b2,a1,1:3     ) + &
                           c3*uuu_t1(b1,a2,1:3     ) + c4*uuu_t1(b2,a2,1:3     )
        ! time interpolation distribution function from father block to son block
        if(n_timeStep==0)LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) = fIn_S1(0:lbmDim)
        if(n_timeStep==0)LBMblks(son)%uuu(zS,yS,xS,1:3     ) = uuu_S_(1:3     )
        if(n_timeStep==1)LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) = 0.5d0*fIn_S2(0:lbmDim)+0.5d0*fIn_S1(0:lbmDim)
        if(n_timeStep==1)LBMblks(son)%uuu(zS,yS,xS,1:3     ) = uuu_S_(1:3     )
        ! LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) = dble(n_timeStep - 0)/dble(m_gridDelta) * fIn_S2(0:lbmDim) + dble(m_gridDelta - n_timeStep)/dble(m_gridDelta) * fIn_S1(0:lbmDim)
    end subroutine

    subroutine fIn_father_to_son(father,son,xS,yS,zS)! Dupius-Chopard method
        implicit none
        integer:: father,son,xS,yS,zS
        real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim),coffe
        ! Guo 2008 P97 6.1.8
        uSqr           = sum(LBMblks(son)%uuu(zS,yS,xS,1:3)**2)
        uxyz(0:lbmDim) = LBMblks(son)%uuu(zS,yS,xS,1) * ee(0:lbmDim,1) + LBMblks(son)%uuu(zS,yS,xS,2) * ee(0:lbmDim,2)+LBMblks(son)%uuu(zS,yS,xS,3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * LBMblks(son)%den(zS,yS,xS) * ( (1.0d0 - 1.5d0 * uSqr) + uxyz(0:lbmDim) * (3.0d0  + 4.5d0 * uxyz(0:lbmDim)) )
        coffe = (LBMblks(son)%tau / LBMblks(father)%tau) / dble(m_gridDelta)
        LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) = fEq(0:lbmDim) + coffe * (LBMblks(son)%fIn(zS,yS,xS,0:lbmDim) - fEq(0:lbmDim))
    end subroutine

    subroutine fIn_son_to_father(father,son,xF,yF,zF)! Dupius-Chopard method
        implicit none
        integer:: father,son,xF,yF,zF
        real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim),coffe
        ! Guo 2008 P97 6.1.8
        uSqr           = sum(LBMblks(father)%uuu(zF,yF,xF,1:3)**2)
        uxyz(0:lbmDim) = LBMblks(father)%uuu(zF,yF,xF,1) * ee(0:lbmDim,1) + LBMblks(father)%uuu(zF,yF,xF,2) * ee(0:lbmDim,2)+LBMblks(father)%uuu(zF,yF,xF,3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * LBMblks(father)%den(zF,yF,xF) * ( (1.0d0 - 1.5d0 * uSqr) + uxyz(0:lbmDim) * (3.0d0  + 4.5d0 * uxyz(0:lbmDim)) )
        coffe = (LBMblks(father)%tau / LBMblks(son)%tau) * dble(m_gridDelta)
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
