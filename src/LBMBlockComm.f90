module LBMBlockComm
    use ConstParams
    use FluidDomain
    implicit none
    !include 'mpif.h'
    type :: CommPair
        integer:: fatherId
        integer:: sonId
        integer:: sds(1:6) ! 0 no need communication; 1(min),-1(max) need communication
        integer:: s(1:6),f(1:6),si(1:6),fi(1:6) ! s son's boundary layer; si son's first inner layer
        integer:: islocal ! local (0) or mpi (1)
    end type CommPair
    type :: blockTreeNode
        integer:: fatherId
        integer:: nsons
        integer:: nt
        integer,allocatable:: sons(:)
        type(CommPair),allocatable:: comm(:)
    end type blockTreeNode
    integer:: m_npairs,m_gridDelta=2! one divides intwo
    integer:: blockTreeRoot ! assume there is only one root
    type(blockTreeNode),allocatable:: blockTree(:)

    contains

    recursive SUBROUTINE build_blocks_comunication(treenode)
        implicit none
        integer,intent(in):: treenode
        integer:: nsons,i,j,sxD,syD,szD,ratio,fId,sId,islocal
        type(CommPair)::pair
        nsons = blockTree(treenode)%nsons
        if(nsons .eq. 0) return
        allocate(blockTree(treenode)%comm(nsons))
        do i=1,nsons
            pair%fatherId = treenode
            pair%sonId = blockTree(treenode)%sons(i)
            fId = pair%fatherId
            sId = pair%sonId
            pair%islocal = 0
            do j=1,6
                if(LBMblks(blockTree(fId)%sons(i))%BndConds(j).eq.BCfluid) then
                    if(mod(j,2).eq.1) then
                        pair%sds(j) = 1
                    else
                        pair%sds(j) = -1
                    endif
                else
                    pair%sds(j) = 0
                endif
                sxD = LBMblks(sId)%xDim
                syD = LBMblks(sId)%yDim
                szD = LBMblks(sId)%zDim
                pair%s([1,3,5]) = 1
                pair%s(2) = sxD
                pair%s(4) = syD
                pair%s(6) = szD
                pair%f(1) = floor((LBMblks(sId)%xmin - LBMblks(fId)%xmin) / LBMblks(fId)%dh + 1.5d0)
                pair%f(3) = floor((LBMblks(sId)%ymin - LBMblks(fId)%ymin) / LBMblks(fId)%dh + 1.5d0)
                pair%f(5) = floor((LBMblks(sId)%zmin - LBMblks(fId)%zmin) / LBMblks(fId)%dh + 1.5d0)
                ratio = floor(LBMblks(fId)%dh / LBMblks(sId)%dh + 0.5d0)
                sxD = (sxD - 1) / ratio
                syD = (syD - 1) / ratio
                szD = (szD - 1) / ratio
                pair%f(2) = pair%f(1) + sxD
                pair%f(4) = pair%f(3) + syD
                pair%f(6) = pair%f(5) + szD
                pair%si(1:6) = pair%s(1:6) + pair%sds(1:6) * ratio
                pair%fi(1:6) = pair%f(1:6) + pair%sds(1:6)
                blockTree(treenode)%comm(i) = pair
            enddo
            call build_blocks_comunication(blockTree(treenode)%sons(i))
        enddo
    end subroutine build_blocks_comunication

    subroutine findremove_blockTreeRoot(iblocks,nb,rootnode)
        implicit none
        integer, intent(inout):: nb,iblocks(1:nb),rootnode
        integer:: fa(1:nb)
        integer:: i, j, tmp,cr,cb
        if(nb.eq.0) return
        ! compare all blocks
        fa = -1
        do i=1,nb-1
            do j=i+1,nb
                tmp = CompareBlocks(iblocks(i), iblocks(j))
                if(tmp.eq.1) then
                    fa(j) = iblocks(i)
                elseif(tmp.eq.-1) then
                    fa(i) = iblocks(j)
                endif
            enddo
        enddo
        ! find all roots
        cr = 0
        cb = 0
        do i=1,nb
            if(fa(i).lt.0) then
                cr = cr + 1
                rootnode = iblocks(i)
            else
                cb = cb + 1
                iblocks(cb) = iblocks(i)
            endif
        enddo
        nb = nb - cr
        if(cr.gt.1) then
            write(*,*) 'Error: there exist more than one block tree root', rootnode
            stop
        endif
    end subroutine

    recursive subroutine array_to_tree(iblocks, nb, rootnode)
        implicit none
        integer, intent(in):: nb, iblocks(1:nb),rootnode
        integer:: nr, roots(1:nb),fa(1:nb),subblock(1:nb)
        integer:: i, j, tmp,cr,cb
        if(nb.eq.0) return
        ! compare all blocks
        fa = -1
        do i=1,nb-1
            do j=i+1,nb
                tmp = CompareBlocks(iblocks(i), iblocks(j))
                if(tmp.eq.1) then
                    fa(j) = iblocks(i)
                elseif(tmp.eq.-1) then
                    fa(i) = iblocks(j)
                endif
            enddo
        enddo
        tmp = 1
        do while(tmp.gt.0)
            tmp = 0
            do i=1,nb
                if(fa(i).gt.0) then
                    if(fa(fa(i)).gt.0 .and. fa(i).ne.fa(fa(i))) then
                        fa(i) = fa(fa(i))
                        tmp = 1
                    endif
                endif
            enddo
        enddo
        ! find all roots
        cr = 0
        do i=1,nb
            if(fa(i).lt.0) then
                cr = cr + 1
                roots(cr) = iblocks(i)
            endif
        enddo
        ! build one layer tree for root node
        blockTree(rootnode)%nsons = cr
        allocate(blockTree(rootnode)%sons(cr))
        do i=1,cr
            blockTree(rootnode)%sons(i) = roots(i)
            blockTree(roots(i))%fatherId = rootnode
        enddo
        ! devide works
        do i=1,cr
            cb = 0
            do j=1,nb
                if(fa(j).eq.roots(i)) then
                    cb = cb + 1
                    subblock(cb) = iblocks(j)
                endif
            enddo
            call array_to_tree(subblock, cb, roots(i))
        enddo
    end subroutine

    subroutine bluid_block_tree()
        implicit none
        integer:: i,f,s,ns(1:m_nblocks),iblocks(1:m_nblocks),nb
        ! initialise blocktree
        allocate(blockTree(1:m_nblocks))
        do i=1,m_nblocks
            blockTree(i)%fatherId = 0
            blockTree(i)%nsons = 0
            iblocks(i) = i
        enddo
        nb = m_nblocks
        call findremove_blockTreeRoot(iblocks,nb,blockTreeRoot)
        call array_to_tree(iblocks,nb,blockTreeRoot)
        call build_blocks_comunication(blockTreeRoot)
        ! allocate father slices in sons node
        call allocate_fIn_uuu(blockTreeRoot)
    endsubroutine bluid_block_tree

    recursive subroutine allocate_fIn_uuu(treenode)
        implicit none
        integer,intent(in):: treenode
        integer:: i,f,s
        do i=1,blockTree(treenode)%nsons
            f = blockTree(treenode)%comm(i)%fatherId
            s = blockTree(treenode)%comm(i)%sonId
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
            call allocate_fIn_uuu(blockTree(treenode)%sons(i))
        enddo
    endsubroutine

    recursive subroutine tree_collision_streaming_IBM_FEM(treenode,time_collision,time_streaming,time_IBM,time_FEM)
        use SolidBody, only: m_nFish,VBodies
        implicit none
        integer:: i, s, treenode, n_timeStep, iFish
        real(8):: time_collision,time_streaming,time_IBM,time_FEM,time_begine2,time_end2
        call LBMblks(treenode)%update_volume_force()
        ! calculate macro quantities for each blocks,must be ahead of collision(Huang Haibo 2024 P162)
        call LBMblks(treenode)%calculate_macro_quantities()
        call LBMblks(treenode)%ResetVolumeForce()
        if (blockTree(treenode)%nsons .eq. 0) then
            do iFish = 1,m_nFish
                if (treenode .eq. VBodies(iFish)%v_carrierFluidId) then
                    call IBM_FEM(treenode,time_IBM,time_FEM,LBMblks(treenode)%blktime)
                endif
            enddo
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
                    call interpolation_father_to_son(blockTree(treenode)%comm(i),n_timeStep)
                enddo
                call deliver_son_to_father(blockTree(treenode)%comm(i))
            enddo
        endif
    endsubroutine tree_collision_streaming_IBM_FEM

    subroutine IBM_FEM(carrierFluidId,time_IBM,time_FEM,time)
        use SolidBody, only: Solver,FSInteraction_force
        use FlowCondition, only: flow
        implicit none
        real(8):: time,dt_solid,time_IBM,time_FEM,time_begine2,time_end2
        integer:: carrierFluidId,isubstep=0
        dt_solid = LBMblks(carrierFluidId)%dh/dble(flow%numsubstep)       !time step of the solid
        call get_now_time(time_begine2)
        call FSInteraction_force(LBMblks(carrierFluidId)%dh,LBMblks(carrierFluidId)%dh,LBMblks(carrierFluidId)%xmin,LBMblks(carrierFluidId)%ymin,LBMblks(carrierFluidId)%zmin, &
                                LBMblks(carrierFluidId)%xDim,LBMblks(carrierFluidId)%yDim,LBMblks(carrierFluidId)%zDim,LBMblks(carrierFluidId)%uuu,LBMblks(carrierFluidId)%force)
        call get_now_time(time_end2)
        time_IBM = time_IBM + (time_end2 - time_begine2)
        call get_now_time(time_begine2)
        do isubstep=1,flow%numsubstep
            call Solver(time,isubstep,LBMblks(carrierFluidId)%dh,dt_solid)
        enddo !do isubstep=1,numsubstep
        call get_now_time(time_end2)
        time_FEM = time_FEM + (time_end2 - time_begine2)
    endsubroutine

    subroutine extract_inner_layer(treenode,time)
        integer:: treenode, time, f, s, i, ns
        integer:: xF, yF, zF, x, y, z
        f = treenode
        ns = blockTree(f)%nsons
        if(ns .ne. 0) then
            do i = 1,ns
                s = blockTree(f)%sons(i)
                if(LBMblks(s)%BndConds(1).eq.BCfluid) then
                    xF = blockTree(f)%comm(i)%f(1)
                    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(y,z)
                    do y = 1,LBMblks(f)%yDim
                    do z = 1,LBMblks(f)%zDim
                        if (time .eq. 1) LBMblks(s)%fIn_Fx1t1(z,y,:) = LBMblks(f)%fIn(z,y,xF,:)
                        if (time .eq. 2) LBMblks(s)%fIn_Fx1t2(z,y,:) = LBMblks(f)%fIn(z,y,xF,:)
                        if (time .eq. 1) LBMblks(s)%uuu_Fx1t1(z,y,:) = LBMblks(f)%uuu(z,y,xF,:)
                    enddo
                    enddo
                    !$OMP END PARALLEL DO
                endif
                if(LBMblks(s)%BndConds(2).eq.BCfluid) then
                    xF = blockTree(f)%comm(i)%f(2)
                    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(y,z)
                    do y = 1,LBMblks(f)%yDim
                    do z = 1,LBMblks(f)%zDim
                        if (time .eq. 1) LBMblks(s)%fIn_Fx2t1(z,y,:) = LBMblks(f)%fIn(z,y,xF,:)
                        if (time .eq. 2) LBMblks(s)%fIn_Fx2t2(z,y,:) = LBMblks(f)%fIn(z,y,xF,:)
                        if (time .eq. 1) LBMblks(s)%uuu_Fx2t1(z,y,:) = LBMblks(f)%uuu(z,y,xF,:)
                    enddo
                    enddo
                    !$OMP END PARALLEL DO
                endif
                if(LBMblks(s)%BndConds(3).eq.BCfluid) then
                    yF = blockTree(f)%comm(i)%f(3)
                    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,z)
                    do x = 1,LBMblks(f)%xDim
                    do z = 1,LBMblks(f)%zDim
                        if (time .eq. 1) LBMblks(s)%fIn_Fy1t1(z,x,:) = LBMblks(f)%fIn(z,yF,x,:)
                        if (time .eq. 2) LBMblks(s)%fIn_Fy1t2(z,x,:) = LBMblks(f)%fIn(z,yF,x,:)
                        if (time .eq. 1) LBMblks(s)%uuu_Fy1t1(z,x,:) = LBMblks(f)%uuu(z,yF,x,:)
                    enddo
                    enddo
                    !$OMP END PARALLEL DO
                endif
                if(LBMblks(s)%BndConds(4).eq.BCfluid) then
                    yF = blockTree(f)%comm(i)%f(4)
                    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,z)
                    do x = 1,LBMblks(f)%xDim
                    do z = 1,LBMblks(f)%zDim
                        if (time .eq. 1) LBMblks(s)%fIn_Fy2t1(z,x,:) = LBMblks(f)%fIn(z,yF,x,:)
                        if (time .eq. 2) LBMblks(s)%fIn_Fy2t2(z,x,:) = LBMblks(f)%fIn(z,yF,x,:)
                        if (time .eq. 1) LBMblks(s)%uuu_Fy2t1(z,x,:) = LBMblks(f)%uuu(z,yF,x,:)
                    enddo
                    enddo
                    !$OMP END PARALLEL DO
                endif
                if(LBMblks(s)%BndConds(5).eq.BCfluid) then
                    zF = blockTree(f)%comm(i)%f(5)
                    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y)
                    do x = 1,LBMblks(f)%xDim
                    do y = 1,LBMblks(f)%yDim
                        if (time .eq. 1) LBMblks(s)%fIn_Fz1t1(y,x,:) = LBMblks(f)%fIn(zF,y,x,:)
                        if (time .eq. 2) LBMblks(s)%fIn_Fz1t2(y,x,:) = LBMblks(f)%fIn(zF,y,x,:)
                        if (time .eq. 1) LBMblks(s)%uuu_Fz1t1(y,x,:) = LBMblks(f)%uuu(zF,y,x,:)
                    enddo
                    enddo
                    !$OMP END PARALLEL DO
                endif
                if(LBMblks(s)%BndConds(6).eq.BCfluid) then
                    zF = blockTree(f)%comm(i)%f(6)
                    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y)
                    do x = 1,LBMblks(f)%xDim
                    do y = 1,LBMblks(f)%yDim
                        if (time .eq. 1) LBMblks(s)%fIn_Fz2t1(y,x,:) = LBMblks(f)%fIn(zF,y,x,:)
                        if (time .eq. 2) LBMblks(s)%fIn_Fz2t2(y,x,:) = LBMblks(f)%fIn(zF,y,x,:)
                        if (time .eq. 1) LBMblks(s)%uuu_Fz2t1(y,x,:) = LBMblks(f)%uuu(zF,y,x,:)
                    enddo
                    enddo
                    !$OMP END PARALLEL DO
                endif
            enddo
        endif
    endsubroutine

    ! verify the parameters of fluid blocks
    recursive SUBROUTINE check_blocks_params(treenode)
        implicit none
        integer,intent(in):: treenode
        integer:: i,nblock
        type(CommPair)::p
        logical:: flag
        real(8)::residual
        flag = .false. ! default flase
        do i=1,blocktree(treenode)%nsons
            p = blocktree(treenode)%comm(i)
            flag =  abs(LBMblks(p%fatherId)%dh - LBMblks(p%sonId)%dh*m_gridDelta).gt.1d-8 .or. &
                    mod(LBMblks(p%sonId)%xDim,m_gridDelta) .ne. 1 .or. &
                    mod(LBMblks(p%sonId)%yDim,m_gridDelta) .ne. 1 .or. &
                    mod(LBMblks(p%sonId)%zDim,m_gridDelta) .ne. 1
            residual = (LBMblks(p%sonId)%xmin-LBMblks(p%fatherId)%xmin)/LBMblks(p%fatherId)%dh + &
                       (LBMblks(p%sonId)%ymin-LBMblks(p%fatherId)%ymin)/LBMblks(p%fatherId)%dh + &
                       (LBMblks(p%sonId)%zmin-LBMblks(p%fatherId)%zmin)/LBMblks(p%fatherId)%dh
            residual = abs(residual - dble(NINT(residual)))
            if(flag .or. residual .gt. 1d-8) then
                write(*,*) 'grid points do not match between fluid blocks',p%fatherId,p%sonId
                stop
            endif
            call check_blocks_params(blocktree(treenode)%sons(i))
        enddo
    END SUBROUTINE

    SUBROUTINE deliver_son_to_father(pair)
        use FluidDomain
        implicit none
        type(CommPair),intent(in):: pair
        integer:: xS,yS,zS,xF,yF,zF
        ! x direction slices
        if(pair%sds(1).eq.1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(yS,zS,xF,yF,zF)
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                xS=                         1+m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,xF,yF,zF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(2).eq.-1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(yS,zS,xF,yF,zF)
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                xS=  LBMblks(pair%sonId)%xDim-m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,xF,yF,zF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        ! y direction slices
        if(pair%sds(3).eq.1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,zS,xF,yF,zF)
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                yS=                         1+m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,xF,yF,zF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(4).eq.-1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,zS,xF,yF,zF)
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  zS=1,LBMblks(pair%sonId)%zDim,m_gridDelta
                yS=  LBMblks(pair%sonId)%yDim-m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,xF,yF,zF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        ! z direction slices
        if(pair%sds(5).eq.1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,yS,xF,yF,zF)
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
                zS=                         1+m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,xF,yF,zF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
        !$OMP END PARALLEL DO
        endif
        if(pair%sds(6).eq.-1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,yS,xF,yF,zF)
            do  xS=1,LBMblks(pair%sonId)%xDim,m_gridDelta
            do  yS=1,LBMblks(pair%sonId)%yDim,m_gridDelta
                zS=  LBMblks(pair%sonId)%zDim-m_gridDelta
                call deliver_grid_distribution(pair%fatherId,pair%sonId,xS,yS,zS,xF,yF,zF)
                call fIn_son_to_father(pair%fatherId,pair%sonId,xF,yF,zF)
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif
    end subroutine 

    SUBROUTINE deliver_grid_distribution(father,son,xS,yS,zS,xF,yF,zF)
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
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(yS,zS)
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
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(2).eq.-1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(yS,zS)
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
            !$OMP END PARALLEL DO
        endif
        ! y direction
        if(pair%sds(3).eq.1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,zS)
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
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(4).eq.-1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,zS)
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
            !$OMP END PARALLEL DO
        endif
        ! z direction
        if(pair%sds(5).eq.1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,yS)
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
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(6).eq.-1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,yS)
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
            !$OMP END PARALLEL DO
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
