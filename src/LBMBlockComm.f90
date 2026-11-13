! SPDX-License-Identifier: GPL-3.0-or-later
!
! FSILBM3D
! Copyright (C) 2025-2026 Ankang Gao and contributors

module LBMBlockComm
    use ConstParams
    use FluidDomain
    use ParallelDomain, only: sum_face, parallel_abort
    implicit none
    !include 'mpif.h'
    type :: CommPair
        integer:: fatherId
        integer:: sonId
        integer:: sds(1:6) ! 0 no need communication; 1(min),-1(max) need communication
        integer:: s(1:6),f(1:6),si(1:6),fi(1:6) ! s son's boundary layer; si son's first inner layer
        integer:: xDimS, xDimF, yDimS, yDimF, zDimS, zDimF
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
        integer:: nsons,i,j,sxD,syD,szD,ratio,fId,sId
        type(CommPair)::pair
        nsons = blockTree(treenode)%nsons
        if(nsons .eq. 0) return
        if(.not.allocated(blockTree(treenode)%comm)) then
            allocate(blockTree(treenode)%comm(nsons))
        endif
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
                if (LBMblks(sId)%periodic_bc(1) .eq. 1) then
                    sxD = sxD - 1
                endif
                if (LBMblks(sId)%periodic_bc(2) .eq. 1) then
                    syD = syD - 1
                endif
                if (LBMblks(sId)%periodic_bc(3) .eq. 1) then
                    szD = szD - 1
                endif
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
                pair%xDimS = LBMblks(sId)%xDim
                pair%xDimF = pair%f(2) - pair%f(1) + 1
                pair%yDimS = LBMblks(sId)%yDim
                pair%yDimF = pair%f(4) - pair%f(3) + 1
                pair%zDimS = LBMblks(sId)%zDim
                pair%zDimF = pair%f(6) - pair%f(5) + 1
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
            write(*,*) 'Error: there exist more than one block tree root, father relationship(-1 is root):', fa
            stop
        endif
    end subroutine

    recursive subroutine array_to_tree(iblocks, nb, rootnode)
        implicit none
        integer, intent(in):: nb, iblocks(1:nb),rootnode
        integer:: roots(1:nb),fa(1:nb),subblock(1:nb)
        integer:: i, j, tmp,cr,cb
        if(nb.eq.0) return
        ! compare all blocks
        fa = -1
        do i=1,nb-1
            do j=i+1,nb
                tmp = CompareBlocks(iblocks(i), iblocks(j))
                if(tmp.eq.1) then
                    fa(j) = i
                elseif(tmp.eq.-1) then
                    fa(i) = j
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
                if(fa(j).gt.0) then
                    if(iblocks(fa(j)).eq.roots(i)) then
                        cb = cb + 1
                        subblock(cb) = iblocks(j)
                    endif
                endif
            enddo
            call array_to_tree(subblock, cb, roots(i))
        enddo
    end subroutine

    subroutine build_block_tree()
        implicit none
        integer:: i,iblocks(1:m_nblocks),nb
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
        call allocate_fIn_tau(blockTreeRoot)
    endsubroutine build_block_tree

    recursive subroutine allocate_fIn_tau(treenode)
        implicit none
        type(CommPair):: pair
        integer,intent(in):: treenode
        integer:: i,f,s,ns,xDimF,yDimF,zDimF
        f  = treenode
        ns = blockTree(f)%nsons
        do i = 1,ns
            s = blockTree(f)%sons(i)
            pair  = blockTree(f)%comm(i)
            xDimF = pair%xDimF
            yDimF = pair%yDimF
            zDimF = pair%zDimF
            if(LBMblks(s)%BndConds(1).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fx1t1(0:lbmDim,zDimF,yDimF))
                allocate(LBMblks(s)%fIn_Fx1t2(0:lbmDim,zDimF,yDimF))
                allocate(LBMblks(s)%tau_Fx1t1(zDimF,yDimF))
                allocate(LBMblks(s)%tau_Fx1t2(zDimF,yDimF))
            endif
            if(LBMblks(s)%BndConds(2).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fx2t1(0:lbmDim,zDimF,yDimF))
                allocate(LBMblks(s)%fIn_Fx2t2(0:lbmDim,zDimF,yDimF))
                allocate(LBMblks(s)%tau_Fx2t1(zDimF,yDimF))
                allocate(LBMblks(s)%tau_Fx2t2(zDimF,yDimF))
            endif
            if(LBMblks(s)%BndConds(3).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fy1t1(0:lbmDim,zDimF,xDimF))
                allocate(LBMblks(s)%fIn_Fy1t2(0:lbmDim,zDimF,xDimF))
                allocate(LBMblks(s)%tau_Fy1t1(zDimF,xDimF))
                allocate(LBMblks(s)%tau_Fy1t2(zDimF,xDimF))
            endif
            if(LBMblks(s)%BndConds(4).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fy2t1(0:lbmDim,zDimF,xDimF))
                allocate(LBMblks(s)%fIn_Fy2t2(0:lbmDim,zDimF,xDimF))
                allocate(LBMblks(s)%tau_Fy2t1(zDimF,xDimF))
                allocate(LBMblks(s)%tau_Fy2t2(zDimF,xDimF))
            endif
            if(LBMblks(s)%BndConds(5).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fz1t1(0:lbmDim,yDimF,xDimF))
                allocate(LBMblks(s)%fIn_Fz1t2(0:lbmDim,yDimF,xDimF))
                allocate(LBMblks(s)%tau_Fz1t1(yDimF,xDimF))
                allocate(LBMblks(s)%tau_Fz1t2(yDimF,xDimF))
            endif
            if(LBMblks(s)%BndConds(6).eq.BCfluid) then
                allocate(LBMblks(s)%fIn_Fz2t1(0:lbmDim,yDimF,xDimF))
                allocate(LBMblks(s)%fIn_Fz2t2(0:lbmDim,yDimF,xDimF))
                allocate(LBMblks(s)%tau_Fz2t1(yDimF,xDimF))
                allocate(LBMblks(s)%tau_Fz2t2(yDimF,xDimF))
            endif
            call allocate_fIn_tau(blockTree(f)%sons(i))
        enddo
    endsubroutine

    recursive subroutine tree_set_boundary_conditions_block(treenode)
        implicit none
        integer:: i, s, treenode
        call set_boundary_conditions_block(treenode)
        ! tree cycle
        if(blockTree(treenode)%nsons.gt.0) then
            do i=1,blockTree(treenode)%nsons
                s = blockTree(treenode)%sons(i)
                call tree_set_boundary_conditions_block(s)
            enddo
        endif
    endsubroutine tree_set_boundary_conditions_block

    recursive subroutine tree_collision_streaming_IBM_FEM(treenode,time_collision,time_streaming,time_IBM,time_FEM)
        implicit none
        integer:: i, s, treenode, n_timeStep
        real(8):: time_collision,time_streaming,time_IBM,time_FEM,time_begine2,time_end2
        call LBMblks(treenode)%update_volume_force()
        ! calculate macro quantities for each blocks,must be ahead of collision(Huang Haibo 2024 P162)
        call LBMblks(treenode)%calculate_macro_quantities()
        call LBMblks(treenode)%ResetVolumeForce()
        call IBM_FEM(treenode,time_IBM,time_FEM,LBMblks(treenode)%blktime)
        call LBMblks(treenode)%add_volume_force()
        ! extract interpolation layer for old time
        call extract_interpolate_layer(treenode,1)
        ! collision
        call get_now_time(time_begine2)
        call collision_block(treenode)
        call get_now_time(time_end2)
        time_collision = time_collision + (time_end2 - time_begine2)
        call halfwayBCset_block(treenode)
        ! streaming
        call get_now_time(time_begine2)
        call streaming_block(treenode)
        call get_now_time(time_end2)
        time_streaming = time_streaming + (time_end2 - time_begine2)
        !set boundary
        call set_boundary_conditions_block(treenode)
        ! extract interpolation layer for new time
        call extract_interpolate_layer(treenode,2)
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
        integer:: carrierFluidId,isubstep
        dt_solid = LBMblks(carrierFluidId)%dh/dble(flow%numsubstep)       !time step of the solid
        call get_now_time(time_begine2)
        associate(b => LBMblks(carrierFluidId))
            call FSInteraction_force(b%carriedBodies,b%dh,b%dh,b%global_xmin,b%ymin,b%zmin, &
                                     b%xDim,b%yDim,b%zDim,b%global_nx,b%first_x,b%den,b%uuu,b%force)
        end associate
        call get_now_time(time_end2)
        time_IBM = time_IBM + (time_end2 - time_begine2)
        call get_now_time(time_begine2)
        do isubstep=1,flow%numsubstep
            call Solver(LBMblks(carrierFluidId)%carriedBodies,time,isubstep,LBMblks(carrierFluidId)%dh,dt_solid)
        enddo !do isubstep=1,numsubstep
        call get_now_time(time_end2)
        time_FEM = time_FEM + (time_end2 - time_begine2)
    endsubroutine

    subroutine extract_interpolate_layer(treenode,time)
        integer, intent(in) :: treenode, time
        integer :: i, face
        type(CommPair) :: pair
        real(8), allocatable :: values(:,:,:)
        do i = 1, blockTree(treenode)%nsons
            pair = blockTree(treenode)%comm(i)
            associate(son => LBMblks(pair%sonId))
                do face = 1, 6
                    if (pair%sds(face) == 0) cycle
                    call capture_face(LBMblks(treenode), pair%f, face, 1, values)
                    select case (face)
                    case (1)
                        call store_face(son%fIn_Fx1t1, son%fIn_Fx1t2, son%tau_Fx1t1, son%tau_Fx1t2)
                    case (2)
                        call store_face(son%fIn_Fx2t1, son%fIn_Fx2t2, son%tau_Fx2t1, son%tau_Fx2t2)
                    case (3)
                        call store_face(son%fIn_Fy1t1, son%fIn_Fy1t2, son%tau_Fy1t1, son%tau_Fy1t2)
                    case (4)
                        call store_face(son%fIn_Fy2t1, son%fIn_Fy2t2, son%tau_Fy2t1, son%tau_Fy2t2)
                    case (5)
                        call store_face(son%fIn_Fz1t1, son%fIn_Fz1t2, son%tau_Fz1t1, son%tau_Fz1t2)
                    case (6)
                        call store_face(son%fIn_Fz2t1, son%fIn_Fz2t2, son%tau_Fz2t1, son%tau_Fz2t2)
                    end select
                enddo
            end associate
        enddo
    contains
        subroutine store_face(f1, f2, tau1, tau2)
            real(8), intent(inout) :: f1(:,:,:), f2(:,:,:), tau1(:,:), tau2(:,:)
            if (time == 1) then
                f1 = values(0:lbmDim,:,:)
                tau1 = values(lbmDim+1,:,:)
            else
                f2 = values(0:lbmDim,:,:)
                tau2 = values(lbmDim+1,:,:)
                f1 = 0.5d0 * (f1 + f2)
                tau1 = 0.5d0 * (tau1 + tau2)
            endif
        end subroutine
    end subroutine

    subroutine face_geometry(bounds, face, step, axes, n)
        integer, intent(in) :: bounds(6), face, step
        integer, intent(out) :: axes(3), n(2)
        integer, parameter :: axis_order(3,3) = reshape([3,2,1, 3,1,2, 2,1,3], [3,3])
        axes = axis_order(:,(face+1)/2)
        n = max(0, (bounds(2*axes(1:2))-bounds(2*axes(1:2)-1)+step)/step)
    end subroutine

    subroutine capture_face(block, bounds, face, step, values)
        type(LBMBlock), intent(in) :: block
        integer, intent(in) :: bounds(6), face, step
        real(8), allocatable, intent(out) :: values(:,:,:)
        integer :: a, b, xyz(3), axes(3), n(2), x, y, z
        call face_geometry(bounds, face, step, axes, n)
        allocate(values(0:lbmDim+1,n(1),n(2)))
        values = 0.0d0
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b,xyz,x,y,z)
        do b = 1, n(2)
            do a = 1, n(1)
                xyz(axes(1:2)) = bounds(2*axes(1:2)-1) + ([a,b]-1)*step
                xyz(axes(3)) = bounds(face)
                if (.not. block%owns_x(xyz(1))) cycle
                x = xyz(1) - block%first_x + 1
                y = xyz(2)
                z = xyz(3)
                values(0:lbmDim,a,b) = block%fIn(z,y,x,0:lbmDim)
                values(lbmDim+1,a,b) = block%tau_all(z,y,x)
            enddo
        enddo
        !$OMP END PARALLEL DO
        call sum_face(values)
    end subroutine

    recursive SUBROUTINE check_blocks_params(treenode)
        implicit none
        integer,intent(in):: treenode
        integer:: i,r(3),bc_pair
        type(CommPair)::p
        logical:: flag
        real(8)::res1,res2,res
        flag = .false. ! default flase
        do i=1,blocktree(treenode)%nsons
            p = blocktree(treenode)%comm(i)
            r = 1
            do bc_pair = 1,3
                if (LBMblks(p%sonId)%periodic_bc(bc_pair) .eq. 1) then
                    r(bc_pair) = 0
                endif
            enddo
            flag =  abs(LBMblks(p%fatherId)%dh - LBMblks(p%sonId)%dh*dble(m_gridDelta)).gt.1d-8 .or. &
                    mod(LBMblks(p%sonId)%xDim,m_gridDelta) .ne. r(1) .or. &
                    mod(LBMblks(p%sonId)%yDim,m_gridDelta) .ne. r(2) .or. &
                    mod(LBMblks(p%sonId)%zDim,m_gridDelta) .ne. r(3)
            res1 =  (LBMblks(p%sonId)%xmin-LBMblks(p%fatherId)%xmin)/LBMblks(p%fatherId)%dh + &
                    (LBMblks(p%sonId)%ymin-LBMblks(p%fatherId)%ymin)/LBMblks(p%fatherId)%dh + &
                    (LBMblks(p%sonId)%zmin-LBMblks(p%fatherId)%zmin)/LBMblks(p%fatherId)%dh
            res2 =  (LBMblks(p%sonId)%xmax-LBMblks(p%fatherId)%xmin)/LBMblks(p%fatherId)%dh + &
                    (LBMblks(p%sonId)%ymax-LBMblks(p%fatherId)%ymin)/LBMblks(p%fatherId)%dh + &
                    (LBMblks(p%sonId)%zmax-LBMblks(p%fatherId)%zmin)/LBMblks(p%fatherId)%dh
            res1 = abs(res1 - dble(NINT(res1)))
            res2 = abs(res2 - dble(NINT(res2)))
            res = res1 + res2
            if(flag .or. res .gt. 1d-8) then
                write(*,*) 'grid points do not match between fluid blocks',LBMBlks(p%fatherId)%ID,LBMBlks(p%sonId)%ID,res1,res2
                write(*,*) 'if son block have periodic boundarys, an even number of grid points is needed. Otherwise an odd number is needed.'
                stop
            endif
            call check_blocks_params(blocktree(treenode)%sons(i))
        enddo
    END SUBROUTINE

    subroutine deliver_son_to_father(pair)
        type(CommPair), intent(in) :: pair
        integer :: face, a, b, xyz(3), axes(3), n(2), x, y, z
        real(8) :: populations(0:lbmDim), coeff, force(3), dh
        real(8), allocatable :: values(:,:,:)
        associate(son => LBMblks(pair%sonId), father => LBMblks(pair%fatherId))
            force = son%volumeForce
            dh = son%dh
            do face = 1, 6
                if (pair%sds(face) == 0) cycle
                call capture_face(son, pair%si, face, m_gridDelta, values)
                call face_geometry(pair%fi, face, 1, axes, n)
                !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b,xyz,x,y,z,populations,coeff)
                do b = 1, size(values,3)
                    do a = 1, size(values,2)
                        xyz(axes(1:2)) = pair%fi(2*axes(1:2)-1) + [a,b] - 1
                        xyz(axes(3)) = pair%fi(face)
                        if (.not. father%owns_x(xyz(1))) cycle
                        x = xyz(1) - father%first_x + 1
                        y = xyz(2)
                        z = xyz(3)
                        coeff = (father%tau_all(z,y,x) / values(lbmDim+1,a,b)) * dble(m_gridDelta)
                        populations = values(0:lbmDim,a,b)
                        call fIn_GridTransform(populations, coeff, force, dh)
                        father%fIn(z,y,x,0:lbmDim) = populations
                    enddo
                enddo
                !$OMP END PARALLEL DO
            enddo
        end associate
    end subroutine

    SUBROUTINE interpolation_father_to_son(pair,n_timeStep)
        implicit none
        type(CommPair),intent(in):: pair
        integer:: n_timeStep,xS,yS,zS,xLocal
        real(8):: coeff,VF(1:3),dh
        integer:: xDimS, xDimF, yDimS, yDimF, zDimS, zDimF
        real(8),allocatable::tmpf(:,:,:)
        real(8),allocatable::tmptau(:,:)
        VF    = LBMblks(pair%fatherId)%volumeForce(1:3)
        dh    = LBMblks(pair%fatherId)%dh
        xDimS = pair%xDimS
        xDimF = pair%xDimF
        yDimS = pair%yDimS
        yDimF = pair%yDimF
        zDimS = pair%zDimS
        zDimF = pair%zDimF
        ! x direction
        if(pair%sds(1).eq.1 .or. pair%sds(2).eq.-1) then
            allocate(tmpf(0:lbmDim,zDimS,yDimS))
            allocate(tmptau(zDimS,yDimS))
        endif
        if(pair%sds(1).eq.1 .and. LBMblks(pair%sonId)%owns_x(pair%s(1))) then
            xS = pair%s(1) - LBMblks(pair%sonId)%first_x + 1
            if(n_timeStep.eq.0) then
                call interpolate_fIn(zDimF,yDimF,LBMblks(pair%sonId)%fIn_Fx1t1,zDimS,yDimS,tmpf)
                call interpolate_tau(zDimF,yDimF,LBMblks(pair%sonId)%tau_Fx1t1,zDimS,yDimS,tmptau)
            else
                call interpolate_fIn(zDimF,yDimF,LBMblks(pair%sonId)%fIn_Fx1t2,zDimS,yDimS,tmpf)
                call interpolate_tau(zDimF,yDimF,LBMblks(pair%sonId)%tau_Fx1t2,zDimS,yDimS,tmptau)
            endif
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(yS,zS,coeff)
            do  yS = 1,yDimS
                do  zS = 1,zDimS
                    coeff = (LBMblks(pair%sonId)%tau_all(zS,yS,xS) / tmptau(zS,yS)) / dble(m_gridDelta)
                    call fIn_GridTransform(tmpf(:,zS,yS), coeff, VF, dh)
                    LBMblks(pair%sonId)%fIn(zS,yS,xS,0:lbmDim) = tmpf(:,zS,yS)
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(2).eq.-1 .and. LBMblks(pair%sonId)%owns_x(pair%s(2))) then
            xS = pair%s(2) - LBMblks(pair%sonId)%first_x + 1
            if(n_timeStep.eq.0) then
                call interpolate_fIn(zDimF,yDimF,LBMblks(pair%sonId)%fIn_Fx2t1,zDimS,yDimS,tmpf)
                call interpolate_tau(zDimF,yDimF,LBMblks(pair%sonId)%tau_Fx2t1,zDimS,yDimS,tmptau)
            else
                call interpolate_fIn(zDimF,yDimF,LBMblks(pair%sonId)%fIn_Fx2t2,zDimS,yDimS,tmpf)
                call interpolate_tau(zDimF,yDimF,LBMblks(pair%sonId)%tau_Fx2t2,zDimS,yDimS,tmptau)
            endif
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(yS,zS,coeff)
            do  yS = 1,yDimS
                do  zS = 1,zDimS
                    coeff = (LBMblks(pair%sonId)%tau_all(zS,yS,xS) / tmptau(zS,yS)) / dble(m_gridDelta)
                    call fIn_GridTransform(tmpf(:,zS,yS), coeff, VF, dh)
                    LBMblks(pair%sonId)%fIn(zS,yS,xS,0:lbmDim) = tmpf(:,zS,yS)
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(allocated(tmpf)) deallocate(tmpf)
        if(allocated(tmptau)) deallocate(tmptau)
        ! y direction
        if(pair%sds(3).eq.1 .or. pair%sds(4).eq.-1) then
            allocate(tmpf(0:lbmDim,zDimS,xDimS))
            allocate(tmptau(zDimS,xDimS))
        endif
        if(pair%sds(3).eq.1) then
            yS = pair%s(3)
            if(n_timeStep.eq.0) then
                call interpolate_fIn(zDimF,xDimF,LBMblks(pair%sonId)%fIn_Fy1t1,zDimS,xDimS,tmpf)
                call interpolate_tau(zDimF,xDimF,LBMblks(pair%sonId)%tau_Fy1t1,zDimS,xDimS,tmptau)
            else
                call interpolate_fIn(zDimF,xDimF,LBMblks(pair%sonId)%fIn_Fy1t2,zDimS,xDimS,tmpf)
                call interpolate_tau(zDimF,xDimF,LBMblks(pair%sonId)%tau_Fy1t2,zDimS,xDimS,tmptau)
            endif
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,zS,xLocal,coeff)
            do xS = LBMblks(pair%sonId)%first_x, LBMblks(pair%sonId)%first_x + LBMblks(pair%sonId)%xDim - 1
                xLocal = xS - LBMblks(pair%sonId)%first_x + 1
                do  zS = 1,zDimS
                    coeff = (LBMblks(pair%sonId)%tau_all(zS,yS,xLocal) / tmptau(zS,xS)) / dble(m_gridDelta)
                    call fIn_GridTransform(tmpf(:,zS,xS), coeff, VF, dh)
                    LBMblks(pair%sonId)%fIn(zS,yS,xLocal,0:lbmDim) = tmpf(:,zS,xS)
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(4).eq.-1) then
            yS = pair%s(4)
            if(n_timeStep.eq.0) then
                call interpolate_fIn(zDimF,xDimF,LBMblks(pair%sonId)%fIn_Fy2t1,zDimS,xDimS,tmpf)
                call interpolate_tau(zDimF,xDimF,LBMblks(pair%sonId)%tau_Fy2t1,zDimS,xDimS,tmptau)
            else
                call interpolate_fIn(zDimF,xDimF,LBMblks(pair%sonId)%fIn_Fy2t2,zDimS,xDimS,tmpf)
                call interpolate_tau(zDimF,xDimF,LBMblks(pair%sonId)%tau_Fy2t2,zDimS,xDimS,tmptau)
            endif
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,zS,xLocal,coeff)
            do xS = LBMblks(pair%sonId)%first_x, LBMblks(pair%sonId)%first_x + LBMblks(pair%sonId)%xDim - 1
                xLocal = xS - LBMblks(pair%sonId)%first_x + 1
                do  zS = 1,zDimS
                    coeff = (LBMblks(pair%sonId)%tau_all(zS,yS,xLocal) / tmptau(zS,xS)) / dble(m_gridDelta)
                    call fIn_GridTransform(tmpf(:,zS,xS), coeff, VF, dh)
                    LBMblks(pair%sonId)%fIn(zS,yS,xLocal,0:lbmDim) = tmpf(:,zS,xS)
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(allocated(tmpf)) deallocate(tmpf)
        if(allocated(tmptau)) deallocate(tmptau)
        ! z direction
        if(pair%sds(5).eq.1 .or. pair%sds(6).eq.-1) then
            allocate(tmpf(0:lbmDim,yDimS,xDimS))
            allocate(tmptau(yDimS,xDimS))
        endif
        if(pair%sds(5).eq.1) then
            zS = pair%s(5)
            if(n_timeStep.eq.0) then
                call interpolate_fIn(yDimF,xDimF,LBMblks(pair%sonId)%fIn_Fz1t1,yDimS,xDimS,tmpf)
                call interpolate_tau(yDimF,xDimF,LBMblks(pair%sonId)%tau_Fz1t1,yDimS,xDimS,tmptau)
            else
                call interpolate_fIn(yDimF,xDimF,LBMblks(pair%sonId)%fIn_Fz1t2,yDimS,xDimS,tmpf)
                call interpolate_tau(yDimF,xDimF,LBMblks(pair%sonId)%tau_Fz1t2,yDimS,xDimS,tmptau)
            endif
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,yS,xLocal,coeff)
            do xS = LBMblks(pair%sonId)%first_x, LBMblks(pair%sonId)%first_x + LBMblks(pair%sonId)%xDim - 1
                xLocal = xS - LBMblks(pair%sonId)%first_x + 1
                do  yS = 1,yDimS
                    coeff = (LBMblks(pair%sonId)%tau_all(zS,yS,xLocal) / tmptau(yS,xS)) / dble(m_gridDelta)
                    call fIn_GridTransform(tmpf(:,yS,xS), coeff, VF, dh)
                    LBMblks(pair%sonId)%fIn(zS,yS,xLocal,0:lbmDim) = tmpf(:,yS,xS)
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(pair%sds(6).eq.-1) then
            zS = pair%s(6)
            if(n_timeStep.eq.0) then
                call interpolate_fIn(yDimF,xDimF,LBMblks(pair%sonId)%fIn_Fz2t1,yDimS,xDimS,tmpf)
                call interpolate_tau(yDimF,xDimF,LBMblks(pair%sonId)%tau_Fz2t1,yDimS,xDimS,tmptau)
            else
                call interpolate_fIn(yDimF,xDimF,LBMblks(pair%sonId)%fIn_Fz2t2,yDimS,xDimS,tmpf)
                call interpolate_tau(yDimF,xDimF,LBMblks(pair%sonId)%tau_Fz2t2,yDimS,xDimS,tmptau)
            endif
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(xS,yS,xLocal,coeff)
            do xS = LBMblks(pair%sonId)%first_x, LBMblks(pair%sonId)%first_x + LBMblks(pair%sonId)%xDim - 1
                xLocal = xS - LBMblks(pair%sonId)%first_x + 1
                do  yS = 1,yDimS
                    coeff = (LBMblks(pair%sonId)%tau_all(zS,yS,xLocal) / tmptau(yS,xS)) / dble(m_gridDelta)
                    call fIn_GridTransform(tmpf(:,yS,xS), coeff, VF, dh)
                    LBMblks(pair%sonId)%fIn(zS,yS,xLocal,0:lbmDim) = tmpf(:,yS,xS)
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(allocated(tmpf)) deallocate(tmpf)
        if(allocated(tmptau)) deallocate(tmptau)
    end subroutine

    subroutine interpolate_fIn(bF,aF,fF,bS,aS,fS) ! space interpolation : 2 for 3rd- and 4th-order, other for linear(lagrange interpolate)
        use FlowCondition, only: flow
        implicit none
        integer,intent(in):: aS,bS,aF,bF
        real(8),intent(in):: fF(0:lbmDim,bF,aF)
        real(8),intent(out):: fS(0:lbmDim,bS,aS) !fine grid values
        integer:: a,b,a1,b1,r(2),aStmp,bStmp
        r = 0
        bStmp = bS
        aStmp = aS
        if (mod(bS,2).eq.0) then
            bStmp = bS - 1
            r(2) = 1
        endif
        if (mod(aS,2).eq.0) then
            aStmp = aS - 1
            r(1) = 1
        endif
        if (flow%interpolateScheme.eq.2) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b,a1,b1)
            do b  = 1,bStmp,2
               b1 = b / 2 + 1
            do a  = 1,aStmp,2
               a1 = a / 2 + 1
                fS(:,b,a) = fF(:,b1,a1)
                if(1.eq.b) then
                    fS(:,b+1,a) = 0.375d0*fF(:,b1,a1) + 0.75d0*fF(:,b1+1,a1) - 0.125d0*fF(:,b1+2,a1)
                else if(b.eq.bStmp-2) then
                    fS(:,b+1,a) = 0.375d0*fF(:,b1+1,a1) + 0.75d0*fF(:,b1,a1) - 0.125d0*fF(:,b1-1,a1)
                else if(b.ne.bStmp) then
                    fS(:,b+1,a) = -0.0625d0*fF(:,b1-1,a1) + 0.5625d0*fF(:,b1,a1) + 0.5625d0*fF(:,b1+1,a1) - 0.0625d0*fF(:,b1+2,a1)
                endif
            enddo
            enddo
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b)
            do b = 1,bStmp
            do a = 2,aStmp,2
                if(2.eq.a) then
                    fS(:,b,a) = 0.375d0*fS(:,b,a-1) + 0.75d0*fS(:,b,a+1) - 0.125d0*fS(:,b,a+3)
                else if(a.eq.aStmp-1) then
                    fS(:,b,a) = 0.375d0*fS(:,b,a+1) + 0.75d0*fS(:,b,a-1) - 0.125d0*fS(:,b,a-3)
                else
                    fS(:,b,a) = -0.0625d0*fS(:,b,a-3) + 0.5625d0*fS(:,b,a-1) + 0.5625d0*fS(:,b,a+1) - 0.0625d0*fS(:,b,a+3)
                endif
            enddo
            enddo
            !$OMP END PARALLEL DO
            if (r(2).eq.1) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a)
                do a = 1,aStmp
                    fS(:,bStmp+1,a) = -0.0625d0*fS(:,bStmp-2,a) + 0.5625d0*fS(:,bStmp,a) + 0.5625d0*fS(:,1,a) - 0.0625d0*fS(:,3,a)
                enddo
                !$OMP END PARALLEL DO
            endif
            if (r(1).eq.1) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(b)
                do b = 1,bStmp
                    fS(:,b,aStmp+1) = -0.0625d0*fS(:,b,aStmp-2) + 0.5625d0*fS(:,b,aStmp) + 0.5625d0*fS(:,b,1) - 0.0625d0*fS(:,b,3)
                enddo
                !$OMP END PARALLEL DO
                if (r(2).eq.1) fS(:,bStmp+1,aStmp+1) = -0.0625d0*fS(:,bStmp+1,aStmp-2) + 0.5625d0*fS(:,bStmp+1,aStmp) + 0.5625d0*fS(:,bStmp+1,1) - 0.0625d0*fS(:,bStmp+1,3)
            endif
        else
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b,a1,b1)
            do b  = 1,bStmp,2
               b1 = b / 2 + 1
            do a  = 1,aStmp,2
               a1 = a / 2 + 1
                fS(:,b,a) = fF(:,b1,a1)
                if(b.lt.bStmp) fS(:,b+1,a) = (fF(:,b1,a1) + fF(:,b1+1,a1))*0.5d0
            enddo
            enddo
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b)
            do b = 1,bStmp
            do a = 2,aStmp,2
                fS(:,b,a) = (fS(:,b,a-1) + fS(:,b,a+1))*0.5d0
            enddo
            enddo
            !$OMP END PARALLEL DO
            if (r(2).eq.1) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a)
                do a = 1,aStmp
                    fS(:,bStmp+1,a) = (fS(:,bStmp,a) + fS(:,1,a))*0.5d0
                enddo
                !$OMP END PARALLEL DO
            endif
            if (r(1).eq.1) then
                !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(b)
                do b = 1,bStmp
                    fS(:,b,aStmp+1) = (fS(:,b,aStmp) + fS(:,b,1))*0.5d0
                enddo
                !$OMP END PARALLEL DO
                if (r(2).eq.1) fS(:,bStmp+1,aStmp+1) = (fS(:,bStmp+1,aStmp) + fS(:,bStmp+1,1))*0.5d0
            endif
        endif
    end subroutine

    subroutine interpolate_tau(bF,aF,fF,bS,aS,fS)
        implicit none
        integer,intent(in):: aS,bS,aF,bF
        real(8),intent(in):: fF(bF,aF)
        real(8),intent(out):: fS(bS,aS) !fine grid values
        integer:: a,b,a1,b1,r(2),aStmp,bStmp
        r = 0
        bStmp = bS
        aStmp = aS
        if (mod(bS,2).eq.0) then
            bStmp = bS - 1
            r(2) = 1
        endif
        if (mod(aS,2).eq.0) then
            aStmp = aS - 1
            r(1) = 1
        endif
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b,a1,b1)
        do b  = 1,bStmp,2
            b1 = b / 2 + 1
        do a  = 1,aStmp,2
            a1 = a / 2 + 1
            fS(b,a) = fF(b1,a1)
            if(b.lt.bStmp) fS(b+1,a) = (fF(b1,a1) + fF(b1+1,a1))*0.5d0
        enddo
        enddo
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a,b)
        do b = 1,bStmp
        do a = 2,aStmp,2
            fS(b,a) = (fS(b,a-1) + fS(b,a+1))*0.5d0
        enddo
        enddo
        !$OMP END PARALLEL DO
        if (r(2).eq.1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(a)
            do a = 1,aStmp
                fS(bStmp+1,a) = (fS(bStmp,a) + fS(1,a))*0.5d0
            enddo
            !$OMP END PARALLEL DO
        endif
        if (r(1).eq.1) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(b)
            do b = 1,bStmp
                fS(b,aStmp+1) = (fS(b,aStmp) + fS(b,1))*0.5d0
            enddo
            !$OMP END PARALLEL DO
            if (r(2).eq.1) fS(bStmp+1,aStmp+1) = (fS(bStmp+1,aStmp) + fS(bStmp+1,1))*0.5d0
        endif
    end subroutine

    subroutine fIn_GridTransform(fIn,coeff,volumeForce,dh)! Dupius-Chopard method
        implicit none
        real(8),intent(in)::coeff,volumeForce(3),dh
        real(8),intent(inout):: fIn(0:lbmDim)
        real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim),den,uuu(3)
        ! Guo 2008 P97 6.1.8
        call cpt_macro()
        uSqr           = sum(uuu(1:3)**2)
        uxyz(0:lbmDim) = uuu(1) * ee(0:lbmDim,1) + uuu(2) * ee(0:lbmDim,2) + uuu(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * den * ( (1.0d0 - 1.5d0 * uSqr) + uxyz(0:lbmDim) * (3.0d0  + 4.5d0 * uxyz(0:lbmDim)) )
        fIn(0:lbmDim)  = fEq(0:lbmDim) + coeff * (fIn(0:lbmDim) - fEq(0:lbmDim))

        contains

        SUBROUTINE cpt_macro()
            implicit none
            den     = (SUM(fIn(0:lbmDim)))
            uuu(1)  = (SUM(fIn(0:lbmDim)*ee(0:lbmDim,1))+0.5d0*volumeForce(1)*dh)/den
            uuu(2)  = (SUM(fIn(0:lbmDim)*ee(0:lbmDim,2))+0.5d0*volumeForce(2)*dh)/den
            uuu(3)  = (SUM(fIn(0:lbmDim)*ee(0:lbmDim,3))+0.5d0*volumeForce(3)*dh)/den
        END SUBROUTINE
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
