module FluidDomain
    use ConstParams
    use FlowCondition
    implicit none
    private
    public:: LBMblks,LBMblksIndex,m_nblocks
    public:: read_fuild_blocks,allocate_fuild_memory_blocks,calculate_macro_quantities_blocks,calculate_turbulent_statistic_blocks,calculate_macro_quantities_iblock,initialise_fuild_blocks, &
             check_is_continue,add_volume_force_blocks,update_volume_force_blocks,write_flow_blocks,set_boundary_conditions_block,collision_block, &
             write_continue_blocks,streaming_block,computeFieldStat_blocks,clear_volume_force, &
             CompareBlocks,FindCarrierFluidBlock,halfwayBCset_block
    integer:: m_nblocks, m_npsize
    type :: LBMBlock
        integer :: ID,iCollidModel,offsetOutput,outputtype
        integer :: xDim,yDim,zDim
        real(8) :: dh,xmin,ymin,zmin,xmax,ymax,zmax
        integer :: BndConds(1:6),periodic_bc(3)
        real(8) :: params(1:10),tau,Omega,Omega2 ! single time, two time relaxation
        real(8) :: M_COLLID(0:lbmDim,0:lbmDim),M_FORCE(0:lbmDim,0:lbmDim) ! multiple time relaxation
        integer, allocatable :: OMPpartition(:),OMPparindex(:),OMPeid(:)
        real(8), allocatable :: OMPedge(:,:,:)
        real(8), allocatable :: fIn(:,:,:,:),uuu(:,:,:,:),force(:,:,:,:),den(:,:,:),uuu_ave(:,:,:,:)
        real(4), allocatable :: OUTtmp(:,:,:,:)
        real(8), allocatable :: fIn_hwx1(:,:,:),fIn_hwx2(:,:,:),fIn_hwy1(:,:,:),fIn_hwy2(:,:,:),fIn_hwz1(:,:,:),fIn_hwz2(:,:,:)
        real(8), allocatable :: fIn_Fx1t1(:,:,:),fIn_Fx1t2(:,:,:),fIn_Fx2t1(:,:,:),fIn_Fx2t2(:,:,:)
        real(8), allocatable :: fIn_Fy1t1(:,:,:),fIn_Fy1t2(:,:,:),fIn_Fy2t1(:,:,:),fIn_Fy2t2(:,:,:)
        real(8), allocatable :: fIn_Fz1t1(:,:,:),fIn_Fz1t2(:,:,:),fIn_Fz2t1(:,:,:),fIn_Fz2t2(:,:,:)
        real(8), allocatable :: tau_Fx1t1(:,:),tau_Fx1t2(:,:),tau_Fx2t1(:,:),tau_Fx2t2(:,:)
        real(8), allocatable :: tau_Fy1t1(:,:),tau_Fy1t2(:,:),tau_Fy2t1(:,:),tau_Fy2t2(:,:)
        real(8), allocatable :: tau_Fz1t1(:,:),tau_Fz1t2(:,:),tau_Fz2t1(:,:),tau_Fz2t2(:,:)
        real(8), allocatable :: tau_all(:,:,:)
        real(8) :: offsetMoveGrid(1:3),volumeForce(3)
        real(8) :: blktime
        integer, allocatable :: carriedBodies(:) ! 0 is the number of carried bodies
    contains
        procedure :: allocate_fluid => allocate_fluid_
        procedure :: Initialise => Initialise_
        procedure :: read_continue => read_continue_
        procedure :: calculate_macro_quantities => calculate_macro_quantities_
        procedure :: calculate_turbulent_statistic => calculate_turbulent_statistic_
        procedure :: collision => collision_
        procedure :: streaming => streaming_
        procedure :: set_boundary_conditions => set_boundary_conditions_
        procedure :: check_periodic_boundary => check_periodic_boundary_
        procedure :: update_volume_force => update_volume_force_
        procedure :: add_volume_force => add_volume_force_
        procedure :: write_flow => write_flow_
        procedure :: write_continue => write_continue_
        procedure :: ComputeFieldStat => ComputeFieldStat_
        procedure :: ResetVolumeForce => ResetVolumeForce_
        procedure :: halfwayBCset => halfwayBCset_
    end type LBMBlock
    type(LBMBlock), allocatable :: LBMblks(:)
    integer,allocatable:: LBMblksIndex(:)
    contains

    SUBROUTINE read_fuild_blocks(filename)
        implicit none
        character(LEN=40),intent(in):: filename
        integer:: iblock
        character(LEN=256):: buffer
        character(LEN=40):: keywordstr
        ! creat LBMblks and read parameters
        open(unit=111, file=filename, status='old', action='read')
        keywordstr = 'FluidBlocks'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*) m_nblocks
        allocate(LBMblks(m_nblocks),LBMblksIndex(m_nblocks))
        do iblock = 1,m_nblocks
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%ID,LBMblks(iblock)%iCollidModel,LBMblks(iblock)%offsetOutput,LBMblks(iblock)%outputtype
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%xDim,LBMblks(iblock)%yDim,LBMblks(iblock)%zDim
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%dh,LBMblks(iblock)%xmin,LBMblks(iblock)%ymin,LBMblks(iblock)%zmin
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%BndConds(1:6)
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%params(1:10)
            if(iblock.lt.m_nblocks) call readequal(111)
            LBMblksIndex(LBMblks(iblock)%ID) = iblock
            ! check bounds
            if (LBMblks(iblock)%xDim.gt.32767 .or. LBMblks(iblock)%yDim.gt.32767 .or. LBMblks(iblock)%zDim.gt.32767) then
                write(*,*) "Grid number exceeds 32767, please try to reduced the grid size."
                stop
            endif
            call LBMblks(iblock)%check_periodic_boundary()
            ! calculate the maximum coordinates
            LBMblks(iblock)%xmax = LBMblks(iblock)%xmin + LBMblks(iblock)%dh*(LBMblks(iblock)%xDim-1)
            LBMblks(iblock)%ymax = LBMblks(iblock)%ymin + LBMblks(iblock)%dh*(LBMblks(iblock)%yDim-1)
            LBMblks(iblock)%zmax = LBMblks(iblock)%zmin + LBMblks(iblock)%dh*(LBMblks(iblock)%zDim-1)
            if (LBMblks(iblock)%periodic_bc(1) .eq. 1) then
                LBMblks(iblock)%xmax = LBMblks(iblock)%xmax + LBMblks(iblock)%dh
            endif
            if (LBMblks(iblock)%periodic_bc(2) .eq. 1) then
                LBMblks(iblock)%ymax = LBMblks(iblock)%ymax + LBMblks(iblock)%dh
            endif
            if (LBMblks(iblock)%periodic_bc(3) .eq. 1) then
                LBMblks(iblock)%zmax = LBMblks(iblock)%zmax + LBMblks(iblock)%dh
            endif
        enddo
        close(111)
    END SUBROUTINE

    subroutine check_periodic_boundary_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer :: i
        this%periodic_bc = 0
        do i = 1,3
            if (this%BndConds(2*i-1) .eq. BCPeriodic .or. this%BndConds(2*i) .eq. BCPeriodic) then
                if (this%BndConds(2*i-1) .eq. this%BndConds(2*i)) then
                    this%periodic_bc(i) = 1
                else
                    write(*,*) 'Stop! Periodic boundaries must apper in pairs: ', this%BndConds(2*i-1), this%BndConds(2*i)
                    stop
                endif
            endif
        enddo
    end subroutine

    ! Check whether the calculation is continued
    SUBROUTINE check_is_continue(step,time,isContinue)
        implicit none
        integer,intent(out):: step
        real(8),intent(out):: time
        integer:: i,j,j2,x,y,z,iblock,indexs(1:6)
        integer:: isContinue,nblocks,index_tmp,idfile=13
        real(8):: xCoord,yCoord,zCoord,coeffs(1:3)
        type(LBMBlock), allocatable :: LBMblks_tmp(:)
        integer, allocatable :: sortdh(:)
        real(8), allocatable :: coor_max(:,:)
        logical:: alive
        inquire(file='./DatContinue/continue', exist=alive)
        if (isContinue==1 .and. alive) then
            write(*,'(A)') '========================================================='
            write(*,'(A)') '=================== Continue computing =================='
            write(*,'(A)') '========================================================='
            ! read continue file
            open(unit=idfile,file='./DatContinue/continue',form='unformatted',status='old',access='stream')
            read(idfile) nblocks,step,time
            allocate(LBMblks_tmp(nblocks),sortdh(nblocks),coor_max(1:3,nblocks))
            do iblock = 1,nblocks
                call LBMblks_tmp(iblock)%read_continue(idfile)
            enddo
            close(idfile)
            ! sort the blocks according to dh
            do iblock = 1,nblocks
                sortdh(iblock) = iblock
            enddo
            do i = 1, nblocks-1
            do j = i+1, nblocks
                if (LBMblks_tmp(i)%dh > LBMblks_tmp(j)%dh) then
                    index_tmp = sortdh(i)
                    sortdh(i) = sortdh(j)
                    sortdh(j) = index_tmp
                end if
            end do
            end do
            ! calculate the max coordinates of each blocks
            do iblock = 1,nblocks
                coor_max(1,iblock) = LBMblks_tmp(iblock)%zmin + (LBMblks_tmp(iblock)%zDim - 1) * LBMblks_tmp(iblock)%dh
                coor_max(2,iblock) = LBMblks_tmp(iblock)%ymin + (LBMblks_tmp(iblock)%yDim - 1) * LBMblks_tmp(iblock)%dh
                coor_max(3,iblock) = LBMblks_tmp(iblock)%xmin + (LBMblks_tmp(iblock)%xDim - 1) * LBMblks_tmp(iblock)%dh
            enddo
            ! interpolate from the continue file
            do i = 1,m_nblocks
                !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,j,j2,xCoord,yCoord,zCoord,indexs,coeffs)
                do x=1,LBMblks(i)%xDim
                    xCoord = LBMblks(i)%xmin + (x - 1) * LBMblks(i)%dh
                do y=1,LBMblks(i)%yDim
                    yCoord = LBMblks(i)%ymin + (y - 1) * LBMblks(i)%dh
                do z=1,LBMblks(i)%zDim
                    zCoord = LBMblks(i)%zmin + (z - 1) * LBMblks(i)%dh
                    ! judge in which fluid block
                    do j = 1,nblocks
                        j2 = sortdh(j)
                        if(zCoord .ge. LBMblks_tmp(j2)%zmin .and. zCoord .le. coor_max(3,j2) .and. &
                           yCoord .ge. LBMblks_tmp(j2)%ymin .and. yCoord .le. coor_max(2,j2) .and. &
                           xCoord .ge. LBMblks_tmp(j2)%xmin .and. xCoord .le. coor_max(1,j2)) then
                            ! calculate the interpolation coefficients of the 8 around points
                            coeffs(1) = (xCoord - LBMblks_tmp(j2)%xmin) / LBMblks_tmp(j2)%dh
                            coeffs(2) = (yCoord - LBMblks_tmp(j2)%ymin) / LBMblks_tmp(j2)%dh
                            coeffs(3) = (zCoord - LBMblks_tmp(j2)%zmin) / LBMblks_tmp(j2)%dh
                            indexs(1) = floor(coeffs(1))
                            indexs(3) = floor(coeffs(2))
                            indexs(5) = floor(coeffs(3))
                            coeffs(1) = coeffs(1) - dble(indexs(1))
                            coeffs(2) = coeffs(2) - dble(indexs(3))
                            coeffs(3) = coeffs(3) - dble(indexs(5))
                            ! calculate the indexs of the 8 around points
                            indexs(1) = indexs(1) + 1  ! x-1
                            indexs(2) = indexs(1) + 1  ! x
                            indexs(3) = indexs(3) + 1  ! y-1
                            indexs(4) = indexs(3) + 1  ! y
                            indexs(5) = indexs(5) + 1  ! z-1
                            indexs(6) = indexs(5) + 1  ! z
                            ! calculate the fIn of the 8 around points
                            LBMblks(i)%fIn(z,y,x,:) = LBMblks_tmp(j2)%fIn(indexs(5),indexs(3),indexs(1),:) * (1 - coeffs(3)) * (1 - coeffs(2)) * (1 - coeffs(1)) + &
                                                      LBMblks_tmp(j2)%fIn(indexs(5),indexs(3),indexs(2),:) * (1 - coeffs(3)) * (1 - coeffs(2)) * coeffs(1) + &
                                                      LBMblks_tmp(j2)%fIn(indexs(5),indexs(4),indexs(1),:) * (1 - coeffs(3)) *      coeffs(2)  * (1 - coeffs(1)) + &
                                                      LBMblks_tmp(j2)%fIn(indexs(5),indexs(4),indexs(2),:) * (1 - coeffs(3)) *      coeffs(2)  * coeffs(1) + &
                                                      LBMblks_tmp(j2)%fIn(indexs(6),indexs(3),indexs(1),:) *      coeffs(3)  * (1 - coeffs(2)) * (1 - coeffs(1)) + &
                                                      LBMblks_tmp(j2)%fIn(indexs(6),indexs(3),indexs(2),:) *      coeffs(3)  * (1 - coeffs(2)) * coeffs(1) + &
                                                      LBMblks_tmp(j2)%fIn(indexs(6),indexs(4),indexs(1),:) *      coeffs(3)  *      coeffs(2)  * (1 - coeffs(1)) + &
                                                      LBMblks_tmp(j2)%fIn(indexs(6),indexs(4),indexs(2),:) *      coeffs(3)  *      coeffs(2)  * coeffs(1)

                            exit
                        endif
                    enddo
                enddo
                enddo
                enddo
                !$OMP END PARALLEL DO
            enddo
            deallocate(LBMblks_tmp,sortdh,coor_max)
        else
            write(*,'(A)') '========================================================='
            write(*,'(A)') '====================== New computing ===================='
            write(*,'(A)') '========================================================='
        endif
    END SUBROUTINE

    !=================================================================================================
    ! Functions to iterate over each block
    SUBROUTINE allocate_fuild_memory_blocks(npsize)
        implicit none
        integer,intent(in):: npsize
        integer:: iblock
        m_npsize = npsize
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%allocate_fluid()
        enddo
    END SUBROUTINE

    SUBROUTINE clear_volume_force()
        implicit none
        integer:: iblock
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%ResetVolumeForce()
        enddo
    END SUBROUTINE

    SUBROUTINE initialise_fuild_blocks(time)
        implicit none
        real(8),intent(in):: time
        integer:: iblock
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%initialise(time)
        enddo
    END SUBROUTINE

    SUBROUTINE write_continue_blocks(step,time)
        implicit none
        real(8):: time
        integer:: step,i,iblock
        integer,parameter::nameLen=10,idfile=13
        character(len=nameLen):: fileName
        write(filename,'(I10)') nint(time/flow%Tref*1d5)
        fileName = adjustr(fileName)
        do  i=1,nameLen
            if(fileName(i:i)==' ') fileName(i:i)='0'
        enddo
        open(idfile,file='./DatContinue/continue'//trim(fileName),form='unformatted',access='stream')
        write(idfile) m_nblocks,step,time
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%write_continue(idfile)
        enddo
        close(idfile)
    END SUBROUTINE

    SUBROUTINE update_volume_force_blocks()
        implicit none
        integer:: iblock
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%update_volume_force()
        enddo
    END SUBROUTINE

    SUBROUTINE add_volume_force_blocks()
        implicit none
        integer:: iblock
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%add_volume_force()
        enddo
    END SUBROUTINE

    SUBROUTINE streaming_block(nblock)
        implicit none
        integer:: nblock
        call LBMblks(nblock)%streaming()
    END SUBROUTINE

    SUBROUTINE collision_block(nblock)
        implicit none
        integer:: nblock
        call LBMblks(nblock)%collision()
    END SUBROUTINE

    SUBROUTINE set_boundary_conditions_block(nblock)
        implicit none
        integer:: nblock
        call LBMblks(nblock)%set_boundary_conditions()
    END SUBROUTINE

    SUBROUTINE halfwayBCset_block(nblock)
        implicit none
        integer:: nblock
        call LBMblks(nblock)%halfwayBCset()
    END SUBROUTINE

    SUBROUTINE calculate_macro_quantities_blocks()
        implicit none
        integer:: iblock
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%calculate_macro_quantities()
        enddo
    END SUBROUTINE

    SUBROUTINE calculate_turbulent_statistic_blocks(step,step_s)
        implicit none
        integer:: iblock,step,step_s
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%calculate_turbulent_statistic(step,step_s)
        enddo
    END SUBROUTINE

    SUBROUTINE calculate_macro_quantities_iblock(iblock)
        implicit none
        integer:: iblock
        call LBMblks(iblock)%calculate_macro_quantities()
    END SUBROUTINE

    SUBROUTINE write_flow_blocks(time)
        implicit none
        real(8),intent(in):: time
        integer:: iblock
        real(8):: waittime,time_begine,time_end
        call get_now_time(time_begine)
        do iblock = 1,m_nblocks
            call mywait()
        enddo
        call get_now_time(time_end)
        waittime = time_end - time_begine
        if(waittime.gt.1.d-1) then
            write(*,'(A,F7.2,A)')'Waiting ', waittime, 's for previous outflow finishing.'
        endif
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%write_flow(time)
        enddo
    END SUBROUTINE

    SUBROUTINE computeFieldStat_blocks()
        implicit none
        integer:: iblock
        do iblock = 1,m_nblocks
            call LBMblks(iblock)%ComputeFieldStat()
        enddo
    END SUBROUTINE

    !=================================================================================================
    ! Functions for single block
    SUBROUTINE allocate_fluid_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: xmin,ymin,zmin,xmax,ymax,zmax
        ! allocate fluid memory
        allocate(this%fIn(this%zDim,this%yDim,this%xDim,0:lbmDim))
        allocate(this%uuu(this%zDim,this%yDim,this%xDim,1:3))
        allocate(this%force(this%zDim,this%yDim,this%xDim,1:3))
        allocate(this%den(this%zDim,this%yDim,this%xDim))
        if (this%outputtype .ge. 2) then
            allocate(this%uuu_ave(this%zDim,this%yDim,this%xDim,1:9))
        endif
        ! allocate output workspace
        xmin = 1 + this%offsetOutput
        ymin = 1 + this%offsetOutput
        zmin = 1 + this%offsetOutput
        xmax = this%xDim - this%offsetOutput
        ymax = this%yDim - this%offsetOutput
        zmax = this%zDim - this%offsetOutput
        if(this%outputtype .ge. 2) then
            allocate(this%outtmp(zmin:zmax,ymin:ymax,xmin:xmax,1:12))
        else
            allocate(this%outtmp(zmin:zmax,ymin:ymax,xmin:xmax,1:3))
        endif
        ! allocate mesh partition
        allocate(this%OMPpartition(1:m_npsize),this%OMPparindex(1:m_npsize+1),this%OMPeid(1:m_npsize))
        allocate(this%OMPedge(1:this%zDim,1:this%yDim, 1:m_npsize))
        call OMPPrePartition(this%xDim, m_npsize, this%OMPpartition, this%OMPparindex)

        contains

        SUBROUTINE OMPPrePartition(xDim, np, partition, parindex)
            implicit none
            integer:: np, xDim
            integer:: partition(1:np), parindex(1:np+1)
            integer:: psize, p, residual
            psize = xDim/np
            residual = xDim - psize * np
            parindex(1) = 1
            parindex(np+1) = xDim + 1
            do p=1,np
                if (p .gt. np-residual) then
                    partition(p) = psize + 1
                else
                    partition(p) = psize
                endif
                if (p .gt. 1) then
                    parindex(p) = parindex(p-1) + partition(p-1)
                endif
            enddo
        END SUBROUTINE OMPPrePartition
    END SUBROUTINE allocate_fluid_

    SUBROUTINE initialise_(this,time)
        implicit none
        class(LBMBlock),intent(inout) :: this
        real(8),intent(in):: time
        ! select the collision model
        call calculate_SRT_params()
        if(this%iCollidModel.eq.2) then
            call calculate_TRT_params(this%params(1))
        elseif(this%iCollidModel.eq.3) then
            call calculate_MRT_params()
        endif
        ! initialize flow information
        this%blktime = time
        call initialise_flow()

        contains

        SUBROUTINE calculate_SRT_params()
            implicit none
            this%tau = flow%nu/(this%dh*Cs2)+0.5d0
            this%Omega =  1.0d0 / this%tau
            allocate(this%tau_all(this%zDim,this%yDim,this%xDim))
            this%tau_all = this%tau
        END SUBROUTINE calculate_SRT_params

        SUBROUTINE calculate_TRT_params(lambda)
            ! lambda = (1/omegap - 0.5) (1/omegam - 0.5)
            implicit none
            real(8):: tmp, lambda
            tmp = (lambda * 4.d0 -1.d0) * this%Omega + 2.D0
            this%Omega2 = 2.d0*(2.d0 - this%Omega) / tmp
        END SUBROUTINE calculate_TRT_params

        SUBROUTINE calculate_MRT_params()
            implicit none
            integer:: I
            real(8):: M_MRT(0:lbmDim,0:lbmDim),M_MRTI(0:lbmDim,0:lbmDim),M(0:lbmDim,0:lbmDim)
            real(8):: S_D(0:lbmDim,0:lbmDim),S(0:lbmDim)
            ! calculate MRTM transformation matrix
            DO  I=0,lbmDim
                M_MRT(0,I)=1.0d0
                M_MRT(1,I)=19.0d0*SUM(ee(I,1:3)**2)-30.0d0
                M_MRT(2,I)=(21.0d0*SUM(ee(I,1:3)**2)**2-53.0d0*SUM(ee(I,1:3)**2)+24.0d0)/2.0d0

                M_MRT(3,I)=ee(I,1)
                M_MRT(5,I)=ee(I,2)
                M_MRT(7,I)=ee(I,3)

                M_MRT(4,I)=(5.0d0*SUM(ee(I,1:3)**2)-9.0d0)*ee(I,1)
                M_MRT(6,I)=(5.0d0*SUM(ee(I,1:3)**2)-9.0d0)*ee(I,2)
                M_MRT(8,I)=(5.0d0*SUM(ee(I,1:3)**2)-9.0d0)*ee(I,3)

                M_MRT(9,I)=3.0d0*ee(I,1)**2-SUM(ee(I,1:3)**2)
                M_MRT(10,I)=(3.0d0*SUM(ee(I,1:3)**2)-5.0d0)*(3.0d0*ee(I,1)**2-SUM(ee(I,1:3)**2))

                M_MRT(11,I)=ee(I,2)**2-ee(I,3)**2
                M_MRT(12,I)=(3.0d0*SUM(ee(I,1:3)**2)-5.0d0)*(ee(I,2)**2-ee(I,3)**2)

                M_MRT(13,I)=ee(I,1)*ee(I,2)
                M_MRT(14,I)=ee(I,2)*ee(I,3)
                M_MRT(15,I)=ee(I,3)*ee(I,1)

                M_MRT(16,I)=(ee(I,2)**2-ee(I,3)**2)*ee(I,1)
                M_MRT(17,I)=(ee(I,3)**2-ee(I,1)**2)*ee(I,2)
                M_MRT(18,I)=(ee(I,1)**2-ee(I,2)**2)*ee(I,3)
            ENDDO
            ! calculate the inverse matrix
            M_MRTI=TRANSPOSE(M_MRT)
            M=MATMUL(M_MRT,M_MRTI)
            DO    I=0,lbmDim
                M_MRTI(0:lbmDim,I)=M_MRTI(0:lbmDim,I)/M(I,I)
            ENDDO
            ! S(0:lbmDim)=Omega ! restore to SRT if S is Omega
            !              0   1  2  3  4  5  6  7  8      9     10      11     12      13        14         15      16  17  18
            S(0:lbmDim)=[  s0,s1,s2,s0,s4,s0,s4,s0,s4,this%Omega,s10,this%Omega,s10,this%Omega,this%Omega,this%Omega,s16,s16,s16]
            ! calculate MRTM collision matrix
            ! IM*S*M
            S_D(0:lbmDim,0:lbmDim)=0.0D0
            DO    i=0,lbmDim
                S_D(i,i)=S(i)
            ENDDO
            this%M_COLLID=MATMUL(MATMUL(M_MRTI,S_D),M_MRT)
            ! calculate MRTM body-force matrix
            ! IM*(I-0.5D0*S)*M=I-0.5*IM*S*M
            S_D(0:lbmDim,0:lbmDim)=0.0D0
            DO    i=0,lbmDim
                S_D(i,i)=1.0d0
            ENDDO
            this%M_FORCE=S_D-0.5*this%M_COLLID
        END SUBROUTINE calculate_MRT_params

        SUBROUTINE initialise_flow()
            implicit none
            real(8):: xCoord,yCoord,zCoord
            integer:: x, y, z
            ! calculating initial flow velocity and the distribution function
            do  x = 1, this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1, this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
            do  z = 1, this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                this%den(z,y,x) = flow%denIn
                call evaluate_velocity(this%blktime,zCoord,yCoord,xCoord,flow%uvwIn(1:SpaceDim),this%uuu(z,y,x,1:SpaceDim),flow%shearRateIn(1:3))
                call calculate_distribution_funcion(this%den(z,y,x),this%uuu(z,y,x,1:SpaceDim),this%fIn(z,y,x,0:lbmDim))
                if(this%outputtype .ge. 2) then
                    this%uuu_ave(z,y,x,:) = 0.0D0
                endif
            enddo
            enddo
            enddo
        END SUBROUTINE initialise_flow
    END SUBROUTINE initialise_

    ! subroutine  initDisturb_(this)
    !    implicit none
    !    class(LBMBlock), intent(inout) :: this
    !    real(8):: xCoord,yCoord,zCoord
    !    integer:: x, y, z
    !    do  x = 1, this%xDim
    !        xCoord = this%xmin + this%dh * (x - 1);
    !    do  y = 1, this%yDim
    !        yCoord = this%ymin + this%dh * (y - 1);
    !    do  z = 1, this%zDim
    !        zCoord = this%zmin + this%dh * (z - 1);
    !        this%uuu(z,y,x,1)=this%uuu(z,y,x,1)+flow%AmplInitDist(1)*flow%Uref*dsin(2.d0*pi*flow%waveInitDist*xCoord)
    !        this%uuu(z,y,x,2)=this%uuu(z,y,x,2)+flow%AmplInitDist(2)*flow%Uref*dsin(2.d0*pi*flow%waveInitDist*yCoord)
    !        this%uuu(z,y,x,3)=this%uuu(z,y,x,3)+flow%AmplInitDist(3)*flow%Uref*dsin(2.d0*pi*flow%waveInitDist*zCoord)
    !    enddo
    !    enddo
    !    enddo
    ! END SUBROUTINE


    SUBROUTINE halfwayBCset_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: x,y,z
        if (this%BndConds(1) .eq. BCstationary_Wall_halfway) then
            do y = 1,this%yDim
            do z = 1,this%zDim
                this%fIn_hwx1([1,7,9,11,13],z,y) = this%fIn(z,y,1,oppo([1,7,9,11,13]))
            enddo
            enddo
        endif
        if (this%BndConds(2) .eq. BCstationary_Wall_halfway) then
            do y = 1,this%yDim
            do z = 1,this%zDim
                this%fIn_hwx2([2,8,10,12,14],z,y) = this%fIn(z,y,this%xDim,oppo([2,8,10,12,14]))
            enddo
            enddo
        endif
        if (this%BndConds(3) .eq. BCstationary_Wall_halfway) then
            do x = 1,this%xDim
            do z = 1,this%zDim
                this%fIn_hwy1([3,7,8,15,17],z,x) = this%fIn(z,1,x,oppo([3,7,8,15,17]))
            enddo
            enddo
        endif
        if (this%BndConds(4) .eq. BCstationary_Wall_halfway) then
            do x = 1,this%xDim
            do z = 1,this%zDim
                this%fIn_hwy2([4,9,10,16,18],z,x) = this%fIn(z,this%yDim,x,oppo([4,9,10,16,18]))
            enddo
            enddo
        endif
        if (this%BndConds(5) .eq. BCstationary_Wall_halfway) then
            do x = 1,this%xDim
            do y = 1,this%yDim
                this%fIn_hwz1([5,11,12,15,16],y,x) = this%fIn(1,y,x,oppo([5,11,12,15,16]))
            enddo
            enddo
        endif
        if (this%BndConds(6) .eq. BCstationary_Wall_halfway) then
            do x = 1,this%xDim
            do y = 1,this%yDim
                this%fIn_hwz2([6,13,14,17,18],y,x) = this%fIn(this%zDim,y,x,oppo([6,13,14,17,18]))
            enddo
            enddo
        endif

    END SUBROUTINE

    SUBROUTINE set_boundary_conditions_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: x, y, z
        real(8):: xCoord,yCoord,zCoord,velocity(1:SpaceDim)
        real(8):: fEq(0:lbmDim),fEqi(0:lbmDim),fTmp(0:lbmDim)
        ! set the x direction (inlet) --------------------------------------------------------------------
        if (this%BndConds(1) .eq. BCEq_DirecletU) then
            ! equilibriun scheme
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,yCoord,this%xmin,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),this%fIn(z,y,1,0:lbmDim))
            enddo
            enddo
        elseif(this%BndConds(1) .eq. BCnEq_DirecletU)then
            ! non-equilibriun extrapoltion scheme
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,yCoord,this%xmin,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                ! equilibriun part
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
                ! non-equilibriun part
                call calculate_distribution_funcion(this%den(z,y,2),this%uuu(z,y,2,1:SpaceDim),fEqi(0:lbmDim))
                ! given equilibriun funciton
                this%fIn(z,y,1,[1,7,9,11,13]) = fEq([1,7,9,11,13]) + (this%fIn(z,y,2,[1,7,9,11,13]) - fEqi([1,7,9,11,13]))
            enddo
            enddo
        elseif(this%BndConds(1) .eq. BCorder1_Extrapolate)then
            this%fIn(:,:,1,[1,7,9,11,13]) = this%fIn(:,:,2,[1,7,9,11,13])
        elseif(this%BndConds(1) .eq. BCorder2_Extrapolate)then
            this%fIn(:,:,1,[1,7,9,11,13]) = 2.0*this%fIn(:,:,2,[1,7,9,11,13]) - this%fIn(:,:,3,[1,7,9,11,13])
        elseif(this%BndConds(1) .eq. BCstationary_Wall)then
            do  y = 1,this%yDim
            do  z = 1,this%zDim
                fTmp([1,7,9,11,13]) = this%fIn(z,y,1,oppo([1,7,9,11,13]))
                this%fIn(z,y,1,[1,7,9,11,13]) = fTmp([1,7,9,11,13])
            enddo
            enddo
        elseif(this%BndConds(1) .eq. BCstationary_Wall_halfway)then
            if (.not.allocated(this%fIn_hwx1)) then
                allocate(this%fIn_hwx1(0:lbmDim,this%zDim,this%yDim))
            else
                do  y = 1,this%yDim
                do  z = 1,this%zDim
                    this%fIn(z,y,1,[1,7,9,11,13]) = this%fIn_hwx1([1,7,9,11,13],z,y)
                enddo
                enddo
            endif
        elseif(this%BndConds(1) .eq. BCmoving_Wall)then
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,yCoord,this%xmin,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call evaluate_moving_wall(flow%denIn,velocity(1:SpaceDim),this%fIn(z,y,1,0:lbmDim),fTmp(0:lbmDim))
                this%fIn(z,y,1,[1,7,9,11,13]) = fTmp([1,7,9,11,13])
            enddo
            enddo
        elseif(this%BndConds(1) .eq. BCSymmetric)then
            do  y = 1,this%yDim
            do  z = 1,this%zDim
                fTmp([1,7,9,11,13]) = this%fIn(z,y,1,[2,8,10,12,14])
                this%fIn(z,y,1,[1,7,9,11,13]) = fTmp([1,7,9,11,13])
            enddo
            enddo
        elseif(this%BndConds(1) .eq. BCPeriodic .or. this%BndConds(1) .eq. BCfluid)then
            ! no need to set
        else
            stop 'inlet (xmin) has no such boundary condition'
        endif
        ! set the x direction (outlet) -------------------------------------------------------------------
        if (this%BndConds(2) .eq. BCEq_DirecletU) then
            ! equilibriun scheme
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,yCoord,this%xmax,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),this%fIn(z,y,this%xDim,0:lbmDim))
            enddo
            enddo
        elseif(this%BndConds(2) .eq. BCnEq_DirecletU)then
            ! non-equilibriun extrapoltion scheme
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,yCoord,this%xmax,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                ! equilibriun part
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
                ! non-equilibriun part
                call calculate_distribution_funcion(this%den(z,y,this%xDim-1),this%uuu(z,y,this%xDim-1,1:SpaceDim),fEqi(0:lbmDim))
                ! given equilibriun funciton
                this%fIn(z,y,this%xDim,[2,8,10,12,14]) = fEq([2,8,10,12,14]) + (this%fIn(z,y,this%xDim-1,[2,8,10,12,14]) - fEqi([2,8,10,12,14]))
            enddo
            enddo
        elseif(this%BndConds(2) .eq. BCorder1_Extrapolate)then
            this%fIn(:,:,this%xDim,[2,8,10,12,14]) = this%fIn(:,:,this%xDim-1,[2,8,10,12,14])
        elseif(this%BndConds(2) .eq. BCorder2_Extrapolate)then
            this%fIn(:,:,this%xDim,[2,8,10,12,14]) = 2.0*this%fIn(:,:,this%xDim-1,[2,8,10,12,14])-this%fIn(:,:,this%xDim-2,[2,8,10,12,14])
        elseif(this%BndConds(2) .eq. BCstationary_Wall)then
            do  y = 1,this%yDim
            do  z = 1,this%zDim
                fTmp([2,8,10,12,14]) = this%fIn(z,y,this%xDim,oppo([2,8,10,12,14]))
                this%fIn(z,y,this%xDim,[2,8,10,12,14]) = fTmp([2,8,10,12,14])
            enddo
            enddo
        elseif(this%BndConds(2) .eq. BCstationary_Wall_halfway)then
            if (.not.allocated(this%fIn_hwx2)) then
                allocate(this%fIn_hwx2(0:lbmDim,this%zDim,this%yDim))
            else
                do  y = 1,this%yDim
                do  z = 1,this%zDim
                    this%fIn(z,y,this%xDim,[2,8,10,12,14]) = this%fIn_hwx2([2,8,10,12,14],z,y)
                enddo
                enddo
            endif
        elseif(this%BndConds(2) .eq. BCmoving_Wall)then
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,yCoord,this%xmax,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call evaluate_moving_wall(flow%denIn,velocity(1:SpaceDim),this%fIn(z,y,this%xDim,0:lbmDim),fTmp(0:lbmDim))
                this%fIn(z,y,this%xDim,[2,8,10,12,14]) = fTmp([2,8,10,12,14])
            enddo
            enddo
        elseif(this%BndConds(2) .eq. BCSymmetric)then
            do  y = 1,this%yDim
            do  z = 1,this%zDim
                fTmp([2,8,10,12,14]) = this%fIn(z,y,this%xDim,[1,7,9,11,13])
                this%fIn(z,y,this%xDim,[2,8,10,12,14]) = fTmp([2,8,10,12,14])
            enddo
            enddo
        elseif(this%BndConds(2) .eq. BCPeriodic .or. this%BndConds(2) .eq. BCfluid)then
            ! no need to set
        else
            stop 'outlet (xmax) has no such boundary condition'
        endif
        ! set the y direction (lower) --------------------------------------------------------------------
        if (this%BndConds(3) .eq. BCEq_DirecletU) then
            ! equilibriun scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,this%ymin,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),this%fIn(z,1,x,0:lbmDim))
            enddo
            enddo
        elseif(this%BndConds(3) .eq. BCnEq_DirecletU)then
            ! non-equilibriun extrapoltion scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,this%ymin,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                ! equilibriun part
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
                ! non-equilibriun part
                call calculate_distribution_funcion(this%den(z,2,x),this%uuu(z,2,x,1:SpaceDim),fEqi(0:lbmDim))
                ! given equilibriun funciton
                this%fIn(z,1,x,[3,7,8,15,17]) = fEq([3,7,8,15,17]) + (this%fIn(z,2,x,[3,7,8,15,17]) - fEqi([3,7,8,15,17]))
            enddo
            enddo
        elseif(this%BndConds(3) .eq. BCorder1_Extrapolate)then
            this%fIn(:,1,:,[3,7,8,15,17]) = this%fIn(:,2,:,[3,7,8,15,17])
        elseif(this%BndConds(3) .eq. BCorder2_Extrapolate)then
            this%fIn(:,1,:,[3,7,8,15,17]) = 2.0*this%fIn(:,2,:,[3,7,8,15,17]) - this%fIn(:,3,:,[3,7,8,15,17])
        elseif(this%BndConds(3) .eq. BCstationary_Wall)then
            do  x = 1,this%xDim
            do  z = 1,this%zDim
                fTmp([3,7,8,15,17]) = this%fIn(z,1,x,oppo([3,7,8,15,17]))
                this%fIn(z,1,x,[3,7,8,15,17]) = fTmp([3,7,8,15,17])
            enddo
            enddo
        elseif(this%BndConds(3) .eq. BCstationary_Wall_halfway)then
            if (.not.allocated(this%fIn_hwy1)) then
                allocate(this%fIn_hwy1(0:lbmDim,this%zDim,this%xDim))
            else
                do  x = 1,this%xDim
                do  z = 1,this%zDim
                    this%fIn(z,1,x,[3,7,8,15,17]) = this%fIn_hwy1([3,7,8,15,17],z,x)
                enddo
                enddo
            endif
        elseif(this%BndConds(3) .eq. BCmoving_Wall)then
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,this%ymin,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call evaluate_moving_wall(flow%denIn,velocity(1:SpaceDim),this%fIn(z,1,x,0:lbmDim),fTmp(0:lbmDim))
                this%fIn(z,1,x,[3,7,8,15,17]) = fTmp([3,7,8,15,17])
            enddo
            enddo
        elseif(this%BndConds(3) .eq. BCSymmetric)then
            do  x = 1,this%xDim
            do  z = 1,this%zDim
                fTmp([3,7,8,15,17]) = this%fIn(z,1,x,[4,9,10,16,18])
                this%fIn(z,1,x,[3,7,8,15,17]) = fTmp([3,7,8,15,17])
            enddo
            enddo
        elseif(this%BndConds(3) .eq. BCPeriodic .or. this%BndConds(3) .eq. BCfluid)then
            ! no need to set
        else
            stop 'lower boundary (ymin) has no such boundary condition'
        endif
        ! set the y direction (higher) -------------------------------------------------------------------
        if (this%BndConds(4) .eq. BCEq_DirecletU) then
            ! equilibriun scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,this%ymax,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),this%fIn(z,this%yDim,x,0:lbmDim))
            enddo
            enddo
        elseif(this%BndConds(4) .eq. BCnEq_DirecletU)then
            ! non-equilibriun extrapoltion scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,this%ymax,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                ! equilibriun part
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
                ! non-equilibriun part
                call calculate_distribution_funcion(this%den(z,this%yDim-1,x),this%uuu(z,this%yDim-1,x,1:SpaceDim),fEqi(0:lbmDim))
                ! given equilibriun funciton
                this%fIn(z,this%yDim,x,[4,9,10,16,18]) = fEq([4,9,10,16,18]) + (this%fIn(z,this%yDim-1,x,[4,9,10,16,18]) - fEqi([4,9,10,16,18]))
            enddo
            enddo
        elseif(this%BndConds(4) .eq. BCorder1_Extrapolate)then
            this%fIn(:,this%yDim,:,[4,9,10,16,18]) = this%fIn(:,this%yDim-1,:,[4,9,10,16,18])
        elseif(this%BndConds(4) .eq. BCorder2_Extrapolate)then
            this%fIn(:,this%yDim,:,[4,9,10,16,18]) = 2.0*this%fIn(:,this%yDim-1,:,[4,9,10,16,18]) - this%fIn(:,this%yDim-2,:,[4,9,10,16,18])
        elseif(this%BndConds(4) .eq. BCstationary_Wall)then
            do  x = 1,this%xDim
            do  z = 1,this%zDim
                fTmp([4,9,10,16,18]) = this%fIn(z,this%yDim,x,oppo([4,9,10,16,18]))
                this%fIn(z,this%yDim,x,[4,9,10,16,18]) = fTmp([4,9,10,16,18])
            enddo
            enddo
        elseif(this%BndConds(4) .eq. BCstationary_Wall_halfway)then
            if (.not.allocated(this%fIn_hwy2)) then
                allocate(this%fIn_hwy2(0:lbmDim,this%zDim,this%xDim))
            else
                do  x = 1,this%xDim
                do  z = 1,this%zDim
                    this%fIn(z,this%yDim,x,[4,9,10,16,18]) = this%fIn_hwy2([4,9,10,16,18],z,x)
                enddo
                enddo
            endif
        elseif(this%BndConds(4) .eq. BCmoving_Wall)then
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  z = 1,this%zDim
                zCoord = this%zmin + this%dh * (z - 1);
                call evaluate_velocity(this%blktime,zCoord,this%ymax,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call evaluate_moving_wall(flow%denIn,velocity(1:SpaceDim),this%fIn(z,this%yDim,x,0:lbmDim),fTmp(0:lbmDim))
                this%fIn(z,this%yDim,x,[4,9,10,16,18]) = fTmp([4,9,10,16,18])
            enddo
            enddo
        elseif(this%BndConds(4) .eq. BCSymmetric)then
            do  x = 1,this%xDim
            do  z = 1,this%zDim
                fTmp([4,9,10,16,18]) = this%fIn(z,this%yDim,x,[3,7,8,15,17])
                this%fIn(z,this%yDim,x,[4,9,10,16,18]) = fTmp([4,9,10,16,18])
            enddo
            enddo
        elseif(this%BndConds(4) .eq. BCPeriodic .or. this%BndConds(4) .eq. BCfluid)then
            ! no need to set
        else
            stop 'higher boundary (ymax) has no such boundary condition'
        endif
        ! set the z direction (front) --------------------------------------------------------------------
        if (this%BndConds(5) .eq. BCEq_DirecletU) then
            ! equilibriun scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
                call evaluate_velocity(this%blktime,this%zmin,yCoord,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),this%fIn(1,y,x,0:lbmDim))
            enddo
            enddo
        elseif(this%BndConds(5) .eq. BCnEq_DirecletU)then
            ! non-equilibriun extrapoltion scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
                call evaluate_velocity(this%blktime,this%zmin,yCoord,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                ! equilibriun part
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
                ! non-equilibriun part
                call calculate_distribution_funcion(this%den(2,y,x),this%uuu(2,y,x,1:SpaceDim),fEqi(0:lbmDim))
                ! given equilibriun funciton
                this%fIn(1,y,x,[5,11,12,15,16]) = fEq([5,11,12,15,16]) + (this%fIn(2,y,x,[5,11,12,15,16]) - fEqi([5,11,12,15,16]))
            enddo
            enddo
        elseif(this%BndConds(5) .eq. BCorder1_Extrapolate)then
            this%fIn(1,:,:,[5,11,12,15,16]) = this%fIn(2,:,:,[5,11,12,15,16])
        elseif(this%BndConds(5) .eq. BCorder2_Extrapolate)then
            this%fIn(1,:,:,[5,11,12,15,16]) = 2.0*this%fIn(2,:,:,[5,11,12,15,16]) - this%fIn(3,:,:,[5,11,12,15,16])
        elseif(this%BndConds(5) .eq. BCstationary_Wall)then
            do  x = 1,this%xDim
            do  y = 1,this%yDim
                fTmp([5,11,12,15,16]) = this%fIn(1,y,x,oppo([5,11,12,15,16]))
                this%fIn(1,y,x,[5,11,12,15,16]) = fTmp([5,11,12,15,16])
            enddo
            enddo
        elseif(this%BndConds(5) .eq. BCstationary_Wall_halfway)then
            if (.not.allocated(this%fIn_hwz1)) then
                allocate(this%fIn_hwz1(0:lbmDim,this%yDim,this%xDim))
            else
                do  x = 1,this%xDim
                do  y = 1,this%yDim
                    this%fIn(1,y,x,[5,11,12,15,16]) = this%fIn_hwz1([5,11,12,15,16],y,x)
                enddo
                enddo
            endif
        elseif(this%BndConds(5) .eq. BCmoving_Wall)then
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
                call evaluate_velocity(this%blktime,this%zmin,yCoord,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call evaluate_moving_wall(flow%denIn,velocity(1:SpaceDim),this%fIn(1,y,x,0:lbmDim),fTmp(0:lbmDim))
                this%fIn(1,y,x,[5,11,12,15,16]) = fTmp([5,11,12,15,16])
            enddo
            enddo
        elseif(this%BndConds(5) .eq. BCSymmetric)then
            do  x = 1,this%xDim
            do  y = 1,this%yDim
                fTmp([5,11,12,15,16]) = this%fIn(1,y,x,[6,13,14,17,18])
                this%fIn(1,y,x,[5,11,12,15,16]) = fTmp([5,11,12,15,16])
            enddo
            enddo
        elseif(this%BndConds(5) .eq. BCPeriodic .or. this%BndConds(5) .eq. BCfluid)then
            ! no need to set
        else
            stop 'front boundary (zmin) has no such boundary condition'
        endif
        ! set the z direction (back) -------------------------------------------------------------------
        if (this%BndConds(6) .eq. BCEq_DirecletU) then
            ! equilibriun scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
                call evaluate_velocity(this%blktime,this%zmax,yCoord,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),this%fIn(this%zDim,y,x,0:lbmDim))
            enddo
            enddo
        elseif(this%BndConds(6) .eq. BCnEq_DirecletU)then
            ! non-equilibriun extrapoltion scheme
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
                call evaluate_velocity(this%blktime,this%zmax,yCoord,xCoord,flow%uvwIn,velocity(1:SpaceDim),flow%shearRateIn(1:3))
                ! equilibriun part
                call calculate_distribution_funcion(flow%denIn,velocity(1:SpaceDim),fEq(0:lbmDim))
                ! non-equilibriun part
                call calculate_distribution_funcion(this%den(this%zDim-1,y,x),this%uuu(this%zDim-1,y,x,1:3),fEqi(0:lbmDim))
                ! given equilibriun funciton
                this%fIn(this%zDim,y,x,[6,13,14,17,18]) = fEq([6,13,14,17,18]) + (this%fIn(this%zDim-1,y,x,[6,13,14,17,18]) - fEqi([6,13,14,17,18]))
            enddo
            enddo
        elseif(this%BndConds(6) .eq. BCorder1_Extrapolate)then
            this%fIn(this%zDim,:,:,[6,13,14,17,18]) = this%fIn(this%zDim-1,:,:,[6,13,14,17,18])
        elseif(this%BndConds(6) .eq. BCorder2_Extrapolate)then
            this%fIn(this%zDim,:,:,[6,13,14,17,18]) = 2.0*this%fIn(this%zDim-1,:,:,[6,13,14,17,18])-this%fIn(this%zDim-2,:,:,[6,13,14,17,18])
        elseif(this%BndConds(6) .eq. BCstationary_Wall)then
            do  x = 1,this%xDim
            do  y = 1,this%yDim
                fTmp([6,13,14,17,18]) = this%fIn(this%zDim,y,x,oppo([6,13,14,17,18]))
                this%fIn(this%zDim,y,x,[6,13,14,17,18]) = fTmp([6,13,14,17,18])
            enddo
            enddo
        elseif(this%BndConds(6) .eq. BCstationary_Wall_halfway)then
            if (.not.allocated(this%fIn_hwz2)) then
                allocate(this%fIn_hwz2(0:lbmDim,this%yDim,this%xDim))
            else
                do  x = 1,this%xDim
                do  y = 1,this%yDim
                    this%fIn(this%zDim,y,x,[6,13,14,17,18]) = this%fIn_hwz2([6,13,14,17,18],y,x)
                enddo
                enddo
            endif
        elseif(this%BndConds(6) .eq. BCmoving_Wall)then
            do  x = 1,this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1,this%yDim
                yCoord = this%ymin + this%dh * (y - 1);
                call evaluate_velocity(this%blktime,this%zmax,yCoord,xCoord,flow%uvwIn(1:SpaceDim),velocity(1:SpaceDim),flow%shearRateIn(1:3))
                call evaluate_moving_wall(flow%denIn,velocity(1:SpaceDim),this%fIn(this%zDim,y,x,0:lbmDim),fTmp(0:lbmDim))
                this%fIn(this%zDim,y,x,[6,13,14,17,18]) = fTmp([6,13,14,17,18])
            enddo
            enddo
        elseif(this%BndConds(6) .eq. BCSymmetric)then
            do  x = 1,this%xDim
            do  y = 1,this%yDim
                fTmp([6,13,14,17,18]) = this%fIn(this%zDim,y,x,[5,11,12,15,16])
                this%fIn(this%zDim,y,x,[6,13,14,17,18]) = fTmp([6,13,14,17,18])
            enddo
            enddo
        elseif(this%BndConds(6) .eq. BCPeriodic .or. this%BndConds(6) .eq. BCfluid)then
            ! no need to set
        else
            stop 'back boundary (zmax) has no such boundary condition'
        endif
    END SUBROUTINE

    SUBROUTINE calculate_macro_quantities_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer::x,y,z
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
        do  x = 1, this%xDim
        do  y = 1, this%yDim
        do  z = 1, this%zDim
            this%den(z,y,x  )  = (SUM(this%fIn(z,y,x,0:lbmDim)))
            this%uuu(z,y,x,1)  = (SUM(this%fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,1))+0.5d0*this%volumeForce(1)*this%dh)/this%den(z,y,x)
            this%uuu(z,y,x,2)  = (SUM(this%fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,2))+0.5d0*this%volumeForce(2)*this%dh)/this%den(z,y,x)
            this%uuu(z,y,x,3)  = (SUM(this%fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,3))+0.5d0*this%volumeForce(3)*this%dh)/this%den(z,y,x)
            !prs(z,y,x)   = Cs2*(den(z,y,x)-denIn)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE

    SUBROUTINE calculate_turbulent_statistic_(this,step,step_s)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer::x,y,z,step,step_s
        real(8)::invStep
        if(step .ge. step_s .and. this%outputtype .ge. 2) then
            invStep = 1 / real(step - step_s + 1)
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
            do  x = 1, this%xDim
            do  y = 1, this%yDim
            do  z = 1, this%zDim
                this%uuu_ave(z,y,x,1) = this%uuu_ave(z,y,x,1) * (1.d0 - invStep) + invStep * this%uuu(z,y,x,1)
                this%uuu_ave(z,y,x,2) = this%uuu_ave(z,y,x,2) * (1.d0 - invStep) + invStep * this%uuu(z,y,x,2)
                this%uuu_ave(z,y,x,3) = this%uuu_ave(z,y,x,3) * (1.d0 - invStep) + invStep * this%uuu(z,y,x,3)
                this%uuu_ave(z,y,x,4) = this%uuu_ave(z,y,x,4) * (1.d0 - invStep) + invStep * (this%uuu(z,y,x,1) - this%uuu_ave(z,y,x,1)) * (this%uuu(z,y,x,1) - this%uuu_ave(z,y,x,1))
                this%uuu_ave(z,y,x,5) = this%uuu_ave(z,y,x,5) * (1.d0 - invStep) + invStep * (this%uuu(z,y,x,2) - this%uuu_ave(z,y,x,2)) * (this%uuu(z,y,x,2) - this%uuu_ave(z,y,x,2))
                this%uuu_ave(z,y,x,6) = this%uuu_ave(z,y,x,6) * (1.d0 - invStep) + invStep * (this%uuu(z,y,x,3) - this%uuu_ave(z,y,x,3)) * (this%uuu(z,y,x,3) - this%uuu_ave(z,y,x,3))
                this%uuu_ave(z,y,x,7) = this%uuu_ave(z,y,x,7) * (1.d0 - invStep) + invStep * (this%uuu(z,y,x,1) - this%uuu_ave(z,y,x,1)) * (this%uuu(z,y,x,2) - this%uuu_ave(z,y,x,2))
                this%uuu_ave(z,y,x,8) = this%uuu_ave(z,y,x,8) * (1.d0 - invStep) + invStep * (this%uuu(z,y,x,1) - this%uuu_ave(z,y,x,1)) * (this%uuu(z,y,x,3) - this%uuu_ave(z,y,x,3))
                this%uuu_ave(z,y,x,9) = this%uuu_ave(z,y,x,9) * (1.d0 - invStep) + invStep * (this%uuu(z,y,x,2) - this%uuu_ave(z,y,x,2)) * (this%uuu(z,y,x,3) - this%uuu_ave(z,y,x,3))
            enddo
            enddo
            enddo
            !$OMP END PARALLEL DO
        endif
    END SUBROUTINE

    SUBROUTINE update_volume_force_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        this%volumeForce(1) = flow%volumeForceIn(1) + flow%volumeForceAmp * dsin(2.d0*pi*flow%volumeForceFreq*this%blktime + flow%volumeForcePhi/180.0d0*pi)
        this%volumeForce(2) = flow%volumeForceIn(2)
        this%volumeForce(3) = flow%volumeForceIn(3)
    END SUBROUTINE

    SUBROUTINE add_volume_force_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: x
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
        do x=1,this%xDim
            this%force(:,:,x,1) = this%force(:,:,x,1) + this%volumeForce(1)
            this%force(:,:,x,2) = this%force(:,:,x,2) + this%volumeForce(2)
            this%force(:,:,x,3) = this%force(:,:,x,3) + this%volumeForce(3)
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE

    SUBROUTINE ResetVolumeForce_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: x
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
        do x=1,this%xDim
            this%force(:,:,x,1) = 0.d0
            this%force(:,:,x,2) = 0.d0
            this%force(:,:,x,3) = 0.d0
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE

    SUBROUTINE collision_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim),Flb(0:lbmDim),dt3,omega,f_1
        integer:: x,y,z
        dt3 = 3.d0*this%dh
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,uSqr,uxyz,fEq,Flb,omega)
        do    x = 1, this%xDim
        do    y = 1, this%yDim
        do    z = 1, this%zDim
            uSqr           = sum(this%uuu(z,y,x,1:3)*this%uuu(z,y,x,1:3))
            uxyz(0:lbmDim) = this%uuu(z,y,x,1) * ee(0:lbmDim,1) + this%uuu(z,y,x,2) * ee(0:lbmDim,2)+this%uuu(z,y,x,3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * this%den(z,y,x) * ( (1.0d0 - 1.5d0 * uSqr) + uxyz(0:lbmDim) * (3.0d0  + 4.5d0 * uxyz(0:lbmDim)) ) - this%fIn(z,y,x,0:lbmDim)
            Flb(0:lbmDim)  = dt3*wt(0:lbmDim)*( &
                 (ee(0:lbmDim,1)-this%uuu(z,y,x,1)+3.d0*uxyz(0:lbmDim)*ee(0:lbmDim,1))*this%force(z,y,x,1) &
                +(ee(0:lbmDim,2)-this%uuu(z,y,x,2)+3.d0*uxyz(0:lbmDim)*ee(0:lbmDim,2))*this%force(z,y,x,2) &
                +(ee(0:lbmDim,3)-this%uuu(z,y,x,3)+3.d0*uxyz(0:lbmDim)*ee(0:lbmDim,3))*this%force(z,y,x,3))
            if(this%iCollidModel==1)then
                ! SRT collision
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + this%Omega*fEq(0:lbmDim) + (1.d0-0.5d0*this%Omega)*Flb(0:lbmDim)
            elseif(this%iCollidModel==2)then
                ! TRT collision
                fEq(0) = this%Omega*fEq(0) + (1.d0-0.5d0*this%Omega)*Flb(0)
                uxyz(positivedirs) = 0.5d0*this%Omega *(fEq(positivedirs)+fEq(negativedirs)) + (0.5d0-0.25d0* this%Omega)*(Flb(positivedirs)+Flb(negativedirs))
                uxyz(negativedirs) = 0.5d0*this%Omega2*(fEq(positivedirs)-fEq(negativedirs)) + (0.5d0-0.25d0*this%Omega2)*(Flb(positivedirs)-Flb(negativedirs))
                fEq(positivedirs) = uxyz(positivedirs) + uxyz(negativedirs)
                fEq(negativedirs) = uxyz(positivedirs) - uxyz(negativedirs)
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + fEq
            elseif(this%iCollidModel==3)then
                ! MRT collision
                this%fIn(z,y,x,0:lbmDim)=this%fIn(z,y,x,0:lbmDim)+MATMUL( this%M_COLLID(0:lbmDim,0:lbmDim), fEq(0:lbmDim) ) + MATMUL( this%M_FORCE(0:lbmDim,0:lbmDim),Flb(0:lbmDim))
            elseif(this%iCollidModel==11)then
                ! SRT collision with LES
                call smag(-fEq(0:lbmDim),this%den(z,y,x),x,y,z,omega)
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + omega*fEq(0:lbmDim) + (1.d0-0.5d0*omega)*Flb(0:lbmDim)
            elseif(this%iCollidModel==12)then
                ! Regularised SRT collision
                call RBGK(-fEq(0:lbmDim),f_1)
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + this%Omega*f_1 + (1.d0-0.5d0*omega)*Flb(0:lbmDim)
            elseif(this%iCollidModel==13)then
                ! SRT collision in ELBM
                call ELBM(fEq(0:lbmDim),this%fIn(z,y,x,0:lbmDim),omega)
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + omega*fEq(0:lbmDim) + (1.d0-0.5d0*omega)*Flb(0:lbmDim)
            elseif(this%iCollidModel==14)then
                ! WALE SRT collision
                call WALE(-fEq(0:lbmDim),this%den(z,y,x),x,y,z,omega)
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + omega*fEq(0:lbmDim) + (1.d0-0.5d0*omega)*Flb(0:lbmDim)
            elseif(this%iCollidModel==15)then
                ! Vremann SRT collision
                call vrem(x,y,z,omega)
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + omega*fEq(0:lbmDim) + (1.d0-0.5d0*omega)*Flb(0:lbmDim)
            endif
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO
        contains
        SUBROUTINE smag(fneq,rho,x0,y0,z0,omega0)
            implicit none
            real(8):: fneq(0:lbmDim),rho,omega0,tau_t
            real(8):: Q11,Q12,Q13,Q22,Q23,Q33
            real(8):: Q
            integer:: x0,y0,z0
            Q11 = fneq(1)+fneq(2)+fneq(7)+fneq(8)+fneq(9)+fneq(10)+fneq(11)+fneq(12)+fneq(13)+fneq(14)
            Q22 = fneq(3)+fneq(4)+fneq(7)+fneq(8)+fneq(9)+fneq(10)+fneq(15)+fneq(16)+fneq(17)+fneq(18)
            Q33 = fneq(5)+fneq(6)+fneq(11)+fneq(12)+fneq(13)+fneq(14)+fneq(15)+fneq(16)+fneq(17)+fneq(18)
            Q12 = fneq(7)-fneq(8)-fneq(9)+fneq(10)
            Q13 = fneq(11)-fneq(12)-fneq(13)+fneq(14)
            Q23 = fneq(15)-fneq(16)-fneq(17)+fneq(18)
            Q = Q11*Q11 + Q22*Q22 + Q33*Q33 + 2.d0*(Q12*Q12 + Q13*Q13 + Q23*Q23)
            tau_t = dsqrt(this%tau*this%tau + CsmagConst*dsqrt(Q)/rho)
            this%tau_all(z0,y0,x0) = 0.5d0 * (this%tau+tau_t)
            omega0 = 2.0d0 / (this%tau+tau_t)
        END SUBROUTINE
        SUBROUTINE RBGK(fneq,f_10)
            implicit none
            real(8):: fneq(0:lbmDim)
            real(8):: Q11,Q12,Q13,Q22,Q23,Q33,f_10
            Q11 = fneq(1)+fneq(2)+fneq(7)+fneq(8)+fneq(9)+fneq(10)+fneq(11)+fneq(12)+fneq(13)+fneq(14)
            Q22 = fneq(3)+fneq(4)+fneq(7)+fneq(8)+fneq(9)+fneq(10)+fneq(15)+fneq(16)+fneq(17)+fneq(18)
            Q33 = fneq(5)+fneq(6)+fneq(11)+fneq(12)+fneq(13)+fneq(14)+fneq(15)+fneq(16)+fneq(17)+fneq(18)
            Q12 = fneq(7)-fneq(8)-fneq(9)+fneq(10)
            Q13 = fneq(11)-fneq(12)-fneq(13)+fneq(14)
            Q23 = fneq(15)-fneq(16)-fneq(17)+fneq(18)
            f_10= 4.5d0*(10.0d0-Cs2)*(Q11+Q22+Q33)
        END SUBROUTINE
        SUBROUTINE ELBM(fneq,fnin,omega0)
            implicit none
            real(8):: fneq(0:lbmDim),fnin(0:lbmDim),omega0
            real(8):: x1(0:lbmDim),x2(0:lbmDim),x3(0:lbmDim)
            real(8):: a,b,c,alpha,beta
            x1(0:lbmDim) = fneq(0:lbmDim) / fnin(0:lbmDim)
            x2(0:lbmDim) = x1(0:lbmDim) * x1(0:lbmDim)
            ! x3(0:lbmDim) = x1(0:lbmDim) * x1(0:lbmDim) * x1(0:lbmDim) * (x1 < 0.0)

            a = sum(fnin(0:lbmDim) * x2(0:lbmDim))
            b = sum(fnin(0:lbmDim) * x3(0:lbmDim))
            c = sum(fnin(0:lbmDim) * (2.0d0 * x2(0:lbmDim) / (2.0d0 + x1(0:lbmDim))))

            alpha = (a - sqrt(a*a - 8.0d0*b*c)) / b / 2.0d0
            beta = flow%timeFlowDelta / (2.0d0 * this%tau + flow%timeFlowDelta)
            omega0 = alpha * beta
        END SUBROUTINE
        SUBROUTINE WALE(fneq,rho,x0,y0,z0,omega0)
            USE, INTRINSIC :: IEEE_ARITHMETIC
            implicit none
            real(8):: fneq(0:lbmDim),rho,omega0
            real(8):: Q11,Q12,Q13,Q22,Q23,Q33
            real(8):: S11,S12,S13,S22,S23,S33
            real(8):: O12,O13,O23,ox,oy,oz
            real(8):: SO11,SO12,SO13,SO22,SO23,SO33
            real(8):: Q,S,O,SO,OP,SdSd
            real(8):: invdh,tau__
            integer:: x0,y0,z0
            invdh = this%dh
            Q11 = fneq(1)+fneq(2)+fneq(7)+fneq(8)+fneq(9)+fneq(10)+fneq(11)+fneq(12)+fneq(13)+fneq(14)
            Q22 = fneq(3)+fneq(4)+fneq(7)+fneq(8)+fneq(9)+fneq(10)+fneq(15)+fneq(16)+fneq(17)+fneq(18)
            Q33 = fneq(5)+fneq(6)+fneq(11)+fneq(12)+fneq(13)+fneq(14)+fneq(15)+fneq(16)+fneq(17)+fneq(18)
            Q12 = fneq(7)-fneq(8)-fneq(9)+fneq(10)
            Q13 = fneq(11)-fneq(12)-fneq(13)+fneq(14)
            Q23 = fneq(15)-fneq(16)-fneq(17)+fneq(18)
            Q = Q11*Q11 + Q22*Q22 + Q33*Q33 + 2.d0*(Q12*Q12 + Q13*Q13 + Q23*Q23)

            tau__ = this%tau_all(z0,y0,x0)
            S11 = -1.5d0*invdh*Q11/(rho*tau__)
            S22 = -1.5d0*invdh*Q22/(rho*tau__)
            S33 = -1.5d0*invdh*Q33/(rho*tau__)
            S12 = -1.5d0*invdh*Q12/(rho*tau__)
            S13 = -1.5d0*invdh*Q13/(rho*tau__)
            S23 = -1.5d0*invdh*Q23/(rho*tau__)
            S = S11*S11 + S22*S22 + S33*S33 + 2.d0*(S12*S12 + S13*S13 + S23*S23)

            ox=0.0d0
            oy=0.0d0
            oz=0.0d0
            ! compute omega_x
            if (y0.gt.1.and.y0.lt.this%yDim) then
                ox = center_diff(this%uuu(z0,y0+1,x0,3),this%uuu(z0,y0-1,x0,3),invdh)
            elseif (y0.eq.1) then
                ox = onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0+1,x0,3),this%uuu(z0,y0+2,x0,3),invdh)
            else
                ox = onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0-1,x0,3),this%uuu(z0,y0-2,x0,3),invdh)
            endif
            if (z0.gt.1.and.z0.lt.this%zDim) then
                ox = 0.5d0*(ox - center_diff(this%uuu(z0+1,y0,x0,2),this%uuu(z0-1,y0,x0,2),invdh))
            elseif (z0.eq.1) then
                ox = 0.5d0*(ox - onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0+1,y0,x0,2),this%uuu(z0+2,y0,x0,2),invdh))
            else
                ox = 0.5d0*(ox - onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0-1,y0,x0,2),this%uuu(z0-2,y0,x0,2),invdh))
            endif
            ! compute omega_y
            if (z0.gt.1.and.z0.lt.this%zDim) then
                oy = center_diff(this%uuu(z0+1,y0,x0,1),this%uuu(z0-1,y0,x0,1),invdh)
            elseif (z0.eq.1) then
                oy = onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0+1,y0,x0,1),this%uuu(z0+2,y0,x0,1),invdh)
            else
                oy = onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0-1,y0,x0,1),this%uuu(z0-2,y0,x0,1),invdh)
            endif
            if (x0.gt.1.and.x0.lt.this%xDim) then
                oy = 0.5d0*(oy - center_diff(this%uuu(z0,y0,x0+1,3),this%uuu(z0,y0,x0-1,3),invdh))
            elseif (x0.eq.1) then
                oy = 0.5d0*(oy - onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0,x0+1,3),this%uuu(z0,y0,x0+2,3),invdh))
            else
                oy = 0.5d0*(oy - onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0,x0-1,3),this%uuu(z0,y0,x0-2,3),invdh))
            endif
            ! compute omega_y
            if (x0.gt.1.and.x0.lt.this%xDim) then
                oz = center_diff(this%uuu(z0,y0,x0+1,2),this%uuu(z0,y0,x0-1,2),invdh)
            elseif (x0.eq.1) then
                oz = onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0,y0,x0+1,2),this%uuu(z0,y0,x0+2,2),invdh)
            else
                oz = onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0,y0,x0-1,2),this%uuu(z0,y0,x0-2,2),invdh)
            endif
            if (y0.gt.1.and.y0.lt.this%yDim) then
                oz = 0.5d0*(oz - center_diff(this%uuu(z0,y0+1,x0,1),this%uuu(z0,y0-1,x0,1),invdh))
            elseif (y0.eq.1) then
                oz = 0.5d0*(oz - onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0,y0+1,x0,1),this%uuu(z0,y0+2,x0,1),invdh))
            else
                oz = 0.5d0*(oz - onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0,y0-1,x0,1),this%uuu(z0,y0-2,x0,1),invdh))
            endif
            O12 = -0.5d0*oz
            O13 =  0.5d0*oy
            O23 = -0.5d0*ox
            ! the boundary need unilateral interpolation

            O = 2.0d0*(O12*O12+O23*O23+O13*O13)

            SO11 =-(0.0d0           + S11*S11*O12*O12 + S11*S11*O13*O13 + &
                    0.0d0           + S12*S12*O12*O12 + S12*S12*O13*O13 + &
                    0.0d0           + S13*S13*O12*O12 + S13*S13*O13*O13)
            SO22 =-(S12*S12*O12*O12 + 0.0d0           + S12*S12*O23*O23 + &
                    S22*S22*O12*O12 + 0.0d0           + S22*S22*O23*O23 + &
                    S23*S23*O12*O12 + 0.0d0           + S23*S23*O23*O23)
            SO33 =-(S13*S13*O13*O13 + S13*S13*O23*O23 + 0.0d0           + &
                    S23*S23*O13*O13 + S23*S23*O23*O23 + 0.0d0           + &
                    S33*S33*O13*O13 + S33*S33*O23*O23 + 0.0d0          )
            SO12 =-(0.0d0           + 0.0d0           + S11*S12*O13*O23 + &
                    0.0d0           + 0.0d0           + S12*S22*O13*O23 + &
                    0.0d0           + 0.0d0           + S13*S23*O13*O23)
            SO13 = (0.0d0           + S11*S13*O12*O23 + 0.0d0           + &
                    0.0d0           + S12*S23*O12*O23 + 0.0d0           + &
                    0.0d0           + S13*S33*O12*O23 + 0.0d0          )
            SO23 =-(S12*S13*O12*O13 + 0.0d0           + 0.0d0           + &
                    S22*S23*O12*O13 + 0.0d0           + 0.0d0           + &
                    S23*S33*O12*O13 + 0.0d0           + 0.0d0          )
            SO = SO11 + SO22 + SO33 + 2.0d0*(SO12 + SO13 + SO23)

            SdSd = (S*S+O*O)/6.0d0+2.0d0*S*O/3.0d0+2.0d0*SO
            OP = SdSd**1.5d0/(S**2.5d0+SdSd**1.25d0)
            if ((.not. IEEE_IS_FINITE(OP)).or.(OP.lt.0.0d0)) then
                OP = 0.0d0
            endif

            tau__ = (flow%nu+CWALEConst*OP*this%dh*this%dh)/(this%dh*Cs2)+0.5d0
            this%tau_all(z0,y0,x0) = tau__
            omega0 = 1.0d0 / (tau__)
        END SUBROUTINE
        FUNCTION center_diff(g1,g2,invdx)
            implicit none
            real(8):: center_diff,g1,g2,invdx
            center_diff = (g1 - g2)*invdx
        END FUNCTION
        FUNCTION onesid_diff(g1,g2,g3,invdx)
            implicit none
            real(8):: onesid_diff,g1,g2,g3,invdx
            onesid_diff = (-3.0d0*g1 + 4.0d0*g2 - g3)*invdx
        END FUNCTION
        SUBROUTINE vrem(x0,y0,z0,omega0)
            USE, INTRINSIC :: IEEE_ARITHMETIC
            implicit none
            real(8):: omega0
            real(8):: a11,a12,a13,a21,a22,a23,a31,a32,a33
            real(8):: b11,b12,b13,b21,b22,b23,b31,b32,b33
            real(8):: aa,bb,OP
            real(8):: invdh,tau__
            integer:: x0,y0,z0
            invdh = this%dh
            if (x0.gt.1.and.x0.lt.this%xDim) then
                a11 = 0.5d0*center_diff(this%uuu(z0,y0,x0+1,1),this%uuu(z0,y0,x0-1,1),invdh)
                a21 = 0.5d0*center_diff(this%uuu(z0,y0,x0+1,2),this%uuu(z0,y0,x0-1,2),invdh)
                a31 = 0.5d0*center_diff(this%uuu(z0,y0,x0+1,3),this%uuu(z0,y0,x0-1,3),invdh)
            elseif (x0.eq.1) then
                a11 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0,y0,x0+1,1),this%uuu(z0,y0,x0+2,1),invdh)
                a21 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0,y0,x0+1,2),this%uuu(z0,y0,x0+2,2),invdh)
                a31 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0,x0+1,3),this%uuu(z0,y0,x0+2,3),invdh)
            else
                a11 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0,y0,x0-1,1),this%uuu(z0,y0,x0-2,1),invdh)
                a21 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0,y0,x0-1,2),this%uuu(z0,y0,x0-2,2),invdh)
                a31 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0,x0-1,3),this%uuu(z0,y0,x0-2,3),invdh)
            endif
            if (y0.gt.1.and.y0.lt.this%yDim) then
                a12 = 0.5d0*center_diff(this%uuu(z0,y0+1,x0,1),this%uuu(z0,y0-1,x0,1),invdh)
                a22 = 0.5d0*center_diff(this%uuu(z0,y0+1,x0,2),this%uuu(z0,y0-1,x0,2),invdh)
                a32 = 0.5d0*center_diff(this%uuu(z0,y0+1,x0,3),this%uuu(z0,y0-1,x0,3),invdh)
            elseif (y0.eq.1) then
                a12 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0,y0+1,x0,1),this%uuu(z0,y0+2,x0,1),invdh)
                a22 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0,y0+1,x0,2),this%uuu(z0,y0+2,x0,2),invdh)
                a32 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0+1,x0,3),this%uuu(z0,y0+2,x0,3),invdh)
            else
                a12 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0,y0-1,x0,1),this%uuu(z0,y0-2,x0,1),invdh)
                a22 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0,y0-1,x0,2),this%uuu(z0,y0-2,x0,2),invdh)
                a32 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0,y0-1,x0,3),this%uuu(z0,y0-2,x0,3),invdh)
            endif
            if (z0.gt.1.and.z0.lt.this%zDim) then
                a13 = 0.5d0*center_diff(this%uuu(z0+1,y0,x0,1),this%uuu(z0-1,y0,x0,1),invdh)
                a23 = 0.5d0*center_diff(this%uuu(z0+1,y0,x0,2),this%uuu(z0-1,y0,x0,2),invdh)
                a33 = 0.5d0*center_diff(this%uuu(z0+1,y0,x0,3),this%uuu(z0-1,y0,x0,3),invdh)
            elseif (z0.eq.1) then
                a13 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0+1,y0,x0,1),this%uuu(z0+2,y0,x0,1),invdh)
                a23 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0+1,y0,x0,2),this%uuu(z0+2,y0,x0,2),invdh)
                a33 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0+1,y0,x0,3),this%uuu(z0+2,y0,x0,3),invdh)
            else
                a13 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,1),this%uuu(z0-1,y0,x0,1),this%uuu(z0-2,y0,x0,1),invdh)
                a23 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,2),this%uuu(z0-1,y0,x0,2),this%uuu(z0-2,y0,x0,2),invdh)
                a33 = 0.5d0*onesid_diff(this%uuu(z0,y0,x0,3),this%uuu(z0-1,y0,x0,3),this%uuu(z0-2,y0,x0,3),invdh)
            endif
            b11 = a11*a11
            b12 = a12*a12
            b13 = a13*a13
            b21 = a21*a21
            b22 = a22*a22
            b23 = a23*a23
            b31 = a31*a31
            b32 = a32*a32
            b33 = a33*a33
            aa = b11+b12+b13+b21+b22+b23+b31+b32+b33
            bb = (b11+b12+b13)*(b21+b22+b23)- &
                 (a11*a21+a12*a22+a13*a23)*(a11*a21+a12*a22+a13*a23)+ &
                 (b11+b12+b13)*(b31+b32+b33)- &
                 (a11*a31+a12*a32+a13*a33)*(a11*a31+a12*a32+a13*a33)+ &
                 (b21+b22+b23)*(b31+b32+b33)- &
                 (a21*a31+a22*a32+a23*a33)*(a21*a31+a22*a32+a23*a33)
            OP = dsqrt(bb/aa)
            if (.not. IEEE_IS_FINITE(OP)) then
                OP = 0.0d0
            endif
            tau__ = (flow%nu+CvremConst*OP*this%dh*this%dh)/(this%dh*Cs2)+0.5d0
            this%tau_all(z0,y0,x0) = tau__
            omega0 = 1.0d0 / (tau__)
        END SUBROUTINE
    END SUBROUTINE collision_

    !0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    !0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0
    !0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1
    !0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1
    SUBROUTINE streaming_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: i
        do  i=0,lbmDim
            call swapzy(this%fIn, ee(i,3), ee(i,2), i, this%zDim, this%yDim, this%xDim, lbmDim)
            call swapx(this%fIn, ee(i,1), i, this%zDim, this%yDim, this%xDim, lbmDim, this%OMPparindex, this%OMPeid, this%OMPedge)
        enddo

        contains

        SUBROUTINE swapzy(f, dz, dy, i, zDim, yDim, xDim, lbmDim)
            implicit none
            integer, intent(in):: dz, dy, i, zDim, yDim, xDim, lbmDim
            real(8), intent(inout):: f(1:zDim,1:yDim,1:xDim,0:lbmDim)
            integer:: z, y, x
            real(8):: temp, tmpz(1:zDim)

            if(dz.eq.0 .and. dy.eq.0) return

            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,temp,tmpz)
            do  x = 1, xDim
                if(dz.eq.1) then
                    do y=1, yDim
                        temp = f(zDim, y, x, i)
                        do z=zDim, 2, -1
                            f(z,y,x, i)=f(z-1,y,x, i)
                        enddo
                        f(1,y,x,i) = temp
                    enddo
                elseif(dz.eq.-1) then
                    do y=1, yDim
                        temp = f(1, y, x, i)
                        do z=1, zDim-1
                            f(z,y,x, i)=f(z+1,y,x, i)
                        enddo
                        f(zDim,y,x,i) = temp
                    enddo
                endif
                if(dy.eq.1) then
                    tmpz = f(:, yDim, x, i)
                    do y=yDim, 2, -1
                        f(:,y,x, i)=f(:,y-1,x, i)
                    enddo
                    f(:,1,x,i) = tmpz
                elseif(dy.eq.-1) then
                    tmpz = f(:, 1, x, i)
                    do y=1, yDim-1
                        f(:,y,x, i)=f(:,y+1,x, i)
                    enddo
                    f(:,yDim,x,i) = tmpz
                endif
            enddo
            !$OMP END PARALLEL DO
        END SUBROUTINE

        SUBROUTINE swapx(f, dx, i, zDim, yDim, xDim, lbmDim, parindex, eid, edge)
            implicit none
            integer, intent(in):: dx, i, zDim, yDim, xDim, lbmDim, parindex(:)
            integer, intent(out):: eid(:)
            real(8), intent(out):: edge(:,:,:)
            real(8), intent(inout):: f(1:zDim, 1:yDim,1:xDim,0:lbmDim)
            integer:: p

            if(dx.eq.0) return

            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
            do  p = 1,m_npsize
                if(dx .eq. -1) then
                    call swapxwAtom(f, edge(:,:,p), eid(p), i, zDim, yDim, xDim, lbmDim, parindex(p), parindex(p+1)-1)
                elseif(dx.eq.1) then
                    call swapxeAtom(f, edge(:,:,p), eid(p), i, zDim, yDim, xDim, lbmDim, parindex(p), parindex(p+1)-1)
                endif
            enddo
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
            do  p = 1,m_npsize
                f(:,:,eid(p),i) = edge(:,:,p)
            enddo
            !$OMP END PARALLEL DO
        END SUBROUTINE

        SUBROUTINE swapxeAtom(f, edge, eid, i, zDim, yDim, xDim, lbmDim, xbgn, xend)
            implicit none
            integer, intent(in):: i, zDim, yDim, xDim, lbmDim, xbgn, xend
            real(8), intent(inout):: f(1:zDim,1:yDim,1:xDim,0:lbmDim)
            real(8), intent(out):: edge(1:zDim,1:yDim)
            integer, intent(out):: eid
            integer:: x
            eid = xend+1
            if(eid .eq. xDim+1) eid = 1
            edge = f(:,:,xend,i)
            do  x = xend,xbgn+1,-1
                f(:,:,x,i) = f(:,:,x-1,i)
            enddo
        endsubroutine

        SUBROUTINE swapxwAtom(f, edge, eid, i, zDim, yDim, xDim, lbmDim, xbgn, xend)
            implicit none
            integer, intent(in):: i, zDim, yDim, xDim, lbmDim, xbgn, xend
            real(8), intent(inout):: f(1:zDim,1:yDim,1:xDim,0:lbmDim)
            real(8), intent(out):: edge(1:zDim,1:yDim)
            integer, intent(out):: eid
            integer:: x
            eid = xbgn-1
            if(eid .eq. 0) eid = xDim
            edge = f(:,:,xbgn,i)
            do  x = xbgn, xend-1
                f(:,:,x,i) = f(:,:,x+1,i)
            enddo
        end subroutine
    END SUBROUTINE streaming_


    SUBROUTINE write_flow_(this,time)
        implicit none
        class(LBMBlock), intent(inout) :: this
        real(8), intent(in):: time
        integer:: x,y,z,pid,i
        real(8)::xmin,ymin,zmin
        integer::nxs,nxe,nys,nye,nzs,nze
        integer,parameter::nameLen=10,blockLen=3,idfile=100
        character (LEN=nameLen):: fileName
        character (LEN=blockLen):: blockName
        real(8):: invUref,invUrefs
        if(this%outputtype .lt. 1) return
        nxs  = 1 + this%offsetOutput
        nys  = 1 + this%offsetOutput
        nzs  = 1 + this%offsetOutput
        nxe  = this%xDim - this%offsetOutput
        nye  = this%yDim - this%offsetOutput
        nze  = this%zDim - this%offsetOutput
        xmin = this%xmin + this%offsetOutput * this%dh
        ymin = this%ymin + this%offsetOutput * this%dh
        zmin = this%zmin + this%offsetOutput * this%dh

        invUref  = 1.d0/flow%Uref
        invUrefs = 1.d0/flow%Uref/flow%Uref

        if(this%outputtype .ne. 2) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
            do x=nxs, nxe
                do y=nys, nye
                    do z=nzs, nze
                        this%outtmp(z,y,x,1) = this%uuu(z,y,x,1)*invUref !u
                        this%outtmp(z,y,x,2) = this%uuu(z,y,x,2)*invUref !v
                        this%outtmp(z,y,x,3) = this%uuu(z,y,x,3)*invUref !w
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        if(this%outputtype .ge. 2) then
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
            do x=nxs, nxe
                do y=nys, nye
                    do z=nzs, nze
                        this%outtmp(z,y,x,4) = this%uuu_ave(z,y,x,1)*invUref  !<u>
                        this%outtmp(z,y,x,5) = this%uuu_ave(z,y,x,2)*invUref  !<v>
                        this%outtmp(z,y,x,6) = this%uuu_ave(z,y,x,3)*invUref  !<w>
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
            do x=nxs, nxe
                do y=nys, nye
                    do z=nzs, nze
                        this%outtmp(z,y,x,7) = this%uuu_ave(z,y,x,4)*invUrefs  !<uu>
                        this%outtmp(z,y,x,8) = this%uuu_ave(z,y,x,5)*invUrefs  !<vv>
                        this%outtmp(z,y,x,9) = this%uuu_ave(z,y,x,6)*invUrefs  !<ww>
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
            do x=nxs, nxe
                do y=nys, nye
                    do z=nzs, nze
                        this%outtmp(z,y,x,10) = this%uuu_ave(z,y,x,7)*invUrefs  !<uv>
                        this%outtmp(z,y,x,11) = this%uuu_ave(z,y,x,8)*invUrefs  !<uw>
                        this%outtmp(z,y,x,12) = this%uuu_ave(z,y,x,9)*invUrefs  !<vw>
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO
        endif
        call myfork(pid)
        if(pid.eq.0) then
            if(this%outputtype .ne. 2) then
                write(fileName,'(I10)') nint(time/flow%Tref*1d5)
                fileName = adjustr(fileName)
                do  i=1,nameLen
                    if(fileName(i:i)==' ')fileName(i:i)='0'
                enddo
                write(blockName,'(I3)') this%ID
                blockName = adjustr(blockName)
                do  i=1,blockLen
                    if(blockName(i:i)==' ')blockName(i:i)='0'
                enddo
                open(idfile,file='./DatFlow/Flow'//trim(fileName)//'_b'//blockName,form='unformatted',access='stream')
                nxe = nxe-nxs+1
                nye = nye-nys+1
                nze = nze-nzs+1
                WRITE(idfile) nxe,nye,nze,this%ID
                WRITE(idfile) xmin,ymin,zmin,this%dh
                write(idfile) this%outtmp(:,:,:,1),this%outtmp(:,:,:,2),this%outtmp(:,:,:,3)
                close(idfile)
            endif
            if(this%outputtype .ge. 2) then
                open(idfile,file='./DatFlow/MeanFlow_b'//blockName,form='unformatted',access='stream')
                WRITE(idfile) nxe,nye,nze,this%ID
                WRITE(idfile) xmin,ymin,zmin,this%dh
                write(idfile) this%outtmp(:,:,:,4),this%outtmp(:,:,:,5),this%outtmp(:,:,:,6)
                write(idfile) this%outtmp(:,:,:,7),this%outtmp(:,:,:,8),this%outtmp(:,:,:,9)
                write(idfile) this%outtmp(:,:,:,10),this%outtmp(:,:,:,11),this%outtmp(:,:,:,12)
                close(idfile)
            endif
            call myexit(0)
        endif
    END SUBROUTINE write_flow_

    subroutine ComputeFieldStat_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: x,y,z,i
        real(8):: invUref, uLinfty(1:3), uL2(1:3), temp
        invUref = 1.d0/flow%Uref
        uLinfty = -1.d0
        uL2 = 0.d0
        do i=1,3
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,temp) &
            !$OMP reduction (+: uL2) reduction(max: uLinfty)
            do x=1, this%xDim
                do y=1, this%yDim
                    do z=1, this%zDim
                        temp = dabs(this%uuu(z,y,x,i)*invUref)
                        uL2(i) = uL2(i) + temp * temp
                        if(temp.gt.uLinfty(i)) uLinfty(i) = temp
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO
            uL2(i) = dsqrt(uL2(i) / (dble(this%xDim) * dble(this%yDim) * dble(this%zDim)))
        enddo
        write(*,'(A,F18.12)')' FIELDSTAT L2 u ', uL2(1)
        write(*,'(A,F18.12)')' FIELDSTAT L2 v ', uL2(2)
        write(*,'(A,F18.12)')' FIELDSTAT L2 w ', uL2(3)
        write(*,'(A,F18.12)')' FIELDSTAT Linfinity u ', uLinfty(1)
        write(*,'(A,F18.12)')' FIELDSTAT Linfinity v ', uLinfty(2)
        write(*,'(A,F18.12)')' FIELDSTAT Linfinity w ', uLinfty(3)
    endsubroutine

    SUBROUTINE write_continue_(this,fID)
        IMPLICIT NONE
        class(LBMBlock), intent(inout) :: this
        integer,intent(in):: fID
        write(fID) this%xmin,this%ymin,this%zmin,this%dh
        write(fID) this%xDim,this%yDim,this%zDim
        write(fID) this%fIn
    ENDSUBROUTINE write_continue_

    SUBROUTINE read_continue_(this,fID)
        IMPLICIT NONE
        class(LBMBlock), intent(inout) :: this
        integer,intent(in) :: fID
        read(fID) this%xmin,this%ymin,this%zmin,this%dh
        read(fID) this%xDim,this%yDim,this%zDim
        read(fID) this%fIn
    ENDSUBROUTINE read_continue_

    SUBROUTINE evaluate_velocity(time,zCoord,yCoord,xCoord,velocityIn,velocityOut,shearRate)
        implicit none
        real(8):: time,xCoord,yCoord,zCoord
        real(8):: velocityIn(1:SpaceDim),velocityOut(1:SpaceDim),shearRate(1:3)
        if (flow%velocityKind .eq. 0) then
            call evaluate_shear_velocity(zCoord,yCoord,xCoord,velocityIn,velocityOut,shearRate)
        elseif (flow%velocityKind .eq. 2) then
            call evaluate_oscillatory_velocity(time,velocityIn,velocityOut,shearRate)
        endif
    END SUBROUTINE

    ! set the shear velocity
    SUBROUTINE evaluate_shear_velocity(zCoord,yCoord,xCoord,velocityIn,velocityOut,shearRate)
        implicit none
        real(8):: zCoord,yCoord,xCoord
        real(8):: velocityIn(1:SpaceDim),velocityOut(1:SpaceDim),shearRate(1:3)
        velocityOut(1) = velocityIn(1) + 0*shearRate(1) + yCoord*shearRate(2) + zCoord*shearRate(3);
        velocityOut(2) = velocityIn(2) + xCoord*shearRate(1) + 0*shearRate(2) + zCoord*shearRate(3);
        velocityOut(3) = velocityIn(3) + xCoord*shearRate(1) + yCoord*shearRate(2) + 0*shearRate(3);
    END SUBROUTINE

    ! set the oscillatory velocity
    SUBROUTINE evaluate_oscillatory_velocity(time,velocityIn,velocityOut,shearRate)
       implicit none
       real(8):: time
       real(8):: velocityIn(1:SpaceDim),velocityOut(1:SpaceDim),shearRate(1:3)
       real(8):: velocityAmp,velocityFreq,velocityPhi
       velocityAmp = shearRate(1)
       velocityFreq = shearRate(2)
       velocityPhi = shearRate(3)
       velocityOut(1) = velocityIn(1) + velocityAmp * dcos(2*pi*velocityFreq*time + velocityPhi/180.0d0*pi)
       velocityOut(2) = velocityIn(2)
       velocityOut(3) = velocityIn(3)
    END SUBROUTINE

    ! calcualte distribution Function
    SUBROUTINE calculate_distribution_funcion(density,velocity,distribution)
        implicit none
        real(8):: distribution(0:lbmDim)
        real(8):: uSqr,uxyz(0:lbmDim),density,velocity(1:SpaceDim)
        uSqr           = sum(velocity(1:3)**2)
        uxyz(0:lbmDim) = velocity(1) * ee(0:lbmDim,1) + velocity(2) * ee(0:lbmDim,2) + velocity(3) * ee(0:lbmDim,3)
        distribution(0:lbmDim)  = wt(0:lbmDim) * density * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
    END SUBROUTINE

    ! set moving wall boundary distribution Function
    SUBROUTINE evaluate_moving_wall(density,velocity,distributionIn,distributionOut)
        implicit none
        real(8):: distributionIn(0:lbmDim),distributionOut(0:lbmDim)
        real(8):: uxyz(0:lbmDim),density,velocity(1:SpaceDim)
        uxyz(0:lbmDim) = velocity(1) * ee(0:lbmDim,1) + velocity(2) * ee(0:lbmDim,2) + velocity(3) * ee(0:lbmDim,3)
        distributionOut(0:lbmDim) = distributionIn(oppo(0:lbmDim)) + 2.0*wt(0:lbmDim)*density*uxyz(0:lbmDim)*3.0
    ENDSUBROUTINE

    function CompareBlocks(i, j)
        ! if block i contains j, return 1
        ! else if block i in block j, return -1
        ! separated, return 0
        ! partial overlaped, output error
        implicit none
        integer:: i, j, CompareBlocks, cnt
        logical:: d(6)
        cnt = 0
        d(1) = LBMblks(i)%xmin.lt.LBMblks(j)%xmin.or.abs(LBMblks(i)%xmin-LBMblks(j)%xmin).lt.MachineTolerace
        d(2) = LBMblks(j)%xmax.lt.LBMblks(i)%xmax.or.abs(LBMblks(j)%xmax-LBMblks(i)%xmax).lt.MachineTolerace
        d(3) = LBMblks(j)%xmin.lt.LBMblks(i)%xmin.or.abs(LBMblks(j)%xmin-LBMblks(i)%xmin).lt.MachineTolerace
        d(4) = LBMblks(i)%xmax.lt.LBMblks(j)%xmax.or.abs(LBMblks(i)%xmax-LBMblks(j)%xmax).lt.MachineTolerace
        d(5) = LBMblks(i)%xmax.lt.LBMblks(j)%xmin
        d(6) = LBMblks(j)%xmax.lt.LBMblks(i)%xmin
        if(d(1) .and. d(2)) then
            cnt = cnt + 1
        else if(d(3) .and. d(4)) then
            cnt = cnt - 1
        else if(d(5) .and. d(6)) then
            CompareBlocks = 0
            return
        endif
        d(1) = LBMblks(i)%ymin.lt.LBMblks(j)%ymin.or.abs(LBMblks(i)%ymin-LBMblks(j)%ymin).lt.MachineTolerace
        d(2) = LBMblks(j)%ymax.lt.LBMblks(i)%ymax.or.abs(LBMblks(j)%ymax-LBMblks(i)%ymax).lt.MachineTolerace
        d(3) = LBMblks(j)%ymin.lt.LBMblks(i)%ymin.or.abs(LBMblks(j)%ymin-LBMblks(i)%ymin).lt.MachineTolerace
        d(4) = LBMblks(i)%ymax.lt.LBMblks(j)%ymax.or.abs(LBMblks(i)%ymax-LBMblks(j)%ymax).lt.MachineTolerace
        d(5) = LBMblks(i)%ymax.lt.LBMblks(j)%ymin
        d(6) = LBMblks(j)%ymax.lt.LBMblks(i)%ymin
        if(d(1) .and. d(2)) then
            cnt = cnt + 1
        else if(d(3) .and. d(4)) then
            cnt = cnt - 1
        else if(d(5) .and. d(6)) then
            CompareBlocks = 0
            return
        endif
        d(1) = LBMblks(i)%zmin.lt.LBMblks(j)%zmin.or.abs(LBMblks(i)%zmin-LBMblks(j)%zmin).lt.MachineTolerace
        d(2) = LBMblks(j)%zmax.lt.LBMblks(i)%zmax.or.abs(LBMblks(j)%zmax-LBMblks(i)%zmax).lt.MachineTolerace
        d(3) = LBMblks(j)%zmin.lt.LBMblks(i)%zmin.or.abs(LBMblks(j)%zmin-LBMblks(i)%zmin).lt.MachineTolerace
        d(4) = LBMblks(i)%zmax.lt.LBMblks(j)%zmax.or.abs(LBMblks(i)%zmax-LBMblks(j)%zmax).lt.MachineTolerace
        d(5) = LBMblks(i)%zmax.lt.LBMblks(j)%zmin
        d(6) = LBMblks(j)%zmax.lt.LBMblks(i)%zmin
        if(d(1) .and. d(2)) then
            cnt = cnt + 1
        else if(d(3) .and. d(4)) then
            cnt = cnt - 1
        else if(d(5) .and. d(6)) then
            CompareBlocks = 0
            return
        endif
        if(cnt.eq.3) then
            CompareBlocks = 1
        else if(cnt.eq.-3) then
            CompareBlocks = -1
        else
            write(*,*) 'Warning, blocks partial overlaps', LBMblks(i)%ID, LBMblks(j)%ID
            write(*,*) LBMblks(i)%xmin, LBMblks(i)%xmax,LBMblks(i)%ymin, LBMblks(i)%ymax,LBMblks(i)%zmin, LBMblks(i)%zmax
            write(*,*) LBMblks(j)%xmin, LBMblks(j)%xmax,LBMblks(j)%ymin, LBMblks(j)%ymax,LBMblks(j)%zmin, LBMblks(j)%zmax
            !stop
        endif
    end function CompareBlocks

    subroutine find_carrier_fluidblock(x, n)
        implicit none
        real(8)::x(1:SpaceDim)
        integer,intent(out):: n
        integer::i
        real(8)::dh
        dh = 1.d10
        n = -1
        do i=1,m_nblocks
            if( LBMblks(i)%xmin.le.x(1) .and. x(1).le.LBMblks(i)%xmax .and. &
                LBMblks(i)%ymin.le.x(2) .and. x(2).le.LBMblks(i)%ymax .and. &
                LBMblks(i)%zmin.le.x(3) .and. x(3).le.LBMblks(i)%zmax) then
                if(LBMblks(i)%dh .lt. dh) then
                    dh = LBMblks(i)%dh
                    n = i
                endif
            endif
        enddo
        if(n.eq.-1) then
            write(*,*) 'Error: carrier fluid block not found', x
            stop
        endif
    end subroutine

    subroutine FindCarrierFluidBlock()
        use SolidBody, only: m_nFish,VBodies
        implicit none
        integer:: iFish, iBlock
        do iFish=1,m_nFish
            call find_carrier_fluidblock(VBodies(iFish)%v_Exyz(1:3,1), VBodies(iFish)%v_carrierFluidId)
        enddo
        do iBlock=1,m_nblocks
            if(.not.allocated(LBMblks(iBlock)%carriedBodies)) then
                allocate(LBMblks(iBlock)%carriedBodies(0:m_nFish))
            endif
            LBMblks(iBlock)%carriedBodies = 0
            do iFish = 1,m_nFish
                if (iBlock .eq. VBodies(iFish)%v_carrierFluidId) then
                    LBMblks(iBlock)%carriedBodies(0) = LBMblks(iBlock)%carriedBodies(0) + 1
                    LBMblks(iBlock)%carriedBodies(LBMblks(iBlock)%carriedBodies(0)) = iFish
                endif
            enddo
        enddo
    end subroutine
end module FluidDomain
