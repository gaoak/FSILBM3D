!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    copyright@ RuNanHua
!    The main program, 3D Lattice Boltzmann Method
!    flow past a flexible plate (uniform flow past a flag or flapping plate, and so on)
!    flexible Plates or shells move in fluid.
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main
    USE omp_lib
    USE ConstParams
    USE FlowCondition
    USE SolidBody
    USE FluidDomain
    USE LBMBlockComm
    implicit none
    character(LEN=40):: parameterFile='inFlow.dat',continueFile='continue.dat',checkFile='check.dat'
    integer:: isubstep=0,step=0
    integer:: n_pairs,n_gridDelta
    real(8):: dt_solid, dt_fluid
    real(8):: time=0.0d0,g(3)=[0,0,0]
    real(8):: time_begine1,time_begine2,time_end1,time_end2
    !==================================================================================================
    ! Read all parameters from input file
    call get_now_time(time_begine1) ! begine time for the preparation before computing
    call read_flow_conditions(parameterFile)
    call read_solid_files(parameterFile)
    call read_fuild_blocks(parameterFile)
    call read_probe_params(parameterFile)
    call read_blocks_comunication(parameterFile)
    !==================================================================================================
    ! Set parallel compute cores
    call omp_set_num_threads(flow%npsize)
    !==================================================================================================
    ! Allocate the memory for simulation
    call allocate_solid_memory(flow%Asfac,flow%Lchod,flow%Lspan,flow%AR)
    call allocate_fuild_memory_blocks(flow%npsize)
    !==================================================================================================
    ! Calculate the tau of each block
    call calculate_blocks_tau()
    !==================================================================================================
    ! Calculate all the reference values
    call calculate_reference_params(flow)
    call set_solidbody_parameters(flow%denIn,flow%uvwIn,LBMblks(1)%BndConds,&
        flow%Aref,flow%Eref,flow%Fref,flow%Lref,flow%Pref,flow%Tref,flow%Uref,flow%ntolLBM,flow%dtolLBM)
    call write_parameter_check_file(checkFile)
    !==================================================================================================
    ! Initialization before simulation
    call initialise_solid_bodies(0.d0, g)
    call initialise_fuild_blocks(flow)
    !==================================================================================================
    ! Determine whether to continue calculating and write output informantion titles
    call check_is_continue(continueFile,step,time,flow%isConCmpt)
    call write_information_titles(m_nFish)
    !==================================================================================================
    ! Update the volumn forces and calculate the macro quantities
    call update_volumn_force_blocks(time)
    call set_boundary_conditions_block()
    call calculate_macro_quantities_blocks()
    !==================================================================================================
    ! Write the initial fluid and solid data
    call write_flow_blocks(time)
    call write_solid_field(time)
    call Write_solid_v_bodies(time)
    call get_now_time(time_end1)             ! end time for the preparation before computing
    write(*,*)'Time for preparation before computing:', (time_end1 - time_begine1)
    write(*,'(A)') '========================================================='
    !==================================================================================================
    dt_fluid = flow%dt                       !time step of the fluid 
    dt_solid = flow%dt/flow%numsubstep       !time step of the solid
    write(*,*) 'Time loop beginning'
    do while(time/flow%Tref < flow%timeSimTotal)
        call get_now_time(time_begine1)
        time = time + dt_fluid
        step = step + 1
        write(*,'(A)') '========================================================='
        write(*,'(A,I6,A,F14.8)')' Steps:',step,'  Time/Tref:',time/flow%Tref
        write(*,'(A)')' --------------------- fluid solver ---------------------'
        ! LBM solver
        if(m_npairs .eq. 0) then ! single block
            call get_now_time(time_begine2)
            call collision_block(1)
            call get_now_time(time_end2)
            write(*,*)'Time for collision step:', (time_end2 - time_begine2)
            call get_now_time(time_begine2)
            call streaming_block(1)
            call get_now_time(time_end2)
            write(*,*)'Time for streaming step:', (time_end2 - time_begine2) 
        elseif(m_npairs .eq. 1) then ! two blocks
            call calculating_public_distribution()
            call collision_block(commpairs(1)%fatherId)
            call streaming_block(commpairs(1)%fatherId)
            do n_gridDelta=1,m_gridDelta
                call blocks_interpolation(1)
                call collision_block(commpairs(1)%sonId)
                call streaming_block(commpairs(1)%sonId)
            enddo
            if(m_npairs .ge. 2) then ! multi-blocks
                do n_pairs=2,m_npairs
                    do n_gridDelta=1,m_gridDelta
                        call blocks_interpolation(n_pairs)
                        call collision_block(commpairs(n_pairs)%sonId)
                        call streaming_block(commpairs(n_pairs)%sonId)
                    enddo
                enddo
            endif
            call set_boundary_conditions_block(commpairs(1)%fatherId)
        endif
        call set_boundary_conditions_block(1)
        call calculate_macro_quantities_blocks()
        write(*,'(A)')' --------------------- solid solver ---------------------'
        call clear_volume_force()
        call get_now_time(time_begine2)
        call FSInteraction_force(dt_fluid,LBMblks(m_containSolidId)%dh,LBMblks(m_containSolidId)%xmin,LBMblks(m_containSolidId)%ymin,LBMblks(m_containSolidId)%zmin, &
                                LBMblks(m_containSolidId)%xDim,LBMblks(m_containSolidId)%yDim,LBMblks(m_containSolidId)%zDim,LBMblks(m_containSolidId)%uuu,LBMblks(m_containSolidId)%force)
        call get_now_time(time_end2)
        write(*,*)'Time   for   IBM   step:', (time_end2 - time_begine2)
        !call get_now_time(time_begine2)
        !call streaming_block()
        !call get_now_time(time_end2)
        !write(*,*)'Time for streaming step:', (time_end2 - time_begine2)
        ! Set fluid boundary conditions
        !call set_boundary_conditions_block()
        !if(m_npairs .ne. 0) call ExchangeFluidInterface()
        !call calculate_macro_quantities_blocks()
        ! Compute volume force exerted on fluids
        !call get_now_time(time_begine2)
        !call clear_volume_force()
        !call FSInteraction_force(dt_fluid,LBMblks(m_containSolidId)%dh,LBMblks(m_containSolidId)%xmin,LBMblks(m_containSolidId)%ymin,LBMblks(m_containSolidId)%zmin, &
        !                         LBMblks(m_containSolidId)%xDim,LBMblks(m_containSolidId)%yDim,LBMblks(m_containSolidId)%zDim,LBMblks(m_containSolidId)%uuu,LBMblks(m_containSolidId)%force)
        !call get_now_time(time_end2)
        !write(*,*)'Time   for   IBM   step:', (time_end2 - time_begine2)
        ! Fluid collision
        !call get_now_time(time_begine2)
        !call collision_block()
        !call get_now_time(time_end2)
        !write(*,*)'Time for collision step:', (time_end2 - time_begine2)
        !IBM solver
        !write(*,'(A)')' --------------------- solid solver ---------------------'
        call get_now_time(time_begine2)
        do isubstep=1,flow%numsubstep
            call Solver(time,isubstep,dt_fluid,dt_solid)
        enddo !do isubstep=1,numsubstep
        call get_now_time(time_end2)
        write(*,*)'Time   for  solid  step:', (time_end2 - time_begine2)
        write(*,'(A)')' ---------------------- write infos ---------------------'
        call get_now_time(time_begine2)
        ! write data for continue computing
        if(DABS(time/flow%Tref-flow%timeContiDelta*NINT(time/flow%Tref/flow%timeContiDelta)) <= 0.5*dt_fluid/flow%Tref)then
            call write_continue_blocks(continueFile,step,time)
        endif
        ! write fluid and soild data
        if((flow%timeWriteBegin .le. time/flow%Tref) .and. (time/flow%Tref .le. flow%timeWriteEnd)) then
            if(DABS(time/flow%Tref-flow%timeBodyDelta*NINT(time/flow%Tref/flow%timeBodyDelta)) <= 0.5*dt_fluid/flow%Tref)then
                call write_solid_field(time)
                call Write_solid_v_bodies(time)
            endif
            if(DABS(time/flow%Tref-flow%timeFlowDelta*NINT(time/flow%Tref/flow%timeFlowDelta)) <= 0.5*dt_fluid/flow%Tref)then
                call write_flow_blocks(time)
            endif
        endif
        ! write processing informations
        if(DABS(time/flow%Tref-flow%timeInfoDelta*NINT(time/flow%Tref/flow%timeInfoDelta)) <= 0.5*dt_fluid/flow%Tref)then
            call write_fluid_information(time,LBMblks(flow%inWhichBlock)%dh,LBMblks(flow%inWhichBlock)%xmin,LBMblks(flow%inWhichBlock)%ymin,LBMblks(flow%inWhichBlock)%zmin, &
                                              LBMblks(flow%inWhichBlock)%xDim,LBMblks(flow%inWhichBlock)%yDim,LBMblks(flow%inWhichBlock)%zDim,LBMblks(flow%inWhichBlock)%uuu)
            call write_solid_information(time,m_nFish)
        endif
        call get_now_time(time_end2)
        write(*,*)'Time  for writing  step:', (time_end2 - time_begine2)
        write(*,'(A)')' ----------------------- one step -----------------------'
        call get_now_time(time_end1)
        write(*,*)'Time   for   one   step:', (time_end1 - time_begine1)
    enddo
    ! write validation informations
    write(*,'(A)') '========================================================='
    call computeFieldStat_blocks()
END PROGRAM main
