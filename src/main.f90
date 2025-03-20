!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    copyright@ An-Kang Gao
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
    character(LEN=40):: parameterFile='inFlow.dat',checkFile='check.dat'
    integer:: step=0,start_ave
    real(8):: dt_fluid
    real(8):: time=0.0d0,start_time=0.0d0,g(3)=[0,0,0]
    real(8):: time_collision,time_streaming,time_IBM,time_FEM,time_begine1,time_begine2,time_end1,time_end2
    !==================================================================================================
    ! Read all parameters from input file
    call get_now_time(time_begine1) ! begine time for the preparation before computing
    call read_flow_conditions(parameterFile)
    call read_solid_files(parameterFile)
    call read_fuild_blocks(parameterFile)
    call read_probe_params(parameterFile)
    call bluid_block_tree()
    !==================================================================================================
    ! Set parallel compute cores
    call omp_set_num_threads(flow%npsize)
    !==================================================================================================
    ! Allocate the memory for simulation
    call allocate_solid_memory(flow%Asfac,flow%Lchod,flow%Lspan,flow%AR)
    call allocate_fuild_memory_blocks(flow%npsize)
    !==================================================================================================
    ! Calculate all the reference values
    call calculate_reference_params(flow)
    call set_solidbody_parameters(flow%denIn,flow%uvwIn,LBMblks(blockTreeRoot)%BndConds,&
        flow%Aref,flow%Eref,flow%Fref,flow%Lref,flow%Pref,flow%Tref,flow%Uref,flow%ntolLBM,flow%dtolLBM)
    !==================================================================================================
    ! Initialization before simulation
    call initialise_solid_bodies(0.d0, g)
    call FindCarrierFluidBlock()
    call initialise_fuild_blocks(time)
    !==================================================================================================
    ! Check blocks number and check the tau of each block
    call check_blocks_params(blockTreeRoot)
    call write_parameter_check_file(checkFile)
    call Write_solid_Check(checkFile)
    !==================================================================================================
    ! Determine whether to continue calculating and write output informantion titles
    call check_is_continue(step,start_time,flow%isConCmpt) ! read dimensionless time in continue file as the computing start time
    call write_information_titles(m_nFish)
    !==================================================================================================
    ! Update the volume forces and calculate the macro quantities
    call update_volume_force_blocks()
    call tree_set_boundary_conditions_block(blockTreeRoot)
    call calculate_macro_quantities_blocks()
    !==================================================================================================
    ! Determine the time and step of continue computing
    dt_fluid = LBMblks(blockTreeRoot)%dh     !time step of the fluid 
    if(flow%isConCmpt .eq. 1) then
        time = start_time * flow%Tref        !calculate the computing start time
    else
        start_time = 0.0d0
        step = 0
    endif
    if(flow%timeWriteBegin .ge. start_time) then
        start_ave = step + nint((flow%timeWriteBegin - start_time) * flow%Tref / dt_fluid)  !the first step for fluid averaging
    else
        start_ave = step
    endif
    write(*,*)'the start step for fluid averaging(if used):', start_ave
    write(*,'(A)') '========================================================='
    !==================================================================================================
    ! Write the initial fluid and solid data
    call write_flow_blocks(time)
    call write_solid_field(time)
    call Write_solid_v_bodies(time)
    call get_now_time(time_end1)             !end time for the preparation before computing
    write(*,*)'Time for preparation before computing:', (time_end1 - time_begine1)
    write(*,'(A)') '========================================================='
    !==================================================================================================
    write(*,*) 'Time loop beginning'
    do while(time/flow%Tref < flow%timeSimTotal)
        call get_now_time(time_begine1)
        time = time + dt_fluid
        step = step + 1
        LBMblks(:)%blktime = time
        write(*,'(A)') '========================================================='
        write(*,'(A,I8,A,F14.8)')' Steps:',step,'  Time/Tref:',time/flow%Tref
        write(*,'(A)')' --------------------- fluid solver ---------------------'
        ! LBM solver
        time_collision = 0.d0
        time_streaming = 0.d0
        time_IBM       = 0.d0
        time_FEM       = 0.d0
        call tree_collision_streaming_IBM_FEM(blockTreeRoot,time_collision,time_streaming,time_IBM,time_FEM)
        call calculate_macro_quantities_blocks()
        call calculate_turbulent_statistic_blocks(step,start_ave)
        write(*,*)'Time for collision step:', time_collision
        write(*,*)'Time for streaming step:', time_streaming
        write(*,*)'Time for       IBM step:', time_IBM
        write(*,'(A)')' --------------------- solid solver ---------------------'
        write(*,*)'Time for     solid step:', time_FEM
        write(*,'(A)')' ---------------------- write info ----------------------'
        call get_now_time(time_begine2)
        ! write data for continue computing
        if(DABS(time/flow%Tref-flow%timeContiDelta*NINT(time/flow%Tref/flow%timeContiDelta)) <= 0.5*dt_fluid/flow%Tref)then
            call write_continue_blocks(step,time / flow%Tref)   ! output dimensionless time
        endif
        ! write fluid and soild data
        if((time/flow%Tref - flow%timeWriteBegin) >= -0.5*dt_fluid/flow%Tref .and. (time/flow%Tref - flow%timeWriteEnd) <= 0.5*dt_fluid/flow%Tref) then
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
                                              LBMblks(flow%inWhichBlock)%xDim,LBMblks(flow%inWhichBlock)%yDim,LBMblks(flow%inWhichBlock)%zDim,LBMblks(flow%inWhichBlock)%uuu,LBMblks(flow%inWhichBlock)%den)
            call write_solid_Information(time,flow%timeInfoDelta,flow%Asfac,flow%solidProbingNum,flow%solidProbingNode)
        endif
        call get_now_time(time_end2)
        write(*,*)'Time for   writing step:', (time_end2 - time_begine2)
        write(*,'(A)')' ----------------------- one step -----------------------'
        call get_now_time(time_end1)
        write(*,*)'Time for       one step:', (time_end1 - time_begine1)
    enddo
    ! write validation informations
    write(*,'(A)') '========================================================='
    call computeFieldStat_blocks()
END PROGRAM main
