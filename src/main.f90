!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    copyright@ RuNanHua
!    The main program, 3D Lattice Boltzmann Method
!    flow past a flexible plate (uniform flow past a flag or flapping plate, and so on)
!    flexible Plates or shells move in fluid.
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main
    USE ConstParams
    USE FlowCondition
    USE SolidBody
    USE FluidDomain
    implicit none
    character(LEN=40):: parameterFile='inflow.dat',continueFile='continue.dat'
    integer:: isubstep=0,step=0
    real(8):: dt_solid, dt_fluid
    real(8):: time=0.0d0,g(3)=[0,0,0]
    real(8):: time_begine1,time_begine2,time_end1,time_end2
    !==================================================================================================
    ! Read all parameters from input file
    call get_now_time(time_begine1) ! begine time for the preparation before computing
    call read_flow_conditions(parameterFile)
    call read_solid_files(parameterFile)
    call read_fuild_blocks(parameterFile)
    !==================================================================================================
    ! Allocate the memory for simulation
    call allocate_solid_memory(flow%Asfac,flow%Lchod,flow%Lspan,flow%AR)
    call allocate_fuild_memory_blocks()
    !==================================================================================================
    ! Calculate all the reference values
    call calculate_reference_params(flow)
    call set_solidbody_parameters(flow%denIn,flow%uvwIn,LBMblks(1)%BndConds,&
        flow%Aref,flow%Eref,flow%Fref,flow%Lref,flow%Pref,flow%Tref,flow%Uref)
    !==================================================================================================
    ! Initialization before simulation
    call initialise_solid_bodies(0.d0, g)
    call initialise_fuild_blocks(flow)
    !==================================================================================================
    ! Determine whether to continue calculating and write output informantion titles
    call check_is_continue(continueFile,step,time,flow%isConCmpt)
    !call write_information_titles()
    !==================================================================================================
    ! Update the volumn forces and calculate the macro quantities
    call update_volumn_force_blocks(time)
    call calculate_macro_quantities_blocks()
    !==================================================================================================
    ! Write the initial fluid and solid data
    call write_flow_blocks(time)
    call write_solid_field(time)
    call Write_solid_v_bodies(time)
    call get_now_time(time_end1)             ! end time for the preparation before computing
    write(*,*)'Time for preparation before computing:', (time_end1 - time_begine1)
    write(*,'(A)') '========================================================'
    !==================================================================================================
    dt_fluid = flow%dt                       !time step of the fluid 
    dt_solid = flow%dt/flow%numsubstep       !time step of the solid
    write(*,*) 'Time loop beginning'
    do while(time/flow%Tref < flow%timeSimTotal)
        call get_now_time(time_begine1)
        time = time + dt_fluid
        step = step + 1
        write(*,'(A)') '========================================================'
        write(*,'(A,I6,A,F15.10)')' Steps:',step,'    Time/Tref:',time/flow%Tref
        write(*,'(A)')' ----------------------fluid solver----------------------'
        ! LBM solver
        call get_now_time(time_begine2)
        CALL streaming_blocks()
        call get_now_time(time_end2)
        write(*,*)'Time for streaming step:', (time_end2 - time_begine2)
        ! Set fluid boundary conditions
        call set_boundary_conditions_blocks()
        call calculate_macro_quantities_blocks()
        ! Compute volume force exerted on fluids
        call get_now_time(time_begine2)
        CALL FSInteraction_force(dt_fluid,LBMblks(1)%dh,LBMblks(1)%xmin,LBMblks(1)%ymin,LBMblks(1)%zmin,LBMblks(1)%xDim,LBMblks(1)%yDim,LBMblks(1)%zDim,LBMblks(1)%uuu,LBMblks(1)%force)
        call get_now_time(time_end2)
        write(*,*)'Time   for   IBM   step:', (time_end2 - time_begine2)
        ! Fluid collision
        call get_now_time(time_begine2)
        call collision_blocks()
        call get_now_time(time_end2)
        write(*,*)'time for collision_step:', (time_end2 - time_begine2)
        !IBM solver
        write(*,'(A)')' ----------------------solid solver----------------------'
        do isubstep=1,flow%numsubstep
            call Solver(time,isubstep,dt_fluid,dt_solid)
        enddo !do isubstep=1,numsubstep
        write(*,'(A)')' --------------------------------------------------------'
        call get_now_time(time_end1)
        write(*,*)'time   for   one   step:', (time_end1 - time_begine1)
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
        !if(DABS(time/flow%Tref-flow%timeInfoDelta*NINT(time/flow%Tref/flow%timeInfoDelta)) <= 0.5*dt_fluid/flow%Tref)then
        !    call write_information()
        !endif
    enddo
    ! write validation informations
    call computeFieldStat_blocks()
END PROGRAM main
