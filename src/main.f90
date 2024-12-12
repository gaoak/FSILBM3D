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
    character(LEN=40):: filename
    integer:: isubstep=0,step=0,iFish,x,y,z
    real(8):: dt_solid, dt_fluid,get_cpu_time
    real(8):: time=0.0d0,g(3)=[0,0,0]
    integer,dimension(8) :: time_begine1,time_begine2,time_end1,time_end2
    !==================================================================================================
    ! Read all parameters from input file
    call date_and_time(VALUES=time_begine1) ! begine time for the preparation before computing
    filename = 'inflow.dat'
    call read_flow_conditions(filename)
    call read_solid_files(filename)
    call read_fuild_blocks(filename)
    !==================================================================================================
    ! Allocate the memory for simulation
    call allocate_solid_memory_blocks(flow%Asfac,flow%Lchod,flow%Lspan,flow%AR)
    call allocate_fuild_memory_blocks()
    !==================================================================================================
    ! Calculate all the reference values
    call calculate_reference_params()
    call set_solidbody_parameters(flow%denIn,flow%uvwIn,LBMblks(1)%BndConds,&
        flow%Aref,flow%Eref,flow%Fref,flow%Lref,flow%Pref,flow%Tref,flow%Uref)
    !==================================================================================================
    ! Initialization before simulation
    call initialise_solid_bodies(0.d0, g)
    call initialise_fuild_blocks(flow)
    !==================================================================================================
    ! Determine whether to continue calculating and write output informantion titles
    call check_is_continue('./DatTemp/conwr.dat',step,time,flow%isConCmpt)
    call write_information_titles()
    !==================================================================================================
    ! Update the volumn forces and calculate the macro quantities
    call update_volumn_force_blocks(time)
    call calculate_macro_quantities_blocks()
    !==================================================================================================
    ! Write the initial fluid and solid data
    call write_flow_blocks(time)
    call write_solid_field(time)
    call Write_solid_v_bodies(time)
    call date_and_time(VALUES=time_end1)              ! end time for the preparation before computing
    write(*,*)'Time for preparation before computing:', (time_end1 - time_begine1)
    write(*,'(A)') '========================================================'
    !==================================================================================================
    dt_fluid = flow%dt                       !time step of the fluid 
    dt_solid = flow%dt/flow%numsubstep       !time step of the solid
    write(*,*) 'Time loop beginning'
    do while(time/flow%Tref < flow%timeSimTotal)
        call date_and_time(VALUES=time_begine1)
        time = time + dt_fluid
        step = step + 1
        write(*,'(A)') '========================================================'
        write(*,'(A,I6,A,F15.10)')' step:',step,'    time/Tref:',time/flow%Tref
        write(*,'(A)')' ----------------------fluid solver----------------------'
        ! LBM solver
        call date_and_time(VALUES=time_begine2)
        CALL streaming_step()
        call date_and_time(VALUES=time_end2)
        write(*,*)'Time for streaming step:', (time_end2 - time_begine2)
        ! Set boundary conditions
        CALL setBndConds()
        CALL calculate_macro_quantities_blocks()
        call date_and_time(VALUES=time_begine2)

        !compute force exerted on fluids
        CALL FSInteraction_force(dt_fluid,LBMblks(1)%dh,LBMblks(1)%xmin,LBMblks(1)%ymin,LBMblks(1)%zmin,LBMblks(1)%xDim,LBMblks(1)%yDim,LBMblks(1)%zDim,LBMblks(1)%uuu,LBMblks(1)%force)
        !compute volume force exerted on fluids
        call date_and_time(VALUES=time_end2)
        write(*,*)'time for IBM step : ',get_cpu_time(time_end2)-get_cpu_time(time_begine2)

        call date_and_time(VALUES=time_begine2)
        CALL collision_step()
        call date_and_time(VALUES=time_end2)
        write(*,*)'time for collision_step:',get_cpu_time(time_end2)-get_cpu_time(time_begine2)
        !******************************************************************************************
        !solve solid
        write(*,'(A)')' ----------------------solid solver----------------------'
        call date_and_time(VALUES=time_begine2)
        do isubstep=1,flow%numsubstep
            call Solver(time,isubstep,dt_fluid,dt_solid)
        enddo !do isubstep=1,numsubstep
        write(*,'(A)')' --------------------------------------------------------'
        call date_and_time(VALUES=time_end1)
        write(*,*)'time for one step:',get_cpu_time(time_end1)-get_cpu_time(time_begine1)
    enddo
    call ComputeFieldStat()
END PROGRAM main
