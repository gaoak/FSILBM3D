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
logical alive
integer:: isubstep=0,step=0,iFish,x,y,z
real(8):: time=0.0d0,g(3)=[0,0,0],time_begine1,time_begine2,time_end1,time_end2
!==================================================================================================
! Read all parameters from input file
call get_cpu_time(time_begine1)          ! begine time for the preparation before computing
call read_flow_conditions('inflow.dat')
call read_solid_files('inflow.dat')
call read_fuild_blocks('inflow.dat')
!==================================================================================================
! Allocate the memory for simulation
call allocate_solid_memory_blocks(flow%Asfac,flow%Lchod,flow%Lspan,flow%AR)
call allocate_fuild_memory()
!==================================================================================================
! Calculate all the reference values
call calculate_reference_params()
call set_solidbody_parameters(flow%dt/dble(flow%numsubstep),LBMblks(1)%dh,flow%denIn,&
    flow%uvwIn,LBMblks(1)%BndConds,LBMblks(1)%zDim,LBMblks(1)%yDim,LBMblks(1)%xDim,&
    flow%Aref,flow%Eref,flow%Fref,flow%Lref,flow%Pref,flow%Tref,flow%Uref)
!==================================================================================================
! Initialization before simulation
call initialise_solid_bodies(0, g)
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
call get_cpu_time(time_end1)              ! end time for the preparation before computing
write(*,*)'Time for preparation before computing:', (time_end1 - time_begine1)
write(*,'(A)') '========================================================'
!==================================================================================================
dt_fluid = flow%dt                       !time step of the fluid 
dt_solid = flow%dt/flow%numsubstep       !time step of the solid
write(*,*) 'Time loop beginning'
do while(time/flow%Tref < flow%timeSimTotal)
    call get_cpu_time(time_begine1)
    time = time + dt_fluid
    step = step + 1
    write(*,'(A)') '========================================================'
    write(*,'(A,I6,A,F15.10)')' step:',step,'    time/Tref:',time/Tref
    write(*,'(A)')' ----------------------fluid solver----------------------'
    ! LBM solver
    call get_cpu_time(time_begine2)
    CALL streaming_step()
    call get_cpu_time(time_end2)
    write(*,*)'Time for streaming step:', (time_end2 - time_begine2)
    ! Set boundary conditions
    if    (iStreamModel==1)then
    xMinBC=boundaryConditions(1)
    xMaxBC=boundaryConditions(2)
    yMinBC=boundaryConditions(3)
    yMaxBC=boundaryConditions(4)
    zMinBC=boundaryConditions(5)
    zMaxBC=boundaryConditions(6)
    CALL set_other_farfld_BCs()
    elseif(iStreamModel==2)then
    CALL set_equilibrium_farfld_BC()
    else
        stop
        write(*,*)'no such type LBMBC'
    endif
    !******************************************************************************************
    if(ismovegrid==1)then ! move fluid grid
        iFish=1
        call cptMove(move(1:3),VBodies(iFish)%rbm%xyzful(NDref,1:3),[xGrid(IXref),yGrid(IYref),zGrid(IZref)],dh,MoveOutputIref(1:3))
        write(*,'(A,3F10.5)')' *Grid:',[xGrid(IXref),yGrid(IYref),zGrid(IZref)]
        write(*,'(A,3F10.5)')' *Body:',VBodies(iFish)%rbm%xyzful(NDref,1:3)
        if(isMoveDimX==1)CALL movGrid(1,move(1))
        if(isMoveDimY==1)CALL movGrid(2,move(2))
        if(isMoveDimZ==1)CALL movGrid(3,move(3))
        !******************************************************************************************
        !DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304,fluid=0, wall=200
        if    (iStreamModel==1)then
        xMinBC=boundaryConditions(1)
        xMaxBC=boundaryConditions(2)
        yMinBC=boundaryConditions(3)
        yMaxBC=boundaryConditions(4)
        zMinBC=boundaryConditions(5)
        zMaxBC=boundaryConditions(6)
        CALL set_other_farfld_BCs()
        elseif(iStreamModel==2)then
        CALL set_equilibrium_farfld_BC()
        else
            stop
            write(*,*)'no such type LBMBC'
        endif
    endif
    !******************************************************************************************
    CALL updateVolumForc()
    CALL calculate_macro_quantities()
    !******************************************************************************************
    !******************************************************************************************
    call date_and_time(VALUES=values0)
    !Pbeta=(1.0d0-dexp(-5.0d0/Pramp*time/Tref))*Pbetatemp
    Pbeta=1.0d0*Pbetatemp

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
    do x=1,xDim
        do y=1,yDim
            do z=1,zDim
                force(z,y,x,1:3) = 0.d0
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO
    !compute force exerted on fluids
    CALL FSInteraction_force(xGrid,yGrid,zGrid,uuu,force)
    !compute volume force exerted on fluids
    CALL addVolumForc()

    call date_and_time(VALUES=values1)
    write(*,*)'time for IBM step : ',get_cpu_time(values1)-get_cpu_time(values0)
    if(time/Tref >begForcDist .and. time/Tref <endForcDist) call forcDisturb() !force disturbance for instability

    call date_and_time(VALUES=values0)
    !call cptForceR(dxmin,dymin,dzmin,nND,nND_max,xyzful,repful,nFish,Lspan)
    call date_and_time(VALUES=values1)
    write(*,*)'time for Lubforce :',get_cpu_time(values1)-get_cpu_time(values0)

    call date_and_time(VALUES=values0)
    CALL collision_step()
    call date_and_time(VALUES=values1)
    write(*,*)'time for collision_step:',get_cpu_time(values1)-get_cpu_time(values0)
    !******************************************************************************************
    !solve solid
    if(isRelease/=1)write(*,'(A)')' ----------------------solid solver----------------------'
    call date_and_time(VALUES=values0)
    subdeltat=deltat/numsubstep
    do isubstep=1,numsubstep
        call Solver(time,isubstep,deltat,subdeltat)
    enddo !do isubstep=1,numsubstep
    call date_and_time(VALUES=values1)
    write(*,*)'time for FEM:',get_cpu_time(values1)-get_cpu_time(values0)
    !----------------------------------------------------------------------
    !******************************************************************************************
    !******************************************************************************************
    !******************************************************************************************
    if(isRelease/=1)then
        do iFish=1,nFish
            write(*,'(A,I5.5)')' Fish number: ', int(VBodies(iFish)%rbm%FishInfo(1))
            write(*,'(A,I5.5,A,E20.10)')' iterFEM = ',int(VBodies(iFish)%rbm%FishInfo(2)),'    dmaxFEM = ',VBodies(iFish)%rbm%FishInfo(3)
        enddo
    endif
    write(*,'(A)')' --------------------------------------------------------'
    call date_and_time(VALUES=values1)
    write(*,*)'max time for FEM  : ',get_cpu_time(values1)-get_cpu_time(values0)
    !******************************************************************************************
    !******************************************************************************************
    !******************************************************************************************
    if(isRelease/=1) write(*,'(A)')' ----------------------post process----------------------'
    if(isRelease/=1)then
        do iFish=1,nFish
            write(*,'(A,I5.5)')' Fish number: ',iFish
            write(*,'(A,3D15.5)')" forceDre: ",sum(VBodies(iFish)%rbm%extful(1:VBodies(iFish)%rbm%nND,1:3),1)/Fref
            write(*,'(A,3D15.5)')" accCentM: ",sum(VBodies(iFish)%rbm%accful(1:VBodies(iFish)%rbm%nND,1:3)*VBodies(iFish)%rbm%mssful(1:VBodies(iFish)%rbm%nND,1:3),1)/sum(VBodies(iFish)%rbm%mssful(1:VBodies(iFish)%rbm%nND,1:3),1)/Aref
            write(*,'(A,3D15.5)')" velCentM: ",sum(VBodies(iFish)%rbm%velful(1:VBodies(iFish)%rbm%nND,1:3)*VBodies(iFish)%rbm%mssful(1:VBodies(iFish)%rbm%nND,1:3),1)/sum(VBodies(iFish)%rbm%mssful(1:VBodies(iFish)%rbm%nND,1:3),1)/Uref
            write(*,'(A,3D15.5)')" xyzCentM: ",sum(VBodies(iFish)%rbm%xyzful(1:VBodies(iFish)%rbm%nND,1:3)*VBodies(iFish)%rbm%mssful(1:VBodies(iFish)%rbm%nND,1:3),1)/sum(VBodies(iFish)%rbm%mssful(1:VBodies(iFish)%rbm%nND,1:3),1)/Lref
        enddo
    endif
    !******************************************************************************************
    !******************************************************************************************
    !******************************************************************************************
    if(DABS(time/Tref-timeOutTemp*NINT(time/Tref/timeOutTemp)) <= 0.5*dt/Tref)then
        CALL write_checkpoint_file()
    endif

    if((timeOutBegin .le. time/Tref) .and. (time/Tref .le. timeOutEnd)) then
        if(DABS(time/Tref-timeOutBody*NINT(time/Tref/timeOutBody)) <= 0.5*dt/Tref)then
            CALL write_solid_field(time)
            call Write_solid_v_bodies(time)
        endif
        if(DABS(time/Tref-timeOutFlow*NINT(time/Tref/timeOutFlow)) <= 0.5*dt/Tref)then
            CALL write_flow_fast()
        endif
    endif

    if(DABS(time/Tref-timeOutInfo*NINT(time/Tref/timeOutInfo)) <= 0.5*dt/Tref)then
        CALL wrtInfo()
    endif
    !******************************************************************************************
    !******************************************************************************************
    !******************************************************************************************
    write(*,'(A)')' --------------------------------------------------------'
    call date_and_time(VALUES=values_e)
    write(*,*)'time for one step:',get_cpu_time(values_e)-get_cpu_time(values_s)
enddo
call ComputeFieldStat
END PROGRAM main
