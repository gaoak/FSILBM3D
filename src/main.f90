!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    copyright@ RuNanHua
!    The main program, 3D Lattice Boltzmann Method
!    flow past a flexible plate (uniform flow past a flag or flapping plate, and so on)
!    flexible Plates or shells move in fluid.
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    PROGRAM main
    USE ConstParams
    USE FluidDomain
    USE FlowCondition
    USE SolidBody
    implicit none
    integer:: isubstep,iFish,x,y,z
    real(8):: Pbetatemp,CPUtime
    logical alive
    integer,dimension(8) :: values0,values1,values_s,values_e


    call read_flow_conditions()
    call initialise_fuild_blocks()
    call read_body_conditions()
    call initialise_solid_bodies()
    call calculate_reference_params()
    call write_params()    


    !call Initialise_Calculate_Solid_params()
    !CALL allocate_fluid_memory()
    !CALL calculate_LB_params()
    !CALL calculate_MRTM_params()
    !call OMP_set_num_threads(npsize)
    !write(*,*)'npsize=', npsize

    if(isRelease==1)then
        write(*,*)'Release'
    endif
    Pbetatemp=Pbeta
    deltat = dt  !set time step of solid deltat the same as fluid time step
!===============================================================================================
    time=0.0d0
    step=0
    CALL Initialise_bodies(npsize,time,zDim,yDim,xDim,dh,g)
    CALL initialize_flow()
    if(ismovegrid==1)then
        iFish = 1
        call cptIref(NDref,IXref,IYref,IZref,VBodies(iFish)%rbm%nND,xDim,yDim,zDim,VBodies(iFish)%rbm%xyzful(1:VBodies(iFish)%rbm%nND,1:3),xGrid,yGrid,zGrid,Xref,Yref,Zref)
        MoveOutputIref=0
    endif
    if(step==0)    CALL wrtInfoTitl()

    inquire(file='./DatTemp/conwr.dat', exist=alive)
    if (isConCmpt==1 .and. alive)then
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'Continue compute!!!!!!!!!!!!!!!!!!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        call read_checkpoint_file()
    else
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'NEW      compute!!!!!!!!!!!!!!!!!!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!===============================================================================================
    endif
!===============================================================================================
    CALL updateVolumForc()
    CALL calculate_macro_quantities()
    CALL write_flow_fast()
    CALL write_solid_field(time)
    call Write_solid_v_bodies(time)
!==================================================================================================
!==================================================================================================
!==================================================================================================
    write(*,*)'time loop'
    do while(time/Tref < timeSimTotl)
        call date_and_time(VALUES=values_s)
        time=time+dt
        step=step+1
        write(*,'(A)')' ======================================================================='
        write(*,'(A,I6,A,F15.10)')' step:',step,'    time/Tref:',time/Tref
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        write(*,'(A)')' ----------------------fluid solver----------------------'
        !solve fluid
        call date_and_time(VALUES=values0)
        CALL streaming_step()
        call date_and_time(VALUES=values1)
        write(*,*)'time for streaming_step:',CPUtime(values1)-CPUtime(values0)
        !******************************************************************************************
        !DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304,fluid=0, wall=200, movingWall=201
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
        write(*,*)'time for IBM step : ',CPUtime(values1)-CPUtime(values0)
        if(time/Tref >begForcDist .and. time/Tref <endForcDist) call forcDisturb() !force disturbance for instability

        call date_and_time(VALUES=values0)
        !call cptForceR(dxmin,dymin,dzmin,nND,nND_max,xyzful,repful,nFish,Lspan)
        call date_and_time(VALUES=values1)
        write(*,*)'time for Lubforce :',CPUtime(values1)-CPUtime(values0)

        call date_and_time(VALUES=values0)
        CALL collision_step()
        call date_and_time(VALUES=values1)
        write(*,*)'time for collision_step:',CPUtime(values1)-CPUtime(values0)
        !******************************************************************************************
        !solve solid
        if(isRelease/=1)write(*,'(A)')' ----------------------solid solver----------------------'
        call date_and_time(VALUES=values0)
        subdeltat=deltat/numsubstep
        do isubstep=1,numsubstep
            call Solver(time,isubstep,deltat,subdeltat)
        enddo !do isubstep=1,numsubstep
        call date_and_time(VALUES=values1)
        write(*,*)'time for FEM:',CPUtime(values1)-CPUtime(values0)
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
        write(*,*)'max time for FEM  : ',CPUtime(values1)-CPUtime(values0)
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
        write(*,*)'time for one step:',CPUtime(values_e)-CPUtime(values_s)
    enddo
    call ComputeFieldStat
    END PROGRAM main
