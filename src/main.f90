!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    copyright@ RuNanHua
!    The main program, 3D Lattice Boltzmann Method
!    flow past a flexible plate (uniform flow past a flag or flapping plate, and so on)
!    flexible Plates or shells move in fluid.
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    PROGRAM main
    USE simParam
    use omp_lib
    USE ImmersedBoundary
    USE FakeBody
    USE FakeBodyspace
    implicit none
    integer:: iND,isubstep,iFish,x,y,z,icount
    real(8), allocatable:: FishInfo(:,:)
    real(8):: Pbetatemp,CPUtime
    logical alive
    !time_and_date
    integer,dimension(8) :: values0,values1,values_s,values_e
    CALL read_file()
    CALL allocate_solid_memory()
    CALL allocate_fluid_memory()
    CALL calculate_LB_params()
    CALL write_params()
    CALL calculate_MRTM_params()
    call Beam_Initial(nFish,nND,xyzful00,XYZ)

    allocate(FishInfo(1:3,1:nFish))

    call OMP_set_num_threads(npsize)
    write(*,*)'npsize=', npsize

    if(isRelease==1)then
        write(*,*)'Release'
    endif
    Pbetatemp=Pbeta
    deltat = dt  !set time step of solid deltat the same as fluid time step
!===============================================================================================
    time=0.0d0
    step=0
    CALL initialize_solid()
    CALL initialize_flow()
    if(ismovegrid==1)then
        iFish = 1
        call cptIref(NDref,IXref,IYref,IZref,nND(iFish),xDim,yDim,zDim,xyzful(1:nND(iFish),1:3,iFish),xGrid,yGrid,zGrid,Xref,Yref,Zref)
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
    CALL write_solid_field(xyzful/Lref,velful/Uref,accful/Aref,extful/Fref,repful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,nFish)
    call Beam_write_solid_fake_field(nFish,Lref,time,Tref)
    if (maxval(Nspan).ne.0) then
        CALL write_solid_span_field(xyzful/Lref,ele,time/Tref,nND,nEL,nND_max,nEL_max,Nspan,dspan,Lref,nFish)
    endif
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
            call cptMove(move(1:3),xyzful(NDref,1:3,iFish),[xGrid(IXref),yGrid(IYref),zGrid(IZref)],dh,MoveOutputIref(1:3))
            write(*,'(A,3F10.5)')' *Grid:',[xGrid(IXref),yGrid(IYref),zGrid(IZref)]
            write(*,'(A,3F10.5)')' *Body:',xyzful(NDref,1:3,iFish)
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

        !package n Fish to one
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iND,icount,iFish)
        do iFish=1,nFish
            if(iFish.eq.1)then
                icount = 0
            elseif(iFish.ge.2)then
                icount = sum(nND(1:iFish-1))
            endif
            do iND=1,nND(iFish)
                xyzful_all(iND+icount,1:6)  =  xyzful(iND,1:6,iFish)
                velful_all(iND+icount,1:6)  =  velful(iND,1:6,iFish)
                extful_all(iND+icount,1:6)  =  extful(iND,1:6,iFish)
            enddo
        enddo
        !$OMP END PARALLEL DO
        !compute force exerted on fluids
        if    (iForce2Body==1)then   !Same force as flow
            if (maxval(isFake_ful) .eq. 1) then
                call Beam_update_package(nFish,nND,xyzful,velful)
                CALL calculate_interaction_force(zDim,yDim,xDim,Beam_nEL,Beam_nND,Beam_ele,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                            Beam_xyzful,Beam_velful,Pbeta,ntolLBM,dtolLBM,force,Beam_extful,isUniformGrid)
                call Beam_unpackage(nFish,nND,extful)       
            elseif    (maxval(Nspan(:)) .eq. 0) then
                CALL calculate_interaction_force(zDim,yDim,xDim,nEL_all,nND_all,ele_all,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                        xyzful_all,velful_all,Pbeta,ntolLBM,dtolLBM,force,extful_all,isUniformGrid)
            else
                CALL calculate_interaction_force_quad(zDim,yDim,xDim,nEL_all,nND_all,ele_all,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                        xyzful_all,velful_all,Pbeta,ntolLBM,dtolLBM,force,extful_all,isUniformGrid,maxval(Nspan(:)),maxval(theta(:)),maxval(dspan(:)),boundaryConditions)
            endif
        elseif(iForce2Body==2)then   !stress force
        endif

        !compute volume force exerted on fluids
        CALL addVolumForc()

        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iND,icount,iFish)
        do iFish=1,nFish
            if(iFish.eq.1)then
                icount = 0
            else
                icount = sum(nND(1:iFish-1))
            endif
            do iND=1,nND(iFish)
                xyzful(iND,1:6,iFish) = xyzful_all(iND + icount,1:6)
                velful(iND,1:6,iFish) = velful_all(iND + icount,1:6)
                extful(iND,1:6,iFish) = extful_all(iND+icount,1:6)
            enddo
        enddo
        !$OMP END PARALLEL DO

        call date_and_time(VALUES=values1)
        write(*,*)'time for IBM step : ',CPUtime(values1)-CPUtime(values0)
        if(time/Tref >begForcDist .and. time/Tref <endForcDist) call forcDisturb() !force disturbance for instability

        call date_and_time(VALUES=values0)
        call cptForceR(dxmin,dymin,dzmin,nND,nND_max,xyzful,repful,nFish,dspan,Nspan)
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
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iND,iFish)
        do iFish=1,nFish
        if(iBodyModel(iFish)==1)then     ! rigid body
            !======================================================
            !prescribed motion
            !------------------------------------------------------
            !translational displacement
            XYZ(1:3,iFish)=XYZo(1:3,iFish)+XYZAmpl(1:3,iFish)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+XYZPhi(1:3,iFish))
            !rotational displacement
            AoA(1:3,iFish)=AoAo(1:3,iFish)+AoAAmpl(1:3,iFish)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+AoAPhi(1:3,iFish))
            call AoAtoTTT(AoA(1:3,iFish),TTTnxt(1:3,1:3,iFish))
            call get_angle_triad(TTT0(1:3,1:3,iFish),TTTnxt(1:3,1:3,iFish),AoAd(1,iFish),AoAd(2,iFish),AoAd(3,iFish))
            !given displacement
            do  iND=1,nND(iFish)
                xyzfulnxt(iND,1:3,iFish)=matmul(TTTnxt(1:3,1:3,iFish),xyzful00(iND,1:3,iFish))+XYZ(1:3,iFish)
                xyzfulnxt(iND,4:6,iFish)=AoAd(1:3,iFish)
            enddo
            xyzful(1:nND(iFish),1:6,iFish)=xyzfulnxt(1:nND(iFish),1:6,iFish)
            !------------------------------------------------------
            !translational velocity
            UVW(1:3,iFish) =-2.0*pi*Freq(iFish)*XYZAmpl(1:3,iFish)*dsin(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+XYZPhi(1:3,iFish))
            !rotational velocity
            WWW1(1:3,iFish)=-2.0*pi*Freq(iFish)*AoAAmpl(1:3,iFish)*dsin(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+AoAPhi(1:3,iFish))
            WWW2(1:3,iFish)=[WWW1(1,iFish)*dcos(AoA(2,iFish))+WWW1(3,iFish),    &
                             WWW1(1,iFish)*dsin(AoA(2,iFish))*dsin(AoA(3,iFish))+WWW1(2,iFish)*dcos(AoA(3,iFish)),   &
                             WWW1(1,iFish)*dsin(AoA(2,iFish))*dcos(AoA(3,iFish))-WWW1(2,iFish)*dsin(AoA(3,iFish))    ]
            WWW3(1:3,iFish)=matmul(TTTnxt(1:3,1:3,iFish),WWW2(1:3,iFish))
            !given velocity
            do  iND=1,nND(iFish)
                velful(iND,1:3,iFish)=[WWW3(2,iFish)*xyzful(iND,3,iFish)-WWW3(3,iFish)*xyzful(iND,2,iFish),    &
                                       WWW3(3,iFish)*xyzful(iND,1,iFish)-WWW3(1,iFish)*xyzful(iND,3,iFish),    &
                                       WWW3(1,iFish)*xyzful(iND,2,iFish)-WWW3(2,iFish)*xyzful(iND,1,iFish)    ]&
                                       + UVW(1:3,iFish)
                velful(iND,4:6,iFish)=WWW3(1:3,iFish)
            enddo
            !-------------------------------------------------------
        elseif(iBodyModel(iFish)==2)then !elastic model
            !translational displacement
            XYZ(1:3,iFish)=XYZo(1:3,iFish)+XYZAmpl(1:3,iFish)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+XYZPhi(1:3,iFish))
            !rotational displacement
            AoA(1:3,iFish)=AoAo(1:3,iFish)+AoAAmpl(1:3,iFish)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+AoAPhi(1:3,iFish))
            call AoAtoTTT(AoA(1:3,iFish),TTTnxt(1:3,1:3,iFish))
            call get_angle_triad(TTT0(1:3,1:3,iFish),TTTnxt(1:3,1:3,iFish),AoAd(1,iFish),AoAd(2,iFish),AoAd(3,iFish))
            !given displacement
            do  iND=1,nND(iFish)
                xyzfulnxt(iND,1:3,iFish)=matmul(TTTnxt(1:3,1:3,iFish),xyzful00(iND,1:3,iFish))+XYZ(1:3,iFish)
                xyzfulnxt(iND,4:6,iFish)=AoAd(1:3,iFish)
            enddo
            !-------------------------------------------------------
            !displacement condition
            vBC(1:nND(iFish),1:6,iFish) = xyzfulnxt(1:nND(iFish),1:6,iFish) - xyzful(1:nND(iFish),1:6,iFish)
            !loading vector
            lodful(1:nND(iFish),1:6,iFish) = extful(1:nND(iFish),1:6,iFish) + grav(1:nND(iFish),1:6,iFish) + repful(1:nND(iFish),1:6,iFish)
            !-----------------------------------------
            CALL structure_solver(jBC(1:nND(iFish),1:6,iFish),vBC(1:nND(iFish),1:6,iFish),ele(1:nEL(iFish),1:5,iFish), &
                                  nloc(1:nND(iFish)*6,iFish),nprof(1:nND(iFish)*6,iFish),nprof2(1:nND(iFish)*6,iFish), &
                                  prop(1:nMT(iFish),1:10,iFish),mss(1:nND(iFish)*6,iFish), &
                                  xyzful0(1:nND(iFish),1:6,iFish),xyzful(1:nND(iFish),1:6,iFish),dspful(1:nND(iFish),1:6,iFish), &
                                  velful(1:nND(iFish),1:6,iFish),accful(1:nND(iFish),1:6,iFish),lodful(1:nND(iFish),1:6,iFish),  &
                                  subdeltat,dampK,dampM,  &
                                  triad_nn(1:3,1:3,1:nND(iFish),iFish),triad_ee(1:3,1:3,1:nEL(iFish),iFish), &
                                  triad_n1(1:3,1:3,1:nEL(iFish),iFish),triad_n2(1:3,1:3,1:nEL(iFish),iFish), &
                                  nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),nBD(iFish),nSTF(iFish),NewmarkGamma,NewmarkBeta,dtolFEM,ntolFEM,    &
                                  nFish,iFish,FishInfo)
        else
            stop 'no define body model'
        endif
        enddo !do iFish=1,nFish
        !$OMP END PARALLEL DO
        enddo !do isubstep=1,numsubstep
        call date_and_time(VALUES=values1)
        write(*,*)'time for FEM:',CPUtime(values1)-CPUtime(values0)
        !----------------------------------------------------------------------
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)then
            do iFish=1,nFish
                write(*,'(A,I5.5)')' Fish number: ', int(FishInfo(1,iFish))
                write(*,'(A,I5.5,A,E20.10)')' iterFEM = ',int(FishInfo(2,iFish)),'    dmaxFEM = ',FishInfo(3,iFish)
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
                write(*,'(A,3D15.5)')" forceDre: ",sum(extful(1:nND(iFish),1:3,iFish),1)/Fref
                write(*,'(A,3D15.5)')" accCentM: ",sum(accful(1:nND(iFish),1:3,iFish)*mssful(1:nND(iFish),1:3,iFish),1)/sum(mssful(1:nND(iFish),1:3,iFish),1)/Aref
                write(*,'(A,3D15.5)')" velCentM: ",sum(velful(1:nND(iFish),1:3,iFish)*mssful(1:nND(iFish),1:3,iFish),1)/sum(mssful(1:nND(iFish),1:3,iFish),1)/Uref
                write(*,'(A,3D15.5)')" xyzCentM: ",sum(xyzful(1:nND(iFish),1:3,iFish)*mssful(1:nND(iFish),1:3,iFish),1)/sum(mssful(1:nND(iFish),1:3,iFish),1)/Lref
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
                CALL write_solid_field(xyzful/Lref,velful/Uref,accful/Aref,extful/Fref,repful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,nFish)
                call Beam_write_solid_fake_field(nFish,Lref,time,Tref)
                if (maxval(Nspan).ne.0) then
                    CALL write_solid_span_field(xyzful/Lref,ele,time/Tref,nND,nEL,nND_max,nEL_max,Nspan,dspan,Lref,nFish)
                endif
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
