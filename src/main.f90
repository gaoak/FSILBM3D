!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    copyright@ RuNanHua
!    The main program, 3D Lattice Boltzmann Method
!    flow past a flexible plate (uniform flow past a flag or flapping plate, and so on)
!    flexible Plates or shells move in fluid.  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    PROGRAM main
    USE simParam
    use omp_lib
    implicit none
    integer:: iND,isubstep
    real(8):: temp(3),Pbetatemp
    logical alive
    !time_and_date
    integer,dimension(8) :: values0,values1,values_s,values_e
    CALL read_file()
    CALL allocate_solid_memory()
    CALL allocate_fluid_memory()
    CALL calculate_LB_params()
    CALL write_params()
    CALL calculate_MRTM_params()
    
    !np=OMP_get_num_procs() 
    
    call OMP_set_num_threads(npsize)
    write(*,*)'npsize=', npsize

    if(isRelease==1)then
        timeOutInfo=timeOutInfo*5
        write(*,*)'Release'
    endif
    Pbetatemp=Pbeta  
    deltat = dt  !set time step of solid deltat the same as fluid time step
!    ===============================================================================================
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
        time=0.0d0
        step=0
        CALL initialize_solid() 
        CALL initialize_flow()                  
!    ===============================================================================================
        if(ismovegrid==1)then
            call cptIref(NDref,IXref,IYref,IZref,nND,xDim,yDim,zDim,xyzful(1:nND,1:3),xGrid,yGrid,zGrid,Xref,Yref,Zref)
        endif  
    endif
!    ===============================================================================================
    if(step==0)    CALL wrtInfoTitl()   !
!    ===============================================================================================
    CALL calculate_macro_quantities()
    CALL write_flow_field(1)
    CALL write_solid_field(xyzful/Lref,velful/Uref,accful/Aref,extful/Fref,ele,time/Tref,nND,nEL,0)
    CALL write_image() 
!==================================================================================================
!==================================================================================================
!==================================================================================================      
    write(*,*)'time loop'
    do while(time/Tref < timeSimTotl)
        call date_and_time(VALUES=values_s)
        time=time+dt
        step=step+1
        write(*,'(A)')' ======================================================================='
        write(*,'(A,I6,A,F15.10)')' step:',step,' time/Tref:',time/Tref
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)write(*,'(A)')' ----------------fluid solver----------------'
        !solve fluid
        call date_and_time(VALUES=values0)        
        CALL streaming_step()
        call date_and_time(VALUES=values1)
        write(*,*)'time for streaming_step:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !******************************************************************************************
        !DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304,fluid=0, wall=200
        if    (iStreamModel==1)then
        xMinBC=DirecletUP
        xMaxBC=DirecletUP
        yMinBC=DirecletUP
        yMaxBC=DirecletUP
        zMinBC=DirecletUP
        zMaxBC=DirecletUP
        CALL set_other_farfld_BCs()
        elseif(iStreamModel==2)then
        CALL set_equilibrium_farfld_BC()
        else
            stop
            write(*,*)'no such type LBMBC'
        endif
        !******************************************************************************************
        if(ismovegrid==1)then ! move fluid grid
            call cptMove(move(1:3),xyzful(NDref,1:3),[xGrid(IXref),yGrid(IYref),zGrid(IZref)],[dx,dy,dz])
            write(*,'(A,3F10.5)')' Grid:',[xGrid(IXref),yGrid(IYref),zGrid(IZref)]
            write(*,'(A,3F10.5)')' Body:',xyzful(NDref,1:3)

            if(isMoveDimX==1)CALL movGrid(1,move(1))
            if(isMoveDimY==1)CALL movGrid(2,move(2))
            if(isMoveDimZ==1)CALL movGrid(3,move(3))
        endif
        !******************************************************************************************
        !DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304,fluid=0, wall=200
        if    (iStreamModel==1)then
        xMinBC=DirecletUP
        xMaxBC=DirecletUP
        yMinBC=DirecletUP
        yMaxBC=DirecletUP
        zMinBC=DirecletUP
        zMaxBC=DirecletUP
        CALL set_other_farfld_BCs()
        elseif(iStreamModel==2)then
        CALL set_equilibrium_farfld_BC()
        else
            stop
            write(*,*)'no such type LBMBC'
        endif
        !******************************************************************************************        
        CALL calculate_macro_quantities()
        !******************************************************************************************
        !******************************************************************************************
        call date_and_time(VALUES=values0)
        Pbeta=(1.0d0-dexp(-5.0d0/Pramp*time/Tref))*Pbetatemp
        CALL calculate_interaction_force(zDim,yDim,xDim,nEL,nND,ele,dx,dy,dz,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                       xyzful,velful,xyzfulIB,Palpha,Pbeta,ntolLBM,dtolLBM,force,extful1)
!        CALL cptStrs(zDim,yDim,xDim,nEL,nND,ele,dh,dx,dy,dz,mu,2.50d0,uuu,prs,xGrid,yGrid,zGrid,xyzful,extful2)
        call date_and_time(VALUES=values1)
        write(*,*)'time for IB:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)

        if(time/Tref >begForcDist .and. time/Tref <endForcDist) call forcDisturb() !force disturbance for instability
        
        call date_and_time(VALUES=values0)  
        CALL collision_step()
        call date_and_time(VALUES=values1)
        write(*,*)'time for collision_step:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !exert force to solid, two types of force: penalty and fluid stress
!        if    (iForce2Body==1)then   !Same force as flow
            extful(1:nND,1:6)  = extful1(1:nND,1:6)
!        elseif(iForce2Body==2)then   !stress force
!            extful(1:nND,1:6)  = extful2(1:nND,1:6)
!        else
!             stop 'no define force to body ' 
!        endif                           
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        !solve solid
        if(isRelease/=1)write(*,'(A)')' ----------------solid solver----------------'
        subdeltat=deltat/numsubstep
        do isubstep=1,numsubstep

        if(iBodyModel==1)then     ! rigid body
            !======================================================
            !prescribed motion
            !------------------------------------------------------
            !translational displacement
            XYZ(1:3)=XYZo(1:3)+XYZAmpl(1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+XYZPhi(1:3))
            !rotational displacement
            AoA(1:3)=AoAo(1:3)+AoAAmpl(1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+AoAPhi(1:3))       
            call AoAtoTTT(AoA,TTTnxt)
            call get_angle_triad(TTT0,TTTnxt,AoAd(1),AoAd(2),AoAd(3))
            !given displacement
            do  iND=1,nND
                xyzfulnxt(iND,1:3)=matmul(TTTnxt,xyzful00(iND,1:3))+XYZ(1:3)
                xyzfulnxt(iND,4:6)=AoAd(1:3)
            enddo
            xyzful(1:nND,1:6)=xyzfulnxt(1:nND,1:6)
            !------------------------------------------------------
            !translational velocity
            UVW(1:3) =-2.0*pi*Freq*XYZAmpl(1:3)*dsin(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+XYZPhi(1:3))
            !rotational velocity
            WWW1(1:3)=-2.0*pi*Freq*AoAAmpl(1:3)*dsin(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+AoAPhi(1:3))
            WWW2(1:3)=[ WWW1(1)*dcos(AoA(2))+WWW1(3),    &                       
                        WWW1(1)*dsin(AoA(2))*dsin(AoA(3))+WWW1(2)*dcos(AoA(3)),     &
                        WWW1(1)*dsin(AoA(2))*dcos(AoA(3))-WWW1(2)*dsin(AoA(3))    &
                      ]
            WWW3(1:3)=matmul(TTTnxt,WWW2(1:3))
            !given velocity
            do  iND=1,nND
                velful(iND,1:3)=[   WWW3(2)*xyzful(iND,3)-WWW3(3)*xyzful(iND,2),    &
                                    WWW3(3)*xyzful(iND,1)-WWW3(1)*xyzful(iND,3),    &
                                    WWW3(1)*xyzful(iND,2)-WWW3(2)*xyzful(iND,1)     &
                                ]+UVW(1:3)
                velful(iND,4:6)=WWW3(1:3)
            enddo
            !-------------------------------------------------------
        elseif(iBodyModel==2)then !elastic model     
            !translational displacement 
            call date_and_time(VALUES=values0)              
            XYZ(1:3)=XYZo(1:3)+XYZAmpl(1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+XYZPhi(1:3))
            !rotational displacement
            AoA(1:3)=AoAo(1:3)+AoAAmpl(1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+AoAPhi(1:3))         
            call AoAtoTTT(AoA,TTTnxt)
            call get_angle_triad(TTT0,TTTnxt,AoAd(1),AoAd(2),AoAd(3))
            !given displacement
            do  iND=1,nND
                xyzfulnxt(iND,1:3)=matmul(TTTnxt,xyzful00(iND,1:3))+XYZ(1:3)
                xyzfulnxt(iND,4:6)=AoAd(1:3)
            enddo
            !-------------------------------------------------------
            !displacement condition
            vBC(   1:nND,1:6) = xyzfulnxt(1:nND,1:6)-xyzful(1:nND,1:6)
            !loading vector
            lodful(1:nND,1:6) = extful(1:nND,1:6)+grav(1:nND,1:6)
            !-----------------------------------------
            CALL structure_solver( jBC,vBC,ele,nloc,nprof,nprof2,prop,mss,xyzful0,xyzful,dspful,velful,accful,lodful,subdeltat,dampK,dampM,  &
                        triad_nn,triad_ee,triad_e0,triad_n1,triad_n2,triad_n3,nND,nEL,nEQ,nMT,nBD,nSTF,NewmarkGamma,NewmarkBeta,dtolFEM,ntolFEM)
            !----------------------------------------------------------------------
            call date_and_time(VALUES=values1)
            write(*,*)'time for FEM:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        else
            stop 'no define body model'
        endif
        enddo
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)write(*,'(A)')' ----------------post process----------------'
        if(isRelease/=1)then
        write(*,'(A,3D15.5)')" forceDirect:",sum(extful(1:nND,1:3),1)/Fref       
        !write(*,'(A,3D15.5)')" forceStress:",sum(extful2(1:nND,1:3),1)/Fref
        write(*,'(A,3D15.5)')" accCentM:",sum(accful(1:nND,1:3)*mssful(1:nND,1:3),1)/sum(mssful(1:nND,1:3),1)/Aref
        write(*,'(A,3D15.5)')" velCentM:",sum(velful(1:nND,1:3)*mssful(1:nND,1:3),1)/sum(mssful(1:nND,1:3),1)/Uref
        write(*,'(A,3D15.5)')" xyzCentM:",sum(xyzful(1:nND,1:3)*mssful(1:nND,1:3),1)/sum(mssful(1:nND,1:3),1)/Lref
        endif
                              
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        if(DABS(time/Tref-timeOutTemp*NINT(time/Tref/timeOutTemp)) <= 0.5*dt/Tref)then
            CALL write_checkpoint_file()
        endif
                          
        if(DABS(time/Tref-timeOutBody*NINT(time/Tref/timeOutBody)) <= 0.5*dt/Tref)then
            if(iBodyModel==2 .or. maxval(dabs(XYZAmpl))>eps .or. maxval(dabs(AoAAmpl))>eps )then
            CALL write_solid_field(xyzful/Lref  ,velful/Uref, accful/Aref,extful/Fref,ele,time/Tref,nND,nEL,0)
            if (dabs(Palpha)>eps) then
            CALL write_solid_field(xyzfulIB/Lref,velful/Uref, accful/Aref,extful/Fref,ele,time/Tref,nND,nEL,1)
            endif
            endif
        endif



        if(DABS(time/Tref-timeOutFlow*NINT(time/Tref/timeOutFlow)) <= 0.5*dt/Tref)then
            CALL write_flow_field(1)
            if(isRelease/=1) then
                CALL write_flow_field(0) 
            endif              
        endif

        !if(DABS(time/Tref-timeOutFlow*NINT(time/Tref/timeOutFlow)) <= 0.5*dt/Tref .and. time/Tref>=timeOutFlBg .and. time/Tref<=timeOutFlEd)then            
        !    CALL write_flow_field(1)
        !endif

        !if(DABS(time/Tref-timeOutBody*NINT(time/Tref/timeOutBody)) <= 0.5*dt/Tref .and. time/Tref>=timeOutFlBg .and. time/Tref<=timeOutFlEd)then
            !CALL write_flow_slice(xDim,yDim,zDim,XGrid,yGrid,zGrid,Lref,Uref,denIn,prs,uuu,time/Tref-timeOutFlBg,3)
        !endif
     
        if(DABS(time/Tref-timeOutInfo*NINT(time/Tref/timeOutInfo)) <= 0.5*dt/Tref)then
            CALL wrtInfo()
        endif
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        call date_and_time(VALUES=values_e)
        write(*,*)'time for total:',(values_e(6)*60.+values_e(7)*1.+values_e(8)*0.001)-(values_s(6)*60.+values_s(7)*1.+values_s(8)*0.001)
    enddo
    END PROGRAM main

