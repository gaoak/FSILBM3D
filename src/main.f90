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
    integer:: iND,isubstep,iFish,x,y,z
    real(8), allocatable:: FishInfo(:,:)
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
    
    allocate(FishInfo(1:nFish,1:3))

    !np=OMP_get_num_procs() 
    
    call OMP_set_num_threads(npsize)
    write(*,*)'npsize=', npsize

    if(isRelease==1)then
        timeOutInfo=timeOutInfo*2
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
            iFish = 1
            call cptIref(NDref,IXref,IYref,IZref,nND(iFish),xDim,yDim,zDim,xyzful(iFish,1:nND(iFish),1:3),xGrid,yGrid,zGrid,Xref,Yref,Zref)
        endif  
    endif
!    ===============================================================================================
    if(step==0)    CALL wrtInfoTitl()
!    ===============================================================================================
    CALL calculate_macro_quantities()
    CALL write_flow_field(1)
    CALL write_solid_field(nFish,xyzful/Lref,velful/Uref,accful/Aref,extful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,0)
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
        write(*,'(A,I6,A,F15.10)')' step:',step,'    time/Tref:',time/Tref
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        write(*,'(A)')' ----------------------fluid solver----------------------'
        !solve fluid
        call date_and_time(VALUES=values0)        
        CALL streaming_step()
        call date_and_time(VALUES=values1)
        write(*,*)'time for streaming:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
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
            call cptMove(move(1:3),xyzful(iFish,NDref,1:3),[xGrid(IXref),yGrid(IYref),zGrid(IZref)],[dx,dy,dz])
            write(*,'(A,3F10.5)')' *Grid:',[xGrid(IXref),yGrid(IYref),zGrid(IZref)]
            write(*,'(A,3F10.5)')' *Body:',xyzful(iFish,NDref,1:3)

            if(isMoveDimX==1)CALL movGrid(1,move(1))
            if(isMoveDimY==1)CALL movGrid(2,move(2))
            if(isMoveDimZ==1)CALL movGrid(3,move(3))
        endif
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
        !******************************************************************************************        
        CALL calculate_macro_quantities()
        !******************************************************************************************
        !******************************************************************************************
        call date_and_time(VALUES=values0)
        Pbeta=(1.0d0-dexp(-5.0d0/Pramp*time/Tref))*Pbetatemp

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
        do iFish=1,nFish
        if(iFish.eq.1)then
            do iND=1,nND(iFish)
               xyzful_all(iND,1:6)   =  xyzful(iFish,iND,1:6)
               velful_all(iND,1:6)   =  velful(iFish,iND,1:6)
               xyzfulIB_all(iND,1:6) =  xyzfulIB(iFish,iND,1:6)
               extful1_all(iND,1:6)  =  extful(iFish,iND,1:6)
               extful2_all(iND,1:6)  =  extful(iFish,iND,1:6)
            enddo
        elseif(iFish.ge.2)then
            do iND=1,nND(iFish)
               xyzful_all(iND+sum(nND(1:iFish-1)),1:6)   =  xyzful(iFish,iND,1:6)
               velful_all(iND+sum(nND(1:iFish-1)),1:6)   =  velful(iFish,iND,1:6)
               xyzfulIB_all(iND+sum(nND(1:iFish-1)),1:6) =  xyzfulIB(iFish,iND,1:6)
               extful1_all(iND+sum(nND(1:iFish-1)),1:6)  =  extful(iFish,iND,1:6)
               extful2_all(iND+sum(nND(1:iFish-1)),1:6)  =  extful(iFish,iND,1:6)
            enddo
        endif    

        enddo
        !compute force exerted on fluids
        if    (iForce2Body==1)then   !Same force as flow
        CALL calculate_interaction_force(zDim,yDim,xDim,nEL_all,nND_all,ele_all,dx,dy,dz,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                       xyzful_all,velful_all,xyzfulIB_all,Palpha,Pbeta,ntolLBM,dtolLBM,force,extful1_all)
        elseif(iForce2Body==2)then   !stress force
        CALL cptStrs(zDim,yDim,xDim,nEL_all,nND_all,ele_all,dh,dx,dy,dz,mu,2.50d0,uuu,prs,xGrid,yGrid,zGrid,xyzful_all,extful2_all)
        else
             stop 'no define force to body ' 
        endif 

        ! unpackage  
        do iFish=1,nFish
            if(iFish.eq.1)then
               do iND=1,nND(iFish)
                 xyzful(iFish,iND,1:6)   =  xyzful_all(iND,1:6) 
                 velful(iFish,iND,1:6)   =  velful_all(iND,1:6)
                 xyzfulIB(iFish,iND,1:6) =  xyzfulIB_all(iND,1:6)
               enddo
            else
               do iND=1,nND(iFish)
                 xyzful(iFish,iND,1:6)   =  xyzful_all(iND + sum(nND(1:iFish-1)),1:6) 
                 velful(iFish,iND,1:6)   =  velful_all(iND + sum(nND(1:iFish-1)),1:6)
                 xyzfulIB(iFish,iND,1:6) =  xyzfulIB_all(iND + sum(nND(1:iFish-1)),1:6)
               enddo
            endif
        enddo 

        call date_and_time(VALUES=values1)
        write(*,*)'time for IBM step : ',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        if(time/Tref >begForcDist .and. time/Tref <endForcDist) call forcDisturb() !force disturbance for instability
        
        call date_and_time(VALUES=values0)  
        call cptForceR(nFish,dxmin,dymin,dzmin,nND,nND_max,nEL,nEL_max,ele,xyzful,repful)
        call date_and_time(VALUES=values1)
        write(*,*)'time for Lubforce :',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)        

        call date_and_time(VALUES=values0)  
        CALL collision_step()
        call date_and_time(VALUES=values1)
        write(*,*)'time for collision: ',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !exert force to solid, two types of force: penalty and fluid stress
        if (iForce2Body==1)then   !Same force as flow
            do iFish=1,nFish
            if(iFish.eq.1)then
                do iND=1,nND(iFish)
                extful(iFish,iND,1:6) = extful1_all(iND,1:6)
                enddo
            else
                do iND=1,nND(iFish)
                extful(iFish,iND,1:6) = extful1_all(iND+sum(nND(1:iFish-1)),1:6) 
                enddo
            endif
            enddo
        elseif(iForce2Body==2)then   !stress force
            do iFish=1,nFish
            if(iFish.eq.1)then
                do iND=1,nND(iFish)
                extful(iFish,iND,1:6) = extful2_all(iND,1:6)
                enddo
            else
                do iND=1,nND(iFish)
                extful(iFish,iND,1:6) = extful2_all(iND+sum(nND(1:iFish-1)),1:6) 
                enddo
            endif
            enddo
        else
             stop 'no define force to body ' 
        endif                           
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        !solve solid
        if(isRelease/=1)write(*,'(A)')' ----------------------solid solver----------------------'
        subdeltat=deltat/numsubstep
        do isubstep=1,numsubstep
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iFish,iND)
        do iFish=1,nFish
        if(iBodyModel(iFish)==1)then     ! rigid body
            !======================================================
            !prescribed motion
            !------------------------------------------------------
            !translational displacement
            XYZ(iFish,1:3)=XYZo(iFish,1:3)+XYZAmpl(iFish,1:3)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+XYZPhi(iFish,1:3))
            !rotational displacement
            AoA(iFish,1:3)=AoAo(iFish,1:3)+AoAAmpl(iFish,1:3)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+AoAPhi(iFish,1:3))       
            call AoAtoTTT(AoA(iFish,1:3),TTTnxt(iFish,1:3,1:3))
            call get_angle_triad(TTT0(iFish,1:3,1:3),TTTnxt(iFish,1:3,1:3),AoAd(iFish,1),AoAd(iFish,2),AoAd(iFish,3))
            !given displacement
            do  iND=1,nND(iFish)
                xyzfulnxt(iFish,iND,1:3)=matmul(TTTnxt(iFish,1:3,1:3),xyzful00(iFish,iND,1:3))+XYZ(iFish,1:3)
                xyzfulnxt(iFish,iND,4:6)=AoAd(iFish,1:3)
            enddo
            xyzful(iFish,1:nND(iFish),1:6)=xyzfulnxt(iFish,1:nND(iFish),1:6)
            !------------------------------------------------------
            !translational velocity
            UVW(iFish,1:3) =-2.0*pi*Freq(iFish)*XYZAmpl(iFish,1:3)*dsin(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+XYZPhi(iFish,1:3))
            !rotational velocity
            WWW1(iFish,1:3)=-2.0*pi*Freq(iFish)*AoAAmpl(iFish,1:3)*dsin(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+AoAPhi(iFish,1:3))
            WWW2(iFish,1:3)=[WWW1(iFish,1)*dcos(AoA(iFish,2))+WWW1(iFish,3),    &                       
                             WWW1(iFish,1)*dsin(AoA(iFish,2))*dsin(AoA(iFish,3))+WWW1(iFish,2)*dcos(AoA(iFish,3)),   &
                             WWW1(iFish,1)*dsin(AoA(iFish,2))*dcos(AoA(iFish,3))-WWW1(iFish,2)*dsin(AoA(iFish,3))    ]
            WWW3(iFish,1:3)=matmul(TTTnxt(iFish,1:3,1:3),WWW2(iFish,1:3))
            !given velocity
            do  iND=1,nND(iFish)
                velful(iFish,iND,1:3)=[WWW3(iFish,2)*xyzful(iFish,iND,3)-WWW3(iFish,3)*xyzful(iFish,iND,2),    &
                                       WWW3(iFish,3)*xyzful(iFish,iND,1)-WWW3(iFish,1)*xyzful(iFish,iND,3),    &
                                       WWW3(iFish,1)*xyzful(iFish,iND,2)-WWW3(iFish,2)*xyzful(iFish,iND,1)    ]&
                                       + UVW(iFish,1:3)
                velful(iFish,iND,4:6)=WWW3(iFish,1:3)
            enddo
            !-------------------------------------------------------
        elseif(iBodyModel(iFish)==2)then !elastic model     
            !translational displacement 
            call date_and_time(VALUES=values0)              
            XYZ(iFish,1:3)=XYZo(iFish,1:3)+XYZAmpl(iFish,1:3)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+XYZPhi(iFish,1:3))
            !rotational displacement
            AoA(iFish,1:3)=AoAo(iFish,1:3)+AoAAmpl(iFish,1:3)*dcos(2.0*pi*Freq(iFish)*(time-deltat+isubstep*subdeltat)+AoAPhi(iFish,1:3))        
            call AoAtoTTT(AoA(iFish,1:3),TTTnxt(iFish,1:3,1:3))
            call get_angle_triad(TTT0(iFish,1:3,1:3),TTTnxt(iFish,1:3,1:3),AoAd(iFish,1),AoAd(iFish,2),AoAd(iFish,3))
            !given displacement
            do  iND=1,nND(iFish)
                xyzfulnxt(iFish,iND,1:3)=matmul(TTTnxt(iFish,1:3,1:3),xyzful00(iFish,iND,1:3))+XYZ(iFish,1:3)
                xyzfulnxt(iFish,iND,4:6)=AoAd(iFish,1:3)
            enddo
            !-------------------------------------------------------
            !displacement condition
            vBC(iFish,1:nND(iFish),1:6) = xyzfulnxt(iFish,1:nND(iFish),1:6) - xyzful(iFish,1:nND(iFish),1:6)
            !loading vector
            lodful(iFish,1:nND(iFish),1:6) = extful(iFish,1:nND(iFish),1:6) + grav(iFish,1:nND(iFish),1:6) + repful(iFish,1:nND(iFish),1:6)
            !-----------------------------------------
            CALL structure_solver(jBC(iFish,1:nND(iFish),1:6),vBC(iFish,1:nND(iFish),1:6),ele(iFish,1:nEL(iFish),1:5), &
                                  nloc(iFish,1:nND(iFish)*6),nprof(iFish,1:nND(iFish)*6),nprof2(iFish,1:nND(iFish)*6), &
                                  prop(iFish,1:nMT(iFish),1:10),mss(iFish,1:nND(iFish)*6), &
                                  xyzful0(iFish,1:nND(iFish),1:6),xyzful(iFish,1:nND(iFish),1:6),dspful(iFish,1:nND(iFish),1:6), &
                                  velful(iFish,1:nND(iFish),1:6),accful(iFish,1:nND(iFish),1:6),lodful(iFish,1:nND(iFish),1:6),  &
                                  subdeltat,dampK,dampM,  &
                                  triad_nn(iFish,1:3,1:3,1:nND(iFish)),triad_ee(iFish,1:3,1:3,1:nEL(iFish)),triad_e0(iFish,1:3,1:3,1:nEL(iFish)), &
                                  triad_n1(iFish,1:3,1:3,1:nEL(iFish)),triad_n2(iFish,1:3,1:3,1:nEL(iFish)),triad_n3(iFish,1:3,1:3,1:nEL(iFish)), &
                                  nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),nBD(iFish),nSTF(iFish),NewmarkGamma,NewmarkBeta,dtolFEM,ntolFEM,    &
                                  nFish,iFish,FishInfo)
        else
            stop 'no define body model'
        endif
        enddo !do iFish=1,nFish
        !$OMP END PARALLEL DO
        enddo !do isubstep=1,numsubstep
        !----------------------------------------------------------------------
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)then
        do iFish=1,nFish
        write(*,'(A,I5.5)')' Fish number: ', int(FishInfo(iFish,1))
        write(*,'(A,I5.5,A,D20.10)')' iterFEM = ',int(FishInfo(iFish,2)),'    dmaxFEM = ',FishInfo(iFish,3)
        enddo
        endif
        write(*,'(A)')' --------------------------------------------------------'
        call date_and_time(VALUES=values1)
        write(*,*)'max time for FEM  : ',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)write(*,'(A)')' ----------------------post process----------------------'
        if(isRelease/=1)then
        do iFish=1,nFish
        write(*,'(A,I5.5)')' Fish number: ',iFish
        write(*,'(A,3D15.5)')" forceDre: ",sum(extful(iFish,1:nND(iFish),1:3),1)/Fref       
        !write(*,'(A,3D15.5)')" forceStress:",sum(extful2(1:nND,1:3),1)/Fref
        write(*,'(A,3D15.5)')" accCentM: ",sum(accful(iFish,1:nND(iFish),1:3)*mssful(iFish,1:nND(iFish),1:3),1)/sum(mssful(iFish,1:nND(iFish),1:3),1)/Aref
        write(*,'(A,3D15.5)')" velCentM: ",sum(velful(iFish,1:nND(iFish),1:3)*mssful(iFish,1:nND(iFish),1:3),1)/sum(mssful(iFish,1:nND(iFish),1:3),1)/Uref
        write(*,'(A,3D15.5)')" xyzCentM: ",sum(xyzful(iFish,1:nND(iFish),1:3)*mssful(iFish,1:nND(iFish),1:3),1)/sum(mssful(iFish,1:nND(iFish),1:3),1)/Lref
        enddo
        endif                
        !******************************************************************************************
        !******************************************************************************************
        !******************************************************************************************
        if(DABS(time/Tref-timeOutTemp*NINT(time/Tref/timeOutTemp)) <= 0.5*dt/Tref)then
            CALL write_checkpoint_file()
        endif
                          
        if(DABS(time/Tref-timeOutBody*NINT(time/Tref/timeOutBody)) <= 0.5*dt/Tref)then
            CALL write_solid_field(nFish,xyzful/Lref  ,velful/Uref,accful/Aref,extful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,0)
            if (dabs(Palpha)>eps) then
            CALL write_solid_field(nFish,xyzfulIB/Lref,velful/Uref,accful/Aref,extful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,1)
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
        write(*,'(A)')' --------------------------------------------------------'
        call date_and_time(VALUES=values_e)
        write(*,*)'time for one step :',(values_e(6)*60.+values_e(7)*1.+values_e(8)*0.001)-(values_s(6)*60.+values_s(7)*1.+values_s(8)*0.001)
    enddo  
    stop 'stop'
    END PROGRAM main
