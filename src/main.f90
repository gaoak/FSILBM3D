!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    PROGRAM main
    USE simParam
    use omp_lib
    implicit none
    integer:: iND,itemp,isubstep,iFish,x,y
    real(8), allocatable:: FishInfo(:,:)
    real(8):: temp(3),Pbetatemp
    real(8),allocatable:: force_temp(:,:,:)
    real, dimension(2) :: tarray
    logical alive
    
    !time_and_date
    integer,dimension(8) :: values0,values1,values_s,values_e
    !call date_and_time(VALUES=values0)
    !write(*,*)values0(:)
    !stop
    !!!!!!!!!
 
    CALL redFile()         
    CALL memFlow()
    CALL memBody()   
    CALL cptPara()
    CALL wrtPara()
    CALL cptMRTM()
    write(*,*) 'Code last updated on 2022/11/17'
    
    allocate(FishInfo(1:nFish,1:3))

    !np=OMP_get_num_procs() 
    call OMP_set_num_threads(npsize)
    write(*,*)'npsize=', npsize
    
    allocate(force_temp(1:yDim,1:xDim,1:SpcDim))
    
    if(isRelease==1)then
        !timeOutInfo=timeOutInfo*5
        timeOutInfo=timeOutInfo*2
        write(*,*)'Release'
    endif 
    Pbetatemp=Pbeta
    deltat = dt
!    ===============================================================================================
    inquire(file='./DatTemp/conwr.dat', exist=alive)
    if (isConCmpt==1 .and. alive)then
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'continue compute!!!!!!!!!!!!!!!!!!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        call redCont()
    else
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'new      compute!!!!!!!!!!!!!!!!!!'
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        
        time=0.0d0
        step=0
        CALL intBody()
        CALL intFlow()
     
        if(ismovegrid==1)then
            iFish=1
            call cptIref(NDref,IXref,IYref,nND(iFish),xDim,yDim,xyzful(iFish,1:nND(iFish),1:2),xGrid,yGrid,Xref,Yref)
        endif
        
        CALL wrtInfoTitl()
        
    endif
!    ===============================================================================================
    if(step==0)  CALL wrtInfoTitl()
!    ===============================================================================================
    CALL cptMacr()
    CALL writeFlow(1)
    CALL writeBody3(nFish,xyzful/Lref,velful/Uref,accful/Aref,extful/Fref,repful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,0)
    CALL writeImag() 
   
    do while(time/Tref<=timeSimTotl)
        call date_and_time(VALUES=values_s)
        time=time+dt
        step=step+1
        write(*,'(A)')' ======================================================================='
        write(*,'(A,I6,A,F15.10)')' Step:',step,' time/Tref:',time/Tref
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1) write(*,'(A)')' --------------------fluid solver--------------------'
        call date_and_time(VALUES=values0)
        !fluid solver
        CALL stream()
        call date_and_time(VALUES=values1)
        write(*,*)'time for streaming:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !******************************************************************************************
        !DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304,MovingWall=305,fluid=0, wall=200, Movingwall=201
        if(iStreamModel==1)then !uniform grid,STLBM
            xMinBC=boundaryConditions(1)
            xMaxBC=boundaryConditions(2)
            yMinBC=boundaryConditions(3)
            yMaxBC=boundaryConditions(4)
            CALL setBund2()
        elseif(iStreamModel==2)then !non-uniform grid,ISLBM
            CALL setBund()
        else
            stop
            write(*,*)'no such type LBMBC'
        endif
        !******************************************************************************************
        if(ismovegrid==1)then
            iFish=1
            call cptMove(move(1:2),xyzful(iFish,NDref,1:2),[xGrid(IXref),yGrid(IYref)],[dxmin,dymin])
            if(isRelease/=1)    write(*,'(A,3F10.5)')'  *Grid:',[xGrid(IXref),yGrid(IYref)]
            if(isRelease/=1)    write(*,'(A,3F10.5)')'  *Body:',xyzful(iFish,NDref,1:2)

            CALL movGrid(iMoveDim,move(iMoveDim))
        endif
        !******************************************************************************************
        !DirecletUP=300,DirecletUU=301,Advection1=302,Advection2=303,Periodic=304,MovingWall=305,fluid=0, wall=200, Movingwall=201
        if(iStreamModel==1)then !uniform grid,STLBM
            xMinBC=boundaryConditions(1)
            xMaxBC=boundaryConditions(2)
            yMinBC=boundaryConditions(3)
            yMaxBC=boundaryConditions(4)
            CALL setBund2()
        elseif(iStreamModel==2)then !non-uniform grid,ISLBM
            CALL setBund()
        else
            stop
            write(*,*)'no such type LBMBC'
        endif
        !******************************************************************************************
        CALL cptMacr()  
        !Pbeta=(1.0d0-dexp(-5.0d0/Pramp*time/Tref))*Pbetatemp
        Pbeta=1.0d0*Pbetatemp
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y)  
         do x=1,xDim
            do y=1,yDim
               force(y,x,1:2)=0.d0
               !force_temp(y,x,1:2)=0.d0
            enddo
          enddo
        !$OMP END PARALLEL DO
        call date_and_time(VALUES=values0)
        ! package n Fish to one 
        do iFish=1,nFish
            
        if(iFish.eq.1)then
            do iND=1,nND(iFish)
               xyzful_all(iND,1:6)  =  xyzful(iFish,iND,1:6)
               velful_all(iND,1:6)  =  velful(iFish,iND,1:6)
               xyzfulIB_all(iND,1:6) = xyzfulIB(iFish,iND,1:6)
               extful_all(iND,1:6)  =  extful(iFish,iND,1:6)
            enddo
        elseif(iFish.ge.2)then
            do iND=1,nND(iFish)
               xyzful_all(iND+ sum(nND(1:iFish-1)),1:6)  =  xyzful(iFish,iND,1:6)
               velful_all(iND+ sum(nND(1:iFish-1)),1:6)  =  velful(iFish,iND,1:6)
               xyzfulIB_all(iND+ sum(nND(1:iFish-1)),1:6) = xyzfulIB(iFish,iND,1:6)
               extful_all(iND+ sum(nND(1:iFish-1)),1:6)  =  extful(iFish,iND,1:6)
            enddo
        endif
        
        enddo        

        !compute force exerted on fluids
        CALL cptForcS(yDim,xDim,nEL_all,nND_all,ele_all,dx,dy,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,  &
                xyzful_all,velful_all,xyzfulIB_all,Palpha,Pbeta,ntolLBM,dtolLBM,force,extful_all)
        
        !unpackage  
        do iFish=1,nFish
            if(iFish.eq.1)then
               do iND=1,nND(iFish)
                 xyzful(iFish,iND,1:6)   =  xyzful_all(iND,1:6) 
                 velful(iFish,iND,1:6)   =  velful_all(iND,1:6)
                 xyzfulIB(iFish,iND,1:6) =  xyzfulIB_all(iND,1:6)
                 extful(iFish,iND,1:6)   =  extful_all(iND,1:6) 
               enddo
            else
               do iND=1,nND(iFish)
                 xyzful(iFish,iND,1:6)   =  xyzful_all(iND+ sum(nND(1:iFish-1)),1:6) 
                 velful(iFish,iND,1:6)   =  velful_all(iND+ sum(nND(1:iFish-1)),1:6)
                 xyzfulIB(iFish,iND,1:6) =  xyzfulIB_all(iND+ sum(nND(1:iFish-1)),1:6)
                 extful(iFish,iND,1:6)   =  extful_all(iND+ sum(nND(1:iFish-1)),1:6) 
               enddo
            endif
        enddo 
        
        call date_and_time(VALUES=values1)
        write(*,*)'time for IBM step : ',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !compute repulsive force in short range. (W.-X. Huang et al. JCP 226 (2007) 2206-2228)
        call date_and_time(VALUES=values0)
        call cptForceR(dxmin,nFish,nND,nND_max,nEL,nEL_max,ele,xyzful,repful)
        call date_and_time(VALUES=values1)
        write(*,*)'time for Lubforce :',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        call date_and_time(VALUES=values0)
        CALL collide()
        call date_and_time(VALUES=values1)
        write(*,*)'time for collision:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)write(*,'(A)')' --------------------solid solver--------------------'
        !using substep for bodies
        subdeltat=deltat/numsubstep
        do isubstep=1,numsubstep
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iFish,iND)
        do iFish=1,nFish
        !======================================================
            if(iBodyModel(iFish)==1)then     ! update xyzful and velful for rigid body
                !prescribed motion
                !------------------------------------------------------
                !translational motion
                call setHeadMotion(iFish,time,deltat,isubstep,subdeltat)
                !XYZ(iFish,1:3)=XYZo(iFish,1:3)+XYZAmpl(iFish,1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+XYZPhi(iFish,1:3))
                !rotational motion
                AoA(iFish,1:3)=AoAo(iFish,1:3)+AoAAmpl(iFish,1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+AoAPhi(iFish,1:3))       
                call AoAtoTTT(AoA(iFish,1:3),TTTnxt(iFish,1:3,1:3))
                call get_angle_triad(TTT0(iFish,1:3,1:3),TTTnxt(iFish,1:3,1:3),AoAd(iFish,1),AoAd(iFish,2),AoAd(iFish,3))
                !prescribed displacement
                do  iND=1,nND(iFish)
                    xyzfulnxt(iFish,iND,1:3)=matmul(TTTnxt(iFish,1:3,1:3),xyzful00(iFish,iND,1:3))+XYZ(iFish,1:3)
                    xyzfulnxt(iFish,iND,4:6)=AoAd(iFish,1:3)
                enddo
                xyzful(iFish,1:nND(iFish),1:6)=xyzfulnxt(iFish,1:nND(iFish),1:6)
                !------------------------------------------------------
                !translational motion
                !UVW(iFish,1:3) =-2.0*pi*Freq*XYZAmpl(iFish,1:3)*dsin(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+XYZPhi(iFish,1:3))
                !rotational motion
                WWW1(iFish,1:3)=-2.0*pi*Freq*AoAAmpl(iFish,1:3)*dsin(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+AoAPhi(iFish,1:3))
                WWW2(iFish,1:3)=[ WWW1(iFish,1)*dcos(AoA(iFish,2))+WWW1(iFish,3),    &                       
                                  WWW1(iFish,1)*dsin(AoA(iFish,2))*dsin(AoA(iFish,3))+WWW1(iFish,2)*dcos(AoA(iFish,3)),   &
                                  WWW1(iFish,1)*dsin(AoA(iFish,2))*dcos(AoA(iFish,3))-WWW1(iFish,2)*dsin(AoA(iFish,3))    &
                                ]
                WWW3(iFish,1:3)=matmul(TTTnxt(iFish,1:3,1:3),WWW2(iFish,1:3))
                !prescribed velocity
                do  iND=1,nND(iFish)
                    velful(iFish,iND,1:3)=[   WWW3(iFish,2)*xyzful(iFish,iND,3)-WWW3(iFish,3)*xyzful(iFish,iND,2),    &
                                              WWW3(iFish,3)*xyzful(iFish,iND,1)-WWW3(iFish,1)*xyzful(iFish,iND,3),    &
                                              WWW3(iFish,1)*xyzful(iFish,iND,2)-WWW3(iFish,2)*xyzful(iFish,iND,1)     &
                                          ]  +UVW(iFish,1:3)
                    velful(iFish,iND,4:6)=WWW3(iFish,1:3)
                enddo
            elseif(iBodyModel(iFish)==2)then !flexible body   
                call date_and_time(VALUES=values0)
                !translational motion
                call setHeadMotion(iFish,time,deltat,isubstep,subdeltat)
                !XYZ(iFish,1:3)=XYZo(iFish,1:3)+XYZAmpl(iFish,1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+XYZPhi(iFish,1:3))
                !rotational motion
                AoA(iFish,1:3)=AoAo(iFish,1:3)+AoAAmpl(iFish,1:3)*dcos(2.0*pi*Freq*(time-deltat+isubstep*subdeltat)+AoAPhi(iFish,1:3))         
                call AoAtoTTT(AoA(iFish,1:3),TTTnxt(iFish,1:3,1:3))
                call get_angle_triad(TTT0(iFish,1:3,1:3),TTTnxt(iFish,1:3,1:3),AoAd(iFish,1),AoAd(iFish,2),AoAd(iFish,3))
                !prescribed displacement
                do  iND=1,nND(iFish)
                    xyzfulnxt(iFish,iND,1:3)=matmul(TTTnxt(iFish,1:3,1:3),xyzful00(iFish,iND,1:3))+XYZ(iFish,1:3)
                    xyzfulnxt(iFish,iND,4:6)=AoAd(iFish,1:3)
                enddo
                !-------------------------------------------------------
                !displacement condition
                vBC(iFish,1:nND(iFish),1:6) = xyzfulnxt(iFish,1:nND(iFish),1:6)-xyzful(iFish,1:nND(iFish),1:6)
                !loading vector
                lodful(iFish,1:nND(iFish),1:6) = extful(iFish,1:nND(iFish),1:6)+grav(iFish,1:nND(iFish),1:6) +repful(iFish,1:nND(iFish),1:6)
                !-----------------------------------------

                CALL solver(jBC(iFish,1:nND(iFish),1:6),vBC(iFish,1:nND(iFish),1:6),ele(iFish,1:nEL(iFish),1:5), &
                            nloc(iFish,1:nND(iFish)*6),nprof(iFish,1:nND(iFish)*6),nprof2(iFish,1:nND(iFish)*6), &
                            prop(iFish,1:nMT(iFish),1:10),mss(iFish,1:nND(iFish)*6), &
                            xyzful0(iFish,1:nND(iFish),1:6),xyzful(iFish,1:nND(iFish),1:6),dspful(iFish,1:nND(iFish),1:6), &
                            velful(iFish,1:nND(iFish),1:6),accful(iFish,1:nND(iFish),1:6),lodful(iFish,1:nND(iFish),1:6),  &
                            subdeltat,dampK,dampM,  &
                            triad_nn(iFish,1:3,1:3,1:nND(iFish)),triad_ee(iFish,1:3,1:3,1:nEL(iFish)),triad_e0(iFish,1:3,1:3,1:nEL(iFish)),  &
                            triad_n1(iFish,1:3,1:3,1:nEL(iFish)),triad_n2(iFish,1:3,1:3,1:nEL(iFish)),triad_n3(iFish,1:3,1:3,1:nEL(iFish)),  &
                            nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),nBD(iFish),maxstiff(iFish),Newmarkdelta,Newmarkalpha,dtolFEM,ntolFEM,&
                            nFish,iFish,FishInfo)
                !----------------------------------------------------------------------

                !call date_and_time(VALUES=values1)
                !write(*,*)'time for FEM:',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
            else
                stop 'no define body model'
            endif
        enddo !do iFish=1,nFish
        !$OMP END PARALLEL DO
        enddo !do isubstep=1,numsubstep
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)then
        do iFish=1,nFish
        write(*,'(A,I5.5)')' Fish number: ', int(FishInfo(iFish,1))
        write(*,'(A,I5.5,A,D20.10)')' iterFEM = ',int(FishInfo(iFish,2)),'    dmaxFEM = ',FishInfo(iFish,3)
        enddo
        endif
        write(*,'(A)')' ----------------------------------------------------'
        call date_and_time(VALUES=values1)
        write(*,*)'max time for FEM  :',(values1(6)*60.+values1(7)*1.+values1(8)*0.001)-(values0(6)*60.+values0(7)*1.+values0(8)*0.001)
        !******************************************************************************************
        !******************************************************************************************
        if(isRelease/=1)write(*,'(A)')' --------------------post process--------------------'
        if(isRelease/=1)then
        do iFish=1,nFish
        write(*,'(A,I5.5)')' Fish number: ',iFish
        write(*,'(A,2E20.10)')" forceDirect:", sum(extful(iFish,1:nND(iFish),1:2),1)/Fref
        write(*,'(A,2E20.10)')" velocityCnt:", sum(velful(iFish,1:nND(iFish),1:2)*mssful(iFish,1:nND(iFish),1:2),1)/sum(mssful(iFish,1:nND(iFish),1:2),1)/Uref
        write(*,'(A,2E20.10)')" coordinaCnt:", sum(xyzful(iFish,1:nND(iFish),1:2)*mssful(iFish,1:nND(iFish),1:2),1)/sum(mssful(iFish,1:nND(iFish),1:2),1)/Lref
        write(*,'(A,2E20.10)')" coordinaHed:", xyzful(iFish,1  ,1)/Lref,xyzful(iFish,1  ,2)/Lref
        write(*,'(A,2E20.10)')" coordinaTal:", xyzful(iFish,nND(iFish),1)/Lref,xyzful(iFish,nND(iFish),2)/Lref
        enddo
        endif                                                  
        !******************************************************************************************
        !****************************************************************************************** 
        if(DABS(time/Tref-timeOutTemp*NINT(time/Tref/timeOutTemp)) <= 0.5*dt/Tref)then
            CALL wrtCont()
        endif
          
        if(DABS(time/Tref-timeOutBody*NINT(time/Tref/timeOutBody)) <= 0.5*dt/Tref)then
            CALL writeBody3(nFish,xyzful/Lref,velful/Uref,accful/Aref,extful/Fref,repful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,0)
            if (dabs(Palpha)>eps) then
            CALL writeBody3(nFish,xyzful/Lref,velful/Uref,accful/Aref,extful/Fref,repful/Fref,ele,time/Tref,nND,nEL,nND_max,nEL_max,1)
            endif
        endif

        if(DABS(time/Tref-timeOutFlow*NINT(time/Tref/timeOutFlow)) <= 0.5*dt/Tref)then
            CALL writeFlow(1)
            if(isRelease/=1) then
                CALL writeFlow(0) 
            endif         
        endif
        
        if(DABS(time/Tref-timeOutInfo*NINT(time/Tref/timeOutInfo)) <= 0.5*dt/Tref)then
            CALL wrtInfo()
        endif
        !******************************************************************************************
        !******************************************************************************************
        write(*,'(A)')' ----------------------------------------------------'
        call date_and_time(VALUES=values_e)
        write(*,*)'time for one step : ',(values_e(6)*60.+values_e(7)*1.+values_e(8)*0.001)-(values_s(6)*60.+values_s(7)*1.+values_s(8)*0.001)
    enddo
    !pause 'stop'
    END PROGRAM main 
