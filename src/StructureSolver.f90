!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Finite element method for solid structure
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    SUBROUTINE structure_solver(  jBC,vBC,ele,nloc,nprof,nprof2,prop,mss,xyzful0,xyzful,dspful,velful,accful,lodExteful,deltat,dampK,dampM,  &
                        triad_nn,triad_ee,triad_e0,triad_n1,triad_n2,triad_n3,nND,nEL,nEQ,nMT,nBD,nSTF,NewmarkGamma,NewmarkBeta,dtol,iterMax)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nBD,nSTF
    integer:: jBC(nND,6),ele(nEL,5),nloc(nEQ),nprof(nEQ),nprof2(nEQ)
    real(8):: vBC(nND,6),xyzful0(nND,6),xyzful(nND,6),prop(nMT,10)
    real(8):: mss(nEQ),dspful(nND,6),velful(nND,6),accful(nND,6),lodExteful(nND,6),deltat
    real(8):: dampK,dampM

    real(8):: triad_nn(3,3,nND),triad_ee(3,3,nEL),triad_e0(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL),triad_n3(3,3,nEL)
!   ----------------------------------------------------------------------------------------------
    real(8):: du(nEQ),ddu(nEQ),lodEffe(nEQ),lodInte(nEQ),lodExte(nEQ),dsp(nEQ),vel(nEQ),acc(nEQ),dspO(nEQ),velO(nEQ),accO(nEQ)
    real(8):: stfEffe(nSTF),stfMatr(nSTF),stfElas(nSTF),stfGeom(nSTF)      
    real(8):: wk1(nEQ),wk2(nEQ) 
    real(8):: NewmarkGamma,NewmarkBeta,AlphaM,AlphaF,RhoInf
    real(8):: ak,a0,a1,a2,a3,a4,a5
    real(8):: beta0,beta,gamma,zi,z0
    real(8):: dsumd,dsumz,dtol,dnorm,geoFRM(nEL),geoPLT(1:9,nEL)
    integer:: i,j,iND,iEQ,iter,iterMax,iloc,ierror,maxramp,iModify
!   -----------------------------------------------------------------------------------------------
!   the generalized a method by Hua, Ru-Nan  not debug!!!!
    RhoInf=1.0d0
    AlphaM=0.0d0 !(2.0d0*RhoInf-1.0d0)/(RhoInf+1.0d0)
    AlphaF=0.0d0 !RhoInf/(RhoInf+1.0d0)

    !NewmarkGamma=0.5d0-AlphaM+AlphaF
    !NewmarkBeta=(1.0d0-AlphaM+AlphaF)**2/4.0d0


    a0 = (1.0d0-AlphaM)/(NewmarkBeta*deltat*deltat)   
    a2 = (1.0d0-AlphaM)/(NewmarkBeta*deltat)
    a3 = (1.0d0-AlphaM)/(NewmarkBeta*2.0d0) - 1.0d0

    ak = 1.0d0-AlphaF
    a1 = ak*NewmarkGamma/(NewmarkBeta*deltat)
    a4 = ak*NewmarkGamma/NewmarkBeta - 1.0d0
    a5 = ak*(NewmarkGamma/NewmarkBeta - 2.0d0)*0.5d0*deltat



    beta0 =1.0
    gamma =1.0
    maxramp =0

    iModify = 1

!   ***********************************************************************************************
    do  i=1,nND
        dsp((i-1)*6+1:(i-1)*6+6)=dspful(i,1:6)
        vel((i-1)*6+1:(i-1)*6+6)=velful(i,1:6)
        acc((i-1)*6+1:(i-1)*6+6)=accful(i,1:6)
        lodExte((i-1)*6+1:(i-1)*6+6)=lodExteful(i,1:6)        
    enddo

    dspO(1:nEQ) = dsp(1:nEQ)
    velO(1:nEQ) = vel(1:nEQ)
    accO(1:nEQ) = acc(1:nEQ)
!   ------------------------------------------------------------------------------
    iter=0
    dnorm=1.0
    do while(dnorm >= dtol .and. iter<= iterMax )
!   ------------------------------------------------------------------------------
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        forces from body stresses
        call body_stress_D( lodInte,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                            ele,prop,triad_n1,triad_n2,triad_n3,triad_ee,triad_e0,triad_nn, nND,nEL,nEQ,nMT,geoFRM,geoPLT &
                          )
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       -------------------------------------------------------------------
        call formstif_s(stfElas,ele,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                        prop,nprof,nloc,triad_e0,triad_ee,nND,nEL,nEQ,nMT,nSTF)
        call formgeom_s(stfGeom,ele,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                        prop,nprof,nloc,triad_e0,triad_ee,nND,nEL,nEQ,nMT,nSTF,geoFRM,geoPLT)
        stfMatr(1:nSTF)     = stfElas(1:nSTF)+ gamma*stfGeom(1:nSTF)
        stfEffe(1:nSTF)     = ak*stfMatr(1:nSTF)
        stfEffe(nloc(1:nEQ))= stfEffe(nloc(1:nEQ)) + a0*mss(1:nEQ) + a1*dampM*mss(1:nEQ)
        stfEffe(1:nSTF)     = stfEffe(1:nSTF)                     + a1*dampK*stfElas(1:nSTF)
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       -------------------------------------------------------------------
        lodEffe(1:nEQ)=+lodExte(1:nEQ)                                                            &            !external force
                       -lodInte(1:nEQ)                                                            &            !internal force
                       -(a0*(dsp(1:nEQ)-dspO(1:nEQ))-a2*vel(1:nEQ)-a3*acc(1:nEQ))*mss(1:nEQ)    &            !inertial force
                       -(a1*(dsp(1:nEQ)-dspO(1:nEQ))-a4*vel(1:nEQ)-a5*acc(1:nEQ))*mss(1:nEQ)*dampM            !mass-ratio dumping force

        if    (dabs(dampK) > 0.0) then  !
            wk1(1:nEQ)= a1*(dsp(1:nEQ)-dspO(1:nEQ)) -a4*vel(1:nEQ) -a5*acc(1:nEQ)
            call AxBCOL(stfElas,nSTF,wk1,wk2,nEQ,nBD,nprof,nprof2,nloc)
            lodEffe(1:nEQ) = lodEffe(1:nEQ) - dampK*wk2(1:nEQ)                                                !stiffness-ratio dumping force
        endif

        if    (dabs(AlphaF) > 0.0) then !
            wk1(1:nEQ)= (dsp(1:nEQ)-dspO(1:nEQ))
            call AxBCOL(stfElas+gamma*stfGeom,nSTF,wk1,wk2,nEQ,nBD,nprof,nprof2,nloc)
            lodEffe(1:nEQ) = lodEffe(1:nEQ) - AlphaF*wk2(1:nEQ)
        endif

!       -------------------------------------------------------------------
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       loading boundary conditions
        do  iND=1,nND
        do  j=1  ,6
            if(jBC(iND,j)>0)then
                iEQ=(iND-1)*6+j
                stfEffe(nloc(iEQ))=stfEffe(nloc(iEQ))*1.0d20
                if(iter==0)then
                    lodEffe(iEQ)=stfEffe(nloc(iEQ))*vBC(iND,j)
                else
                    lodEffe(iEQ)=stfEffe(nloc(iEQ))*0.0
                endif
            endif
        enddo
        enddo
                           
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       -------------------------------------------------------------------
!       solve equation
        call uduCOL_D(stfEffe,nSTF,nEQ,nBD,ierror,nprof,nloc)
        if    (ierror.eq.0) then
            write(*,*)'@@ ERROR: zero diagonal term !!!'
            return
        endif
        call bakCOL_D(stfEffe,nSTF,lodEffe,nEQ,nBD,du,ierror,nprof,nprof2,nloc)
!       -------------------------------------------------------------------
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if    (iter <= maxramp) then
            zi=2.0**(iter)
            z0=2.0**(maxramp)
            beta=zi/z0*beta0
        else
            beta=1.0d0*beta0
        endif

        ddu(1:nEQ)= beta*du(1:nEQ)
        dsp(1:nEQ)= dsp(1:nEQ) + ddu(1:nEQ)
        

        do  i=1,nND
            dspful(i,1:6)=dsp((i-1)*6+1:(i-1)*6+6)
        enddo

        xyzful(1:nND,1:6)=xyzful0(1:nND,1:6)+dspful(1:nND,1:6)

        call update_triad_D(ele,ddu,triad_nn,triad_n1,triad_n2,triad_n3,nND,nEL,nEQ)

        call make_triad_ee(ele,xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3),triad_ee,triad_n1,triad_n2,triad_n3,nND,nEL)

!       -------------------------------------------------------------------
!        test for convergence
        !if(iter==0)dsumd=dsqrt(sum((du(1:nEQ)*beta)**2))
        if(iter==0)dsumd=dabs(maxval((du(1:nEQ)*beta)**2))
        
        dsumz=dsqrt(sum(dsp(1:nEQ)**2))
        !if (dsumz < dtol/10.0) dsumz=dtol/10.0
        !dnorm=dsumd/dsumz

        !dnorm=dabs(maxval((du(1:nEQ)*beta)**2))
        if(iter==0)then
            dnorm=1.0
        else
            dnorm=dabs(maxval((du(1:nEQ)*beta)**2))/dsumd
        endif

        iter=iter+1
        if(iter>=100) write(*,*)'iter=',iter,'dnorm=',dnorm
        write(*,*)'iter=',iter,'dnorm=',dnorm 
    enddo

    write(*,'(A,I5,A,D20.10)')' iterFEM=',iter,' dmaxFEM   =',dnorm


    acc(1:nEQ)  = 1.0d0/(NewmarkBeta*deltat)*( (dsp(1:nEQ)-dspO(1:nEQ))/deltat -velO(1:nEQ) ) - (1.0d0/(NewmarkBeta*2.0d0) - 1.0d0)*accO(1:nEQ)

    vel(1:nEQ)  = velO(1:nEQ) + ((1.0d0 - NewmarkGamma)*accO(1:nEQ) +NewmarkGamma*acc(1:nEQ))*deltat

    do  i=1,nND
        velful(i,1:6)=vel((i-1)*6+1:(i-1)*6+6)
        accful(i,1:6)=acc((i-1)*6+1:(i-1)*6+6)        
    enddo
  
    return
    ENDSUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ structural DaTafile
    subroutine read_structural_datafile(jBC,ele,nloc,nprof,nprof2,xyzful0,prop,nND,nEL,nEQ,nMT,nBD,nSTF,idat)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nBD,nSTF,idat
    integer:: ele(nEL,5),jBC(nND,6),nloc(nND*6),nprof(nND*6),nprof2(nND*6)
    real(8):: xyzful0(nND,6),prop(nMT,10)
!   ---------------------------------------------------------------------------
    integer:: i,j,nm1,nm2,material,nbc,node,nmp,ii,ind,ibandh,iend,ibandv,ji1
    character*50 title,endin
!   -----------------------------------------------------------------------------------------------
!   READ  node
    read(idat,*)  nND
    do    i= 1, nND
        read(idat,*) node,xyzful0(node,1),xyzful0(node,2),xyzful0(node,3)
    enddo
    read(idat,'(1a50)') endin  
!   -----------------------------------------------------------------------------------------------
!   READ elem data
    read(idat,*) nEL
    do  i= 1, nEL
        read(idat,*) j,ele(j,1:5)
    enddo
    read(idat,'(1a50)') endin

!   -----------------------------------------------------------------------------------------------
!    READ  bcs  default is 0=free
    jBC(1:nND,1:6) = 0
    read(idat,*)  nbc
    do  i=1,nbc
        read(idat,*)node,jBC(node,1),jBC(node,2),jBC(node,3), &
                         jBC(node,4),jBC(node,5),jBC(node,6)   
       
    enddo
    read(idat,'(1a50)') endin  
!   -----------------------------------------------------------------------------------------------
!   READ element material properties
    read(idat,*) nMT
    do    i= 1, nMT
        read(idat,*) nmp,prop(nmp,1:8)
    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------

    nEQ=nND*6

    call max_band_width(ele,nprof,nND,nEL,nEQ,nBD)

!   nprof2
    do i=1,nEQ
       ibandh=1
       iend=i+nBD-1
       if (iend .gt. nEQ) iend=nEQ
       do j=i+1,iend
            ibandv=nprof(j)
            ji1=j-i+1
            if  (ibandv .ge. ji1) then
                ibandh = ji1
            endif
       enddo
       nprof2(i)=ibandh
    enddo

    nloc(1)=1
    do    i=1,nEQ-1
        nloc(i+1) = nloc(i) + nprof(i)
    enddo
    nSTF=nloc(nEQ)+nprof(nEQ)

    !write(*,'(3(A,1x,I8,2x))')'nND=',nND,'nEL=',nEL,'nEQ=',nEQ
    !write(*,'(3(A,1x,I8,2x))')'nMT=',nMT,'nBD=',nBD,'nSTF=',nSTF
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   MAX Band Width calculation
    subroutine max_band_width(ele,nprof,nND,nEL,nEQ,nBD)
    implicit none
    integer:: nBD,nEQ,nEL,nND
    integer:: ele(nEL,5),nprof(nEQ) 

    integer:: ipv(18)
    integer:: i,j,n,ihfbnd,idof,jdof,kdof,ieqn1,ieqn2,jband,ieq,ieq2

    nprof(1:nEQ)=0

    ihfbnd=0
    do  n=1,nEL
        idof=(ele(n,1)-1)*6
        jdof=(ele(n,2)-1)*6
        kdof=(ele(n,3)-1)*6
        do  i=1,6
            ipv(i   )=idof+i
            ipv(i+6 )=jdof+i
            ipv(i+12)=kdof+i
        enddo

        do  i=1,18
            ieqn1=ipv(i)          
            do  j=i,18
                ieqn2=ipv(j)
                ihfbnd = max0(ihfbnd,iabs(ieqn1-ieqn2))
                jband=abs(ieqn1-ieqn2)+1
                ieq=max(ieqn1,ieqn2)
                if  (jband .gt. nprof(ieq)) then
                    nprof(ieq)=jband
                endif
                ieq2=min(ieqn1,ieqn2)
!               if  (jband .gt. nprof2(ieq2)) then
!                   nprof2(ieq2)=jband
!               endif
            enddo
        enddo

    enddo
    nBD=ihfbnd+1

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE write_coord(xyzful,velful,accful,ele,time,nND,nEL)
    implicit none
    integer:: nND,nEL
    real(8):: xyzful(nND,6),velful(nND,6),accful(nND,6)
    integer:: ele(nEL,5)
    real(8):: time
!   -------------------------------------------------------

    integer:: i,j,ireduc
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName
    !==================================================================================================        
    integer::    nv
    integer,parameter:: namLen=40,idfile=100,numVar=6
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=['x','y','z','u','v','w'] 
    !==================================================================================================

    write(fileName,'(F7.3)') time
    fileName = adjustr(fileName)
    DO  I=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    END DO

    OPEN(idfile,FILE='./DatBody/Body'//trim(filename)//'.plt',form='unformatted')

!   I. The header section.
!   =============================================        
    write(idfile) "#!TDV101"    
    write(idfile) 1
    call dumpstring(title,idfile)
    write(idfile) numVar
    do  nv=1,numVar
        call dumpstring(varname(nv),idfile)
    enddo
!   ---------------------------------------------
!   zone head for ZONE1
!   ---------------------------------------------
    write(idfile) ZONEMARKER              
    call dumpstring(zonename,idfile)
    write(idfile) -1,2,1,0,0,nND,nEL,0,0,0,0

!   =============================================
    write(idfile) EOHMARKER
!   =============================================
!   II. Data section
!   =============================================
    write(idfile) Zonemarker
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    do    i=1,nND
        write(idfile)   real(xyzful(i,1:3)),real(velful(i,1:3))                            
    enddo
    do    i=1,nEL
        write(idfile) ele(i,1),ele(i,2),ele(i,3)
    enddo
    close(idfile)
    ENDSUBROUTINE write_coord

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    set angle of initial triads
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE init_triad_D(ele,xord0,yord0,zord0,triad_nn,triad_n1,triad_n2,triad_n3,triad_ee,triad_e0,nND,nEL)
    implicit none
    integer:: nND,nEL         
    integer:: ele(nEL,5)
    real(8):: xord0(nND),yord0(nND),zord0(nND)
    real(8):: triad_nn(3,3,nND),triad_ee(3,3,nEL),triad_e0(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL),triad_n3(3,3,nEL)
!
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,dx,dy,dz,xl0
    real(8):: xll1,xmm1,xnn1, xll2,xmm2,xnn2, xll3,xmm3,xnn3
    real(8):: axy,ayz,azx,area,dd

    integer:: i,j,n,i1,j1,k1,nELt

!   For each node, set the triad  to global system
    do  n=1,nND

        do i=1,3
        do j=1,3
            triad_nn(i,j,n)=0.0
        enddo
        enddo
        triad_nn(1,1,n)=1.0
        triad_nn(2,2,n)=1.0
        triad_nn(3,3,n)=1.0

    enddo
!
!   For each element, calculate the current orientation triad
    do  n=1,nEL

        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
!
        if    (nELt == 2) then
!           frame
            dx = xord0(j1) - xord0(i1)
            dy = yord0(j1) - yord0(i1)
            dz = zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll1=dx/xl0
            xmm1=dy/xl0
            xnn1=dz/xl0
            dd=dsqrt(xll1*xll1+xmm1*xmm1)
            triad_n1(1,1,n)=xll1
            triad_n1(2,1,n)=xmm1
            triad_n1(3,1,n)=xnn1

            if    (dd .lt. 0.001) then
                triad_n1(1,2,n)=0.0
                triad_n1(2,2,n)=1.0
                triad_n1(3,2,n)=0.0
                triad_n1(1,3,n)=-xnn1
                triad_n1(2,3,n)=0.00
                triad_n1(3,3,n)=0.00
            else
                triad_n1(1,2,n)=-xmm1/dd
                triad_n1(2,2,n)=+xll1/dd
                triad_n1(3,2,n)=0.0

                triad_n1(1,3,n)=-xll1*xnn1/dd
                triad_n1(2,3,n)=-xmm1*xnn1/dd                                    
                triad_n1(3,3,n)= dd 
            endif
!            all element triads have same initial orientation
            triad_n2(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_ee(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_e0(1:3,1:3,n)=triad_n1(1:3,1:3,n)
        elseif (nELt == 3) then
!           plate
            x1=xord0(i1)
            x2=xord0(j1)
            x3=xord0(k1)
            y1=yord0(i1)
            y2=yord0(j1)
            y3=yord0(k1)
            z1=zord0(i1)
            z2=zord0(j1)
            z3=zord0(k1)

            dx = xord0(j1) - xord0(i1)
            dy = yord0(j1) - yord0(i1)
            dz = zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll1=dx/xl0
            xmm1=dy/xl0
            xnn1=dz/xl0
            !tengent1 vector (one edge) (Hua)
            triad_n1(1,1,n)=xll1
            triad_n1(2,1,n)=xmm1
            triad_n1(3,1,n)=xnn1
!
!            determine vector area
            axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
            ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
            azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
            area=dsqrt( axy*axy + ayz*ayz + azx*azx)
            xll3=ayz/area
            xmm3=azx/area
            xnn3=axy/area
            !normal vector    (Hua)
            triad_n1(1,3,n)=xll3
            triad_n1(2,3,n)=xmm3
            triad_n1(3,3,n)=xnn3

            !cross product            
            xll2=xmm3*xnn1 - xnn3*xmm1
            xmm2=xnn3*xll1 - xll3*xnn1
            xnn2=xll3*xmm1 - xmm3*xll1
            !tengent2 vector (Hua)
            triad_n1(1,2,n)=xll2
            triad_n1(2,2,n)=xmm2
            triad_n1(3,2,n)=xnn2
!
!           all element triads have same initial orientation

            triad_n2(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_n3(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_ee(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_e0(1:3,1:3,n)=triad_n1(1:3,1:3,n)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
    ENDSUBROUTINE init_triad_D


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    update angle of  triads
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE update_triad_D(ele,ddut,triad_nn,triad_n1,triad_n2,triad_n3,nND,nEL,nEQ)
    implicit none
    integer:: nND,nEL,nEQ         
    integer:: ele(nEL,5)
    real(8):: ddut(nEQ)
    real(8):: triad_nn(3,3,nND),triad_n1(3,3,nEL),triad_n2(3,3,nEL),triad_n3(3,3,nEL)
   
    real(8):: ddutful(nND,6),rr(3,3),sumd(3)
    real(8):: dtx1,dty1,dtz1
    integer:: i,j,k,m,n,i1,j1,k1,node,jdof,ireduc,nELt

    do    i=1,nND*6
        node=(i+5)/6
        jdof=i-(node-1)*6
        ireduc=i
        ddutful(node,jdof)=ddut(ireduc)
    enddo
!
!   For each node, set the triad_nn = [R]triad_nn
    do    n=1,nND
        dtx1=ddutful(n,4)
        dty1=ddutful(n,5)
        dtz1=ddutful(n,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_nn(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_nn(1:3,1:3,n))
    enddo

!   For each element, calculate the current orientation triad
    do  n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
!       n1 node
        dtx1=ddutful(i1,4)
        dty1=ddutful(i1,5)
        dtz1=ddutful(i1,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_n1(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_n1(1:3,1:3,n))
!       n2 node
        dtx1=ddutful(j1,4)
        dty1=ddutful(j1,5)
        dtz1=ddutful(j1,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_n2(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_n2(1:3,1:3,n))
!
        if (nELt == 2) then
!           frame
        elseif (nELt == 3) then
!            plate
!           n3 node
            dtx1=ddutful(k1,4)
            dty1=ddutful(k1,5)
            dtz1=ddutful(k1,6)
            call finite_rot(dtx1,dty1,dtz1,rr)
            triad_n3(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_n3(1:3,1:3,n))
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
    ENDSUBROUTINE update_triad_D

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    get orientation of element
!    copyright@ RuNanHua 
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE make_triad_ee(ele,xord,yord,zord,triad_ee,triad_n1,triad_n2,triad_n3,nND,nEL)
    implicit none
    integer:: nND,nEL
    integer:: ele(nEL,5)
    real(8):: xord(nND), yord(nND), zord(nND)
!
    real(8):: triad_aa(3,3)
    real(8):: triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL),triad_n3(3,3,nEL)
    real(8):: rr(3,3),tx,ty,tz
    real(8):: triad_11(3,3),triad_22(3,3)
    real(8):: sumd(3)
!
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,dx,dy,dz,xl0
    real(8):: xll1,xmm1,xnn1, xll2,xmm2,xnn2, xll3,xmm3,xnn3
    real(8):: axy,ayz,azx,area
    
    real(8):: xll,xmm,xnn,dd,r2e1,r3e1
    integer:: i,j,k,m,n,i1,j1,k1,nELt
!
!

!   For each element, calculate the current orientation triad
    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
!
        if    (nELt == 2) then
!            frame
            dx = xord(j1) - xord(i1)
            dy = yord(j1) - yord(i1)
            dz = zord(j1) - zord(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl0
            xmm=dy/xl0
            xnn=dz/xl0
            dd =dsqrt(xll*xll+xmm*xmm)
            do    i=1,3
            do    j=1,3
                triad_ee(i,j,n)=0.0
            enddo
            enddo
            triad_ee(1,1,n)=xll
            triad_ee(2,1,n)=xmm
            triad_ee(3,1,n)=xnn
!
!            get angle between two triads
            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_n1(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo
            call get_angle_triad( triad_11,triad_22,tx,ty,tz)
!
!           rotate n1 to intermediate
            tx=tx/2.0
            ty=ty/2.0
            tz=tz/2.0
            call finite_rot(tx,ty,tz,rr)
            triad_aa(1:3,1:3)=matmul(rr(1:3,1:3),triad_n1(1:3,1:3,n))
!
!
!           vectors e2 e3
            r2e1 = 0.0
            r3e1 = 0.0
            do    k=1,3
                r2e1 = r2e1 + triad_aa(k,2)*triad_ee(k,1,n)
                r3e1 = r3e1 + triad_aa(k,3)*triad_ee(k,1,n)
            enddo
            do    j=1,3
                triad_ee(j,2,n)=triad_aa(j,2) - r2e1*(triad_aa(j,1)+triad_ee(j,1,n))/2.0
                triad_ee(j,3,n)=triad_aa(j,3) - r3e1*(triad_aa(j,1)+triad_ee(j,1,n))/2.0
            enddo
!
        elseif (nELt == 3) then
!           plate
            x1=xord(i1)
            x2=xord(j1)
            x3=xord(k1)
            y1=yord(i1)
            y2=yord(j1)
            y3=yord(k1)
            z1=zord(i1)
            z2=zord(j1)
            z3=zord(k1)
            dx = xord(j1) - xord(i1)
            dy = yord(j1) - yord(i1)
            dz = zord(j1) - zord(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll1=dx/xl0
            xmm1=dy/xl0
            xnn1=dz/xl0
            triad_ee(1,1,n)=xll1
            triad_ee(2,1,n)=xmm1
            triad_ee(3,1,n)=xnn1
!
!            determine vector area
            axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
            ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
            azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
            area=dsqrt( axy*axy + ayz*ayz + azx*azx)
            xll3=ayz/area
            xmm3=azx/area
            xnn3=axy/area

            triad_ee(1,3,n)=xll3
            triad_ee(2,3,n)=xmm3
            triad_ee(3,3,n)=xnn3
!
            xll2=xmm3*xnn1 - xnn3*xmm1
            xmm2=xnn3*xll1 - xll3*xnn1
            xnn2=xll3*xmm1 - xmm3*xll1

            triad_ee(1,2,n)=xll2
            triad_ee(2,2,n)=xmm2
            triad_ee(3,2,n)=xnn2
!
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo
!   end of loop over elements
    return
    ENDSUBROUTINE make_triad_ee

!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE finite_rot(t1,t2,t3,rr)
    implicit none
    real(8):: rr(3,3),rr1(3,3),rr2(3,3),rr3(3,3)
    real(8):: t1,t2,t3,tt,ss,cc,c1,c2

    integer:: i,j
!
    tt=dsqrt( t1**2 + t2**2 + t3**2 )
    ss=dsin(tt)
    cc=dcos(tt)
!
    rr1(1,1)=1.0
    rr1(1,2)=0.0
    rr1(1,3)=0.0
    rr1(2,1)=0.0
    rr1(2,2)=1.0
    rr1(2,3)=0.0
    rr1(3,1)=0.0
    rr1(3,2)=0.0
    rr1(3,3)=1.0
!
    rr2(1,1)=0.0
    rr2(1,2)=-t3
    rr2(1,3)= t2
    rr2(2,1)= t3
    rr2(2,2)=0.0
    rr2(2,3)=-t1
    rr2(3,1)=-t2
    rr2(3,2)= t1
    rr2(3,3)=0
!
    rr3(1,1)=-t3*t3-t2*t2
    rr3(1,2)= t2*t1
    rr3(1,3)= t3*t1
    rr3(2,1)= t1*t2
    rr3(2,2)=-t3*t3-t1*t1
    rr3(2,3)= t3*t2
    rr3(3,1)= t1*t3
    rr3(3,2)= t2*t3
    rr3(3,3)=-t2*t2-t1*t1
!
    do    i=1,3
    do    j=1,3
        if    (tt .lt. 1.0e-10) then
            c1=1.0
!!!!        c2=1.0
            c2=0.5
        else
            c1 = ss/tt
            c2 = (1.0-cc)/tt**2
        endif
        rr(i,j) = rr1(i,j) + rr2(i,j)*c1 + rr3(i,j)*c2
    enddo
    enddo
!
    return
    ENDSUBROUTINE finite_rot


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    GET ANGLE of between TRIADs
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE get_angle_triad(triad_n1,triad_n2,tx,ty,tz)
    implicit none
    real(8):: triad_n1(3,3),triad_n2(3,3)
    real(8):: rr(3,3)
    real(8):: tx,ty,tz, dtx,dty,dtz,c1,tt,sint
    integer:: i,j,k
!
!   get angle between two triads
    do    i=1,3
    do    j=1,3
        rr(i,j)=0.0
        do    k=1,3
            rr(i,j)=rr(i,j) + triad_n2(i,k)*triad_n1(j,k)
        enddo
    enddo
    enddo

    dtx = (rr(3,2)-rr(2,3))/2.0
    dty = (rr(1,3)-rr(3,1))/2.0
    dtz = (rr(2,1)-rr(1,2))/2.0

    c1=1.0
    sint = dsqrt(dtx*dtx+dty*dty+dtz*dtz)

    if (sint .gt. 1.0) sint=1.0
    tt = dasin(sint)
    if ( sint .lt. 1.0e-6) then
         c1=1.0
    else
         c1 = tt/sint
    endif

    tx=c1*dtx
    ty=c1*dty
    tz=c1*dtz

    return
    ENDSUBROUTINE get_angle_triad
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE global_to_local(triad,tx,ty,tz,tx2,ty2,tz2)
    implicit none
    real(8):: triad(3,3)
    real(8):: tx,ty,tz,tx2,ty2,tz2
    tx2 = triad(1,1)*tx+triad(2,1)*ty+triad(3,1)*tz
    ty2 = triad(1,2)*tx+triad(2,2)*ty+triad(3,2)*tz
    tz2 = triad(1,3)*tx+triad(2,3)*ty+triad(3,3)*tz
    return
    ENDSUBROUTINE global_to_local

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Nodal loads due to body stresses
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE body_stress_D(    gforce,xord0,yord0,zord0,xord,yord,zord, &
                                ele,prop,triad_n1,triad_n2,triad_n3,triad_ee,triad_e0,triad_nn, &
                                nND,nEL,nEQ,nMT,geoFRM,geoPLT &
                            )
    implicit none
    integer:: nMT,nEL,nND,nEQ
    integer:: ele(nEL,5)
    real(8):: xord0(nND), yord0(nND), zord0(nND)
    real(8):: xord(nND), yord(nND), zord(nND)
!
    real(8):: prop(nMT,10),geoFRM(nEL),geoPLT(1:9,nEL)
 
    real(8):: strain(18),stress(18), force(18),forceb(18)
    real(8):: gforce(nEQ)
    real(8):: ek9(9,9)
    real(8):: ekb12(12,12)
    real(8):: ekb(18,18)
!
    real(8):: triad_nn(3,3,nND)
    real(8):: triad_ee(3,3,nEL),triad_e0(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL),triad_n3(3,3,nEL)

    real(8):: triad_00(3,3),triad_11(3,3),triad_22(3,3)
    real(8):: rr(3,3)
    real(8):: ub(18),dl,temp(6)
!
    real(8):: dx0,dy0,dz0,du,dv,dw
    real(8):: dx,dy,dz,xl0
    real(8):: tx1,tx2,tx3,ty1,ty2,ty3,tz1,tz2,tz3,tx,ty,tz
    real(8):: xyz012(3),xyz013(3),xyzb012(3),xyzb013(3)
    real(8):: xyz12(3),xyz13(3),xyzb12(3),xyzb13(3)
    real(8):: xyz11(3),xyzb11(3)
    real(8):: xb01,xb02,xb03,yb01,yb02,yb03
    real(8):: xb1,xb2,xb3,yb1,yb2,yb3
    real(8):: uub(6)


    real(8):: e0,g0,a0,b0,r0,zix0,ziy0,ziz0,xl,fxx,t0,pl0,zip0,zia0,zib0,alpha,beta
    integer:: i,j,k,n,i1,j1,k1,mat,nELt

!    pi=4.0*datan(1.0d0)

    gforce(1:nEQ)=0.0
  
    !rewind(igeoFRM)
    !rewind(igeoPLT)


!   For each element, calculate the nodal forces
    do    n=1,nEL

        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)        
        nELt= ele(n,4)
        mat = ele(n,5)
        if    ( nELt == 2) then
!           frame
            e0  =prop(mat,1)
            g0  =prop(mat,2)
            a0  =prop(mat,3)
            r0  =prop(mat,4)
            b0  =prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)

            dx0 = xord0(j1) - xord0(i1)
            dy0 = yord0(j1) - yord0(i1)
            dz0 = zord0(j1) - zord0(i1)
            xl0 = dsqrt(dx0*dx0+dy0*dy0+dz0*dz0)
!
!           orientation
            du = (xord(j1)-xord0(j1))-(xord(i1)-xord0(i1))
            dv = (yord(j1)-yord0(j1))-(yord(i1)-yord0(i1))
            dw = (zord(j1)-zord0(j1))-(zord(i1)-zord0(i1))

            dx = dx0 + du
            dy = dy0 + dv
            dz = dz0 + dw
            xl =dsqrt(dx*dx+dy*dy+dz*dz)
!
            dl = ( (2*dx0+du)*du +(2*dy0+dv)*dv +(2*dz0+dw)*dw )/ (xl+xl0)
!           get twisting angles
            do    i=1,3
            do    j=1,3
                triad_00(i,j)=triad_ee(i,j,n)
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n1(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)
!
            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

!            non-zero ty1 tz1 u2 tx2 ty2 tz2
            ub(1)=0.0
            ub(2)=0.0
            ub(3)=0.0
            ub(4)=tx1
            ub(5)=ty1
            ub(6)=tz1
!
            ub(7)=dl
            ub(8)=0.0
            ub(9)=0.0
            ub(10)=tx2
            ub(11)=ty2
            ub(12)=tz2
!
!
!           get current stiffness in local coords. use L0

            call elmstfFRM_D(xl0,zix0,ziy0,ziz0,a0,e0,g0,ekb12,nELt )
!
!           compute axial force
            fxx=dl*e0*a0/xl0
!           save local force for geo stiff
            !write(igeoFRM,'(D25.15)') fxx
            geoFRM(n)=fxx
!           nodal forces in local coords
!           {F}=[k]{u}
            forceb(1:12) =matmul(ekb12(1:12,1:12),ub(1:12))
            forceb(13:18)=0.0

!           transform to global
            do  i=1,3
            do  j=1,3
                rr(i,j)=triad_ee(i,j,n)
            enddo
            enddo
            do  i=1,18
                force(i)=0.0
            enddo
            do    i=1,3
            do    k=1,3
                force(0+i) = force(0+i) + rr(i,k)*forceb(0+k)
                force(3+i) = force(3+i) + rr(i,k)*forceb(3+k)
                force(6+i) = force(6+i) + rr(i,k)*forceb(6+k)
                force(9+i) = force(9+i) + rr(i,k)*forceb(9+k)
            enddo
            enddo

            call assembFOR(nEQ,nND,gforce,force,i1,j1,k1)
!
        elseif (nELt == 3) then
!           plate
            e0=prop(mat,1)
            g0=prop(mat,2)
            t0=prop(mat,3)
            r0=prop(mat,4)
            pl0=prop(mat,5)
            zip0=prop(mat,6)
            zia0=prop(mat,7)
            zib0=prop(mat,8)
!
!           original lengths of triangle
            xyz012(1) = xord0(j1)-xord0(i1)
            xyz012(2) = yord0(j1)-yord0(i1)
            xyz012(3) = zord0(j1)-zord0(i1)
            xyz013(1) = xord0(k1)-xord0(i1)
            xyz013(2) = yord0(k1)-yord0(i1)
            xyz013(3) = zord0(k1)-zord0(i1)
!
!           use element triad to rotate to local x=[Rt]x
            do    i=1,3
            do    j=1,3
                rr(i,j)=triad_e0(i,j,n)
            enddo
            enddo
            do    i=1,3
                xyzb012(i)=0.0
                xyzb013(i)=0.0
                do    k=1,3
                    xyzb012(i) = xyzb012(i) + rr(k,i)*xyz012(k)
                    xyzb013(i) = xyzb013(i) + rr(k,i)*xyz013(k)
                enddo
            enddo
            xb01=0.0
            yb01=0.0
            xb02=xyzb012(1)
            yb02=xyzb012(2)
            xb03=xyzb013(1)
            yb03=xyzb013(2)

!           current lengths of triangle
            xyz11(1) = xord(i1)-xord0(i1)
            xyz11(2) = yord(i1)-yord0(i1)
            xyz11(3) = zord(i1)-zord0(i1)
            xyz12(1) = xord(j1)-xord0(i1)
            xyz12(2) = yord(j1)-yord0(i1)
            xyz12(3) = zord(j1)-zord0(i1)
            xyz13(1) = xord(k1)-xord0(i1)
            xyz13(2) = yord(k1)-yord0(i1)
            xyz13(3) = zord(k1)-zord0(i1)
!
!           use element triad to rotate to local x=[Rt]x
            do    i=1,3
            do    j=1,3
                rr(i,j)=triad_ee(i,j,n)
            enddo
            enddo
            do    i=1,3
                xyzb11(i)=0.0
                xyzb12(i)=0.0
                xyzb13(i)=0.0
                do    k=1,3
                    xyzb11(i) = xyzb11(i) + rr(k,i)*xyz11(k)
                    xyzb12(i) = xyzb12(i) + rr(k,i)*xyz12(k)
                    xyzb13(i) = xyzb13(i) + rr(k,i)*xyz13(k)
                enddo
            enddo

            xb1=xyzb11(1)
            yb1=xyzb11(2)
            xb2=xyzb12(1)
            yb2=xyzb12(2)
            xb3=xyzb13(1)
            yb3=xyzb13(2)

!           get angles ref to orig from triads
            do    i=1,3
            do    j=1,3
!2001!!!        triad_00(i,j)=triad_e0(i,j,n)
                triad_00(i,j)=triad_ee(i,j,n)
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n1(i,j,n)
            enddo
            enddo

            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)

            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo

            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n3(i,j,n)
            enddo
            enddo

            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx3,ty3,tz3)
!
!           local DoF
            ub(1)=xb1-xb01
            ub(2)=yb1-yb01
            ub(3)=0.0
            ub(4)=tx1
            ub(5)=ty1
            ub(6)=tz1
            ub(7)=xb2-xb02
            ub(8)=yb2-yb02
            ub(9)=0.0
            ub(10)=tx2
            ub(11)=ty2
            ub(12)=tz2
            ub(13)=xb3-xb03
            ub(14)=yb3-yb03
            ub(15)=0.0
            ub(16)=tx3
            ub(17)=ty3
            ub(18)=tz3

!           get forces
            do    k=1,18
                force(k)=0.0
            enddo
!           Calculate LOCAL stiffness, use orig coords
            do    i=1,18
            do    j=1,18
                ekb(i,j)=0.0
            enddo
            enddo


!            membrane
            alpha=zia0
            beta =zib0
            call elmstfMRT_D(e0,g0,t0,pl0,alpha,beta,xb01,xb02,xb03,yb01,yb02,yb03,ekb,ek9)

!           flexure
            call elmstfDKT_D(e0,g0,t0,zip0, xb01,xb02,xb03,yb01,yb02,yb03,ekb,ek9)

!            nodal forces in local coords
!           {F}=[k]{u}
            forceb(1:18) =matmul(ekb(1:18,1:18),ub(1:18))
!           save local force for geo stiff           
            !write(igeoPLT,'(9(D25.15,1x))') forceb(1),forceb(2),forceb(3),forceb(7),forceb(8),forceb(9),forceb(13),forceb(14),forceb(15)
            geoPLT(1:9,n)=[forceb(1),forceb(2),forceb(3),forceb(7),forceb(8),forceb(9),forceb(13),forceb(14),forceb(15)]
!           transform to global
            rr(1:3,1:3)=triad_ee(1:3,1:3,n)

            force(1:18)=0

            do    i=1,3
            do    k=1,3
                force(0+i) = force(0+i) + rr(i,k)*forceb(0+k)
                force(3+i) = force(3+i) + rr(i,k)*forceb(3+k)
                force(6+i) = force(6+i) + rr(i,k)*forceb(6+k)
                force(9+i) = force(9+i) + rr(i,k)*forceb(9+k)
                force(12+i)= force(12+i)+ rr(i,k)*forceb(12+k)
                force(15+i)= force(15+i)+ rr(i,k)*forceb(15+k)
            enddo
            enddo

            call assembFOR(nEQ,nND,gforce,force,i1,j1,k1)

        else
            write(*,*)'not this nELt:',nELt
            stop
        endif

    enddo
    
    return
    ENDSUBROUTINE body_stress_D

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    FORM MASS matrix: only lumped
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine formmass_D(ele,xord,yord,zord,prop,mss,nND,nEL,nEQ,nMT,alphaf,alpham,alphap)
    implicit none                  
    integer:: nND,nEL,nEQ,nMT
    integer:: ele(nEL,5)
    real(8):: prop(nMT,10)
!
    real(8):: xord(nND), yord(nND), zord(nND)
    real(8):: mss(nEQ), em(18,18), emb(18,18)
    real(8):: em12(12,12)
 
    integer:: i,i1,j1,k1,ii,jj,mat,nELt,n
    real(8):: a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl,xll,xmm,xnn,tt0,rh0,dxb,dyb,dzb
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,axy,ayz,azx,area,xb1,yb1,zb1,xb2,yb2,zb2,xb3,yb3,zb3

    !rotational inertial factors
    real(8):: alphaf,alpham,alphap

!   zero array before assembling

    mss(1:nEQ)=0.0
!
!   form the element form matrix, and assemble
    do  i=1,nEL
        i1  = ele(i,1)
        j1  = ele(i,2)
        k1  = ele(i,3)        
        nELt= ele(i,4)
        mat = ele(i,5)
        if (nELt == 2) then
            a0=prop(mat,3)
            r0=prop(mat,4)
            b0=prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl
            call elmmasFRM_D(r0,a0,xl,zix0,em12,alphaf)
            call trans3d_D(xll,xmm,xnn,em12,b0)
            em(1:18,1:18)=0.0
            em(1:12,1:12)=em12(1:12,1:12)
            call assembLUM(nND,nEQ,mss,em,i1,j1,k1)

        elseif(nELt == 3) then
            tt0=prop(mat,3)
            rh0=prop(mat,4)
            x1=xord(i1)
            x2=xord(j1)
            x3=xord(k1)
            y1=yord(i1)
            y2=yord(j1)
            y3=yord(k1)
            z1=zord(i1)
            z2=zord(j1)
            z3=zord(k1)
!
!           determine vector area  A=0.5*v1Xv2
            axy =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
            ayz =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.
            azx =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.
            area=dsqrt( axy*axy + ayz*ayz + azx*azx)
            xll=ayz/area
            xmm=azx/area
            xnn=axy/area
!
!           transform element coords to local X-Y
            xb1=x1
            yb1=y1
            zb1=z1
            call rotvec_D(xll,xmm,xnn,(x2-x1),(y2-y1),(z2-z1),dxb,dyb,dzb)
            xb2=x1+dxb
            yb2=y1+dyb
            zb2=z1+dzb
            call rotvec_D(xll,xmm,xnn,(x3-x1),(y3-y1),(z3-z1),dxb,dyb,dzb)
            xb3=x1+dxb
            yb3=y1+dyb
            zb3=z1+dzb

!           Calculate the mass matrix and assemble
            emb(1:18,1:18)=0.0        
            call elmmasCST_D(rh0,area,tt0,emb)
            call elmmasMRT_D(rh0,area,tt0,emb,alpham)
            call elmmasPLT_D(rh0,area,tt0,emb,alphap)

            call rotate_D(xll,xmm,xnn,emb,em)
            call assembLUM(nND,nEQ,mss,em,i1,j1,k1)

        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo    
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    FORM STIFfness matrix  [K]
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine formstif_s(stf,ele,xord0,yord0,zord0,xord,yord,zord,prop,nprof,nloc,triad_e0,triad_ee,nND,nEL,nEQ,nMT,nSTF)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nSTF  
    integer:: ele(nEL,5),nprof(nEQ), nloc(nEQ) 
    real(8):: prop(nMT,10)
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)
    real(8):: stf(nSTF)
    real(8):: ek(18,18), ekb(18,18),ek12(12,12),ek9(9,9)
      
    real(8):: triad_e0(3,3,nEL),triad_ee(3,3,nEL),rr(3,3)
    
    real(8):: xyz12(3),xyz13(3),xyzb12(3),xyzb13(3)
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3
    real(8):: xb1,xb2,xb3,yb1,yb2,yb3

    integer:: i,j,k,ipress,n,i1,j1,k1,mat,nELt,ii,jj
    real(8):: e0,g0,a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl0,xl,xll,xmm,xnn,xl9,area,t0,zip0,zia0,zib0,pl0,alpha,beta
    real(8):: temp

!   initialize [K]  to zero
    stf(1:nSTF)=0.0

!   form each element matrix, and assemble
    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)       
        nELt= ele(n,4)
        mat = ele(n,5)
!
        if    (nELt == 2) then
!           frame
!           material props not change
            e0=prop(mat,1)
            g0=prop(mat,2)
            a0=prop(mat,3)
            r0=prop(mat,4)
            b0=prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)
            dx= xord0(j1) - xord0(i1)
            dy= yord0(j1) - yord0(i1)
            dz= zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
!
!           orientation
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl

!           Calculate the stiffness matrix and assemble
            xl9= xl0
            call elmstfFRM_D(xl9,zix0,ziy0,ziz0,a0,e0,g0,ek12, nELt)
!           use trans3d to transform from local to global stifness
            call trans3d_D(xll,xmm,xnn,ek12,b0)
!
!           expand to [18x18]

            ek(1:18,1:18)=0.0

            do    j=1,12
            do    k=1,12
                ek(j,k)=ek12(j,k)
            enddo
            enddo

            call assembCOL(nSTF,nND,nEQ,stf,ek,i1,j1,k1,nloc)
!
        elseif (nELt == 3) then
            e0=prop(mat,1)
            g0=prop(mat,2)
            t0=prop(mat,3)
            r0=prop(mat,4)
            pl0=prop(mat,5)
            zip0=prop(mat,6)
            zia0=prop(mat,7)
            zib0=prop(mat,8)

            x1=xord0(i1)
            x2=xord0(j1)
            x3=xord0(k1)
            y1=yord0(i1)
            y2=yord0(j1)
            y3=yord0(k1)
            z1=zord0(i1)
            z2=zord0(j1)
            z3=zord0(k1)
!
!           use element triad to rotate coords to local x=Rtx
            rr(1:3,1:3)=triad_e0(1:3,1:3,n)

            xyz12(1)=x2-x1
            xyz12(2)=y2-y1
            xyz12(3)=z2-z1
            xyz13(1)=x3-x1
            xyz13(2)=y3-y1
            xyz13(3)=z3-z1
            do    i=1,3
                xyzb12(i)=0.0
                xyzb13(i)=0.0
                do    k=1,3
                    xyzb12(i)=xyzb12(i)+rr(k,i)*xyz12(k)
                    xyzb13(i)=xyzb13(i)+rr(k,i)*xyz13(k)
                enddo
            enddo

            xb1=0.0
            yb1=0.0
            xb2=xyzb12(1)
            yb2=xyzb12(2)
            xb3=xyzb13(1)
            yb3=xyzb13(2)

!           Calculate the LOCAL stiffness matrix
            ekb(1:18,1:18)=0.0
            alpha=zia0
            beta =zib0
            call elmstfMRT_D(e0,g0,t0,pl0,alpha,beta,xb1,xb2,xb3,yb1,yb2,yb3,ekb,ek9)
            call elmstfDKT_D(e0,g0,t0,zip0,xb1,xb2,xb3,yb1,yb2,yb3,ekb,ek9)
            rr(1:3,1:3)=triad_ee(1:3,1:3,n)
            call rot_mat_LG(rr,ekb,ek)
            call assembCOL(nSTF,nND,nEQ,stf,ek,i1,j1,k1,nloc)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    FORM GEOMetric stiffness matrices  [KGx],[KGy],[KGxy]
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine formgeom_s(geo,ele,xord0,yord0,zord0,xord,yord,zord,prop,nprof,nloc,triad_e0,triad_ee,nND,nEL,nEQ,nMT,nSTF,geoFRM,geoPLT)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nSTF
    integer:: ele(nEL,5),nprof(nEQ), nloc(nEQ)
    real(8):: prop(nMT,10),geoFRM(nEL),geoPLT(1:9,nEL)
    real(8):: eg12(12,12),eg(18,18)
!
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)

    real(8):: geo(nSTF),egm(18,18),eg6(6,6)
    real(8):: x0b1,x0b2,x0b3,y0b1,y0b2,y0b3
    real(8):: xb1,xb2,xb3,yb1,yb2,yb3,zb1,zb2,zb3
    real(8):: x1,x2,x3,y1,y2,y3,x10,y10,z10
!
    real(8):: triad_e0(3,3,nEL),triad_ee(3,3,nEL),rr(3,3)
    real(8):: triad_11(3,3),triad_22(3,3)
    real(8):: xd1,xd2,xd3,yd1,yd2,yd3,zd1,zd2,zd3
    real(8):: xyz12(3),xyz13(3),xyzb12(3),xyzb13(3)
    real(8):: xyz11(3),xyzb11(3)

    integer:: n,i,j,k,i1,j1,k1,mat,nELt,ii,jj
    real(8):: t0,e0,g0,a0,b0,dx,dy,dz,xl0,xl,xll,xmm,xnn,sxx,s,xl9,z1,z2,z3,x01,y01,z01,forceb(9)

!   initialize [G] to zero
    geo(1:nSTF)=0.0

    !rewind(igeoFRM)
    !rewind(igeoPLT)

    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)       
        nELt= ele(n,4)
        mat = ele(n,5)

        t0=prop(mat,3)
!
        if (nELt == 2) then
!           constit relation not change
            e0=prop(mat,1)
            g0=prop(mat,2)
            a0=prop(mat,3)
            b0=prop(mat,5)
            dx= xord0(j1) - xord0(i1)
            dy= yord0(j1) - yord0(i1)
            dz= zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
!           orientation
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl
!
!           Calculate the stiffness matrix and assemble
            !read(igeoFRM,*) sxx
            sxx=geoFRM(n)
            s=sxx
            xl9= xl0
            call elmgeomFRM_D(xl9,eg12,s)
!           use trans3d to transform from local to global stifness
            call trans3d_D(xll,xmm,xnn,eg12,b0)

!           expand to [18x18]
            eg(1:18,1:18)=0.0
            eg(1:12,1:12)=eg12(1:12,1:12)
            call assembCOL(nSTF,nND,nEQ,geo,eg ,i1,j1,k1,nloc)

        elseif (nELt == 3) then
!           use element triad to rotate coords to local x=Rtx

            rr(1:3,1:3)=triad_e0(1:3,1:3,n)

            xd1=xord0(i1)
            xd2=xord0(j1)
            xd3=xord0(k1)
            yd1=yord0(i1)
            yd2=yord0(j1)
            yd3=yord0(k1)
            zd1=zord0(i1)
            zd2=zord0(j1)
            zd3=zord0(k1)
!
            xyz12(1)=xd2-xd1
            xyz12(2)=yd2-yd1
            xyz12(3)=zd2-zd1
            xyz13(1)=xd3-xd1
            xyz13(2)=yd3-yd1
            xyz13(3)=zd3-zd1
            do  i=1,3
                xyzb12(i)=0.0
                xyzb13(i)=0.0
                do  k=1,3
                    xyzb12(i)=xyzb12(i)+rr(k,i)*xyz12(k)
                    xyzb13(i)=xyzb13(i)+rr(k,i)*xyz13(k)
                enddo
            enddo

            x0b1=0.0
            y0b1=0.0
            x0b2=xyzb12(1)
            y0b2=xyzb12(2)
            x0b3=xyzb13(1)
            y0b3=xyzb13(2)
!
!           membrane form
            x01=xord0(i1)
            y01=yord0(i1)
            x1=xord(i1)
            x2=xord(j1)
            x3=xord(k1)
            y1=yord(i1)
            y2=yord(j1)
            y3=yord(k1)
!
            x10=xord0(i1)
            y10=yord0(i1)
            z10=zord0(i1)
            x1=xord(i1)
            x2=xord(j1)
            x3=xord(k1)
            y1=yord(i1)
            y2=yord(j1)
            y3=yord(k1)
            z1=zord(i1)
            z2=zord(j1)
            z3=zord(k1)
!

            rr(1:3,1:3)=triad_ee(1:3,1:3,n)

            xyz11(1)=x1-x10
            xyz11(2)=y1-y10
            xyz11(3)=z1-z10
            xyz12(1)=x2-x10
            xyz12(2)=y2-y10
            xyz12(3)=z2-z10
            xyz13(1)=x3-x10
            xyz13(2)=y3-y10
            xyz13(3)=z3-z10
            do  i=1,3
                xyzb11(i)=0.0
                xyzb12(i)=0.0
                xyzb13(i)=0.0
                do  k=1,3
                    xyzb11(i)=xyzb11(i)+rr(k,i)*xyz11(k)
                    xyzb12(i)=xyzb12(i)+rr(k,i)*xyz12(k)
                    xyzb13(i)=xyzb13(i)+rr(k,i)*xyz13(k)
                enddo
            enddo

            xb1=xyzb11(1)
            yb1=xyzb11(2)
            zb1=xyzb11(3)
            xb2=xyzb12(1)
            yb2=xyzb12(2)
            zb2=xyzb12(3)
            xb3=xyzb13(1)
            yb3=xyzb13(2)
            zb3=xyzb13(3)
!

            triad_11(1:3,1:3) = triad_e0(1:3,1:3,n)
            triad_22(1:3,1:3) = triad_ee(1:3,1:3,n)


            !read(igeoPLT,*) forceb(1:9)
            forceb(1:9)=geoPLT(1:9,n)
            call elmgeomCST_s(x1,x2,x3,y1,y2,y3,x0b1,x0b2,x0b3,y0b1,y0b2,y0b3, &
                                    xb1,xb2,xb3,yb1,yb2,yb3,zb1,zb2,zb3, &
                                    t0,eg6,egm,forceb,triad_11,triad_22,x10,y10)

            call assembCOL(nSTF,nND,nEQ,geo,egm,i1,j1,k1,nloc)

        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent MASs matrix for the FRaMe
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmmasFRM_D(rho,area,length,zix,em,alphaf)
    implicit none
    real(8):: rho, area, length,zix
    real(8):: em(12,12)
    real(8):: alphaf,roal

    em(1:12,1:12) = 0.0
    roal = rho*area*length/2.0
    em(1,1)     = roal
    em(2,2)     = roal
    em(3,3)     = roal
    em(4,4)     = roal*zix/area
    em(5,5)     = roal*length*length*alphaf/48
    em(6,6)     = roal*length*length*alphaf/48
    em(7,7)     = em(1,1)
    em(8,8)     = em(2,2)
    em(9,9)     = em(3,3)
    em(10,10)   = em(4,4)
    em(11,11)   = em(5,5)
    em(12,12)   = em(6,6)
    return
    end
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ELeMent MASs matrix for Constant Strain Triangle
    subroutine elmmasCST_D(rho, area, th, em)
    implicit none              
    real(8):: rho,area,th,em(18,18)
    real(8):: emb(9,9)   
    real(8):: roat
    integer:: inew(9),i,j,ii,jj
!
    emb(1:9,1:9) = 0.0

    roat = rho*area*th/3.0
    emb(1,1) = roat
    emb(2,2) = roat
    emb(3,3) = roat
    emb(4,4) = roat
    emb(5,5) = roat
    emb(6,6) = roat
    emb(7,7) = roat
    emb(8,8) = roat
    emb(9,9) = roat

    inew(1)=1
    inew(2)=2
    inew(3)=3
    inew(4)=7
    inew(5)=8
    inew(6)=9
    inew(7)=13
    inew(8)=14
    inew(9)=15

    do    i=1,9
        ii=inew(i)
        do    j=1,9
            jj=inew(j)
            em(ii,jj) = emb(i,j)
        enddo
    enddo
!
    return
    end
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent MASs matrix for Moment Rotation Triangle
    subroutine elmmasMRT_D(rho,area,th,em ,alpham)
    implicit none  
    real(8):: alpham,rho,area,th,em(18,18),emb(9,9)
    real(8):: roat,rad,alpha
    integer:: inew(9),i,j,ii,jj

    emb(1:9,1:9) = 0.0


    roat    = rho*area*th/3.0
    rad     = dsqrt(area/3.0)
    alpha   = area/6.0
    rad     = dsqrt(area/3.12)
    alpha   = rad*rad/40
    alpha   = 20*alpha*alpham

    emb(1,1) = roat
    emb(2,2) = roat
    emb(3,3) = roat*alpha
    emb(4,4) = roat
    emb(5,5) = roat
    emb(6,6) = roat*alpha
    emb(7,7) = roat
    emb(8,8) = roat
    emb(9,9) = roat*alpha

    inew(1)=1
    inew(2)=2
    inew(3)=6
    inew(4)=7
    inew(5)=8
    inew(6)=12
    inew(7)=13
    inew(8)=14
    inew(9)=18

    do    i=1,9
        ii=inew(i)
        do    j=1,9
            jj=inew(j)
            em(ii,jj) = emb(i,j)
        enddo
    enddo
!
    return
    end
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ELeMent MASs matrix for PLaTe: triangle
    subroutine elmmasPLT_D(rho, area, th, em ,alphap)
    implicit none
    real(8):: rho, area, th,alphap
    real(8):: em(18,18), emb(9,9)
    integer:: inew(9)
    real(8):: roat,alpha,rad
    integer:: i,ii,j,jj
!
    emb(1:9,1:9) = 0.0

    roat    = rho*area*th/3.0
    alpha   = 1.0e-6
    alpha   = 100*area/12.0
    rad     = dsqrt(area/3.12)
    alpha   = rad*rad/40
    alpha   = alpha*alphap

    emb(1,1) = roat
    emb(2,2) = roat*alpha
    emb(3,3) = roat*alpha
    emb(4,4) = roat
    emb(5,5) = roat*alpha
    emb(6,6) = roat*alpha
    emb(7,7) = roat
    emb(8,8) = roat*alpha
    emb(9,9) = roat*alpha

!   assign to full element matrix
    inew(1)=3
    inew(2)=4
    inew(3)=5
    inew(4)=9
    inew(5)=10
    inew(6)=11
    inew(7)=15
    inew(8)=16
    inew(9)=17
    do    i=1,9
        ii=inew(i)
        do    j=1,9
            jj=inew(j)
            em(ii,jj) = emb(i,j)
        enddo
    enddo
!
!
    return
    end
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent STiFfness for Discrete Kirchhoff Triangle: CMP pp332
!    copyright@ RuNanHua 
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmstfDKT_D(e0,g0,tt0,zip0,x1,x2,x3,y1,y2,y3,ek,ekb)
    implicit none
    real(8):: e0,g0,tt0,zip0,x1,x2,x3,y1,y2,y3
    real(8):: ekb(9,9),ek(18,18)

    real(8):: gg(10,9),b(3),c(3),als(3),px(3,3),d(3,3),dd(9,9),qq(9,9),pp(3,3),pt(2,3),rs(2,3),q(3)
    real(8):: sum,inew(9)
    integer:: kod(2,9)
    data kod/1,1,2,3,3,2,4,4,5,6,6,5,7,7,8,9,9,8/
    data pp /12.0,4.0,4.0,4.0,2.0,1.0,4.0,1.0,2.0/

    integer:: i,j,k,k1,k2,ii,jj,l
    real(8):: dt,det,znu,r1,r2,r3,p1,p2,p3
!
    znu = e0/(2.0*g0) - 1.0
    dt=(e0*tt0**3)/12.0/(1-znu*znu)
    dt=(e0*zip0   )    /(1-znu*znu)

    d(1,1) = dt
    d(1,2) = dt*znu
    d(1,3) = 0.0
    d(2,1) = dt*znu
    d(2,2) = dt
    d(2,3) = 0.0
    d(3,1) = 0.0
    d(3,2) = 0.0
    d(3,3) = dt*(1.0-znu)/2.0
!
    b(1) = y2-y3
    b(2) = y3-y1
    b(3) = y1-y2
    c(1) = x3-x2
    c(2) = x1-x3
    c(3) = x2-x1
    det  = 24.0*(b(1)*c(2)-b(2)*c(1))
!
    do  i=1,3
    do  j=1,3
         px(i,j) = pp(i,j)/det
    enddo
    enddo
!
    do  i=1,3
    do  j=1,3
        do  k1=1,3
            ii=(i-1)*3+k1
            do  k2=1,3
                jj=(j-1)*3+k2
                dd(ii,jj)=d(i,j)*px(k1,k2)
            enddo
        enddo
    enddo
    enddo
!
    do  i=1,3
        als(i)  = b(i)*b(i)+c(i)*c(i)
        pt(1,i) = 6.0*c(i)/als(i)
        pt(2,i) = 6.0*b(i)/als(i)
        rs(1,i) = 3.0*c(i)*c(i)/als(i)
        rs(2,i) = 3.0*b(i)*b(i)/als(i)
        q(i)    = 3.0*b(i)*c(i)/als(i)
    enddo
!
    do  i=1,10
    do  j=1,9
        gg(i,j) = 0.0
    enddo
    enddo
!
    do  i=1,2
        ii=(i-1)*5
        p1=pt(i,1)
        p2=pt(i,2)
        p3=pt(i,3)
        r1=rs(i,1)
        r2=rs(i,2)
        r3=rs(i,3)
        gg(ii+1,kod(i,1)) =  p3
        gg(ii+2,kod(i,1)) = -p2
        gg(ii+3,kod(i,1)) = -p3
        gg(ii+4,kod(i,1)) =  p2-p3
        gg(ii+5,kod(i,1)) =  p2
        gg(ii+1,kod(i,2)) = -q(3)
        gg(ii+2,kod(i,2)) = -q(2)
        gg(ii+3,kod(i,2)) =  q(3)
        gg(ii+4,kod(i,2)) =  q(2)+q(3)
        gg(ii+5,kod(i,2)) =  q(2)
        gg(ii+1,kod(i,3)) = -1.0-r3
        gg(ii+2,kod(i,3)) = -1.0-r2
        gg(ii+3,kod(i,3)) =  r3
        gg(ii+4,kod(i,3)) =  r2+r3
        gg(ii+5,kod(i,3)) =  r2
        gg(ii+1,kod(i,4)) = -p3
        gg(ii+3,kod(i,4)) =  p3
        gg(ii+4,kod(i,4)) =  p1+p3
        gg(ii+1,kod(i,5)) = -q(3)
        gg(ii+3,kod(i,5)) =  q(3)
        gg(ii+4,kod(i,5)) =  q(3)-q(1)
        gg(ii+1,kod(i,6)) =  1.0-r3
        gg(ii+3,kod(i,6)) =  r3
        gg(ii+4,kod(i,6)) =  r3-r1
        gg(ii+2,kod(i,7)) =  p2
        gg(ii+4,kod(i,7)) = -p1-p2
        gg(ii+5,kod(i,7)) = -p2
        gg(ii+2,kod(i,8)) = -q(2)
        gg(ii+4,kod(i,8)) =  q(2)-q(1)
        gg(ii+5,kod(i,8)) =  q(2)
        gg(ii+2,kod(i,9)) =  1.0-r2
        gg(ii+4,kod(i,9)) =  r2-r1
        gg(ii+5,kod(i,9)) =  r2
    enddo
!
    do  i=1,9
        qq(1,i) =     b(2)*gg(1,i) +     b(3)*gg(2,i)
        qq(2,i) = 2.0*b(2)*gg(3,i) +     b(3)*gg(4,i)
        qq(3,i) =     b(2)*gg(4,i) + 2.0*b(3)*gg(5,i)
        qq(4,i) =    -c(2)*gg(6,i) -     c(3)*gg(7,i)
        qq(5,i) =-2.0*c(2)*gg(8,i) -     c(3)*gg(9,i)
        qq(6,i) =-    c(2)*gg(9,i) - 2.0*c(3)*gg(10,i)
        qq(7,i) =     c(2)*gg(1,i) +     c(3)*gg(2 ,i) -    b(2)*gg(6,i) -     b(3)*gg(7,i)
        qq(8,i) = 2.0*c(2)*gg(3,i) +     c(3)*gg(4 ,i) -2.0*b(2)*gg(8,i) -     b(3)*gg(9,i)
        qq(9,i) =     c(2)*gg(4,i) + 2.0*c(3)*gg(5 ,i) -    b(2)*gg(9,i) - 2.0*b(3)*gg(10,i)
    enddo
!
    do  i=1,9
    do  j=1,9
        gg(i,j) = 0.0
        do  k=1,9
            gg(i,j) = gg(i,j) + dd(i,k)*qq(k,j)
        enddo
    enddo
    enddo
!
    do  L=1,9
    do  j=L,9
        sum = 0.0
        do  k=1,9
            sum = sum  + qq(k,L)*gg(k,j)
        enddo
        ekb(L,j)=sum
        ekb(j,L)=sum
    enddo
    enddo

!    assign to full element matrix[18X18]
    inew(1)=3
    inew(2)=4
    inew(3)=5
    inew(4)=9
    inew(5)=10
    inew(6)=11
    inew(7)=15
    inew(8)=16
    inew(9)=17
    do    i=1,9
        ii=inew(i)
        do    j=1,9
            jj=inew(j)
            ek(ii,jj) = ekb(i,j)
        enddo
    enddo

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent STiFfness for Constant Strain Triangle
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmstfCST_D(e0,g0,tt0,pl0,x1,x2,x3,y1,y2,y3,ek,ek9)
    implicit none
    integer:: inew(9)
    real(8):: e0,g0,tt0,pl0,x1,x2,x3,y1,y2,y3
    real(8):: ek(18,18),ekb(6,6),b(3,6),c(3,3),etemp(3,6),ek9(9,9)

    integer:: i,j,k,ii,jj
    real(8):: znu0,xkap,a,a2,att0
!
    ekb(1:6,1:6)=0.0
    etemp(1:3,1:6)=0.0
!   tag
    znu0=(e0/g0/2.0)-1.0
!   plane strain
    if (pl0 .lt. 0) xkap=3.0-4.0*znu0
!   plane stress
    if (pl0 .gt. 0) xkap=(3.0-znu0)/(1.0+znu0)
!
!   {s}=[C]{e}
    c(1,1)=g0/(xkap-1.0)*(xkap+1.0)
    c(1,2)=g0/(xkap-1.0)*(3.0-xkap)
    c(1,3)=0.0
    c(2,1)=c(1,2)
    c(2,2)=c(1,1)
    c(3,2)=0.0
    c(3,1)=0.0
    c(3,2)=0.0
    c(3,3)=g0
!
!   {e}=[B]{u}
    a=((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.
    a2=a*2.0
    b(1,1)=(y2-y3)/a2
    b(1,2)=0.0
    b(1,3)=(y3-y1)/a2
    b(1,4)=0.0
    b(1,5)=(y1-y2)/a2
    b(1,6)=0.0
    b(2,1)=0.0
    b(2,2)=(x3-x2)/a2
    b(2,3)=0.0
    b(2,4)=(x1-x3)/a2
    b(2,5)=0.0
    b(2,6)=(x2-x1)/a2
    b(3,1)=b(2,2)
    b(3,2)=b(1,1)
    b(3,3)=b(2,4)
    b(3,4)=b(1,3)
    b(3,5)=b(2,6)
    b(3,6)=b(1,5)
!
!   [k] = At[B][C][B]
    att0=a*tt0
    do  i=1,3
    do  j=1,6
    do  k=1,3
        etemp(i,j)=etemp(i,j)+c(i,k)*b(k,j)
    enddo
    enddo
    enddo

    do  i=1,6
    do  j=1,6
    do  k=1,3
        ekb(i,j)=ekb(i,j)+b(k,i)*etemp(k,j)*att0
    enddo
    enddo
    enddo
!
!   assign [ekb] to [18x18]
    inew(1)=1
    inew(2)=2
    inew(3)=7
    inew(4)=8
    inew(5)=13
    inew(6)=14
    do  i=1,6
        ii=inew(i)
        do  j=1,6
            jj=inew(j)
            ek(ii,jj) = ekb(i,j)
        enddo
    enddo

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent STiFfness for FRaMe
!    calculates the element stiffness matrices.
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmstfFRM_D(length, ix, iy, iz,area, emod, gmod, ek ,nELt)
    implicit none
    real(8):: area, length, ix, iy, iz, emod, gmod
    real(8):: ek(12,12)
    integer:: nELt

    integer:: i,j
    real(8):: emlen,emlen2,emlen3
!
!   initialize all ek elements to zero
    ek(1:12,1:12)=0.0
!
!   STIFFNESS matrix in local coordinates
!
    emlen  = emod/length
    emlen2 = emlen/length
    emlen3 = emlen2/length
    if  (nELt == 2) then
        ek(1,1)   =   area*emlen
        ek(2,2)   =   12.0*emlen3*iz
        ek(3,3)   =   12.0*emlen3*iy
        ek(4,4)   =   gmod*ix/length
        ek(5,5)   =   4.0*emlen*iy
        ek(6,6)   =   4.0*emlen*iz
!
        ek(2,6)   =   6.0*emlen2*iz
        ek(3,5)   =  -6.0*emlen2*iy
    else
        write(*,*)'not this nELt:',nELt
        stop
    endif
!
    ek(7,7)   =   ek(1,1)
    ek(8,8)   =   ek(2,2)
    ek(9,9)   =   ek(3,3)
    ek(10,10) =   ek(4,4)
    ek(11,11) =   ek(5,5)
    ek(12,12) =   ek(6,6)
!
    ek(1,7)   =   -ek(1,1)
    ek(2,8)   =   -ek(2,2)
    ek(2,12)  =    ek(2,6)
    ek(3,9)   =   -ek(3,3)
    ek(3,11)  =    ek(3,5)
    ek(4,10)  =   -ek(4,4)
    ek(5,9)   =   -ek(3,5)
    ek(5,11)  =    ek(5,5)/2.0
    ek(6,8)   =   -ek(2,6)
    ek(6,12)  =    ek(6,6)/2.0
!
    ek(8,12)  =   -ek(2,6)
    ek(9,11)  =   -ek(3,5)
!
!
!   impose the symmetry
    do  i= 1, 12
    do  j= i, 12
        ek(j,i) = ek(i,j)
    enddo
    enddo
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent STiFFness for Membrane with Rotation Triangle
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmstfMRT_D(e0,g0,tt0,pl0,alpha0,beta0,x1,x2,x3,y1,y2,y3,ek,ekb)
    implicit none
    real(8):: e0,g0,tt0,pl0,alpha0,beta0,x1,x2,x3,y1,y2,y3
    real(8):: ek(18,18),ekb(9,9)
    real(8):: xt(3), yt(3),dmt(3,3)
    real(8):: alpha,f,beta,fbeta,xkap,znu0
    integer:: inew(9)
    integer:: lst(12)=[1,2,4,5,7,8,3,6,9,10,11,12]
    integer:: i,j,m,ii,jj
!
    znu0=(e0/g0/2.0)-1.0
!   plane stress
    if (pl0 .ge. 0) xkap=(3.0-znu0)/(1.0+znu0)
!   plane strain
    if (pl0 .lt. 0) xkap=3.0-4.0*znu0
!
!   {s}=[C]{e}
    dmt(1,1)=g0*tt0/(xkap-1.0)*(xkap+1.0)
    dmt(1,2)=g0*tt0/(xkap-1.0)*(3.0-xkap)
    dmt(1,3)=0.0
    dmt(2,1)=dmt(1,2)
    dmt(2,2)=dmt(1,1)
    dmt(3,2)=0.0
    dmt(3,1)=0.0
    dmt(3,2)=0.0
    dmt(3,3)=g0*tt0

    xt(1)=x1
    xt(2)=x2
    xt(3)=x3
    yt(1)=y1
    yt(2)=y2
    yt(3)=y3
    f=1.0
    alpha=alpha0
    beta =beta0
!

    ekb(1:9,1:9)=0.0

    m=9
    call sm3mb_D(xt,yt,dmt,alpha,f,lst,ekb)
!
    if    (beta .gt. 0.0) then
        fbeta=f*beta
        call sm3mh_D(xt,yt,dmt,fbeta,lst,ekb)
    endif

!   assign [ekb] to [18x18]
    inew(1)=1
    inew(2)=2
    inew(3)=6
    inew(4)=7
    inew(5)=8
    inew(6)=12
    inew(7)=13
    inew(8)=14
    inew(9)=18
    do    i=1,9
        ii=inew(i)
        do    j=1,9
            jj=inew(j)
            ek(ii,jj) = ekb(i,j)
        enddo
    enddo
!
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent GEOMetric stiffness for Constant Strain Triangle
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmgeomCST_s(x1,x2,x3,y1,y2,y3,x0b1,x0b2,x0b3,y0b1,y0b2,y0b3,xb1,xb2,xb3,yb1,yb2,yb3,zb1,zb2,zb3, &
                            tt0,ekg6,ekg18,forceb,triad_11,triad_22,x10,y10)
    implicit none    
    real(8):: x1,x2,x3,y1,y2,y3,x0b1,x0b2,x0b3,y0b1,y0b2,y0b3
    real(8):: xb1,xb2,xb3,yb1,yb2,yb3,zb1,zb2,zb3
    real(8):: tt0,ekg6(6,6),ekg18(18,18),forceb(9)
    real(8):: triad_11(3,3),triad_22(3,3),x10,y10

    real(8):: vv02t(3,9),vv02(9,3),ssf(9,3),ssx(3,9)
    real(8):: sk1(9,9),sk2(9,9),sk3(9,9),skb18(18,18)
    real(8):: wk1(9,9),wk2(9,9)
    real(8):: fb1x,fb1y,fb1z,fb2x,fb2y,fb2z,fb3x,fb3y,fb3z,a2,b1,b2,b3,c1,c2,c3
    integer:: inew(9)
    integer:: i,j,k,ii,jj,kk

    fb1x=forceb(1)
    fb1y=forceb(2)
    fb1z=forceb(3)
    fb2x=forceb(4)
    fb2y=forceb(5)
    fb2z=forceb(6)
    fb3x=forceb(7)
    fb3y=forceb(8)
    fb3z=forceb(9)

!   3-D form
    a2=(y0b1-y0b2)*(x0b3-x0b2) + (x0b2-x0b1)*(y0b3-y0b2)
    b1=(y0b2-y0b3)/a2
    b2=(y0b3-y0b1)/a2
    b3=(y0b1-y0b2)/a2
    c1=(x0b3-x0b2)/a2
    c2=(x0b1-x0b3)/a2
    c3=(x0b2-x0b1)/a2
!
!   3-D form

    vv02t(1:3,1:9)=0.0
    ssf(1:9,1:3)=0.0
    ssx(1:3,1:9)=0.0


    vv02t(1,3)= 2*c1
    vv02t(2,3)=-2*b1
    vv02t(3,1)=-  c1
    vv02t(3,2)=   b1
    vv02t(1,6)= 2*c2
    vv02t(2,6)=-2*b2
    vv02t(3,4)=-  c2
    vv02t(3,5)=   b2
    vv02t(1,9)= 2*c3
    vv02t(2,9)=-2*b3
    vv02t(3,7)=-  c3
    vv02t(3,8)=   b3
    do    i=1,3
    do    j=1,9
        vv02(j,i) = vv02t(i,j)
    enddo
    enddo
!
!   fb1z=0
!   fb2z=0
!   fb3z=0
    ssf(1,2) = -fb1z
    ssf(1,3) =  fb1y
    ssf(2,1) =  fb1z
    ssf(2,3) = -fb1x
    ssf(3,1) = -fb1y
    ssf(3,2) =  fb1x
    ssf(4,2) = -fb2z
    ssf(4,3) =  fb2y
    ssf(5,1) =  fb2z
    ssf(5,3) = -fb2x
    ssf(6,1) = -fb2y
    ssf(6,2) =  fb2x
    ssf(7,2) = -fb3z
    ssf(7,3) =  fb3y
    ssf(8,1) =  fb3z
    ssf(8,3) = -fb3x
    ssf(9,1) = -fb3y
    ssf(9,2) =  fb3x
!
    ssx(1,2) = -zb1
    ssx(1,3) =  yb1
    ssx(2,1) =  zb1
    ssx(2,3) = -xb1
    ssx(3,1) = -yb1
    ssx(3,2) =  xb1
    ssx(1,5) = -zb2
    ssx(1,6) =  yb2
    ssx(2,4) =  zb2
    ssx(2,6) = -xb2
    ssx(3,4) = -yb2
    ssx(3,5) =  xb2
    ssx(1,8) = -zb3
    ssx(1,9) =  yb3
    ssx(2,7) =  zb3
    ssx(2,9) = -xb3
    ssx(3,7) = -yb3
    ssx(3,8) =  xb3
!
!   [S(F)][Vt]
    do  i=1,9
    do  j=1,9
        sk1(i,j) = 0.0
        do  k=1,3
            sk1(i,j) = sk1(i,j) + ssf(i,k)*vv02t(k,j)
        enddo
        sk1(i,j) = -sk1(i,j)/2.0
    enddo
    enddo

    do  i=1,9
    do  j=1,9
        sk2(i,j) =  sk1(j,i)
    enddo
    enddo
!
!   [S(x)][S(F)]
    do  i=1,3
    do  j=1,3
        wk1(i,j) = 0.0
        do  k=1,9
            wk1(i,j) = wk1(i,j) + ssx(i,k)*ssf(k,j)
        enddo
    enddo
    enddo
!
    do  jj=1,3
    do  i=1,3
    do  j=1,3
        wk2(i,j) = 0.0
        do  k=1,3
            ii=(jj-1)*3+i
            kk=(jj-1)*3+k
            wk2(i,j) = wk2(i,j) + ssx(i,kk)*ssf(kk,j)
        enddo
    enddo
    enddo
    enddo

!   [S(x)][S(F)][Vt]
    do  i=1,3
    do  j=1,9
        wk2(i,j) = 0.0
        do  k=1,3
            wk2(i,j) = wk2(i,j) + wk1(i,k)*vv02t(k,j)
        enddo
    enddo
    enddo
!   [V][S(x)][S(F)][Vt]
    do  i=1,9
    do  j=1,9
        sk3(i,j) = 0.0
        do  k=1,3
            sk3(i,j) = sk3(i,j) + vv02(i,k)*wk2(k,j)
        enddo
        sk3(i,j) = sk3(i,j)/4.0
    enddo
    enddo

!   symmetrize sk3
    do  i=1,9
    do  j=1,9
        wk2(i,j) = (sk3(i,j) + sk3(j,i))/2.0
    enddo
    enddo

    do  i=1,9
    do  j=1,9
        sk3(i,j) = wk2(i,j)
    enddo
    enddo

!   assign [skb] to [18x18]

    skb18(1:18,1:18) = 0.0


    inew(1)=1
    inew(2)=2
    inew(3)=3
    inew(4)=7
    inew(5)=8
    inew(6)=9
    inew(7)=13
    inew(8)=14
    inew(9)=15
    do    i=1,9
        ii=inew(i)
        do    j=1,9
            jj=inew(j)
            skb18(ii,jj) = sk1(i,j) + sk2(i,j) + sk3(i,j)
        enddo
    enddo

    call rot_mat_LG(triad_22,skb18,ekg18)

    return
    end
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent GEOMetric stiffness matrix for a FRaMe
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmgeomFRM_D(length,eg,s)
    implicit none
    real(8):: length,eg(12,12),s
    real(8):: emlenz,alpha,beta
    integer:: i,j
!
!   initialize all eg elements to zero
    eg(1:12,1:12)=0.0
!   Stiffness matrix in local coordinates
!** if (s .gt. 200.) s=200.
!   emlenz  =   s/(30.0*length)
!***alpha   =   s*1.0e-0
!   beta    =   1.0
    alpha   =   s*1.0e-3
    emlenz  =   s/length


    eg(1,1)   =   alpha
    eg(2,2)   =   emlenz
    eg(3,3)   =   emlenz
    eg(4,4)   =   alpha
    eg(5,5)   =   0.0
    eg(6,6)   =   0.0

!
    eg(7,7)   =   eg(1,1)
    eg(8,8)   =   eg(2,2)
    eg(9,9)   =   eg(3,3)
    eg(10,10) =   eg(4,4)
    eg(11,11) =   eg(5,5)
    eg(12,12) =   eg(6,6)
!
    eg(1,7)   =   -eg(1,1)
    eg(2,8)   =   -eg(2,2)
    eg(2,12)  =    eg(2,6)
    eg(3,9)   =   -eg(3,3)
    eg(3,11)  =    eg(3,5)
    eg(4,10)  =   -eg(4,4)
    eg(5,9)   =   -eg(3,5)
    eg(5,11)  =   -eg(5,5)/4.0
    eg(6,8)   =   -eg(2,6)
    eg(6,12)  =   -eg(6,6)/4.0
!
    eg(8,12)  =   -eg(2,6)
    eg(9,11)  =   -eg(3,5)
!
!   impose the symmetry
    do  i= 1, 12
    do  j= i, 12
        eg(j,i) = eg(i,j)
    enddo
    enddo
!
!   check diagonal terms
    do    i=1,12
        if (abs(eg(i,i)) .lt. 1.0e-12) eg(i,i)=1.0e-12
    enddo
!
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ASSEMBle element matrices in COLumn form
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine assembCOL(nSTF,nND,nEQ,aa,a,i1,j1,k1,nloc)    
    implicit none
    integer:: nSTF,nND,nEQ            
    real(8):: aa(nSTF),a(18,18)
    integer:: idof(18)
    integer:: i1,j1,k1,nloc(nEQ)
    integer:: i,j,imax,jmax,ieqn1,ieqn2,jband,iloc    
!
!   Set idof to posible DoF number of each nodes
    if    (j1 .gt. 0  .AND. k1 .gt. 0) then
        imax=18
        jmax=18
    elseif (j1 .eq. 0  .AND. k1 .eq. 0) then
        imax=6
        jmax=6
    endif

    do  i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
!
!   Store the values for individual array in global array
    do  i= 1, imax
        ieqn1 = idof(i)
        do    j= i, jmax
            ieqn2 = idof(j)          
            if    (ieqn1 .gt. ieqn2) then
                    jband= (ieqn1-ieqn2)+1
                    iloc = nloc(ieqn1)
                    aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                   aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
            else
                    jband= (ieqn2-ieqn1)+1
                    iloc = nloc(ieqn2)
                    aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                   aa(ieqn1,jband) = aa(ieqn1,jband) + a(i,j)
            endif
        enddo
    enddo
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ASSEMBle LUMped element matrices
    subroutine assembLUM(nND,nEQ,aa,a,i1,j1,k1)
    implicit none
    integer:: nND,nEQ
    real(8):: aa(nEQ),a(18,18)          
    integer:: i1,j1,k1     
    integer:: idof(18),i,ieqn1
!
!   Set idof to posible DoF number of each nodes
    do  i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
!
!   Store the values for individual array in global array
    do  i= 1, 18
        ieqn1 = idof(i)
        aa(ieqn1) = aa(ieqn1) + a(i,i)
    enddo
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ASSEMBle consistent FORce for pressurized flex plate
    subroutine assembFOR(nEQ,nND, aa, a,i1, j1,k1)
    implicit none
    integer:: nEQ,nND   
    real(8):: aa(nEQ  ),a(18)
    integer:: i1,j1,k1
    integer:: idof(18),i,ieqn1
!
!   Set idof to posible DoF number of each nodes
    do    i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
!
!   Store the values for individual array in global array
    do  i= 1, 18
        ieqn1 = idof(i)
        aa(ieqn1) = aa(ieqn1) + a(i)
    enddo
!
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   TRANSformation of matrix in 3D. L->G
    subroutine trans3d_D(l,m,n,ek,beta)
    implicit none
    real(8):: ek(12,12),rt(3,3),r(3,3),ktemp(12,12)
    real(8):: m,n,l,beta,pi,sb,cb,d
    integer:: i,j,k,j1,j2,ii,jj,in,jn

    pi=4.0*datan(1.0d0)
!
    sb=dsin(beta*pi/180)
    cb=dcos(beta*pi/180)
    d=dsqrt(1.0-n**2)
!   if (abs(l).ge. 0.995 .and. abs(beta).le. 0.01) return
    if (abs(n).gt.0.995) then
           r(1,1)  =  0.0
           r(1,2)  =  0.0
           r(1,3)  =  n
           r(2,1)  = -n*sb
           r(2,2)  =  cb
           r(2,3)  =  0.0
           r(3,1)  = -n*cb
           r(3,2)  = -sb
           r(3,3)  =  0.0
    else
           r(1,1)  =  l
           r(1,2)  =  m
           r(1,3)  =  n
           if (abs(beta) .le. .01) then
              r(2,1)  =  -m/d
              r(2,2)  =  l/d
              r(2,3)  =  0.0
              r(3,1)  =  -l*n/d
              r(3,2)  =  -m*n/d
              r(3,3)  =  d
           else
              r(2,1)  =  -(m*cb+l*n*sb)/d
              r(2,2)  =  (l*cb-m*n*sb)/d
              r(2,3)  =  d*sb
              r(3,1)  =  (m*sb-l*n*cb)/d
              r(3,2)  =  -(l*sb+m*n*cb)/d
              r(3,3)  =  d*cb
           endif
    endif
!
    do  in=1,3
    do  jn=1,3
        rt(jn,in)=r(in,jn)
    enddo
    enddo
!    take [Rtrans][K][R] using the nature of [R] to speed computation.
!    k is sectioned off into 3x3s then multiplied [rtrans][k][r]
!
    do  i=0,3
    do  j=0,3
        do    k=1,3
        do    ii=1,3
            j1=i*3
            j2=j*3
            ktemp(j1+k,j2+ii)=0.0
            do     jj=1,3
            ktemp(j1+k,j2+ii)=ktemp(j1+k,j2+ii)+ek(j1+k,j2+jj)*r(jj,ii)
            enddo
        enddo
        enddo
        do  k=1,3
        do  ii=1,3
            ek(j1+k,j2+ii)=0.0
            do  jj=1,3
                ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rt(k,jj)*ktemp(j1+jj,j2+ii)
            enddo
        enddo
        enddo
    enddo
    enddo

    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   AxB product for COLumn storage
    subroutine AxBCOL(matrix,nSTF,vecin,vecout, nEQ,nBD,nprof,nprof2,nloc)
    implicit none
    integer:: nSTF,nEQ,nBD
    integer:: nprof(nEQ),nprof2(nEQ),nloc(nEQ)
    real(8):: vecout(nEQ),matrix(nSTF),vecin(nEQ)
    real(8):: val,valmat
    integer:: i,j,io,is,jlim,iloc
!
    do    i=1,nEQ
        jlim=max(1,(i-nprof(i)+1))
        do  j=jlim,i
            is=i
            io=i-j+1
            if (io .gt. nprof(is)) cycle
            iloc = nloc(is) + io -1
            valmat=matrix(iloc)
            val = vecin(j)
            vecout(i)=vecout(i) + val*valmat
        enddo

        jlim=min(nprof2(i),(nEQ-i+1))
        do  j=2,jlim
            is=i+j-1
            io=j
            if (io .gt. nprof(is)) cycle
            iloc = nloc(is) + io -1
            valmat=matrix(iloc)
            val = vecin(i+j-1)
            vecout(i)=vecout(i) + val*valmat
        enddo
    enddo


    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    makes a 3-D coordinate transformations.
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ROTate VECtor
    subroutine rotvec_D(xll,xmm,xnn,dx,dy,dz,dxb,dyb,dzb)
    implicit none
    real(8):: xll,xmm,xnn,dx,dy,dz,dxb,dyb,dzb
    real(8):: r(3,3)
    real(8):: ddd
!
    ddd=dsqrt(1.0-xnn**2)
    if    (abs(xnn) .gt. 0.9999) then
        r(1,1)  = +xnn
        r(1,2)  = +0.0
        r(1,3)  = -0.0
        r(2,1)  =  -0.0
        r(2,2)  =  1.0
        r(2,3)  =  0.0
        r(3,1)  =  0.0
        r(3,2)  =  0.0
        r(3,3)  =  xnn
        dxb = r(1,1)*dx
        dyb = r(2,2)*dy
        dzb = r(3,3)*dz
        return
    elseif (abs(xll) .gt. 0.9999) then
        r(1,1)  = +0.0
        r(1,2)  = +0.0
        r(1,3)  = -1.0
        r(2,1)  = -0.0
        r(2,2)  = +xll
        r(2,3)  =  0.0
        r(3,1)  = +xll
        r(3,2)  =  0.0
        r(3,3)  =  0.0
        dxb = r(1,3)*dz
        dyb = r(2,2)*dy
        dzb = r(3,1)*dx
        return
    else
        r(1,1)  = +xll*xnn/ddd
        r(1,2)  = +xmm*xnn/ddd
        r(1,3)  = -ddd
        r(2,1)  =  -xmm/ddd
        r(2,2)  =  +xll/ddd
        r(2,3)  =  0.0
        r(3,1)  =  xll
        r(3,2)  =  xmm
        r(3,3)  =  xnn
    endif

    dxb = r(1,1)*dx + r(1,2)*dy + r(1,3)*dz
    dyb = r(2,1)*dx + r(2,2)*dy + r(2,3)*dz
    dzb = r(3,1)*dx + r(3,2)*dy + r(3,3)*dz

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ROTate stiffness MATrix Local to Global
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rot_mat_LG(rr,ekb,ek)
    implicit none
    real(8):: ek(18,18),ekb(18,18)
    real(8):: rt(3,3),rr(3,3),ektemp(18,18)
    integer:: i,j,k,j1,j2,ii,kk,in,jn
!
!   take [R][k][Rt] using the nature of [R] to speed computation.
!   [k] is sectioned off into 6 3x3s then multiplied [rtrans][k][r]
!
!   get transpose
    do  in=1,3
    do  jn=1,3
        rt(jn,in)=rr(in,jn)
    enddo
    enddo
!
    do  i=0,5
    do  j=0,5
        j1=i*3
        j2=j*3
!
!       [k][R]
        do  k=1,3
            do  ii=1,3
                ektemp(j1+k,j2+ii)=0.0
                do  kk=1,3
                    ektemp(j1+k,j2+ii)=ektemp(j1+k,j2+ii)+ekb(j1+k,j2+kk)*rt(kk,ii)

                enddo
            enddo
        enddo
!

!        [R][k]
        do  k=1,3
            do  ii=1,3
                ek(j1+k,j2+ii)=0.0
                do  kk=1,3
                    ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rr(k,kk)*ektemp(j1+kk,j2+ii)

                enddo
            enddo
        enddo
    enddo
    enddo
!

    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ROTATE stiffness matrix
!   makes a 2-D coordinate transformations.
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rotate_D(xll,xmm,xnn,ekb,ek)   
    implicit none
    real(8):: xll,xmm,xnn
    real(8):: ek(18,18),ekb(18,18)
    real(8):: rt(3,3),r(3,3),ektemp(18,18)
    real(8):: ddd
    integer:: i,j,k,ii,jj,j1,j2,in,jn
!
    ddd=dsqrt(1.0-xnn**2)
    if    (abs(xnn) .gt. 0.9999) then
        r(1,1)  = +xnn
        r(1,2)  = +0.0
        r(1,3)  = -0.0
        r(2,1)  =  -0.0
        r(2,2)  =   1.0
        r(2,3)  =   0.0
        r(3,1)  =  0.0
        r(3,2)  =  0.0
        r(3,3)  =  xnn
    elseif (abs(xll) .gt. 0.9999) then
        r(1,1)  = +0.0
        r(1,2)  = +0.0
        r(1,3)  = -1.0
        r(2,1)  =  -0.0
        r(2,2)  =  +xll
        r(2,3)  =   0.0
        r(3,1)  = +xll
        r(3,2)  =  0.0
        r(3,3)  =  0.0
    else
        r(1,1)  = +xll*xnn/ddd
        r(1,2)  = +xmm*xnn/ddd
        r(1,3)  = -ddd
        r(2,1)  =  -xmm/ddd
        r(2,2)  =  +xll/ddd
        r(2,3)  =  0.0
        r(3,1)  =  xll
        r(3,2)  =  xmm
        r(3,3)  =  xnn
    endif
!
!  take [Rtrans][k][R] using the nature of [R] to speed computation.
!  [k] is sectioned off into 6 3x3s then multiplied [rtrans][k][r]
!
    if    (abs(xnn) .gt. 0.9999) then
        do    i=0,5
        do    j=0,5
            ii=i*3
            jj=j*3
            ek(ii+1,jj+1) = ekb(ii+1,jj+1)
            ek(ii+1,jj+2) = ekb(ii+1,jj+2)*r(1,1)
            ek(ii+1,jj+3) = ekb(ii+1,jj+3)
            ek(ii+2,jj+1) = ekb(ii+2,jj+1)*r(1,1)
            ek(ii+2,jj+2) = ekb(ii+2,jj+2)
            ek(ii+2,jj+3) = ekb(ii+2,jj+3)*r(1,1)
            ek(ii+3,jj+1) = ekb(ii+3,jj+1)
            ek(ii+3,jj+2) = ekb(ii+3,jj+2)*r(1,1)
            ek(ii+3,jj+3) = ekb(ii+3,jj+3)
        enddo
        enddo
        return
    elseif (abs(xll) .gt. 0.9999) then
        do    i=0,5
        do    j=0,5
            ii=i*3
            jj=j*3
            ek(ii+1,jj+1) = ekb(ii+3,jj+3)
            ek(ii+1,jj+2) = ekb(ii+3,jj+2)
            ek(ii+1,jj+3) =-ekb(ii+3,jj+1)*r(2,2)
            ek(ii+2,jj+1) = ekb(ii+2,jj+3)
            ek(ii+2,jj+2) = ekb(ii+2,jj+2)
            ek(ii+2,jj+3) =-ekb(ii+2,jj+1)*r(2,2)
            ek(ii+3,jj+1) =-ekb(ii+1,jj+3)*r(2,2)
            ek(ii+3,jj+2) =-ekb(ii+1,jj+2)*r(2,2)
            ek(ii+3,jj+3) = ekb(ii+1,jj+1)
        enddo
        enddo
        return
    endif
!
!   get transpose
    do    in=1,3
    do    jn=1,3
        rt(jn,in)=r(in,jn)
    enddo
    enddo
!
    do  i=0,5
    do  j=0,5
        j1=i*3
        j2=j*3
!
!        [k][R]
        do  k=1,3
        do  ii=1,3
            ektemp(j1+k,j2+ii)=0.0
            do  jj=1,3
                ektemp(j1+k,j2+ii)=ektemp(j1+k,j2+ii)+ekb(j1+k,j2+jj)*r(jj,ii)
            enddo
        enddo
        enddo
!
!       [Rtrans][k]
        do    k=1,3
        do  ii=1,3
            ek(j1+k,j2+ii)=0.0
            do  jj=1,3
                ek(j1+k,j2+ii)=ek(j1+k,j2+ii) +rt(k,jj)*ektemp(j1+jj,j2+ii)
            enddo
        enddo
        enddo

    enddo
    enddo
!
    return
    end

 
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   from Bergan and Felippa, CMAME, 1985
!   basic stiffness
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sm3mb_D(x,y,dm,alpha,f,ls,sm)
    implicit none
    real(8):: x(3),y(3),dm(3,3),alpha,f,sm(9,9)
    integer:: ls(9)

    real(8):: p(9,3)
    real(8):: area2,c
    real(8):: d11,d12,d13,d22,d23,d33
    real(8):: x21,x32,x13,y21,y32,y13
    real(8):: x12,x23,x31,y12,y23,y31
    real(8):: s1,s2,s3

    integer:: i,j,k,l,n
!
    x21 = x(2) - x(1)
    x32 = x(3) - x(2)
    x13 = x(1) - x(3)
    x12 = -x21
    x23 = -x32
    x31 = -x13
!
    y21 = y(2) - y(1)
    y32 = y(3) - y(2)
    y13 = y(1) - y(3)
    y12 = -y21
    y23 = -y32
    y31 = -y13
!
    area2  = y21*x13 - x21*y13
    p(1,1) = y23
    p(2,1) = 0.0
    p(3,1) = y31
    p(4,1) = 0.0
    p(5,1) = y12
    p(6,1) = 0.0
    p(1,2) = 0.0
    p(2,2) = x32
    p(3,2) = 0.0
    p(4,2) = x13
    p(5,2) = 0.0
    p(6,2) = x21
    p(1,3) = x32
    p(2,3) = y23
    p(3,3) = x13
    p(4,3) = y31
    p(5,3) = x21
    p(6,3) = y12
    n=6
!
    if (alpha .ne. 0.0) then
        p(7,1) = y23 * (y13-y21)*alpha/6.0
        p(7,2) = x32 * (x31-x12)*alpha/6.0
        p(7,3) =       (x31*y13-x12*y21)*alpha/3.0
        p(8,1) = y31 * (y21-y32)*alpha/6.0
        p(8,2) = x13 * (x12-x23)*alpha/6.0
        p(8,3) =       (x12*y21-x23*y32)*alpha/3.0
        p(9,1) = y12 * (y32-y13)*alpha/6.0
        p(9,2) = x21 * (x23-x31)*alpha/6.0
        p(9,3) =       (x23*y32-x31*y13)*alpha/3.0
        n=9
    endif
!
    c= 0.5*f/area2
    d11 = c*dm(1,1)
    d22 = c*dm(2,2)
    d33 = c*dm(3,3)
    d12 = c*dm(1,2)
    d13 = c*dm(1,3)
    d23 = c*dm(2,3)

    do  j=1,n
        l  = ls(j)
        s1 = d11*p(j,1) + d12*p(j,2) + d13*p(j,3)
        s2 = d12*p(j,1) + d22*p(j,2) + d23*p(j,3)
        s3 = d13*p(j,1) + d23*p(j,2) + d33*p(j,3)
        do  i=1,j
            k = ls(i)
            sm(k,l) = sm(k,l) + (s1*p(i,1) + s2*p(i,2) + s3*p(i,3))
            sm(l,k) = sm(k,l)
        enddo
    enddo
!
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   higher stiffness
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sm3mh_D(x,y,dm,f,ls,sm)
    implicit none
    real(8):: x(3),y(3),dm(3,3),f,sm(9,9)
    integer:: ls(9)

    real(8):: xc(3), yc(3), xm(3),ym(3)
    real(8):: gt(9,9),hh(3,9)
    real(8):: sqh(3,3), qx(3,3),qy(3,3)
    real(8):: area,area2,c
    real(8):: a1j, a2j, a3j, b1j, b2j,b3j
    real(8):: d11,d12,d13,d22,d23,d33,jxx,jxy,jyy
    real(8):: x0, y0, xi,yi
    real(8):: cj, sj,dl,dx,dy
    real(8):: s1,s2,s3,s4,s5,s6
    integer:: iperm(9)
    real(8):: gti(9,9),wk(9,9)
    integer:: i,j,k,l,n
!
!
    area2 = (y(2) -y(1))*(x(1)-x(3)) - (x(2)-x(1))*(y(1)-y(3))
!
    x0 = (x(1)+x(2)+x(3))/3.0
    y0 = (y(1)+y(2)+y(3))/3.0
    area=0.5*area2
    c = 1.0/dsqrt(area)
!
    xc(1) = c*(x(1)-x0)
    xc(2) = c*(x(2)-x0)
    xc(3) = c*(x(3)-x0)
    yc(1) = c*(y(1)-y0)
    yc(2) = c*(y(2)-y0)
    yc(3) = c*(y(3)-y0)
    xm(1) = 0.5*(xc(2)+xc(3))
    xm(2) = 0.5*(xc(3)+xc(1))
    xm(3) = 0.5*(xc(1)+xc(2))
    ym(1) = 0.5*(yc(2)+yc(3))
    ym(2) = 0.5*(yc(3)+yc(1))
    ym(3) = 0.5*(yc(1)+yc(2))
!
!   form G^T in GT and initialize HH
!
    do  i=1,9
        do  j=1,6
            gt(j,i) = 0
        enddo
        hh(1,i)=0
        hh(2,i)=0
        hh(3,i)=0
    enddo
!
    d11 = f*dm(1,1)
    d22 = f*dm(2,2)
    d33 = f*dm(3,3)
    d12 = f*dm(1,2)
    d13 = f*dm(1,3)
    d23 = f*dm(2,3)
    jxx = -2.0*(xc(1)*xc(2) + xc(2)*xc(3) + xc(3)*xc(1))/3.0
    jxy =      (xc(1)*yc(1) + xc(2)*yc(2) + xc(3)*yc(3))/3.0
    jyy = -2.0*(yc(1)*yc(2) + yc(2)*yc(3) + yc(3)*yc(1))/3.0
    do    j=1,3
        dx = xm(j) - xc(j)
        dy = ym(j) - yc(j)
        dl = dsqrt(dx*dx + dy*dy)
        cj = dx/dl
        sj = dy/dl
!
!   !!!a2j b2j different than paper
        a1j = -0.5*sj*cj**2
        a2j =  0.5*cj**3
        b2j = -0.5*sj**3
        b3j =  0.5*sj**2*cj
        a3j = -(b2j + a1j + a1j)
        b1j = -(b3j + b3j + a2j)
!
        gt(1,2*j-1) =   1.
        gt(2,2*j  ) =   1.
        gt(3,2*j-1) =  -yc(j)
        gt(3,2*j  ) =  xc(j)
        gt(3,  j+6) =   c
        gt(4,2*j-1) =  xc(j)
        gt(6,2*j-1) =  yc(j)
        gt(5,2*j  ) =  yc(j)
        gt(6,2*j  ) =  xc(j)
        hh(j,j+6)   =   1.
        qx(j,1)     =   a1j
        qx(j,2)     =   b2j
        qx(j,3)     =  -2.0*b3j
        qy(j,1)     =   a2j
        qy(j,2)     =   b3j
        qy(j,3)     = -2.0*a1j
        s1 =         d11*qx(j,1)    + d12*qx(j,2)   + d13*qx(j,3)
        s2 =         d12*qx(j,1)    + d22*qx(j,2)   + d23*qx(j,3)
        s3 =         d13*qx(j,1)    + d23*qx(j,2)   + d33*qx(j,3)
        s4 =         d11*qy(j,1)    + d12*qy(j,2)   + d13*qy(j,3)
        s5 =         d12*qy(j,1)    + d22*qy(j,2)   + d23*qy(j,3)
        s6 =         d13*qy(j,1)    + d23*qy(j,2)   + d33*qy(j,3)
        do  i = 1,3
            xi = xc(i)
            yi = yc(i)
            gt(j+6,2*i-1) =    a1j*xi*xi + 2.*a2j*xi*yi + a3j*yi*yi
            gt(j+6,2*i)   =    b1j*xi*xi + 2.*b2j*xi*yi + b3j*yi*yi
            gt(j+6,i+6)   =   -c*(cj*xi+sj*yi)
        enddo
        do  i=1,j
            sqh(i,j) = jxx*( qx(i,1)*s1+qx(i,2)*s2+qx(i,3)*s3) &
                     + jxy*( qx(i,1)*s4+qx(i,2)*s5+qx(i,3)*s6+qy(i,1)*s1+qy(i,2)*s2+qy(i,3)*s3) &
                     + jyy*( qy(i,1)*s4+qy(i,2)*s5+qy(i,3)*s6)
        enddo
    enddo
!
!   Factor G' and backsolve to obtain H
!   Form physical stiffness and add to incoming SM
    do  i=1,9
    do  j=1,9
        gti(i,j)=gt(i,j)
    enddo
    enddo

    call ainver(gt,9,iperm,wk)

    do  i=1,3
    do  j=1,9
        hh(i,j)=gt(j,i+6)
    enddo
    enddo

    do  j = 1,9
        l  = ls(j)
        s1 = sqh(1,1)*hh(1,j) + sqh(1,2)*hh(2,j) + sqh(1,3)*hh(3,j)
        s2 = sqh(1,2)*hh(1,j) + sqh(2,2)*hh(2,j) + sqh(2,3)*hh(3,j)
        s3 = sqh(1,3)*hh(1,j) + sqh(2,3)*hh(2,j) + sqh(3,3)*hh(3,j)
        do  i = 1,j
            k = ls(i)
            sm(k,l) = sm(k,l) + (s1*hh(1,i) + s2*hh(2,i) + s3*hh(3,i))
            sm(l,k) = sm(k,l)
        enddo
    enddo

    return
    end



!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   UDU decomposition of COLumn profiled system
    subroutine uduCOL_D(a,maxstore,nEQ,nBD,imult,nprof,nloc)
    implicit none
    integer:: maxstore,nEQ,nBD,imult
    real(8):: a(maxstore)    
    integer:: nprof(nEQ), nloc(nEQ)

    real(8):: temp,sum
    integer:: i,j,k,j2,j3,im1,jm1,is,io,iloc,jcol,iloci,ilocj,iloc1

    if    (a(1) .eq. 0.0d0) then
        imult=0
        return
    endif
!
    if    (nEQ .eq. 1) then
        imult=1
        return
    endif
!

    do  j=2,nEQ
        jm1=j-1
        j2=j-nprof(j)+1
        if    (j2.lt.1) then
            j2=1
        endif
!
!       off-diagonal terms
        if    (jm1.eq.1) then
            is=j
            io=1
            iloc = nloc(is) + io - 1
            sum=a(iloc)
!           sum=a(j,1)
        else
            do  i=j2+1,jm1
                im1=i-1
                is=j
                io=j-i+1
                iloc = nloc(is) + io - 1
                sum=a(iloc   )
!               sum=a(i,j-i+1)
!
                j3=i-nprof(i)+1
                jcol=j3
                if    (j3 .lt. j2) then
                    jcol=j2
                endif
!               do    k=j2,im1
                do  k=jcol,im1
                    is=i
                    io=i-k+1
                    iloci = nloc(is) + io - 1
                    is=j
                    io=j-k+1
                    ilocj = nloc(is) + io - 1
                    sum=sum-a(iloci  )*a(ilocj  )
!                   sum=sum-a(k,i-k+1)*a(k,j-k+1)
                    imult=imult+1
                enddo
                a(iloc   )=sum
!               a(i,j-i+1)=sum
            enddo
            is=j
            io=1
            iloc = nloc(is) + io - 1
            sum=a(iloc   )
!           sum=a(j,1)
        endif
!
!        diagonal terms
        do  k=j2,jm1
            is=j
            io=j-k+1
            ilocj = nloc(is) + io - 1
            is=k
            io=1
            iloc1 = nloc(is) + io - 1
            temp=a(ilocj  )/a(iloc1)
            sum=sum-temp*a(ilocj  )
            a(ilocj  )=temp
!           temp=a(k,j-k+1)/a(k,1)
!           sum=sum-temp*a(k,j-k+1)
!           a(k,j-k+1)=temp
            imult=imult+2
        enddo

        if (sum.eq.0.0d0) then
            imult=0
            return
        endif
        is=j
        io=1
        iloc = nloc(is) + io - 1
        a(iloc   ) = sum
!       a(j,1)=sum
    enddo
!
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   BAcK solver of COLumn profiled system
    subroutine bakCOL_D(a,maxstore,b,nEQ,nBD,wk,imult,nprof,nprof2,nloc)
    implicit none
    integer:: maxstore,nEQ,nBD,imult
    real(8):: a(maxstore),b(nEQ),wk(nEQ)    
    integer:: nprof(nEQ), nprof2(nEQ),nloc(nEQ)

    real(8):: sum
    integer:: i,j,k,i1,k2,jb,jbb,km1,is,io,iloc

!
!   forward substitutions
    do    i=1,nEQ
!        j=i-nBD+1
        j=i-nprof(i)+1
        if    (i.le.nprof(i) ) then
            j=1
        endif
        jb=i-nBD+1
        jbb=jb
        if    (i.le.nBD    ) then
            jbb=1
        endif
        sum=b(i)
        km1=i-1
        if    (j.gt.km1) then
            wk(i)=sum
        else
            do  k=j,km1
                is=i
                io=i-k+1
                iloc = nloc(is) + io - 1
                sum=sum-a(iloc   )*wk(k)
!               sum=sum-a(k,i-k+1)*wk(k)
                imult=imult+1
            enddo
            wk(i)=sum
        endif
    enddo
!
!   middle terms
    do  i=1,nEQ
        is=i
        io=1
        iloc = nloc(is) + io - 1
        wk(i)=wk(i)/a(iloc )
!       wk(i)=wk(i)/a(i,1)
        imult=imult+1
    enddo
!
!   backward substitution
    do  i1=1,nEQ
        i=nEQ-i1+1
        j=i+nprof2(i) -1
        if    (j.gt.nEQ) then
            j=nEQ
        endif
        jb=i+nBD-1
        jbb=jb
        if    (jb.gt.nEQ) then
            jbb=nEQ
        endif
        sum=wk(i)
        k2=i+1
        if    (k2.gt.j) then
            wk(i)=sum
        else
            do  k=k2,j
                is=k
                io=k-i+1
                if (io .gt. nprof(is)) cycle
                iloc = nloc(is) + io - 1
                sum=sum-a(iloc   )*wk(k)
!               sum=sum-a(i,k-i+1)*wk(k)
                imult=imult+1
            enddo

            wk(i)=sum
        endif
    enddo
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ainver(a,n,indx,yn)
    implicit none
    integer:: n
    real(8):: a(n,n),yn(n,n)
    integer:: indx(n)
    real(8):: d
    integer:: i,j,np

    do    i = 1,n
        do    j = 1,n
            yn(i,j) = 0.0
        enddo
        yn(i,i) = 1.0
    enddo
!
    np=n
    call ludcmp(a,n,np,indx,d)
!
    do  j = 1,n
        call lubksb(a,n,np,indx,yn(1,j))
    enddo
!
    do  j = 1,n
    do  i = 1,n
        a(i,j) = yn(i,j)
    enddo
    enddo
!
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ludcmp(a,n,np,indx,d)
    implicit none
    integer:: n,np
    integer,parameter:: nmax=100
    real(8),parameter:: tiny=1.0e-20
    real(8):: a(np,np), vv(nmax)
    real(8):: d,sum,aamax,dum
    integer:: indx(n)
    integer:: i,j,k,imax
!
    d = 1.0
    do    i = 1,n
        aamax = 0.0
        do  j = 1,n
            if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
        enddo
        if (aamax .eq. 0.0) pause 'Singular Matrix'
        vv(i) = 1.0/aamax
    enddo

    do  j = 1,n
        do  i = 1,j-1
            sum = a(i,j)
            do   k = 1,i-1
               sum = sum - a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
        enddo

        aamax = 0.0
        do  i = j,n
            sum = a(i,j)
            do    k = 1,j-1
                sum = sum - a(i,k)*a(k,j)
            enddo
            a(i,j) = sum
            dum = vv(i)*abs(sum)
            if (dum .ge. aamax) then
               imax = i
               aamax = dum
            endif
        enddo

        if    (j .ne. imax) then
            do  k = 1,n
                dum = a(imax,k)
                a(imax,k) = a(j,k)
                a(j,k) = dum
            enddo
            d = -d
            vv(imax) = vv(j)
        endif

        indx(j) = imax
        if (a(j,j) .eq. 0.0) a(j,j) = tiny
        if (j .ne. n) then
            dum = 1.0/a(j,j)
            do    i = j+1,n
                a(i,j) = a(i,j)*dum
            enddo
         endif

    enddo
!
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lubksb(a,n,np,indx,b)
    implicit none
    integer:: n,np
    real(8):: a(np,np), b(n)
    real(8):: sum
    integer:: indx(n)
    integer:: i,j,ll,ii
!
    ii = 0
    do  i = 1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if    (ii .ne. 0) then
            do  j = ii , i-1
                sum = sum - a(i,j)*b(j)
            enddo
        else if (sum .ne. 0.0) then
            ii = i
        endif
        b(i) = sum
    enddo
    do  i = n,1,-1
        sum = b(i)
        if    (i .lt. n) then
            do  j = i+1,n
                sum = sum - a(i,j)*b(j)
            enddo
        endif
        b(i) = sum/a(i,i)
    enddo
    return
    end


