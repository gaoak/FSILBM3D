!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_file()
    USE simParam
    implicit none
    integer:: i
    open(unit=111,file='inFlow.dat') 
    read(111,*)     !======================================
    read(111,*)     npsize
    read(111,*)     isRelease
    read(111,*)     isConCmpt,  iCollidModel
    read(111,*)     iFlapRef ,  iKB
    read(111,*)     timeSimTotl,timeOutTemp
    read(111,*)     timeOutFlEd,timeOutFlBg
    read(111,*)     timeOutFlow,timeOutBody,timeOutInfo
    read(111,*)     Palpha,Pbeta,Pramp 
    read(111,*)     !======================================
    read(111,*)     uuuIn(1:3)
    read(111,*)     denIn     
    read(111,*)     Re             
    read(111,*)     !======================================
    read(111,*)     LBmeshName
    read(111,*)     dt    
    read(111,*)     !======================================
    read(111,*)     iBC 
    read(111,*)     !======================================
    read(111,*)     dtolLBM,ntolLBM
    read(111,*)     !======================================
    read(111,*)     FEmeshName
    read(111,*)     iBodyModel,iForce2Body,iChordDirection
    read(111,*)     isMotionGiven(1:3)
    read(111,*)     isMotionGiven(4:6)
    read(111,*)     !======================================
    read(111,*)     denR,  psR
    if(iKB==0) read(111,*)     EmR, tcR  
    if(iKB==1) read(111,*)     KB,    KS
    read(111,*)     !======================================
    read(111,*)     numsubstep
    read(111,*)     dampK,dampM
    read(111,*)     NewmarkGamma,NewmarkBeta
    read(111,*)     alphaf,alpham,alphap
    read(111,*)     dtolFEM,ntolFEM  
    read(111,*)     !======================================
    read(111,*)     Freq,   St
    read(111,*)     XYZo(1:3)
    read(111,*)     XYZAmpl(1:3)
    read(111,*)     XYZPhi(1:3)
    read(111,*)     AoAo(1:3)
    read(111,*)     AoAAmpl(1:3)
    read(111,*)     AoAPhi(1:3)
    read(111,*)     !======================================
    read(111,*)     Frod(1:3)
    read(111,*)     !======================================
    read(111,*)     isMoveGrid,      numOutput
    read(111,*)     isMoveDimX,      isMoveOutputX
    read(111,*)     isMoveDimY,      isMoveOutputY
    read(111,*)     isMoveDimZ,      isMoveOutputZ
    read(111,*)     Xref,Yref,Zref
    read(111,*)     !======================================
    read(111,*)     waveInitDist,AmplInitDist(1:SpcDim)
    read(111,*)     FreqForcDist,AmplForcDist(1:SpcDim)
    read(111,*)     begForcDist,endForcDist
    read(111,*)     posiForcDist(1:SpcDim)
    read(111,*)     !======================================
    read(111,*)     numSampFlow
    allocate(SampFlowPint(1:numSampFlow,1:3))
    do  i=1,     numSampFlow
        read(111,*) SampFlowPint(i,1:3)
    enddo
    read(111,*)     !======================================
    read(111,*)     numSampBody
    allocate(SampBodyNode(1:numSampBody))
    read(111,*) SampBodyNode(1:numSampBody)
    read(111,*)     !======================================

    close(111)
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    allocate memory for fluid simulation
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE allocate_fluid_memory()
    USE simParam
    USE PartitionXDim
    implicit none
    integer,parameter::idFile=111
    integer:: x,y,z,temp,allocatemsg
     
!   read initial grid**********************************************************************************
    open(unit=idFile,file=trim(adjustl(LBmeshName)))
    read(idFile,*)xDim
    allocate(xGrid0(xDim),xGrid(xDim),dx(xDim))
    do x=1,xDim
    read(idFile,*)temp,xGrid0(x)   
    enddo
    read(idFile,*)yDim
    allocate(yGrid0(yDim),yGrid(yDim),dy(yDim))
    do y=1,yDim
    read(idFile,*)temp,yGrid0(y)   
    enddo
    read(idFile,*)zDim
    allocate(zGrid0(zDim),zGrid(zDim),dz(zDim))
    do z=1,zDim
    read(idFile,*)temp,zGrid0(z)   
    enddo
    close(idFile)

!   allocate fluid memory**********************************************************************************
    allocate(fIn(zDim,yDim,xDim,0:18),fInTemp(zDim,yDim,xDim))
    allocate(uuu(zDim,yDim,xDim,1:3),force(zDim,yDim,xDim,1:3))
    allocate(den(zDim,yDim,xDim),prs(zDim,yDim,xDim),image(zDim,yDim,xDim))

!   calculate grid size**********************************************************************************    
    dx(1)=dabs(xGrid0(2)-xGrid0(1))
    do  x=2,xDim-1
        dx(x)=dabs(0.5*(xGrid0(x)+xGrid0(x+1))-0.5*(xGrid0(x-1)+xGrid0(x)))
    enddo
    dx(xDim)=dabs(xGrid0(xDim)-xGrid0(xDim-1))

    dy(1)=dabs(yGrid0(2)-yGrid0(1))
    do  y=2,yDim-1
        dy(y)=dabs(0.5*(yGrid0(y)+yGrid0(y+1))-0.5*(yGrid0(y-1)+yGrid0(y)))
    enddo
    dy(yDim)=dabs(yGrid0(yDim)-yGrid0(yDim-1))

    dz(1)=dabs(zGrid0(2)-zGrid0(1))
    do  z=2,zDim-1
        dz(z)=dabs(0.5*(zGrid0(z)+zGrid0(z+1))-0.5*(zGrid0(z-1)+zGrid0(z)))
    enddo
    dz(zDim)=dabs(zGrid0(zDim)-zGrid0(zDim-1))

    dxmin=dabs(minval(dx(1:xDim)))
    dymin=dabs(minval(dy(1:yDim)))
    dzmin=dabs(minval(dz(1:zDim)))

    dxmax=dabs(maxval(dx(1:xDim)))
    dymax=dabs(maxval(dy(1:yDim)))
    dzmax=dabs(maxval(dz(1:zDim)))

    cptxmin=xGrid0(1)
    cptxmax=xGrid0(xDim)
    cptymin=yGrid0(1)
    cptymax=yGrid0(yDim)
    cptzmin=zGrid0(1)
    cptzmax=zGrid0(zDim)

    if(dabs(dxmin-dymin)>eps) then
    write(*,*)'dxmin=',dxmin
    write(*,*)'dymin=',dymin
    write(*,*)'dxmin/=dymin'
    stop
    endif
    if(dabs(dxmin-dzmin)>eps) then
    write(*,*)'dxmin=',dxmin
    write(*,*)'dzmin=',dzmin
    write(*,*)'dxmin/=dzmin'
    stop
    endif
    if(dabs(dzmin-dymin)>eps) then
    write(*,*)'dzmin=',dzmin
    write(*,*)'dymin=',dymin
    write(*,*)'dzmin/=dymin'
    stop
    endif
    dh=dxmin

    !!!!!!!!!!!!! mannual partition in x-direction
    npsize_copy = npsize
    xDim_copy = xDim
    allocate(partition(1:npsize),parindex(1:npsize+1),eid(1:npsize), STAT=allocatemsg)
    if(allocatemsg .ne. 0) write(*, *) 'Allocation of partition in swapx failed'
    allocate(edge(1:zDim,1:yDim, 1:npsize), STAT=allocatemsg)
    if(allocatemsg .ne. 0) write(*, *) 'Allocation of edge in swapx failed'
    call OMPPartition(xDim_copy, npsize_copy, partition, parindex)
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Allocate memory for solid simulation
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE allocate_solid_memory()
    USE simParam
    implicit none
    integer:: iND
    real(8):: xCT,yCT,zCT,xl,xr,yl,yr,zl,zr
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az

    write(*,'(A)') '=============================================================================='
    open(unit=idat, file = trim(adjustl(FEmeshName)))
    rewind(idat)
    read(idat,*)  
    read(idat,*)nND,nEL,nMT 
    read(idat,*)             
!   ===============================================================================================
    allocate( ele(nEL,5),nloc(nND*6),nprof(nND*6),nprof2(nND*6),jBC(nND,6))
    allocate( xyzful00(nND,6),xyzful0(nND,6),mssful(nND,6),lodful(nND,6),extful(nND,6),extful1(nND,6),streI(nND),bendO(nND))
    allocate( extful2(nND,6),grav(nND,6),vBC(nND,6),mss(nND*6),prop(nMT,10),areaElem00(nEL),areaElem(nEL))

    allocate( xyzful(nND,6),xyzfulIB(nND,6),xyzfulnxt(nND,6),dspful(nND,6),velful(nND,6),accful(nND,6)) 
    allocate( triad_nn(3,3,nND),triad_ee(3,3,nEL),triad_e0(3,3,nEL) )
    allocate( triad_n1(3,3,nEL),triad_n2(3,3,nEL),triad_n3(3,3,nEL) ) 
            
!   ===============================================================================================     
    call read_structural_datafile(jBC,ele,nloc,nprof,nprof2,xyzful00,prop,nND,nEL,nEQ,nMT,nBD,nSTF,idat)
    close(idat)
!    ===============================================================================================
    !calculate area
    call cptArea(areaElem00,nND,nEL,ele,xyzful00)
    Asfac=sum(areaElem00(:))
    elmax=maxval(areaElem00(:))
    elmin=minval(areaElem00(:))
    streI(1:nND)=0.0d0
    bendO(1:nND)=0.0d0
    !calculate spanwise length, chord length, aspect ratio
!    if(iChordDirection==1)then      
!        Lchod   = maxval(xyzful00(:,1))-minval(xyzful00(:,1))
!        Lspan   = maxval(xyzful00(:,2))-minval(xyzful00(:,2))
!    elseif(iChordDirection==2)then
!        Lchod   = maxval(xyzful00(:,2))-minval(xyzful00(:,2))
!        Lspan   = maxval(xyzful00(:,1))-minval(xyzful00(:,1))        
!    else
!        stop 'no define chordwise'
!    endif

!    if((Lchod-1.0d0)<=1.0d-2)Lchod=1.0d0
!    if((Lspan-1.0d0)<=1.0d-2)Lspan=1.0d0
    Lchod=1.0d0 
    Lspan=1.0d0
    AR     = Lspan**2/Asfac

    !calculate central position, central node
    xCT=sum(xyzful00(:,1))/nND
    yCT=sum(xyzful00(:,2))/nND
    zCT=sum(xyzful00(:,3))/nND


    !calculate leading-edge central node, three trailing nodes, the nearest neighbouring node
    xl=minval(xyzful00(:,1))
    xr=maxval(xyzful00(:,1))
    yl=minval(xyzful00(:,2))
    yr=maxval(xyzful00(:,2))
    zl=minval(xyzful00(:,3))
    zr=maxval(xyzful00(:,3))

    NDct   =minloc(dsqrt((xyzful00(1:nND,1)-xCT)**2+(xyzful00(1:nND,2)-        yCT)**2),1)
    
    if(trim(adjustl(FEmeshName))=="fi30.dat")then 
        NDtl(1)=513
        NDtl(2)=221
        NDtl(3)=536
        NDhd(1)=2
        NDhd(2)=108
        NDhd(3)=392
    elseif(trim(adjustl(FEmeshName))=="fi45.dat")then
        NDtl(1)=502
        NDtl(2)=309
        NDtl(3)=530
        NDhd(1)=2
        NDhd(2)=109
        NDhd(3)=417
    elseif(trim(adjustl(FEmeshName))=="fi60.dat")then
        NDtl(1)=472
        NDtl(2)=342
        NDtl(3)=524
        NDhd(1)=2
        NDhd(2)=112
        NDhd(3)=417
    elseif(trim(adjustl(FEmeshName))=="fi75.dat")then
        NDtl(1)=469
        NDtl(2)=445
        NDtl(3)=533
        NDhd(1)=2
        NDhd(2)=110
        NDhd(3)=421
    elseif(trim(adjustl(FEmeshName))=="fi90.dat")then
    NDtl(1)=minloc(dsqrt((xyzful00(1:nND,1)-xr )**2+(xyzful00(1:nND,2)-         yl)**2),1)
    NDtl(2)=minloc(dsqrt((xyzful00(1:nND,1)-xr )**2+(xyzful00(1:nND,2)-0.5*(yl+yr))**2),1)
    NDtl(3)=minloc(dsqrt((xyzful00(1:nND,1)-xr )**2+(xyzful00(1:nND,2)-         yr)**2),1)
    
    NDtl(4)=minloc(dsqrt((xyzful00(1:nND,1)-0.5*(xl+xr) )**2+(xyzful00(1:nND,2)-yl)**2),1)
    NDtl(5)=minloc(dsqrt((xyzful00(1:nND,1)-0.5*(xl+xr) )**2+(xyzful00(1:nND,2)-yr)**2),1)
    
    NDhd(1)=minloc(dsqrt((xyzful00(1:nND,1)-xl )**2+(xyzful00(1:nND,2)-         yl)**2),1)
    NDhd(2)=minloc(dsqrt((xyzful00(1:nND,1)-xl )**2+(xyzful00(1:nND,2)-0.5*(yl+yr))**2),1)
    NDhd(3)=minloc(dsqrt((xyzful00(1:nND,1)-xl )**2+(xyzful00(1:nND,2)-         yr)**2),1)
    elseif(trim(adjustl(FEmeshName))=="fi105.dat")then
        NDtl(1)=369
        NDtl(2)=491
        NDtl(3)=530
        NDhd(1)=2
        NDhd(2)=109
        NDhd(3)=425
    elseif(trim(adjustl(FEmeshName))=="fi120.dat")then
        NDtl(1)=330
        NDtl(2)=523
        NDtl(3)=534
        NDhd(1)=2
        NDhd(2)=112
        NDhd(3)=422
    elseif(trim(adjustl(FEmeshName))=="fi135.dat")then
        NDtl(1)=237
        NDtl(2)=528
        NDtl(3)=509
        NDhd(1)=2
        NDhd(2)=109
        NDhd(3)=414
    elseif(trim(adjustl(FEmeshName))=="fi150.dat")then
        NDtl(1)=150
        NDtl(2)=515
        NDtl(3)=441
        NDhd(1)=2
        NDhd(2)=112
        NDhd(3)=400
    endif

!   loading boundary type*******************************************************************************
    do    iND=1,nND
        if(jBC(iND,1)==1)jBC(iND,1:6)=isMotionGiven(1:6)
    enddo

    END SUBROUTINE 

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    initialize flow field
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initialize_flow()
    USE simParam
    implicit none
    real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim)
    integer:: x, y, z

!   grid coordinate***************************************************************************************         
    xGrid(1:xDim)=xGrid0(1:xDim)
    yGrid(1:yDim)=yGrid0(1:yDim)
    zGrid(1:zDim)=zGrid0(1:zDim)
!   macro quantities***************************************************************************************    
    uuu(1:zDim,1:yDim,1:xDim,1) = uuuIn(1)
    uuu(1:zDim,1:yDim,1:xDim,2) = uuuIn(2)
    uuu(1:zDim,1:yDim,1:xDim,3) = uuuIn(3)
    den(1:zDim,1:yDim,1:xDim)   = denIn
    prs(1:zDim,1:yDim,1:xDim)   = Cs2*(den(1:zDim,1:yDim,1:xDim)-denIn)
!   initial disturbance***************************************************************************************
    call initDisturb()
!   distribution function***************************************************************************************
    do  x = 1, xDim
    do  y = 1, yDim
    do  z = 1, zDim
        uSqr       = sum(uuu(z,y,x,1:3)**2)
        uxyz(0:lbmDim) = uuu(z,y,x,1) * ee(0:lbmDim,1) + uuu(z,y,x,2) * ee(0:lbmDim,2)+uuu(z,y,x,3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)= wt(0:lbmDim) * den(z,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        fIn(z,y,x,0:lbmDim)=fEq(0:lbmDim)
    enddo
    enddo
    enddo
!   grid type***************************************************************************************
    image(1:zDim,1:yDim,1:xDim)  = fluid
    image(1:zDim,1:yDim,1     )  = xMinBC
    image(1:zDim,1:yDim,xDim  )  = xMaxBC
    image(1:zDim,     1,1:xDim)  = yMinBC   
    image(1:zDim,  yDim,1:xDim)  = yMaxBC
    image(1     ,1:yDim,1:xDim)  = zMinBC
    image(zDim  ,1:yDim,1:xDim)  = zMaxBC
   
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    initialize solid field
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initialize_solid()
    USE simParam
    implicit none
    integer:: iND

    TTT00(:,:)=0.0d0
    TTT00(1,1)=1.0d0
    TTT00(2,2)=1.0d0
    TTT00(3,3)=1.0d0

    XYZ(1:3)=XYZo(1:3)+XYZAmpl(1:3)*dcos(2.0*pi*Freq*time+XYZPhi(1:3))
    AoA(1:3)=AoAo(1:3)+AoAAmpl(1:3)*dcos(2.0*pi*Freq*time+AoAPhi(1:3))
    call AoAtoTTT(AoA,TTT0)
    call AoAtoTTT(AoA,TTTnow)
    call get_angle_triad(TTT0,TTTnow,AoAd(1),AoAd(2),AoAd(3))
    do  iND=1,nND
            xyzful0(iND,1:3)=matmul(TTT0,xyzful00(iND,1:3))+XYZ(1:3)
            xyzful0(iND,4:6)=AoAd(1:3)
    enddo

    xyzful(1:nND,1:6)=xyzful0(1:nND,1:6)     
    velful(1:nND,1:6)=0.0
                
    dspful(1:nND,1:6)=0.0
    accful(1:nND,1:6)=0.0

    xyzfulIB(1:nND,1:6)=xyzful(1:nND,1:6)

    CALL formmass_D(  ele,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),prop,mss,nND,nEL,nEQ,nMT,alphaf,alpham,alphap ) 
    do  iND = 1, nND
        mssful(iND,1:6)=mss((iND-1)*6+1:(iND-1)*6+6)
        grav(iND,1:6) = mssful(iND,1:6)*[g(1),g(2),g(3),0.0d0,0.0d0,0.0d0]
    enddo
    CALL init_triad_D(ele,xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3),triad_nn,triad_n1,triad_n2,triad_n3,triad_ee,triad_e0,nND,nEL) 

    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    calculate Bolztman parameters from flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE calculate_LB_params()
    USE simParam
    implicit none
    integer:: nt
!   angle to radian
    AoAo(1:3)   = AoAo(1:3)/180.0*pi
    AoAAmpl(1:3)= AoAAmpl(1:3)/180.0*pi
    AoAPhi(1:3) = AoAPhi(1:3)/180.0*pi
!   reference values: length, velocity, time
    Lref    = Lchod
    if(maxval(dabs(uuuIn(1:3)))>0.0001d0)then !with incoming flow
        Uref    = maxval(dabs(uuuIn(1:3)))        
        if    (MAXVAL(dabs(xyzAmpl(1:3)))>0.0001d0)then
            !Freq    = St*Uref/(2.0*MAXVAL(dabs(xyzAmpl(1:3))))
            Tref    = 1.0d0/Freq
            St      =  2.0*MAXVAL(dabs(xyzAmpl(1:3)))*Freq/maxval(dabs(uuuIn(1:3))) ! St=2fA/U
            if(iFlapRef==1) Uref    = Lref*Freq           
        elseif(dabs(MAXVAL(AoAAmpl(1:3)))>0.0001d0)then
            Tref    = 1.0d0/Freq            
        else
            Tref    = Lref/Uref            
        endif
        g(1:3)=Frod(1:3)*Uref**2/Lref
    else                                    !without incoming flow
        if    (MAXVAL(dabs(xyzAmpl(1:3)))>0.0001d0)then !heaving 
        !Uref    = 2.0*pi*MAXVAL(dabs(xyzAmpl(1:3)))*Freq  !maximum velocity
        !Uref    = 4.0*   MAXVAL(dabs(xyzAmpl(1:3)))*Freq  !mean velocity
        !Uref    =        MAXVAL(dabs(xyzAmpl(1:3)))*Freq  !mean velocity
        Uref    = Lref*Freq
        Tref    = 1.0d0/Freq
        St      = 1000.0
        g(1:3)=Frod(1:3)*Uref**2/Lref       
        elseif(MAXVAL(dabs(AoAAmpl(1:3)))>0.0001d0)then !picthing
        !Uref = 2.0*pi*MAXVAL(dabs(AoAAmpl(1:3))*[maxval(dabs(xyzful00(:,2))),maxval(dabs(xyzful00(:,1))),maxval(dabs(xyzful00(:,3)))])*Freq  !maximum velocity
        Uref  = 4.0*   MAXVAL(dabs(AoAAmpl(1:3))*[maxval(dabs(xyzful00(:,2))),maxval(dabs(xyzful00(:,1))),maxval(dabs(xyzful00(:,3)))])*Freq  !mean velocity
        !Uref    =     MAXVAL(dabs(AoAAmpl(1:3))*[maxval(dabs(xyzful00(:,2))),maxval(dabs(xyzful00(:,1))),maxval(dabs(xyzful00(:,3)))])*Freq  !mean velocity
        Tref    = 1.0d0/Freq
        St      = 1000.0
        g(1:3)=Frod(1:3)*Uref**2/Lref 
        elseif(maxval(dabs(Frod(1:3)))>0.0001d0)then    !gravity
            Uref=maxval(dabs(Frod(1:3)))
            g(1:3)=(Frod(1:3))**2/(Lref*denR)
            if(Frod(1)<0.000001d0) g(1)=-g(1)
            if(Frod(2)<0.000001d0) g(2)=-g(2)
            if(Frod(3)<0.000001d0) g(3)=-g(3)
            Tref    = Lref/Uref
        else
            stop 'no flow/moving body/gravity'
        endif
    endif

    uMax=maxval([maxval(dabs(uuuIn(1:3))),2.0*pi*MAXVAL(dabs(xyzAmpl(1:3)))*Freq, &
    2.0*pi*MAXVAL(dabs(AoAAmpl(1:3))*[maxval(dabs(xyzful00(:,2))),maxval(dabs(xyzful00(:,1))),maxval(dabs(xyzful00(:,3)))])*Freq])
    
    
!   calculate viscosity, LBM relexation time
    ratio  =  dt/dh
    if(ratio>1.0d0+eps)then
        write(*,*)'dt >  dhmin !!!!!!!!!!'
        write(*,*)'dt <= dhmin (we use streching mesh for LBM)'
        stop
    endif
    
    !for uniform grid, advection length equals grid size
    if(dabs(dt/dh-1.0d0)<eps .and. dabs(dxmax/dh-1.0d0)<eps .and.  dabs(dymax/dh-1.0d0)<eps .and.  dabs(dzmax/dh-1.0d0)<eps)then 
        iStreamModel=1
        write(*,*)'    uniform grid,STLBM'
    else
        iStreamModel=2
        write(*,*)'non-uniform grid,ISLBM'
    endif



    Cs2   =  (1/dsqrt(3.0d0))**2
    nu    =  Uref * Lref/ Re
    Mu    =  nu*denIn
    tau   =  nu/(dt*Cs2)+0.5d0
    Omega =  1.0d0 / tau
    
    Aref=Uref/Tref
    Fref=0.5*denIn*Uref**2*Asfac
    Eref=denIn*Uref**2*Lref**2*Lref
    Pref=denIn*Uref**2*Lref**2*Uref
    !calculate material parameters
    nt=ele(1,4)  
    if(iKB==0)then
        prop(1:nMT,1) = EmR*denIn*Uref**2
        prop(1:nMT,2) = prop(1:nMT,1)/2.0d0/(1.0+psR)
        Lthck= tcR*Lref
        prop(1:nMT,3) = tcR*Lref
        prop(1:nMT,4) = denR*Lref*denIn/prop(1:nMT,3)        
        if    (nt==2)then   !frame
        prop(1:nMT,7) = prop(1:nMT,3)**3/12.0d0
        prop(1:nMT,8) = prop(1:nMT,3)**3/12.0d0
        elseif(nt==3)then   !plate
        prop(1:nMT,6) = prop(1:nMT,3)**3/12.0d0
        else
        endif
        KB=prop(nMT,1)*prop(nMT,6)/(denIn*Uref**2*Lref**3)
        KS=prop(nMT,1)*prop(nMT,3)/(denIn*Uref**2*Lref)
    endif
    
    if(iKB==1)then
        prop(1:nMT,3) = dsqrt(KB/KS*12.0d0)*Lref
        prop(1:nMT,4) = denR*Lref*denIn/prop(1:nMT,3)
        if    (nt==2)then   !frame
        prop(1:nMT,1) = KS*denIn*Uref**2*Lref/prop(1:nMT,3)
        prop(1:nMT,2) = prop(1:nMT,1)/2.0d0/(1.0d0+psR)
        prop(1:nMT,7) = prop(1:nMT,3)**3/12.0d0
        prop(1:nMT,8) = prop(1:nMT,3)**3/12.0d0
        elseif(nt==3)then   !plate
        prop(1:nMT,1) = KS*denIn*Uref**2*Lref/prop(1:nMT,3)
        prop(1:nMT,2) = prop(1:nMT,1)/2.0d0/(1.0d0+psR)
        prop(1:nMT,6) = prop(1:nMT,3)**3/12.0d0
        else
        endif
        EmR = prop(nMT,1)/(denIn*Uref**2) 
        tcR = prop(nMT,3)/Lref
        Lthck=prop(nMT,3)
    endif


    END SUBROUTINE
     
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    calculate multiple relexasion time model parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE calculate_MRTM_params()
    USE simParam
    implicit none
!   ===============================================================================================
    integer:: I 
    real(8):: M_MRT(0:lbmDim,0:lbmDim),M_MRTI(0:lbmDim,0:lbmDim),M(0:lbmDim,0:lbmDim)  
    real(8):: S_D(0:lbmDim,0:lbmDim),S(0:lbmDim)   
!   ======================================================= 
!   calculate MRTM transformation matrix
    DO    I=0,lbmDim
        M_MRT(0,I)=1
        M_MRT(1,I)=19*SUM(ee(I,1:3)**2)-30
        M_MRT(2,I)=(21*SUM(ee(I,1:3)**2)**2-53*SUM(ee(I,1:3)**2)+24)/2.0

        M_MRT(3,I)=ee(I,1)
        M_MRT(5,I)=ee(I,2)
        M_MRT(7,I)=ee(I,3)

        M_MRT(4,I)=(5*SUM(ee(I,1:3)**2)-9)*ee(I,1)
        M_MRT(6,I)=(5*SUM(ee(I,1:3)**2)-9)*ee(I,2)
        M_MRT(8,I)=(5*SUM(ee(I,1:3)**2)-9)*ee(I,3)

        M_MRT(9,I)=3*ee(I,1)**2-SUM(ee(I,1:3)**2)
        M_MRT(11,I)=ee(I,2)**2-ee(I,3)**2

        M_MRT(13,I)=ee(I,1)*ee(I,2)
        M_MRT(14,I)=ee(I,2)*ee(I,3)
        M_MRT(15,I)=ee(I,3)*ee(I,1)

        M_MRT(10,I)=(3*SUM(ee(I,1:3)**2)-5)*(3*ee(I,1)**2-SUM(ee(I,1:3)**2))
        M_MRT(12,I)=(3*SUM(ee(I,1:3)**2)-5)*(ee(I,2)**2-ee(I,3)**2)

        M_MRT(16,I)=(ee(I,2)**2-ee(I,3)**2)*ee(I,1)
        M_MRT(17,I)=(ee(I,3)**2-ee(I,1)**2)*ee(I,2)
        M_MRT(18,I)=(ee(I,1)**2-ee(I,2)**2)*ee(I,3)
    ENDDO
!   calculate the inverse matrix
    M_MRTI=TRANSPOSE(M_MRT)    
    M=MATMUL(M_MRT,M_MRTI)   
    DO    I=0,lbmDim
        M_MRTI(0:lbmDim,I)=M_MRTI(0:lbmDim,I)/M(I,I)
    ENDDO

!   ----------------------------------------------------------
    !S(0:lbmDim)=Omega ! restore to SRT if S is Omega
    !              0   1  2  3  4  5  6  7  8  9     10  11    12  13    14    15    16  17  18
    S(0:lbmDim)=[  s0,s1,s2,s0,s4,s0,s4,s0,s4,Omega,s10,Omega,s10,Omega,Omega,Omega,s16,s16,s16]
    !=====================
!   calculate MRTM collision matrix
    !IM*S*M
    S_D(0:lbmDim,0:lbmDim)=0.0D0
    DO    i=0,lbmDim
        S_D(i,i)=S(i)
    ENDDO    
    M_COLLID=MATMUL(MATMUL(M_MRTI,S_D),M_MRT)
    !=====================
!   calculate MRTM body-force matrix
    !IM*(I-0.5D0*S)*M=I-0.5*IM*S*M
    S_D(0:lbmDim,0:lbmDim)=0.0D0
    DO    i=0,lbmDim
        S_D(i,i)=1.0d0
    ENDDO
    M_FORCE=S_D-0.5*M_COLLID

    END SUBROUTINE



!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write check point file for restarting simulation
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_checkpoint_file()
    USE simParam
    IMPLICIT NONE
    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='replace')
    write(13) step,time
    write(13) fIn,xGrid,yGrid,zGrid,image
    write(13) IXref,IYref,IZref,NDref 
    write(13) xyzful0,xyzful,xyzfulIB,dspful,velful,accful,extful,mss,mssful,grav
    write(13) triad_nn,triad_ee,triad_e0
    write(13) triad_n1,triad_n2,triad_n3
    write(13) UPre,UNow,Et,Ek,Ep,Es,Eb,Ew 
    close(13)
    ENDSUBROUTINE


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read check point file for restarting simulation
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_checkpoint_file()
    USE simParam
    IMPLICIT NONE
    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='old')
    read(13) step,time 
    read(13) fIn,xGrid,yGrid,zGrid,image
    read(13) IXref,IYref,IZref,NDref  
    read(13) xyzful0,xyzful,xyzfulIB,dspful,velful,accful,extful,mss,mssful,grav
    read(13) triad_nn,triad_ee,triad_e0
    read(13) triad_n1,triad_n2,triad_n3
    read(13) UPre,UNow,Et,Ek,Ep,Es,Eb,Ew 
    close(13)
    ENDSUBROUTINE