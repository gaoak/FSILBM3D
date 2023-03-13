!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE redFile()
    USE simParam
    implicit none
    integer, allocatable:: Fishnum(:),SampBodyNodeKind(:)
    integer:: i,temp,Fishkind,iFish,itemp,iFishkind,FishnumBegin,FishnumLast,stepnum,iBodyModelKind
    real(8):: denRKind,psRKind,EmRKind,tcRKind,KBKind,KSKind,isMotionGivenKind(6)
    real(8):: waittingTimeKind(2),XYZointial(3),dXYZo(3),XYZAmplKind(3),XYZPhiKind(3),AoAoKind(3),AoAAmplKind(3),AoAPhiKind(3)
    character(LEN=40):: FEmeshNameKind
    open(unit=111,file='inFlow.dat') 
    call readequal(111)     !======================================
    read(111,*)     npsize
    read(111,*)     isRelease
    read(111,*)     isConCmpt,  iCollidModel,  isVelInit
    read(111,*)     RefVelocity , RefLength , isKB
    read(111,*)     timeSimTotl,timeOutTemp
    read(111,*)     timeOutFlow,timeOutBody,timeOutInfo
    read(111,*)     timeOutFlEd,timeOutFlBg
    read(111,*)     Palpha, Pbeta, Pramp
    call readequal(111)     !======================================
    read(111,*)     uIn(1:SpcDim)
    read(111,*)     shearRateIn(1:SpcDim)
    read(111,*)     boundaryConditions(1:2),VelocityKind
    if(VelocityKind==1) then
        ParabolicParameter(1:2) = shearRateIn(1:2)
    elseif(VelocityKind==2) then
        VelocityAmp = shearRateIn(1)
        VelocityFreq = shearRateIn(2)
    endif
    read(111,*)     boundaryConditions(3),MovingKind1,MovingVel1,MovingFreq1   
    read(111,*)     boundaryConditions(4),MovingKind2,MovingVel2,MovingFreq2          
    call readequal(111)     !====================================== 
    read(111,*)     Re   
    read(111,*)     dt
    read(111,*)     Frod(1:3)
    call readequal(111)     !====================================== 
    read(111,*)     iBC
    read(111,*)     LBmeshName
    read(111,*)     denIn    
    read(111,*)     dtolLBM,ntolLBM      
    call readequal(111)     !======================================
    read(111,*)     isMoveGrid,     iMoveDim,      isMoveOutput
    read(111,*)     Xref,      Yref,   Zref
    read(111,*)     AmplInitDist(1:2)
    read(111,*)     FreqInitDist(1:2)
    call readequal(111)     !======================================
    read(111,*)     nFish,Fishkind
    call readequal(111)     !======================================
    
    allocate( FEmeshName(1:nFish) )
    allocate( iBodyModel(1:nFish) )
    allocate( isMotionGiven(1:nFish,1:DOFDim) )

    allocate( Fishnum(1:Fishkind+1) )
    Fishnum(1)=1
    FishnumBegin=0
    FishnumLast=0

    allocate( denR(1:nFish), tcR(1:nFish), EmR(1:nFish), psR(1:nFish), KB(1:nFish), KS(1:nFish) )

    do iFishkind=1,Fishkind
    read(111,*)     Fishnum(iFishkind+1), FEmeshNameKind  ! the length of array which records each kind fish series numbers is iFishkind+1
    read(111,*)     iBodyModelKind  !1 for rigid body; 2 for flexible body
    read(111,*)     isMotionGivenKind(1:3)
    read(111,*)     isMotionGivenKind(4:6)
    read(111,*)     denRKind, psRKind
    if(isKB==0) read(111,*)     EmRKind, tcRKind 
    if(isKB==1) read(111,*)     KBKind, KSKind

    FishnumBegin=FishnumBegin+Fishnum(iFishkind)
    FishnumLast=FishnumLast+Fishnum(iFishkind+1)

    do iFish=FishnumBegin,FishnumLast
    itemp = iFish
    FEmeshName(iFish)=FEmeshNameKind
    iBodyModel(iFish)=iBodyModelKind
    isMotionGiven(iFish,1:6)=isMotionGivenKind(1:6)
    denR(iFish)=denRKind
    psR(iFish)=psRKind

    if(isKB==0) then
    EmR(iFish)=EmRKind
    tcR(iFish)=tcRKind
    elseif(isKB==1) then
    KB(iFish)=KBKind
    KS(iFish)=KSKind
    endif

    enddo
    enddo
    
    call readequal(111)     !======================================
    read(111,*)     Freq,   St
    call readequal(111)     !======================================
    allocate( waittingTime(1:nFish,1:2) )
    allocate( XYZo(1:nFish,1:3),XYZAmpl(1:nFish,1:3),XYZPhi(1:nFish,1:3) )
    allocate( AoAo(1:nFish,1:3),AoAAmpl(1:nFish,1:3),AoAPhi(1:nFish,1:3) )
    
    FishnumBegin=0
    FishnumLast=0

    do iFishkind=1,Fishkind
    read(111,*)     itemp
    read(111,*)     waittingTimeKind(1:2)
    read(111,*)     XYZointial(1:3)
    read(111,*)     dXYZo(1:3)
    read(111,*)     XYZAmplKind(1:3)
    read(111,*)     XYZPhiKind(1:3)
    read(111,*)     AoAoKind(1:3)
    read(111,*)     AoAAmplKind(1:3)
    read(111,*)     AoAPhiKind(1:3)

    FishnumBegin=FishnumBegin+Fishnum(iFishkind)
    FishnumLast=FishnumLast+Fishnum(iFishkind+1)
    stepnum=0

    do iFish=FishnumBegin,FishnumLast
    waittingTime(iFish,1:2)=waittingTimeKind(1:2)
    XYZo(iFish,1:3)=XYZointial(1:3)+dXYZo(1:3)*stepnum
    XYZAmpl(iFish,1:3)=XYZAmplKind(1:3)
    AoAo(iFish,1:3)=AoAoKind(1:3)
    AoAAmpl(iFish,1:3)=AoAAmplKind(1:3)
    AoAPhi(iFish,1:3)=AoAPhiKind(1:3)
    stepnum=stepnum+1
    enddo

    enddo

    call readequal(111)     !======================================
    read(111,*)     numsubstep
    read(111,*)     dampK,dampM
    read(111,*)     Newmarkdelta,Newmarkalpha
    read(111,*)     alphaf,alpham,alphap
    read(111,*)     dtolFEM,ntolFEM 
    
    call readequal(111)     !======================================
    read(111,*)     numSampFlow
    allocate(SampFlowPint(1:numSampFlow,1:2))
    do  i=1,     numSampFlow
        read(111,*) SampFlowPint(i,1:2)
    enddo
    
    call readequal(111)     !======================================
    read(111,*)     numSampBody
    allocate(SampBodyNodeKind(1:numSampBody))
    allocate(SampBodyNode(1:nFish,1:numSampBody))
    read(111,*)     SampBodyNodeKind(1:numSampBody)
    do iFish=1,nFish
    SampBodyNode(iFish,1:numSampBody)=SampBodyNodeKind(1:numSampBody)
    enddo
    call readequal(111)     !======================================
    close(111)
    END SUBROUTINE

    SUBROUTINE readequal(ifile)
    implicit none
    integer ifile
    character (40):: linestring
    linestring = 'null'
    do while(linestring(1:1)/='=')
        read(ifile, *) linestring
    enddo
    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    SUBROUTINE memFlow()
    USE simParam
    implicit none
    integer,parameter::idFile=111
    integer:: x,y,temp

!   read fluid mesh**********************************************************************************
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
    close(idFile)
!   Define variables**********************************************************************************
    allocate(image(yDim,xDim))
    allocate(fIn(1:yDim,1:xDim,0:lbmDim),fInTemp(1:yDim,1:xDim))
    allocate(uuu(1:yDim,1:xDim,1:SpcDim),force(1:yDim,1:xDim,1:SpcDim))
    allocate(prs(1:yDim,1:xDim),den(1:yDim,1:xDim),vor(1:yDim,1:xDim))

!   compute grid size********************************************************************************* 
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

    dxmin=dabs(minval(dx(1:xDim)))
    dymin=dabs(minval(dy(1:yDim)))
    dxmax=dabs(maxval(dx(1:xDim)))
    dymax=dabs(maxval(dy(1:yDim)))

    cptxmin=xGrid0(1)
    cptxmax=xGrid0(xDim)
    cptymin=yGrid0(1)
    cptymax=yGrid0(yDim)

    if(dabs(dxmin-dymin)>eps) then
    write(*,*)'dxmin/=dymin'
    stop
    endif
    dh=dxmin

    END SUBROUTINE


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE memBody()
    USE simParam
    implicit none
    integer:: iND,iEL,jND,iFish,itemp
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az,lll
    
    write(*,'(A)') '=============================================================================='
    
    allocate( nND(1:nFish),nEL(1:nFish),nMT(1:nFish),nEQ(1:nFish),nBD(1:nFish),maxstiff(1:nFish)  )
    !open(idat,file='BMeshNum.dat')
    !do iFish=1,nFish
    !        read(idat,*) itemp,nND(iFish),nEL(iFish),nMT(iFish) 
    !enddo
    !close(idat)
    do iFish=1,nFish
    open(unit=idat, file = trim(adjustl(FEmeshName(iFish))))
    rewind(idat)
    read(idat,*)  
    read(idat,*)nND(iFish),nEL(iFish),nMT(iFish)
    read(idat,*)  
    close(idat)
    enddo

    nND_max=maxval(nND(:))
    nEL_max=maxval(nEL(:))
    nMT_max=maxval(nMT(:))
    nEQ_max=nND_max*6
    nEQ(:)=nND(:)*6
    
    nEL_all = sum(nEL(:)) 
    nND_all = sum(nND(:))                                      
    
    allocate( ele_all(nEL_all,5),xyzful_all(nND_all,6),velful_all(nND_all,6),xyzfulIB_all(nND_all,6),extful_all(nND_all,6) )
    
    allocate( ele(1:nFish,nEL_max,5),nloc(1:nFish,nND_max*6),nprof(1:nFish,nND_max*6),nprof2(1:nFish,nND_max*6),jBC(1:nFish,nND_max,6))
    allocate( xyzful00(1:nFish,nND_max,6),xyzful0(1:nFish,nND_max,6),mssful(1:nFish,nND_max,6),lodful(1:nFish,nND_max,6), &
              extful(1:nFish,nND_max,6),extful1(1:nFish,nND_max,6),extful2(1:nFish,nND_max,6),grav(1:nFish,nND_max,6),repful(1:nFish,nND_max,6),    &
              vBC(1:nFish,nND_max,6),mss(1:nFish,nND_max*6),prop(1:nFish,nMT_max,10))

    allocate( xyzful(1:nFish,nND_max,6),xyzfulIB(1:nFish,nND_max,6),xyzfulnxt(1:nFish,nND_max,6),dspful(1:nFish,nND_max,6),velful(1:nFish,nND_max,6),accful(1:nFish,nND_max,6)) 
    allocate( triad_nn(1:nFish,3,3,nND_max),triad_ee(1:nFish,3,3,nEL_max),triad_e0(1:nFish,3,3,nEL_max) )
    allocate( triad_n1(1:nFish,3,3,nEL_max),triad_n2(1:nFish,3,3,nEL_max),triad_n3(1:nFish,3,3,nEL_max) )
    
    extful(:,:,1:6)=0.d0
    extful1(:,:,1:6)=0.d0
    extful2(:,:,1:6)=0.d0
    repful(:,:,:)=0.d0
     
!   ===============================================================================================  
    do iFish=1,nFish
    open(unit=idat, file = trim(adjustl(FEmeshName(iFish))))
    rewind(idat)
    read(idat,*)  
    read(idat,*)nND(iFish),nEL(iFish),nMT(iFish)
    read(idat,*)             
    !write(*,*)'read FEMeshFile ',iFish      
!   ===============================================================================================
    call readdt(jBC(iFish,1:nND(iFish),1:6),ele(iFish,1:nEL(iFish),1:5),nloc(iFish,1:nND(iFish)*6),nprof(iFish,1:nND(iFish)*6), &
                nprof2(iFish,1:nND(iFish)*6),xyzful00(iFish,1:nND(iFish),1:6),prop(iFish,1:nMT(iFish),1:10), &
                nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),nBD(iFish),maxstiff(iFish),idat)
    close(idat)

    write(*,*)'read FEMeshFile ',iFish,' end' 
    enddo
    
    !package
    do iFish=1,nFish
        if(iFish.eq.1)then
           do iEL=1,nEL(iFish)
              ele_all(iEL,1:3)=ele(iFish,iEL,1:3) 
              ele_all(iEL,4:5)=ele(iFish,iEL,4:5) 
           enddo
        else
        do iEL=1,nEL(iFish)
            ele_all(iEL+sum(nEL(1:iFish-1)),1:3)=ele(iFish,iEL,1:3)+ sum(nND(1:iFish-1))
            ele_all(iEL+sum(nEL(1:iFish-1)),4:5)=ele(iFish,iEL,4:5) 
        enddo 
        endif
    enddo
    
    allocate(Lchod(1:nFish)) 
    Lchod(:)=0.0d0
    
    !iFish=1 !maximum body diameter
    do  iFish=1,nFish
    do    iND=1,nND(iFish)
    do    jND=iND+1,nND(iFish)
        lll=dsqrt(sum((xyzful00(iFish,iND,1:2)-xyzful00(iFish,jND,1:2))**2)) !two points distance
        if(lll>Lchod(iFish))then
            Lchod(iFish)=lll
        endif
    enddo        
    enddo
    if(dabs(Lchod(iFish)-1.0d0)<=1.0d-2) Lchod(iFish)=1.0d0
    enddo
    
!    Lleng=0.0d0
!    do    iEL=1,nEL
!        Lleng=Lleng+dsqrt(sum((xyzful00(ele(iEL,1),1:2)-xyzful00(ele(iEL,2),1:2))**2))
!    enddo

!   Determine the degree of freedom of each point***************************************************
    do  iFish=1,nFish
    do    iND=1,nND(iFish)
        if(jBC(iFish,iND,1)==1) jBC(iFish,iND,1:6)=isMotionGiven(iFish,1:6)
    enddo
    enddo
   
    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    SUBROUTINE intFlow()
    USE simParam
    implicit none
    integer:: x,y
    real(8):: uSqr,uxy(0:lbmDim),fEq(0:lbmDim)
    real(8):: vel(1:SpcDim)

!   copy fluid mesh***************************************************************************************
    xGrid(1:xDim)=xGrid0(1:xDim)
    yGrid(1:yDim)=yGrid0(1:yDim)
!   compute macro velocity
    do  x = 1, xDim
    do  y = 1, yDim
        if(isVelInit==0) then
            vel(1) = 0.0
            vel(2) = 0.0
        else
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(x), yGrid(y), vel)
            elseif(VelocityKind==1) then
                call evaluateParabolicVelocity(xGrid(x), yGrid(y), vel)
            elseif(VelocityKind==2) then
                vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
                vel(2)=uIn(2)
            endif
        endif
        uuu(y, x, 1) = vel(1)
        uuu(y, x, 2) = vel(2)
    enddo
    enddo
    den(1:yDim,1:xDim)   = denIn
    prs(1:yDim,1:xDim)   = Cs2*(den(1:yDim,1:xDim)-denIn)

    call initDisturb(xDim,yDim,xGrid,yGrid,AmplInitDist,FreqInitDist,uuu)
    !set fIn as equlibrium state
    do  x = 1, xDim
    do  y = 1, yDim
        uSqr          = sum(uuu(y,x,1:spcDim)**2)
        uxy(0:lbmDim) = uuu(y,x,1) * ee(0:lbmDim,1) + uuu(y,x,2) * ee(0:lbmDim,2)
        fEq(0:lbmDim) = wt(0:lbmDim) * den(y,x) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)
        fIn(y,x,0:lbmDim)=fEq(0:lbmDim)
    enddo
    enddo
!   Record boundary conditions and flow field******************************************************
    image          =  fluid
    image(:,1)     =  xMinBC
    image(:,xDim)  =  xMaxBC
    image(1,:)     =  yMinBC
    image(yDim,:)  =  yMaxBC
    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE intBody()
    USE simParam
    implicit none
    integer:: iND,iFish,iCount
    allocate(TTT00(1:nFish,1:3,1:3),TTT0(1:nFish,1:3,1:3),TTTnow(1:nFish,1:3,1:3),TTTnxt(1:nFish,1:3,1:3))
    allocate(XYZ(1:nFish,1:3),XYZd(1:nFish,1:3),UVW(1:nFish,1:3) )
    allocate(AoA(1:nFish,1:3),AoAd(1:nFish,1:3),WWW1(1:nFish,1:3),WWW2(1:nFish,1:3),WWW3(1:nFish,1:3) )
    iCount = 0
    do iFish=1,nFish
     
    TTT00(iFish,:,:)=0.0d0
    TTT00(iFish,1,1)=1.0d0
    TTT00(iFish,2,2)=1.0d0
    TTT00(iFish,3,3)=1.0d0
    
    XYZ(iFish,1:3)=XYZo(iFish,1:3)+XYZAmpl(iFish,1:3)*dcos(2.0*pi*Freq*time+XYZPhi(iFish,1:3))
    AoA(iFish,1:3)=AoAo(iFish,1:3)+AoAAmpl(iFish,1:3)*dcos(2.0*pi*Freq*time+AoAPhi(iFish,1:3))
    
    call AoAtoTTT(AoA(iFish,1:3),TTT0(iFish,1:3,1:3))
    call AoAtoTTT(AoA(iFish,1:3),TTTnow(iFish,1:3,1:3))
    call get_angle_triad(TTT0(iFish,1:3,1:3),TTTnow(iFish,1:3,1:3),AoAd(iFish,1),AoAd(iFish,2),AoAd(iFish,3))
    
    do  iND=1,nND(iFish)
        xyzful0(iFish,iND,1:3)=matmul(TTT0(iFish,1:3,1:3),xyzful00(iFish,iND,1:3))+XYZ(iFish,1:3)
        xyzful0(iFish,iND,4:6)=AoAd(iFish,1:3)
    enddo

    xyzful(iFish,1:nND(iFish),1:6)=xyzful0(iFish,1:nND(iFish),1:6)     
    velful(iFish,1:nND(iFish),1:6)=0.0
                
    dspful(iFish,1:nND(iFish),1:6)=0.0
    accful(iFish,1:nND(iFish),1:6)=0.0

    xyzfulIB(iFish,1:nND(iFish),1:6)=xyzful(iFish,1:nND(iFish),1:6)
 
    do iND=1,nND(iFish)
        xyzful_all(iND+iCount,1:6)=xyzful(iFish,iND,1:6)
        velful_all(iND+iCount,1:6)=velful(iFish,iND,1:6)
        xyzfulIB_all(iND+iCount,1:6)=xyzfulIB(iFish,iND,1:6)
        extful_all(iND+iCount,1:6)=0.d0
    enddo
  

    CALL formmass_D(  ele(iFish,1:nEL(iFish),1:5),xyzful0(iFish,1:nND(iFish),1),xyzful0(iFish,1:nND(iFish),2),xyzful0(iFish,1:nND(iFish),3), &
                     prop(iFish,1:nMT(iFish),1:10),mss(iFish,1:nND(iFish)*6),nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),alphaf,alpham,alphap ) 
    do  iND = 1, nND(iFish)
        mssful(iFish,iND,1:6)=mss(iFish,(iND-1)*6+1:(iND-1)*6+6)
        grav(iFish,iND,1:6) = mssful(iFish,iND,1:6)*[g(1),g(2),g(3),0.0d0,0.0d0,0.0d0]
    enddo
    
    CALL init_triad_D( ele(iFish,1:nEL(iFish),1:5),xyzful(iFish,1:nND(iFish),1),xyzful(iFish,1:nND(iFish),2),xyzful(iFish,1:nND(iFish),3), &
                      triad_nn(iFish,1:3,1:3,1:nND(iFish)),triad_n1(iFish,1:3,1:3,1:nEL(iFish)),triad_n2(iFish,1:3,1:3,1:nEL(iFish)),triad_n3(iFish,1:3,1:3,1:nEL(iFish)), &
                      triad_ee(iFish,1:3,1:3,1:nEL(iFish)),triad_e0(iFish,1:3,1:3,1:nEL(iFish)),nND(iFish),nEL(iFish) ) 
  
    
    iCount = iCount + nND(iFish)
    enddo  !do iFish=1,nFish
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptPara()
    USE simParam
    implicit none
    integer:: iFish
    !angle degrees to radians
    AoAo(:,1:3)=AoAo(:,1:3)/180.0*pi
    AoAAmpl(:,1:3)=AoAAmpl(:,1:3)/180.0*pi
    AoAPhi(:,1:3)=AoAPhi(:,1:3)/180.0*pi
    XYZPhi(:,1:3)=XYZPhi(:,1:3)/180.0*pi
    !use the body with minimum diameter as reference body
    iFish=minloc(abs(Lchod(1:nFish)),1)
    if(iFish.lt.0 .and. iFish.gt.nFish)then
        write(*,*)'bad iFish in cptPara, stop', iFish
        stop
    endif
    
    if(RefLength==0) then
        Lref = Lchod(iFish) !body diameter as the reference length
    else
        Lref = 1.d0
    endif

    if(RefVelocity==0) then
        Uref = dabs(uIn(1))
    elseif(RefVelocity==1) then
        Uref = dabs(uIn(2))
    elseif(RefVelocity==2) then
        Uref = dsqrt(uIn(1)**2 + uIn(2)**2)  
    elseif(RefVelocity==3) then
        Uref = abs(max(MovingVel1,MovingVel2))  !Moving Wall Velocity
    elseif(RefVelocity==4) then
        Uref = abs(VelocityAmp)  !Velocity Amplitude
    elseif(RefVelocity==10) then
        Uref = Lref * Freq
    elseif(RefVelocity==11) then
        Uref = 2.d0*pi*Freq*MAXVAL(dabs(xyzAmpl(1:nFish,1:3)))
    elseif(RefVelocity==12) then
        Uref = 2.d0*pi*Freq*MAXVAL(dabs(xyzAmpl(1:nFish,1:3)))*2.D0 !Park 2017 pof
    else
        Uref = 1.d0
    endif
    Tref = Lref / Uref

    St = Lref*Freq/Uref
    g(1:3)=Frod(1:3)*Uref**2/Lref
    
    uMax=maxval([maxval(dabs(uIn(1:2))),2.0*pi*MAXVAL(dabs(xyzAmpl(:,1:3)))*Freq,2.0*pi*dabs(MAX(MovingVel1*MovingFreq1,MovingVel2*MovingFreq1)/Tref),&
    2.0*pi*MAXVAL(dabs(AoAAmpl(iFish,1:3))*[maxval(dabs(xyzful00(iFish,:,2))),maxval(dabs(xyzful00(iFish,:,1))),maxval(dabs(xyzful00(iFish,:,3)))])*Freq])


    !time step should be smaller than spatial grid size
    ratio  =  dt/dh
    if(ratio>1.0d0+eps)then
        write(*,*)'dt>dh'
        stop
    endif
    
    !test if the grid is uniform
    if(dabs(dt/dh-1.0d0)<eps .and. dabs(dxmax/dh-1.0d0)<eps .and.  dabs(dymax/dh-1.0d0)<eps )then 
        iStreamModel=1
        write(*,*)'uniform grid,STLBM'
    else
        iStreamModel=2
        write(*,*)'non-uniform grid,ISLBM'
    endif

    Cs2     =  (1/dsqrt(3.0d0))**2
    Nu      =  Uref * Lref / Re
    Mu      =  DenIn*Nu
    Tau     =  3.0d0*Nu/dt+0.5d0 
    Omega   =  1.0d0 / Tau

    Aref=Uref/Tref
    Fref=0.5*denIn*Uref**2*Lref
    Eref=denIn*Uref**2*Lref*Lref
    Pref=denIn*Uref**2*Lref*Uref  
    !Calculate the parameters of the bodies
    do iFish=1,nFish
        
    if(isKB==0)then
        prop(iFish,1:nMT(iFish),1) = EmR(iFish)*denIn*Uref**2
        prop(iFish,1:nMT(iFish),2) = prop(iFish,1:nMT(iFish),1)/2.0d0/(1.0+psR(iFish))
        Lthck= tcR(iFish)*Lref
        prop(iFish,1:nMT(iFish),3) = tcR(iFish)*Lref
        prop(iFish,1:nMT(iFish),4) = denR(iFish)*Lref*denIn/prop(iFish,1:nMT(iFish),3)        
        prop(iFish,1:nMT(iFish),7) = prop(iFish,1:nMT(iFish),3)**3/12.0d0
        prop(iFish,1:nMT(iFish),8) = prop(iFish,1:nMT(iFish),3)**3/12.0d0

        KB(iFish)=prop(iFish,nMT(iFish),1)*prop(iFish,nMT(iFish),6)/(denIn*Uref**2*Lref**3)
        KS(iFish)=prop(iFish,nMT(iFish),1)*prop(iFish,nMT(iFish),3)/(denIn*Uref**2*Lref)
    endif
    
    if(isKB==1)then
        prop(iFish,1:nMT(iFish),3) = dsqrt(KB(iFish)/KS(iFish)*12.0)*Lref ! h[eight] of rectangluar-section beam, width b = 1
        prop(iFish,1:nMT(iFish),4) = denR(iFish)*Lref*denIn/prop(iFish,1:nMT(iFish),3) ! real density per volume
        prop(iFish,1:nMT(iFish),1) = KS(iFish)*denIn*Uref**2*Lref/prop(iFish,1:nMT(iFish),3) ! E
        prop(iFish,1:nMT(iFish),2) = prop(iFish,1:nMT(iFish),1)/2.0d0/(1.0d0+psR(iFish)) ! G, shear modulus
        prop(iFish,1:nMT(iFish),7) = prop(iFish,1:nMT(iFish),3)**3/12.0d0 ! I, the second cross-section moment, width = 1
        prop(iFish,1:nMT(iFish),8) = prop(iFish,1:nMT(iFish),3)**3/12.0d0 ! I

        EmR(iFish) = prop(iFish,nMT(iFish),1)/(denIn*Uref**2)      
        tcR(iFish) = prop(iFish,nMT(iFish),3)/Lref
        Lthck=prop(iFish,nMT(iFish),3)
    endif
    
    enddo

    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptMRTM()
    USE simParam
    implicit none
    integer:: I
    real(8):: M_MRT(0:lbmDim,0:lbmDim),M_MRTI(0:lbmDim,0:lbmDim),M(0:lbmDim,0:lbmDim)
    real(8):: S_D(0:lbmDim,0:lbmDim),S(0:lbmDim)

!   ======================================================= 
    DO  I=0,lbmDim
        M_MRT(0,I)= 1.0d0
        M_MRT(1,I)=-4.0d0+3.0d0*SUM(ee(I,1:2)**2)
        M_MRT(2,I)=(9.0d0*SUM(ee(I,1:2)**2)**2-21.0d0*SUM(ee(I,1:2)**2)+8.0d0)/2.0

        M_MRT(3,I)=ee(I,1)
        M_MRT(5,I)=ee(I,2)
        
        M_MRT(4,I)=(3.0d0*SUM(ee(I,1:2)**2)-5.0d0)*ee(I,1)
        M_MRT(6,I)=(3.0d0*SUM(ee(I,1:2)**2)-5.0d0)*ee(I,2)

        M_MRT(7,I)=ee(I,1)**2-ee(I,2)**2
        M_MRT(8,I)=ee(I,1)*ee(I,2)
    enddo
    M_MRTI=TRANSPOSE(M_MRT)    
    M=MATMUL(M_MRT,M_MRTI)   
    DO  I=0,lbmDim
        M_MRTI(0:lbmDim,I)=M_MRTI(0:lbmDim,I)/M(I,I)
    ENDDO
!==============================================================================
    S(0:lbmDim)=[s0,s1,s2,s0,s4,s0,s4,Omega,Omega]
    !IM*S*M
    S_D(0:lbmDim,0:lbmDim)=0.0D0
    DO    i=0,lbmDim
        S_D(i,i)=S(i)
    ENDDO
    !M_COLLID1=MATMUL(M_MRTI,S_D)    
    M_COLLID=MATMUL(MATMUL(M_MRTI,S_D),M_MRT)
    !=====================
    !IM*(I-0.5D0*S)*M=I-0.5*IM*S*M
    S_D(0:lbmDim,0:lbmDim)=0.0D0
    DO    i=0,lbmDim
        S_D(i,i)=1.0d0
    ENDDO
    M_FORCE=S_D-0.5*M_COLLID
    END SUBROUTINE
    
    
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write temp data for continue computation
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wrtCont()
    USE simParam
    IMPLICIT NONE
    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='replace')
    write(13) Step,time
    write(13) fIn,xGrid,yGrid,image
    write(13) IXref,IYref,NDref 
    write(13) xyzful0,xyzful,xyzfulIB,dspful,velful,accful,extful,mss,mssful,grav
    write(13) triad_nn,triad_ee,triad_e0
    write(13) triad_n1,triad_n2,triad_n3
    write(13) UPre,UNow,Et,Ek,Ep,Es,Eb,Ew 
    close(13)
    ENDSUBROUTINE


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     read temp data for continue computation      
!   
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE redCont()
    USE simParam
    IMPLICIT NONE
      
    allocate(TTT00(1:nFish,1:3,1:3),TTT0(1:nFish,1:3,1:3),TTTnow(1:nFish,1:3,1:3),TTTnxt(1:nFish,1:3,1:3))
    allocate(XYZ(1:nFish,1:3),XYZd(1:nFish,1:3),UVW(1:nFish,1:3) )
    allocate(AoA(1:nFish,1:3),AoAd(1:nFish,1:3),WWW1(1:nFish,1:3),WWW2(1:nFish,1:3),WWW3(1:nFish,1:3) )
    
    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='old')
    read(13) Step,time 
    read(13) fIn,xGrid,yGrid,image
    read(13) IXref,IYref,NDref  
    read(13) xyzful0,xyzful,xyzfulIB,dspful,velful,accful,extful,mss,mssful,grav
    read(13) triad_nn,triad_ee,triad_e0
    read(13) triad_n1,triad_n2,triad_n3
    read(13) UPre,UNow,Et,Ek,Ep,Es,Eb,Ew
    close(13)
    ENDSUBROUTINE

