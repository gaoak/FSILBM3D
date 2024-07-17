!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_file()
    USE simParam
    implicit none
    real(8):: iXYZ(1:3),dXYZ(1:3)
    integer:: i,iFish,iKind,FishKind,Order0
    integer:: FishOrder1,FishOrder2,LineX,LineY
    integer, allocatable:: FishNum(:),NumX(:),NumY(:)
    character(LEN=40):: nFEmeshName
    integer:: niBodyModel,nisMotionGiven(1:6)
    real(8):: ndenR,npsR,nEmR,ntcR,nKB,nKS,nFreq,nSt
    real(8):: nXYZAmpl(1:3),nXYZPhi(1:3),nAoAo(1:3),nAoAAmpl(1:3),nAoAPhi(1:3)
    open(unit=111,file='inFlow.dat') 
    read(111,*)     !======================================
    read(111,*)     npsize
    read(111,*)     isRelease
    read(111,*)     isConCmpt,  iCollidModel
    read(111,*)     RefVelocity, iKB
    read(111,*)     timeSimTotl,timeOutTemp
    read(111,*)     timeOutFlow,timeOutBody,timeOutInfo
    read(111,*)     timeOutFlEd,timeOutFlBg
    read(111,*)     Palpha,Pbeta,Pramp 
    read(111,*)     !======================================
    read(111,*)     uuuIn(1:3)
    read(111,*)     shearRateIn(1:3)     
    read(111,*)     boundaryConditions(1:6)
    read(111,*)     !======================================
    read(111,*)     Re     
    read(111,*)     dt      
    read(111,*)     Frod(1:3)            !Gravity
    read(111,*)     !======================================
    read(111,*)     iBC 
    read(111,*)     LBmeshName
    read(111,*)     denIn   
    read(111,*)     dtolLBM,ntolLBM      !velocity iteration
    read(111,*)     !======================================
    read(111,*)     nFish,FishKind
    read(111,*)     iForce2Body,iChordDirection

    allocate(FEmeshName(1:nFish),iBodyModel(1:nFish),isMotionGiven(1:nFish,1:DOFDim))
    allocate(denR(1:nFish),EmR(1:nFish),tcR(1:nFish),psR(1:nFish),KB(1:nFish),KS(1:nFish))
    allocate(FishNum(1:(FishKind+1)),NumX(1:FishKind),NumY(1:FishKind))
    FishNum(1)=1
    FishOrder1=0
    FishOrder2=0

    read(111,*)     !======================================
    do iKind=1,FishKind
        read(111,*)     FishNum(iKind+1),NumX(iKind),NumY(iKind)
        read(111,*)     niBodyModel, nFEmeshName
        read(111,*)     nisMotionGiven(1:3)
        read(111,*)     nisMotionGiven(4:6)
        read(111,*)     ndenR, npsR
        if(iKB==0) read(111,*)     nEmR, ntcR
        if(iKB==1) read(111,*)     nKB, nKS      
        FishOrder1=FishOrder1+FishNum(iKind  )
        FishOrder2=FishOrder2+FishNum(iKind+1)
        do iFish=FishOrder1,FishOrder2
            iBodyModel(iFish)=niBodyModel
            FEmeshName(iFish)=nFEmeshName
            isMotionGiven(iFish,1:6)=nisMotionGiven(1:6)
            denR(iFish)= ndenR
            psR(iFish) = npsR
            if(iKB==0) then
            EmR(iFish) =nEmR
            tcR(iFish) =ntcR
            elseif(iKB==1) then
            KB(iFish)  =nKB
            KS(iFish)  =nKS
            endif
        enddo
    enddo

    read(111,*)     !======================================
    read(111,*)     numsubstep
    read(111,*)     dampK,dampM
    read(111,*)     NewmarkGamma,NewmarkBeta
    read(111,*)     alphaf,alpham,alphap
    read(111,*)     dtolFEM,ntolFEM  
    read(111,*)     !======================================

    allocate(Freq(1:nFish),St(1:nFish))
    allocate(XYZo(1:nFish,1:3),XYZAmpl(1:nFish,1:3),XYZPhi(1:nFish,1:3))
    allocate(AoAo(1:nFish,1:3),AoAAmpl(1:nFish,1:3),AoAPhi(1:nFish,1:3))
    FishOrder1 =0
    FishOrder2 =0

    do iKind=1,FishKind
        read(111,*)     nFreq, nSt
        read(111,*)     iXYZ(1:3)
        read(111,*)     dXYZ(1:3)
        read(111,*)     nXYZAmpl(1:3)
        read(111,*)     nXYZPhi(1:3)
        read(111,*)     nAoAo(1:3)
        read(111,*)     nAoAAmpl(1:3)
        read(111,*)     nAoAPhi(1:3)
        FishOrder1=FishOrder1+FishNum(iKind  )
        FishOrder2=FishOrder2+FishNum(iKind+1)
        do iFish=FishOrder1,FishOrder2
            Freq(iFish)=nFreq
            St(iFish)  =nSt
            XYZAmpl(iFish,1:3)=nXYZAmpl(1:3)
            XYZPhi(iFish,1:3) =nXYZPhi(1:3)
            AoAo(iFish,1:3)   =nAoAo(1:3)
            AoAAmpl(iFish,1:3)=nAoAAmpl(1:3)
            AoAPhi(iFish,1:3) =nAoAPhi(1:3)
            ! initial position distribution
            Order0 = iFish - FishOrder1
            LineX  = mod(Order0,NumX(iKind))
            LineY  = Order0/NumX(iKind)
            XYZo(iFish,1) = iXYZ(1) + dXYZ(1) * LineX
            XYZo(iFish,2) = iXYZ(2) + dXYZ(2) * LineY
            XYZo(iFish,3) = iXYZ(3)
        enddo
    enddo
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

    allocate(SampFlowPint(1:numSampFlow,1:3))

    read(111,*)     numSampFlow

    do i=1,numSampFlow
        read(111,*) SampFlowPint(i,1:3)
    enddo
    read(111,*)     !======================================
    read(111,*)     numSampBody

    allocate(SampBodyNode(1:nFish,1:numSampBody))
    read(111,*)     SampBodyNode(1,1:numSampBody)
    do iFish=1,nFish
    SampBodyNode(iFish,1:numSampBody)=SampBodyNode(1,1:numSampBody)
    enddo
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
    integer:: iEL,iND,iFish,maxN(1)
    real(8):: xCT,yCT,zCT,xl,xr,yl,yr,zl,zr
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az
    allocate(nND(1:nFish),nEL(1:nFish),nMT(1:nFish),nEQ(1:nFish),nBD(1:nFish),nSTF(1:nFish))

    write(*,'(A)') '=============================================================================='
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

!   ===============================================================================================
    allocate( ele_all(nEL_all,5),xyzful_all(nND_all,6),velful_all(nND_all,6),xyzfulIB_all(nND_all,6),extful1_all(nND_all,6),extful2_all(nND_all,6))

    allocate( ele(1:nFish,nEL_max,5),xyzful00(1:nFish,nND_max,6),xyzful0(1:nFish,nND_max,6),mssful(1:nFish,nND_max,6),lodful(1:nFish,nND_max,6), &
              extful(1:nFish,nND_max,6),repful(1:nFish,nND_max,1:6),extful1(1:nFish,nND_max,6),extful2(1:nFish,nND_max,6),nloc(1:nFish,nND_max*6),nprof(1:nFish,nND_max*6), &
              nprof2(1:nFish,nND_max*6),jBC(1:nFish,nND_max,6),streI(1:nFish,nND_max),bendO(1:nFish,nND_max))
    allocate( grav(1:nFish,nND_max,6),vBC(1:nFish,nND_max,6),mss(1:nFish,nND_max*6),prop(1:nFish,nMT_max,10),areaElem00(1:nFish,nEL_max),areaElem(1:nFish,nEL_max))

    allocate( xyzful(1:nFish,nND_max,6),xyzfulIB(1:nFish,nND_max,6),xyzfulnxt(1:nFish,nND_max,6),dspful(1:nFish,nND_max,6),velful(1:nFish,nND_max,6),accful(1:nFish,nND_max,6)) 
    allocate( triad_nn(1:nFish,3,3,nND_max),triad_ee(1:nFish,3,3,nEL_max),triad_e0(1:nFish,3,3,nEL_max) )
    allocate( triad_n1(1:nFish,3,3,nEL_max),triad_n2(1:nFish,3,3,nEL_max),triad_n3(1:nFish,3,3,nEL_max) ) 

    repful(:,:,1:6) =0.d0
    extful1(:,:,1:6)=0.d0
    extful2(:,:,1:6)=0.d0

!   ===============================================================================================     
    do iFish=1,nFish
    open(unit=idat, file = trim(adjustl(FEmeshName(iFish))))
    rewind(idat)
    read(idat,*)  
    read(idat,*)nND(iFish),nEL(iFish),nMT(iFish)
    read(idat,*) 

    call read_structural_datafile(jBC(iFish,1:nND(iFish),1:6),ele(iFish,1:nEL(iFish),1:5),nloc(iFish,1:nND(iFish)*6),nprof(iFish,1:nND(iFish)*6), &
                                      nprof2(iFish,1:nND(iFish)*6),xyzful00(iFish,1:nND(iFish),1:6),prop(iFish,1:nMT(iFish),1:10),nND(iFish), &
                                      nEL(iFish),nEQ(iFish),nMT(iFish),nBD(iFish),nSTF(iFish),idat)
    close(idat)

    write(*,*)'read FEMeshFile ',iFish,' end' 
    enddo   
!    ===============================================================================================
!    package
    do iFish=1,nFish
        if(iFish.eq.1)then
            do iEL=1,nEL(iFish)
                ele_all(iEL,1:3)=ele(iFish,iEL,1:3) 
                ele_all(iEL,4:5)=ele(iFish,iEL,4:5) 
            enddo
        else
            do iEL=1,nEL(iFish)
                ele_all(iEL+sum(nEL(1:iFish-1)),1:3)=ele(iFish,iEL,1:3)+sum(nND(1:iFish-1))
                ele_all(iEL+sum(nEL(1:iFish-1)),4:5)=ele(iFish,iEL,4:5) 
            enddo 
        endif
    enddo  
!   ===============================================================================================
!   calculate area
    allocate(NDtl(1:nFish,1:5),NDhd(1:nFish,1:3),NDct(1:nFish),elmax(1:nFish),elmin(1:nFish))
    allocate(nAsfac(1:nFish),nLchod(1:nFish),nLspan(1:nFish))
    do iFish=1,nFish
        call cptArea(areaElem00(iFish,1:nEL(iFish)),nND(iFish),nEL(iFish),ele(iFish,1:nEL(iFish),1:5),xyzful00(iFish,1:nND(iFish),1:6))
        nAsfac(iFish)=sum(areaElem00(iFish,1:nEL(iFish)))
        elmax(iFish)=maxval(areaElem00(iFish,1:nEL(iFish)))
        elmin(iFish)=minval(areaElem00(iFish,1:nEL(iFish)))
        streI(iFish,1:nND(iFish))=0.0d0
        bendO(iFish,1:nND(iFish))=0.0d0
    enddo
    !calculate spanwise length, chord length, aspect ratio
    if(iChordDirection==1)then     
        do iFish=1,nFish
        nLchod(iFish) = maxval(xyzful00(iFish,:,1))-minval(xyzful00(iFish,:,1))
        nLspan(iFish) = maxval(xyzful00(iFish,:,2))-minval(xyzful00(iFish,:,2))
        enddo
    elseif(iChordDirection==2)then
        do iFish=1,nFish
        nLchod(iFish) = maxval(xyzful00(iFish,:,2))-minval(xyzful00(iFish,:,2))
        nLspan(iFish) = maxval(xyzful00(iFish,:,1))-minval(xyzful00(iFish,:,1))  
        enddo      
    else
        stop 'no define chordwise'
    endif
    !Use the plate with the largest area as the reference plate
    maxN  = maxloc(nAsfac)
    Asfac = nAsfac(maxN(1))
    Lchod = nLchod(maxN(1))
    Lspan = nLspan(maxN(1))

    if((Lchod-1.0d0)<=1.0d-2)Lchod=1.0d0
    if((Lspan-1.0d0)<=1.0d-2)Lspan=1.0d0

    AR    = Lspan**2/Asfac
!    Lchod=1.0d0 
!    Lspan=1.0d0

    !calculate central position, central node
    do iFish=1,nFish
    xCT=sum(xyzful00(iFish,:,1))/nND(iFish)
    yCT=sum(xyzful00(iFish,:,2))/nND(iFish)
    zCT=sum(xyzful00(iFish,:,3))/nND(iFish)

    !calculate leading-edge central node, three trailing nodes, the nearest neighbouring node
    xl=minval(xyzful00(iFish,:,1))
    xr=maxval(xyzful00(iFish,:,1))
    yl=minval(xyzful00(iFish,:,2))
    yr=maxval(xyzful00(iFish,:,2))
    zl=minval(xyzful00(iFish,:,3))
    zr=maxval(xyzful00(iFish,:,3))
    ! center point
    NDct(iFish)   =minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-xCT)**2+(xyzful00(iFish,1:nND(iFish),2)-yCT)**2),1)
    ! edge points
    NDtl(iFish,1)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-xr )**2+(xyzful00(iFish,1:nND(iFish),2)-         yl)**2),1)
    NDtl(iFish,2)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-xr )**2+(xyzful00(iFish,1:nND(iFish),2)-0.5*(yl+yr))**2),1)
    NDtl(iFish,3)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-xr )**2+(xyzful00(iFish,1:nND(iFish),2)-         yr)**2),1)
    NDtl(iFish,4)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-0.5*(xl+xr) )**2+(xyzful00(iFish,1:nND(iFish),2)-yl)**2),1)
    NDtl(iFish,5)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-0.5*(xl+xr) )**2+(xyzful00(iFish,1:nND(iFish),2)-yr)**2),1)
    
    NDhd(iFish,1)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-xl )**2+(xyzful00(iFish,1:nND(iFish),2)-         yl)**2),1)
    NDhd(iFish,2)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-xl )**2+(xyzful00(iFish,1:nND(iFish),2)-0.5*(yl+yr))**2),1)
    NDhd(iFish,3)=minloc(dsqrt((xyzful00(iFish,1:nND(iFish),1)-xl )**2+(xyzful00(iFish,1:nND(iFish),2)-         yr)**2),1)

    ! if(trim(adjustl(FEmeshName))=="fi30.dat")then 
    !     NDtl(1)=513
    !     NDtl(2)=221
    !     NDtl(3)=536
    !     NDhd(1)=2
    !     NDhd(2)=108
    !     NDhd(3)=392
    ! elseif(trim(adjustl(FEmeshName))=="fi45.dat")then
    !     NDtl(1)=502
    !     NDtl(2)=309
    !     NDtl(3)=530
    !     NDhd(1)=2
    !     NDhd(2)=109
    !     NDhd(3)=417
    ! elseif(trim(adjustl(FEmeshName))=="fi60.dat")then
    !     NDtl(1)=472
    !     NDtl(2)=342
    !     NDtl(3)=524
    !     NDhd(1)=2
    !     NDhd(2)=112
    !     NDhd(3)=417
    ! elseif(trim(adjustl(FEmeshName))=="fi75.dat")then
    !     NDtl(1)=469
    !     NDtl(2)=445
    !     NDtl(3)=533
    !     NDhd(1)=2
    !     NDhd(2)=110
    !     NDhd(3)=421
    ! elseif(trim(adjustl(FEmeshName))=="fi90.dat")then
    ! NDtl(1)=minloc(dsqrt((xyzful00(1:nND,1)-xr )**2+(xyzful00(1:nND,2)-         yl)**2),1)
    ! NDtl(2)=minloc(dsqrt((xyzful00(1:nND,1)-xr )**2+(xyzful00(1:nND,2)-0.5*(yl+yr))**2),1)
    ! NDtl(3)=minloc(dsqrt((xyzful00(1:nND,1)-xr )**2+(xyzful00(1:nND,2)-         yr)**2),1)
    
    ! NDtl(4)=minloc(dsqrt((xyzful00(1:nND,1)-0.5*(xl+xr) )**2+(xyzful00(1:nND,2)-yl)**2),1)
    ! NDtl(5)=minloc(dsqrt((xyzful00(1:nND,1)-0.5*(xl+xr) )**2+(xyzful00(1:nND,2)-yr)**2),1)
    
    ! NDhd(1)=minloc(dsqrt((xyzful00(1:nND,1)-xl )**2+(xyzful00(1:nND,2)-         yl)**2),1)
    ! NDhd(2)=minloc(dsqrt((xyzful00(1:nND,1)-xl )**2+(xyzful00(1:nND,2)-0.5*(yl+yr))**2),1)
    ! NDhd(3)=minloc(dsqrt((xyzful00(1:nND,1)-xl )**2+(xyzful00(1:nND,2)-         yr)**2),1)
    ! elseif(trim(adjustl(FEmeshName))=="fi105.dat")then
    !     NDtl(1)=369
    !     NDtl(2)=491
    !     NDtl(3)=530
    !     NDhd(1)=2
    !     NDhd(2)=109
    !     NDhd(3)=425
    ! elseif(trim(adjustl(FEmeshName))=="fi120.dat")then
    !     NDtl(1)=330
    !     NDtl(2)=523
    !     NDtl(3)=534
    !     NDhd(1)=2
    !     NDhd(2)=112
    !     NDhd(3)=422
    ! elseif(trim(adjustl(FEmeshName))=="fi135.dat")then
    !     NDtl(1)=237
    !     NDtl(2)=528
    !     NDtl(3)=509
    !     NDhd(1)=2
    !     NDhd(2)=109
    !     NDhd(3)=414
    ! elseif(trim(adjustl(FEmeshName))=="fi150.dat")then
    !     NDtl(1)=150
    !     NDtl(2)=515
    !     NDtl(3)=441
    !     NDhd(1)=2
    !     NDhd(2)=112
    !     NDhd(3)=400
    ! endif
    enddo
!   loading boundary type*******************************************************************************
    do  iFish=1,nFish
    do    iND=1,nND(iFish)
        if(jBC(iFish,iND,1)==1) jBC(iFish,iND,1:6)=isMotionGiven(iFish,1:6)
    enddo
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
    real(8):: vel(1:SpcDim)
    integer:: x, y, z

!   grid coordinate***************************************************************************************         
    xGrid(1:xDim)=xGrid0(1:xDim)
    yGrid(1:yDim)=yGrid0(1:yDim)
    zGrid(1:zDim)=zGrid0(1:zDim)
!   macro quantities***************************************************************************************    
    do  x = 1, xDim
    do  y = 1, yDim
    do  z = 1, zDim
        call evaluateShearVelocity(xGrid(x),yGrid(y),zGrid(z), vel)
        uuu(z, y, x, 1) = vel(1)
        uuu(z, y, x, 2) = vel(2)
        uuu(z, y, x, 3) = vel(3)
    enddo
    enddo
    enddo
    !uuu(1:zDim,1:yDim,1:xDim,1) = uuuIn(1)
    !uuu(1:zDim,1:yDim,1:xDim,2) = uuuIn(2)
    !uuu(1:zDim,1:yDim,1:xDim,3) = uuuIn(3)
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

        do iND=1,nND(iFish)
            xyzful0(iFish,iND,1:3)=matmul(TTT0(iFish,1:3,1:3),xyzful00(iFish,iND,1:3))+XYZ(iFish,1:3)
            xyzful0(iFish,iND,4:6)=AoAd(iFish,1:3)
        enddo

        xyzful(iFish,1:nND(iFish),1:6)=xyzful0(iFish,1:nND(iFish),1:6)     
        velful(iFish,1:nND(iFish),1:6)=0.0
                    
        dspful(iFish,1:nND(iFish),1:6)=0.0
        accful(iFish,1:nND(iFish),1:6)=0.0

        xyzfulIB(iFish,1:nND(iFish),1:6)=xyzful(iFish,1:nND(iFish),1:6)

        do iND=1,nND(iFish)
            xyzful_all(iND+iCount,1:6)   =xyzful(iFish,iND,1:6)
            velful_all(iND+iCount,1:6)   =velful(iFish,iND,1:6)
            xyzfulIB_all(iND+iCount,1:6) =xyzfulIB(iFish,iND,1:6)
            extful1_all(iND+iCount,1:6)  =0.d0
            extful2_all(iND+iCount,1:6)  =0.d0
        enddo

        CALL formmass_D(ele(iFish,1:nEL(iFish),1:5),xyzful0(iFish,1:nND(iFish),1),xyzful0(iFish,1:nND(iFish),2),xyzful0(iFish,1:nND(iFish),3), &
                        prop(iFish,1:nMT(iFish),1:10),mss(iFish,1:nND(iFish)*6),nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),alphaf,alpham,alphap) 
        
        do iND = 1, nND(iFish)
            mssful(iFish,iND,1:6)= mss(iFish,(iND-1)*6+1:(iND-1)*6+6)
            grav(iFish,iND,1:6)  = mssful(iFish,iND,1:6)*[g(1),g(2),g(3),0.0d0,0.0d0,0.0d0]
        enddo

        CALL init_triad_D(ele(iFish,1:nEL(iFish),1:5),xyzful(iFish,1:nND(iFish),1),xyzful(iFish,1:nND(iFish),2),xyzful(iFish,1:nND(iFish),3), &
                          triad_nn(iFish,1:3,1:3,1:nND(iFish)),triad_n1(iFish,1:3,1:3,1:nEL(iFish)),triad_n2(iFish,1:3,1:3,1:nEL(iFish)),triad_n3(iFish,1:3,1:3,1:nEL(iFish)), &
                          triad_ee(iFish,1:3,1:3,1:nEL(iFish)),triad_e0(iFish,1:3,1:3,1:nEL(iFish)),nND(iFish),nEL(iFish)) 
        
        iCount = iCount + nND(iFish)
    enddo
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    calculate Bolztman parameters from flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE calculate_LB_params()
    USE simParam
    implicit none
    integer:: nt(1:nFish),iFish
    real(8):: nUref(1:nFish)
!   angle to radian
    AoAo(:,1:3)=AoAo(:,1:3)/180.0*pi
    AoAAmpl(:,1:3)=AoAAmpl(:,1:3)/180.0*pi
    AoAPhi(:,1:3)=AoAPhi(:,1:3)/180.0*pi
    XYZPhi(:,1:3)=XYZPhi(:,1:3)/180.0*pi
!   reference values: length, velocity, time
    Lref  = Lchod

    if(RefVelocity==0) then
        Uref = dabs(uuuIn(1))
    elseif(RefVelocity==1) then
        Uref = dabs(uuuIn(2))
    elseif(RefVelocity==2) then
        Uref = dabs(uuuIn(3))
    elseif(RefVelocity==3) then
        Uref = dsqrt(uuuIn(1)**2 + uuuIn(2)**2 + uuuIn(3)**2)  
    elseif(RefVelocity==10) then
        Uref = Lref * MAXVAL(Freq(1:nFish))
    elseif(RefVelocity==11) then
        do iFish=1,nFish
        nUref(iFish)=2.d0*pi*Freq(iFish)*MAXVAL(dabs(xyzAmpl(iFish,1:3)))
        enddo
        Uref = MAXVAL(nUref(1:nFish))
    elseif(RefVelocity==12) then
        do iFish=1,nFish
        nUref(iFish)=2.d0*pi*Freq(iFish)*MAXVAL(dabs(xyzAmpl(iFish,1:3)))*2.D0 !Park 2017 pof
        enddo
        Uref = MAXVAL(nUref(1:nFish))
    else
        Uref = 1.d0
    endif   

    Tref = Lref / Uref
    St = Lref * Freq / Uref
    g(1:3)=Frod(1:3) * Uref ** 2/Lref

    ! if(maxval(dabs(uuuIn(1:3)))>0.0001d0)then !with incoming flow
    !     Uref    = maxval(dabs(uuuIn(1:3)))        
    !     if    (MAXVAL(dabs(xyzAmpl(1:3)))>0.0001d0)then
    !         !Freq    = St*Uref/(2.0*MAXVAL(dabs(xyzAmpl(1:3))))
    !         Tref    = 1.0d0/Freq
    !         St      =  2.0*MAXVAL(dabs(xyzAmpl(1:3)))*Freq/maxval(dabs(uuuIn(1:3))) ! St=2fA/U
    !         if(iFlapRef==1) Uref    = Lref*Freq           
    !     elseif(dabs(MAXVAL(AoAAmpl(1:3)))>0.0001d0)then
    !         Tref    = 1.0d0/Freq            
    !     else
    !         Tref    = Lref/Uref            
    !     endif
    !     g(1:3)=Frod(1:3)*Uref**2/Lref
    ! else                                    !without incoming flow
    !     if    (MAXVAL(dabs(xyzAmpl(1:3)))>0.0001d0)then !heaving 
    !     !Uref    = 2.0*pi*MAXVAL(dabs(xyzAmpl(1:3)))*Freq  !maximum velocity
    !     !Uref    = 4.0*   MAXVAL(dabs(xyzAmpl(1:3)))*Freq  !mean velocity
    !     !Uref    =        MAXVAL(dabs(xyzAmpl(1:3)))*Freq  !mean velocity
    !     Uref    = Lref*Freq
    !     Tref    = 1.0d0/Freq
    !     St      = 1000.0
    !     g(1:3)=Frod(1:3)*Uref**2/Lref       
    !     elseif(MAXVAL(dabs(AoAAmpl(1:3)))>0.0001d0)then !picthing
    !     !Uref = 2.0*pi*MAXVAL(dabs(AoAAmpl(1:3))*[maxval(dabs(xyzful00(:,2))),maxval(dabs(xyzful00(:,1))),maxval(dabs(xyzful00(:,3)))])*Freq  !maximum velocity
    !     Uref  = 4.0*   MAXVAL(dabs(AoAAmpl(1:3))*[maxval(dabs(xyzful00(:,2))),maxval(dabs(xyzful00(:,1))),maxval(dabs(xyzful00(:,3)))])*Freq  !mean velocity
    !     !Uref    =     MAXVAL(dabs(AoAAmpl(1:3))*[maxval(dabs(xyzful00(:,2))),maxval(dabs(xyzful00(:,1))),maxval(dabs(xyzful00(:,3)))])*Freq  !mean velocity
    !     Tref    = 1.0d0/Freq
    !     St      = 1000.0
    !     g(1:3)=Frod(1:3)*Uref**2/Lref 
    !     elseif(maxval(dabs(Frod(1:3)))>0.0001d0)then    !gravity
    !         Uref=maxval(dabs(Frod(1:3)))
    !         g(1:3)=(Frod(1:3))**2/(Lref*denR)
    !         if(Frod(1)<0.000001d0) g(1)=-g(1)
    !         if(Frod(2)<0.000001d0) g(2)=-g(2)
    !         if(Frod(3)<0.000001d0) g(3)=-g(3)
    !         Tref    = Lref/Uref
    !     else
    !         stop 'no flow/moving body/gravity'
    !     endif
    ! endif

    uMax=maxval([maxval(dabs(uuuIn(1:3))),2.0*pi*MAXVAL(dabs(xyzAmpl(iFish,1:3)))*Freq, &
    2.0*pi*MAXVAL(dabs(AoAAmpl(iFish,1:3))*[maxval(dabs(xyzful00(iFish,:,2))),maxval(dabs(xyzful00(iFish,:,1))),maxval(dabs(xyzful00(iFish,:,3)))])*Freq])
    
    
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
        write(*,*)'uniform grid,STLBM'
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
    do iFish=1,nFish

    nt(iFish)=ele(iFish,1,4)  
    if(iKB==0)then
        prop(iFish,1:nMT(iFish),1) = EmR(iFish)*denIn*Uref**2
        prop(iFish,1:nMT(iFish),2) = prop(iFish,1:nMT(iFish),1)/2.0d0/(1.0+psR(iFish))
        Lthck= tcR(iFish)*Lref
        prop(iFish,1:nMT(iFish),3) = tcR(iFish)*Lref
        prop(iFish,1:nMT(iFish),4) = denR(iFish)*Lref*denIn/prop(iFish,1:nMT(iFish),3)        
        if    (nt(iFish)==2)then   !frame
        prop(iFish,1:nMT(iFish),7) = prop(iFish,1:nMT(iFish),3)**3/12.0d0
        prop(iFish,1:nMT(iFish),8) = prop(iFish,1:nMT(iFish),3)**3/12.0d0
        elseif(nt(iFish)==3)then   !plate
        prop(iFish,1:nMT(iFish),6) = prop(iFish,1:nMT(iFish),3)**3/12.0d0
        else
        endif
        KB=prop(iFish,nMT(iFish),1)*prop(iFish,nMT(iFish),6)/(denIn*Uref**2*Lref**3)
        KS=prop(iFish,nMT(iFish),1)*prop(iFish,nMT(iFish),3)/(denIn*Uref**2*Lref)
    endif
    
    if(iKB==1)then
        prop(iFish,1:nMT(iFish),3) = dsqrt(KB(iFish)/KS(iFish)*12.0d0)*Lref
        prop(iFish,1:nMT(iFish),4) = denR(iFish)*Lref*denIn/prop(iFish,1:nMT(iFish),3)
        if    (nt(iFish)==2)then   !frame
        prop(iFish,1:nMT(iFish),1) = KS(iFish)*denIn*Uref**2*Lref/prop(iFish,1:nMT(iFish),3)
        prop(iFish,1:nMT(iFish),2) = prop(iFish,1:nMT(iFish),1)/2.0d0/(1.0d0+psR(iFish))
        prop(iFish,1:nMT(iFish),7) = prop(iFish,1:nMT(iFish),3)**3/12.0d0
        prop(iFish,1:nMT(iFish),8) = prop(iFish,1:nMT(iFish),3)**3/12.0d0
        elseif(nt(iFish)==3)then   !plate
        prop(iFish,1:nMT(iFish),1) = KS(iFish)*denIn*Uref**2*Lref/prop(iFish,1:nMT(iFish),3)
        prop(iFish,1:nMT(iFish),2) = prop(iFish,1:nMT(iFish),1)/2.0d0/(1.0d0+psR(iFish))
        prop(iFish,1:nMT(iFish),6) = prop(iFish,1:nMT(iFish),3)**3/12.0d0
        else
        endif
        EmR(iFish) = prop(iFish,nMT(iFish),1)/(denIn*Uref**2) 
        tcR(iFish) = prop(iFish,nMT(iFish),3)/Lref
        Lthck=prop(iFish,nMT(iFish),3)
    endif
    enddo

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
    allocate(TTT00(1:nFish,1:3,1:3),TTT0(1:nFish,1:3,1:3),TTTnow(1:nFish,1:3,1:3),TTTnxt(1:nFish,1:3,1:3))
    allocate(XYZ(1:nFish,1:3),XYZd(1:nFish,1:3),UVW(1:nFish,1:3) )
    allocate(AoA(1:nFish,1:3),AoAd(1:nFish,1:3),WWW1(1:nFish,1:3),WWW2(1:nFish,1:3),WWW3(1:nFish,1:3) )
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