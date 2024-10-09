!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_file()
    USE simParam
    USE ImmersedBoundary
    implicit none
    real(8):: iXYZ(1:3),dXYZ(1:3)
    integer:: i,iFish,iKind,FishKind,Order0
    integer:: FishOrder1,FishOrder2,LineX,LineY,LineZ
    integer, allocatable:: FishNum(:),NumX(:),NumY(:)
    character(LEN=40):: nFEmeshName
    integer:: niBodyModel,nisMotionGiven(1:6),nNspan
    real(8):: ndenR,npsR,nEmR,ntcR,nKB,nKS,nFreq,nSt,ndspan,ntheta
    real(8):: nXYZAmpl(1:3),nXYZPhi(1:3),nAoAo(1:3),nAoAAmpl(1:3),nAoAPhi(1:3)
    open(unit=111,file='inFlow.dat')
    call readequal(111)
    read(111,*)     npsize
    read(111,*)     isRelease,isConCmpt
    read(111,*)     iCollidModel,iKB
    read(111,*)     RefVelocity,Uref
    read(111,*)     timeSimTotl,timeOutTemp
    read(111,*)     timeOutFlow,timeOutBody,timeOutInfo
    read(111,*)     timeOutBegin,timeOutEnd
    read(111,*)     Palpha,Pbeta,Pramp
    if(timeOutBegin.gt.timeOutEnd) then
        write(*,*) 'The start time for output should be earlier than the end time.'
        stop
    endif
    call readequal(111)
    read(111,*)     uuuIn(1:3)
    read(111,*)     shearRateIn(1:3)
    read(111,*)     boundaryConditions(1:6)
    read(111,*)     VelocityKind
    if(VelocityKind==2) then
        VelocityAmp = shearRateIn(1)
        VelocityFreq = shearRateIn(2)
        VelocityPhi = shearRateIn(3)
    endif
    read(111,*)     MovingKind1,MovingVel1,MovingFreq1
    read(111,*)     MovingKind2,MovingVel2,MovingFreq2
    call readequal(111)
    read(111,*)     VolumeForceIn(1:SpcDim)
    read(111,*)     VolumeForceAmp,VolumeForceFreq,VolumeForcePhi
    call readequal(111)
    read(111,*)     Re,dt
    read(111,*)     RefTime,Tref
    read(111,*)     Frod(1:3)            !Gravity
    call readequal(111)
    read(111,*)     iBC
    read(111,*)     LBmeshName
    read(111,*)     denIn
    read(111,*)     dtolLBM,ntolLBM      !velocity iteration
    call readequal(111)
    read(111,*)     nFish,FishKind,iForce2Body
    if(nFish.lt.FishKind) then
        write(*, *) "Fish kind is more than fish number ", FishKind, nFish
    endif
    if(nFish.eq.0) iForce2Body = 0

    if(nFish>0) then
        allocate(FEmeshName(1:nFish),iBodyModel(1:nFish),isMotionGiven(1:DOFDim,1:nFish))
        allocate(denR(1:nFish),EmR(1:nFish),tcR(1:nFish),psR(1:nFish),KB(1:nFish),KS(1:nFish))
        allocate(dspan(1:nFish),theta(1:nFish),Nspan(1:nFish))
        allocate(FishNum(1:(FishKind+1)),NumX(1:FishKind),NumY(1:FishKind))
        FishNum(1)=1
        FishOrder1=0
        FishOrder2=0
    endif

    call readequal(111)
    do iKind=1,FishKind
        read(111,*)     FishNum(iKind+1),NumX(iKind),NumY(iKind)
        read(111,*)     niBodyModel, nFEmeshName
        read(111,*)     nisMotionGiven(1:3)
        read(111,*)     nisMotionGiven(4:6)
        read(111,*)     ndenR, npsR
        if(iKB==0) read(111,*)     nEmR, ntcR
        if(iKB==1) read(111,*)     nKB, nKS
        read(111,*)     ndspan,ntheta,nNspan
        FishOrder1=FishOrder1+FishNum(iKind  )
        FishOrder2=FishOrder2+FishNum(iKind+1)
        do iFish=FishOrder1,FishOrder2
            iBodyModel(iFish)=niBodyModel
            FEmeshName(iFish)=nFEmeshName
            isMotionGiven(1:6,iFish)=nisMotionGiven(1:6)
            denR(iFish)= ndenR
            psR(iFish) = npsR
            if(iKB==0) then
            EmR(iFish) = nEmR
            tcR(iFish) = ntcR
            elseif(iKB==1) then
            KB(iFish)  = nKB
            KS(iFish)  = nKS
            endif
            dspan(iFish) = ndspan
            theta(iFish) = ntheta
            Nspan(iFish) = nNspan
        enddo
    enddo

    call readequal(111)
    read(111,*)     numsubstep
    read(111,*)     dampK,dampM
    read(111,*)     NewmarkGamma,NewmarkBeta
    read(111,*)     alphaf,alpham,alphap
    read(111,*)     dtolFEM,ntolFEM
    call readequal(111)

    if(nFish>0) then
        allocate(Freq(1:nFish),St(1:nFish))
        allocate(XYZo(1:3,1:nFish),XYZAmpl(1:3,1:nFish),XYZPhi(1:3,1:nFish))
        allocate(AoAo(1:3,1:nFish),AoAAmpl(1:3,1:nFish),AoAPhi(1:3,1:nFish))
        FishOrder1 =0
        FishOrder2 =0
    endif

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
            XYZAmpl(1:3,iFish)=nXYZAmpl(1:3)
            XYZPhi(1:3,iFish) =nXYZPhi(1:3)
            AoAo(1:3,iFish)   =nAoAo(1:3)
            AoAAmpl(1:3,iFish)=nAoAAmpl(1:3)
            AoAPhi(1:3,iFish) =nAoAPhi(1:3)
            ! initial position distribution
            Order0 = iFish - FishOrder1
            LineX  = mod(Order0,NumX(iKind))
            LineY  = mod(Order0/NumX(iKind),NumY(iKind))
            LineZ  = Order0/(NumX(iKind)*NumY(iKind))
            XYZo(1,iFish) = iXYZ(1) + dXYZ(1) * LineX
            XYZo(2,iFish) = iXYZ(2) + dXYZ(2) * LineY
            XYZo(3,iFish) = iXYZ(3) + dXYZ(3) * LineZ
        enddo
    enddo
    call readequal(111)
    read(111,*)     isMoveGrid,offsetOutput
    read(111,*)     isMoveDimX,isMoveOutputX
    read(111,*)     isMoveDimY,isMoveOutputY
    read(111,*)     isMoveDimZ,isMoveOutputZ
    read(111,*)     Xref,Yref,Zref
    if(nFish .eq. 0 .and. isMoveGrid .eq. 1) then
        write(*, *) 'No fish, resetting move grid as false'
        isMoveGrid = 0
    endif
    call readequal(111)
    read(111,*)     waveInitDist,AmplInitDist(1:SpcDim)
    read(111,*)     FreqForcDist,AmplForcDist(1:SpcDim)
    read(111,*)     begForcDist,endForcDist
    read(111,*)     posiForcDist(1:SpcDim)
    call readequal(111)

    read(111,*)     numSampFlow,isFluidOutput
    allocate(SampFlowPint(1:numSampFlow,1:3))
    do i=1,numSampFlow
        read(111,*) SampFlowPint(i,1:3)
    enddo
    call readequal(111)
    if(nFish>0) then
        read(111,*)     numSampBody,isBodyOutput
        allocate(SampBodyNode(1:numSampBody,1:nFish))
        read(111,*)     SampBodyNode(1,1:numSampBody)
        do iFish=1,nFish
        SampBodyNode(1:numSampBody,iFish)=SampBodyNode(1,1:numSampBody)
        enddo
    endif
    call readequal(111)
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
    allocate(den(zDim,yDim,xDim),prs(zDim,yDim,xDim))

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
    ! allocate workspace for flow field output
    call initOutFlowWorkspace
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
    USE ImmersedBoundary
    implicit none
    integer:: iEL,iND,iFish
    real(8):: lentemp
    if(nFish==0) return
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
    allocate( ele_all(nEL_all,5),xyzful_all(nND_all,6),velful_all(nND_all,6),extful_all(nND_all,6))

    allocate( ele(nEL_max,5,1:nFish),xyzful00(nND_max,6,1:nFish),xyzful0(nND_max,6,1:nFish),mssful(nND_max,6,1:nFish),lodful(nND_max,6,1:nFish), &
              extful(nND_max,6,1:nFish),repful(nND_max,1:6,1:nFish),extful1(nND_max,6,1:nFish),extful2(nND_max,6,1:nFish),nloc(nND_max*6,1:nFish),nprof(nND_max*6,1:nFish), &
              nprof2(nND_max*6,1:nFish),jBC(nND_max,6,1:nFish))
    allocate( grav(nND_max,6,1:nFish),vBC(nND_max,6,1:nFish),mss(nND_max*6,1:nFish),prop(nMT_max,10,1:nFish),areaElem00(nEL_max,1:nFish),areaElem(nEL_max,1:nFish))

    allocate( xyzful(nND_max,6,1:nFish),xyzfulnxt(nND_max,6,1:nFish),dspful(nND_max,6,1:nFish),velful(nND_max,6,1:nFish),accful(nND_max,6,1:nFish))
    allocate( triad_nn(3,3,nND_max,1:nFish),triad_ee(3,3,nEL_max,1:nFish),triad_e0(3,3,nEL_max,1:nFish) )
    allocate( triad_n1(3,3,nEL_max,1:nFish),triad_n2(3,3,nEL_max,1:nFish),triad_n3(3,3,nEL_max,1:nFish) )

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

    call read_structural_datafile(jBC(1:nND(iFish),1:6,iFish),ele(1:nEL(iFish),1:5,iFish),nloc(1:nND(iFish)*6,iFish),nprof(1:nND(iFish)*6,iFish), &
                                      nprof2(1:nND(iFish)*6,iFish),xyzful00(1:nND(iFish),1:6,iFish),prop(1:nMT(iFish),1:10,iFish),nND(iFish), &
                                      nEL(iFish),nEQ(iFish),nMT(iFish),nBD(iFish),nSTF(iFish),idat)
    close(idat)
    if (maxval(Nspan).gt.0 .and. maxval(dabs(prop(1:nMT(iFish),5,iFish))).gt.1d-6) then
        write(*,*) 'Extruded body should have zero rotation angle, gamma', prop(1:nMT(iFish),5,iFish)
        stop
    endif
    write(*,*)'read FEMeshFile ',iFish,' end'
    enddo
!    ===============================================================================================
!    package
    do iFish=1,nFish
        if(iFish.eq.1)then
            do iEL=1,nEL(iFish)
                ele_all(iEL,1:3)=ele(iEL,1:3,iFish)
                ele_all(iEL,4:5)=ele(iEL,4:5,iFish)
            enddo
        else
            do iEL=1,nEL(iFish)
                ele_all(iEL+sum(nEL(1:iFish-1)),1:3)=ele(iEL,1:3,iFish)+sum(nND(1:iFish-1))
                ele_all(iEL+sum(nEL(1:iFish-1)),4:5)=ele(iEL,4:5,iFish)
            enddo
        endif
    enddo
!   ===============================================================================================
    iFish = 1
    Lchod   = maxval(xyzful00(:,1,iFish))-minval(xyzful00(:,1,iFish))
    lentemp = maxval(xyzful00(:,2,iFish))-minval(xyzful00(:,2,iFish))
    if(lentemp .gt. Lchod) Lchod = lentemp
    if (maxval(Nspan).gt.0) then
        Lspan = maxval(dspan*Nspan)
    else
        Lspan = maxval(xyzful00(:,3,iFish))-minval(xyzful00(:,3,iFish))
    endif
    Asfac = Lchod * Lspan
    AR    = Lspan / Lchod

!   loading boundary type*******************************************************************************
    do  iFish=1,nFish
    do    iND=1,nND(iFish)
        if(jBC(iND,1,iFish)==1) jBC(iND,1:6,iFish)=isMotionGiven(1:6,iFish)
    enddo
    enddo
    END SUBROUTINE allocate_solid_memory

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
        if(VelocityKind==0) then
            call evaluateShearVelocity(xGrid(x),yGrid(y),zGrid(z), vel)
        elseif(VelocityKind==2) then
            call evaluateOscillatoryVelocity(vel)
        endif
        uuu(z, y, x, 1) = vel(1)
        uuu(z, y, x, 2) = vel(2)
        uuu(z, y, x, 3) = vel(3)
    enddo
    enddo
    enddo
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
    END SUBROUTINE initialize_flow

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    initialize solid field
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initialize_solid()
    USE simParam
    USE ImmersedBoundary
    implicit none
    integer:: iND,iFish,iCount
    if(nFish.eq.0) return
    allocate(TTT00(1:3,1:3,1:nFish),TTT0(1:3,1:3,1:nFish),TTTnow(1:3,1:3,1:nFish),TTTnxt(1:3,1:3,1:nFish))
    allocate(XYZ(1:3,1:nFish),XYZd(1:3,1:nFish),UVW(1:3,1:nFish) )
    allocate(AoA(1:3,1:nFish),AoAd(1:3,1:nFish),WWW1(1:3,1:nFish),WWW2(1:3,1:nFish),WWW3(1:3,1:nFish) )
    iCount = 0
    do iFish=1,nFish
        TTT00(:,:,iFish)=0.0d0
        TTT00(1,1,iFish)=1.0d0
        TTT00(2,2,iFish)=1.0d0
        TTT00(3,3,iFish)=1.0d0

        XYZ(1:3,iFish)=XYZo(1:3,iFish)+XYZAmpl(1:3,iFish)*dcos(2.0*pi*Freq(iFish)*time+XYZPhi(1:3,iFish))
        AoA(1:3,iFish)=AoAo(1:3,iFish)+AoAAmpl(1:3,iFish)*dcos(2.0*pi*Freq(iFish)*time+AoAPhi(1:3,iFish))

        call AoAtoTTT(AoA(1:3,iFish),TTT0(1:3,1:3,iFish))
        call AoAtoTTT(AoA(1:3,iFish),TTTnow(1:3,1:3,iFish))
        call get_angle_triad(TTT0(1:3,1:3,iFish),TTTnow(1:3,1:3,iFish),AoAd(1,iFish),AoAd(2,iFish),AoAd(3,iFish))

        do iND=1,nND(iFish)
            xyzful0(iND,1:3,iFish)=matmul(TTT0(1:3,1:3,iFish),xyzful00(iND,1:3,iFish))+XYZ(1:3,iFish)
            xyzful0(iND,4:6,iFish)=AoAd(1:3,iFish)
        enddo

        xyzful(1:nND(iFish),1:6,iFish)=xyzful0(1:nND(iFish),1:6,iFish)
        velful(1:nND(iFish),1:6,iFish)=0.0

        dspful(1:nND(iFish),1:6,iFish)=0.0
        accful(1:nND(iFish),1:6,iFish)=0.0

        do iND=1,nND(iFish)
            xyzful_all(iND+iCount,1:6)   =xyzful(iND,1:6,iFish)
            velful_all(iND+iCount,1:6)   =velful(iND,1:6,iFish)
            extful_all(iND+iCount,1:6)  =0.d0
        enddo

        CALL formmass_D(ele(1:nEL(iFish),1:5,iFish),xyzful0(1:nND(iFish),1,iFish),xyzful0(1:nND(iFish),2,iFish),xyzful0(1:nND(iFish),3,iFish), &
                        prop(1:nMT(iFish),1:10,iFish),mss(1:nND(iFish)*6,iFish),nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),alphaf)

        do iND = 1, nND(iFish)
            mssful(iND,1:6,iFish)= mss((iND-1)*6+1:(iND-1)*6+6,iFish)
            grav(iND,1:6,iFish)  = mssful(iND,1:6,iFish)*[g(1),g(2),g(3),0.0d0,0.0d0,0.0d0]
        enddo

        CALL init_triad_D(ele(1:nEL(iFish),1:5,iFish),xyzful(1:nND(iFish),1,iFish),xyzful(1:nND(iFish),2,iFish),xyzful(1:nND(iFish),3,iFish), &
                          triad_nn(1:3,1:3,1:nND(iFish),iFish),triad_n1(1:3,1:3,1:nEL(iFish),iFish),triad_n2(1:3,1:3,1:nEL(iFish),iFish), &
                          triad_ee(1:3,1:3,1:nEL(iFish),iFish),triad_e0(1:3,1:3,1:nEL(iFish),iFish),nND(iFish),nEL(iFish))

        iCount = iCount + nND(iFish)
    enddo
    call initializexyzIB
    END SUBROUTINE initialize_solid

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    calculate Bolztman parameters from flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE calculate_LB_params()
    USE simParam
    implicit none
    integer:: nt(1:nFish),iFish
    real(8):: nUref(1:nFish)
!   reference values: length, velocity, time
    if(nFish.eq.0) then
        Lref = 1.d0
    else
        Lref  = Lchod
    endif

    if(RefVelocity==0) then
        Uref = dabs(uuuIn(1))
    elseif(RefVelocity==1) then
        Uref = dabs(uuuIn(2))
    elseif(RefVelocity==2) then
        Uref = dabs(uuuIn(3))
    elseif(RefVelocity==3) then
        Uref = dsqrt(uuuIn(1)**2 + uuuIn(2)**2 + uuuIn(3)**2)
    elseif(RefVelocity==4) then
        Uref = dabs(VelocityAmp)  !Velocity Amplitude
    elseif(RefVelocity==10) then
        Uref = Lref * MAXVAL(Freq(1:nFish))
    elseif(RefVelocity==11) then
        do iFish=1,nFish
        nUref(iFish)=2.d0*pi*Freq(iFish)*MAXVAL(dabs(xyzAmpl(1:3,iFish)))
        enddo
        Uref = MAXVAL(nUref(1:nFish))
    elseif(RefVelocity==12) then
        do iFish=1,nFish
        nUref(iFish)=2.d0*pi*Freq(iFish)*MAXVAL(dabs(xyzAmpl(1:3,iFish)))*2.D0 !Park 2017 pof
        enddo
        Uref = MAXVAL(nUref(1:nFish))
    !else
        !Uref = 1.0d0
    endif

    if(RefTime==0) then
        Tref = Lref / Uref
    elseif(RefTime==1) then
        Tref = 1 / maxval(Freq(:))
    !else
    endif

    do iFish=1,nFish
        St(iFish) = Lref * Freq(iFish) / Uref
    enddo
    g(1:3)=Frod(1:3) * Uref ** 2/Lref
    uMax = 0.
    do iFish=1,nFish
        ! angle to radian
        AoAo(1:3,iFish)=AoAo(1:3,iFish)/180.0*pi
        AoAAmpl(1:3,iFish)=AoAAmpl(1:3,iFish)/180.0*pi
        AoAPhi(1:3,iFish)=AoAPhi(1:3,iFish)/180.0*pi
        XYZPhi(1:3,iFish)=XYZPhi(1:3,iFish)/180.0*pi
        uMax=maxval([uMax, maxval(dabs(uuuIn(1:3))),2.0*pi*MAXVAL(dabs(xyzAmpl(1:3,iFish)))*Freq(iFish), &
            2.0*pi*MAXVAL(dabs(AoAAmpl(1:3,iFish))*[maxval(dabs(xyzful00(:,2,iFish))), &
            maxval(dabs(xyzful00(:,1,iFish))),maxval(dabs(xyzful00(:,3,iFish)))])*Freq(iFish)])
    enddo

!   calculate viscosity, LBM relexation time
    ratio  =  dt/dh
    if(ratio>1.0d0+eps)then
        write(*,*)'dt >  dhmin !!!!!!!!!!'
        write(*,*)'dt <= dhmin (we use streching mesh for LBM)'
        stop
    endif

    !for uniform grid, advection length equals grid size
    isUniformGrid(1) = dabs(dxmax/dh-1.0d0)<eps
    isUniformGrid(2) = dabs(dymax/dh-1.0d0)<eps
    isUniformGrid(3) = dabs(dzmax/dh-1.0d0)<eps
    if(dabs(dt/dh-1.0d0)<eps .and. isUniformGrid(1) .and. isUniformGrid(2) .and. isUniformGrid(3))then
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
    Eref=denIn*Uref**2*Asfac*Lref
    Pref=denIn*Uref**2*Asfac*Uref
    !calculate material parameters
    do iFish=1,nFish
    nt(iFish)=ele(1,4,iFish)
    if(iKB==0)then
        prop(1:nMT(iFish),1,iFish) = EmR(iFish)*denIn*Uref**2
        prop(1:nMT(iFish),2,iFish) = prop(1:nMT(iFish),1,iFish)/2.0d0/(1.0+psR(iFish))
        Lthck= tcR(iFish)*Lref
        prop(1:nMT(iFish),3,iFish) = tcR(iFish)*Lref
        prop(1:nMT(iFish),4,iFish) = denR(iFish)*Lref*denIn/prop(1:nMT(iFish),3,iFish)
        if    (nt(iFish)==2)then   !frame
        prop(1:nMT(iFish),7,iFish) = prop(1:nMT(iFish),3,iFish)**3/12.0d0
        prop(1:nMT(iFish),8,iFish) = prop(1:nMT(iFish),3,iFish)**3/12.0d0
        elseif(nt(iFish)==3)then   !plate
        prop(1:nMT(iFish),6,iFish) = prop(1:nMT(iFish),3,iFish)**3/12.0d0
        else
        endif
        KB=prop(nMT(iFish),1,iFish)*prop(nMT(iFish),6,iFish)/(denIn*Uref**2*Lref**3)
        KS=prop(nMT(iFish),1,iFish)*prop(nMT(iFish),3,iFish)/(denIn*Uref**2*Lref)
    endif

    if(iKB==1)then
        prop(1:nMT(iFish),3,iFish) = dsqrt(KB(iFish)/KS(iFish)*12.0d0)*Lref
        prop(1:nMT(iFish),4,iFish) = denR(iFish)*Lref*denIn/prop(1:nMT(iFish),3,iFish)
        if    (nt(iFish)==2)then   !frame
        prop(1:nMT(iFish),1,iFish) = KS(iFish)*denIn*Uref**2*Lref/prop(1:nMT(iFish),3,iFish)
        prop(1:nMT(iFish),2,iFish) = prop(1:nMT(iFish),1,iFish)/2.0d0/(1.0d0+psR(iFish))
        prop(1:nMT(iFish),7,iFish) = prop(1:nMT(iFish),3,iFish)**3/12.0d0
        prop(1:nMT(iFish),8,iFish) = prop(1:nMT(iFish),3,iFish)**3/12.0d0
        elseif(nt(iFish)==3)then   !plate
        prop(1:nMT(iFish),1,iFish) = KS(iFish)*denIn*Uref**2*Lref/prop(1:nMT(iFish),3,iFish)
        prop(1:nMT(iFish),2,iFish) = prop(1:nMT(iFish),1,iFish)/2.0d0/(1.0d0+psR(iFish))
        prop(1:nMT(iFish),6,iFish) = prop(1:nMT(iFish),3,iFish)**3/12.0d0
        else
        endif
        EmR(iFish) = prop(nMT(iFish),1,iFish)/(denIn*Uref**2)
        tcR(iFish) = prop(nMT(iFish),3,iFish)/Lref
        Lthck=prop(nMT(iFish),3,iFish)
    endif
    enddo

    END SUBROUTINE calculate_LB_params

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

    END SUBROUTINE calculate_MRTM_params



!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write check point file for restarting simulation
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_checkpoint_file()
    USE simParam
    USE ImmersedBoundary
    IMPLICIT NONE
    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='replace')
    write(13) step,time
    write(13) fIn,xGrid,yGrid,zGrid
    write(13) nFish, nND_max
    write(13) IXref,IYref,IZref,NDref
    write(13) xyzful0,xyzful,dspful,velful,accful,extful,mss,mssful,grav
    if(Palpha.gt.0.d0) write(13) xyzfulIB
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
    USE ImmersedBoundary
    IMPLICIT NONE
    integer:: tmpnfish, tmpND_max

    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='old')
    read(13) step,time
    read(13) fIn,xGrid,yGrid,zGrid
    read(13) tmpnfish, tmpND_max
    if((tmpnfish .eq. nFish) .and. (tmpND_max .eq. nND_max)) then
        read(13) IXref,IYref,IZref,NDref
        read(13) xyzful0,xyzful,dspful,velful,accful,extful,mss,mssful,grav
        if(Palpha.gt.0.d0) read(13) xyzfulIB
        read(13) triad_nn,triad_ee,triad_e0
        read(13) triad_n1,triad_n2,triad_n3
        read(13) UPre,UNow,Et,Ek,Ep,Es,Eb,Ew
    endif
    close(13)
    ENDSUBROUTINE
