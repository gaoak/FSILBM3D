!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_file()
    USE simParam
    USE SolidSolver
    implicit none
    real(8):: iXYZ(1:3),dXYZ(1:3)
    integer:: i,iFish,iKind,FishKind,Order0
    integer:: FishOrder1,FishOrder2,LineX,LineY,LineZ
    integer, allocatable:: FishNum(:),NumX(:),NumY(:)
    character(LEN=40):: tmpFEmeshName
    integer:: niBodyModel,nisMotionGiven(1:6)
    real(8):: ndenR,npsR,nEmR,ntcR,nKB,nKS,nFreq,nSt
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
        allocate(FEmeshName(1:nFish),iBodyModel(1:nFish))
        allocate(FishNum(1:(FishKind+1)),NumX(1:FishKind),NumY(1:FishKind))
        FishNum(1)=1
        FishOrder1=0
        FishOrder2=0
    endif

    call readequal(111)
    do iKind=1,FishKind
        read(111,*)     FishNum(iKind+1),NumX(iKind),NumY(iKind)
        read(111,*)     niBodyModel, tmpFEmeshName
        read(111,*)     nisMotionGiven(1:3)
        read(111,*)     nisMotionGiven(4:6)
        read(111,*)     ndenR, npsR
        if(iKB==0) read(111,*)     nEmR, ntcR
        if(iKB==1) read(111,*)     nKB, nKS
        FishOrder1=FishOrder1+FishNum(iKind  )
        FishOrder2=FishOrder2+FishNum(iKind+1)
        do iFish=FishOrder1,FishOrder2
            iBodyModel(iFish)=niBodyModel
            FEmeshName(iFish)=tmpFEmeshName
            BeamInfo(iFish)%isMotionGiven(1:6)=nisMotionGiven(1:6)
            BeamInfo(iFish)%denR= ndenR
            BeamInfo(iFish)%psR = npsR
            if(iKB==0) then
            BeamInfo(iFish)%EmR = nEmR
            BeamInfo(iFish)%tcR = ntcR
            elseif(iKB==1) then
            BeamInfo(iFish)%KB  = nKB
            BeamInfo(iFish)%KS  = nKS
            endif
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
            BeamInfo(iFish)%Freq=nFreq
            BeamInfo(iFish)%St  =nSt
            BeamInfo(iFish)%XYZAmpl(1:3)=nXYZAmpl(1:3)
            BeamInfo(iFish)%XYZPhi(1:3) =nXYZPhi(1:3)
            BeamInfo(iFish)%AoAo(1:3)   =nAoAo(1:3)
            BeamInfo(iFish)%AoAAmpl(1:3)=nAoAAmpl(1:3)
            BeamInfo(iFish)%AoAPhi(1:3) =nAoAPhi(1:3)
            ! initial position distribution
            Order0 = iFish - FishOrder1
            LineX  = mod(Order0,NumX(iKind))
            LineY  = mod(Order0/NumX(iKind),NumY(iKind))
            LineZ  = Order0/(NumX(iKind)*NumY(iKind))
            BeamInfo(iFish)%XYZo(1) = iXYZ(1) + dXYZ(1) * LineX
            BeamInfo(iFish)%XYZo(2) = iXYZ(2) + dXYZ(2) * LineY
            BeamInfo(iFish)%XYZo(3) = iXYZ(3) + dXYZ(3) * LineZ
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
        read(111,*)     SampBodyNode(1:numSampBody,1)
        do iFish=1,nFish
        SampBodyNode(1:numSampBody,iFish)=SampBodyNode(1:numSampBody,1)
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
!    calculate Bolztman parameters from flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE calculate_LB_params()
    USE simParam
    USE SolidSolver
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
        Uref = Lref * MAXVAL(BeamInfo(:)%Freq)
    elseif(RefVelocity==11) then
        do iFish=1,nFish
        nUref(iFish)=2.d0*pi*BeamInfo(iFish)%Freq*MAXVAL(dabs(BeamInfo(iFish)%xyzAmpl(1:3)))
        enddo
        Uref = MAXVAL(nUref(1:nFish))
    elseif(RefVelocity==12) then
        do iFish=1,nFish
        nUref(iFish)=2.d0*pi*BeamInfo(iFish)%Freq*MAXVAL(dabs(BeamInfo(iFish)%xyzAmpl(1:3)))*2.D0 !Park 2017 pof
        enddo
        Uref = MAXVAL(nUref(1:nFish))
    !else
        !Uref = 1.0d0
    endif

    if(RefTime==0) then
        Tref = Lref / Uref
    elseif(RefTime==1) then
        Tref = 1 / maxval(BeamInfo(:)%Freq)
    !else
    endif

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
    Eref=0.5*denIn*Uref**2*Asfac*Lref
    Pref=0.5*denIn*Uref**2*Asfac*Uref

    g(1:3)=Frod(1:3) * Uref ** 2/Lref
    uMax = 0.
    do iFish = 1,nFish
        call BeamInfo(iFish)%calculate_angle_material
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
    USE SolidSolver
    IMPLICIT NONE
    integer:: iFish
    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='replace')
    write(13) step,time
    write(13) fIn,xGrid,yGrid,zGrid
    write(13) nFish
    write(13) IXref,IYref,IZref,NDref
    do iFish = 1,nFish
        call BeamInfo(iFish)%write_solid_temp(13)
    enddo
    write(13) UPre,UNow
    close(13)
    ENDSUBROUTINE


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read check point file for restarting simulation
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_checkpoint_file()
    USE simParam
    USE SolidSolver
    IMPLICIT NONE
    integer:: tmpnfish,iFish

    open(unit=13,file='./DatTemp/conwr.dat',form='unformatted',status='old')
    read(13) step,time
    read(13) fIn,xGrid,yGrid,zGrid
    read(13) tmpnfish
    if((tmpnfish .eq. nFish)) then
        read(13) IXref,IYref,IZref,NDref
        do iFish = 1,nFish
        call BeamInfo(iFish)%read_solid_temp(13)
    enddo
        read(13) UPre,UNow
    endif
    close(13)
    ENDSUBROUTINE
