!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read flow parameters
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_file()
    USE simParam
    USE SolidBody
    implicit none
    real(8):: iXYZ(1:3),dXYZ(1:3)
    integer:: i,iFish,iKind,FishKind,Order0
    integer:: FishOrder1,FishOrder2,LineX,LineY,LineZ
    integer, allocatable:: FishNum(:),NumX(:),NumY(:)
    character(LEN=40):: tmpFEmeshName
    integer:: niBodyModel,niBodyType,nisMotionGiven(1:6)
    real(8):: ndenR,npsR,nEmR,ntcR,nKB,nKS,nFreq,nSt
    real(8):: nXYZAmpl(1:3),nXYZPhi(1:3),nAoAo(1:3),nAoAAmpl(1:3),nAoAPhi(1:3)
    character (LEN=40), allocatable:: FEmeshName(:)
    integer, allocatable:: isMotionGiven(:,:)
    real(8), allocatable:: denR(:),KB(:),KS(:),EmR(:),psR(:),tcR(:),St(:)
    real(8), allocatable:: Freq(:)
    real(8), allocatable:: XYZo(:,:),XYZAmpl(:,:),XYZPhi(:,:)
    real(8), allocatable:: AoAo(:,:),AoAAmpl(:,:),AoAPhi(:,:)
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
        allocate(FEmeshName(1:nFish),iBodyModel(1:nFish),iBodyType(1:nFish),isMotionGiven(1:DOFDim,1:nFish))
        allocate(denR(1:nFish),EmR(1:nFish),tcR(1:nFish),psR(1:nFish),KB(1:nFish),KS(1:nFish))
        allocate(FishNum(1:(FishKind+1)),NumX(1:FishKind),NumY(1:FishKind))
        FishNum(1)=1
        FishOrder1=0
        FishOrder2=0
    endif

    call readequal(111)
    do iKind=1,FishKind
        read(111,*)     FishNum(iKind+1),NumX(iKind),NumY(iKind)
        read(111,*)     niBodyModel,niBodyType,tmpFEmeshName
        read(111,*)     nisMotionGiven(1:3)
        read(111,*)     nisMotionGiven(4:6)
        read(111,*)     ndenR, npsR
        if(iKB==0) read(111,*)     nEmR, ntcR
        if(iKB==1) read(111,*)     nKB, nKS
        FishOrder1=FishOrder1+FishNum(iKind  )
        FishOrder2=FishOrder2+FishNum(iKind+1)
        do iFish=FishOrder1,FishOrder2
            iBodyModel(iFish)=niBodyModel
            iBodyType(iFish)=niBodyType
            FEmeshName(iFish)=tmpFEmeshName
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
        read(111,*)     SampBodyNode(1:numSampBody,1)
        do iFish=1,nFish
        SampBodyNode(1:numSampBody,iFish)=SampBodyNode(1:numSampBody,1)
        enddo
    endif
    call readequal(111)
    close(111)

    call read_solid_file(nFish,FEmeshName,iBodyModel,iBodyType,isMotionGiven,denR,KB,KS,EmR,psR,tcR,St, &
                         Freq,XYZo,XYZAmpl,XYZPhi,AoAo,AoAAmpl,AoAPhi, &
                         ntolLBM,dtolLBM,Pbeta,dt,denIn,uuuIn,boundaryConditions, &
                         dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,dtolFEM,ntolFEM,iForce2Body,iKB)
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
