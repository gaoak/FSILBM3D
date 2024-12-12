module FlowCondition
    implicit none
    private
    public :: flow
    type :: FlowCondType
        integer :: isRelease,isConCmpt,numsubstep
        real(8) :: timeSimTotal,timeWriteBegin,timeWriteEnd,timeWriteFlow,timeWriteBody,timeWriteInfo
        real(8) :: Re,dt,denIn,nu,Mu,dtolLBM
        integer :: TrefType,UrefType,ntolLBM
        real(8) :: uvwIn(1:3),shearRateIn(1:3)
        real(8) :: volumeForceIn(1:3),volumeForceAmp,volumeForceFreq,volumeForcePhi
        real(8) :: Uref,Lref,Tref
        real(8) :: Aref,Fref,Eref,Pref
        !real(8) :: velocityType,velocityAmp,velocityFreq,velocityPhi
        real(8):: Asfac,Lchod,Lspan,AR ! fish area, length
    end type FlowCondType
    type(FlowCondType) :: flow
    contains

    SUBROUTINE read_flow_conditions(filename)
        ! read computing core and block numbers
        implicit none
        character(LEN=40),intent(in):: filename
        open(unit=111, file=filename, status='old', action='read')
        call found_keyword(111,'FlowCondition')
        call readNextData(111, buffer)
        read(111,*)    flow%isRelease,flow%isConCmpt,flow%numsubstep
        call readNextData(111, buffer)
        read(111,*)    flow%timeSimTotal,flow%timeWriteBegin,flow%timeWriteEnd
        call readNextData(111, buffer)
        read(111,*)    flow%timeWriteFlow,flow%timeWriteBody,flow%timeWriteInfo
        call readNextData(111, buffer)
        read(111,*)    flow%Re,flow%dt,flow%denIn
        call readNextData(111, buffer)
        read(111,*)    flow%uvwIn(1:3)
        call readNextData(111, buffer)
        read(111,*)    flow%shearRateIn(1:3)
        call readNextData(111, buffer)
        read(111,*)    flow%volumeForceIn(1:3)
        call readNextData(111, buffer)
        read(111,*)    flow%volumeForceAmp,flow%volumeForceFreq,flow%volumeForcePhi
        call readNextData(111, buffer)
        read(111,*)    flow%TrefType,flow%Tref
        call readNextData(111, buffer)
        read(111,*)    flow%UrefType,flow%Uref
        call readNextData(111, buffer)
        read(111,*)    flow%ntolLBM,flow%dtolLBM
        close(111)
    END SUBROUTINE

    SUBROUTINE calculate_reference_params()
        USE ConstParams
        USE SolidBody
        implicit none
        integer:: iFish
        real(8):: nUref(1:m_nFish)
        ! reference length
        if(m_nFish.eq.0) then
            flow%Lref = 1.d0
        else
            flow%Lref  = flow%Lchod
        endif
        ! reference velocity
        if(flow%UrefType==0) then
            flow%Uref = dabs(flow%uvwIn(1))
        elseif(flow%UrefType==1) then
            flow%Uref = dabs(flow%uvwIn(2))
        elseif(flow%UrefType==2) then
            flow%Uref = dabs(flow%uvwIn(3))
        elseif(flow%UrefType==3) then
            flow%Uref = dsqrt(flow%uvwIn(1)**2 + flow%uvwIn(2)**2 + flow%uvwIn(3)**2)
        !elseif(flow%UrefType==4) then
        !    Uref = dabs(VelocityAmp)  !Velocity Amplitude
        elseif(flow%UrefType==5) then
            flow%Uref = flow%Lref * MAXVAL(VBodies(:)%rbm%Freq)
        elseif(flow%UrefType==6) then
            do iFish=1,m_nFish
            nUref(iFish)=2.d0*pi*VBodies(iFish)%rbm%Freq*MAXVAL(dabs(VBodies(iFish)%rbm%xyzAmpl(1:3)))
            enddo
            flow%Uref = MAXVAL(nUref(:))
        elseif(flow%UrefType==7) then
            do iFish=1,m_nFish
            nUref(iFish)=2.d0*pi*VBodies(iFish)%rbm%Freq*MAXVAL(dabs(VBodies(iFish)%rbm%xyzAmpl(1:3)))*2.D0 !Park 2017 pof
            enddo
            flow%Uref = MAXVAL(nUref(1:m_nFish))
        else
            write(*,*) 'use input reference velocity'
        endif
        ! reference time
        if(flow%TrefType==0) then
            flow%Tref = flow%Lref / flow%Uref
        elseif(flow%TrefType==1) then
            flow%Tref = 1 / maxval(VBodies(:)%rbm%Freq)
        else
            write(*,*) 'use input reference time'
        endif
        ! reference acceleration, force, energy, power
        flow%Aref = flow%Uref/flow%Tref
        flow%Fref = 0.5*flow%denIn*flow%Uref**2*flow%Asfac
        flow%Eref = 0.5*flow%denIn*flow%Uref**2*flow%Asfac*flow%Lref
        flow%Pref = 0.5*flow%denIn*flow%Uref**2*flow%Asfac*flow%Uref
        ! fluid viscosity
        flow%nu =  flow%Uref*flow%Lref/flow%Re
        flow%Mu =  flow%nu*flow%denIn
    END SUBROUTINE

end module FlowCondition