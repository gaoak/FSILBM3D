module FlowCondition
    implicit none
    private
    public :: FlowCondType,flow
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
        character(LEN=256):: buffer
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

end module FlowCondition