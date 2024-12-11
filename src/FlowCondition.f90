module FlowCondition
    implicit none
    private
    public :: conditions
    type :: FlowCondType
        integer :: TrefType,UrefType
        real(8) :: Re,dt,denIn,nu,Mu
        real(8) :: uvwIn(1:3),shearRateIn(1:3)
        real(8) :: volumeForceIn(1:3),volumeForceAmp,volumeForceFreq,volumeForcePhi
        real(8) :: Uref,Lref,Tref
        integer :: isRelease,isConCmpt
        real(8) :: timeWriteBegin,timeWriteEnd,timeWriteFlow,timeWriteBody,timeWriteInfo
        real(8) :: VelocityAmp,VelocityFreq,VelocityPhi
    end type FlowCondType
    type(FlowCondType) :: conditions
    contains

    SUBROUTINE read_FlowCondition()
        ! read computing core and block numbers
        implicit none
        open(unit=111, file='inFlow.dat', status='old', action='read')
        call found_keyword(111,'FlowCondition')
        read(111,*)    conditions%Re,conditions%dt,conditions%denIn
        read(111,*)    conditions%uvwIn(1:3)
        read(111,*)    conditions%shearRateIn(1:3)
        read(111,*)    conditions%volumeForceIn(1:3)
        read(111,*)    conditions%volumeForceAmp,conditions%volumeForceFreq,conditions%volumeForcePhi
        read(111,*)    conditions%TrefType,conditions%Tref
        read(111,*)    conditions%UrefType,conditions%Uref
        read(111,*)    conditions%isRelease,conditions%isConCmpt
        read(111,*)    conditions%timeWriteBegin,conditions%timeWriteEnd
        read(111,*)    conditions%timeWriteFlow,conditions%timeWriteBody,conditions%timeWriteInfo

        conditions%nu = conditions%Re !needed

        close(111)
    END SUBROUTINE

end module FlowCondition