module FlowCondition
    implicit none
    private
    public:: flow
    type :: FlowCondType
        real(8):: denIn,Uref,Lref,Tref,nu,Mu
        real(8):: uuuIn(1:3),shearRateIn(1:3),VelocityAmp,VelocityFreq,VelocityPhi
        real(8):: VolumeForce(1:3),VolumeForceAmp,VolumeForceFreq,VolumeForcePhi,VolumeForceIn(1:3)
    end type FlowCondType
    type(FlowCondType):: flow

    contains

    subroutine read_flowCondition()
        implicit none
        ! to do
        ! read inflow.dat FlowCondition to flow
    end subroutine

end module FlowCondition