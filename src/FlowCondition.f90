module FlowCondition
    implicit none
    private
    public :: FlowCondType,flow
    public :: read_flow_conditions,write_parameter_check_file
    type :: FlowCondType
        integer :: isConCmpt,numsubstep
        real(8) :: timeSimTotal,timeContiDelta,timeWriteBegin,timeWriteEnd,timeFlowDelta,timeBodyDelta,timeInfoDelta
        real(8) :: Re,dt,denIn,nu,Mu,dtolLBM
        integer :: TrefType,UrefType,ntolLBM
        real(8) :: uvwIn(1:3),shearRateIn(1:3)
        real(8) :: volumeForceIn(1:3),volumeForceAmp,volumeForceFreq,volumeForcePhi
        real(8) :: Uref,Lref,Tref
        real(8) :: Aref,Fref,Eref,Pref
        real(8) :: Asfac,Lchod,Lspan,AR 
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
        read(111,*)    flow%isConCmpt,flow%numsubstep
        call readNextData(111, buffer)
        read(111,*)    flow%timeSimTotal,flow%timeContiDelta
        call readNextData(111, buffer)
        read(111,*)    flow%timeWriteBegin,flow%timeWriteEnd
        call readNextData(111, buffer)
        read(111,*)    flow%timeFlowDelta,flow%timeBodyDelta,flow%timeInfoDelta
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

    SUBROUTINE  write_parameter_check_file(filename)
        implicit none
        character(LEN=40):: filename
        open(111,file=filename)
        write(111,'(A      )')'===================================================================='
        write(111,'(A,F20.10)')'Re   =', flow%Re
        write(111,'(A,F20.10)')'dt   =', flow%dt
        write(111,'(A,F20.10)')'den  =', flow%denIn
        write(111,'(A      )')'===================================================================='
        write(111,'(A,F20.10)')'Nu   =', flow%Nu
        write(111,'(A,F20.10)')'Mu   =', flow%Mu
        write(111,'(A      )')'===================================================================='
        write(111,'(A,F20.10)')'Lref =', flow%Lref
        write(111,'(A,F20.10)')'Uref =', flow%Uref
        write(111,'(A,F20.10)')'Tref =', flow%Tref
        write(111,'(A,F20.10)')'Aref =', flow%Aref
        write(111,'(A,F20.10)')'Pref =', flow%Pref
        write(111,'(A,F20.10)')'Eref =', flow%Eref
        write(111,'(A,F20.10)')'Fref =', flow%Fref
        write(111,'(A      )')'===================================================================='
        close(111)
    END SUBROUTINE
end module FlowCondition