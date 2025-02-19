module FlowCondition
    implicit none
    private
    public :: FlowCondType,flow
    public :: read_flow_conditions,read_probe_params,write_information_titles,write_fluid_information
    type :: FlowCondType
        integer :: isConCmpt,numsubstep,npsize
        real(8) :: timeSimTotal,timeContiDelta,timeWriteBegin,timeWriteEnd,timeFlowDelta,timeBodyDelta,timeInfoDelta
        real(8) :: Re,denIn,nu,Mu,dtolLBM
        integer :: TrefType,UrefType,ntolLBM
        integer :: velocityKind,interpolateScheme
        real(8) :: uvwIn(1:3),shearRateIn(1:3)
        real(8) :: volumeForceIn(1:3),volumeForceAmp,volumeForceFreq,volumeForcePhi
        real(8) :: Uref,Lref,Tref
        real(8) :: Aref,Fref,Eref,Pref
        real(8) :: Asfac,Lchod,Lspan,AR 
        integer :: fluidProbingNum,inWhichBlock,solidProbingNum
        integer, allocatable :: solidProbingNode(:)
        real(8), allocatable :: fluidProbingCoords(:,:)
        real(8) :: AmplInitDist(1:3),waveInitDist
    end type FlowCondType
    type(FlowCondType) :: flow

    contains

    SUBROUTINE read_flow_conditions(filename)
        implicit none
        character(LEN=40),intent(in):: filename
        character(LEN=256):: buffer
        character(LEN=40):: keywordstr
        open(unit=111, file=filename, status='old', action='read')
        keywordstr = 'Parallel'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*)    flow%npsize
        rewind(111)
        keywordstr = 'FlowCondition'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*)    flow%isConCmpt,flow%numsubstep
        call readNextData(111, buffer)
        read(buffer,*)    flow%timeSimTotal,flow%timeContiDelta
        call readNextData(111, buffer)
        read(buffer,*)    flow%timeWriteBegin,flow%timeWriteEnd
        call readNextData(111, buffer)
        read(buffer,*)    flow%timeFlowDelta,flow%timeBodyDelta,flow%timeInfoDelta
        call readNextData(111, buffer)
        read(buffer,*)    flow%Re,flow%denIn
        call readNextData(111, buffer)
        read(buffer,*)    flow%uvwIn(1:3)
        call readNextData(111, buffer)
        read(buffer,*)    flow%shearRateIn(1:3),flow%velocityKind
        call readNextData(111, buffer)
        read(buffer,*)    flow%volumeForceIn(1:3)
        call readNextData(111, buffer)
        read(buffer,*)    flow%volumeForceAmp,flow%volumeForceFreq,flow%volumeForcePhi
        call readNextData(111, buffer)
        read(buffer,*)    flow%TrefType,flow%Tref
        call readNextData(111, buffer)
        read(buffer,*)    flow%UrefType,flow%Uref
        call readNextData(111, buffer)
        read(buffer,*)    flow%ntolLBM,flow%dtolLBM
        call readNextData(111, buffer)
        read(buffer,*)    flow%interpolateScheme
        close(111)
        ! flow%denIn is not 1
        if(abs(flow%denIn-1.d0).gt.1e-6) then
            write(*,*) 'Warning, denIn is not 1, ', flow%denIn
        endif
    END SUBROUTINE

    SUBROUTINE read_probe_params(filename)
        implicit none
        integer:: i
        character(LEN=40),intent(in):: filename
        character(LEN=256):: buffer
        character(LEN=40):: keywordstr
        ! read fluid probes
        open(unit=111, file=filename, status='old', action='read')
        keywordstr = 'ProbingFluid'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*)    flow%fluidProbingNum,flow%inWhichBlock
        if (flow%fluidProbingNum .ne. 0) then
            allocate(flow%fluidProbingCoords(flow%fluidProbingNum,1:3))
            do i=1,flow%fluidProbingNum
                call readNextData(111, buffer)
                read(buffer,*)    flow%fluidProbingCoords(i,1:3)
            enddo
        endif
        close(111)
        ! read solid probes
        open(unit=111, file=filename, status='old', action='read')
        keywordstr = 'ProbingSolid'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*)    flow%solidProbingNum
        if (flow%solidProbingNum .ne. 0) then
            allocate(flow%solidProbingNode(flow%solidProbingNum))
            do i=1,flow%solidProbingNum
                call readNextData(111, buffer)
                read(buffer,*)    flow%solidProbingNode(i)
            enddo
        endif
        close(111)
    END SUBROUTINE

    ! SUBROUTINE read_distb_params(filename)
    !     implicit none
    !     character(LEN=40),intent(in):: filename
    !     character(LEN=256):: buffer
    !     character(LEN=40):: keywordstr
    !     ! read fluid disturb
    !     open(unit=111, file=filename, status='old', action='read')
    !     keywordstr = 'Disturb'
    !     call found_keyword(111,keywordstr)
    !     call readNextData(111, buffer)
    !     read(buffer,*)    flow%waveInitDist,flow%AmplInitDist(1:3)
    !     close(111)
    ! END SUBROUTINE

    SUBROUTINE write_information_titles(nFish)
        implicit none
        integer:: i,iFish,nFish
        integer,parameter::nameLen=3
        character (LEN=nameLen):: fishNum,probeNum
        ! fish informations
        do iFish=1,nFish
            ! get fish numbers
            write(fishNum,'(I3)') iFish
            fishNum = adjustr(fishNum)
            do  i=1,nameLen
                    if(fishNum(i:i)==' ') fishNum(i:i)='0'
            enddo
            ! write forces title
            open(111,file='./DatInfo/FishForce_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t"  "Fx"  "Fy"  "Fz"'
            close(111)
            ! write begin information titles
            open(111,file='./DatInfo/FishNodeBegin_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az"'
            close(111)
            ! write end information titles
            open(111,file='./DatInfo/FishNodeEnd_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az"'
            close(111)
            ! write center information titles
            open(111,file='./DatInfo/FishNodeCenter_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az"'
            close(111)
            ! write mean information titles
            open(111,file='./DatInfo/FishNodeMean_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az"'
            close(111)
            ! write angular information titles
            open(111,file='./DatInfo/FishAngular_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t"  "AoA"  "Ty-Hy"  "Hy"  "Ty"'
            close(111) 
            ! write power title
            open(111,file='./DatInfo/FishPower_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t" "Ptot" "Paero" "Piner" "Pax" "Pay" "Paz" "Pix" "Piy" "Piz"'
            close(111)
            ! write area title
            ! open(111,file='./DatInfo/FishArea_'//trim(fishNum)//'.plt')
            ! write(111,*) 'variables= "t"  "Area"'
            ! close(111)
            ! write energy title
            open(111,file='./DatInfo/FishEnergy_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t","Es","Eb","Ep","Ek","Ew","Et"'
            close(111)
            ! write solid probing title
            do  i=1,flow%solidProbingNum
                write(probeNum,'(I3.3)') i
                open(111,file='./DatInfo/FishProbes_'//trim(fishNum)//'_'//trim(probeNum)//'.plt')
                write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az"'
                close(111)
            enddo
        enddo
        ! write fluid probing title
        do  i=1,flow%fluidProbingNum
            write(probeNum,'(I3.3)') i
            open(111,file='./DatInfo/FluidProbes_'//trim(probeNum)//'.plt')
            write(111,*) 'variables= "t"  "u"  "v"  "w" '
            close(111)
        enddo
        ! ! write max vel of fluid title
        ! open(111,file='./DatInfo/MaMax.plt')
        ! write(111,*)'variables= "t"  "MaMax"  '
        ! close(111)
        ! ! write convergence title
        ! open(111,file='./DatInfo/Converg.plt')
        ! write(111,*)'variables= "t"  "Convergence"  '
        ! close(111)
    END SUBROUTINE

    SUBROUTINE write_fluid_information(time,dh,xmin,ymin,zmin,xDim,yDim,zDim,velocityIn)
        implicit none
        integer:: i,j,xDim,yDim,zDim
        real(8):: time,dh,xmin,ymin,zmin,xmax,ymax,zmax
        real(8):: velocityIn(zDim,yDim,xDim,1:3),velocityOut(1:3)
        integer,parameter::nameLen=3
        character (LEN=nameLen):: probeNum
        ! write fluid probing information
        do  i=1,flow%fluidProbingNum
            xmax = xmin + dh * (xDim - 1)
            ymax = ymin + dh * (yDim - 1)
            zmax = zmin + dh * (zDim - 1)
            if(flow%fluidProbingCoords(i,1) .lt. xmin .or. flow%fluidProbingCoords(i,1) .gt. xmax .or. & 
               flow%fluidProbingCoords(i,2) .lt. ymin .or. flow%fluidProbingCoords(i,2) .gt. ymax .or. &
               flow%fluidProbingCoords(i,3) .lt. zmin .or. flow%fluidProbingCoords(i,3) .gt. zmax) then 
                write(*,'(A,I2,A)') 'fluid probe', i ,' is not in selected block'
                stop
            endif
            ! velocity interpolation
            do j=1,3
                call grid_value_interpolation(dh,xmin,ymin,zmin,xDim,yDim,zDim,flow%fluidProbingCoords(i,1:3),velocityIn(1:zDim,1:yDim,1:xDim,j),velocityOut(j))
            enddo
            ! write file
            write(probeNum,'(I3.3)') i
            open(111,file='./DatInfo/FluidProbes_'//trim(probeNum)//'.plt',position='append')
            write(111,'(4E20.10)') time/flow%Tref,velocityOut(1:3)/flow%Uref
            close(111)
        enddo
        END SUBROUTINE

end module FlowCondition
