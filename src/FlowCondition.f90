module FlowCondition
    implicit none
    private
    public :: FlowCondType,flow
    public :: read_flow_conditions,write_parameter_check_file,read_probe_params,write_information_titles,write_fluid_information,write_solid_information
    type :: FlowCondType
        integer :: isConCmpt,numsubstep,npsize
        real(8) :: timeSimTotal,timeContiDelta,timeWriteBegin,timeWriteEnd,timeFlowDelta,timeBodyDelta,timeInfoDelta
        real(8) :: Re,denIn,nu,Mu,dtolLBM
        integer :: TrefType,UrefType,ntolLBM
        integer :: velocityKind
        real(8) :: uvwIn(1:3),shearRateIn(1:3)
        real(8) :: volumeForceIn(1:3),volumeForceAmp,volumeForceFreq,volumeForcePhi
        real(8) :: Uref,Lref,Tref
        real(8) :: Aref,Fref,Eref,Pref
        real(8) :: Asfac,Lchod,Lspan,AR 
        integer :: fluidProbingNum,inWhichBlock,solidProbingNum
        integer, allocatable :: solidProbingNode(:)
        real(8), allocatable :: fluidProbingCoords(:,:)
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
        close(111)
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

    SUBROUTINE  write_parameter_check_file(filename)
        implicit none
        character(LEN=40):: filename
        open(111,file=filename)
        write(111,'(A      )')'===================================================================='
        write(111,'(A,F20.10)')'Re   =', flow%Re
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
            ! write average information titles
            open(111,file='./DatInfo/FishMeanInfo_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v" '
            close(111)
            ! write power title
            open(111,file='./DatInfo/FishPower_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t" "Ptot" "Paero" "Piner" "Pax" "Pay" "Paz" "Pix" "Piy" "Piz"'
            close(111)
            ! write energy title
            open(111,file='./DatInfo/FishEnergy_'//trim(fishNum)//'.plt')
            write(111,*) 'variables= "t","Es","Eb","Ep","Ek","Ew","Et"'
            close(111)
            ! probing informations
            do  i=1,flow%solidProbingNum
                write(probeNum,'(I3.3)') i
                open(111,file='./DatInfo/FishProbes_'//trim(fishNum)//'_'//trim(probeNum)//'.plt')
                write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w" '
                close(111)
            enddo
        enddo
        ! write fluid probing title
        do  i=1,flow%fluidProbingNum
            write(probeNum,'(I3.3)') i
            open(111,file='./DatInfo/FluidProbes_'//trim(probeNum)//'.plt')
            write(111,*) 'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w" '
            close(111)
        enddo
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
            write(111,'(4E20.10)') time/flow%Tref,flow%fluidProbingCoords(i,1:3)/flow%Lref,velocityOut(1:3)/flow%Uref
            close(111)
        enddo
        END SUBROUTINE

        SUBROUTINE write_solid_information(time,nFish)
            implicit none
            integer:: i,iFish,nFish
            real(8):: time
            integer,parameter::nameLen=3
            character (LEN=nameLen):: fishNum,probeNum
        
            do iFish=1,nFish
                ! get fish numbers
                write(fishNum,'(I3)') iFish
                fishNum = adjustr(fishNum)
                do  i=1,nameLen
                        if(fishNum(i:i)==' ') fishNum(i:i)='0'
                enddo
                ! write forces
                open(111,file='./DatInfo/FishForce_'//trim(fishNum)//'.plt',position='append')
                write(111,*) 'test'
                !write(111,'(4E20.10)') time/flow%Tref,sum(extful(1:nND(iFish),1:3,iFish),1)/Fref
                close(111)
                ! write average information
                open(111,file='./DatInfo/FishMeanInfo_'//trim(fishNum)//'.plt',position='append')
                write(111,*) 'test'
                !write(111,'(10E20.10)')time/flow%Tref,sum(xyzful(1:nND(iFish),1:3,iFish)*mssful(1:nND(iFish),1:3,iFish),1)/sum(mssful(1:nND(iFish),1:3,iFish),1)/Lref, &
                !                                      sum(velful(1:nND(iFish),1:3,iFish)*mssful(1:nND(iFish),1:3,iFish),1)/sum(mssful(1:nND(iFish),1:3,iFish),1)/Uref
                close(111)
                ! write power
                open(111,file='./DatInfo/FishPower_'//trim(fishNum)//'.plt',position='append')
                !Pax=sum(extful(1:nND(iFish),1,iFish)*velful(1:nND(iFish),1,iFish))/Pref
                !Pay=sum(extful(1:nND(iFish),2,iFish)*velful(1:nND(iFish),2,iFish))/Pref
                !Paz=sum(extful(1:nND(iFish),3,iFish)*velful(1:nND(iFish),3,iFish))/Pref
                !Pix=-sum(mssful(1:nND(iFish),1,iFish)*accful(1:nND(iFish),1,iFish)*velful(1:nND(iFish),1,iFish))/Pref
                !Piy=-sum(mssful(1:nND(iFish),2,iFish)*accful(1:nND(iFish),2,iFish)*velful(1:nND(iFish),2,iFish))/Pref
                !Piz=-sum(mssful(1:nND(iFish),3,iFish)*accful(1:nND(iFish),3,iFish)*velful(1:nND(iFish),3,iFish))/Pref
                !Paero=Pax+Pay+Paz
                !Piner=Pix+Piy+Piz
                !Ptot=Paero+Piner
                !write(111,'(10E20.10)')time/flow%Tref,Ptot,Paero,Piner,Pax,Pay,Paz,Pix,Piy,Piz
                write(111,*) 'test'
                close(111)
                ! write energy title
                open(111,file='./DatInfo/FishEnergy_'//trim(fishNum)//'.plt', position='append')
                !call strain_energy_D(strainEnergy(1:nEL(iFish),1:2,iFish),xyzful0(1:nND(iFish),1,iFish),xyzful0(1:nND(iFish),2,iFish),xyzful0(1:nND(iFish),3,iFish), &
                !                        xyzful(1:nND(iFish),1,iFish), xyzful(1:nND(iFish),2,iFish), xyzful(1:nND(iFish),3,iFish),ele(1:nEL(iFish),1:5,iFish), prop(1:nMT(iFish),1:10,iFish), &
                !                        triad_n1(1:3,1:3,1:nEL(iFish),iFish),triad_n2(1:3,1:3,1:nEL(iFish),iFish), &
                !                        triad_ee(1:3,1:3,1:nEL(iFish),iFish), &
                !                        nND(iFish),nEL(iFish),nMT(iFish))
                !EEE(1)=sum(strainEnergy(1:nEL(iFish),1,iFish))
                !EEE(2)=sum(strainEnergy(1:nEL(iFish),2,iFish))
                !Es=EEE(1)/Eref
                !Eb=EEE(2)/Eref
                !Ep=Es+Eb
                !Ew=Ew+Paero*timeOutInfo
                !Ek=0.5*sum(mssful(1:nND(iFish),1:6,iFish)*velful(1:nND(iFish),1:6,iFish)*velful(1:nND(iFish),1:6,iFish))/Eref
                !Et=Ek+Ep
                !write(111,'(7E20.10)')time/flow%Tref,Es,Eb,Ep,Ek,Ew,Et
                write(111,*) 'test'
                close(111)
                
                ! probing informations
                do  i=1,flow%solidProbingNum
                    write(probeNum,'(I3.3)') i
                    open(111,file='./DatInfo/FishProbes_'//trim(fishNum)//'_'//trim(probeNum)//'.plt',position='append')
                !    write(111,'(5E20.10)') time/flow%Tref,  Pressure/(0.5*denIn*Uref**2), velocity(1:3)/Uref
                    write(111,*) 'test'
                    close(111)
                enddo
            enddo   
        END SUBROUTINE

end module FlowCondition
