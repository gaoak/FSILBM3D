!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write files' header
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_information_titles()
    implicit none
    integer:: i,iFish
    integer,parameter::nameLen=4
    character (LEN=nameLen):: fileName,Nodename

    do iFish=1,nFish

        write(fileName,'(I4)') iFish
        fileName = adjustr(fileName)
        do  i=1,nameLen
             if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo

        if    (iForce2Body==1)then   !Same force as flow
        open(111,file='./DatInfo/ForceDirect_'//trim(filename)//'.plt')
        write(111,*)'variables= "t"  "Fx"  "Fy"  "Fz"'
        close(111)
        elseif(iForce2Body==2)then   !stress force
        open(111,file='./DatInfo/ForceStress_'//trim(filename)//'.plt')
        write(111,*)'variables= "t"  "Fx"  "Fy"  "Fz"'
        close(111)
        endif

        !===============================================================================
        open(111,file='./DatInfo/SampBodyNodeBegin_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)
        !===============================================================================
        open(111,file='./DatInfo/SampBodyNodeEnd_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)

        open(111,file='./DatInfo/SampBodyNodeCenter_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "u"  "v"  "ax"  "ay" '
        close(111)

        open(111,file='./DatInfo/SampBodyMean_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)

        open(111,file='./DatInfo/SampBodyAngular_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "AoA"  "Ty-Hy"  "Hy"  "Ty"'
        close(111)

        open(111,file='./DatInfo/Power_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t" "Ptot" "Paero" "Piner" "Pax" "Pay" "Paz" "Pix" "Piy" "Piz"'
        close(111)

        open(111,file='./DatInfo/Area_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "Area"  '
        close(111)

        open(111,file='./DatInfo/Energy_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t","Es","Eb","Ep","Ek","Ew","Et"'
        close(111)
    enddo

    if(isBodyOutput==1)then
        do iFish=1,nFish
            write(fileName,'(I4)') iFish
            fileName = adjustr(fileName)
            do  i=1,nameLen
                if(fileName(i:i)==' ')fileName(i:i)='0'
            enddo
            do  i=1,numSampBody
                write(Nodename,'(I4.4)') SampBodyNode(i,iFish)
                open(111,file='./DatInfo/SampBodyNode_'//trim(fileName)//'_'//trim(Nodename)//'.plt')
                write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
                close(111)
            enddo
        enddo !iFish
    endif

    if(isFluidOutput==1)then
    do  i=1,numSampFlow
        write(Nodename,'(I4.4)') i
        open(111,file='./DatInfo/SampFlowPint_'//trim(Nodename)//'.plt')
        write(111,*)'variables= "t"  "p" "u"  "v"  "w" '
        close(111)
    enddo
    endif
    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     write data
!    copyright@ RuNanHua
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wrtInfo()
    USE SolidBody
    implicit none
    integer:: i,z,y,x,zbgn,ybgn,xbgn,zend,yend,xend
    real(8):: weightm,velocity(1:3),Pressure
    integer,parameter::nameLen=4
    character (LEN=nameLen):: fileName

    call Write_solid_Data(111,time,timeOutInfo,Asfac)

    if(isBodyOutput==1)then
        call Write_SampBodyNode(111,time,numSampBody,SampBodyNode)
    endif

    if(isFluidOutput==1)then
    do  i=1,numSampFlow
        x=minloc(dabs(SampFlowPint(i,1)-xGrid(1:xDim)),1)
        if(SampFlowPint(i,1)-xGrid(x)>0.0d0)then
                xbgn=x; xend=x+1
        else
                xbgn=x-1;xend=x
        endif
        if(xend.eq.xDim+1) then
            xbgn = xDim - 1
            xend = xDim
        else if(xbgn.eq.0) then
            xbgn = 1
            xend = 2
        endif
        y=minloc(dabs(SampFlowPint(i,2)-yGrid(1:yDim)),1)
        if(SampFlowPint(i,2)-yGrid(y)>0.0d0)then
                ybgn=y; yend=y+1
        else
                ybgn=y-1;yend=y
        endif
        if(yend.eq.yDim+1) then
            ybgn = yDim - 1
            yend = yDim
        else if(ybgn.eq.0) then
            ybgn = 1
            yend = 2
        endif
        z=minloc(dabs(SampFlowPint(i,3)-zGrid(1:zDim)),1)
        if(SampFlowPint(i,3)-zGrid(z)>0.0d0)then
                zbgn=z; zend=z+1
        else
                zbgn=z-1;zend=z
        endif
        if(zend.eq.zDim+1) then
            zbgn = zDim - 1
            zend = zDim
        else if(zbgn.eq.0) then
            zbgn = 1
            zend = 2
        endif

        velocity(1:3)=0.0d0
        Pressure=0.0d0

        do    x=xbgn,xend
        do    y=ybgn,yend
        do    z=zbgn,zend
            weightm=(dx(x)-dabs(xGrid(x)-SampFlowPint(i,1)))*(dy(y)-dabs(yGrid(y)-SampFlowPint(i,2)))*(dz(z)-dabs(zGrid(z)-SampFlowPint(i,3)))
            velocity(1:3)=velocity(1:3)+uuu(z,y,x,1:3)*weightm/(dx(x)*dy(y)*dz(z))
            Pressure=Pressure+prs(z,y,x)*weightm/(dx(x)*dy(y)*dz(z))
        enddo
        enddo
        enddo

        write(fileName,'(I4.4)') i
        open(111,file='./DatInfo/SampFlowPint'//trim(fileName)//'.plt',position='append')
        write(111,'(5E20.10)') time/Tref,  Pressure/(0.5*denIn*Uref**2), velocity(1:3)/Uref
        close(111)
    enddo
    endif
    END SUBROUTINE

    subroutine cptArea(areaElem,nND,nEL,ele,xyzful)
    implicit none
    integer:: nND,nEL,ele(nEL,5)
    real(8):: areaElem(nEL),xyzful(nND,6)
    integer:: i,j,k,nt,iEL
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az
        do  iEL=1,nEL
        i =ele(iEL,1)
        j =ele(iEL,2)
        k =ele(iEL,3)
        nt=ele(iEL,4)

        x1=xyzful(i,1)
        x2=xyzful(j,1)
        x3=xyzful(k,1)
        y1=xyzful(i,2)
        y2=xyzful(j,2)
        y3=xyzful(k,2)
        z1=xyzful(i,3)
        z2=xyzful(j,3)
        z3=xyzful(k,3)

        if(nt==2)then
        ax = x2-x1
        ay = y2-y1
        az = z2-z1
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az)
        endif

        if(nt==3)then
        ax =(z1-z2)*(y3-y2) + (y2-y1)*(z3-z2)
        ay =(x1-x2)*(z3-z2) + (z2-z1)*(x3-x2)
        az =(y1-y2)*(x3-x2) + (x2-x1)*(y3-y2)
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az) * 0.5d0
        endif

    enddo
    endsubroutine
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AoAtoTTT(AoA,TTT)
    implicit none
    real(8), parameter:: pi=3.141562653589793d0,eps=1.0d-5
    real(8)::AoA(3),TTT(3,3),rrx(3,3),rry(3,3),rrz(3,3),vcos,vsin
    TTT(:,:)=0.0d0
    TTT(1,1)=1.0d0
    TTT(2,2)=1.0d0
    TTT(3,3)=1.0d0

    vcos=dcos(AoA(1))
    vsin=dsin(AoA(1))
    if(dabs(AoA(1))<eps)then
    vcos=1.0d0
    vsin=0.0d0
    endif
    if(dabs(AoA(1)-0.5d0*pi)<eps)then
    vcos=0.0d0
    vsin=1.0d0
    endif
    if(dabs(AoA(1)+0.5d0*pi)<eps)then
    vcos=0.0d0
    vsin=-1.0d0
    endif

    rrx(1:3,1:3)=reshape([  1.0d0,0.0d0,0.0d0,  &
                            0.0d0,vcos,vsin, &
                            0.0d0,-vsin,vcos],[3,3])
    vcos=dcos(AoA(2))
    vsin=dsin(AoA(2))
    if(dabs(AoA(2))<eps)then
    vcos=1.0d0
    vsin=0.0d0
    endif
    if(dabs(AoA(2)-0.5d0*pi)<eps)then
    vcos=0.0d0
    vsin=1.0d0
    endif
    if(dabs(AoA(2)+0.5d0*pi)<eps)then
    vcos=0.0d0
    vsin=-1.0d0
    endif

    rry(1:3,1:3)=reshape([  vcos,0.0d0,-vsin,  &
                            0.0d0,1.0d0,0.0d0, &
                           vsin,0.0d0,vcos],[3,3])

    vcos=dcos(AoA(3))
    vsin=dsin(AoA(3))
    if(dabs(AoA(3))<eps)then
    vcos=1.0d0
    vsin=0.0d0
    endif
    if(dabs(AoA(3)-0.5d0*pi)<eps)then
    vcos=0.0d0
    vsin=1.0d0
    endif
    if(dabs(AoA(3)+0.5d0*pi)<eps)then
    vcos=0.0d0
    vsin=-1.0d0
    endif

    rrz(1:3,1:3)=reshape([ vcos,vsin,0.0d0, &
                          -vsin,vcos,0.0d0, &
                           0.0d0,0.0d0,1.0d0],[3,3])
    TTT=matmul(rrz,TTT)
    TTT=matmul(rry,TTT)
    TTT=matmul(rrx,TTT)

    end subroutine
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write  string to unit=idfile  in binary format
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dumpstring(instring,idfile)
    implicit none
    character(40) instring
    integer:: nascii,ii,len,idfile
!
    len=LEN_TRIM(instring)

    do    ii=1,len
        nascii=ICHAR(instring(ii:ii))
        write(idfile) nascii
    enddo

    write(idfile) 0
!
    return
    endsubroutine dumpstring

    FUNCTION CPUtime(values)
        IMPLICIT NONE
        real(8)::CPUtime
        integer,dimension(8) :: values
        CPUtime = dble(values(6))*60.d0+dble(values(7))*1.d0+dble(values(8))*0.001d0
    ENDFUNCTION

    ! found keyword in inflow.dat for next parameters read
    SUBROUTINE found_keyword(fileID,keyword)
        implicit none
        integer fileID
        character :: keyword
        character(len=50) :: readString
        readString = 'null'
        do while(.not.EOF(fileID))
            read(fileID, *) readString
            if (index(LOWCASE(readString), LOWCASE(keyword)) .GT. 0) then
                exit
            endif
        enddo
        if (index(LOWCASE(readString), LOWCASE(keyword)) .EQ. 0) then
            write(*,*) 'the parameters in inflow.dat do not exist.'
            stop
        endif
    END SUBROUTINE

    subroutine readNextData(ifile, buffer)
        implicit none
        character(LEN=256):: buffer
        integer:: ifile, IOstatus
        do while(.true.)
            read(ifile, *, IOSTAT=IOstatus) buffer
            if (IOstatus.ne.0) then
                write(*, *) 'end of file encounter in readNextData', ifile
                stop
            endif
            adjustl(buffer)
            if(buffer(1:1).ne.'#') then
                exit
            endif
        enddo
    end subroutine readNextData

    SUBROUTINE readequal(ifile)
        implicit none
        integer:: IOstatus, ifile
        character (40):: buffer
        do while(.true.)
            read(ifile, *, IOSTAT=IOstatus) buffer
            if (IOstatus.ne.0) then
                write(*, *) 'end of file encounter in readequal', ifile
                stop
            endif
            adjustl(buffer)
            if(buffer(1:1).eq.'=') then
                exit
            endif
        enddo
    END SUBROUTINE

    SUBROUTINE check_is_continue(filename,step,time,isContinue)
        implicit none
        integer:: isContinue
        character(LEN=40),intent(in):: filename
        inquire(file=filename, exist=alive)
        if (isContinue==1 .and. alive)then
            write(*,*)  '=============================================='
            write(*,*)  '==============Continue computing=============='
            write(*,*)  '=============================================='
            call read_continue_blocks(filename,step,time)
        else
            write(*,*)  '=============================================='
            write(*,*)  '=================New computing================'
            write(*,*)  '=============================================='
        endif
    END SUBROUTINE