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
!   copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dumpstring(instring,idfile)
    implicit none
    character(40) instring
    integer:: nascii,ii,len,idfile

    len=LEN_TRIM(instring)

    do    ii=1,len
        nascii=ICHAR(instring(ii:ii))
        write(idfile) nascii
    enddo

    write(idfile) 0

    return
    endsubroutine dumpstring

    ! get the time right now
    SUBROUTINE get_now_time(now_time) 
        IMPLICIT NONE
        real(8)::now_time
        integer,dimension(8) :: cpu_time
        call date_and_time(VALUES=cpu_time)
        now_time = cpu_time(6)*60.d0 + cpu_time(7)*1.d0 + cpu_time(8)*0.001d0
    END SUBROUTINE

    ! Found keyword in inflow.dat for next parameters read
    SUBROUTINE found_keyword(fileID,keyword)
        implicit none
        integer:: fileID, IOstatus
        character :: keyword
        character(len=50) :: readString
        readString = 'null'
        IOstatus = 0
        call to_lowercase(keyword)
        do while(IOstatus.eq.0)
            read(fileID, *, IOSTAT=IOstatus) readString
            call to_lowercase(readString)
            readString = adjustl(readString)
            if (index(readString, keyword) .GT. 0) then
                exit
            endif
        enddo
        if (index(readString, keyword) .EQ. 0) then
            write(*,*) 'the parameters in inflow.dat do not exist.'
            stop
        endif
    END SUBROUTINE

    ! Convert an uppercase string to lowercase
    SUBROUTINE to_lowercase(string)
        implicit none
        character:: string
        integer :: i
        do i = 1, len_trim(string)
            if (iachar(string(i:i)) .ge. iachar('A') .and. iachar(string(i:i)).le. iachar('Z')) then
                string(i:i) = achar(iachar(string(i:i)) + iachar('a') - iachar('A'))
            end if
        end do
    END SUBROUTINE

    SUBROUTINE readNextData(ifile, buffer)
        implicit none
        character(LEN=256):: buffer
        integer:: ifile, IOstatus
        do while(.true.)
            read(ifile, *, IOSTAT=IOstatus) buffer
            if (IOstatus.ne.0) then
                write(*, *) 'end of file encounter in readNextData', ifile
                stop
            endif
            buffer = adjustl(buffer)
            if(buffer(1:1).ne.'#') then
                exit
            endif
        enddo
    END SUBROUTINE readNextData

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
            buffer = adjustl(buffer)
            if(buffer(1:1).eq.'=') then
                exit
            endif
        enddo
    END SUBROUTINE