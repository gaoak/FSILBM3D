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
        character(LEN=40) :: keyword
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
            write(*,*) trim(keyword)//' is not found in inFlow.dat'
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
            read(ifile, '(a)', IOSTAT=IOstatus) buffer
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
    
    SUBROUTINE grid_value_interpolation(dh,xmin,ymin,zmin,xDim,yDim,zDim,Coord,valueIn,valueOut)
        implicit none
        integer:: xDim,yDim,zDim
        integer:: x1,y1,z1,x2,y2,z2
        real(8):: dh,xmin,ymin,zmin,Coord(3)
        real(8):: dx1,dy1,dz1,dx2,dy2,dz2
        real(8):: coffe,c1,c2,c3,c4,c5,c6,c7,c8
        real(8):: valueIn(zDim,yDim,xDim),valueOut
        !real(8):: interCoord(3)
        ! coordinate number of the surrounding grid points
        x1 = FLOOR((Coord(1)-xmin)/dh + 1)
        y1 = FLOOR((Coord(2)-ymin)/dh + 1)
        z1 = FLOOR((Coord(3)-zmin)/dh + 1)
        x2 = x1 + 1
        y2 = y1 + 1
        z2 = z1 + 1
        ! coordinate difference between the point and the surrounding grid points
        dx1 = Coord(1) - (xmin + dh*(x1 - 1))
        dy1 = Coord(2) - (ymin + dh*(y1 - 1))
        dz1 = Coord(3) - (zmin + dh*(z1 - 1))
        dx2 = dh - dx1
        dy2 = dh - dy1
        dz2 = dh - dz1
        ! interpolation coefficient
        coffe = 1.0d0/(dh*dh*dh)
        c1 = dx2*dy2*dz2*coffe
        c2 = dx2*dy2*dz1*coffe
        c3 = dx2*dy1*dz2*coffe
        c4 = dx1*dy2*dz2*coffe
        c5 = dx2*dy1*dz1*coffe
        c6 = dx1*dy2*dz1*coffe
        c7 = dx1*dy1*dz2*coffe
        c8 = dx1*dy1*dz1*coffe
        ! interpolation
        valueOut = c1*valueIn(z1,y1,x1) + c2*valueIn(z2,y1,x1) + c3*valueIn(z1,y2,x1) + c4*valueIn(z1,y1,x2) + &
                   c5*valueIn(z2,y2,x1) + c6*valueIn(z2,y1,x2) + c7*valueIn(z1,y2,x2) + c8*valueIn(z2,y2,x2)
        ! test interpolation coordinate
        !interCoord(1) = c1*(xmin + dh*(x1 - 1)) + c2*(xmin + dh*(x1 - 1)) + c3*(xmin + dh*(x1 - 1)) + c4*(xmin + dh*(x2 - 1)) + &
        !                c5*(xmin + dh*(x1 - 1)) + c6*(xmin + dh*(x2 - 1)) + c7*(xmin + dh*(x2 - 1)) + c8*(xmin + dh*(x2 - 1))
        !interCoord(2) = c1*(ymin + dh*(y1 - 1)) + c2*(ymin + dh*(y1 - 1)) + c3*(ymin + dh*(y2 - 1)) + c4*(ymin + dh*(y1 - 1)) + &
        !                c5*(ymin + dh*(y2 - 1)) + c6*(ymin + dh*(y1 - 1)) + c7*(ymin + dh*(y2 - 1)) + c8*(ymin + dh*(y2 - 1))
        !interCoord(3) = c1*(zmin + dh*(z1 - 1)) + c2*(zmin + dh*(z2 - 1)) + c3*(zmin + dh*(z1 - 1)) + c4*(zmin + dh*(z1 - 1)) + &
        !                c5*(zmin + dh*(z2 - 1)) + c6*(zmin + dh*(z2 - 1)) + c7*(zmin + dh*(z1 - 1)) + c8*(zmin + dh*(z2 - 1))
    END SUBROUTINE

    SUBROUTINE  write_parameter_check_file(filename)
        use FlowCondition
        use FluidDomain
        use LBMBlockComm
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
        call write_parameter_blocks(blockTreeRoot)
        close(111)

        contains

        recursive subroutine write_parameter_blocks(treenode)
            implicit none
            integer:: treenode,i,j,s
            character(len=4):: IDstr
            write(IDstr,'(I4.4)')LBMblks(treenode)%ID
            write(111,'(A,A,A  )')'============================ BlockID = ',IDstr,' ============================='
            write(111,'(A,F20.10)')'Tau  =', LBMblks(treenode)%tau
            write(111,'(A,F20.10)')'Omega =', LBMblks(treenode)%Omega
            write(111,'(A,A    )')'---------------------------------------------------------------------------'
            write(IDstr,'(I4.4)')LBMblks(treenode)%carriedBodies(0)
            write(111,'(A,A    )') 'carryBodies : ', IDstr
            write(111,'(A      )', advance='no') 'bodyIDs :'
            do i=1,LBMblks(treenode)%carriedBodies(0)
                write(IDstr,'(I4.4)')LBMblks(treenode)%carriedBodies(i)
                write(111,'(A,A    )', advance='no') '    ', IDstr
            enddo
            write(111,'(A      )') ''
            write(111,'(A      )', advance='no') 'sonBlock :'
            do j=1,blockTree(treenode)%nsons
                s = blockTree(treenode)%comm(j)%sonId
                write(IDstr,'(I4.4)')s
                write(111,'(A,A    )', advance='no') '    ', IDstr
                write(111,'(A      )') ''
                call write_parameter_blocks(s)
            enddo
            end subroutine
    END SUBROUTINE