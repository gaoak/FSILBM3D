!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    write field, tecplot binary format
!    copyright@ RuNanHua       
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_flow_field(isT)
    USE simParam
    implicit none
    integer:: x,y,z,i,isT
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName
    !==================================================================================================        
    integer::    nv
    integer,parameter:: namLen=40,idfile=100,numVar=7 
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=['x','y','z','p','u','v','w']  

    write(fileName,'(F10.5)') time/Tref
    fileName = adjustr(fileName)
    do  i=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    enddo
    if(isT==0)then
        open(idfile,file='./DatFlow/Flow.plt',form='unformatted',access='stream')
    else
        open(idfile,file='./DatFlow/Flow'//trim(fileName)//'.plt',form='unformatted',access='stream')
    endif
    
    
!   =============================================================================================
    write(idfile) "#!TDV101"    
    write(idfile) 1
    call dumpstring(title,idfile)

    write(idfile) numVar
    do  nv=1,numVar
        call dumpstring(varname(nv),idfile)
    enddo
!   ===============================================================================================
    write(idfile) ZONEMARKER              
    call dumpstring(zonename,idfile)
    write(idfile) -1,0,1,0,0,zDim-2*numOutput,ydim-2*numOutput,xDim-2*numOutput,0

    write(idfile) ZONEMARKER              
    call dumpstring(zonename,idfile)
    write(idfile) -1,2,1,0,0,nND,nEL,0,0,0,0

    write(idfile) EOHMARKER

    write(idfile) ZONEMARKER
!    variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    
    if(isMoveGrid==1 .and. (isMoveOutputX==1 .or. isMoveOutputY==1 .or. isMoveOutputZ==1))then 

        if    ((isMoveDimX==1 .and. isMoveOutputX==1) .and. (isMoveDimY/=1 .and. isMoveOutputY/=1) .and. (isMoveDimZ/=1.and. isMoveOutputZ/=1))then
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput
            write(idfile)real((xGrid(x)-xyzful(NDref,1))/Lref),real(yGrid(y)/Lref),real(zGrid(z)/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                     real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
        elseif((isMoveDimX/=1 .and. isMoveOutputX/=1) .and. (isMoveDimY==1 .and. isMoveOutputY==1) .and. (isMoveDimZ/=1.and. isMoveOutputZ/=1))then
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput
            write(idfile)real(xGrid(x)/Lref),real((yGrid(y)-xyzful(NDref,2))/Lref),real(zGrid(z)/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                     real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
        elseif((isMoveDimX/=1 .and. isMoveOutputX/=1) .and. (isMoveDimY/=1 .and. isMoveOutputY/=1) .and. (isMoveDimZ==1.and. isMoveOutputZ==1))then
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput
            write(idfile)real(xGrid(x)/Lref),real(yGrid(y)/Lref),real((zGrid(z)-xyzful(NDref,3))/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                     real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
        elseif((isMoveDimX==1 .and. isMoveOutputX==1) .and. (isMoveDimY==1 .and. isMoveOutputY==1) .and. (isMoveDimZ/=1.and. isMoveOutputZ/=1))then
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput
            write(idfile)real((xGrid(x)-xyzful(NDref,1))/Lref),real((yGrid(y)-xyzful(NDref,2))/Lref),real(zGrid(z)/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                     real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
        elseif((isMoveDimX==1 .and. isMoveOutputX==1) .and. (isMoveDimY/=1 .and. isMoveOutputY/=1) .and. (isMoveDimZ==1.and. isMoveOutputZ==1))then
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput
            write(idfile)real((xGrid(x)-xyzful(NDref,1))/Lref),real(yGrid(y)/Lref),real((zGrid(z)-xyzful(NDref,3))/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                     real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
        elseif((isMoveDimX/=1 .and. isMoveOutputX/=1) .and. (isMoveDimY==1 .and. isMoveOutputY==1) .and. (isMoveDimZ==1.and. isMoveOutputZ==1))then
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput
            write(idfile)real(xGrid(x)/Lref),real((yGrid(y)-xyzful(NDref,2))/Lref),real((zGrid(z)-xyzful(NDref,3))/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                     real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
        else
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput      
            write(idfile)real(xGrid(x)/Lref),real(yGrid(y)/Lref),real(zGrid(z)/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                        real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
        endif

    else
            do    x=1+numOutput, xDim-numOutput
            do    y=1+numOutput, yDim-numOutput
            do    z=1+numOutput, zDim-numOutput
            write(idfile)real(xGrid(x)/Lref),real(yGrid(y)/Lref),real(zGrid(z)/Lref),real(prs(z,y,x)/(0.5*denIn*Uref**2)),    &
                        real(uuu(z,y,x,1:3)/Uref)
            enddo
            enddo
            enddo
    endif

    write(idfile) ZONEMARKER
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    if(isMoveGrid==1 .and. (isMoveOutputX==1 .or. isMoveOutputY==1 .or. isMoveOutputZ==1)) then
        
            if    ((isMoveDimX==1 .and. isMoveOutputX==1) .and. (isMoveDimY/=1 .and. isMoveOutputY/=1) .and. (isMoveDimZ/=1 .and. isMoveOutputZ/=1))then
                do    i=1,nND
                write(idfile)   real((xyzful(i,1)-xyzful(NDref,1))/Lref),real(xyzful(i,2)/Lref),real(xyzful(i,3)/Lref),real(0.0), &
                            real(velful(i,1:3)/Uref)
                enddo
            elseif((isMoveDimX/=1 .and. isMoveOutputX/=1) .and. (isMoveDimY==1 .and. isMoveOutputY==1) .and. (isMoveDimZ/=1 .and. isMoveOutputZ/=1))then
                do    i=1,nND
                write(idfile)   real(xyzful(i,1)/Lref),real((xyzful(i,2)-xyzful(NDref,2))/Lref),real(xyzful(i,3)/Lref),real(0.0), &
                            real(velful(i,1:3)/Uref)
                enddo
            elseif((isMoveDimX/=1 .and. isMoveOutputX/=1) .and. (isMoveDimY/=1 .and. isMoveOutputY/=1) .and. (isMoveDimZ==1 .and. isMoveOutputZ==1))then
                do    i=1,nND
                write(idfile)   real(xyzful(i,1)/Lref),real(xyzful(i,2)/Lref),real((xyzful(i,3)-xyzful(NDref,3))/Lref),real(0.0), &
                            real(velful(i,1:3)/Uref)
                enddo
            elseif((isMoveDimX==1 .and. isMoveOutputX==1) .and. (isMoveDimY==1 .and. isMoveOutputY==1) .and. (isMoveDimZ/=1 .and. isMoveOutputZ/=1))then
                do    i=1,nND
                write(idfile)   real((xyzful(i,1)-xyzful(NDref,1))/Lref),real((xyzful(i,2)-xyzful(NDref,2))/Lref),real(xyzful(i,3)/Lref),real(0.0), &
                            real(velful(i,1:3)/Uref)
                enddo
            elseif((isMoveDimX==1 .and. isMoveOutputX==1) .and. (isMoveDimY/=1 .and. isMoveOutputY/=1) .and. (isMoveDimZ==1 .and. isMoveOutputZ==1))then
                do    i=1,nND
                write(idfile)   real((xyzful(i,1)-xyzful(NDref,1))/Lref),real(xyzful(i,2)/Lref),real((xyzful(i,3)-xyzful(NDref,3))/Lref),real(0.0), &
                            real(velful(i,1:3)/Uref)
                enddo
            elseif((isMoveDimX/=1 .and. isMoveOutputX/=1) .and. (isMoveDimY==1 .and. isMoveOutputY==1) .and. (isMoveDimZ==1 .and. isMoveOutputZ==1))then
                do    i=1,nND
                write(idfile)   real(xyzful(i,1)/Lref),real((xyzful(i,2)-xyzful(NDref,2))/Lref),real((xyzful(i,3)-xyzful(NDref,3))/Lref),real(0.0), &
                            real(velful(i,1:3)/Uref)
                enddo
            else
                do    i=1,nND
                write(idfile)   real(xyzful(i,1:3)/Lref),real(0.0),real(velful(i,1:3)/Uref)                           
                enddo
            endif                         

    else
            do    i=1,nND
            write(idfile)   real(xyzful(i,1:3)/Lref),real(0.0),real(velful(i,1:3)/Uref)                           
            enddo
    endif

    do    i=1,nEL
        write(idfile) ele(i,1),ele(i,2),ele(i,3)
    enddo
    close(idfile)


    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot binary format
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE write_solid_field(xyzful,velful,accful,extful,ele,time,nND,nEL,isIB)
    implicit none
    integer:: nND,nEL,isIB
    real(8):: xyzful(nND,6),velful(nND,6),accful(nND,6),extful(1:nND,1:6) !,streI(1:nND),bendO(1:nND)
    integer:: ele(nEL,5)
    real(8):: time
!   -------------------------------------------------------

    integer:: i,j,ireduc
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName
    !==================================================================================================        
    integer::    nv
    integer,parameter:: namLen=40,idfile=100,numVar=12
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=[character(namLen)::'x','y','z','u','v','w','ax','ay','az','fx','fy','fz'] 
    !==================================================================================================

    write(fileName,'(F10.5)') time
    fileName = adjustr(fileName)
    DO  I=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    END DO

    if(isIB==0)then
        OPEN(idfile,FILE='./DatBody/Body'//trim(filename)//'.plt',form='unformatted')
    elseif(isIB==1)then
        OPEN(idfile,FILE='./DatBodyIB/Body'//trim(filename)//'.plt',form='unformatted')
    else
    endif

!   I. The header section.
!   =============================================        
    write(idfile) "#!TDV101"    
    write(idfile) 1
    call dumpstring(title,idfile)
    write(idfile) numVar
    do  nv=1,numVar
        call dumpstring(varname(nv),idfile)
    enddo
!   ---------------------------------------------
!   zone head for ZONE1
!   ---------------------------------------------
    write(idfile) ZONEMARKER              
    call dumpstring(zonename,idfile)
    write(idfile) -1,2,1,0,0,nND,nEL,0,0,0,0

!   =============================================
    write(idfile) EOHMARKER
!   =============================================
!   II. Data section
!   =============================================
    write(idfile) Zonemarker
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    do    i=1,nND
        write(idfile)   real(xyzful(i,1:3)),real(velful(i,1:3)),real(accful(i,1:3)),real(extful(i,1:3)) !,real(streI(i)),real(bendO(i))                                                  
    enddo
    do    i=1,nEL
        write(idfile) ele(i,1),ele(i,2),ele(i,3)
    enddo
    close(idfile)
!   =============================================
    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    write flow slice, tecplot binary
!    copyright@ RuNanHua  
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_flow_slice(xDim,yDim,zDim,XGrid,yGrid,zGrid,Lref,Uref,dIn,p,u,time,islic)
    implicit none
    integer:: x,y,z,i,isT
    integer:: xDim,yDim,zDim,islic
    real(8):: Lref,Uref,dIn
    real(8):: XGrid(xDim),yGrid(yDim),zGrid(zDim),p(zDim,yDim,xDim),u(zDim,yDim,xDim,3),time
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName
    !==================================================================================================        
    integer::    nv
    integer,parameter:: namLen=40,idfile=100,numVar=7 
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=['x','y','z','p','u','v','w']  

    write(fileName,'(F10.5)') time
    fileName = adjustr(fileName)
    do  i=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    enddo
    if    (islic==1)then    !yz plane
        open(idfile,file='./DatOthe/SlcX'//trim(fileName)//'.plt') !,form='BINARY')
        write(idfile,*)'VARIABLES = "x" "y" "z" "p" "u" "v" "w"'
        write(idfile,*)'ZONE I=',zDim,', J=',yDim,', K=',1,',F=POINT'
        x=(xDim+1)/2
        do  y=1,yDim
        do  z=1,zDim
            write(idfile,'(7D20.10)')xGrid(x)/Lref,yGrid(y)/Lref,zGrid(z)/Lref,p(z,y,x)/(0.5*dIn*Uref**2),u(z,y,x,1:3)/Uref
        enddo
        enddo
        close(idfile)
    elseif(islic==2)then    !zx plane
        open(idfile,file='./DatOthe/SlcY'//trim(fileName)//'.plt') !,form='BINARY')
        write(idfile,*)'VARIABLES = "x" "y" "z" "p" "u" "v" "w"'
        write(idfile,*)'ZONE I=',zDim,', J=',1,', K=',xDim,',F=POINT'
        y=(yDim+1)/2
        do  x=1,xDim
        do  z=1,zDim
            write(idfile,'(7D20.10)')xGrid(x)/Lref,yGrid(y)/Lref,zGrid(z)/Lref,p(z,y,x)/(0.5*dIn*Uref**2),u(z,y,x,1:3)/Uref
        enddo
        enddo
        close(idfile)   
    elseif(islic==3)then    !zx plane
        open(idfile,file='./DatOthe/SlcZ'//trim(fileName)//'.plt') !,form='BINARY')
        write(idfile,*)'VARIABLES = "x" "y" "z" "p" "u" "v" "w"'
        write(idfile,*)'ZONE I=',1,', J=',yDim,', K=',xDim,',F=POINT'
        z=(zDim+1)/2
        do  x=1,xDim
        do  y=1,yDim
            write(idfile,'(7D20.10)')xGrid(x)/Lref,yGrid(y)/Lref,zGrid(z)/Lref,p(z,y,x)/(0.5*dIn*Uref**2),u(z,y,x,1:3)/Uref
        enddo
        enddo
        close(idfile)
    else
    endif
    

    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_image()
    USE simParam
    implicit none
    integer:: x,y,z
    !==================================================================================================
    integer:: iimage        
    integer:: nv
    integer,parameter:: namLen=40,idfile=100,numVar=4
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=[character(namLen)::'x','y','z','image'] 
    !==================================================================================================
!   =============================================================================================
    open(idfile,file='./DatOthe/Image.plt',form='unformatted',access='stream')        
!   =============================================================================================
    write(idfile) "#!TDV101"    
    write(idfile) 1
    call dumpstring(title,idfile)
    write(idfile) numVar
    do  nv=1,numVar
        call dumpstring(varname(nv),idfile)
    enddo
!   ===============================================================================================
    write(idfile) ZONEMARKER              
    call dumpstring(zonename,idfile)
    write(idfile) -1,0,1,0,0,zDim,ydim,xDim,0
!   =============================================
    write(idfile) EOHMARKER
    write(idfile) ZONEMARKER
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    do    x=1, xDim
    do    y=1, yDim
    do    z=1, zDim
        write(idfile)    real(xGrid(x)/Lref),real(yGrid(y)/Lref),real(zGrid(z)/Lref),real(image(z,y,x))
    enddo
    enddo
    enddo
    close(idfile)
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write parameters for checking
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE  write_params()
    USE simParam
    implicit none
    integer:: iMT,i
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileNameR,fileNameM,fileNameS,fileNameK
    write(fileNameR,'(F10.5)') Re
    write(fileNameM,'(F10.5)') DenR
    write(fileNameK,'(F10.5)') KB
    write(fileNameS,'(F10.5)') KS

    fileNameR = adjustr(fileNameR)
    do  i=1,nameLen
        if(fileNameR(i:i)==' ')fileNameR(i:i)='0'
    enddo
    fileNameM = adjustr(fileNameM)
    do  i=1,nameLen
        if(fileNameM(i:i)==' ')fileNameM(i:i)='0'
    enddo
    fileNameS = adjustr(fileNameS)
    do  i=1,nameLen
        if(fileNameS(i:i)==' ')fileNameS(i:i)='0'
    enddo
    fileNameK = adjustr(fileNameK)
    do  i=1,nameLen
        if(fileNameK(i:i)==' ')fileNameK(i:i)='0'
    enddo
    open(111,file='./R'//trim(fileNameR)//'M'//trim(fileNameM)//'S'//trim(fileNameS)//'K'//trim(fileNameK))
    close(111)

        open(111,file='./Check.dat')
        write(111,'(A      )')'===================================='
        if    (iBodyModel==1)then
            write(111,'(A      )')'This is a RIGID    body problem'
        elseif(iBodyModel==2)then
            write(111,'(A      )')'This is a FLRXIBLE body problem'
        endif
        if    (isMotionGiven(1)==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in X-direction freely'
        endif
        if    (isMotionGiven(2)==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in Y-direction freely'
        endif
        if    (isMotionGiven(3)==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in Z-direction freely'
        endif
        write(111,'(A      )')'===================================='
        write(111,'(3(A,1x,I8,2x))')'nND=',nND,'nEL=',nEL,'nEQ=',nEQ
        write(111,'(3(A,1x,I8,2x))')'nMT=',nMT,'nBD=',nBD,'nSTF=',nSTF
        write(111,'(A      )')'===================================='
        write(111,'(A,3I20.10)') 'xDim*yDim*zDim          :',xDim*yDim*zDim
        write(111,'(A,3I20.10)') 'xDim,yDim,zDim          :',xDim, yDim, zDim
        write(111,'(A,3F20.10)') 'dh,dt,ratio             :',dh,     dt, ratio
        write(111,'(A,3F20.10)') 'dxmin,dymin,dzmin       :',dxmin, dymin, dzmin
        write(111,'(A,3F20.10)') 'dxmax,dymax,dzmax       :',dxmax, dymax, dzmax
        write(111,'(A,3F20.10)') 'cptxMin,cptyMin,cptzMin :',cptxMin, cptyMin, cptzMin
        write(111,'(A,3F20.10)') 'cptxMax,cptyMax,cptzMax :',cptxMax, cptyMax, cptzMax
        write(111,'(A,2F20.10)') 'elmin,elmax             :',elmin, elmax
        write(111,'(A      )')'===================================='
        write(111,'(A,F20.10)')'Re   =',Re  
        write(111,'(A,F20.10)')'KB   =',KB
        write(111,'(A,F20.10)')'KS   =',KS
        write(111,'(A,F20.10)')'EmR  =',EmR
        write(111,'(A,F20.10)')'tcR  =',tcR        
        write(111,'(A,F20.10)')'denR =',denR
        write(111,'(A,F20.10)')'Ampl =',maxval(dabs(XYZAmpl(1:3)))
        write(111,'(A,F20.10)')'AR   =',AR
        write(111,'(A      )')'===================================='
        write(111,'(A,F20.10)')'Asfac=',Asfac
        write(111,'(A,F20.10)')'Lchod=',Lchod
        write(111,'(A,F20.10)')'Lspan=',Lspan       
        write(111,'(A,F20.10)')'psR  =',psR
        write(111,'(A,F20.10)')'uMax =',uMax
        write(111,'(A,F20.10)')'Lref =',Lref
        write(111,'(A,F20.10)')'Uref =',Uref
        write(111,'(A,F20.10)')'Tref =',Tref
        write(111,'(A,F20.10)')'Pref =',Pref
        write(111,'(A,F20.10)')'Eref =',Eref
        write(111,'(A,F20.10)')'Fref =',Fref
        write(111,'(A,F20.10)')'Aref =',Aref
        write(111,'(A,F20.10)')'mxMa =',uMax/dsqrt(Cs2)
        write(111,'(A,F20.10)')'St   =',St
        write(111,'(A,F20.10)')'Nu   =',Nu
        write(111,'(A,F20.10)')'Mu   =',Mu
        write(111,'(A,F20.10)')'Tau  =',Tau
        write(111,'(A,F20.10)')'Omega=',Omega
        write(111,'(A      )')'===================================='
        do iMT=1,nMT
        write(111,'(A      )')'================='
        write(111,'(A,I5.5 )')'MT:',iMT
        write(111,'(A,D20.10 )')'E    =',prop(1,1)
        write(111,'(A,D20.10 )')'G    =',prop(1,2)
        write(111,'(A,D20.10 )')'h    =',prop(1,3)
        write(111,'(A,D20.10 )')'rho  =',prop(iMT,4)
        write(111,'(A,D20.10 )')'gamma=',prop(iMT,5)
        write(111,'(A,D20.10 )')'Ip   =',prop(iMT,6)
        write(111,'(A,D20.10 )')'alpha=',prop(iMT,7)
        write(111,'(A,D20.10 )')'beta =',prop(iMT,8)
        enddo
        write(111,'(A      )')'====================================' 
        write(111,'(A,3I20.10)')'T1,T2,T3:',NDtl(1:3) 
        write(111,'(A,3I20.10)')'H1,H2,H3:',NDhd(1:3)
        write(111,'(A,1I20.10)')'      CT:',NDct 
        write(111,'(A      )')'===================================='     
        close(111)
    END SUBROUTINE



!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    compute strain energy
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE strain_energy_D(  strainEnergy, xord0,yord0,zord0,xord,yord,zord, &
                                 ele,prop,triad_n1,triad_n2,triad_n3,triad_ee,triad_e0,triad_nn, &
                                 nND,nEL,nMT &
                               )
    implicit none
    integer:: nMT,nEL,nND
    integer:: ele(nEL,5)
    real(8):: xord0(nND), yord0(nND), zord0(nND), strainEnergy(nEL,2)
    real(8):: xord(nND), yord(nND), zord(nND)
!
    real(8):: prop(nMT,10)
 
    real(8):: ek9(9,9),ekb12(12,12),ekb12Strech(12,12),ekb12BendTor(12,12)
    real(8):: ekb(18,18),ekbInplane(18,18),ekbOutplane(18,18)
!
    real(8):: triad_nn(3,3,nND),triad_ee(3,3,nEL),triad_e0(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL),triad_n3(3,3,nEL)

    real(8):: triad_00(3,3),triad_11(3,3),triad_22(3,3),rr(3,3)
    real(8):: ub(18),dl,temp(6)
!
    real(8):: dx0,dy0,dz0,du,dv,dw,dx,dy,dz,xl0
    real(8):: tx1,tx2,tx3,ty1,ty2,ty3,tz1,tz2,tz3,tx,ty,tz
    real(8):: xyz012(3),xyz013(3),xyzb012(3),xyzb013(3)
    real(8):: xyz12(3),xyz13(3),xyzb12(3),xyzb13(3),xyz11(3),xyzb11(3)
    real(8):: xb01,xb02,xb03,yb01,yb02,yb03,xb1,xb2,xb3,yb1,yb2,yb3
    real(8):: uub(6)

    real(8):: e0,g0,a0,b0,r0,zix0,ziy0,ziz0,xl,fxx,t0,pl0,zip0,zia0,zib0,alpha,beta
    integer:: i,j,k,n,i1,j1,k1,mat,nELt

!   For each element, calculate the nodal forces
    do    n=1,nEL

        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)        
        nELt= ele(n,4)
        mat = ele(n,5)
        if    ( nELt == 2) then
!           frame
            e0  =prop(mat,1)
            g0  =prop(mat,2)
            a0  =prop(mat,3)
            r0  =prop(mat,4)
            b0  =prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)

            dx0 = xord0(j1) - xord0(i1)
            dy0 = yord0(j1) - yord0(i1)
            dz0 = zord0(j1) - zord0(i1)
            xl0 = dsqrt(dx0*dx0+dy0*dy0+dz0*dz0)
!
!           orientation
            du = (xord(j1)-xord0(j1))-(xord(i1)-xord0(i1))
            dv = (yord(j1)-yord0(j1))-(yord(i1)-yord0(i1))
            dw = (zord(j1)-zord0(j1))-(zord(i1)-zord0(i1))

            dx = dx0 + du
            dy = dy0 + dv
            dz = dz0 + dw
            xl =dsqrt(dx*dx+dy*dy+dz*dz)
!
            dl = ( (2*dx0+du)*du +(2*dy0+dv)*dv +(2*dz0+dw)*dw )/ (xl+xl0)
!           get twisting angles
            do    i=1,3
            do    j=1,3
                triad_00(i,j)=triad_ee(i,j,n)
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n1(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)
!
            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

!            non-zero ty1 tz1 u2 tx2 ty2 tz2
            ub(1)=0.0
            ub(2)=0.0
            ub(3)=0.0
            ub(4)=tx1
            ub(5)=ty1
            ub(6)=tz1
!
            ub(7)=dl
            ub(8)=0.0
            ub(9)=0.0
            ub(10)=tx2
            ub(11)=ty2
            ub(12)=tz2

!
!
!           get current stiffness in local coords. use L0

            call elmstfFRM_D(xl0,zix0,ziy0,ziz0,a0,e0,g0,ekb12,nELt )

            ekb12Strech(1:12,1:12)=0.0d0
            ekb12Strech(1,1)=ekb12(1,1)
            ekb12Strech(1,7)=ekb12(1,7)
            ekb12Strech(7,7)=ekb12(7,7)
            ekb12Strech(7,1)=ekb12(7,1)

            ekb12BendTor(1:12,1:12)=ekb12(1:12,1:12)
            ekb12BendTor(1,1)=0.0d0
            ekb12BendTor(1,7)=0.0d0
            ekb12BendTor(7,7)=0.0d0
            ekb12BendTor(7,1)=0.0d0

!           nodal forces in local coords
            strainEnergy(n,1)=0.5*sum(matmul(ekb12Strech(1:12,1:12),ub(1:12))*ub(1:12))
            strainEnergy(n,2)=0.5*sum(matmul(ekb12BendTor(1:12,1:12),ub(1:12))*ub(1:12))

        elseif (nELt == 3) then
!           plate
            e0=prop(mat,1)
            g0=prop(mat,2)
            t0=prop(mat,3)
            r0=prop(mat,4)
            pl0=prop(mat,5)
            zip0=prop(mat,6)
            zia0=prop(mat,7)
            zib0=prop(mat,8)
!
!           original lengths of triangle
            xyz012(1) = xord0(j1)-xord0(i1)
            xyz012(2) = yord0(j1)-yord0(i1)
            xyz012(3) = zord0(j1)-zord0(i1)
            xyz013(1) = xord0(k1)-xord0(i1)
            xyz013(2) = yord0(k1)-yord0(i1)
            xyz013(3) = zord0(k1)-zord0(i1)
!
!           use element triad to rotate to local x=[Rt]x
            do    i=1,3
            do    j=1,3
                rr(i,j)=triad_e0(i,j,n)
            enddo
            enddo
            do    i=1,3
                xyzb012(i)=0.0
                xyzb013(i)=0.0
                do    k=1,3
                    xyzb012(i) = xyzb012(i) + rr(k,i)*xyz012(k)
                    xyzb013(i) = xyzb013(i) + rr(k,i)*xyz013(k)
                enddo
            enddo
            xb01=0.0
            yb01=0.0
            xb02=xyzb012(1)
            yb02=xyzb012(2)
            xb03=xyzb013(1)
            yb03=xyzb013(2)

!           current lengths of triangle
            xyz11(1) = xord(i1)-xord0(i1)
            xyz11(2) = yord(i1)-yord0(i1)
            xyz11(3) = zord(i1)-zord0(i1)
            xyz12(1) = xord(j1)-xord0(i1)
            xyz12(2) = yord(j1)-yord0(i1)
            xyz12(3) = zord(j1)-zord0(i1)
            xyz13(1) = xord(k1)-xord0(i1)
            xyz13(2) = yord(k1)-yord0(i1)
            xyz13(3) = zord(k1)-zord0(i1)
!
!           use element triad to rotate to local x=[Rt]x
            do    i=1,3
            do    j=1,3
                rr(i,j)=triad_ee(i,j,n)
            enddo
            enddo
            do    i=1,3
                xyzb11(i)=0.0
                xyzb12(i)=0.0
                xyzb13(i)=0.0
                do    k=1,3
                    xyzb11(i) = xyzb11(i) + rr(k,i)*xyz11(k)
                    xyzb12(i) = xyzb12(i) + rr(k,i)*xyz12(k)
                    xyzb13(i) = xyzb13(i) + rr(k,i)*xyz13(k)
                enddo
            enddo

            xb1=xyzb11(1)
            yb1=xyzb11(2)
            xb2=xyzb12(1)
            yb2=xyzb12(2)
            xb3=xyzb13(1)
            yb3=xyzb13(2)

!           get angles ref to orig from triads
            do    i=1,3
            do    j=1,3
!2001!!!        triad_00(i,j)=triad_e0(i,j,n)
                triad_00(i,j)=triad_ee(i,j,n)
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n1(i,j,n)
            enddo
            enddo

            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)

            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo

            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n3(i,j,n)
            enddo
            enddo

            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx3,ty3,tz3)
!
!           local DoF
            ub(1)=xb1-xb01
            ub(2)=yb1-yb01
            ub(3)=0.0
            ub(4)=tx1
            ub(5)=ty1
            ub(6)=tz1

            ub(7)=xb2-xb02
            ub(8)=yb2-yb02
            ub(9)=0.0
            ub(10)=tx2
            ub(11)=ty2
            ub(12)=tz2

            ub(13)=xb3-xb03
            ub(14)=yb3-yb03
            ub(15)=0.0
            ub(16)=tx3
            ub(17)=ty3
            ub(18)=tz3

!           Calculate LOCAL stiffness, use orig coords
!            membrane
            ekbInplane(1:18,1:18)=0.0
            alpha=zia0
            beta =zib0
            call elmstfMRT_D(e0,g0,t0,pl0,alpha,beta,xb01,xb02,xb03,yb01,yb02,yb03,ekbInplane,ek9)
!           flexure
            ekbOutplane(1:18,1:18)=0.0
            call elmstfDKT_D(e0,g0,t0,zip0,          xb01,xb02,xb03,yb01,yb02,yb03,ekbOutplane,ek9)

            ekb=ekbInplane+ekbOutplane

            strainEnergy(n,1)=0.5*sum(matmul(ekbInplane(1:18,1:18),ub(1:18))*ub(1:18))
            strainEnergy(n,2)=0.5*sum(matmul(ekbOutplane(1:18,1:18),ub(1:18))*ub(1:18))
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif

    enddo


    
    return
    ENDSUBROUTINE strain_energy_D