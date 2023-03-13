!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE writeFlow(isT)
    USE simParam
    implicit none
    integer:: x,y,i,isT,iFish
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName
    !==================================================================================================        
    integer::    nv
    integer,parameter:: namLen=40,idfile=100,numVar=6 
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=[character(namLen)::'x','y','p','u','v','vort']  

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

    DO  x=2, xDim-1
    DO  y=2, yDim-1
        vor(y,x)=(uuu(y,x+1,2)-uuu(y,x-1,2))/(2.0*dx(x))-(uuu(y+1,x,1)-uuu(y-1,x,1))/(2.0*dy(y))
    ENDDO
    ENDDO
    DO    y=2, yDim-1
        vor(y,1)=2*vor(y,2)-VOR(y,3)
        vor(y,xDim)=2*vor(y,xDim-1)-vor(y,xDim-2)
    ENDDO
    DO    x=1,xDim
        vor(1,x)=2*vor(2,x)-vor(3,x)
        vor(yDim,x)=2*vor(yDim-1,x)-vor(yDim-2,x)
    ENDDO

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
    write(idfile) -1,0,1,0,0,ydim-4,xDim-4,1,0
    
    do iFish=1,nFish
       write(idfile) ZONEMARKER              
       call dumpstring(zonename,idfile)
       write(idfile) -1,1,1,0,0,nND(iFish),nEL(iFish),0,0,0,0
    enddo
 
    write(idfile) EOHMARKER
    write(idfile) ZONEMARKER
!    variable data format, 1=Float, 2=Double, 3=LongInt, 4=ShortInt, 5=Byte, 6=Bit
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    if(isMoveGrid==1 .and. isMoveOutput==1)then
        do  x=3, xDim-2
        do  y=3, yDim-2
            if    (iMoveDim==1)then
            write(idfile) real((xGrid(x)-xyzful(1,NDref,1))/Lref),real(yGrid(y)/Lref),real(prs(y,x)/(0.5*denIn*Uref**2)), &
                          real(uuu(y,x,1:2)/Uref),real(vor(y,x)/(Uref/Lref))
            elseif(iMoveDim==2)then
            write(idfile) real(xGrid(x)/Lref),real((yGrid(y)-xyzful(1,NDref,2))/Lref),real(prs(y,x)/(0.5*denIn*Uref**2)), &
                          real(uuu(y,x,1:2)/Uref),real(vor(y,x)/(Uref/Lref))
            else
            endif
        enddo
        enddo
    else
        do  x=3, xDim-2
        do  y=3, yDim-2
            write(idfile) real(xGrid(x)/Lref),real(yGrid(y)/Lref),real(prs(y,x)/(0.5*denIn*Uref**2)), &
                        real(uuu(y,x,1:2)/Uref),real(vor(y,x)/(Uref/Lref))
        enddo
        enddo
    endif
    
    do iFish=1,nFish        
    write(idfile) ZONEMARKER
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    if(isMoveGrid==1 .and. isMoveOutput==1) then
        do    i=1,nND(iFish)
        if    (iMoveDim==1)then
        write(idfile)real((xyzful(iFish,i,1)-xyzful(iFish,NDref,1))/Lref), real(xyzful(iFish,i,2)/Lref),real(0.0d0),real(velful(iFish,i,1:2)/Uref),real(0.0d0)
        elseif(iMoveDim==2)then
        write(idfile)real(xyzful(iFish,i,1)/Lref), real((xyzful(iFish,i,2)-xyzful(iFish,NDref,2))/Lref),real(0.0d0),real(velful(iFish,i,1:2)/Uref),real(0.0d0)
        else
        endif                        
        enddo
    else
        do    i=1,nND(iFish)
        write(idfile)   real(xyzful(iFish,i,1:2)/Lref),real(0.0d0),real(velful(iFish,i,1:2)/Uref),real(0.0d0)                        
        enddo
    endif
    do    i=1,nEL(iFish)
        write(idfile) ele(iFish,i,1),ele(iFish,i,2) !,ele(iFish,i,3)
    enddo
    
    enddo
    
    close(idfile)

    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!    Runan Hua, Oct. 2010
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE wrtBody(xyzful,velful,accful,extful,ele,time,nND,nEL,isIB)        
    implicit none
    integer:: nND,nEL,isIB
    real(8):: xyzful(nND,6),velful(nND,6),accful(1:nND,1:6),extful(1:nND,1:6)
    integer:: ele(nEL,5)
    real(8):: time
!   -------------------------------------------------------

    integer:: i,j,ireduc,iFish
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName
    !==================================================================================================        
    integer::    nv
    integer,parameter:: namLen=40,idfile=100,numVar=8
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=[character(namLen)::'x','y','u','v','ax','ay','fx','fy'] 
    !==================================================================================================

    write(fileName,'(F10.5)') time
    fileName = adjustr(fileName)
    DO  I=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    END DO

    if(isIB==0)then
        OPEN(idfile,FILE='./DatBody/Body'//trim(filename)//'.plt',form='unformatted',access='stream')
    elseif(isIB==1)then
        OPEN(idfile,FILE='./DatBodyIB/Body'//trim(filename)//'.plt',form='unformatted',access='stream')
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
    write(idfile) -1,1,1,0,0,nND,nEL,0,0,0,0

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
        write(idfile)   real(xyzful(i,1:2)),real(velful(i,1:2)),real(accful(i,1:2)),real(extful(i,1:2))                                                   
    enddo
    do    i=1,nEL
        write(idfile) ele(i,1),ele(i,2) !,ele(i,3)
    enddo
    close(idfile)
!   =============================================
    END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    Runan Hua, Oct. 2010
!   ZeRui Peng, Mar. 2015
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    SUBROUTINE wrtBody2(nFish,xyzful,velful,accful,extful,repful,ele,time,nND,nEL,nND_max,nEL_max,isIB)
    implicit none
    integer:: nFish,isIB,nND_max,nEL_max
    integer:: nND(nFish),nEL(nFish)
    
    real(8):: xyzful(nFish,nND_max,6),velful(nFish,nND_max,6),accful(nFish,1:nND_max,1:6)
    real(8):: extful(nFish,1:nND_max,1:6),repful(nFish,1:nND_max,1:6)
    integer:: ele(nFish,nEL_max,5)
    real(8):: time

    !local
    integer:: i,j,ireduc,iFish
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName
    !==================================================================================================        
    integer::    nv
    integer,parameter:: namLen=50,idfile=100,numVar=10
    real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
    character(namLen):: ZoneName='ZONE 1',title="Binary File.",    &
                        varname(numVar)=[character(namLen)::'x','y','u','v','ax','ay','fxi','fyi','fxr','fyr'] 
    !==================================================================================================

    write(fileName,'(F9.5)') time
    fileName = adjustr(fileName)
    DO  I=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    END DO

    if(isIB==0)then
        OPEN(idfile,FILE='./DatBody/Body'//trim(filename)//'.plt',form='unformatted',access='stream')
    elseif(isIB==1)then
        OPEN(idfile,FILE='./DatBodyIB/Body'//trim(filename)//'.plt',form='unformatted',access='stream')
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
    do iFish=1,nFish
       write(idfile) ZONEMARKER              
       call dumpstring(zonename,idfile)
       write(idfile) -1,1,1,0,0,nND(iFish),nEL(iFish),0,0,0,0
    enddo
!   =============================================
    write(idfile) EOHMARKER
!   =============================================
!   II. Data section
!   =============================================
    do iFish=1,nFish

    write(idfile) Zonemarker
    do  nv=1,numVar
        write(idfile) 1                                 
    enddo
    write(idfile) 0,-1
    do    i=1,nND(iFish)
        write(idfile)   real(xyzful(iFish,i,1:2)),real(velful(iFish,i,1:2)),real(accful(iFish,i,1:2)),real(extful(iFish,i,1:2)),  &
                        real(repful(iFish,i,1:2))                                                 
    enddo
    do    i=1,nEL(iFish)
        write(idfile) ele(iFish,i,1),ele(iFish,i,2) !,ele(iFish,i,3)
    enddo

    enddo !iFish=1,nFish
    close(idfile)
!   =============================================
    END SUBROUTINE

    !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    Runan Hua, Oct. 2010
!    Ankang Gao, Jul. 2022, write reduced ASCII file
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    SUBROUTINE writeBody3(nFish,xyzful,velful,accful,extful,repful,ele,time,nND,nEL,nND_max,nEL_max,isIB)
        implicit none
        integer:: nFish,isIB,nND_max,nEL_max
        integer:: nND(nFish),nEL(nFish)
        
        real(8):: xyzful(nFish,nND_max,6),velful(nFish,nND_max,6),accful(nFish,1:nND_max,1:6)
        real(8):: extful(nFish,1:nND_max,1:6),repful(nFish,1:nND_max,1:6)
        integer:: ele(nFish,nEL_max,5)
        real(8):: time
    
        !local
        integer:: i,j,ireduc,iFish
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName, idstr
        !==================================================================================================        
        integer::    nv
        integer,parameter:: namLen=50,idfile=100,numVar=10
        real(4),parameter:: ZONEMARKER=299.0,EOHMARKER =357.0
        character(namLen):: varname(numVar)=[character(namLen)::'x','y','u','v','ax','ay','fxi','fyi','fxr','fyr'] 
        !==================================================================================================
    
        write(fileName,'(F9.5)') time
        fileName = adjustr(fileName)
        DO  I=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        END DO

        do iFish = 1, nFish
            write(idstr, '(I3.3)') iFish ! assume iFish < 1000
            if(isIB==0)then
                OPEN(idfile,FILE='./DatBody/Body'//trim(idstr)//'_'//trim(filename)//'.dat')
            elseif(isIB==1)then
                OPEN(idfile,FILE='./DatBodyIB/Body'//trim(idstr)//'_'//trim(filename)//'.dat')
            else
            endif
            write(idfile, '(A)', advance='no') 'variables = '
            do i=1,numVar-1
                write(idfile, '(2A)', advance='no') trim(varname(i)), ','
            enddo
            write(idfile, '(A)') varname(numVar)
            do i=1,nND(iFish)
                write(idfile, '(10E28.18 )')   real(xyzful(iFish,i,1:2)),real(velful(iFish,i,1:2)),real(accful(iFish,i,1:2)),real(extful(iFish,i,1:2)),  &
                            real(repful(iFish,i,1:2))  
            enddo
            close(idfile)
        enddo
    !   =============================================
        END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE  wrtPara()  !check.dat
    USE simParam
    implicit none
    integer:: iMT,i,iFish
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileNameR,fileNameM,fileNameS,fileNameK
    write(fileNameR,'(F10.5)') Re 
    write(fileNameM,'(F10.5)') DenR(1)
    write(fileNameK,'(F10.5)') KB(1)
    write(fileNameS,'(F10.5)') KS(1)

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
        write(111,'(A      )')'========================================================================'
        if    (maxval(iBodyModel(:))==1)then
            write(111,'(A      )')'This is a RIGID body problem'
        elseif(minval(iBodyModel(:))==2)then
            write(111,'(A      )')'This is a FLRXIBLE body problem'
        else
            write(111,'(A      )')'This is a FLRXIBLE And RIGID body problem'
        endif
        if    (minval(isMotionGiven(:,1))==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in X-direction freely'
        endif
        if    (minval(isMotionGiven(:,2))==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in Y-direction freely'
        endif

        write(111,'(A      )')'===================================================================='
        write(111,'(A,I20.10)')'number of fish is',nFish

        write(111,'(A      )')'===================================================================='
        write(111,'(A,2I20.10)') 'xDim,   yDim    :',xDim,yDim
        write(111,'(A,2F20.10)') 'cptxMin,cptyMin :',cptxMin,cptyMin
        write(111,'(A,2F20.10)') 'cptxMax,cptyMax :',cptxMax,cptyMax
        write(111,'(A,2F20.10)') 'dxmin,  dymin   :',dxmin,dymin
        write(111,'(A,2F20.10)') 'dxmax,  dymax   :',dxmax,dymax
        write(111,'(A,2F20.10)') 'dh,     dt      :',dh,dt
        write(111,'(A,1F20.10)') 'ratio           :',ratio

        write(111,'(A      )')'===================================================================='
        write(111,'(A,F20.10)')'Re   =',Re 

        write(111,'(A      )')'===================================================================='
        write(111,'(A,F20.10)')'Freq =',Freq
        write(111,'(A,F20.10)')'Ampl =',maxval(dabs(XYZAmpl(1:nFish,1:3)))
        write(111,'(A,F20.10)')'Lref =',Lref
        write(111,'(A,F20.10)')'Uref =',Uref
        write(111,'(A,F20.10)')'Tref =',Tref

        write(111,'(A      )')'===================================================================='
        !write(111,'(A,F20.10)')'Lleng=',Lleng
        write(111,'(A,F20.10)')'Pref =',Pref
        write(111,'(A,F20.10)')'Eref =',Eref
        write(111,'(A,F20.10)')'Fref =',Fref
        write(111,'(A,F20.10)')'Aref =',Aref
        write(111,'(A,F20.10)')'St   =',St
        write(111,'(A,F20.10)')'uMax =',uMax
        write(111,'(A,F20.10)')'maxMa=',uMax/dsqrt(Cs2)
        write(111,'(A,F20.10)')'Tau  =',Tau
        write(111,'(A,F20.10)')'Omega=',Omega
        write(111,'(A,F20.10)')'Nu   =',Nu
        write(111,'(A,F20.10)')'Mu   =',Mu

        do iFish=1,nFish
        write(111,'(A      )')'===================================================================='
        write(111,'(A,I20.10)')'Fish number is',iFish
        write(111,'(A,1F20.10)')'Lchod       =',Lchod(iFish)
        write(111,'(A,2F20.10)')'denR,psR    =',denR(iFish),psR(iFish)
        write(111,'(A,2F20.10)')'KB,  KS     =',KB(iFish),KS(iFish)
        write(111,'(A,2F20.10)')'EmR, tcR    =',EmR(iFish),tcR(iFish) 
        write(111,'(A,1x,3F20.10,2x)')'XYZo(1:3)   =',XYZo(iFish,1:3)
        write(111,'(A,1x,3F20.10,2x)')'XYZAmpl(1:3)=',XYZAmpl(iFish,1:3) 
        write(111,'(A,1x,3F20.10,2x)')'XYZPhi(1:3) =',XYZPhi(iFish,1:3) 
        write(111,'(A,1x,3F20.10,2x)')'AoAo(1:3)   =',AoAo(iFish,1:3)
        write(111,'(A,1x,3F20.10,2x)')'AoAAmpl(1:3)=',AoAAmpl(iFish,1:3) 
        write(111,'(A,1x,3F20.10,2x)')'AoAPhi(1:3) =',AoAPhi(iFish,1:3) 
        write(111,'(3(A,1x,I20.10,2x))')'nND=',nND(iFish),'nEL=',nEL(iFish),'nEQ=',nEQ(iFish)
        write(111,'(3(A,1x,I20.10,2x))')'nMT=',nMT(iFish),'nBD=',nBD(iFish),'maxstiff=',maxstiff(iFish)
        enddo      
        
        do iFish=1,nFish
        do iMT=1,nMT(iFish)
        write(111,'(A      )')'===================================================================='
        write(111,'(A,I5.5 )')'iFish:',iFish
        write(111,'(A,I5.5 )')'MT:',iMT
        write(111,'(A,E20.10)')'E    =',prop(iFish,iMT,1)
        write(111,'(A,E20.10)')'G    =',prop(iFish,iMT,2)
        write(111,'(A,E20.10)')'A    =',prop(iFish,iMT,3)
        write(111,'(A,E20.10)')'rho  =',prop(iFish,iMT,4)
        write(111,'(A,E20.10)')'gamma=',prop(iFish,iMT,5)
        write(111,'(A,E20.10)')'Ix   =',prop(iFish,iMT,6)
        write(111,'(A,E20.10)')'Iy   =',prop(iFish,iMT,7)
        write(111,'(A,E20.10)')'Iz   =',prop(iFish,iMT,8)
        enddo
        enddo
        write(111,'(A      )')'===================================================================='      
        close(111)
    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE writeImag()
    USE simParam
    implicit none
    integer:: x,y
    open(14,file='./DatOthe/Image.plt')
    write(14,*)'variables="x","y","image"'
    WRITE(14,*) 'ZONE T="ZONE 1"'
    WRITE(14,*) 'I=',yDim, ' J=',xDim
    WRITE(14,*) 'F=POINT'
    do  x=1, xDim
    do  y=1, yDim
        write(14,'(2f20.10,i10)') xGrid(x),yGrid(y),image(y,x)
    enddo
    enddo
    close(14)

    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   compute strain energy
!
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
!           membrane
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
