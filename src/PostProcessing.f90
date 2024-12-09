!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write parameters for checking
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE  write_params()
        USE SolidBody
        implicit none
        integer:: i
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileNameR
        write(fileNameR,'(F10.5)') Re

        fileNameR = adjustr(fileNameR)
        do  i=1,nameLen
            if(fileNameR(i:i)==' ')fileNameR(i:i)='0'
        enddo

        open(111,file='./Check.dat')
        write(111,'(A      )')'===================================='
        if    (maxval(iBodyModel(:)).eq.1)then
            write(111,'(A      )')'This is a RIGID    body problem'
        elseif(maxval(iBodyModel(:)).eq.2)then
            write(111,'(A      )')'This is a FLRXIBLE body problem'
        else
            write(111,'(A      )')'This is a FLRXIBLE And RIGID body problem'
        endif
        if    (minval(VBodies(:)%rbm%isMotionGiven(1))==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in X-direction freely'
        endif
        if    (minval(VBodies(:)%rbm%isMotionGiven(2))==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in Y-direction freely'
        endif
        if    (minval(VBodies(:)%rbm%isMotionGiven(3))==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in Z-direction freely'
        endif
        write(111,'(A      )')'===================================================================='
        write(111,'(A,I20.10)')'number of fish is',nFish
        write(111,'(A      )')'===================================='
        write(111,'(A,3I20.10)') 'xDim,yDim,zDim          :',xDim, yDim, zDim
        write(111,'(A,3F20.10)') 'dh,dt,ratio             :',dh, dt, ratio
        write(111,'(A,3F20.10)') 'dxmin,dymin,dzmin       :',dxmin, dymin, dzmin
        write(111,'(A,3F20.10)') 'dxmax,dymax,dzmax       :',dxmax, dymax, dzmax
        write(111,'(A,3F20.10)') 'cptxMin,cptyMin,cptzMin :',cptxMin, cptyMin, cptzMin
        write(111,'(A,3F20.10)') 'cptxMax,cptyMax,cptzMax :',cptxMax, cptyMax, cptzMax
        write(111,'(A,2F20.10)') 'elmin,elmax             :',minval(VBodies(:)%rbm%elmin), maxval(VBodies(:)%rbm%elmax)
        write(111,'(A      )')'===================================='
        write(111,'(A,F20.10)')'Re   =',Re
        write(111,'(A,F20.10)')'Lref =',Lref
        write(111,'(A,F20.10)')'Uref =',Uref
        write(111,'(A,F20.10)')'Tref =',Tref
        write(111,'(A      )')'===================================='
        write(111,'(A,F20.10)')'Freq =',maxval(VBodies(:)%rbm%Freq)
        write(111,'(A,F20.10)')'Ampl =',maxval([dabs(VBodies(:)%rbm%XYZAmpl(1)),dabs(VBodies(:)%rbm%XYZAmpl(2)),dabs(VBodies(:)%rbm%XYZAmpl(3))])
        write(111,'(A,F20.10)')'Lchod=',Lchod
        write(111,'(A,F20.10)')'Lspan=',Lspan
        write(111,'(A,F20.10)')'Asfac=',Asfac
        write(111,'(A,F20.10)')'AR   =',AR
        write(111,'(A      )')'===================================='
        write(111,'(A,F20.10)')'denIn=',denIn
        write(111,'(A,F20.10)')'Aref =',Aref
        write(111,'(A,F20.10)')'Pref =',Pref
        write(111,'(A,F20.10)')'Eref =',Eref
        write(111,'(A,F20.10)')'Fref =',Fref
        write(111,'(A,F20.10)')'uMax =',uMax
        write(111,'(A,F20.10)')'maxMa=',uMax/dsqrt(Cs2)
        write(111,'(A,F20.10)')'Tau  =',Tau
        write(111,'(A,F20.10)')'Omega=',Omega
        write(111,'(A,F20.10)')'Nu   =',Nu
        write(111,'(A,F20.10)')'Mu   =',Mu

        call Write_solid_Check(111)
        close(111)
    END SUBROUTINE
