!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_solid_field(xyzful,velful,accful,extful,repful,ele,time,nND,nEL,nND_max,nEL_max,nFish)
    implicit none
    integer:: nND_max,nEL_max,nFish
    integer:: ele(nEL_max,5,nFish),nND(nFish),nEL(nFish)
    real(8):: xyzful(nND_max,6,nFish),velful(nND_max,6,nFish),accful(nND_max,6,nFish),extful(1:nND_max,1:6,nFish),repful(1:nND_max,1:6,nFish)
    real(8):: time
!   -------------------------------------------------------

    integer:: i,iFish,ElmType
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName,idstr
    !==================================================================================================
    integer,parameter:: namLen=40,idfile=100,numVar=15
    character(namLen):: varname(numVar)=[character(namLen)::'x','y','z','u','v','w','ax','ay','az','fxi','fyi','fzi','fxr','fyr','fzr']
    !==================================================================================================

    write(fileName,'(I10)') nint(time*1d5)
    fileName = adjustr(fileName)
    DO  I=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    END DO

    do iFish=1,nFish

        ElmType = ele(1,4,iFish)

        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        OPEN(idfile,FILE='./DatBody/Body'//trim(idstr)//'_'//trim(fileName)//'.dat')
        
        ! Write header information
        write(idfile, '(A)') 'TITLE    = "ASCII File."'
        write(idfile, '(A)', advance='no') 'variables= '
        do i=1,numVar-1
            write(idfile, '(3A)', advance='no') '"', trim(varname(i)), '" '
        enddo
        write(idfile, '(A)') varname(numVar)

        write(idfile, '(A)') 'ZONE    T= "ZONE 1"'
        write(idfile, '(A)') ' STRANDID=0, SOLUTIONTIME=0'
        write(idfile, '(A,I8,A,I8,A)', advance='no') ' Nodes=',nND(iFish),', Elements=',nEL(iFish),', ZONETYPE='
        if(ElmType.eq.2) then
            write(idfile, '(A)') 'FELINESEG'
        elseif (ElmType.eq.3) then
            write(idfile, '(A)') 'FETRIANGLE'
        endif
        write(idfile, '(A)') ' DATAPACKING=POINT'
        write(idfile, '(A)', advance='no') ' DT=('
        do i=1,numVar-1
            write(idfile, '(A)', advance='no') 'SINGLE '
        enddo
        write(idfile, '(A)') 'SINGLE )'

        ! Write node data
        do i=1,nND(iFish)
            write(idfile, '(10E28.18 )')   real(xyzful(i,1:3,iFish)),real(velful(i,1:3,iFish)),real(accful(i,1:3,iFish)),real(extful(i,1:3,iFish)),real(repful(i,1:3,iFish))
        enddo

        ! Write element data
        if(ElmType.eq.2) then
            do i = 1, nEL(iFish)
                write(idfile, *) ele(i,1,iFish),ele(i,2,iFish)
            enddo
        elseif (ElmType.eq.3) then
            do i = 1, nEL(iFish)
                write(idfile, *) ele(i,1,iFish),ele(i,2,iFish),ele(i,3,iFish)
            enddo
        endif

        close(idfile)
    enddo
!   =============================================
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot binary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_solidIB_field(xyzfulIB,ele,time,nND,nEL,nND_max,nEL_max,Nspan,nFish)
    implicit none
    integer,intent(in):: nND_max,nEL_max,nFish,Nspan(nFish)
    integer,intent(in):: ele(nEL_max,5,nFish),nND(nFish),nEL(nFish)
    real(8),intent(in):: xyzfulIB(1:maxval(Nspan)+1,nFish,nND_max,6)
    real(8),intent(in):: time
    !   -------------------------------------------------------
    integer:: i,j,iFish,Nspanpts,ElmType,off1, off2
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName,idstr
    integer,parameter:: idfile=100
    !==========================================================================
    write(fileName,'(I10)') nint(time*1d5)
    fileName = adjustr(fileName)
    DO  I=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    END DO
    do iFish=1,nFish
        if(Nspan(iFish).eq.0) then
            Nspanpts = 1
            ElmType = ele(1,4,iFish)
        else
            Nspanpts = Nspan(iFish) + 1
            ElmType = 4
        endif
        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        OPEN(idfile,FILE='./DatBodyIB/Body'//trim(idstr)//'_'//trim(filename)//'.dat')
        !   I. The header section.
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        write(idfile, '(A,I7,A,I7,A)', advance='no') 'ZONE N=',nND(iFish)*Nspanpts,', E=',nEL(iFish)*Nspan(iFish),', DATAPACKING=POINT, ZONETYPE='
        if(ElmType.eq.2) then
            write(idfile, '(A)') 'FELINESEG'
        elseif (ElmType.eq.3) then
            write(idfile, '(A)') 'FETRIANGLE'
        elseif(ElmType.eq.4) then
            write(idfile, '(A)') 'FEQUADRILATERAL'
        endif
        do  i=1,nND(iFish)
            do j=1,Nspanpts
                write(idfile, *)   xyzfulIB(j,iFish,i,1:3)
            enddo
        enddo
        do  i=1,nEL(iFish)
            if(ElmType.eq.2) then
                write(idfile, *) ele(i,1,iFish),ele(i,2,iFish)
            elseif(ElmType.eq.3) then
                write(idfile, *) ele(i,1,iFish),ele(i,2,iFish),ele(i,3,iFish)
            else
                do j=1,Nspan(iFish)
                    off1 = (ele(i,1,iFish)-1) * Nspanpts + j
                    off2 = (ele(i,2,iFish)-1) * Nspanpts + j
                    write(idfile, *) off1, off1+1, off2+1, off2
                enddo
            endif
        enddo
        close(idfile)
    enddo
    !   =============================================
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot binary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE write_solid_span_field(xyzful,ele,time,nND,nEL,nND_max,nEL_max,Nspan,dspan,Lref,nFish)
    implicit none
    integer,intent(in):: nND_max,nEL_max,nFish,Nspan(nFish)
    real(8),intent(in):: dspan(nFish)
    integer,intent(in):: ele(nEL_max,5,nFish),nND(nFish),nEL(nFish)
    real(8),intent(in):: xyzful(nND_max,6,nFish)
    real(8),intent(in):: time
    !   -------------------------------------------------------
    integer:: i,j,iFish,Nspanpts,ElmType,off1, off2
    real(8):: Lspan,Lref
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileName,idstr
    integer,parameter:: idfile=100
    !==========================================================================
    write(fileName,'(I10)') nint(time*1d5)
    fileName = adjustr(fileName)
    DO  I=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    END DO
    do iFish=1,nFish
        if(Nspan(iFish).eq.0) then
            Nspanpts = 1
            ElmType = ele(1,4,iFish)
        else
            Nspanpts = Nspan(iFish) + 1
            ElmType = 4
        endif
        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        OPEN(idfile, FILE='./DatBodySpan/BodySpan'//trim(idstr)//'_'//trim(filename)//'.dat')
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        write(idfile, '(A,I7,A,I7,A)', advance='no') 'ZONE N=',nND(iFish)*Nspanpts,', E=',nEL(iFish)*Nspan(iFish),', DATAPACKING=POINT, ZONETYPE='
        if(ElmType.eq.2) then
            write(idfile, '(A)') 'FELINESEG'
        elseif (ElmType.eq.3) then
            write(idfile, '(A)') 'FETRIANGLE'
        elseif(ElmType.eq.4) then
            write(idfile, '(A)') 'FEQUADRILATERAL'
        endif
        do  i=1,nND(iFish)
            do j=1,Nspanpts
                Lspan = (j-1)*dspan(iFish)
                write(idfile, *)  xyzful(i,1:2,iFish),xyzful(i,3,iFish)+Lspan/Lref
            enddo
        enddo
        do  i=1,nEL(iFish)
            if(ElmType.eq.2) then
                write(idfile, *) ele(i,1,iFish),ele(i,2,iFish)
            elseif(ElmType.eq.3) then
                write(idfile, *) ele(i,1,iFish),ele(i,2,iFish),ele(i,3,iFish)
            else
                do j=1,Nspan(iFish)
                    off1 = (ele(i,1,iFish)-1) * Nspanpts + j
                    off2 = (ele(i,2,iFish)-1) * Nspanpts + j
                    write(idfile, *) off1, off1+1, off2+1, off2
                enddo
            endif
        enddo
        close(idfile)
    enddo
    !   =============================================
    END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write parameters for checking
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE  write_params()
    USE simParam
    USE ImmersedBoundary
    implicit none
    integer:: iMT,i,iFish
    integer,parameter::nameLen=10
    character (LEN=nameLen):: fileNameR,fileNameM,fileNameS,fileNameK
    write(fileNameR,'(F10.5)') Re

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
        if    (maxval(iBodyModel(:))==1)then
            write(111,'(A      )')'This is a RIGID    body problem'
        elseif(maxval(iBodyModel(:))==2)then
            write(111,'(A      )')'This is a FLRXIBLE body problem'
        else
            write(111,'(A      )')'This is a FLRXIBLE And RIGID body problem'
        endif
        if    (minval(isMotionGiven(1,:))==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in X-direction freely'
        endif
        if    (minval(isMotionGiven(2,:))==0)then
            write(111,'(A      )')'This FLRXIBLE body can move in Y-direction freely'
        endif
        if    (minval(isMotionGiven(3,:))==0)then
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
        write(111,'(A,2F20.10)') 'elmin,elmax             :',minval(elmin(1:nFish)), maxval(elmax(1:nFish))
        write(111,'(A      )')'===================================='
        write(111,'(A,F20.10)')'Re   =',Re
        write(111,'(A,F20.10)')'Lref =',Lref
        write(111,'(A,F20.10)')'Uref =',Uref
        write(111,'(A,F20.10)')'Tref =',Tref
        write(111,'(A      )')'===================================='
        write(111,'(A,F20.10)')'Freq =',maxval(Freq(:))
        write(111,'(A,F20.10)')'Ampl =',maxval(dabs(XYZAmpl(1:3,1:nFish)))
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

        do iFish=1,nFish
        write(111,'(A      )')'===================================='
        write(111,'(A,I20.10)')'Fish number is',iFish
        write(111,'(A      )')'===================================='
        write(111,'(A,2F20.10)')'Freq,St     =',Freq(iFish),St(iFish)
        write(111,'(A,2F20.10)')'denR,psR    =',denR(iFish),psR(iFish)
        write(111,'(A,2F20.10)')'KB,  KS     =',KB(iFish),KS(iFish)
        write(111,'(A,2F20.10)')'EmR, tcR    =',EmR(iFish),tcR(iFish)
        write(111,'(A,2F20.10)')'theta, span    =',theta(iFish),dspan(iFish) * Nspan(iFish)
        write(111,'(A,1x,3F20.10,2x)')'XYZo(1:3)   =',XYZo(1:3,iFish)
        write(111,'(A,1x,3F20.10,2x)')'XYZAmpl(1:3)=',XYZAmpl(1:3,iFish)
        write(111,'(A,1x,3F20.10,2x)')'XYZPhi(1:3) =',XYZPhi(1:3,iFish)
        write(111,'(A,1x,3F20.10,2x)')'AoAo(1:3)   =',AoAo(1:3,iFish)
        write(111,'(A,1x,3F20.10,2x)')'AoAAmpl(1:3)=',AoAAmpl(1:3,iFish)
        write(111,'(A,1x,3F20.10,2x)')'AoAPhi(1:3) =',AoAPhi(1:3,iFish)
        write(111,'(3(A,1x,I8,2x))')'nND=',nND(iFish),'nEL=',nEL(iFish),'nEQ=',nEQ(iFish)
        write(111,'(3(A,1x,I8,2x))')'nMT=',nMT(iFish),'nBD=',nBD(iFish),'nSTF=',nSTF(iFish)
        enddo

        do iFish=1,nFish
        do iMT=1,nMT(iFish)
        write(111,'(A      )')'===================================='
        write(111,'(A,I5.5 )')'Fish number is',iFish
        write(111,'(A,I5.5 )')'MT:',iMT
        write(111,'(A,E20.10 )')'E    =',prop(iMT,1,iFish)
        write(111,'(A,E20.10 )')'G    =',prop(iMT,2,iFish)
        write(111,'(A,E20.10 )')'h    =',prop(iMT,3,iFish)
        write(111,'(A,E20.10 )')'rho  =',prop(iMT,4,iFish)
        write(111,'(A,E20.10 )')'gamma=',prop(iMT,5,iFish)
        write(111,'(A,E20.10 )')'Ip   =',prop(iMT,6,iFish)
        write(111,'(A,E20.10 )')'alpha=',prop(iMT,7,iFish)
        write(111,'(A,E20.10 )')'beta =',prop(iMT,8,iFish)
        enddo
        enddo
        close(111)
    END SUBROUTINE



!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    compute strain energy
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE strain_energy_D(  strainEnergy, xord0,yord0,zord0,xord,yord,zord, &
                                 ele,prop,triad_n1,triad_n2,triad_ee, &
                                 nND,nEL,nMT &
                               )
    implicit none
    integer:: nMT,nEL,nND
    integer:: ele(nEL,5)
    real(8):: xord0(nND), yord0(nND), zord0(nND), strainEnergy(nEL,2)
    real(8):: xord(nND), yord(nND), zord(nND)
!
    real(8):: prop(nMT,10)

    real(8):: ekb12(12,12),ekb12Strech(12,12),ekb12BendTor(12,12)
!
    real(8):: triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)

    real(8):: triad_00(3,3),triad_11(3,3),triad_22(3,3)
    real(8):: ub(18),dl
!
    real(8):: dx0,dy0,dz0,du,dv,dw,dx,dy,dz,xl0
    real(8):: tx1,tx2,ty1,ty2,tz1,tz2,tx,ty,tz

    real(8):: e0,g0,a0,b0,r0,zix0,ziy0,ziz0,xl
    integer:: i,j,n,i1,j1,k1,mat,nELt

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
            ub(1)=0.0d0
            ub(2)=0.0d0
            ub(3)=0.0d0
            ub(4)=tx1
            ub(5)=ty1
            ub(6)=tz1
!
            ub(7)=dl
            ub(8)=0.0d0
            ub(9)=0.0d0
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
            strainEnergy(n,1)=0.5d0*sum(matmul(ekb12Strech(1:12,1:12),ub(1:12))*ub(1:12))
            strainEnergy(n,2)=0.5d0*sum(matmul(ekb12BendTor(1:12,1:12),ub(1:12))*ub(1:12))
        else
            !write(*,*)'not this nELt:',nELt
            !stop
            strainEnergy(n,:) = 0.0d0
        endif

    enddo



    return
    ENDSUBROUTINE strain_energy_D

SUBROUTINE write_flow_fast()
USE simParam
USE OutFlowWorkspace
implicit none
integer:: x,y,z,pid,i
integer::xmin,xmax,ymin,ymax,zmin,zmax
integer,parameter::nameLen=10,idfile=100
character (LEN=nameLen):: fileName
real(8):: invUref
real(8):: CPUtime, waittime
integer,dimension(8) :: values0,values1
call date_and_time(VALUES=values0)
call mywait()
call date_and_time(VALUES=values1)
waittime = CPUtime(values1)-CPUtime(values0)
if(waittime.gt.1.d-1) then
    write(*,'(A,F7.2,A)')'Waiting ', waittime, 's for previous outflow finishing.'
endif
xmin = 1 + offsetOutput
ymin = 1 + offsetOutput
zmin = 1 + offsetOutput
xmax = xDim - offsetOutput
ymax = yDim - offsetOutput
zmax = zDim - offsetOutput
invUref = 1.d0/Uref

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
do x=xmin, xmax
    do y=ymin, ymax
        do z=zmin, zmax
            oututmp(z,y,x) = uuu(z,y,x,1)*invUref
        enddo
    enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
do x=xmin, xmax
    do y=ymin, ymax
        do z=zmin, zmax
            outvtmp(z,y,x) = uuu(z,y,x,2)*invUref
        enddo
    enddo
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
do x=xmin, xmax
    do y=ymin, ymax
        do z=zmin, zmax
            outwtmp(z,y,x) = uuu(z,y,x,3)*invUref
        enddo
    enddo
enddo
!$OMP END PARALLEL DO
offsetMoveGrid=0.0
if(isMoveGrid==1)then
    if(isMoveDimX==1) offsetMoveGrid(1) = dh*dble(MoveOutputIref(1))
    if(isMoveDimY==1) offsetMoveGrid(2) = dh*dble(MoveOutputIref(2))
    if(isMoveDimZ==1) offsetMoveGrid(3) = dh*dble(MoveOutputIref(3))
endif
call myfork(pid)
if(pid.eq.0) then
    write(fileName,'(I10)') nint(time/Tref*1d5)
    fileName = adjustr(fileName)
    do  i=1,nameLen
        if(fileName(i:i)==' ')fileName(i:i)='0'
    enddo
    open(idfile,file='./DatFlow/Flow'//trim(fileName)//'.plt',form='unformatted',access='stream')
    WRITE(idfile) xmin,xmax,ymin,ymax,zmin,zmax
    WRITE(idfile) offsetMoveGrid(1:3)
    write(idfile) oututmp,outvtmp,outwtmp
    close(idfile)
    call myexit(0)
endif
END SUBROUTINE

subroutine ComputeFieldStat
USE simParam
USE OutFlowWorkspace
implicit none
integer:: x,y,z,i
real(8):: invUref, uLinfty(1:3), uL2(1:3), temp
invUref = 1.d0/Uref
uLinfty = -1.d0
uL2 = 0.d0
do i=1,3
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,temp) &
    !$OMP reduction (+: uL2) reduction(max: uLinfty)
    do x=1, xDim
        do y=1, yDim
            do z=1, zDim
                temp = dabs(uuu(z,y,x,i)*invUref)
                uL2(i) = uL2(i) + temp * temp
                if(temp.gt.uLinfty(i)) uLinfty(i) = temp
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO
    uL2(i) = dsqrt(uL2(i) / dble(xDim * yDim * zDim))
enddo
write(*,'(A,F18.12)')'FIELDSTAT L2 u ', uL2(1)
write(*,'(A,F18.12)')'FIELDSTAT L2 v ', uL2(2)
write(*,'(A,F18.12)')'FIELDSTAT L2 w ', uL2(3)
write(*,'(A,F18.12)')'FIELDSTAT Linfinity u ', uLinfty(1)
write(*,'(A,F18.12)')'FIELDSTAT Linfinity v ', uLinfty(2)
write(*,'(A,F18.12)')'FIELDSTAT Linfinity w ', uLinfty(3)
endsubroutine ComputeFieldStat

subroutine initOutFlowWorkspace
use simParam
use OutFlowWorkspace
implicit none
integer::xmin,xmax,ymin,ymax,zmin,zmax
xmin = 1 + offsetOutput
ymin = 1 + offsetOutput
zmin = 1 + offsetOutput
xmax = xDim - offsetOutput
ymax = yDim - offsetOutput
zmax = zDim - offsetOutput
allocate( oututmp(zmin:zmax,ymin:ymax,xmin:xmax),outvtmp(zmin:zmax,ymin:ymax,xmin:xmax),outwtmp(zmin:zmax,ymin:ymax,xmin:xmax) )
endsubroutine initOutFlowWorkspace
