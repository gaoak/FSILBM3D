!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write parameters for checking
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE  write_params()
    USE simParam
    USE SolidBody
    implicit none
    integer:: i,iFish
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
