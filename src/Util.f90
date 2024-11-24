!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write files' header
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wrtInfoTitl()
    USE simParam
    USE SolidBody
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
        open(111,file='./DatInfo/ForceDirect'//trim(filename)//'.plt')
        write(111,*)'variables= "t"  "Fx"  "Fy"  "Fz"'
        close(111)
        elseif(iForce2Body==2)then   !stress force
        open(111,file='./DatInfo/ForceStress'//trim(filename)//'.plt')
        write(111,*)'variables= "t"  "Fx"  "Fy"  "Fz"'
        close(111)
        endif

        !===============================================================================
        open(111,file='./DatInfo/SampBodyNodeBegin'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)
        !===============================================================================
        write(Nodename,'(I4.4)') VBodies(iFish)%rbm%nNd
        open(111,file='./DatInfo/SampBodyNodeEnd'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)

        open(111,file='./DatInfo/SampBodyMean'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)

        open(111,file='./DatInfo/SampBodyAngular'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "AoA"  "Ty-Hy"  "Hy"  "Ty"'
        close(111)

        open(111,file='./DatInfo/Power'//trim(fileName)//'.plt')
        write(111,*)'variables= "t" "Ptot" "Paero" "Piner" "Pax" "Pay" "Paz" "Pix" "Piy" "Piz"'
        close(111)

        !open(111,file='./DatInfo/Area'//trim(fileName)//'.plt')
        !write(111,*)'variables= "t"  "Area"  '
        !close(111)

        open(111,file='./DatInfo/Energy'//trim(fileName)//'.plt')
        write(111,*)'variables= "t","Es","Eb","Ep","Ek","Ew","Et"'
        close(111)
    enddo

    !open(111,file='./DatInfo/MaxValBody.plt')
    !write(111,*)'variables= "t"  "xyzMax" "velMax" "accMax"  '
    !close(111)

    open(111,file='./DatInfo/MaMax.plt')
    write(111,*)'variables= "t"  "MaMax"  '
    close(111)

    open(111,file='./DatInfo/Converg.plt')
    write(111,*)'variables= "t"  "Convergence"  '
    close(111)

    if(isBodyOutput==1)then
    do  i=1,numSampBody
        write(Nodename,'(I4.4)') SampBodyNode(i,iFish)
        open(111,file='./DatInfo/SampBodyNode'//trim(fileName)//'_'//trim(Nodename)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)
    enddo
    endif

    if(isFluidOutput==1)then
    do  i=1,numSampFlow
        write(Nodename,'(I4.4)') i
        open(111,file='./DatInfo/SampFlowPint'//trim(Nodename)//'.plt')
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
    USE simParam
    USE SolidBody
    implicit none
    integer:: i,z,y,x,zbgn,ybgn,xbgn,zend,yend,xend,iFish
    real(8):: convergence,MaMax,weightm,velocity(1:3),Pressure
    integer,parameter::nameLen=4
    character (LEN=nameLen):: fileName,Nodename

    call Write_solid_Data(111,time,timeOutInfo,Asfac)

    UNow=sum(dsqrt( uuu(:,:,:,1)**2+uuu(:,:,:,2)**2+uuu(:,:,:,3)**2))
    convergence=dabs(UNow-UPre)/UNow
    UPre=UNow
    open(111,file='./DatInfo/Converg.plt',position='append')
    write(111,'(2E20.10)')time/Tref,convergence
    close(111)

    MaMax=MaxVal(dsqrt( uuu(:,:,:,1)**2+uuu(:,:,:,2)**2+uuu(:,:,:,3)**2))/dsqrt(Cs2)
    open(111,file='./DatInfo/MaMax.plt',position='append')
    write(111,'(2E20.10)')time/Tref,MaMax
    close(111)

    !open(111,file='./DatInfo/MaxValBody.plt',position='append')
    !write(111,'(4E20.10)')time/Tref,maxval(dsqrt(xyzful(1:nND_max,1,1:nFish)**2+xyzful(1:nND_max,2,1:nFish)**2+xyzful(1:nND_max,3,1:nFish)**2))/Lref, &
    !                                maxval(dsqrt(velful(1:nND_max,1,1:nFish)**2+velful(1:nND_max,2,1:nFish)**2+velful(1:nND_max,3,1:nFish)**2))/Uref, &
    !                                maxval(dsqrt(accful(1:nND_max,1,1:nFish)**2+accful(1:nND_max,2,1:nFish)**2+accful(1:nND_max,3,1:nFish)**2))/Aref
    !close(111)

    if(isBodyOutput==1)then
    do  i=1,numSampBody
        write(Nodename,'(I4.4)') SampBodyNode(i,iFish)
        open(111,file='./DatInfo/SampBodyNode'//trim(fileName)//'_'//trim(Nodename)//'.plt',position='append')
        write(111,'(10E20.10)')time/Tref, VBodies(iFish)%rbm%xyzful(SampBodyNode(i,iFish),1:3)/Lref,VBodies(iFish)%rbm%velful(SampBodyNode(i,iFish),1:3)/Uref,VBodies(iFish)%rbm%accful(SampBodyNode(i,iFish),1:3)/Aref
        close(111)
    enddo
    endif

    if(isFluidOutput==1)then
    do  i=1,numSampFlow
        x=minloc(dabs(SampFlowPint(i,1)-xGrid(1:xDim)),1)
        if(SampFlowPint(i,1)-xGrid(x)>0.0d0)then
                xbgn=x; xend=x+1
        else
                xbgn=x-1;xend=x
        endif
        y=minloc(dabs(SampFlowPint(i,2)-yGrid(1:yDim)),1)
        if(SampFlowPint(i,2)-yGrid(y)>0.0d0)then
                ybgn=y; yend=y+1
        else
                ybgn=y-1;yend=y
        endif
        z=minloc(dabs(SampFlowPint(i,3)-zGrid(1:zDim)),1)
        if(SampFlowPint(i,3)-zGrid(z)>0.0d0)then
                zbgn=z; zend=z+1
        else
                zbgn=z-1;zend=z
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

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE movGrid(dim,direction)
    USE simParam
    implicit none
    integer:: dim,direction

    if    (dim==1)then      !x
        if(direction<0)       then      !move to -x
            fIn(:,:,2:xDim,:)=fIn(:,:,1:xDim-1,:)
            xGrid(:)=xGrid(:)-dx
        elseif(direction>0)    then     !move to +x
            fIn(:,:,1:xDim-1,:)=fIn(:,:,2:xDim,:)
            xGrid(:)=xGrid(:)+dx
        else
        endif
    elseif(dim==2)then      !y
        if(direction<0)        then     !move to -y
            fIn(:,2:yDim,:,:)=fIn(:,1:yDim-1,:,:)
            yGrid(:)=yGrid(:)-dy
        elseif(direction>0)    then     !move to +y
            fIn(:,1:yDim-1,:,:)=fIn(:,2:yDim,:,:)
            yGrid(:)=yGrid(:)+dy
        else
        endif
    elseif(dim==3)then      !z
        if(direction<0)        then     !move to -z
            fIn(2:zDim,:,:,:)=fIn(1:zDim-1,:,:,:)
            zGrid(:)=zGrid(:)-dz
        elseif(direction>0)    then     !move to +z
            fIn(1:zDim-1,:,:,:)=fIn(2:zDim,:,:,:)
            zGrid(:)=zGrid(:)+dz
        else
        endif
    else
    endif
    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Counting the number of steps the grid has moved
!    Determine if the grd is moving left or right
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  cptMove(move,xB,xG,d,MoveOutputIref)
    implicit none
    integer:: move(3),i,MoveOutputIref(3)
    real(8):: xB(3),xG(3),d
    do  i=1,3
        if    (xB(i)-xG(i)> d)then
            move(i)=1
            MoveOutputIref(i)=MoveOutputIref(i)+1
        elseif(xB(i)-xG(i)<-d)then
            move(i)=-1
            MoveOutputIref(i)=MoveOutputIref(i)-1
        else
            move(i)=0
        endif
    enddo
    end subroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    find reference point, moving grid
!    copyright@ RuNanHua
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cptIref(NDref,IXref,IYref,IZref,nND,xDim,yDim,zDim,xyzful,xGrid,yGrid,zGrid,Xref,Yref,Zref)
    implicit none
    integer:: NDref,IXref,IYref,IZref,nND,xDim,yDim,zDim
    real(8):: xyzful(1:nND,1:3),xGrid(xDim),yGrid(yDim),zGrid(zDim),Xref,Yref,Zref

    NDref=minloc(dsqrt((xyzful(1:nND,1)-Xref)**2+(xyzful(1:nND,2)-Yref)**2+(xyzful(1:nND,3)-Zref)**2),1)
    IXref=minloc(dabs(xyzful(NDref,1)-xGrid(1:xDim)),1)
    IYref=minloc(dabs(xyzful(NDref,2)-yGrid(1:yDim)),1)
    IZref=minloc(dabs(xyzful(NDref,3)-zGrid(1:zDim)),1)

    end subroutine
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  initDisturb()
    USE simParam
    implicit none

    integer:: x, y,z
    do  z = 1, zDim
    do  y = 1, yDim
    do  x = 1, xDim
        uuu(z,y,x,1)=uuu(z,y,x,1)+AmplInitDist(1)*Uref*dsin(2.0*pi*waveInitDist*xGrid(x))
        uuu(z,y,x,2)=uuu(z,y,x,2)+AmplInitDist(2)*Uref*dsin(2.0*pi*waveInitDist*xGrid(x))
        uuu(z,y,x,3)=uuu(z,y,x,3)+AmplInitDist(3)*Uref*dsin(2.0*pi*waveInitDist*xGrid(x))
    enddo
    enddo
    enddo

    end

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  forcDisturb()
    USE simParam
    implicit none
    integer:: x, y, z,xbgn,xend,ybgn,yend,zbgn,zend
    real(8):: rx,ry,rz

    xbgn    = minloc(dabs(posiForcDist(1)-xGrid(1:xDim)),1)-3
    xend    = minloc(dabs(posiForcDist(1)-xGrid(1:xDim)),1)+4
    ybgn    = minloc(dabs(posiForcDist(2)-yGrid(1:yDim)),1)-3
    yend    = minloc(dabs(posiForcDist(2)-yGrid(1:yDim)),1)+4
    zbgn    = minloc(dabs(posiForcDist(3)-zGrid(1:zDim)),1)-3
    zend    = minloc(dabs(posiForcDist(3)-zGrid(1:zDim)),1)+4
    do    x=xbgn,xend
    do    y=ybgn,yend
    do    z=zbgn,zend
        rx=(posiForcDist(1)-xGrid(x))/dx(x)
        ry=(posiForcDist(2)-yGrid(y))/dy(y)
        rz=(posiForcDist(3)-zGrid(z))/dz(z)

        force(z,y,x,1)=force(z,y,x,1)+AmplforcDist(1)*(0.5*denIn*Uref**2*Asfac)*dsin(2.0*pi*FreqforcDist*time/Tref)*Phi(rx)*Phi(ry)*Phi(rz)/(dx(x)*dy(y)*dz(z))
        force(z,y,x,2)=force(z,y,x,2)+AmplforcDist(2)*(0.5*denIn*Uref**2*Asfac)*dsin(2.0*pi*FreqforcDist*time/Tref)*Phi(rx)*Phi(ry)*Phi(rz)/(dx(x)*dy(y)*dz(z))
        force(z,y,x,3)=force(z,y,x,3)+AmplforcDist(3)*(0.5*denIn*Uref**2*Asfac)*dsin(2.0*pi*FreqforcDist*time/Tref)*Phi(rx)*Phi(ry)*Phi(rz)/(dx(x)*dy(y)*dz(z))
    enddo
    enddo
    enddo
    contains
    FUNCTION Phi(x)
        IMPLICIT NONE
        real(8)::Phi,x,r
        r=dabs(x)
        if(r<1.0d0)then
            Phi=(3.d0-2.d0*r+dsqrt( 1.d0+4.d0*r*(1.d0-r)))*0.125d0
        elseif(r<2.0d0)then
            Phi=(5.d0-2.d0*r-dsqrt(-7.d0+4.d0*r*(3.d0-r)))*0.125d0
        else
            Phi=0.0d0
        endif
    ENDFUNCTION Phi
    end

    SUBROUTINE OMPPartition(xDim, np, partition, parindex)
        implicit none
        integer:: np, xDim
        integer:: partition(1:np), parindex(1:np+1)
        integer:: psize, p, residual
        psize = xDim/np
        residual = xDim - psize * np
        parindex(1) = 1
        parindex(np+1) = xDim + 1
        do p=1,np
            if (p .gt. np-residual) then
                partition(p) = psize + 1
            else
                partition(p) = psize
            endif
            if (p .gt. 1) then
                parindex(p) = parindex(p-1) + partition(p-1)
            endif
        enddo
    endsubroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   lxguang 2023.02 Add Shear flow velocity boundary
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE evaluateShearVelocity(x, y, z, vel)
    USE simParam
    real(8):: x, y, z, vel(1:SpcDim)
    vel(1) = uuuIn(1) + 0 * shearRateIn(1) + y * shearRateIn(2) + z * shearRateIn(3)
    vel(2) = uuuIn(2) + x * shearRateIn(1) + 0 * shearRateIn(2) + z * shearRateIn(3)
    vel(3) = uuuIn(3) + x * shearRateIn(1) + y * shearRateIn(2) + 0 * shearRateIn(3)
    END SUBROUTINE

    SUBROUTINE evaluateOscillatoryVelocity(vel)
        USE simParam
        real(8):: vel(1:SpcDim)
        vel(1) = uuuIn(1) + VelocityAmp * dcos(2*pi*VelocityFreq*time + VelocityPhi/180.0*pi)
        vel(2) = uuuIn(2)
        vel(3) = uuuIn(3)
    END SUBROUTINE

    SUBROUTINE updateVolumForc()
        USE simParam
        implicit none
        VolumeForce(1) = VolumeForceIn(1) + VolumeForceAmp * dsin(2.d0 * pi * VolumeForceFreq * time + VolumeForcePhi/180.0*pi)
        VolumeForce(2) = VolumeForceIn(2)
        VolumeForce(3) = VolumeForceIn(3)
    END SUBROUTINE

    SUBROUTINE addVolumForc()
        USE simParam
        implicit none
        integer:: x
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
        do x=1,xDim
            force(:,:,x,1) = force(:,:,x,1) + VolumeForce(1)
            force(:,:,x,2) = force(:,:,x,2) + VolumeForce(2)
            force(:,:,x,3) = force(:,:,x,3) + VolumeForce(3)
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE

    FUNCTION CPUtime(values)
        IMPLICIT NONE
        real(8)::CPUtime
        integer,dimension(8) :: values
        CPUtime = dble(values(6))*60.d0+dble(values(7))*1.d0+dble(values(8))*0.001d0
    ENDFUNCTION