!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Calculate the repulsive force between solids
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptForceR(dxmin,dymin,dzmin,nND,nND_max,xyzful,repful,nFish)
    USE ImmersedBoundary
    implicit none
    integer:: nFish,nND_max,nND(1:nFish)
    real(8):: dxmin,dymin,dzmin
    real(8):: xyzful(1:nND_max,1:6,1:nFish),repful(1:nND_max,1:6,1:nFish)
    !local
    integer:: iND,jND,iFish,jFish
    real(8):: delta_h,r(1:3),ds(1:3),phi_r(1:3),SpanLength
    real(8):: minx,miny,maxx,maxy,minz,maxz
    real(8):: xmin(1:nFish),xmax(1:nFish),ymin(1:nFish),ymax(1:nFish),zmin(1:nFish),zmax(1:nFish)

    repful(:,:,:)=0.d0
    ds(1)=dxmin
    ds(2)=dymin
    ds(3)=dzmin
    SpanLength=dspan*Nspan
    if(Nspan.eq.0)then
        SpanLength=1.0d0
    endif

    do iFish=1,nFish
        xmin(iFish) = minval(xyzful(1:nND(iFish),1,iFish))-dxmin*2.5d0
        xmax(iFish) = maxval(xyzful(1:nND(iFish),1,iFish))+dxmin*2.5d0
        ymin(iFish) = minval(xyzful(1:nND(iFish),2,iFish))-dymin*2.5d0
        ymax(iFish) = maxval(xyzful(1:nND(iFish),2,iFish))+dymin*2.5d0
        zmin(iFish) = minval(xyzful(1:nND(iFish),3,iFish))-dzmin*2.5d0
        zmax(iFish) = maxval(xyzful(1:nND(iFish),3,iFish))+dzmin*2.5d0
    enddo

    do iFish=1,nFish
        do jFish=iFish+1,nFish
            minx = max(xmin(iFish),xmin(jFish))
            miny = max(ymin(iFish),ymin(jFish))
            minz = max(zmin(iFish),zmin(jFish))
            maxx = min(xmax(iFish),xmax(jFish))
            maxy = min(ymax(iFish),ymax(jFish))
            maxz = min(zmax(iFish),zmax(jFish))
            if((minx > maxx).or.(miny > maxy).or.(minz > maxz)) then
                cycle
            endif
            ! overlapping regin [minx, maxx] X [miny, maxy] X [minz, maxz]
            do iND=1,nND(iFish)
                if ( (minx>xyzful(iND,1,iFish)) .or. (xyzful(iND,1,iFish)>maxx)  &
                .or. (miny>xyzful(iND,2,iFish)) .or. (xyzful(iND,2,iFish)>maxy)  &
                .or. (minz>xyzful(iND,3,iFish)) .or. (xyzful(iND,3,iFish)>maxz)) then
                    cycle !point iND not in the overpalling region
                endif
                do jND=1,nND(jFish)
                    if ( (minx>xyzful(jND,1,jFish)) .or. (xyzful(jND,1,jFish)>maxx)  &
                    .or. (miny>xyzful(jND,2,jFish)) .or. (xyzful(jND,2,jFish)>maxy)  &
                    .or. (minz>xyzful(jND,3,jFish)) .or. (xyzful(jND,3,jFish)>maxz)) then
                        cycle !point jND not in the overpalling region
                    endif
                    r(1)=(xyzful(iND,1,iFish)-xyzful(jND,1,jFish))/dxmin
                    r(2)=(xyzful(iND,2,iFish)-xyzful(jND,2,jFish))/dymin
                    r(3)=(xyzful(iND,3,iFish)-xyzful(jND,3,jFish))/dzmin
                    call get_smoother_phi_r(r,phi_r)
                    delta_h=phi_r(1)*phi_r(2)*phi_r(3)/dxmin/dymin/dzmin/sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
                    repful(iND,1:3,iFish)=repful(iND,1:3,iFish) + delta_h*r(1:3)*ds(1:3)*SpanLength ! force
                    repful(jND,1:3,jFish)=repful(jND,1:3,jFish) - delta_h*r(1:3)*ds(1:3)*SpanLength ! reaction force
                enddo !jND=1,nND(jFish)
            enddo !iND=1,nND(iFish)
        enddo !jFish=iFish+1,nFish
    enddo !iFish=1,nFish

    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    (displacement, velocity) spring, penalty method
!    calculate force at element center, distribute force to three nodes
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_interaction_force(zDim,yDim,xDim,nEL,nND,ele,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                                xyzful,velful,xyzfulIB,Palpha,Pbeta,ntolLBM,dtolLBM,force,extful,isUniformGrid)
USE, INTRINSIC :: IEEE_ARITHMETIC
IMPLICIT NONE
integer,intent(in):: zDim,yDim,xDim,nEL,nND,ele(nEL,5),ntolLBM
real(8),intent(in):: dh,Uref,denIn,dtolLBM,dt,Palpha,Pbeta
real(8),intent(in):: den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
logical,intent(in):: isUniformGrid(1:3)
real(8),intent(in):: xyzful(nND,6),velful(nND,6)
real(8),intent(inout)::xyzfulIB(nND,6),uuu(zDim,yDim,xDim,1:3)
real(8),intent(out)::extful(nND,6),force(zDim,yDim,xDim,1:3)
!==================================================================================================
integer:: i,j,k,x,y,z,iEL,nt,iterLBM,iND
real(8):: rx,ry,rz,Phi,dmaxLBM,dsum,invdh,forcetemp(1:3)
real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az
real(8):: forceElemTemp(nEL,3)
real(8):: forceElem(nEL,3),forceNode(nND,3),velfulIB(nND,3)
real(8):: posElem(nEL,3),velElem(nEL,3),velElemIB(nEL,3),areaElem(nEL)
real(8),allocatable::posElemIB(:,:)
!==================================================================================================
!   compute velocity and displacement at IB nodes
invdh = 1.D0/dh
if(Palpha.gt.0.d0) then
    allocate(posElemIB(1:nEL,1:3))
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iND,i,j,k,x,y,z,rx,ry,rz)
    do  iND=1,nND
        call my_minloc(xyzful(iND,1), xGrid, xDim, isUniformGrid(1), i)
        call my_minloc(xyzful(iND,2), yGrid, yDim, isUniformGrid(2), j)
        call my_minloc(xyzful(iND,3), zGrid, zDim, isUniformGrid(3), k)
        ! xbgn             xxx             xend
        !   -1   0   1   2
        ! interpolate fluid velocity to body nodes
        velfulIB(iND,1:3)=0.0d0
        do x=-1+i,2+i
            rx=Phi((xyzful(iND,1)-xGrid(x))*invdh)
            do y=-1+j,2+j
                ry=Phi((xyzful(iND,2)-yGrid(y))*invdh)
                do z=-1+k,2+k
                    rz=Phi((xyzful(iND,3)-zGrid(z))*invdh)
                    velfulIB(iND,1:3)=velfulIB(iND,1:3)+uuu(z,y,x,1:3)*rx*ry*rz
                enddo
            enddo
        enddo
        xyzfulIB(iND,1:3)=xyzfulIB(iND,1:3)+velfulIB(iND,1:3)*dt
    enddo
    !$OMP END PARALLEL DO
endif
!==================================================================================================
!   compute displacement, velocity, area at surface element center
do  iEL=1,nEL
    i=ele(iEL,1)
    j=ele(iEL,2)
    k=ele(iEL,3)
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
        posElem(iEL,1:3)=(xyzful(i,1:3)+xyzful(j,1:3))/2.0d0
        velElem(iEL,1:3)=(velful(i,1:3)+velful(j,1:3))/2.0d0
        if(Palpha.gt.0.d0) then
            posElemIB(iEL,1:3)=(xyzfulIB(i,1:3)+xyzfulIB(j,1:3))/2.0d0
        endif
        ax =(x1-x2)
        ay =(y1-y2)
        az =(z1-z2)
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az)

    elseif(nt==3)then
        posElem(iEL,1:3)=(xyzful(i,1:3)+xyzful(j,1:3)+xyzful(k,1:3))/3.0d0
        velElem(iEL,1:3)=(velful(i,1:3)+velful(j,1:3)+velful(k,1:3))/3.0d0
        if(Palpha.gt.0.d0) then
            posElemIB(iEL,1:3)=(xyzfulIB(i,1:3)+xyzfulIB(j,1:3)+xyzfulIB(k,1:3))/3.0d0
        endif
        ax =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.0d0
        ay =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.0d0
        az =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.0d0
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az)
    else
            write(*,*)'cell type is not defined'
    endif
enddo

!**************************************************************************************************
!**************************************************************************************************
forceElem(1:nEL,1:3)=0.0d0
dmaxLBM=1.0d0
iterLBM=0
!   ***********************************************************************************************
do  while( iterLBM<ntolLBM .and. dmaxLBM>dtolLBM)
!   ***********************************************************************************************
!   compute the velocity of IB nodes at element center
    do  iEL=1,nEL
        call my_minloc(posElem(iEL,1), xGrid, xDim, isUniformGrid(1), i)
        call my_minloc(posElem(iEL,2), yGrid, yDim, isUniformGrid(2), j)
        call my_minloc(posElem(iEL,3), zGrid, zDim, isUniformGrid(3), k)
        velElemIB(iEL,1:3)=0.0d0
        do x=-1+i,2+i
            rx=Phi((posElem(iEL,1)-xGrid(x))*invdh)
            do y=-1+j,2+j
                ry=Phi((posElem(iEL,2)-yGrid(y))*invdh)
                do z=-1+k,2+k
                    rz=Phi((posElem(iEL,3)-zGrid(z))*invdh)
                    velElemIB(iEL,1:3)=velElemIB(iEL,1:3)+uuu(z,y,x,1:3)*rx*ry*rz
                enddo
            enddo
        enddo
    enddo
!   ***********************************************************************************************
!   calculate interaction force
    do  iEL=1,nEL
        if(ele(iEL,4)==3) then
            forceElemTemp(iEL,1:3) = -Pbeta* 2.0d0*denIn*(velElem(iEL,1:3)-velElemIB(iEL,1:3))/dt*areaElem(iEL)*dh
            if(Palpha.gt.0.d0) then
                forceElemTemp(iEL,1:3) = forceElemTemp(iEL,1:3) - Palpha*2.0d0*denIn*(posElem(iEL,1:3)-posElemIB(iEL,1:3))/dt*areaElem(iEL)*dh
            endif
        else
            forceElemTemp(iEL,1:3) = 0.0d0
        endif
        if ((.not. IEEE_IS_FINITE(forceElemTemp(iEL,1))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(iEL,2))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(iEL,3)))) then
            write(*, *) 'Nan found in forceElemTemp', forceElemTemp
            stop
        endif
    enddo
!   ***********************************************************************************************
!   calculate Eulerian body force
    do    iEL=1,nEL
        call my_minloc(posElem(iEL,1), xGrid, xDim, isUniformGrid(1), i)
        call my_minloc(posElem(iEL,2), yGrid, yDim, isUniformGrid(2), j)
        call my_minloc(posElem(iEL,3), zGrid, zDim, isUniformGrid(3), k)
        do x=-1+i,2+i
            rx=Phi((posElem(iEL,1)-xGrid(x))*invdh)
            do y=-1+j,2+j
                ry=Phi((posElem(iEL,2)-yGrid(y))*invdh)
                do z=-1+k,2+k
                    rz=Phi((posElem(iEL,3)-zGrid(z))*invdh)
                    forcetemp(1:3) = -forceElemTemp(iEL,1:3)*rx*ry*rz*invdh*invdh*invdh
                    ! update velocity
                    uuu(z,y,x,1:3)  = uuu(z,y,x,1:3)+0.5d0*dt*forceTemp(1:3)/den(z,y,x)
                    force(z,y,x,1:3)=force(z,y,x,1:3) + forcetemp(1:3)
                enddo
            enddo
        enddo
    enddo
    forceElem(1:nEL,1:3) = forceElem(1:nEL,1:3)+forceElemTemp(1:nEL,1:3)
!   convergence test
    if(iterLBM==0)then
        dsum=0.0d0
        do iEL=1,nEL
        dsum=dsum+dsqrt(sum((velElem(iEL,1:3)-velElemIB(iEL,1:3))**2))
        enddo
    endif
    dsum=Uref*nEL

    dmaxLBM=0.0d0
    do iEL=1,nEL
        dmaxLBM=dmaxLBM+dsqrt(sum((velElem(iEL,1:3)-velElemIB(iEL,1:3))**2))
    enddo
    dmaxLBM=dmaxLBM/dsum
    iterLBM=iterLBM+1
!   ***********************************************************************************************
enddo
!write(*,'(A,I5,A,D20.10)')' iterLBM =',iterLBM,'    dmaxLBM =',dmaxLBM
!**************************************************************************************************
!**************************************************************************************************
!   element force to nodal force
forceNode(1:nND,1:3)=0.0d0
do    iEL=1,nEL
    i=ele(iEL,1)
    j=ele(iEL,2)
    k=ele(iEL,3)
    nt=ele(iEL,4)
    forceNode(i,1:3)=forceNode(i,1:3)+forceElem(iEl,1:3)/3.0d0
    forceNode(j,1:3)=forceNode(j,1:3)+forceElem(iEl,1:3)/3.0d0
    forceNode(k,1:3)=forceNode(k,1:3)+forceElem(iEl,1:3)/3.0d0
enddo

extful(1:nND,1:3) = forceNode(1:nND,1:3)
extful(1:nND,4:6) = 0.0d0
if(Palpha.gt.0.d0) then
    deallocate(posElemIB)
endif
END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    (displacement, velocity) spring, penalty method
!    calculate force at element center, distribute force to four nodes
!    solid should be in uniform grid
!    period body should span full from 1 to xDim+1
!    symmetric body should not exceed domain
!    the first body point should be in the domain or on the boundary
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_interaction_force_quad(zDim,yDim,xDim,nEL,nND,ele,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
    xyzful,velful,xyzfulIB,Palpha,Pbeta,ntolLBM,dtolLBM,force,extful,isUniformGrid,Nspan,dspan,boundaryConditions)
USE, INTRINSIC :: IEEE_ARITHMETIC
use BoundCondParams
IMPLICIT NONE
integer,intent(in):: zDim,yDim,xDim,nEL,nND,ele(nEL,5),ntolLBM,Nspan
real(8),intent(in):: dh,Uref,denIn,dtolLBM,dt,Palpha,Pbeta,dspan
real(8),intent(in):: den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
logical,intent(in):: isUniformGrid(1:3)
real(8),intent(in):: xyzful(nND,6),velful(nND,6)
integer,intent(in):: boundaryConditions(1:6)
real(8),intent(inout)::xyzfulIB(1:Nspan+1,nND,6),uuu(zDim,yDim,xDim,1:3)
real(8),intent(out)::extful(nND,6),force(zDim,yDim,xDim,1:3)
!==================================================================================================
integer:: i,j,k,x,y,z,s,iEL,nt,iterLBM,iND
integer:: ix(-1:2),jy(-1:2),kz(-1:2)
real(8):: rx(-1:2),ry(-1:2),rz(-1:2),Phi,dmaxLBM,dsum,invdh,forceTemp(1:3)
real(8):: x1,x2,y1,y2,z1,ax,ay
real(8):: velfulIB(1:Nspan+1,nND,3)
real(8):: forceElem(1:Nspan,nEL,3),forceElemTemp(1:Nspan,nEL,3),areaElem(nEL)
real(8):: posElem(1:Nspan,nEL,3),velElem(1:Nspan,nEL,3),velElemIB(1:Nspan,nEL,3)
real(8),allocatable::posElemIB(:,:,:)
real(8)::x0,y0,z0,detx,dety,detz
integer::i0,j0,k0
!==================================================================================================
invdh = 1.D0/dh
call my_minloc(xyzful(1,1), xGrid, xDim, isUniformGrid(1), i0)
call my_minloc(xyzful(1,2), yGrid, yDim, isUniformGrid(2), j0)
call my_minloc(xyzful(1,3), zGrid, zDim, isUniformGrid(3), k0)
x0 = xGrid(i0)
y0 = yGrid(j0)
z0 = zGrid(k0)
!   compute velocity and displacement at IB nodes
if(Palpha.gt.0.d0) then
    allocate(posElemIB(1:Nspan,1:nEL,1:3))
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iND,i,j,k,x,y,z,s,rx,ry,rz,z1,detx,dety,detz,ix,jy,kz)
    do iND=1,nND
        call minloc_fast(xyzful(iND,1), x0, i0, invdh, i, detx)
        call minloc_fast(xyzful(iND,2), y0, j0, invdh, j, dety)
        do x=-1,2
            rx(x)=Phi(dble(x)-detx)
        enddo
        do y=-1,2
            ry(y)=Phi(dble(y)-dety)
        enddo
        call trimedindex(i, xDim, ix, boundaryConditions(1:2))
        call trimedindex(j, yDim, jy, boundaryConditions(3:4))
        do s=1,Nspan+1
            z1 = xyzful(iND,3)+dspan * (s - 1)
            call minloc_fast(z1, z0, k0, invdh, k, detz)
            do z=-1,2
                rz(z)=Phi(dble(z)-detz)
            enddo
            call trimedindex(k, zDim, kz, boundaryConditions(5:6))
            ! interpolate fluid velocity to body nodes
            velfulIB(s,iND,1:3)=0.0d0
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        velfulIB(s,iND,1:3)=velfulIB(s,iND,1:3)+uuu(kz(z),jy(y),ix(x),1:3)*rx(x)*ry(y)*rz(z)
                    enddo
                enddo
            enddo
        enddo
        xyzfulIB(:,iND,1:3)=xyzfulIB(:,iND,1:3)+velfulIB(:,iND,1:3)*dt
    enddo
    !$OMP END PARALLEL DO
endif

!==================================================================================================
!   compute displacement, velocity, area at surface element center
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,s,nt,x1,x2,y1,y2,ax,ay)
do  iEL=1,nEL
    i=ele(iEL,1)
    j=ele(iEL,2)
    nt=ele(iEL,4)

    x1=xyzful(i,1)
    x2=xyzful(j,1)
    y1=xyzful(i,2)
    y2=xyzful(j,2)
    if(nt/=2) write(*,*) 'only support line segments'
    do s=1,Nspan
        posElem(s,iEL,1)=(x1+x2)*0.5d0
        posElem(s,iEL,2)=(y1+y2)*0.5d0
        posElem(s,iEL,3)=xyzful(i,3) + dspan*(s-0.5)
        velElem(s,iEL,1:2)=(velful(i,1:2)+velful(j,1:2))*0.5d0
        velElem(s,iEL,3)=0.d0
        if(Palpha.gt.0.d0) then
            posElemIB(s,iEL,1:3)=(xyzfulIB(s,i,1:3)+xyzfulIB(s,j,1:3)+xyzfulIB(s+1,i,1:3)+xyzfulIB(s+1,j,1:3))*0.25d0
        endif
    enddo
    ax =(x1-x2)
    ay =(y1-y2)
    areaElem(iEL)=dsqrt( ax*ax + ay*ay) * dspan
enddo
!$OMP END PARALLEL DO

!**************************************************************************************************
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
do x = 1, xDim
    force(:,:,x,1)=0.0d0
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
do x = 1, xDim
    force(:,:,x,2)=0.0d0
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
do x = 1, xDim
    force(:,:,x,3)=0.0d0
enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL)
do iEL = 1, nEL
    forceElem(:,iEL,1:3)=0.0d0
enddo
!$OMP END PARALLEL DO

!***********************************************************************************************
dmaxLBM=1.0d0
iterLBM=0
do  while( iterLBM<ntolLBM .and. dmaxLBM>dtolLBM)
    !***********************************************************************************************
    ! compute the velocity of IB nodes at element center
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,k,x,y,z,s,rx,ry,rz,detx,dety,detz,ix,jy,kz)
    do  iEL=1,nEL
        call minloc_fast(posElem(1,iEL,1), x0, i0, invdh, i, detx)
        call minloc_fast(posElem(1,iEL,2), y0, j0, invdh, j, dety)
        do x=-1,2
            rx(x)=Phi(dble(x)-detx)
        enddo
        do y=-1,2
            ry(y)=Phi(dble(y)-dety)
        enddo
        call trimedindex(i, xDim, ix, boundaryConditions(1:2))
        call trimedindex(j, yDim, jy, boundaryConditions(3:4))
        do s=1,Nspan
            call minloc_fast(posElem(s,iEL,3), z0, k0, invdh, k, detz)
            do z=-1,2
                rz(z)=Phi(dble(z)-detz)
            enddo
            call trimedindex(k, zDim, kz, boundaryConditions(5:6))
            velElemIB(s,iEL,1:3)=0.0d0
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        velElemIB(s,iEL,1:3)=velElemIB(s,iEL,1:3)+uuu(kz(z),jy(y),ix(x),1:3)*rx(x)*ry(y)*rz(z)
                    enddo
                enddo
            enddo
        enddo
    enddo
    !$OMP END PARALLEL DO
    !***********************************************************************************************
    ! calculate interaction force
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,s)
    do  iEL=1,nEL
        do s=1,Nspan
            forceElemTemp(s,iEL,1:3) = -Pbeta* 2.0d0*denIn*(velElem(s,iEL,1:3)-velElemIB(s,iEL,1:3))/dt*areaElem(iEL)*dh
            if(Palpha.gt.0.d0) then
                forceElemTemp(s,iEL,1:3) = forceElemTemp(s,iEL,1:3) - Palpha*2.0d0*denIn*(posElem(s,iEL,1:3)-posElemIB(s,iEL,1:3))/dt*areaElem(iEL)*dh
            endif
            if ((.not. IEEE_IS_FINITE(forceElemTemp(s,iEL,1))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(s,iEL,2))) .or. (.not. IEEE_IS_FINITE(forceElemTemp(s,iEL,3)))) then
                write(*, *) 'Nan found in forceElemTemp', forceElemTemp
                stop
            endif
        enddo
    enddo
    !$OMP END PARALLEL DO
    !***********************************************************************************************
    ! calculate Eulerian body force
    ! no parallel to avoid write conflict to forceTemp
    do iEL=1,nEL
        call minloc_fast(posElem(1,iEL,1), x0, i0, invdh, i, detx)
        call minloc_fast(posElem(1,iEL,2), y0, j0, invdh, j, dety)
        do x=-1,2
            rx(x)=Phi(dble(x)-detx)*invdh
        enddo
        do y=-1,2
            ry(y)=Phi(dble(y)-dety)*invdh
        enddo
        call trimedindex(i, xDim, ix, boundaryConditions(1:2))
        call trimedindex(j, yDim, jy, boundaryConditions(3:4))
        do s=1,Nspan
            call minloc_fast(posElem(s,iEL,3), z0, k0, invdh, k, detz)
            do z=-1,2
                rz(z)=Phi(dble(z)-detz)*invdh
            enddo
            call trimedindex(k, zDim, kz, boundaryConditions(5:6))
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        forceTemp(1:3) = -forceElemTemp(s,iEL,1:3)*rx(x)*ry(y)*rz(z)
                        ! update velocity
                        uuu(kz(z),jy(y),ix(x),1:3)  = uuu(kz(z),jy(y),ix(x),1:3)+0.5d0*dt*forceTemp(1:3)/den(kz(z),jy(y),ix(x))
                        force(kz(z),jy(y),ix(x),1:3) = force(kz(z),jy(y),ix(x),1:3) + forceTemp(1:3)
                    enddo
                enddo
            enddo
        enddo
    enddo
    forceElem(:,1:nEL,1:3) = forceElem(:,1:nEL,1:3)+forceElemTemp(:,1:nEL,1:3)
    ! convergence test
    dsum=Uref*nEL*Nspan
    dmaxLBM=0.0d0
    do iEL=1,nEL
        do s=1,Nspan
            dmaxLBM=dmaxLBM+dsqrt(sum((velElem(s,iEL,1:3)-velElemIB(s,iEL,1:3))**2))
        enddo
    enddo
    dmaxLBM=dmaxLBM/dsum
    iterLBM=iterLBM+1
!***********************************************************************************************
enddo
!write(*,'(A,I5,A,D20.10)')' iterLBM =',iterLBM,'    dmaxLBM =',dmaxLBM
!**************************************************************************************************
if(.false. .and. boundaryConditions(5).eq.symmetric) then
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
    do x=1,xDim
        force(1,:,x,1) = force(1,:,x,1) * 2.d0
        force(1,:,x,2) = force(1,:,x,2) * 2.d0
        force(1,:,x,3) = 0.d0
    enddo
    !$OMP END PARALLEL DO
endif
if(.false. .and. boundaryConditions(6).eq.symmetric) then
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
    do x=1,xDim
        force(zDim,:,x,1) = force(zDim,:,x,1) * 2.d0
        force(zDim,:,x,2) = force(zDim,:,x,2) * 2.d0
        force(zDim,:,x,3) = 0.d0
    enddo
    !$OMP END PARALLEL DO
endif
!**************************************************************************************************
!   element force to nodal force
extful(1:nND,1:6)=0.0d0
do iEL=1,nEL
    i=ele(iEL,1)
    j=ele(iEL,2)
    do s=1,Nspan
        extful(i,1:2)=extful(i,1:2)+forceElem(s,iEl,1:2)*0.5d0
        extful(j,1:2)=extful(j,1:2)+forceElem(s,iEl,1:2)*0.5d0
    enddo
enddo
if(Palpha.gt.0.d0) then
    deallocate(posElemIB)
endif
END SUBROUTINE calculate_interaction_force_quad

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    one dimensional delta function
!    copyright@ RuNanHua
!    xbgn    xxx     xend
!     -1      0   1   2
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

    subroutine get_phi_r(r,phi_r)
        implicit none
        real(8):: r(1:3),phi_r(1:3)
        real(8):: rr
        integer:: i

        do i=1,3
           rr=dabs(r(i))
           if(rr<1.0d0)then
                Phi_r(i)=(3.0d0-2.0d0*rr+dsqrt(1.0d0+4.0d0*rr-4.0d0*rr*rr))/8.0d0
           elseif(rr<2.0d0)then
                Phi_r(i)=(5.0d0-2.0d0*rr-dsqrt(-7.0d0+12.0d0*rr-4.0d0*rr*rr))/8.0d0
           else
                Phi_r(i)=0.0d0
           endif
        enddo

    end subroutine get_phi_r

    subroutine get_smoother_phi_r(r,phi_r)
        implicit none
        real(8):: r(1:3),phi_r(1:3)
        real(8):: rr , pi
        integer:: i
        pi=4.0d0*datan(1.0d0)

        do i=1,3
           rr=dabs(r(i))
           if(rr.le.0.5d0)then
                phi_r(i)=(12.0d0+pi-rr*rr*8.0d0)/32.0d0
           elseif(rr.le.1.5d0)then
                phi_r(i)=(2.0d0+(1-rr)*dsqrt(-2.0d0+8.0d0*rr-4.0d0*rr*rr)-dasin(dsqrt(2.0d0)*(rr-1)))/8.0d0
            elseif(rr.le.2.5d0)then
                phi_r(i)=(68.0d0-pi-48*rr+8*rr*rr+4*(rr-2)*dsqrt(-14.0d0+16.0d0*rr-4.0d0*rr*rr)+4*dasin(dsqrt(2.0d0)*(rr-2)))/64.0d0
           else
                phi_r(i)=0.0d0
           endif
        enddo

    end subroutine get_smoother_phi_r

subroutine initializexyzIB
USE simParam
USE ImmersedBoundary
implicit none
integer:: s,iFish,Nspanpts
if(Palpha.gt.0.d0) then
    Nspanpts = Nspan + 1
    allocate(xyzfulIB_all(1:Nspanpts,nND_all,6),xyzfulIB(1:Nspanpts,nND_max,6,1:nFish))
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(s,iFish)
    do iFish=1,nFish
        do s=1,Nspanpts
            xyzfulIB(s,1:nND(iFish),1:6,iFish)=xyzful(1:nND(iFish),1:6,iFish)
            xyzfulIB(s,1:nND(iFish),3,iFish)=xyzful(1:nND(iFish),3,iFish)+(s-1)*dspan
        enddo
    enddo
    !$OMP END PARALLEL DO
    call packxyzIB
endif
endsubroutine initializexyzIB

subroutine packxyzIB
USE simParam
USE ImmersedBoundary
implicit none
integer:: iND,iFish,icount
if(Palpha.gt.0.d0) then
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iND,icount,iFish)
    do iFish=1,nFish
        if(iFish.eq.1)then
            icount = 0
        elseif(iFish.ge.2)then
            icount = sum(nND(1:iFish-1))
        endif
        do iND=1,nND(iFish)
            xyzfulIB_all(:,iND+icount,1:6) = xyzfulIB(:,iND,1:6,iFish)
        enddo
    enddo
    !$OMP END PARALLEL DO
endif
endsubroutine packxyzIB

subroutine unpackxyzIB
USE simParam
USE ImmersedBoundary
implicit none
integer:: iND,iFish,icount
if(Palpha.gt.0.d0) then
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iND,icount,iFish)
    do iFish=1,nFish
        if(iFish.eq.1)then
            icount = 0
        else
            icount = sum(nND(1:iFish-1))
        endif
        do iND=1,nND(iFish)
            xyzfulIB(:,iND,1:6,iFish) =  xyzfulIB_all(:,iND + icount,1:6)
        enddo
    enddo
    !$OMP END PARALLEL DO
endif
endsubroutine unpackxyzIB
