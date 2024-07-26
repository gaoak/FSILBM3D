!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Calculate the repulsive force between solids
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptForceR(nFish,dxmin,dymin,dzmin,nND,nND_max,nEL,nEL_max,ele,xyzful,repful)
    USE ImmersedBoundary
    implicit none
    integer:: nFish,nEL_max,nND_max
    real(8):: dxmin,dymin,dzmin
    integer:: nND(1:nFish),nEL(1:nFish),ele(1:nFish,1:nEL_max,1:5)
    real(8):: xyzful(1:nFish,1:nND_max,1:6),repful(1:nFish,1:nND_max,1:6)
    !local
    integer:: iND,jND,iFish,jFish
    real(8):: delta_h,Phi,r(1:3),ds(1:3),phi_r(1:3),Lspan
    real(8):: minx,miny,maxx,maxy,minz,maxz
    real(8):: xmin(1:nFish),xmax(1:nFish),ymin(1:nFish),ymax(1:nFish),zmin(1:nFish),zmax(1:nFish)

    repful(:,:,:)=0.d0
    ds(1)=dxmin
    ds(2)=dymin
    ds(3)=dzmin
    Lspan=dspan*Nspan
    if(Nspan.eq.0)then
        Lspan=1.0d0
    endif
    
    do iFish=1,nFish
        xmin(iFish) = minval(xyzful(iFish,1:nND(iFish),1))-dxmin*1.5d0
        xmax(iFish) = maxval(xyzful(iFish,1:nND(iFish),1))+dxmin*1.5d0
        ymin(iFish) = minval(xyzful(iFish,1:nND(iFish),2))-dymin*1.5d0
        ymax(iFish) = maxval(xyzful(iFish,1:nND(iFish),2))+dymin*1.5d0 
    enddo

    do iFish=1,nFish
        do jFish=iFish+1,nFish
            minx = max(xmin(iFish),xmin(jFish))
            miny = max(ymin(iFish),ymin(jFish))
            maxx = min(xmax(iFish),xmax(jFish))
            maxy = min(ymax(iFish),ymax(jFish))
            if((minx > maxx).or.(miny > maxy)) then
                cycle
            endif
            ! overlapping regin [minx, maxx] X [miny, maxy] X [minz, maxz]
            do iND=1,nND(iFish)
                if ( (minx>xyzful(iFish,iND,1)) .or. (xyzful(iFish,iND,1)>maxx)  &
                .or. (miny>xyzful(iFish,iND,2)) .or. (xyzful(iFish,iND,2)>maxy)) then
                    cycle !point iND not in the overpalling region
                endif
                do jND=1,nND(jFish)
                    if ( (minx>xyzful(jFish,jND,1)) .or. (xyzful(jFish,jND,1)>maxx)  &
                    .or. (miny>xyzful(jFish,jND,2)) .or. (xyzful(jFish,jND,2)>maxy)) then
                        cycle !point jND not in the overpalling region
                    endif
                    r(1)=(xyzful(iFish,iND,1)-xyzful(jFish,jND,1))/dxmin
                    r(2)=(xyzful(iFish,iND,2)-xyzful(jFish,jND,2))/dymin
                    call get_phi_r(r,phi_r)
                    delta_h=phi_r(1)*phi_r(2)/dxmin/dymin/dsqrt(r(1)*r(1)+r(2)*r(2))
                    repful(iFish,iND,1:2)=repful(iFish,iND,1:2) + delta_h*r(1:2)*ds(1:2)*Lspan ! force
                    repful(jFish,jND,1:2)=repful(jFish,jND,1:2) - delta_h*r(1:2)*ds(1:2)*Lspan ! reaction force
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
SUBROUTINE calculate_interaction_force(zDim,yDim,xDim,nEL,nND,ele,dx,dy,dz,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
                                xyzful,velful,xyzfulIB,Palpha,Pbeta,ntolLBM,dtolLBM,force,extful,isUniformGrid)    
USE, INTRINSIC :: IEEE_ARITHMETIC
IMPLICIT NONE
integer,intent(in):: zDim,yDim,xDim,nEL,nND,ele(nEL,5),ntolLBM
real(8),intent(in):: dz(zDim),dy(yDim),dx(xDim),dh,Uref,denIn,dtolLBM,dt,Palpha,Pbeta
real(8),intent(in):: den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
logical,intent(in):: isUniformGrid(1:3)
real(8),intent(in):: xyzful(nND,6),velful(nND,6)
real(8),intent(inout)::xyzfulIB(nND,6),uuu(zDim,yDim,xDim,1:3)
real(8),intent(out)::extful(nND,6),force(zDim,yDim,xDim,1:3)
!==================================================================================================
integer:: i,j,k,x,y,z,xbgn,ybgn,zbgn,xend,yend,zend,iEL,nt,iterLBM,iND
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
        velfulIB(iND,1:3)=0.0
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
        ax =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.0
        ay =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.0
        az =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.0
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az)
    else
            write(*,*)'cell type is not defined'
    endif
enddo

!**************************************************************************************************
!**************************************************************************************************
forceElem(1:nEL,1:3)=0.0d0
dmaxLBM=1.0
iterLBM=0
!   ***********************************************************************************************
do  while( iterLBM<ntolLBM .and. dmaxLBM>dtolLBM)  
!   ***********************************************************************************************
!   compute the velocity of IB nodes at element center    
    do  iEL=1,nEL
        call my_minloc(posElem(iEL,1), xGrid, xDim, isUniformGrid(1), i)
        call my_minloc(posElem(iEL,2), yGrid, yDim, isUniformGrid(2), j)
        call my_minloc(posElem(iEL,3), zGrid, zDim, isUniformGrid(3), k)
        velElemIB(iEL,1:3)=0.0
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
            forceElemTemp(iEL,1:3) = -Pbeta* 2.0*denIn*(velElem(iEL,1:3)-velElemIB(iEL,1:3))/dt*areaElem(iEL)*dh
            if(Palpha.gt.0.d0) then
                forceElemTemp(iEL,1:3) = forceElemTemp(iEL,1:3) - Palpha*2.0*denIn*(posElem(iEL,1:3)-posElemIB(iEL,1:3))/dt*areaElem(iEL)*dh
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
                    uuu(z,y,x,1:3)  = uuu(z,y,x,1:3)+0.5*dt*forceTemp(1:3)/den(z,y,x)
                    force(z,y,x,1:3)=force(z,y,x,1:3) + forcetemp(1:3)
                enddo
            enddo
        enddo
    enddo
    forceElem(1:nEL,1:3) = forceElem(1:nEL,1:3)+forceElemTemp(1:nEL,1:3)
!   convergence test
    if(iterLBM==0)then
        dsum=0.0
        do iEL=1,nEL
        dsum=dsum+dsqrt(sum((velElem(iEL,1:3)-velElemIB(iEL,1:3))**2))
        enddo
    endif
    dsum=Uref*nEL
    
    dmaxLBM=0.0
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
forceNode(1:nND,1:3)=0.0
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
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_interaction_force_quad(zDim,yDim,xDim,nEL,nND,ele,dx,dy,dz,dh,Uref,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
    xyzful,velful,xyzfulIB,Palpha,Pbeta,ntolLBM,dtolLBM,force,extful,isUniformGrid,Nspan,dspan)
USE, INTRINSIC :: IEEE_ARITHMETIC
IMPLICIT NONE
integer,intent(in):: zDim,yDim,xDim,nEL,nND,ele(nEL,5),ntolLBM,Nspan
real(8),intent(in):: dz(zDim),dy(yDim),dx(xDim),dh,Uref,denIn,dtolLBM,dt,Palpha,Pbeta,dspan
real(8),intent(in):: den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
logical,intent(in):: isUniformGrid(1:3)
real(8),intent(in):: xyzful(nND,6),velful(nND,6)
real(8),intent(inout)::xyzfulIB(1:Nspan+1,nND,6),uuu(zDim,yDim,xDim,1:3)
real(8),intent(out)::extful(nND,6),force(zDim,yDim,xDim,1:3)
!==================================================================================================
integer:: i,j,k,x,y,z,s,iEL,nt,iterLBM,iND,ix,iy,iz
real(8):: rx(-1:2),ry(-1:2),rz(-1:2),Phi,dmaxLBM,dsum,invdh,forceTemp(1:3)
real(8):: x1,x2,y1,y2,z1,ax,ay
real(8):: forceNode(nND,3),velfulIB(1:Nspan+1,nND,3)
real(8):: forceElem(1:Nspan,nEL,3),forceElemTemp(1:Nspan,nEL,3),areaElem(nEL)
real(8):: posElem(1:Nspan,nEL,3),velElem(1:Nspan,nEL,3),velElemIB(1:Nspan,nEL,3)
real(8),allocatable::posElemIB(:,:,:)
!==================================================================================================
!   compute velocity and displacement at IB nodes
invdh = 1.D0/dh
if(Palpha.gt.0.d0) then
    allocate(posElemIB(1:Nspan,1:nEL,1:3))
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iND,i,j,k,x,y,z,s,rx,ry,rz,z1)
    do iND=1,nND
        call my_minloc(xyzful(iND,1), xGrid, xDim, isUniformGrid(1), i)
        call my_minloc(xyzful(iND,2), yGrid, yDim, isUniformGrid(2), j)
        do x=-1+i,2+i
            rx(x-i)=Phi((xyzful(iND,1)-xGrid(x))*invdh)
        enddo
        do y=-1+j,2+j
            ry(y-j)=Phi((xyzful(iND,2)-yGrid(y))*invdh)
        enddo
        do s=1,Nspan+1
            z1 = xyzful(iND,3)+dspan * (s - 1)
            call my_minloc(z1, zGrid, zDim, isUniformGrid(3), k)
            do z=-1+k,2+k
                rz(z-k)=Phi((z1-zGrid(z))*invdh)
            enddo
            ! interpolate fluid velocity to body nodes
            velfulIB(s,iND,1:3)=0.0
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        velfulIB(s,iND,1:3)=velfulIB(s,iND,1:3)+uuu(z+k,y+j,x+i,1:3)*rx(x)*ry(y)*rz(z)
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
        posElem(s,iEL,3)=xyzful(i,3) + dspan* (s-0.5)
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
dmaxLBM=1.0
iterLBM=0
do  while( iterLBM<ntolLBM .and. dmaxLBM>dtolLBM)  
    !***********************************************************************************************
    ! compute the velocity of IB nodes at element center
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,k,x,y,z,rx,ry,rz)
    do  iEL=1,nEL
        call my_minloc(posElem(1,iEL,1), xGrid, xDim, isUniformGrid(1), i)
        call my_minloc(posElem(1,iEL,2), yGrid, yDim, isUniformGrid(2), j)
        if(i.lt.2 .or. i.gt.xDim-2 .or. j.lt.2 .or. j.gt.yDim-2) then
            write(*, *) 'Point too close to boundary posElem(1,iEL,1:2): ', posElem(1,iEL,1:2), iEL
            stop
        endif
        do x=-1+i,2+i
            rx(x-i)=Phi((posElem(1,iEL,1)-xGrid(x))*invdh)
        enddo
        do y=-1+j,2+j
            ry(y-j)=Phi((posElem(1,iEL,2)-yGrid(y))*invdh)
        enddo
        do s=1,Nspan
            call my_minloc(posElem(s,iEL,3), zGrid, zDim, isUniformGrid(3), k)
            if(k.lt.2 .or. k.gt.zDim-2) then
                write(*, *) 'Point too close to boundary posElem(s,iEL,3): ', posElem(s,iEL,3), s, iEL
                stop
            endif
            do z=-1+k,2+k
                rz(z-k)=Phi((posElem(s,iEL,3)-zGrid(z))*invdh)
            enddo
            velElemIB(s,iEL,1:3)=0.0
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        velElemIB(s,iEL,1:3)=velElemIB(s,iEL,1:3)+uuu(z+k,y+j,x+i,1:3)*rx(x)*ry(y)*rz(z)
                    enddo
                enddo
            enddo
        enddo
    enddo
    !***********************************************************************************************
    ! calculate interaction force
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,s)
    do  iEL=1,nEL
        do s=1,Nspan
            forceElemTemp(s,iEL,1:3) = -Pbeta* 2.0*denIn*(velElem(s,iEL,1:3)-velElemIB(s,iEL,1:3))/dt*areaElem(iEL)*dh
            if(Palpha.gt.0.d0) then
                forceElemTemp(s,iEL,1:3) = forceElemTemp(s,iEL,1:3) - Palpha*2.0*denIn*(posElem(s,iEL,1:3)-posElemIB(s,iEL,1:3))/dt*areaElem(iEL)*dh
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
        call my_minloc(posElem(1,iEL,1), xGrid, xDim, isUniformGrid(1), i)
        call my_minloc(posElem(1,iEL,2), yGrid, yDim, isUniformGrid(2), j)
        do x=-1+i,2+i
            rx(x-i)=Phi((posElem(1,iEL,1)-xGrid(x))*invdh)
        enddo
        do y=-1+j,2+j
            ry(y-j)=Phi((posElem(1,iEL,2)-yGrid(y))*invdh)
        enddo
        do s=1,Nspan
            call my_minloc(posElem(s,iEL,3), zGrid, zDim, isUniformGrid(3), k)
            do z=-1+k,2+k
                rz(z-k)=Phi((posElem(s,iEL,3)-zGrid(z))*invdh)
            enddo
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        forceTemp(1:3) = -forceElemTemp(s,iEL,1:3)*rx(x)*ry(y)*rz(z)*invdh*invdh*invdh
                        ! update velocity
                        ix = i + x
                        iy = j + y
                        iz = k + z
                        uuu(iz,iy,ix,1:3)  = uuu(iz,iy,ix,1:3)+0.5*dt*forceTemp(1:3)/den(iz,iy,ix)
                        force(iz,iy,ix,1:3) = force(iz,iy,ix,1) + forceTemp(1:3)
                    enddo
                enddo
            enddo
        enddo
    enddo
    forceElem(:,1:nEL,1:3) = forceElem(:,1:nEL,1:3)+forceElemTemp(:,1:nEL,1:3)
    ! convergence test
    dsum=Uref*nEL
    dmaxLBM=0.0
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
!**************************************************************************************************
!   element force to nodal force
forceNode(1:nND,1:3)=0.0
do iEL=1,nEL
    i=ele(iEL,1)
    j=ele(iEL,2)
    do s=1,Nspan
        forceNode(i,1:3)=forceNode(i,1:3)+forceElem(s,iEl,1:3)*0.5d0
        forceNode(j,1:3)=forceNode(j,1:3)+forceElem(s,iEl,1:3)*0.5d0
    enddo
enddo
extful(1:nND,1:2) = forceNode(1:nND,1:2)
extful(1:nND,3:6) = 0.0d0
if(Palpha.gt.0.d0) then
    deallocate(posElemIB)
endif
END SUBROUTINE calculate_interaction_force_quad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    �����������Ӧ������
!   ��������������Ӧ��
!    copyright@ RuNanHua 
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptStrs(zDim,yDim,xDim,nEL,nND,ele,dh,dx,dy,dz,mu,rr,u,p,xGrid,yGrid,zGrid,xyzful,extful)    
    IMPLICIT NONE
    integer,intent(in):: zDim,yDim,xDim,nEL,nND,ele(nEL,5)
    real(8):: dh,dx(xDim),dy(yDim),dz(zDim),mu,rr
    real(8):: u(zDim,yDim,xDim,1:3),p(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
    real(8):: xyzful(nND,6),extful(nND,6)
!   ===============================================================================================
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,cx,cy,cz,px(2),py(2),pz(2)
    real(8):: nx(2),ny(2),nz(2),ax(2),ay(2),az(2),area,ll1,ll2,w1,w2
    real(8):: Sxx1,Sxy1,Sxz1,Syy1,Syz1,Szz1,Pre1,Sxx2,Sxy2,Sxz2,Syy2,Syz2,Szz2,Pre2
    real(8):: SxxTemp,SxyTemp,SxzTemp,SyyTemp,SyzTemp,SzzTemp
    integer:: i,j,k,m,nELt,iEL
    integer:: x,y,z,xbgn,ybgn,zbgn,xend,yend,zend
    real(8):: weightm,forceElem(nEL,3),forceElemTemp(nEL,3),forceNode(nND,3)

!   ��������λ��ƽ��
!
    ll1=0.5*rr*dh
    ll2=1.0*rr*dh

    w1=ll2/(ll2-ll1)
    w2=1.0d0-w1
!   ***********************************************************************************************
!   ���㵥Ԫ��
    forceElem(1:nEL,1:3)=0.0d0
    do    iEL=1,nEL
        
        i  = ele(iEL,1)
        j  = ele(iEL,2)
        k  = ele(iEL,3)
        nELt= ele(iEL,4)
!
        if    (nELt == 2) then
!            frame
            forceElem(iEL,1:3)=[0.0,0.0,0.0]
        elseif (nELt == 3) then
!           plate
            x1=xyzful(i,1)
            x2=xyzful(j,1)
            x3=xyzful(k,1)
            y1=xyzful(i,2)
            y2=xyzful(j,2)
            y3=xyzful(k,2)
            z1=xyzful(i,3)
            z2=xyzful(j,3)
            z3=xyzful(k,3)

            cx=(x1+x2+x3)/3.0
            cy=(y1+y2+y3)/3.0
            cz=(z1+z2+z3)/3.0

!            determine vector area         
            ax(1) =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.0
            ay(1) =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.0
            az(1) =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.0
            area=dsqrt( ax(1)*ax(1) + ay(1)*ay(1) + az(1)*az(1))
            nx(1)=ax(1)/area
            ny(1)=ay(1)/area
            nz(1)=az(1)/area

            
            ax(2) =-ax(1)
            ay(2) =-ay(1)
            az(2) =-az(1)

            nx(2)=-nx(1)
            ny(2)=-ny(1)
            nz(2)=-nz(1)

!==================================================================================================
!   interpolation of flow variables by inverse-distance-weighted summation           
            forceElem(iEL,1:3)=0.0d0
            do  m=1,2
                !******************************************
                !******************************************
                px(m)=cx+nx(m)*ll1
                py(m)=cy+ny(m)*ll1
                pz(m)=cz+nz(m)*ll1
                
                x=minloc(dabs(px(m)-xGrid(1:xDim)),1)
                if(px(m)-xGrid(x)>0.0d0)then
                xbgn=x; xend=x+1
                else
                xbgn=x-1;xend=x
                endif              
                y=minloc(dabs(py(m)-yGrid(1:yDim)),1)
                if(py(m)-yGrid(y)>0.0d0)then
                ybgn=y; yend=y+1
                else
                ybgn=y-1;yend=y
                endif                                
                z=minloc(dabs(pz(m)-zGrid(1:zDim)),1)
                if(pz(m)-zGrid(z)>0.0d0)then
                zbgn=z; zend=z+1
                else
                zbgn=z-1;zend=z
                endif

                Sxx1=0.d0
                Syy1=0.d0
                Szz1=0.d0
                Sxy1=0.d0
                Sxz1=0.d0
                Syz1=0.d0
                Pre1=0.d0

                do    x=xbgn,xend
                do    y=ybgn,yend
                do    z=zbgn,zend

                weightm=(dx(x)-dabs(xGrid(x)-px(m)))*(dy(y)-dabs(yGrid(y)-py(m)))*(dz(z)-dabs(zGrid(z)-pz(m)))
!               stress on the interpolation supported node
                SxxTemp=Mu*(u(z,y,x+1,1)-u(z,y,x-1,1))/dx(x)      !  x- normal stress             
                SyyTemp=Mu*(u(z,y+1,x,2)-u(z,y-1,x,2))/dy(y)      !  y- normal stress
                SzzTemp=Mu*(u(z+1,y,x,3)-u(z-1,y,x,3))/dz(z)      !  z- normal stress

                SxyTemp=0.5d0*Mu*(  (u(z,y+1,x,1)-u(z,y-1,x,1))/dy(y)  + (u(z,y,x+1,2)-u(z,y,x-1,2))/dx(x)   )  ! cross stress
                SxzTemp=0.5d0*Mu*(  (u(z+1,y,x,1)-u(z-1,y,x,1))/dz(z)  + (u(z,y,x+1,3)-u(z,y,x-1,3))/dx(x)   )  ! cross stress
                SyzTemp=0.5d0*Mu*(  (u(z+1,y,x,2)-u(z-1,y,x,2))/dz(z)  + (u(z,y+1,x,3)-u(z,y-1,x,3))/dy(y)   )  ! cross stress

                Sxx1=Sxx1+SxxTemp*weightm/(dx(x)*dy(y)*dz(z))
                Syy1=Syy1+SyyTemp*weightm/(dx(x)*dy(y)*dz(z))
                Szz1=Szz1+SzzTemp*weightm/(dx(x)*dy(y)*dz(z))
                Sxy1=Sxy1+SxyTemp*weightm/(dx(x)*dy(y)*dz(z))
                Sxz1=Sxz1+SxzTemp*weightm/(dx(x)*dy(y)*dz(z))
                Syz1=Syz1+SyzTemp*weightm/(dx(x)*dy(y)*dz(z))
                Pre1=Pre1+p(z,y,x)*weightm/(dx(x)*dy(y)*dz(z))

                enddo
                enddo
                enddo
                !******************************************
                !******************************************
                px(m)=cx+nx(m)*ll2
                py(m)=cy+ny(m)*ll2
                pz(m)=cz+nz(m)*ll2

                x=minloc(dabs(px(m)-xGrid(1:xDim)),1)
                if(px(m)-xGrid(x)>0.0d0)then
                xbgn=x; xend=x+1
                else
                xbgn=x-1;xend=x
                endif              
                y=minloc(dabs(py(m)-yGrid(1:yDim)),1)
                if(py(m)-yGrid(y)>0.0d0)then
                ybgn=y; yend=y+1
                else
                ybgn=y-1;yend=y
                endif                                
                z=minloc(dabs(pz(m)-zGrid(1:zDim)),1)
                if(pz(m)-zGrid(z)>0.0d0)then
                zbgn=z; zend=z+1
                else
                zbgn=z-1;zend=z
                endif

                Sxx2=0.d0
                Syy2=0.d0
                Szz2=0.d0
                Sxy2=0.d0
                Sxz2=0.d0
                Syz2=0.d0
                Pre2=0.d0

                do    x=xbgn,xend
                do    y=ybgn,yend
                do    z=zbgn,zend

                weightm=(dx(x)-dabs(xGrid(x)-px(m)))*(dy(y)-dabs(yGrid(y)-py(m)))*(dz(z)-dabs(zGrid(z)-pz(m)))
!               stress on the interpolation supported node
                SxxTemp=Mu*(u(z,y,x+1,1)-u(z,y,x-1,1))/dx(x)      !  x- normal stress             
                SyyTemp=Mu*(u(z,y+1,x,2)-u(z,y-1,x,2))/dy(y)      !  y- normal stress
                SzzTemp=Mu*(u(z+1,y,x,3)-u(z-1,y,x,3))/dz(z)      !  z- normal stress

                SxyTemp=0.5d0*Mu*(  (u(z,y+1,x,1)-u(z,y-1,x,1))/dy(y)  + (u(z,y,x+1,2)-u(z,y,x-1,2))/dx(x)   )  ! cross stress
                SxzTemp=0.5d0*Mu*(  (u(z+1,y,x,1)-u(z-1,y,x,1))/dz(z)  + (u(z,y,x+1,3)-u(z,y,x-1,3))/dx(x)   )  ! cross stress
                SyzTemp=0.5d0*Mu*(  (u(z+1,y,x,2)-u(z-1,y,x,2))/dz(z)  + (u(z,y+1,x,3)-u(z,y-1,x,3))/dy(y)   )  ! cross stress

                Sxx2=Sxx2+SxxTemp*weightm/(dx(x)*dy(y)*dz(z))
                Syy2=Syy2+SyyTemp*weightm/(dx(x)*dy(y)*dz(z))
                Szz2=Szz2+SzzTemp*weightm/(dx(x)*dy(y)*dz(z))
                Sxy2=Sxy2+SxyTemp*weightm/(dx(x)*dy(y)*dz(z))
                Sxz2=Sxz2+SxzTemp*weightm/(dx(x)*dy(y)*dz(z))
                Syz2=Syz2+SyzTemp*weightm/(dx(x)*dy(y)*dz(z))
                Pre2=Pre2+p(z,y,x)*weightm/(dx(x)*dy(y)*dz(z))

                enddo
                enddo
                enddo
                !******************************************
                !******************************************

                forceElem(iEL,1)=forceElem(iEL,1)-Pre2*ax(m)+(w1*Sxx1+w2*Sxx2)*ax(m)+(w1*Sxy1+w2*Sxy2)*ay(m)+(w1*Sxz1+w2*Sxz2)*az(m)
                forceElem(iEL,2)=forceElem(iEL,2)-Pre2*ay(m)+(w1*Sxy1+w2*Sxy2)*ax(m)+(w1*Syy1+w2*Syy2)*ay(m)+(w1*Syz1+w2*Syz2)*az(m)
                forceElem(iEL,3)=forceElem(iEL,3)-Pre2*az(m)+(w1*Sxz1+w2*Sxz2)*ax(m)+(w1*Syz1+w2*Syz2)*ay(m)+(w1*Szz1+w2*Szz2)*az(m)
            enddo
            
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif       
    enddo
!   ***********************************************************************************************
!   element force to nodal force
    forceNode(1:nND,1:3)=0.0d0
    do    iEL=1,nEL       
        i  = ele(iEL,1)
        j  = ele(iEL,2)
        k  = ele(iEL,3)
        nELt= ele(iEL,4)
        forceNode(i,1:3)=forceNode(i,1:3)+forceElem(iEL,1:3)/3.0d0
        forceNode(j,1:3)=forceNode(j,1:3)+forceElem(iEL,1:3)/3.0d0
        forceNode(k,1:3)=forceNode(k,1:3)+forceElem(iEL,1:3)/3.0d0
    enddo

    extful(1:nND,1:3) = forceNode(1:nND,1:3)
    extful(1:nND,4:6) = 0.0d0
    END SUBROUTINE



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

subroutine initializexyzIB
USE simParam
USE ImmersedBoundary
implicit none
integer:: iND,s,iFish,Nspanpts
if(Palpha.gt.0.d0) then
    Nspanpts = Nspan + 1
    allocate(xyzfulIB_all(1:Nspanpts,nND_all,6),xyzfulIB(1:Nspanpts,1:nFish,nND_max,6))
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iFish,s)
    do iFish=1,nFish
        do s=1,Nspanpts
            xyzfulIB(s,iFish,1:nND(iFish),1:6)=xyzful(iFish,1:nND(iFish),1:6)
            xyzfulIB(s,iFish,1:nND(iFish),3)=xyzful(iFish,1:nND(iFish),3)+(s-1)*dspan
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
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iFish,iND,icount)
    do iFish=1,nFish
        if(iFish.eq.1)then
            icount = 0
        elseif(iFish.ge.2)then
            icount = sum(nND(1:iFish-1))
        endif
        do iND=1,nND(iFish)
            xyzfulIB_all(:,iND+icount,1:6) = xyzfulIB(:,iFish,iND,1:6)
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
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iFish,iND,icount)
    do iFish=1,nFish
        if(iFish.eq.1)then
            icount = 0
        else
            icount = sum(nND(1:iFish-1))
        endif
        do iND=1,nND(iFish)
            xyzfulIB(:,iFish,iND,1:6) =  xyzfulIB_all(:,iND + icount,1:6)
        enddo
    enddo
    !$OMP END PARALLEL DO
endif
endsubroutine unpackxyzIB