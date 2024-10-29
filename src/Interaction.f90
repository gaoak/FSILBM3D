!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Calculate the repulsive force between solids
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptForceR(dxmin,dymin,dzmin,nND,nND_max,xyzful,repful,nFish,dspan,Nspan)
    implicit none
    integer:: nFish,nND_max,nND(1:nFish),Nspan(1:nFish)
    real(8):: dxmin,dymin,dzmin
    real(8):: xyzful(1:nND_max,1:6,1:nFish),repful(1:nND_max,1:6,1:nFish),dspan(1:nFish)
    !local
    integer:: iND,jND,iFish,jFish
    real(8):: delta_h,r(1:3),ds(1:3),phi_r(1:3),SpanLength
    real(8):: minx,miny,maxx,maxy,minz,maxz
    real(8):: xmin(1:nFish),xmax(1:nFish),ymin(1:nFish),ymax(1:nFish),zmin(1:nFish),zmax(1:nFish)

    repful(:,:,:)=0.d0
    ds(1)=dxmin
    ds(2)=dymin
    ds(3)=dzmin

    do iFish=1,nFish
        xmin(iFish) = minval(xyzful(1:nND(iFish),1,iFish))-dxmin*2.5d0
        xmax(iFish) = maxval(xyzful(1:nND(iFish),1,iFish))+dxmin*2.5d0
        ymin(iFish) = minval(xyzful(1:nND(iFish),2,iFish))-dymin*2.5d0
        ymax(iFish) = maxval(xyzful(1:nND(iFish),2,iFish))+dymin*2.5d0
        zmin(iFish) = minval(xyzful(1:nND(iFish),3,iFish))-dzmin*2.5d0
        zmax(iFish) = maxval(xyzful(1:nND(iFish),3,iFish))+dzmin*2.5d0
    enddo

    do iFish=1,nFish
        SpanLength=dspan(iFish)*Nspan(iFish)
        if(Nspan(iFish).eq.0)then
            SpanLength=1.0d0
        endif
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
SUBROUTINE calculate_interaction_force()
    USE simParam
    IMPLICIT NONE
    !==================================================================================================
    integer:: iFish
    integer:: i,j,k,iEL,nt,iterLBM
    real(8):: dmaxLBM,dsum
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az
    real(8):: forceElem(nEL_max,3,nFish),forceNode(nND_max,3,nFish),areaElemTemp(nEL_max,nFish)
    real(8):: posElem(nEL_max,3,nFish),velElem(nEL_max,3,nFish),velElemIB(nEL_max,3,nFish)
    !==================================================================================================
    !   compute displacement, velocity, area at surface element center
    do iFish=1,nFish
        do  iEL=1,nEL(iFish)
            i=ele(iEL,1,iFish)
            j=ele(iEL,2,iFish)
            k=ele(iEL,3,iFish)
            nt=ele(iEL,4,iFish)

            x1=xyzful(i,1,iFish)
            x2=xyzful(j,1,iFish)
            x3=xyzful(k,1,iFish)
            y1=xyzful(i,2,iFish)
            y2=xyzful(j,2,iFish)
            y3=xyzful(k,2,iFish)
            z1=xyzful(i,3,iFish)
            z2=xyzful(j,3,iFish)
            z3=xyzful(k,3,iFish)

            if(nt==2)then
                posElem(iEL,1:3,iFish)=(xyzful(i,1:3,iFish)+xyzful(j,1:3,iFish))/2.0d0
                velElem(iEL,1:3,iFish)=(velful(i,1:3,iFish)+velful(j,1:3,iFish))/2.0d0
                ax =(x1-x2)
                ay =(y1-y2)
                az =(z1-z2)
                areaElemTemp(iEL,iFish)=dsqrt( ax*ax + ay*ay + az*az)

            elseif(nt==3)then
                posElem(iEL,1:3,iFish)=(xyzful(i,1:3,iFish)+xyzful(j,1:3,iFish)+xyzful(k,1:3,iFish))/3.0d0
                velElem(iEL,1:3,iFish)=(velful(i,1:3,iFish)+velful(j,1:3,iFish)+velful(k,1:3,iFish))/3.0d0
                ax =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.0d0
                ay =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.0d0
                az =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.0d0
                areaElemTemp(iEL,iFish)=dsqrt( ax*ax + ay*ay + az*az)
            else
                    write(*,*)'cell type is not defined'
            endif
        enddo
    enddo

    !**************************************************************************************************
    !**************************************************************************************************
    forceElem(1:nEL_max,1:3,1:nFish)=0.0d0
    dmaxLBM=1.0d0
    iterLBM=0
    !   ***********************************************************************************************
    do  while( iterLBM<ntolLBM .and. dmaxLBM>dtolLBM)

        dmaxLBM=0.0d0
        dsum=0.0d0

        do iFish=1,nFish

            call calculate_interaction_force_core(zDim,yDim,xDim,nEL(iFish),ele(1:nEL(iFish),1:5,iFish),dh,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
            Pbeta,force,isUniformGrid,posElem(1:nEL(iFish),1:3,iFish),velElem(1:nEL(iFish),1:3,iFish), &
            areaElemTemp(1:nEL(iFish),iFish),forceElem(1:nEL(iFish),1:3,iFish),velElemIB(1:nEL(iFish),1:3,iFish))

            dsum=dsum+Uref*nEL(iFish)

            do iEL=1,nEL(iFish)
                dmaxLBM=dmaxLBM+dsqrt(sum((velElem(iEL,1:3,iFish)-velElemIB(iEL,1:3,iFish))**2))
            enddo
    !   ***********************************************************************************************
        enddo
        dmaxLBM=dmaxLBM/dsum
        iterLBM=iterLBM+1
    enddo
    !write(*,'(A,I5,A,D20.10)')' iterLBM =',iterLBM,'    dmaxLBM =',dmaxLBM
    !**************************************************************************************************
    !**************************************************************************************************
    !   element force to nodal force
    forceNode(1:nND_max,1:3,1:nFish)=0.0d0
    extful(1:nND_max,1:6,1:nFish)=0.0d0
    do iFish=1,nFish
        do    iEL=1,nEL(iFish)
            i=ele(iEL,1,iFish)
            j=ele(iEL,2,iFish)
            k=ele(iEL,3,iFish)
            nt=ele(iEL,4,iFish)
            forceNode(i,1:3,iFish)=forceNode(i,1:3,iFish)+forceElem(iEL,1:3,iFish)/3.0d0
            forceNode(j,1:3,iFish)=forceNode(j,1:3,iFish)+forceElem(iEL,1:3,iFish)/3.0d0
            forceNode(k,1:3,iFish)=forceNode(k,1:3,iFish)+forceElem(iEL,1:3,iFish)/3.0d0
        enddo
        extful(1:nND(iFish),1:3,iFish) = forceNode(1:nND(iFish),1:3,iFish)
        extful(1:nND(iFish),4:6,iFish) = 0.0d0
    enddo
END SUBROUTINE

SUBROUTINE calculate_interaction_force_core(zDim,yDim,xDim,nEL,ele,dh,denIn,dt,uuu,den,xGrid,yGrid,zGrid,  &
    Pbeta,force,isUniformGrid,posElem,velElem,areaElem,forceElem,velElemIB)
    USE, INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    integer,intent(in):: zDim,yDim,xDim,nEL,ele(nEL,5)
    real(8),intent(in):: dh,denIn,dt,Pbeta
    real(8),intent(in):: den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
    logical,intent(in):: isUniformGrid(1:3)
    real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
    real(8),intent(out)::force(zDim,yDim,xDim,1:3)
    !==================================================================================================
    real(8),intent(in):: posElem(nEL,3),velElem(nEL,3)
    real(8),intent(in):: areaElem(nEL)
    real(8),intent(inout)::forceElem(nEL,3)
    real(8),intent(out)::velElemIB(nEL,3)
    !==================================================================================================
    integer:: i,j,k,x,y,z,iEL
    real(8):: rx,ry,rz,Phi,invdh,forcetemp(1:3)
    real(8):: forceElemTemp(nEL,3)
    !==================================================================================================
    !   compute velocity and displacement at IB nodes
    invdh = 1.D0/dh
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
END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    (displacement, velocity) spring, penalty method
!    calculate force at element center, distribute force to four nodes
!    solid should be in uniform grid
!    period body should span full from 1 to xDim+1
!    symmetric body should not exceed domain
!    the first body point should be in the domain or on the boundary
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_interaction_force_quad()
    USE simParam
    USE ImmersedBoundary
    use BoundCondParams
    IMPLICIT NONE
    !==================================================================================================
    integer:: iFish
    integer:: i,j,x,s,iEL,nt,iterLBM
    real(8):: dmaxLBM,dsum
    real(8):: x1,x2,y1,y2,ax,ay
    real(8):: forceElem(1:maxval(Nspan),nEL_max,3,nFish),areaElemTemp(nEL_max,nFish)
    real(8):: posElem(1:maxval(Nspan),nEL_max,3,nFish),velElem(1:maxval(Nspan),nEL_max,3,nFish),velElemIB(1:maxval(Nspan),nEL_max,3,nFish)
    real(8)::x0(nFish),y0(nFish),z0(nFish)
    integer::i0(nFish),j0(nFish),k0(nFish)
    !==================================================================================================
    do iFish=1,nFish
        call my_minloc(xyzful(1,1,iFish), xGrid, xDim, isUniformGrid(1), i0(iFish))
        call my_minloc(xyzful(1,2,iFish), yGrid, yDim, isUniformGrid(2), j0(iFish))
        call my_minloc(xyzful(1,3,iFish), zGrid, zDim, isUniformGrid(3), k0(iFish))
        x0(iFish) = xGrid(i0(iFish))
        y0(iFish) = yGrid(j0(iFish))
        z0(iFish) = zGrid(k0(iFish))
        !==================================================================================================
        !   compute displacement, velocity, area at surface element center
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,s,nt,x1,x2,y1,y2,ax,ay)
        do  iEL=1,nEL(iFish)
            i=ele(iEL,1,iFish)
            j=ele(iEL,2,iFish)
            nt=ele(iEL,4,iFish)

            x1=xyzful(i,1,iFish)
            x2=xyzful(j,1,iFish)
            y1=xyzful(i,2,iFish)
            y2=xyzful(j,2,iFish)
            if(nt/=2) write(*,*) 'only support line segments'
            do s=1,Nspan(iFish)
                posElem(s,iEL,1,iFish)=(x1+x2)*0.5d0
                posElem(s,iEL,2,iFish)=(y1+y2)*0.5d0
                posElem(s,iEL,3,iFish)=xyzful(i,3,iFish) + dspan(iFish)*(s-0.5)
                velElem(s,iEL,1:2,iFish)=(velful(i,1:2,iFish)+velful(j,1:2,iFish))*0.5d0
                velElem(s,iEL,3,iFish)=0.d0
            enddo
            ax =(x1-x2)
            ay =(y1-y2)
            areaElemTemp(iEL,iFish)=dsqrt( ax*ax + ay*ay) * dspan(iFish)
        enddo
        !$OMP END PARALLEL DO

        !**************************************************************************************************
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL)
        do iEL = 1, nEL(iFish)
            forceElem(:,iEL,1:3,iFish)=0.0d0
        enddo
        !$OMP END PARALLEL DO
    enddo

    !***********************************************************************************************
    dmaxLBM=1.0d0
    iterLBM=0
    do  while( iterLBM<ntolLBM .and. dmaxLBM>dtolLBM)

        dmaxLBM=0.0d0
        dsum=0.0d0

        do iFish=1,nFish

            call calculate_interaction_force_quad_core(zDim,yDim,xDim,nEL(iFish),dh,denIn,dt,uuu,den,  &
            Pbeta,force,Nspan(iFish),theta(iFish),boundaryConditions,i0(iFish),j0(iFish),k0(iFish),x0(iFish),y0(iFish),z0(iFish), &
            posElem(1:Nspan(iFish),1:nEL(iFish),1:3,iFish),velElem(1:Nspan(iFish),1:nEL(iFish),1:3,iFish), &
            areaElemTemp(1:nEL(iFish),iFish),forceElem(1:Nspan(iFish),1:nEL(iFish),1:3,iFish),velElemIB(1:Nspan(iFish),1:nEL(iFish),1:3,iFish))

        ! convergence test
        dsum=dsum+Uref*nEL(iFish)*Nspan(iFish)
        
        do iEL=1,nEL(iFish)
            do s=1,Nspan(iFish)
                dmaxLBM=dmaxLBM+dsqrt(sum((velElem(s,iEL,1:3,iFish)-velElemIB(s,iEL,1:3,iFish))**2))
            enddo
        enddo
    !***********************************************************************************************
        enddo
        dmaxLBM=dmaxLBM/dsum
        iterLBM=iterLBM+1
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
    extful(1:nND_max,1:6,1:nFish)=0.0d0
    do iFish=1,nFish
        do iEL=1,nEL(iFish)
            i=ele(iEL,1,iFish)
            j=ele(iEL,2,iFish)
            do s=1,Nspan(iFish)
                extful(i,1:2,iFish)=extful(i,1:2,iFish)+forceElem(s,iEl,1:2,iFish)*0.5d0
                extful(j,1:2,iFish)=extful(j,1:2,iFish)+forceElem(s,iEl,1:2,iFish)*0.5d0
            enddo
        enddo
    enddo
END SUBROUTINE calculate_interaction_force_quad

SUBROUTINE calculate_interaction_force_quad_core(zDim,yDim,xDim,nEL,dh,denIn,dt,uuu,den,  &
    Pbeta,force,Nspan,theta,boundaryConditions,i0,j0,k0,x0,y0,z0,posElem,velElem,areaElem,forceElem,velElemIB)
    USE, INTRINSIC :: IEEE_ARITHMETIC
    IMPLICIT NONE
    integer,intent(in):: zDim,yDim,xDim,nEL,Nspan
    real(8),intent(in):: dh,denIn,dt,Pbeta,theta
    real(8),intent(in):: den(zDim,yDim,xDim)
    integer,intent(in):: boundaryConditions(1:6)
    real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
    real(8),intent(out)::force(zDim,yDim,xDim,1:3)
    !==================================================================================================
    integer,intent(in):: i0,j0,k0
    real(8),intent(in):: x0,y0,z0
    real(8),intent(in):: posElem(1:Nspan,nEL,3),velElem(1:Nspan,nEL,3)
    real(8),intent(in):: areaElem(nEL)
    real(8),intent(inout)::forceElem(1:Nspan,nEL,3)
    real(8),intent(out)::velElemIB(1:Nspan,nEL,3)
    !==================================================================================================
    integer:: i,j,k,x,y,z,s,iEL
    integer:: ix(-1:2),jy(-1:2),kz(-1:2)
    real(8):: rx(-1:2),ry(-1:2),rz(-1:2),Phi,invdh,forceTemp(1:3)
    real(8):: forceElemTemp(1:Nspan,nEL,3)
    real(8)::detx,dety,detz
    !==================================================================================================
    !   compute velocity and displacement at IB nodes
    invdh = 1.D0/dh
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