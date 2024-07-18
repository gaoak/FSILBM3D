!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Calculate the repulsive force between solids
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptForceR(nFish,dxmin,dymin,dzmin,nND,nND_max,nEL,nEL_max,ele,xyzful,repful)
    implicit none
    integer:: nFish,nEL_max,nND_max
    real(8):: dxmin,dymin,dzmin
    integer:: nND(1:nFish),nEL(1:nFish),ele(1:nFish,1:nEL_max,1:5)
    real(8):: xyzful(1:nFish,1:nND_max,1:6),repful(1:nFish,1:nND_max,1:6)
    !local
    integer:: iND,jND,iFish,jFish
    real(8):: delta_h,Phi,r(1:3),ds(1:3)
    real(8):: minx,miny,maxx,maxy,minz,maxz
    real(8):: xmin(1:nFish),xmax(1:nFish),ymin(1:nFish),ymax(1:nFish),zmin(1:nFish),zmax(1:nFish)

    repful(:,:,:)=0.d0
    ds(1)=dxmin
    ds(2)=dymin
    ds(3)=dzmin
    
    do iFish=1,nFish
        xmin(iFish) = minval(xyzful(iFish,1:nND(iFish),1))-dxmin*3.d0
        xmax(iFish) = maxval(xyzful(iFish,1:nND(iFish),1))+dxmin*3.d0
        ymin(iFish) = minval(xyzful(iFish,1:nND(iFish),2))-dymin*3.d0
        ymax(iFish) = maxval(xyzful(iFish,1:nND(iFish),2))+dymin*3.d0
        zmin(iFish) = minval(xyzful(iFish,1:nND(iFish),3))-dzmin*3.d0
        zmax(iFish) = maxval(xyzful(iFish,1:nND(iFish),3))+dzmin*3.d0   
    enddo

    do iFish=2,nFish
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
                if ( (minx>xyzful(iFish,iND,1)) .or. (xyzful(iFish,iND,1)>maxx)  &
                .or. (miny>xyzful(iFish,iND,2)) .or. (xyzful(iFish,iND,2)>maxy)  &
                .or. (minz>xyzful(iFish,iND,3)) .or. (xyzful(iFish,iND,3)>maxz)) then
                    cycle !point iND not in the overpalling region
                endif
                do jND=1,nND(jFish)
                    if ( (minx>xyzful(jFish,jND,1)) .or. (xyzful(jFish,jND,1)>maxx)  &
                    .or. (miny>xyzful(jFish,jND,2)) .or. (xyzful(jFish,jND,2)>maxy)  &
                    .or. (minz>xyzful(jFish,jND,3)) .or. (xyzful(jFish,jND,3)>maxz)) then
                        cycle !point jND not in the overpalling region
                    endif
                    r(1)=(xyzful(iFish,iND,1)-xyzful(jFish,jND,1))/dxmin
                    r(2)=(xyzful(iFish,iND,2)-xyzful(jFish,jND,2))/dymin
                    r(3)=(xyzful(iFish,iND,3)-xyzful(jFish,jND,3))/dzmin
                    delta_h=Phi(r(1))*Phi(r(2))*Phi(r(3))/dxmin/dymin/dzmin/sqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
                    repful(iFish,iND,1:3)=repful(iFish,iND,1:3) + delta_h*r(1:3)*ds(1:3) ! force
                    repful(jFish,jND,1:3)=repful(jFish,jND,1:3) - delta_h*r(1:3)*ds(1:3) ! reaction force
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
    IMPLICIT NONE
    integer:: zDim,yDim,xDim,nEL,nND,ele(nEL,5),ntolLBM
    real(8):: dz(zDim),dy(yDim),dx(xDim),dh,Uref,denIn,dtolLBM,dt,Palpha,Pbeta
    real(8):: force(zDim,yDim,xDim,1:3),uuu(zDim,yDim,xDim,1:3),den(zDim,yDim,xDim),xGrid(xDim),yGrid(yDim),zGrid(zDim)
    real(8):: xyzful(nND,6),xyzfulIB(nND,6),velful(nND,6),extful(nND,6)
    logical,intent(in):: isUniformGrid(1:3)
!==================================================================================================
    integer:: i,j,k,x,y,z,xbgn,ybgn,zbgn,xend,yend,zend,iEL,nt,iterLBM,iND
    real(8):: rx,ry,rz,Phi,dmaxLBM,dsum,invdh
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az
    real(8):: forceTemp(zDim,yDim,xDim,1:3),forceElemTemp(nEL,3)
    real(8):: forceElem(nEL,3),forceNode(nND,3),velfulIB(nND,3)
    real(8):: posElem(nEL,3),posElemIB(nEL,3),velElem(nEL,3),velElemIB(nEL,3),areaElem(nEL)
!==================================================================================================
!   compute velocity and displacement at IB nodes
    invdh = 1.D0/dh
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

        posElemIB(iEL,1:3)=(xyzfulIB(i,1:3)+xyzfulIB(j,1:3))/2.0d0

        ax =(x1-x2)
        ay =(y1-y2)
        az =(z1-z2)
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az)

        elseif(nt==3)then
        posElem(iEL,1:3)=(xyzful(i,1:3)+xyzful(j,1:3)+xyzful(k,1:3))/3.0d0
        velElem(iEL,1:3)=(velful(i,1:3)+velful(j,1:3)+velful(k,1:3))/3.0d0

        posElemIB(iEL,1:3)=(xyzfulIB(i,1:3)+xyzfulIB(j,1:3)+xyzfulIB(k,1:3))/3.0d0

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
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)  
    do  x = 1, xDim
    do  y = 1, yDim
    do  z = 1, zDim
        force(z,y,x,1:3)=0.0d0
     enddo
     enddo
     enddo
!$OMP END PARALLEL DO
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
        if(ele(iEL,4)==2)then
            forceElemTemp(iEL,1:3) = 0.0d0
        elseif(ele(iEL,4)==3)then
            forceElemTemp(iEL,1:3) = -Palpha*2.0*denIn*(posElem(iEL,1:3)-posElemIB(iEL,1:3))/dt*areaElem(iEL)*dh  &
                                     -Pbeta* 2.0*denIn*(velElem(iEL,1:3)-velElemIB(iEL,1:3))/dt*areaElem(iEL)*dh
        else
        endif
    enddo
!   ***********************************************************************************************
!   calculate Eulerian body force
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)  
    do  x = 1, xDim
    do  y = 1, yDim
    do  z = 1, zDim
        forceTemp(z,y,x,1:3)=0.0d0
    enddo
    enddo
    enddo
    !$OMP END PARALLEL DO
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
                    forceTemp(z,y,x,1:3)=forceTemp(z,y,x,1:3)-forceElemTemp(iEL,1:3)*rx*ry*rz*invdh*invdh*invdh
                enddo
            enddo
        enddo
    enddo
!   ***********************************************************************************************
!   update velocity
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)  
    do  x = 1, xDim
    do  y = 1, yDim
    do  z = 1, zDim         
        uuu(z,y,x,1:3)  = uuu(z,y,x,1:3)+0.5*dt*forceTemp(z,y,x,1:3)/den(z,y,x)
        force(z,y,x,1:3) = force(z,y,x,1:3) + forceTemp(z,y,x,1:3)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
!    force(1:zDim,1:yDim,1:xDim,1:3)=force(1:zDim,1:yDim,1:xDim,1:3)+forceTemp(1:zDim,1:yDim,1:xDim,1:3)
    forceElem(1:nEL,1:3) = forceElem(1:nEL,1:3)+forceElemTemp(1:nEL,1:3)   
!   ***********************************************************************************************
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

    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    �����������Ӧ������
!   ��������������Ӧ��
!    copyright@ RuNanHua 
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptStrs(zDim,yDim,xDim,nEL,nND,ele,dh,dx,dy,dz,mu,rr,u,p,xGrid,yGrid,zGrid,xyzful,extful)    
    IMPLICIT NONE
    integer:: zDim,yDim,xDim,nEL,nND,ele(nEL,5)
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
!   ����ڵ���
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
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION Phi(x)
    IMPLICIT NONE
    real(8)::Phi,x,r

    r=dabs(x)

    if(r<1.0d0)then
        Phi=(3.0-2.0*r+dsqrt(1.0+4.0*r-4.0*r*r))/8.0
    elseif(r<2.0d0)then
        Phi=(5.0-2.0*r-dsqrt(-7.0+12.0*r-4.0*r*r))/8.0
    else
        Phi=0.0d0
    endif
    ENDFUNCTION