!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Calculate the repulsive force between solids
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptForceR(dxmin,dymin,dzmin,nND,nND_max,xyzful,repful,nFish,Lspan)
    implicit none
    integer:: nFish,nND_max,nND(1:nFish)
    real(8):: dxmin,dymin,dzmin
    real(8):: xyzful(1:nND_max,1:6,1:nFish),repful(1:nND_max,1:6,1:nFish)
    !local
    integer:: iND,jND,iFish,jFish
    real(8):: delta_h,r(1:3),ds(1:3),phi_r(1:3),Lspan(nFish)
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
                    repful(iND,1:3,iFish)=repful(iND,1:3,iFish) + delta_h*r(1:3)*ds(1:3)*Lspan(iFish) ! force
                    repful(jND,1:3,jFish)=repful(jND,1:3,jFish) - delta_h*r(1:3)*ds(1:3)*Lspan(iFish) ! reaction force
                enddo !jND=1,nND(jFish)
            enddo !iND=1,nND(iFish)
        enddo !jFish=iFish+1,nFish
    enddo !iFish=1,nFish

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