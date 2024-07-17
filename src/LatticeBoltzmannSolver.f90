!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Calculate macro quantities
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE calculate_macro_quantities()
    USE simParam
    implicit none
    integer::x,y,z,ig
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z) 
    do  x = 1, xDim
    do  y = 1, yDim
    do  z = 1, zDim
        den(z,y,x  )  = SUM(fIn(z,y,x,0:lbmDim)) 
        uuu(z,y,x,1)  = (SUM(fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,1))+0.5d0*VolumeForce(1)*dt)/den(z,y,x)
        uuu(z,y,x,2)  = (SUM(fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,2))+0.5d0*VolumeForce(2)*dt)/den(z,y,x)
        uuu(z,y,x,3)  = (SUM(fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,3))+0.5d0*VolumeForce(3)*dt)/den(z,y,x)    
        prs(z,y,x)   = Cs2*(den(z,y,x)-denIn)
    enddo
    enddo
    enddo
    !$OMP END PARALLEL DO
    !prs(1:zDim,1:yDim,1:xDim)   = Cs2*(den(1:zDim,1:yDim,1:xDim)-denIn)
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    collision step SRT or MRT
!    collision model: single relexation time or multiple relexation time
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE collision_step()
    USE simParam
    implicit none 
    real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim),m(0:lbmDim),mEq(0:lbmDim),Flb(0:lbmDim)
    integer:: x,y,z

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,uSqr,uxyz,fEq,Flb)
    do    x = 1, xDim
    do    y = 1, yDim
    do    z = 1, zDim
        uSqr           = sum(uuu(z,y,x,1:3)**2)
        uxyz(0:lbmDim) = uuu(z,y,x,1) * ee(0:lbmDim,1) + uuu(z,y,x,2) * ee(0:lbmDim,2)+uuu(z,y,x,3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * den(z,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        Flb(0:lbmDim)  = dt*wt(0:lbmDim)*( &
                          (3.0*(ee(0:lbmDim,1)-uuu(z,y,x,1))+9.0*(ee(0:lbmDim,1)*uuu(z,y,x,1)+ee(0:lbmDim,2)*uuu(z,y,x,2)+ee(0:lbmDim,3)*uuu(z,y,x,3))*ee(0:lbmDim,1))*force(z,y,x,1) &
                         +(3.0*(ee(0:lbmDim,2)-uuu(z,y,x,2))+9.0*(ee(0:lbmDim,1)*uuu(z,y,x,1)+ee(0:lbmDim,2)*uuu(z,y,x,2)+ee(0:lbmDim,3)*uuu(z,y,x,3))*ee(0:lbmDim,2))*force(z,y,x,2) &
                         +(3.0*(ee(0:lbmDim,3)-uuu(z,y,x,3))+9.0*(ee(0:lbmDim,1)*uuu(z,y,x,1)+ee(0:lbmDim,2)*uuu(z,y,x,2)+ee(0:lbmDim,3)*uuu(z,y,x,3))*ee(0:lbmDim,3))*force(z,y,x,3) &
                                         )

        if    (iCollidModel==1)then
        !SRT collision
        fIn(z,y,x,0:lbmDim) = fIn(z,y,x,0:lbmDim) + Omega * (fEq(0:lbmDim)-fIn(z,y,x,0:lbmDim))+(1.0-0.5*Omega)*Flb(0:lbmDim)
        elseif(iCollidModel==2)then
        !MRT collision
        fIn(z,y,x,0:lbmDim)=fIn(z,y,x,0:lbmDim)+MATMUL( M_COLLID(0:lbmDim,0:lbmDim), fEq(0:lbmDim)-fIn(z,y,x,0:lbmDim) ) &
                                               +MATMUL( M_FORCE(0:lbmDim,0:lbmDim),Flb(0:lbmDim))
        else
            write(*,*)' collision_step Model is not defined'
        endif
    enddo
    enddo
    enddo
    !$OMP END PARALLEL DO
    END SUBROUTINE

    !0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    !0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0
    !0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1
    !0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1
    SUBROUTINE streams()
        USE simParam
        implicit none
        integer:: i
        integer:: strmDir(0:lbmDim,1:3)

        !stream direction
        strmDir(0:lbmDim,1)=-ee(0:lbmDim,3)
        strmDir(0:lbmDim,2)=-ee(0:lbmDim,2)
        strmDir(0:lbmDim,3)=-ee(0:lbmDim,1)

        do  i=0,lbmDim
            call swapzy(fIn, strmDir(i,3), strmDir(i,2), i, zDim, yDim, xDim, lbmDim)
            call swapx(fIn, strmDir(i,1), i, zDim, yDim, xDim, lbmDim)
        enddo
    END SUBROUTINE

    SUBROUTINE swapzy(f, dz, dy, i, zDim, yDim, xDim, lbmDim)
        implicit none
        integer, intent(in):: dz, dy, i, zDim, yDim, xDim, lbmDim
        real(8), intent(inout):: f(1:zDim,1:yDim,1:xDim,0:lbmDim)
        integer:: z, y, x
        real(8):: temp, tmpz(1:zDim)

        if(dz.eq.0 .and. dy.eq.0) return

        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,temp,tmpz)
        do  x = 1, xDim
            if(dz.eq.1) then
                do y=1, yDim
                    temp = f(zDim, y, x, i)
                    do z=zDim, 2, -1
                        f(z,y,x, i)=f(z-1,y,x, i)
                    enddo
                    f(1,y,x,i) = temp
                enddo
            elseif(dz.eq.-1) then
                do y=1, yDim
                    temp = f(1, y, x, i)
                    do z=1, zDim-1
                        f(z,y,x, i)=f(z+1,y,x, i)
                    enddo
                    f(zDim,y,x,i) = temp
                enddo
            endif
            if(dy.eq.1) then
                tmpz = f(:, yDim, x, i)
                do y=yDim, 2, -1
                    f(:,y,x, i)=f(:,y-1,x, i)
                enddo
                f(:,1,x,i) = tmpz
            elseif(dy.eq.-1) then
                tmpz = f(:, 1, x, i)
                do y=1, yDim-1
                    f(:,y,x, i)=f(:,y+1,x, i)
                enddo
                f(:,yDim,x,i) = tmpz
            endif
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE

    SUBROUTINE swapx(f, dx, i, zDim, yDim, xDim, lbmDim)
        USE PartitionXDim
        implicit none
        integer, intent(in):: dx, i, zDim, yDim, xDim, lbmDim
        real(8), intent(inout):: f(1:zDim, 1:yDim,1:xDim,0:lbmDim)
        integer:: p

        if(dx.eq.0) return

        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
        do  p = 1,npsize_copy
            if(dx .eq. -1) then
                call swapxwAtom(f, edge(:,:,p), eid(p), i, zDim, yDim, xDim, lbmDim, parindex(p), parindex(p+1)-1)
            elseif(dx.eq.1) then
                call swapxeAtom(f, edge(:,:,p), eid(p), i, zDim, yDim, xDim, lbmDim, parindex(p), parindex(p+1)-1)
            endif
        enddo
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
        do  p = 1,npsize_copy
            f(:,:,eid(p),i) = edge(:,:,p)
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE

    SUBROUTINE swapxeAtom(f, edge, eid, i, zDim, yDim, xDim, lbmDim, xbgn, xend)
        implicit none
        integer, intent(in):: i, zDim, yDim, xDim, lbmDim, xbgn, xend
        real(8), intent(inout):: f(1:zDim,1:yDim,1:xDim,0:lbmDim)
        real(8), intent(out):: edge(1:zDim,1:yDim)
        integer, intent(out):: eid
        integer:: x
        eid = xend+1
        if(eid .eq. xDim+1) eid = 1
        edge = f(:,:,xend,i)
        do  x = xend,xbgn+1,-1
            f(:,:,x,i) = f(:,:,x-1,i)
        enddo
    endsubroutine

    SUBROUTINE swapxwAtom(f, edge, eid, i, zDim, yDim, xDim, lbmDim, xbgn, xend)
        implicit none
        integer, intent(in):: i, zDim, yDim, xDim, lbmDim, xbgn, xend
        real(8), intent(inout):: f(1:zDim,1:yDim,1:xDim,0:lbmDim)
        real(8), intent(out):: edge(1:zDim,1:yDim)
        integer, intent(out):: eid
        integer:: x
        eid = xbgn-1
        if(eid .eq. 0) eid = xDim
        edge = f(:,:,xbgn,i)
        do  x = xbgn, xend-1
            f(:,:,x,i) = f(:,:,x+1,i)
        enddo
    endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Advection model: uniform grid advection, interpolation on the non-uniform grid
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE streaming_step()
    USE simParam
    implicit none
    integer:: x,y,z,i,k
    logical:: upwind,center,outer
    integer:: strmDir(0:lbmDim,1:3)

    !stream direction
    strmDir(0:lbmDim,1)=-ee(0:lbmDim,3)
    strmDir(0:lbmDim,2)=-ee(0:lbmDim,2)
    strmDir(0:lbmDim,3)=-ee(0:lbmDim,1)

!============================================
    if    (iStreamModel==1) then   !Advection uniform grid
!============================================
        call streams()
!============================================
    elseif(iStreamModel==2)then    !Interpolation
!============================================
         
        do  i=1,lbmDim
        fInTemp(1:zDim,1:yDim,1:xDim)=fIn(1:zDim,1:yDim,1:xDim,i)
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,upwind,center,outer,upxc0,upxcm,upxcmm,upyc0,upycm,upycmm,upzc0,upzcm,upzcmm,cnxc0,cnxcm,cnxcp,cnyc0,cnycm,cnycp,cnzc0,cnzcm,cnzcp)
        do  x=2,xDim-1
        do  y=2,yDim-1
        do  z=2,zDim-1
            !set logical flag
            upwind = (x>=3).and.(x<=xDim-2).and.(y>=3).and.(y<=yDim-2) .and. (z>=3) .and. (z<=zDim-2)                    
            center = (x>=2).and.(x<=xDim-1).and.(y>=2).and.(y<=yDim-1) .and. (z>=2) .and. (z<=zDim-1) .and. (.not.upwind) 
            outer  = .not.(upwind .or. center)

        !******************************************************************
        !******************************************************************
        if(upwind) then      
            !2nd-order upwind interpolation
            !2nd-order upwind interpolation coefficient
            if(ee(i,1)/=0 )then
            upxc0  = (xGrid(x)-ratio*dh*ee(i,1)-xGrid(x-ee(i,1)))*(xGrid(x)-ratio*dh*ee(i,1)-xGrid(x-2*ee(i,1)))/((xGrid(x          )-xGrid(x-ee(i,1)))*(xGrid(x          )-xGrid(x-2*ee(i,1))))
            upxcm  = (xGrid(x)-ratio*dh*ee(i,1)-xGrid(x        ))*(xGrid(x)-ratio*dh*ee(i,1)-xGrid(x-2*ee(i,1)))/((xGrid(x-  ee(i,1))-xGrid(x        ))*(xGrid(x-  ee(i,1))-xGrid(x-2*ee(i,1))))
            upxcmm = (xGrid(x)-ratio*dh*ee(i,1)-xGrid(x        ))*(xGrid(x)-ratio*dh*ee(i,1)-xGrid(x-  ee(i,1)))/((xGrid(x-2*ee(i,1))-xGrid(x        ))*(xGrid(x-2*ee(i,1))-xGrid(x-  ee(i,1))))
            else
            upxc0  = 1.0d0/3.0d0
            upxcm  = 1.0d0/3.0d0
            upxcmm = 1.0d0/3.0d0
            endif

            if(ee(i,2)/=0)then
            upyc0  = (yGrid(y)-ratio*dh*ee(i,2)-yGrid(y-ee(i,2)))*(yGrid(y)-ratio*dh*ee(i,2)-yGrid(y-2*ee(i,2)))/((yGrid(y          )-yGrid(y-ee(i,2)))*(yGrid(y          )-yGrid(y-2*ee(i,2))))
            upycm  = (yGrid(y)-ratio*dh*ee(i,2)-yGrid(y        ))*(yGrid(y)-ratio*dh*ee(i,2)-yGrid(y-2*ee(i,2)))/((yGrid(y-  ee(i,2))-yGrid(y        ))*(yGrid(y-  ee(i,2))-yGrid(y-2*ee(i,2))))
            upycmm = (yGrid(y)-ratio*dh*ee(i,2)-yGrid(y        ))*(yGrid(y)-ratio*dh*ee(i,2)-yGrid(y-  ee(i,2)))/((yGrid(y-2*ee(i,2))-yGrid(y        ))*(yGrid(y-2*ee(i,2))-yGrid(y-  ee(i,2))))
            else
            upyc0  = 1.0d0/3.0d0
            upycm  = 1.0d0/3.0d0
            upycmm = 1.0d0/3.0d0
            endif

            if(ee(i,3)/=0 )then
            upzc0  = (zGrid(z)-ratio*dh*ee(i,3)-zGrid(z-ee(i,3)))*(zGrid(z)-ratio*dh*ee(i,3)-zGrid(z-2*ee(i,3)))/((zGrid(z          )-zGrid(z-ee(i,3)))*(zGrid(z          )-zGrid(z-2*ee(i,3))))
            upzcm  = (zGrid(z)-ratio*dh*ee(i,3)-zGrid(z        ))*(zGrid(z)-ratio*dh*ee(i,3)-zGrid(z-2*ee(i,3)))/((zGrid(z-  ee(i,3))-zGrid(z        ))*(zGrid(z-  ee(i,3))-zGrid(z-2*ee(i,3))))
            upzcmm = (zGrid(z)-ratio*dh*ee(i,3)-zGrid(z        ))*(zGrid(z)-ratio*dh*ee(i,3)-zGrid(z-  ee(i,3)))/((zGrid(z-2*ee(i,3))-zGrid(z        ))*(zGrid(z-2*ee(i,3))-zGrid(z-  ee(i,3))))
            else
            upzc0  = 1.0d0/3.0d0
            upzcm  = 1.0d0/3.0d0
            upzcmm = 1.0d0/3.0d0
            endif

            fIn(z,y,x,i)= upzc0 *(upxc0 *( upyc0 *fInTemp(z, y, x          )+upycm *fInTemp(z, y-ee(i,2), x          )+upycmm *fInTemp(z, y-2*ee(i,2), x          ))+&
                                  upxcm *( upyc0 *fInTemp(z, y, x-  ee(i,1))+upycm *fInTemp(z, y-ee(i,2), x-  ee(i,1))+upycmm *fInTemp(z, y-2*ee(i,2), x-  ee(i,1)))+&
                                  upxcmm*( upyc0 *fInTemp(z, y, x-2*ee(i,1))+upycm *fInTemp(z, y-ee(i,2), x-2*ee(i,1))+upycmm *fInTemp(z, y-2*ee(i,2), x-2*ee(i,1))) &
                                 ) +&
                          upzcm *(upxc0 *( upyc0 *fInTemp(z-ee(i,3), y, x          )+upycm *fInTemp(z-ee(i,3), y-ee(i,2), x          )+upycmm *fInTemp(z-ee(i,3), y-2*ee(i,2), x          ))+&
                                  upxcm *( upyc0 *fInTemp(z-ee(i,3), y, x-  ee(i,1))+upycm *fInTemp(z-ee(i,3), y-ee(i,2), x-  ee(i,1))+upycmm *fInTemp(z-ee(i,3), y-2*ee(i,2), x-  ee(i,1)))+&
                                  upxcmm*( upyc0 *fInTemp(z-ee(i,3), y, x-2*ee(i,1))+upycm *fInTemp(z-ee(i,3), y-ee(i,2), x-2*ee(i,1))+upycmm *fInTemp(z-ee(i,3), y-2*ee(i,2), x-2*ee(i,1))) &
                                 ) +& 
                          upzcmm*(upxc0 *( upyc0 *fInTemp(z-2*ee(i,3), y, x          )+upycm *fInTemp(z-2*ee(i,3), y-ee(i,2), x          )+upycmm *fInTemp(z-2*ee(i,3), y-2*ee(i,2), x          ))+&
                                  upxcm *( upyc0 *fInTemp(z-2*ee(i,3), y, x-  ee(i,1))+upycm *fInTemp(z-2*ee(i,3), y-ee(i,2), x-  ee(i,1))+upycmm *fInTemp(z-2*ee(i,3), y-2*ee(i,2), x-  ee(i,1)))+&
                                  upxcmm*( upyc0 *fInTemp(z-2*ee(i,3), y, x-2*ee(i,1))+upycm *fInTemp(z-2*ee(i,3), y-ee(i,2), x-2*ee(i,1))+upycmm *fInTemp(z-2*ee(i,3), y-2*ee(i,2), x-2*ee(i,1))) &
                                 )            
        endif           
        !******************************************************************
        !******************************************************************
        if(center)then 
            !2nd-order central difference
            !2n-order central difference coefficient
            if(ee(i,1)/=0)then
            cnxcm = (xGrid(x)-ratio*dh*ee(i,1)-xGrid(x        ))*(xGrid(x)-ratio*dh*ee(i,1)-xGrid(x+ee(i,1)))/((xGrid(x-ee(i,1))-xGrid(x        ))*(xGrid(x-ee(i,1))-xGrid(x+ee(i,1))))
            cnxc0 = (xGrid(x)-ratio*dh*ee(i,1)-xGrid(x-ee(i,1)))*(xGrid(x)-ratio*dh*ee(i,1)-xGrid(x+ee(i,1)))/((xGrid(x        )-xGrid(x-ee(i,1)))*(xGrid(x        )-xGrid(x+ee(i,1))))
            cnxcp = (xGrid(x)-ratio*dh*ee(i,1)-xGrid(x-ee(i,1)))*(xGrid(x)-ratio*dh*ee(i,1)-xGrid(x        ))/((xGrid(x+ee(i,1))-xGrid(x-ee(i,1)))*(xGrid(x+ee(i,1))-xGrid(x        )))
            else
            cnxcm = 1.0d0/3.0d0
            cnxc0 = 1.0d0/3.0d0
            cnxcp = 1.0d0/3.0d0
            endif
            
            if(ee(i,2)/=0)then
            cnycm = (yGrid(y)-ratio*dh*ee(i,2)-yGrid(y        ))*(yGrid(y)-ratio*dh*ee(i,2)-yGrid(y+ee(i,2)))/((yGrid(y-ee(i,2))-yGrid(y        ))*(yGrid(y-ee(i,2))-yGrid(y+ee(i,2))))
            cnyc0 = (yGrid(y)-ratio*dh*ee(i,2)-yGrid(y-ee(i,2)))*(yGrid(y)-ratio*dh*ee(i,2)-yGrid(y+ee(i,2)))/((yGrid(y        )-yGrid(y-ee(i,2)))*(yGrid(y        )-yGrid(y+ee(i,2))))
            cnycp = (yGrid(y)-ratio*dh*ee(i,2)-yGrid(y-ee(i,2)))*(yGrid(y)-ratio*dh*ee(i,2)-yGrid(y        ))/((yGrid(y+ee(i,2))-yGrid(y-ee(i,2)))*(yGrid(y+ee(i,2))-yGrid(y        )))
            else
            cnycm = 1.0d0/3.0d0
            cnyc0 = 1.0d0/3.0d0
            cnycp = 1.0d0/3.0d0    
            endif

            if(ee(i,3)/=0)then
            cnzcm = (zGrid(z)-ratio*dh*ee(i,3)-zGrid(z        ))*(zGrid(z)-ratio*dh*ee(i,3)-zGrid(z+ee(i,3)))/((zGrid(z-ee(i,3))-zGrid(z        ))*(zGrid(z-ee(i,3))-zGrid(z+ee(i,3))))
            cnzc0 = (zGrid(z)-ratio*dh*ee(i,3)-zGrid(z-ee(i,3)))*(zGrid(z)-ratio*dh*ee(i,3)-zGrid(z+ee(i,3)))/((zGrid(z        )-zGrid(z-ee(i,3)))*(zGrid(z        )-zGrid(z+ee(i,3))))
            cnzcp = (zGrid(z)-ratio*dh*ee(i,3)-zGrid(z-ee(i,3)))*(zGrid(z)-ratio*dh*ee(i,3)-zGrid(z        ))/((zGrid(z+ee(i,3))-zGrid(z-ee(i,3)))*(zGrid(z+ee(i,3))-zGrid(z        )))
            else
            cnzcm = 1.0d0/3.0d0
            cnzc0 = 1.0d0/3.0d0
            cnzcp = 1.0d0/3.0d0    
            endif
                
            fIn(z,y,x,i)= cnzc0 *(cnxc0 *(cnyc0*fInTemp(z,y, x        )+cnycm*fInTemp(z,y-ee(i,2),x        )+cnycp*fInTemp(z,y+ee(i,2),x        ))+&
                                  cnxcm *(cnyc0*fInTemp(z,y, x-ee(i,1))+cnycm*fInTemp(z,y-ee(i,2),x-ee(i,1))+cnycp*fInTemp(z,y+ee(i,2),x-ee(i,1)))+&
                                  cnxcp *(cnyc0*fInTemp(z,y, x+ee(i,1))+cnycm*fInTemp(z,y-ee(i,2),x+ee(i,1))+cnycp*fInTemp(z,y+ee(i,2),x+ee(i,1))) &
                                 )+ &
                          cnzcm *(cnxc0 *(cnyc0*fInTemp(z-ee(i,3),y, x        )+cnycm*fInTemp(z-ee(i,3),y-ee(i,2),x        )+cnycp*fInTemp(z-ee(i,3),y+ee(i,2),x        ))+&
                                  cnxcm *(cnyc0*fInTemp(z-ee(i,3),y, x-ee(i,1))+cnycm*fInTemp(z-ee(i,3),y-ee(i,2),x-ee(i,1))+cnycp*fInTemp(z-ee(i,3),y+ee(i,2),x-ee(i,1)))+&
                                  cnxcp *(cnyc0*fInTemp(z-ee(i,3),y, x+ee(i,1))+cnycm*fInTemp(z-ee(i,3),y-ee(i,2),x+ee(i,1))+cnycp*fInTemp(z-ee(i,3),y+ee(i,2),x+ee(i,1))) &
                                 )+ &
                          cnzcp *(cnxc0 *(cnyc0*fInTemp(z+ee(i,3),y, x        )+cnycm*fInTemp(z+ee(i,3),y-ee(i,2),x        )+cnycp*fInTemp(z+ee(i,3),y+ee(i,2),x        ))+&
                                  cnxcm *(cnyc0*fInTemp(z+ee(i,3),y, x-ee(i,1))+cnycm*fInTemp(z+ee(i,3),y-ee(i,2),x-ee(i,1))+cnycp*fInTemp(z+ee(i,3),y+ee(i,2),x-ee(i,1)))+&
                                  cnxcp *(cnyc0*fInTemp(z+ee(i,3),y, x+ee(i,1))+cnycm*fInTemp(z+ee(i,3),y-ee(i,2),x+ee(i,1))+cnycp*fInTemp(z+ee(i,3),y+ee(i,2),x+ee(i,1))) &
                                 )
         endif
            
         !******************************************************************
         !******************************************************************
         if(outer)then
         endif

        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO
        enddo
!============================================
    else
    endif
!============================================
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    far field boundary set     
!    set all far-field boundary conditions as equilibrium function
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE set_equilibrium_farfld_BC()
    USE simParam
    implicit none
    integer:: i,k
    integer:: x, y, z,ig
    logical:: upwind,center,outer
    real(8):: uSqr ,uxyz(0:lbmDim) ,fEq(0:lbmDim), vel(1:SpcDim)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,upwind,center,outer,vel)
    do  x=1,xDim 
    do  y=1,yDim 
    do  z=1,zDim 
        !set logical flag
        if(iBC==1)then
            if      ( outer  ) then               !outer most layer
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(x), yGrid(y), zGrid(z), vel)
            elseif(VelocityKind==2) then
                call evaluateOscillatoryVelocity(vel)
            endif
            uSqr           = sum(vel(1:3)**2)
            uxyz(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)+vel(3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
            fIn(z,y,x,0:lbmDim)= fEq(0:lbmDim)
            endif
        elseif(iBC==2)then
            if      ( outer .or. center ) then    !second outer-most layer
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(x), yGrid(y), zGrid(z), vel)
            elseif(VelocityKind==2) then
                call evaluateOscillatoryVelocity(vel)
            endif
            uSqr           = sum(vel(1:3)**2)
            uxyz(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)+vel(3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
            fIn(z,y,x,0:lbmDim)= fEq(0:lbmDim)
            endif
        else
            write(*,*)'BC is not defined!'
            stop
        endif
    enddo
    enddo
    enddo  
    !$OMP END PARALLEL DO     
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    far field boundary set     
!    set far-field boundary conditions, other types
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE set_other_farfld_BCs()
    USE simParam
    implicit none
    integer:: i,k
    integer:: x, y, z,ig
    real(8):: uSqr ,uxyz(0:lbmDim) ,fEq(0:lbmDim)
    real(8):: uSqri,uxyzi(0:lbmDim),fEqi(0:lbmDim)
    real(8):: fTmp(0:lbmDim), vel(1:SpcDim)   
!    ---------------------------------------
!    ---------------------------------------

    !=======x-direction
    if(xMinBC==DirecletUP)then
        do  y = 1, yDim
        do  z = 1, zDim
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(1), yGrid(y), zGrid(z), vel)
            elseif(VelocityKind==2) then
                call evaluateOscillatoryVelocity(vel)
            endif
            uSqr           = sum(vel(1:3)**2)
            uxyz(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)+vel(3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)  
            fIn(z,y,1,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(xMinBC==DirecletUU)then
        do  y = 1, yDim
        do  z = 1, zDim
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(1), yGrid(y), zGrid(z), vel)
            elseif(VelocityKind==2) then
                call evaluateOscillatoryVelocity(vel)
            endif
            uSqr           = sum(vel(1:3)**2)
            uxyz(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)+vel(3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * den(z,y,2) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

            uSqri           = sum(uuu(z,y,2,1:3)**2)
            uxyzi(0:lbmDim) = uuu(z,y,2,1) * ee(0:lbmDim,1) + uuu(z,y,2,2) * ee(0:lbmDim,2)+uuu(z,y,2,3) * ee(0:lbmDim,3)
            fEqi(0:lbmDim)  = wt(0:lbmDim) * den(z,y,2) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

            fIn(z,y,1,[1,7,9,11,13]) =fEq([1,7,9,11,13])+ (fIn(z,y,2,[1,7,9,11,13])-fEqi([1,7,9,11,13]))
        enddo
        enddo
    elseif(xMinBC==Advection1)then
        fIn(:,:,1,[1,7,9,11,13]) = fIn(:,:,2,[1,7,9,11,13])
    elseif(xMinBC==Advection2)then
        fIn(:,:,1,[1,7,9,11,13]) = 2.0*fIn(:,:,2,[1,7,9,11,13])-fIn(:,:,3,[1,7,9,11,13])
    elseif(xMinBC==Periodic .or. xMinBC==fluid )then
        ! no need set
    elseif(xMinBC==wall) then
        do  y = 1, yDim
        do  z = 1, zDim
            fTmp([1,7,9,11,13]) = fIn(z,y,1,oppo([1,7,9,11,13]))
            fIn(z,y,1,[1,7,9,11,13]) = fTmp([1,7,9,11,13])
        enddo
        enddo
    else
        stop 'no define BC'
    endif

    if(xMaxBC==DirecletUP)then
        do  y = 1, yDim
        do  z = 1, zDim
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(xDim), yGrid(y), zGrid(z), vel)
            elseif(VelocityKind==2) then
                call evaluateOscillatoryVelocity(vel)
            endif
            uSqr           = sum(vel(1:3)**2)
            uxyz(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)+vel(3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
            fIn(z,y,1,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(xMaxBC==DirecletUU)then
        do  y = 1, yDim
        do  z = 1, zDim
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(xDim), yGrid(y), zGrid(z), vel)
            elseif(VelocityKind==2) then
                call evaluateOscillatoryVelocity(vel)
            endif
            uSqr           = sum(vel(1:3)**2)
            uxyz(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)+vel(3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * den(z,y,xDim-1) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

            uSqri           = sum(uuu(z,y,xDim-1,1:3)**2)
            uxyzi(0:lbmDim) = uuu(z,y,xDim-1,1) * ee(0:lbmDim,1) + uuu(z,y,xDim-1,2) * ee(0:lbmDim,2)+uuu(z,y,xDim-1,3) * ee(0:lbmDim,3)
            fEqi(0:lbmDim)  = wt(0:lbmDim) * den(z,y,xDim-1) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

            fIn(z,y,xDim,[2,8,10,12,14]) =fEq([2,8,10,12,14])+ (fIn(z,y,xDim-1,[2,8,10,12,14])-fEqi([2,8,10,12,14]))
        enddo
        enddo
    elseif(xMaxBC==Advection1)then
        fIn(:,:,xDim,[2,8,10,12,14]) = fIn(:,:,xDim-1,[2,8,10,12,14])
    elseif(xMaxBC==Advection2)then
        fIn(:,:,xDim,[2,8,10,12,14]) = 2.0*fIn(:,:,xDim-1,[2,8,10,12,14])-fIn(:,:,xDim-2,[2,8,10,12,14])
    elseif(xMaxBC==Periodic .or. xMaxBC==fluid)then
        ! no need set
    elseif(xMaxBC==wall) then
        do  y = 1, yDim
        do  z = 1, zDim
            fTmp([2,8,10,12,14]) = fIn(z,y,xDim,oppo([2,8,10,12,14]))
            fIn(z,y,xDim,[2,8,10,12,14]) = fTmp([2,8,10,12,14])
        enddo
        enddo
    else
        stop 'no define BC'
    endif

    !=======y-direction
    if(yMinBC==DirecletUP)then
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        do  x = 1, xDim
        do  z = 1, zDim
        fIn(z,1,x,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(yMinBC==DirecletUU)then
        do  x = 1, xDim
        do  z = 1, zDim
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * den(z,2,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        uSqri           = sum(uuu(z,2,x,1:3)**2)
        uxyzi(0:lbmDim) = uuu(z,2,x,1) * ee(0:lbmDim,1) + uuu(z,2,x,2) * ee(0:lbmDim,2)+uuu(z,2,x,3) * ee(0:lbmDim,3)
        fEqi(0:lbmDim)  = wt(0:lbmDim) * den(z,2,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        fIn(z,1,x,[3,7,8,15,17]) =fEq([3,7,8,15,17])+ (fIn(z,2,x,[3,7,8,15,17])-fEqi([3,7,8,15,17]))
        enddo
        enddo
    elseif(yMinBC==Advection1)then
        fIn(:,1,:,[3,7,8,15,17]) = fIn(:,2,:,[3,7,8,15,17])
    elseif(yMinBC==Advection2)then
        fIn(:,1,:,[3,7,8,15,17]) = 2.0*fIn(:,2,:,[3,7,8,15,17])-fIn(:,3,:,[3,7,8,15,17])
    elseif(yMinBC==Periodic .or. yMinBC==fluid)then
        ! no need set
    elseif(yMinBC==wall) then
        do  x = 1, xDim
        do  z = 1, zDim
            fTmp([3,7,8,15,17]) = fIn(z,1,x,oppo([3,7,8,15,17]))
            fIn(z,1,x,[3,7,8,15,17]) = fTmp([3,7,8,15,17])
        enddo
        enddo
    else
        stop 'no define BC'
    endif

    if(yMaxBC==DirecletUP)then
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        do  x = 1, xDim
        do  z = 1, zDim
        fIn(z,yDim,x,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(yMaxBC==DirecletUU)then
        do  x = 1, xDim
        do  z = 1, zDim
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * den(z,yDim-1,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        uSqri           = sum(uuu(z,yDim-1,x,1:3)**2)
        uxyzi(0:lbmDim) = uuu(z,yDim-1,x,1) * ee(0:lbmDim,1) + uuu(z,yDim-1,x,2) * ee(0:lbmDim,2)+uuu(z,yDim-1,x,3) * ee(0:lbmDim,3)
        fEqi(0:lbmDim)  = wt(0:lbmDim) * den(z,yDim-1,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        fIn(z,yDim,x,[4,9,10,16,18]) =fEq([4,9,10,16,18])+ (fIn(z,yDim-1,x,[4,9,10,16,18])-fEqi([4,9,10,16,18]))
        enddo
        enddo
    elseif(yMaxBC==Advection1)then
        fIn(:,yDim,:,[4,9,10,16,18]) = fIn(:,yDim-1,:,[4,9,10,16,18])
    elseif(yMaxBC==Advection2)then
        fIn(:,yDim,:,[4,9,10,16,18]) = 2.0*fIn(:,yDim-1,:,[4,9,10,16,18])-fIn(:,yDim-2,:,[4,9,10,16,18])
    elseif(yMaxBC==Periodic .or. yMaxBC==fluid)then
        ! no need set
    elseif(yMaxBC==wall) then
        do  x = 1, xDim
        do  z = 1, zDim
            fTmp([4,9,10,16,18]) = fIn(z,yDim,x,oppo([4,9,10,16,18]))
            fIn(z,yDim,x,[4,9,10,16,18]) = fTmp([4,9,10,16,18])
        enddo
        enddo
    else
        stop 'no define BC'
    endif

    !=======z-direction
    if(zMinBC==DirecletUP)then
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        do  x = 1, xDim
        do  y = 1, yDim
        fIn(1,y,x,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(zMinBC==DirecletUU)then
        do  x = 1, xDim
        do  y = 1, yDim
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * den(2,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        uSqri           = sum(uuu(2,y,x,1:3)**2)
        uxyzi(0:lbmDim) = uuu(2,y,x,1) * ee(0:lbmDim,1) + uuu(2,y,x,2) * ee(0:lbmDim,2)+uuu(2,y,x,3) * ee(0:lbmDim,3)
        fEqi(0:lbmDim)  = wt(0:lbmDim) * den(1,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        fIn(1,y,x,[5,11,12,15,16]) =fEq([5,11,12,15,16])+ (fIn(2,y,x,[5,11,12,15,16])-fEqi([5,11,12,15,16]))
        enddo
        enddo
    elseif(zMinBC==Advection1)then
        fIn(1,:,:,[5,11,12,15,16]) = fIn(2,:,:,[5,11,12,15,16])
    elseif(zMinBC==Advection2)then
        fIn(1,:,:,[5,11,12,15,16]) = 2.0*fIn(2,:,:,[5,11,12,15,16])-fIn(3,:,:,[5,11,12,15,16])
    elseif(zMinBC==Periodic .or. zMinBC==fluid)then
        ! no need set
    elseif(zMinBC==wall) then
        do  x = 1, xDim
        do  y = 1, yDim
            fTmp([5,11,12,15,16]) = fIn(1,y,x,oppo([5,11,12,15,16]))
            fIn(1,y,x,[5,11,12,15,16]) = fTmp([5,11,12,15,16])
        enddo
        enddo
    elseif(zMinBC==movingWall) then
        do  x = 1, xDim
        do  y = 1, yDim
            if(MovingKind1==0) then
                if(VelocityKind==0) then
                    call evaluateShearVelocity(xGrid(x),yGrid(y),zGrid(1),vel)
                elseif(VelocityKind==2) then
                    call evaluateOscillatoryVelocity(vel)
                endif
            elseif(MovingKind1==1) then
                vel(1)=MovingVel1
                vel(2)=0.0d0
                vel(3)=0.0d0
            elseif(MovingKind1==2) then
                vel(1)=MovingVel1 * dcos(2*pi*MovingFreq1*time)
                vel(2)=0.0d0
                vel(3)=0.0d0             
            endif
            fTmp(5) = fIn(1,y,x,oppo(5)) + 2.0 * wt(5) * denIn * (ee(5,1) * vel(1) + ee(5,2) * vel(2) + ee(5,3) * vel(3)) * 3.0
            fTmp(11) = fIn(1,y,x,oppo(11)) + 2.0 * wt(11) * denIn * (ee(11,1) * vel(1) + ee(11,2) * vel(2) + ee(11,3) * vel(3)) * 3.0
            fTmp(12) = fIn(1,y,x,oppo(12)) + 2.0 * wt(12) * denIn * (ee(12,1) * vel(1) + ee(12,2) * vel(2) + ee(12,3) * vel(3)) * 3.0
            fTmp(15) = fIn(1,y,x,oppo(15)) + 2.0 * wt(15) * denIn * (ee(15,1) * vel(1) + ee(15,2) * vel(2) + ee(15,3) * vel(3)) * 3.0
            fTmp(16) = fIn(1,y,x,oppo(16)) + 2.0 * wt(16) * denIn * (ee(16,1) * vel(1) + ee(16,2) * vel(2) + ee(16,3) * vel(3)) * 3.0
            fIn(1,y,x,[5,11,12,15,16]) = fTmp([5,11,12,15,16])
        enddo
        enddo
    else
        stop 'no define BC'
    endif

    if(zMaxBC==DirecletUP)then
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        do  x = 1, xDim
        do  y = 1, yDim
        fIn(zDim,y,x,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(zMaxBC==DirecletUU)then
        do  x = 1, xDim
        do  y = 1, yDim
        uSqr           = sum(uuuIn(1:3)**2)
        uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * den(zDim-1,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        uSqri           = sum(uuu(zDim-1,y,x,1:3)**2)
        uxyzi(0:lbmDim) = uuu(zDim-1,y,x,1) * ee(0:lbmDim,1) + uuu(zDim-1,y,x,2) * ee(0:lbmDim,2)+uuu(zDim-1,y,x,3) * ee(0:lbmDim,3)
        fEqi(0:lbmDim)  = wt(0:lbmDim) * den(zDim-1,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)

        fIn(zDim,y,x,[6,13,14,17,18]) =fEq([6,13,14,17,18])+ (fIn(zDim-1,y,x,[6,13,14,17,18])-fEqi([6,13,14,17,18]))
        enddo
        enddo
    elseif(zMaxBC==Advection1)then
        fIn(zDim,:,:,[6,13,14,17,18]) = fIn(zDim-1,:,:,[6,13,14,17,18])
    elseif(zMaxBC==Advection2)then
        fIn(zDim,:,:,[6,13,14,17,18]) = 2.0*fIn(zDim-1,:,:,[6,13,14,17,18])-fIn(zDim-2,:,:,[6,13,14,17,18])
    elseif(zMaxBC==Periodic .or. zMaxBC==fluid)then
        ! no need set
    elseif(zMaxBC==wall) then
        do  x = 1, xDim
        do  y = 1, yDim
            fTmp([6,13,14,17,18]) = fIn(zDim,y,x,oppo([6,13,14,17,18]))
            fIn(zDim,y,x,[6,13,14,17,18]) = fTmp([6,13,14,17,18])
        enddo
        enddo
    elseif(zMaxBC==movingWall) then
        do  x = 1, xDim
        do  y = 1, yDim
            if(MovingKind2==0) then
                if(VelocityKind==0) then
                    call evaluateShearVelocity(xGrid(x),yGrid(y),zGrid(zDim),vel)
                elseif(VelocityKind==2) then
                    call evaluateOscillatoryVelocity(vel)
                endif
            elseif(MovingKind2==1) then
                vel(1)=MovingVel2
                vel(2)=0.0d0
                vel(3)=0.0d0
            elseif(MovingKind2==2) then
                vel(1)=MovingVel2 * dcos(2*pi*MovingFreq2*time)
                vel(2)=0.0d0
                vel(3)=0.0d0
            endif
            fTmp(6) = fIn(zDim,y,x,oppo(6)) + 2.0 * wt(6) * denIn * (ee(6,1) * vel(1) + ee(6,2) * vel(2) + ee(6,3) * vel(3)) * 3.0
            fTmp(13) = fIn(zDim,y,x,oppo(13)) + 2.0 * wt(13) * denIn * (ee(13,1) * vel(1) + ee(13,2) * vel(2) + ee(13,3) * vel(3)) * 3.0
            fTmp(14) = fIn(zDim,y,x,oppo(14)) + 2.0 * wt(14) * denIn * (ee(14,1) * vel(1) + ee(14,2) * vel(2) + ee(14,3) * vel(3)) * 3.0
            fTmp(17) = fIn(zDim,y,x,oppo(17)) + 2.0 * wt(17) * denIn * (ee(17,1) * vel(1) + ee(17,2) * vel(2) + ee(17,3) * vel(3)) * 3.0
            fTmp(18) = fIn(zDim,y,x,oppo(18)) + 2.0 * wt(18) * denIn * (ee(18,1) * vel(1) + ee(18,2) * vel(2) + ee(18,3) * vel(3)) * 3.0
            fIn(zDim,y,x,[6,13,14,17,18]) = fTmp([6,13,14,17,18])
        enddo
        enddo
    else
        stop 'no define BC'
    endif
!   ==========================
    END SUBROUTINE
