!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptMacr()
    USE simParam
    implicit none       
    integer:: x,y
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y) 
    do    x = 1, xDim
    do    y = 1, yDim
        den(y,x)    = SUM(fIn(y,x,:))
        uuu(y,x,1)  = SUM(fIn(y,x,:)*ee(:,1))/den(y,x)
        uuu(y,x,2)  = SUM(fIn(y,x,:)*ee(:,2))/den(y,x)  
        prs(y,x)    = Cs2*(den(y,x)-denIn)
    enddo
    enddo
    !$OMP END PARALLEL DO
    
    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE collide()
    USE simParam
    implicit none
    real(8):: uSqr,uxy(0:lbmDim),fEq(0:lbmDim),Flb(0:lbmDim)
    integer:: x,y
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,uSqr,uxy,fEq,Flb)
    do    x = 1, xDim 
    do    y = 1, yDim 
        uSqr          = sum(uuu(y,x,1:spcDim)**2) 
        uxy(0:lbmDim) = uuu(y,x,1) * ee(0:lbmDim,1) + uuu(y,x,2) * ee(0:lbmDim,2)
        fEq(0:lbmDim) = wt(0:lbmDim) * den(y,x) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)
        Flb(0:lbmDim) =  dt*wt(0:lbmDim)*( &
                (3.0*(ee(0:lbmDim,1)-uuu(y,x,1))+9.0*(ee(0:lbmDim,1)*uuu(y,x,1)+ee(0:lbmDim,2)*uuu(y,x,2))*ee(0:lbmDim,1))*force(y,x,1) &
               +(3.0*(ee(0:lbmDim,2)-uuu(y,x,2))+9.0*(ee(0:lbmDim,1)*uuu(y,x,1)+ee(0:lbmDim,2)*uuu(y,x,2))*ee(0:lbmDim,2))*force(y,x,2) &
                                         )
        if    (iCollidModel==1)then
        !single relaxation model
        fIn(y,x,0:lbmDim) = fIn(y,x,0:lbmDim) + Omega * (fEq(0:lbmDim)-fIn(y,x,0:lbmDim))+(1.0-0.5*Omega)*Flb(0:lbmDim)
        elseif(iCollidModel==2)then
        !multiple relaxation model
        fIn(y,x,0:lbmDim) = fIn(y,x,0:lbmDim) + MATMUL( M_COLLID(0:lbmDim,0:lbmDim), fEq(0:lbmDim)-fIn(y,x,0:lbmDim) )  &
                                              + MATMUL( M_FORCE(0:lbmDim,0:lbmDim), Flb(0:lbmDim))
        else
        endif
        
    enddo
    enddo
    !$OMP END PARALLEL DO
    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE streams()
    USE simParam
    implicit none
    real(8):: f_hlp(1:yDim,1:xDim,0:lbmDim)
    integer  x,y,x_e,x_w,y_n,y_s
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,y_n,x_e,y_s,x_w)
    do  x = 1, xDim
    do  y = 1, yDim
          y_n = mod(y,yDim) + 1
          x_e = mod(x,xDim) + 1
          y_s = yDim - mod(yDim + 1 - y, yDim)
          x_w = xDim - mod(xDim + 1 - x, xDim)
          
          f_hlp(y,x_e,1) = fIn(y,x,1)
          f_hlp(y_n,x,2) = fIn(y,x,2)
          f_hlp(y,x_w,3) = fIn(y,x,3)
          f_hlp(y_s,x,4) = fIn(y,x,4)
          f_hlp(y_n,x_e,5) = fIn(y,x,5)
          f_hlp(y_n,x_w,6) = fIn(y,x,6)
          f_hlp(y_s,x_w,7) = fIn(y,x,7)
          f_hlp(y_s,x_e,8) = fIn(y,x,8)

    enddo
    enddo
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y)   
    do  x = 1, xDim
    do  y = 1, yDim
             fIn(y,x,1) = f_hlp(y,x,1)
             fIn(y,x,2) = f_hlp(y,x,2)
             fIn(y,x,3) = f_hlp(y,x,3)
             fIn(y,x,4) = f_hlp(y,x,4)
             fIn(y,x,5) = f_hlp(y,x,5)
             fIn(y,x,6) = f_hlp(y,x,6)
             fIn(y,x,7) = f_hlp(y,x,7)
             fIn(y,x,8) = f_hlp(y,x,8)
             
   enddo
   enddo
   !$OMP END PARALLEL DO
   END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE stream()
    USE simParam
    implicit none
    integer:: x, y,i,k
    logical:: upwind,center,outer
    integer:: strmDir(0:lbmDim,1:2)
    !note: derections of strmDir(:) & physical spcace, PengGL, Fortran95, Page 571, CSHIFT usage
    strmDir(0:lbmDim,1) =-ee(0:lbmDim,2)
    strmDir(0:lbmDim,2) =-ee(0:lbmDim,1) 
!============================================
    if(iStreamModel==1) then   
!============================================
        call streams()
!============================================
    elseif(iStreamModel==2)then    ! interpolation
!============================================             
        do  i = 1, lbmDim
        fInTemp(1:yDim,1:xDim)=fIn(1:yDim,1:xDim,i)
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,upwind,center,outer,upxc0,upxcm,upxcmm,upyc0,upycm,upycmm,cnxc0,cnxcm,cnxcp,cnyc0,cnycm,cnycp)
        do  x = 1, xDim
        do  y = 1, yDim
            !point localtion
            upwind =  (x>=3).and.(x<=xDim-2).and.(y>=3).and.(y<=yDim-2)                    
            center = ((x>=2).and.(x<=xDim-1).and.(y>=2).and.(y<=yDim-1)) .and. (.not.upwind) 
            outer  = .not.(upwind .or. center)
            !******************************************************************
            if(upwind) then
            !2nd-order upwind scheme coefficients
            if(ee(i,1)/=0)then
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

            fIn(y,x,i)= upxc0 *( upyc0*fInTemp(y, x          )+upycm*fInTemp(y-ee(i,2), x          )+upycmm*fInTemp(y-2*ee(i,2), x          ))+&
                        upxcm *( upyc0*fInTemp(y, x-  ee(i,1))+upycm*fInTemp(y-ee(i,2), x-  ee(i,1))+upycmm*fInTemp(y-2*ee(i,2), x-  ee(i,1)))+&
                        upxcmm*( upyc0*fInTemp(y, x-2*ee(i,1))+upycm*fInTemp(y-ee(i,2), x-2*ee(i,1))+upycmm*fInTemp(y-2*ee(i,2), x-2*ee(i,1)))     
        
            endif
            !******************************************************************
            if(center)then 
            !2nd-order center scheme coefficients
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
                
            fIn(y,x,i)= cnxc0*(cnyc0*fInTemp(y, x        )+cnycm*fInTemp(y-ee(i,2),x        )+cnycp*fInTemp(y+ee(i,2),x        ))+&
                        cnxcm*(cnyc0*fInTemp(y, x-ee(i,1))+cnycm*fInTemp(y-ee(i,2),x-ee(i,1))+cnycp*fInTemp(y+ee(i,2),x-ee(i,1)))+&
                        cnxcp*(cnyc0*fInTemp(y, x+ee(i,1))+cnycm*fInTemp(y-ee(i,2),x+ee(i,1))+cnycp*fInTemp(y+ee(i,2),x+ee(i,1)))
            endif

            !******************************************************************
            !******************************************************************
            if(outer)then 
            endif

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
!    
!  
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE setBund()  !Boundary set for uniform grid
    USE simParam
    implicit none  
    integer:: x, y
    logical:: upwind,center,outer
    real(8):: uSqr,uxy(0:lbmDim),fEq(0:lbmDim)
    real(8):: vel(1:SpcDim)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,upwind,center,outer,vel)
    do  x = 1, xDim
    do  y = 1, yDim
        !set logic
        upwind = (x>=3).and.(x<=xDim-2).and.(y>=3).and.(y<=yDim-2)                    
        center = ((x>=2).and.(x<=xDim-1).and.(y>=2).and.(y<=yDim-1)) .and. (.not.upwind) 
        outer  = .not.(upwind .or. center)
        if (((iBC==1).and.(outer)) .or. ((iBC==2).and.( outer .or. center )))  then
            !out most layer, or second out most layer
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(x), yGrid(y), vel)
            elseif(VelocityKind==1) then
                call evaluateParabolicVelocity(xGrid(x), yGrid(y), vel)
            elseif(VelocityKind==2) then
                vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
                vel(2)=uIn(2)
            endif
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim) = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)
            fIn(y,x,0:lbmDim) = fEq(0:lbmDim)
        endif
    enddo
    enddo
    !$OMP END PARALLEL DO

    END SUBROUTINE
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   far field boundary set     
!   
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE setBund2()  !Boundary set for non-uniform grid
    USE simParam
    implicit none
    integer:: x, y
    real(8):: uSqr ,uxy(0:lbmDim) ,fEq(0:lbmDim)
    real(8):: uSqri,uxyi(0:lbmDim),fEqi(0:lbmDim)
    real(8):: fTmp(0:lbmDim), vel(1:SpcDim)
!    ---------------------------------------
!    ---------------------------------------

    !=======x-direction
    if(xMinBC==DirecletUP)then
        do  y = 1, yDim
        if(VelocityKind==0) then
            call evaluateShearVelocity(xGrid(1), yGrid(y), vel)
        elseif(VelocityKind==1) then
            call evaluateParabolicVelocity(xGrid(1), yGrid(y), vel)
        elseif(VelocityKind==2) then
            vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
            vel(2)=uIn(2)
        endif
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim) = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)
            fIn(y,1,0:lbmDim)   = fEq(0:lbmDim)
        enddo
    elseif(xMinBC==DirecletUU)then
        do  y = 1, yDim
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(1), yGrid(y), vel)
            elseif(VelocityKind==1) then
                call evaluateParabolicVelocity(xGrid(1), yGrid(y), vel)
            elseif(VelocityKind==2) then
            vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
            vel(2)=uIn(2)
            endif
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim) = wt(0:lbmDim) * den(y,2) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            uSqri           = sum(uuu(y,2,1:2)**2)
            uxyi(0:lbmDim) = uuu(y,2,1) * ee(0:lbmDim,1) + uuu(y,2,2) * ee(0:lbmDim,2)
            fEqi(0:lbmDim)  = wt(0:lbmDim) * den(y,2) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            fIn(y,1,[1,5,8]) =fEq([1,5,8])+ (fIn(y,2,[1,5,8])-fEqi([1,5,8]))
        enddo
    elseif(xMinBC==Advection1)then
        fIn(:,1,[1,5,8]) = fIn(:,2,[1,5,8])
    elseif(xMinBC==Advection2)then
        fIn(:,1,[1,5,8]) = 2.0*fIn(:,2,[1,5,8])-fIn(:,3,[1,5,8])
    elseif(xMinBC==Periodic .or. xMinBC==fluid )then
        ! no need set
    elseif(xMinBC==wall) then
        do  y = 1, yDim
            fTmp([1,5,8]) = fIn(y,1,oppo([1,5,8]))
            fIn(y,1,[1,5,8]) = fTmp([1,5,8])
        enddo
    else
        stop 'no define BC'
    endif

    if(xMaxBC==DirecletUP)then
        do  y = 1, yDim
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(xDim), yGrid(y), vel)
            elseif(VelocityKind==1) then
                call evaluateParabolicVelocity(xGrid(xDim), yGrid(y), vel)
            elseif(VelocityKind==2) then
                vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
                vel(2)=uIn(2)
            endif
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim) = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)
            fIn(y,xDim,0:lbmDim)   = fEq(0:lbmDim)
        enddo
    elseif(xMaxBC==DirecletUU)then
        do  y = 1, yDim
            if(VelocityKind==0) then
                call evaluateShearVelocity(xGrid(xDim), yGrid(y), vel)
            elseif(VelocityKind==1) then
                call evaluateParabolicVelocity(xGrid(xDim), yGrid(y), vel)
            elseif(VelocityKind==2) then
                vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
                vel(2)=uIn(2)
            endif        
            call evaluateShearVelocity(xGrid(xDim), yGrid(y), vel)
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim) = wt(0:lbmDim) * den(y,xDim-1) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            uSqri           = sum(uuu(y,xDim-1,1:2)**2)
            uxyi(0:lbmDim) = uuu(y,xDim-1,1) * ee(0:lbmDim,1) + uuu(y,xDim-1,2) * ee(0:lbmDim,2)
            fEqi(0:lbmDim)  = wt(0:lbmDim) * den(y,xDim-1) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            fIn(y,xDim,[3,6,7]) =fEq([3,6,7])+ (fIn(y,xDim-1,[3,6,7])-fEqi([3,6,7]))
        enddo
    elseif(xMaxBC==Advection1)then
        fIn(:,xDim,[3,6,7]) = fIn(:,xDim-1,[3,6,7])
    elseif(xMaxBC==Advection2)then
        fIn(:,xDim,[3,6,7]) = 2.0*fIn(:,xDim-1,[3,6,7])-fIn(:,xDim-2,[3,6,7])
    elseif(xMaxBC==Periodic .or. xMaxBC==fluid)then
        ! no need set
    elseif(xMaxBC==wall) then
        do  y = 1, yDim
            fTmp([3,6,7]) = fIn(y,xDim,oppo([3,6,7]))
            fIn(y,xDim,[3,6,7]) = fTmp([3,6,7])
        enddo
    else
        stop 'no define BC'
    endif

    !=======y-direction
    if(yMinBC==DirecletUP)then
        do  x = 1, xDim
            call evaluateShearVelocity(xGrid(x), yGrid(1), vel)
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)
            fIn(1,x,0:lbmDim)   = fEq(0:lbmDim)
        enddo
    elseif(yMinBC==DirecletUU)then
        do  x = 1, xDim
            call evaluateShearVelocity(xGrid(x), yGrid(1), vel)
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim)  = wt(0:lbmDim) * den(2,x) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            uSqri           = sum(uuu(2,x,1:2)**2)
            uxyi(0:lbmDim) = uuu(2,x,1) * ee(0:lbmDim,1) + uuu(2,x,2) * ee(0:lbmDim,2)
            fEqi(0:lbmDim)  = wt(0:lbmDim) * den(2,x) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            fIn(1,x,[2,5,6]) =fEq([2,5,6])+ (fIn(2,x,[2,5,6])-fEqi([2,5,6]))
        enddo
    elseif(yMinBC==Advection1)then
        fIn(1,:,[2,5,6]) = fIn(2,:,[2,5,6])
    elseif(yMinBC==Advection2)then
        fIn(1,:,[2,5,6]) = 2.0*fIn(2,:,[2,5,6])-fIn(3,:,[2,5,6])
    elseif(yMinBC==Periodic .or. yMinBC==fluid)then
        ! no need set
    elseif(yMinBC==wall) then
        do  x = 1, xDim
            fTmp([2,5,6]) = fIn(1,x,oppo([2,5,6]))
            fIn(1,x,[2,5,6]) = fTmp([2,5,6])
        enddo
    ! xianguang Luo, July,2022
    elseif(yMinBC==MovingWall) then
        do  x = 1, xDim
            if(MovingKind1==0) then
                if(VelocityKind==0) then
                    call evaluateShearVelocity(xGrid(x), yGrid(1), vel)
                elseif(VelocityKind==1) then
                    call evaluateParabolicVelocity(xGrid(x), yGrid(1), vel)
                elseif(VelocityKind==2) then
                    vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
                    vel(2)=uIn(2)
                endif
            elseif(MovingKind1==1) then
                vel(1)=MovingVel1
                vel(2)=0.0d0
            elseif(MovingKind1==2) then
                vel(1)=MovingVel1 * dsin(2*pi*MovingFreq1*time/Tref)
                vel(2)=0.0d0               
            endif
            fTmp(2) = fIn(1,x,oppo(2)) + 2.0*wt(2)*denIn*(ee(2,1)*vel(1)+ee(2,2)*vel(2))*3.0
            fTmp(5) = fIn(1,x,oppo(5)) + 2.0*wt(5)*denIn*(ee(5,1)*vel(1)+ee(5,2)*vel(2))*3.0
            fTmp(6) = fIn(1,x,oppo(6)) + 2.0*wt(6)*denIn*(ee(6,1)*vel(1)+ee(6,2)*vel(2))*3.0
            fIn(1,x,[2,5,6]) = fTmp([2,5,6])
        enddo
    else
        stop 'no define BC'
    endif

    if(yMaxBC==DirecletUP)then
        do  x = 1, xDim
            call evaluateShearVelocity(xGrid(x), yGrid(yDim), vel)
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)
            fIn(yDim,x,0:lbmDim)   = fEq(0:lbmDim)
        enddo
    elseif(yMaxBC==DirecletUU)then
        do  x = 1, xDim
            call evaluateShearVelocity(xGrid(x), yGrid(yDim), vel)
            uSqr   = vel(1)**2 + vel(2)**2
            uxy(0:lbmDim) = vel(1) * ee(0:lbmDim,1) + vel(2) * ee(0:lbmDim,2)
            fEq(0:lbmDim)  = wt(0:lbmDim) * den(yDim-1,x) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            uSqri           = sum(uuu(yDim-1,x,1:2)**2)
            uxyi(0:lbmDim) = uuu(yDim-1,x,1) * ee(0:lbmDim,1) + uuu(yDim-1,x,2) * ee(0:lbmDim,2)
            fEqi(0:lbmDim)  = wt(0:lbmDim) * den(yDim-1,x) * (1.0d0 + 3.0d0 * uxy(0:lbmDim) + 4.5d0 * uxy(0:lbmDim) * uxy(0:lbmDim) - 1.5d0 * uSqr)

            fIn(yDim,x,[4,7,8]) =fEq([4,7,8])+ (fIn(yDim-1,x,[4,7,8])-fEqi([4,7,8]))
        enddo
    elseif(yMaxBC==Advection1)then
        fIn(yDim,:,[4,7,8]) = fIn(yDim-1,:,[4,7,8])
    elseif(yMaxBC==Advection2)then
        fIn(yDim,:,[4,7,8]) = 2.0*fIn(yDim-1,:,[4,7,8])-fIn(yDim-2,:,[4,7,8])
    elseif(yMaxBC==Periodic .or. yMaxBC==fluid)then
        ! no need set
    elseif(yMaxBC==wall) then
        do  x = 1, xDim
            fTmp([4,7,8]) = fIn(yDim,x,oppo([4,7,8]))
            fIn(yDim,x,[4,7,8]) = fTmp([4,7,8])
        enddo
    ! xianguang Luo, July,2022
    elseif(yMaxBC==MovingWall) then
        do  x = 1, xDim
        call evaluateShearVelocity(xGrid(x), yGrid(yDim), vel)
            if(MovingKind2==0) then
                if(VelocityKind==0) then
                    call evaluateShearVelocity(xGrid(x), yGrid(yDim), vel)
                elseif(VelocityKind==1) then
                    call evaluateParabolicVelocity(xGrid(x), yGrid(yDim), vel)
                elseif(VelocityKind==2) then
                    vel(1)=uIn(1) + VelocityAmp * dsin(2*pi*VelocityFreq*time/Tref)
                    vel(2)=uIn(2)
                endif
            elseif(MovingKind2==1) then
                vel(1)=MovingVel2
                vel(2)=0.0d0
            elseif(MovingKind2==2) then
                vel(1)=MovingVel2 * dsin(2*pi*MovingFreq2*time/Tref)
                vel(2)=0.0d0
            endif
            fTmp(4) = fIn(yDim,x,oppo(4)) + 2*wt(4)*denIn*(ee(4,1)*vel(1)+ee(4,2)*vel(2))*3.0
            fTmp(7) = fIn(yDim,x,oppo(7)) + 2*wt(7)*denIn*(ee(7,1)*vel(1)+ee(7,2)*vel(2))*3.0
            fTmp(8) = fIn(yDim,x,oppo(8)) + 2*wt(8)*denIn*(ee(8,1)*vel(1)+ee(8,2)*vel(2))*3.0
            fIn(yDim,x,[4,7,8]) = fTmp([4,7,8])
        enddo
    else
        stop 'no define BC'
    endif
!   ==========================
    END SUBROUTINE
