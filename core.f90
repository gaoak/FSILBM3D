!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	计算宏观量：密度和速度
!	copyright@ RuNanHua 
!	版权所有，华如南（中国科大近代力学系）  
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE cptMacr()
    USE simParam
    implicit none
	integer::x,y,z,ig
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z) 
    do  x = 1, xDim
    do  y = 1, yDim
	do  z = 1, zDim
        den(z,y,x  )  = SUM(fIn(z,y,x,0:lbmDim)) 
        uuu(z,y,x,1)  = SUM(fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,1))/den(z,y,x)
        uuu(z,y,x,2)  = SUM(fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,2))/den(z,y,x)
	    uuu(z,y,x,3)  = SUM(fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,3))/den(z,y,x)	
        prs(z,y,x)   = Cs2*(den(z,y,x)-denIn)
	enddo
    enddo
    enddo
    !$OMP END PARALLEL DO
    !prs(1:zDim,1:yDim,1:xDim)   = Cs2*(den(1:zDim,1:yDim,1:xDim)-denIn)
    END SUBROUTINE

!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	collision step SRT or MRT
!   碰撞模型：单松弛或多松弛
!	copyright@ RuNanHua 
!	版权所有，华如南（中国科大近代力学系）
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE collide()
    USE simParam
    implicit none 
    real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim),m(0:lbmDim),mEq(0:lbmDim),Flb(0:lbmDim)
    integer:: x,y,z

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,uSqr,uxyz,fEq,Flb)
    do	x = 1, xDim
    do	y = 1, yDim
    do	z = 1, zDim
        uSqr           = sum(uuu(z,y,x,1:3)**2)
		uxyz(0:lbmDim) = uuu(z,y,x,1) * ee(0:lbmDim,1) + uuu(z,y,x,2) * ee(0:lbmDim,2)+uuu(z,y,x,3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * den(z,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        Flb(0:lbmDim)  = dt*wt(0:lbmDim)*( &
                          (3.0*(ee(0:lbmDim,1)-uuu(z,y,x,1))+9.0*(ee(0:lbmDim,1)*uuu(z,y,x,1)+ee(0:lbmDim,2)*uuu(z,y,x,2)+ee(0:lbmDim,3)*uuu(z,y,x,3))*ee(0:lbmDim,1))*force(z,y,x,1) &
						 +(3.0*(ee(0:lbmDim,2)-uuu(z,y,x,2))+9.0*(ee(0:lbmDim,1)*uuu(z,y,x,1)+ee(0:lbmDim,2)*uuu(z,y,x,2)+ee(0:lbmDim,3)*uuu(z,y,x,3))*ee(0:lbmDim,2))*force(z,y,x,2) &
                         +(3.0*(ee(0:lbmDim,3)-uuu(z,y,x,3))+9.0*(ee(0:lbmDim,1)*uuu(z,y,x,1)+ee(0:lbmDim,2)*uuu(z,y,x,2)+ee(0:lbmDim,3)*uuu(z,y,x,3))*ee(0:lbmDim,3))*force(z,y,x,3) &
						                 )

        if    (iCollidModel==1)then
        !本空间单松弛碰撞
		fIn(z,y,x,0:lbmDim) = fIn(z,y,x,0:lbmDim) + Omega * (fEq(0:lbmDim)-fIn(z,y,x,0:lbmDim))+(1.0-0.5*Omega)*Flb(0:lbmDim)
        elseif(iCollidModel==2)then
        !矩空间多松弛碰撞
        fIn(z,y,x,0:lbmDim)=fIn(z,y,x,0:lbmDim)+MATMUL( M_COLLID(0:lbmDim,0:lbmDim), fEq(0:lbmDim)-fIn(z,y,x,0:lbmDim) ) &
                                               +MATMUL( M_FORCE(0:lbmDim,0:lbmDim),Flb(0:lbmDim))
        else
            write(*,*)' Collide Model is not defined'
        endif
    enddo
    enddo
    enddo
    !$OMP END PARALLEL DO
    END SUBROUTINE

!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	对流模型：均匀网格的正常对流和非均匀网格的插值对流 
!	copyright@ RuNanHua 
!	版权所有，华如南（中国科大近代力学系）   
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE streams()
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
    if    (iStreamModel==1) then   !对流
!============================================
        do  i=0,lbmDim
        do  k=1,3
            if(strmDir(i,k)/=0)then
                fIn(1:zDim,1:yDim,1:xDim, i)=cshift(fIn(1:zDim,1:yDim,1:xDim, i),shift=strmDir(i,k),dim=k)
            endif
        enddo
        enddo
!============================================
    elseif(iStreamModel==2)then    !插值
!============================================
         
        do  i=1,lbmDim
        fInTemp(1:zDim,1:yDim,1:xDim)=fIn(1:zDim,1:yDim,1:xDim,i)
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,upwind,center,outer,upxc0,upxcm,upxcmm,upyc0,upycm,upycmm,upzc0,upzcm,upzcmm,cnxc0,cnxcm,cnxcp,cnyc0,cnycm,cnycp,cnzc0,cnzcm,cnzcp)
        do  x=2,xDim-1
        do  y=2,yDim-1
        do  z=2,zDim-1
            !设置逻辑判断
            upwind = (x>=3).and.(x<=xDim-2).and.(y>=3).and.(y<=yDim-2) .and. (z>=3) .and. (z<=zDim-2)                    
            center = (x>=2).and.(x<=xDim-1).and.(y>=2).and.(y<=yDim-1) .and. (z>=2) .and. (z<=zDim-1) .and. (.not.upwind) 
            outer  = .not.(upwind .or. center)

        !******************************************************************
        !******************************************************************
        if(upwind) then      
            !二阶迎风插值
            !二阶迎风插值系数
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
            !二阶中心插值
            !二阶中心插值系数
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

!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   far field boundary set	 
!   设置远场边界条件，全部设置为平衡函数
!	copyright@ RuNanHua 
!	版权所有，华如南（中国科大近代力学系） 
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE setBund()
    USE simParam
    implicit none
    integer:: i,k
    integer:: x, y, z,ig
    logical:: upwind,center,outer
    real(8):: uSqr ,uxyz(0:lbmDim) ,fEq(0:lbmDim)
    uSqr           = sum(uuuIn(1:3)**2)
    uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
    fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,upwind,center,outer)
    do  x=1,xDim 
    do  y=1,yDim 
    do  z=1,zDim 
        !设置逻辑判断
        upwind = (x>=3).and.(x<=xDim-2).and.(y>=3).and.(y<=yDim-2) .and. (z>=3) .and. (z<=zDim-2)                    
        center = (x>=2).and.(x<=xDim-1).and.(y>=2).and.(y<=yDim-1) .and. (z>=2) .and. (z<=zDim-1) .and. (.not.upwind) 
        outer  = .not.(upwind .or. center)
        if(iBC==1)then
            if	  ( outer  ) then               !最外一层            
            fIn(z,y,x,0:lbmDim)= fEq(0:lbmDim)
            endif
        elseif(iBC==2)then
            if	  ( outer .or. center ) then    !最外两层
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

!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   far field boundary set	 
!   设置远场边界条件，各种条件
!	copyright@ RuNanHua 
!	版权所有，华如南（中国科大近代力学系）
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE setBund2()
    USE simParam
    implicit none
    integer:: i,k
    integer:: x, y, z,ig
    real(8):: uSqr ,uxyz(0:lbmDim) ,fEq(0:lbmDim)
    real(8):: uSqri,uxyzi(0:lbmDim),fEqi(0:lbmDim)
    real(8):: fTmp(0:lbmDim)   
!	---------------------------------------
!	---------------------------------------

    !=======x-direction
    if(xMinBC==DirecletUP)then
        uSqr           = sum(uuuIn(1:3)**2)
		uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        do  y = 1, yDim
	    do  z = 1, zDim
		fIn(z,y,1,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(xMinBC==DirecletUU)then
        do  y = 1, yDim
	    do  z = 1, zDim
        uSqr           = sum(uuuIn(1:3)**2)
		uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
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
        uSqr           = sum(uuuIn(1:3)**2)
		uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
        fEq(0:lbmDim)  = wt(0:lbmDim) * denIn * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
        do  y = 1, yDim
	    do  z = 1, zDim
		fIn(z,y,1,0:lbmDim)   = fEq(0:lbmDim)
        enddo
        enddo
    elseif(xMaxBC==DirecletUU)then
        do  y = 1, yDim
	    do  z = 1, zDim
        uSqr           = sum(uuuIn(1:3)**2)
		uxyz(0:lbmDim) = uuuIn(1) * ee(0:lbmDim,1) + uuuIn(2) * ee(0:lbmDim,2)+uuuIn(3) * ee(0:lbmDim,3)
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
    else
        stop 'no define BC'
    endif
!   ==========================
    END SUBROUTINE