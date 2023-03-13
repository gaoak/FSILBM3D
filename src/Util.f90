!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!     
!        
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wrtInfoTitl()
    USE simParam
    implicit none
    integer::   i,iFish
    integer,parameter::nameLen=4
    character (LEN=nameLen):: fileName,Nodename

    !open(111,file='./DatInfo/ForceStress.plt')
    !write(111,*)'variables= "t"  "Fx"  "Fy"  '
    !close(111)
        
    do iFish=1,nFish

        write(fileName,'(I4)') iFish
        fileName = adjustr(fileName)
        do  i=1,nameLen
             if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo
        
        open(111,file='./DatInfo/ForceDirect_'//trim(filename)//'.plt')
        write(111,*)'variables= "t"  "Fx"  "Fy"  '
        close(111)

        open(111,file='./DatInfo/SampBodyNodeBegin_'//trim(filename)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "u"  "v"  "ax"  "ay" '
        close(111)

        write(Nodename,'(I4.4)') nNd(iFish)
        open(111,file='./DatInfo/SampBodyNodeEnd_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "u"  "v"  "ax"  "ay" '
        close(111)

        open(111,file='./DatInfo/SampBodyNodeCenter_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "u"  "v"  "ax"  "ay" '
        close(111)
        open(111,file='./DatInfo/SampBodyMean_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "u"  "v"  "ax"  "ay" '
        close(111)

        do  i=1,numSampBody
        write(Nodename,'(I4.4)') SampBodyNode(iFish,i)
        open(111,file='./DatInfo/SampBodyNode'//trim(Nodename)//'_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "u"  "v"  "ax"  "ay" '
        close(111)
        enddo

        do  i=1,numSampFlow
        write(Nodename,'(I4.4)') i
        open(111,file='./DatInfo/SampFlowPint'//trim(Nodename)//'.plt')
        write(111,*)'variables= "t"  "p" "u"  "v"  '
        close(111)
        enddo

        open(111,file='./DatInfo/SampBodyAngular_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "AoA"  "Ty-Hy"  "Hy"  "Ty"'
        close(111)
        open(111,file='./DatInfo/Power_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t" "Ptot" "Paero" "Piner" "Pax" "Pay"  "Pix" "Piy" '
        close(111)
        open(111,file='./DatInfo/Converg.plt')
        write(111,*)'variables= "t"  "Convergence"  '
        close(111)
        open(111,file='./DatInfo/MaMax.plt')
        write(111,*)'variables= "t"  "MaMax"  '
        close(111)
        open(111,file='./DatInfo/MaxValBody_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "xyzMax" "velMax" "accMax"  '
        close(111)

        open(111,file='./DatInfo/Length_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "Length"  '
        close(111)

        open(111,file='./DatInfo/Energy_'//trim(fileName)//'.plt')
        write(111,*)'variables= "t","Es","Eb","Ep","Ek","Ew","Et"'
        close(111)
        
    enddo
    END SUBROUTINE
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!     
!        
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wrtInfo()
    USE simParam
    implicit none
    integer:: iEL,i,z,y,x,zbgn,ybgn,xbgn,zend,yend,xend,iFish
    real(8):: EEE(2),strainEnergy(nEL_max,2)
    real(8):: convergence,MaMax,Length,weightm,velocity(1:3),Pressure
    real(8):: Ptot,Paero,Piner,Pax,Pay,Pix,Piy
    integer,parameter::nameLen=4
    character (LEN=nameLen):: fileName,Nodename
      
    do iFish=1,nFish
        write(fileName,'(I4)') iFish
        fileName = adjustr(fileName)
        do  i=1,nameLen
             if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo  
        !open(111,file='./DatInfo/ForceStress.plt',position='append')
        !write(111,'(4E20.10)')time/Tref,sum(extful2(1:nND,1:2),1)/Fref                
        !close(111)
    !===============================================================================    
        open(111,file='./DatInfo/ForceDirect_'//trim(fileName)//'.plt',position='append')
        write(111,'(4E20.10)')time/Tref,sum(extful(iFish,1:nND(iFish),1:2),1)/Fref  
        close(111)
    !===============================================================================
        open(111,file='./DatInfo/SampBodyNodeBegin_'//trim(fileName)//'.plt',position='append')
        write(111,'(7E20.10)')time/Tref,xyzful(iFish,1,1:2)/Lref,velful(iFish,1,1:2)/Uref,accful(iFish,1,1:2)/Aref                     
        close(111)
        !===
        write(Nodename,'(I4.4)') nNd(iFish)
        open(111,file='./DatInfo/SampBodyNodeEnd_'//trim(fileName)//'.plt',position='append')
        write(111,'(7E20.10)')time/Tref,xyzful(iFish,nND(iFish),1:2)/Lref,velful(iFish,nND(iFish),1:2)/Uref,accful(iFish,nND(iFish),1:2)/Aref                  
        close(111)
      
    !===============================================================================    
        open(111,file='./DatInfo/SampBodyNodeCenter_'//trim(fileName)//'.plt',position='append')
        write(111,'(7E20.10)')time/Tref,xyzful(iFish,(nND(iFish)+1)/2,1:2)/Lref,velful(iFish,(nND(iFish)+1)/2,1:2)/Uref,accful(iFish,(nND(iFish)+1)/2,1:2)/Aref                     
        close(111)

        open(111,file='./DatInfo/SampBodyMean_'//trim(fileName)//'.plt',position='append')
        write(111,'(7E20.10)')time/Tref,sum(xyzful(iFish,1:nND(iFish),1:2)*mssful(iFish,1:nND(iFish),1:2),1)/sum(mssful(iFish,1:nND(iFish),1:2),1)/Lref, &
                                        sum(velful(iFish,1:nND(iFish),1:2)*mssful(iFish,1:nND(iFish),1:2),1)/sum(mssful(iFish,1:nND(iFish),1:2),1)/Uref, &
                                        sum(accful(iFish,1:nND(iFish),1:2)*mssful(iFish,1:nND(iFish),1:2),1)/sum(mssful(iFish,1:nND(iFish),1:2),1)/Aref                        
        close(111)

    !===============================================================================
        do  i=1,numSampBody 
           write(Nodename,'(I4.4)') SampBodyNode(iFish,i)
           open(111,file='./DatInfo/SampBodyNode'//trim(Nodename)//'_'//trim(fileName)//'.plt',position='append')
           write(111,'(10E20.10)')time/Tref, xyzful(iFish,SampBodyNode(iFish,i),1:2)/Lref,velful(iFish,SampBodyNode(iFish,i),1:2)/Uref,accful(iFish,SampBodyNode(iFish,i),1:2)/Aref 
           close(111)
        enddo
    !===============================================================================
        open(111,file='./DatInfo/SampBodyAngular_'//trim(fileName)//'.plt',position='append')
           write(111,'(5E20.10)')time/Tref,datan((xyzful(iFish,nND(iFish),2)-xyzful(iFish,1,2))/(xyzful(iFish,nND(iFish),1)-xyzful(iFish,1,1))),    &
                                        xyzful(iFish,nND(iFish),2)/Lref-xyzful(iFish,1,2)/Lref,xyzful(iFish,1,2)/Lref,xyzful(iFish,nND(iFish),2)/Lref
        close(111)
    !===============================================================================
        Pax=sum(extful(iFish,1:nND(iFish),1)*velful(iFish,1:nND(iFish),1))/Pref
        Pay=sum(extful(iFish,1:nND(iFish),2)*velful(iFish,1:nND(iFish),2))/Pref
        Paero=Pax+Pay
        Pix=-sum(mssful(iFish,1:nND(iFish),1)*accful(iFish,1:nND(iFish),1)*velful(iFish,1:nND(iFish),1))/Pref
        Piy=-sum(mssful(iFish,1:nND(iFish),2)*accful(iFish,1:nND(iFish),2)*velful(iFish,1:nND(iFish),2))/Pref
    
        Piner=Pix+Piy
        Ptot=Paero+Piner
        open(111,file='./DatInfo/Power_'//trim(fileName)//'.plt',position='append')
           write(111,'(8E20.10)')time/Tref,Ptot,Paero,Piner,Pax,Pay,Pix,Piy
                                                                        
        close(111)
    !===============================================================================
        open(111,file='./DatInfo/MaxValBody_'//trim(fileName)//'.plt',position='append')
        write(111,'(4E20.10)')time/Tref,maxval(dsqrt(xyzful(iFish,1:nND(iFish),1)**2+xyzful(iFish,1:nND(iFish),2)**2))/Lref, &
                                        maxval(dsqrt(velful(iFish,1:nND(iFish),1)**2+velful(iFish,1:nND(iFish),2)**2))/Uref, &
                                        maxval(dsqrt(accful(iFish,1:nND(iFish),1)**2+accful(iFish,1:nND(iFish),2)**2))/Aref
        close(111)
    !===============================================================================
        Length=0.0d0
        do    iEL=1,nEL(iFish)
            Length=Length+dsqrt(sum((xyzful(iFish,ele(iFish,iEL,1),1:2)-xyzful(iFish,ele(iFish,iEL,2),1:2))**2))
        enddo
        open(111,file='./DatInfo/Length_'//trim(fileName)//'.plt',position='append')
        write(111,'(2E20.10)')time/Tref, Length
        close(111)
    
    !===============================================================================
        call strain_energy_D(    strainEnergy(1:nEL(iFish),1:2),xyzful0(iFish,1:nND(iFish),1),xyzful0(iFish,1:nND(iFish),2),xyzful0(iFish,1:nND(iFish),3), &
                                               xyzful(iFish,1:nND(iFish),1), xyzful(iFish,1:nND(iFish),2), xyzful(iFish,1:nND(iFish),3), &
                                                  ele(iFish,1:nEL(iFish),1:5), prop(iFish,1:nMT(iFish),1:10), &
                                 triad_n1(iFish,1:3,1:3,1:nEL(iFish)),triad_n2(iFish,1:3,1:3,1:nEL(iFish)),triad_n3(iFish,1:3,1:3,1:nEL(iFish)), &
                                 triad_ee(iFish,1:3,1:3,1:nEL(iFish)),triad_e0(iFish,1:3,1:3,1:nEL(iFish)),triad_nn(iFish,1:3,1:3,1:nND(iFish)), &
                                 nND(iFish),nEL(iFish),nMT(iFish) &
                            )
        EEE(1)=sum(strainEnergy(1:nEL(iFish),1))
        EEE(2)=sum(strainEnergy(1:nEL(iFish),2))
        Es=EEE(1)/Eref
        Eb=EEE(2)/Eref
        Ep=Es+Eb
        Ew=Ew+Paero*timeOutInfo
        Ek=0.5*sum(mssful(iFish,1:nND(iFish),1:6)*velful(iFish,1:nND(iFish),1:6)*velful(iFish,1:nND(iFish),1:6))/Eref

        !write(*,*)'Ekt:', 0.5*sum(mssful(1:nND,1:3)*velful(1:nND,1:3)*velful(1:nND,1:3))/Eref
        !write(*,*)'Ekr:', 0.5*sum(mssful(1:nND,4:6)*velful(1:nND,4:6)*velful(1:nND,4:6))/Eref

        Et=Ek+Ep
        open(111,file='./DatInfo/Energy_'//trim(fileName)//'.plt',position='append')
        write(111,'(7E20.10)')time/Tref,Es,Eb,Ep,Ek,Ew,Et
        close(111)
        
    enddo !iFish       
        
    !Write data of monitoring points of fluid
        !===============================================================================
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


        velocity(1:2)=0.0d0
        Pressure=0.0d0

        do    x=xbgn,xend
        do    y=ybgn,yend
            weightm=(dx(x)-dabs(xGrid(x)-SampFlowPint(i,1)))*(dy(y)-dabs(yGrid(y)-SampFlowPint(i,2)))
            velocity(1:2)=velocity(1:2)+uuu(y,x,1:2)*weightm/(dx(x)*dy(y))
            Pressure=Pressure+prs(y,x)*weightm/(dx(x)*dy(y))
        enddo
        enddo

        write(fileName,'(I4.4)') i
        open(111,file='./DatInfo/SampFlowPint'//trim(fileName)//'.plt',position='append')
        write(111,'(4E20.10)') time/Tref,  Pressure/(0.5*denIn*Uref**2), velocity(1:2)/Uref
        close(111)
        enddo
        !============================================================================================== 
        MaMax=MaxVal(dsqrt( uuu(1:yDim,1:xDim,1)**2+uuu(1:yDim,1:xDim,2)**2))/dsqrt(Cs2)
        open (111,file='./DatInfo/MaMax.plt',position='append')
        write(111,'(2E20.10)')time/Tref,MaMax             
        close(111)
        !============================================================================================== 
        UNow=sum(dsqrt( uuu(1:yDim,1:xDim,1)**2+uuu(1:yDim,1:xDim,2)**2))          
        convergence=dabs(UNow-UPre)/UNow
        UPre=UNow
        open (111,file='./DatInfo/Converg.plt',position='append')
        write(111,'(2E20.10)')time/Tref,convergence             
        close(111)
        !============================================================================================== 
    
     
    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!     
!        
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE movGrid(dim,direction)
    USE simParam
    implicit none
    integer:: dim,direction
    
    if    (dim==1)then      !x
        if(direction<0)       then      !move to -x
            fIn(1:yDim,2:xDim,0:LBMDim)=fIn(1:yDim,1:xDim-1,0:LBMDim)
            xGrid(1:xDim)=xGrid(1:xDim)-dx(1:xDim)
        elseif(direction>0)    then     !move to +x
            fIn(1:yDim,1:xDim-1,0:LBMDim)=fIn(1:yDim,2:xDim,0:LBMDim)
            xGrid(1:xDim)=xGrid(1:xDim)+dx(1:xDim)
        else
        endif
    elseif(dim==2)then      !y
        if(direction<0)        then     !move to -y
            fIn(2:yDim,1:xDim,0:LBMDim)=fIn(1:yDim-1,1:xDim,0:LBMDim)
            yGrid(1:yDim)=yGrid(:)-dy(:)
        elseif(direction>0)    then     !move to +y
            fIn(1:yDim-1,1:xDim,0:LBMDim)=fIn(2:yDim,1:xDim,0:LBMDim)
            yGrid(1:yDim)=yGrid(1:yDim)+dy(1:yDim)
        else
        endif
    else
    endif
    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!            
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  cptMove(move,xB,xG,d)  !Judge moving direction of moving grid
    implicit none
    integer:: move(2),i
    real(8):: xB(2),xG(2),d(2)
    do  i=1,2
        if    (xB(i)-xG(i)> d(i))then
            move(i)=1
        elseif(xB(i)-xG(i)<-d(i))then
            move(i)=-1
        else
            move(i)=0
        endif
    enddo
    end subroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!            
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cptIref(NDref,IXref,IYref,nND,xDim,yDim,xyzful,xGrid,yGrid,Xref,Yref)
    implicit none
    integer:: NDref,IXref,IYref,nND,xDim,yDim,iND,i
    real(8):: distance,xyzful(1:nND,1:2),xGrid(xDim),yGrid(yDim),Xref,Yref

    NDref=minloc(dsqrt((xyzful(1:nND,1)-Xref)**2+(xyzful(1:nND,2)-Yref)**2),1)                         
    IXref=minloc(dabs(xyzful(NDref,1)-xGrid(1:xDim)),1)
    IYref=minloc(dabs(xyzful(NDref,2)-yGrid(1:yDim)),1)
    end subroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!   Coordinate transformation of rotation       
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine AoAtoTTT(AoA,TTT)
    real(8)::AoA(3),TTT(3,3),rrx(3,3),rry(3,3),rrz(3,3)
    TTT(:,:)=0.0d0
    TTT(1,1)=1.0d0
    TTT(2,2)=1.0d0
    TTT(3,3)=1.0d0
    rrx(1:3,1:3)=reshape([  1.0d0,0.0d0,0.0d0,  &
                            0.0d0,dcos(AoA(1)),dsin(AoA(1)), &
                            0.0d0,-dsin(AoA(1)),dcos(AoA(1))],[3,3])

    rry(1:3,1:3)=reshape([  dcos(AoA(2)),0.0d0,-dsin(AoA(2)),  &
                            0.0d0       ,1.0d0,       0.0d0, &
                           dsin(AoA(2)),0.0d0,dcos(AoA(2))],[3,3])

    rrz(1:3,1:3)=reshape([ dcos(AoA(3)),dsin(AoA(3)),0.0d0, &
                          -dsin(AoA(3)),dcos(AoA(3)),0.0d0, &
                           0.0d0       ,0.0d0       ,1.0d0],[3,3])
    TTT=matmul(rrz,TTT)
    TTT=matmul(rry,TTT)
    TTT=matmul(rrx,TTT)

    end subroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write  string to unit=idfile  in binary format
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
!   Sinusoidal velocity noise of feed flow field
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  initDisturb(xDim,yDim,xGrid,yGrid,AmplInitDist,FreqInitDist,uuu)
    implicit none
    real(8), parameter:: Pi=3.141592653589793d0
    integer:: xDim, yDim
    real(8):: AmplInitDist(1:2),FreqInitDist(1:2), uuu(yDim,xDim,2),xGrid(1:xDim),yGrid(1:yDim)

    real(8):: r
    integer:: x, y,z

    do  y = 1, yDim
    do  x = 1, xDim
        !call random_number(r)
        !uuu(z,y,x,1)=uuu(z,y,x,1)+(r-0.5d0)*AmplInitDist(1)
        !call random_number(r)
        !uuu(z,y,x,2)=uuu(z,y,x,2)+(r-0.5d0)*AmplInitDist(2)
        uuu(y,x,1)=uuu(y,x,1)+AmplInitDist(1)*dsin(2.0*pi*FreqInitDist(1)*xGrid(x))*dsin(2.0*pi*FreqInitDist(2)*yGrid(y))
        uuu(y,x,2)=uuu(y,x,2)+AmplInitDist(2)*dsin(2.0*pi*FreqInitDist(1)*xGrid(x))*dsin(2.0*pi*FreqInitDist(2)*yGrid(y))
    enddo
    enddo


    endsubroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Shear flow velocity boundary
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE evaluateShearVelocity(x, y, vel)
        USE simParam
        real(8):: x, y, vel(1:SpcDim)
        vel(1) = uIn(1) + y * shearRateIn(1)
        vel(2) = uIn(2) - x * shearRateIn(2)
    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   lxguang 2022.11 Add Parabolic Velocity Inlet Boundary
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE evaluateParabolicVelocity(x, y, vel)
        USE simParam
        real(8):: x, y, vel(1:SpcDim)
        vel(1) = uIn(1) - ParabolicParameter(1) * (y - yGrid(1)) * (y -  yGrid(yDim)) -uIn(1)
        vel(2) = uIn(2) - ParabolicParameter(2) * (x - xGrid(1)) * (x -  xGrid(yDim)) -uIn(2)
    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   lxguang 2022.09 Add flapping board intermittent movement mode
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE setHeadMotion(ifish,ttt,deltatt,iisubstep,subdeltatt)
        USE simParam
        implicit none
        integer:: ifish
        real(8):: tTol, judgetime, ttt, sumt, deltatt, iisubstep, subdeltatt
        tTol = 1.0 + sum(waittingTime(ifish,:))
        judgetime = (ttt*Freq/tTol-floor(ttt*Freq/tTol))*tTol
        !Judge the position of the current time in the Period
        if (judgetime <= 0.5)then
            sumt = judgetime/Freq-deltatt+iisubstep*subdeltatt
            XYZ(ifish,1:3) = XYZo(ifish,1:3)+XYZAmpl(ifish,1:3)*dcos(2.0*pi*Freq*sumt+XYZPhi(ifish,1:3))
            UVW(ifish,1:3) = -2.0*pi*Freq*XYZAmpl(ifish,1:3)*dsin(2.0*pi*Freq*sumt+XYZPhi(ifish,1:3))
        elseif (judgetime <= waittingTime(ifish,1) + 0.5)then
            XYZ(ifish,1:3) = XYZo(ifish,1:3)-XYZAmpl(ifish,1:3)
            UVW(ifish,1:3) = 0.0
        elseif (judgetime <= waittingTime(ifish,1) + 1.0)then
            sumt = (judgetime-waittingTime(ifish,1))/Freq-deltatt+iisubstep*subdeltatt
            XYZ(ifish,1:3) = XYZo(ifish,1:3)+XYZAmpl(ifish,1:3)*dcos(2.0*pi*Freq*sumt+XYZPhi(ifish,1:3))
            UVW(ifish,1:3) = -2.0*pi*Freq*XYZAmpl(ifish,1:3)*dsin(2.0*pi*Freq*sumt+XYZPhi(ifish,1:3))
        else
            XYZ(ifish,1:3) = XYZo(ifish,1:3)+XYZAmpl(ifish,1:3)
            UVW(ifish,1:3) = 0.0
        endif
    END SUBROUTINE