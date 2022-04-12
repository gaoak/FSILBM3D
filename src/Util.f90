!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    write files' header
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wrtInfoTitl()
    USE simParam
    implicit none
    integer::   i
    integer,parameter::nameLen=4
    character (LEN=nameLen):: fileName

        open(111,file='./DatInfo/ForceStress.plt')
        write(111,*)'variables= "t"  "Fx"  "Fy"  "Fz"'
        close(111)
        open(111,file='./DatInfo/ForceDirect.plt')
        write(111,*)'variables= "t"  "Fx"  "Fy"  "Fz"'
        close(111)

    do  i=1,5
        write(fileName,'(I4.4)') NDtl(i)
        open(111,file='./DatInfo/SampBodyNodeB'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)
    enddo 
    do  i=1,3   
        write(fileName,'(I4.4)') NDhd(i)
        open(111,file='./DatInfo/SampBodyNodeB'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)
    enddo

        open(111,file='./DatInfo/SampBodyCentP.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)
        open(111,file='./DatInfo/SampBodyCentM.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)

        open(111,file='./DatInfo/SampBodyAngular1.plt')
        write(111,*)'variables= "t"  "AoA"  "Ty-Hy"  "Hy"  "Ty"'
        close(111)
        open(111,file='./DatInfo/SampBodyAngular2.plt')
        write(111,*)'variables= "t"  "AoA"  "Ty-Hy"  "Hy"  "Ty"'
        close(111)
        open(111,file='./DatInfo/SampBodyAngular3.plt')
        write(111,*)'variables= "t"  "AoA"  "Ty-Hy"  "Hy"  "Ty"'
        close(111)
    do  i=1,numSampBody
        write(fileName,'(I4.4)') SampBodyNode(i)
        open(111,file='./DatInfo/SampBodyNode'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "x"  "y"  "z"  "u"  "v"  "w"  "ax"  "ay"  "az" '
        close(111)
    enddo
    do  i=1,numSampFlow
        write(fileName,'(I4.4)') i
        open(111,file='./DatInfo/SampFlowPint'//trim(fileName)//'.plt')
        write(111,*)'variables= "t"  "p" "u"  "v"  "w" '
        close(111)
    enddo

        open(111,file='./DatInfo/Power.plt')
        write(111,*)'variables= "t" "Ptot" "Paero" "Piner" "Pax" "Pay" "Paz" "Pix" "Piy" "Piz"'
        close(111)

        open(111,file='./DatInfo/Converg.plt')
        write(111,*)'variables= "t"  "Convergence"  '
        close(111)

        open(111,file='./DatInfo/MaMax.plt')
        write(111,*)'variables= "t"  "MaMax"  '
        close(111)

        open(111,file='./DatInfo/MaxValBody.plt')
        write(111,*)'variables= "t"  "xyzMax" "velMax" "accMax"  '
        close(111)

        open(111,file='./DatInfo/Area.plt')
        write(111,*)'variables= "t"  "Area"  '
        close(111)

        open(111,file='./DatInfo/Energy.plt')
        write(111,*)'variables= "t","Es","Eb","Ep","Ek","Ew","Et"'
        close(111)

    END SUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!     write data
!    copyright@ RuNanHua 
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��     
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE wrtInfo()
    USE simParam
    implicit none
    real(8):: EEE(2),strainEnergy(nEL,2)
    integer:: i,j,k,nt,iEL,z,y,x,zbgn,ybgn,xbgn,zend,yend,xend
    real(8):: convergence,MaMax,weightm,velocity(1:3),Pressure
    real(8):: Ptot,Paero,Piner,Pax,Pay,Paz,Pix,Piy,Piz
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3,ax,ay,az
    integer,parameter::nameLen=4
    character (LEN=nameLen):: fileName

        open(111,file='./DatInfo/ForceDirect.plt',position='append')
        write(111,'(4D20.10)')time/Tref,sum(extful1(1:nND,1:3),1)/Fref
                                        
        close(111)

        open(111,file='./DatInfo/ForceStress.plt',position='append')
        write(111,'(4D20.10)')time/Tref,sum(extful2(1:nND,1:3),1)/Fref
                                        
        close(111)

    !==============================================================================================       
        !===
    do  i=1,5
        write(fileName,'(I4.4)') NDtl(i)
        open(111,file='./DatInfo/SampBodyNodeB'//trim(fileName)//'.plt',position='append')
        write(111,'(10D20.10)')time/Tref,xyzful(NDtl(i),1:3)/Lref,velful(NDtl(i),1:3)/Uref,accful(NDtl(i),1:3)/Aref                  
        close(111)
    enddo
        !===
    do  i=1,3
        write(fileName,'(I4.4)') NDhd(i)
        open(111,file='./DatInfo/SampBodyNodeB'//trim(fileName)//'.plt',position='append')
        write(111,'(10D20.10)')time/Tref,xyzful(NDhd(i),1:3)/Lref,velful(NDhd(i),1:3)/Uref,accful(NDhd(i),1:3)/Aref                     
        close(111)
    enddo
        !===
        open(111,file='./DatInfo/SampBodyCentP.plt',position='append')
        write(111,'(10D20.10)')time/Tref,xyzful(NDct,1:3)/Lref,velful(NDct,1:3)/Uref,accful(NDct,1:3)/Aref                     
        close(111)

        open(111,file='./DatInfo/SampBodyCentM.plt',position='append')
        write(111,'(10D20.10)')time/Tref,sum(xyzful(1:nND,1:3)*mssful(1:nND,1:3),1)/sum(mssful(1:nND,1:3),1)/Lref, &                                                                          
                                         sum(velful(1:nND,1:3)*mssful(1:nND,1:3),1)/sum(mssful(1:nND,1:3),1)/Uref, &
                                         sum(accful(1:nND,1:3)*mssful(1:nND,1:3),1)/sum(mssful(1:nND,1:3),1)/Aref                                                             
        close(111)

        open(111,file='./DatInfo/SampBodyAngular1.plt',position='append')
        write(111,'(5D20.10)')time/Tref,datan((xyzful(NDtl(1),3)-xyzful(NDhd(1),3))/(xyzful(NDtl(1),1)-xyzful(NDhd(1),1))),    &
                                        xyzful(NDtl(1),3)/Lref-xyzful(NDhd(1),3)/Lref,xyzful(NDhd(1),3)/Lref,xyzful(NDtl(1),3)/Lref
        close(111)
        open(111,file='./DatInfo/SampBodyAngular2.plt',position='append')
        write(111,'(5D20.10)')time/Tref,datan((xyzful(NDtl(2),3)-xyzful(NDhd(2),3))/(xyzful(NDtl(2),1)-xyzful(NDhd(2),1))),    &
                                        xyzful(NDtl(2),3)/Lref-xyzful(NDhd(2),3)/Lref,xyzful(NDhd(2),3)/Lref,xyzful(NDtl(2),3)/Lref
        close(111)
        open(111,file='./DatInfo/SampBodyAngular3.plt',position='append')
        write(111,'(5D20.10)')time/Tref,datan((xyzful(NDtl(3),3)-xyzful(NDhd(3),3))/(xyzful(NDtl(3),1)-xyzful(NDhd(3),1))),    &
                                        xyzful(NDtl(3),3)/Lref-xyzful(NDhd(3),3)/Lref,xyzful(NDhd(3),3)/Lref,xyzful(NDtl(3),3)/Lref
        close(111)
   
    do  i=1,numSampBody
        write(fileName,'(I4.4)') SampBodyNode(i)
        open(111,file='./DatInfo/SampBodyNode'//trim(fileName)//'.plt',position='append')
        write(111,'(10D20.10)')time/Tref, xyzful(SampBodyNode(i),1:3)/Lref,velful(SampBodyNode(i),1:3)/Uref,accful(SampBodyNode(i),1:3)/Aref 
        close(111)
    enddo
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
        write(111,'(5D20.10)') time/Tref,  Pressure/(0.5*denIn*Uref**2), velocity(1:3)/Uref
        close(111)

    enddo
    !==============================================================================================  

        Pax=sum(extful(1:nND,1)*velful(1:nND,1))/Pref
        Pay=sum(extful(1:nND,2)*velful(1:nND,2))/Pref
        Paz=sum(extful(1:nND,3)*velful(1:nND,3))/Pref
        Paero=Pax+Pay+Paz
        !write(*,*)'Pt:',Paero
        !write(*,*)'Pr:',(sum(extful(1:nND,4)*velful(1:nND,4))+sum(extful(1:nND,5)*velful(1:nND,5))+sum(extful(1:nND,6)*velful(1:nND,6)))/Pref
        Pix=-sum(mssful(1:nND,1)*accful(1:nND,1)*velful(1:nND,1))/Pref
        Piy=-sum(mssful(1:nND,2)*accful(1:nND,2)*velful(1:nND,2))/Pref
        Piz=-sum(mssful(1:nND,3)*accful(1:nND,3)*velful(1:nND,3))/Pref
        Piner=Pix+Piy+Piz
        Ptot=Paero+Piner
        open(111,file='./DatInfo/Power.plt',position='append')
        write(111,'(10D20.10)')time/Tref,Ptot,Paero,Piner,Pax,Pay,Paz,Pix,Piy,Piz                                                                      
        close(111)


        UNow=sum(dsqrt( uuu(:,:,:,1)**2+uuu(:,:,:,2)**2+uuu(:,:,:,3)**2))
        convergence=dabs(UNow-UPre)/UNow 
        UPre=UNow
        open(111,file='./DatInfo/Converg.plt',position='append')
        write(111,'(2D20.10)')time/Tref,convergence             
        close(111)

        MaMax=MaxVal(dsqrt( uuu(:,:,:,1)**2+uuu(:,:,:,2)**2+uuu(:,:,:,3)**2))/dsqrt(Cs2)
        open(111,file='./DatInfo/MaMax.plt',position='append')
        write(111,'(2D20.10)')time/Tref,MaMax             
        close(111)

        open(111,file='./DatInfo/MaxValBody.plt',position='append')
        write(111,'(4D20.10)')time/Tref,maxval(dsqrt(xyzful(1:nND,1)**2+xyzful(1:nND,2)**2+xyzful(1:nND,3)**2))/Lref, &
                                        maxval(dsqrt(velful(1:nND,1)**2+velful(1:nND,2)**2+velful(1:nND,3)**2))/Uref, &
                                        maxval(dsqrt(accful(1:nND,1)**2+accful(1:nND,2)**2+accful(1:nND,3)**2))/Aref
        close(111)

        !�������
        call cptArea(areaElem,nND,nEL,ele,xyzful)
        open(111,file='./DatInfo/Area.plt',position='append')
        write(111,'(2D20.10)')time/Tref,sum(areaElem(:))/Asfac
        close(111)
        
        call strain_energy_D(  strainEnergy,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3),&
                                 ele,prop,triad_n1,triad_n2,triad_n3,triad_ee,triad_e0,triad_nn, &
                                 nND,nEL,nMT &
                               )
        EEE(1)=sum(strainEnergy(1:nEL,1))
        EEE(2)=sum(strainEnergy(1:nEL,2))
        Es=EEE(1)/Eref
        Eb=EEE(2)/Eref
        Ep=Es+Eb
        Ew=Ew+Paero*timeOutInfo
        Ek=0.5*sum(mssful(1:nND,1:6)*velful(1:nND,1:6)*velful(1:nND,1:6))/Eref
        
        !write(*,*)'Ekt:', 0.5*sum(mssful(1:nND,1:3)*velful(1:nND,1:3)*velful(1:nND,1:3))/Eref
        !write(*,*)'Ekr:', 0.5*sum(mssful(1:nND,4:6)*velful(1:nND,4:6)*velful(1:nND,4:6))/Eref

        Et=Ek+Ep
        open(111,file='./DatInfo/Energy.plt', position='append')
        write(111,'(7D20.10)')time/Tref,Es,Eb,Ep,Ek,Ew,Et
        close(111)

        streI(1:nND)=0.0d0
        bendO(1:nND)=0.0d0
        do  iEL=1,nEL
            streI(ele(iEL,1:3))=streI(ele(iEL,1:3))+strainEnergy(iEL,1)/3
            bendO(ele(iEL,1:3))=bendO(ele(iEL,1:3))+strainEnergy(iEL,2)/3           
        enddo



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
!    copyright@ RuNanHua 
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��            
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  cptMove(move,xB,xG,d)
    implicit none
    integer:: move(3),i
    real(8):: xB(3),xG(3),d(3)
    do  i=1,3
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
    integer:: i,j,k,nt,iND,iEL
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
        ax =(x2-x1)*(x2-x1)
        ay =(y2-y1)*(y2-y1)
        az =(z2-z1)*(z2-z1)
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az)
        endif

        if(nt==3)then        
        x1=xyzful(i,1)
        x2=xyzful(j,1)
        x3=xyzful(k,1)
        y1=xyzful(i,2)
        y2=xyzful(j,2)
        y3=xyzful(k,2)
        z1=xyzful(i,3)
        z2=xyzful(j,3)
        z3=xyzful(k,3)

        ax =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))/2.0
        ay =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))/2.0
        az =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))/2.0
        areaElem(iEL)=dsqrt( ax*ax + ay*ay + az*az)
        endif

    enddo
    endsubroutine
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!
!    copyright@ RuNanHua 
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��            
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
!    ��Ȩ���У������ϣ��й��ƴ������ѧϵ��
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
    real(8):: rx,ry,rz,Phi
    
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
    end 

