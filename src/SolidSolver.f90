module SolidSolver
    implicit none
    private
    INTEGER, PARAMETER:: SpcDim = 3
    integer, parameter:: idat=12
    integer:: nFish!,nND_max,nEL_max,nMT_max,nEQ_max,NDref
    real(8):: dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,alpham,alphap
    real(8):: dtolFEM
    integer:: ntolFEM
    real(8):: g(3)
    real(8):: deltaT
    public :: BeamSolver
    type :: BeamSolver
        real(8), allocatable :: r_Lspan(:)
        real(8), allocatable :: r_Rspan(:)
        integer, allocatable :: r_Nspan(:)
        real(8) :: r_dirc(3)
        real(8):: Freq
        real(8), allocatable:: elmax(:),elmin(:)
        real(8), allocatable:: XYZ(:),XYZo(:),XYZAmpl(:),XYZPhi(:),XYZd(:),UVW(:)
        real(8), allocatable:: AoA(:),AoAo(:),AoAAmpl(:),AoAPhi(:),AoAd(:),WWW1(:),WWW2(:),WWW3(:)
        real(8), allocatable:: TTT00(:,:),TTT0(:,:),TTTnow(:,:),TTTnxt(:,:)
        integer:: nND,nEL,nEQ,nMT,nBD,nSTF
        integer, allocatable:: isMotionGiven(:)
        integer, allocatable:: ele(:,:),jBC(:,:),nloc(:),nprof(:),nprof2(:)
        real(8), allocatable:: xyzful00(:,:),mssful(:,:),vBC(:,:),prop(:,:),mss(:),areaElem00(:),areaElem(:)
        real(8), allocatable:: lodful(:,:),repful(:,:),extful(:,:),extful1(:,:),extful2(:,:),grav(:,:)

        real(8), allocatable:: xyzful0(:,:),xyzfulnxt(:,:),dspful(:,:),accful(:,:)
        real(8), allocatable:: xyzful(:,:),velful(:,:)
        real(8), allocatable:: triad_nn(:,:,:),triad_ee(:,:,:),triad_e0(:,:,:)
        real(8), allocatable:: triad_n1(:,:,:),triad_n2(:,:,:),triad_n3(:,:,:)
        real(8), allocatable:: FishInfo(:)
    contains
        procedure :: Allocate_solid => Allocate_solid_
        procedure :: Initialise => Initialise_
        procedure :: structure => structure_
    end type BeamSolver
  contains
    SUBROUTINE Allocate_solid_(this,FEmeshName,nAsfac,nLchod,lentemp)
        implicit none
        class(BeamSolver), intent(inout) :: this
        real(8):: nAsfac,nLchod,lentemp
        integer:: iND
        character (LEN=40):: FEmeshName
        
        write(*,'(A)') '=============================================================================='
        open(unit=idat, file = trim(adjustl(FEmeshName)))
        rewind(idat)
        read(idat,*)
        read(idat,*)this%nND,this%nEL,this%nMT
        read(idat,*)

        this%nEQ=this%nND*6
    
    !   ===============================================================================================
        allocate( this%ele(this%nEL,5),this%xyzful00(this%nND,6),this%xyzful0(this%nND,6),this%mssful(this%nND,6),this%lodful(this%nND,6), &
                  this%extful(this%nND,6),this%repful(this%nND,1:6),this%extful1(this%nND,6),this%extful2(this%nND,6),this%nloc(this%nND*6),this%nprof(this%nND*6), &
                  this%nprof2(this%nND*6),this%jBC(this%nND,6))
        allocate( this%grav(this%nND,6),this%vBC(this%nND,6),this%mss(this%nND*6),this%prop(this%nMT,10),this%areaElem00(this%nEL),this%areaElem(this%nEL))
    
        allocate( this%xyzful(this%nND,6),this%xyzfulnxt(this%nND,6),this%dspful(this%nND,6),this%velful(this%nND,6),this%accful(this%nND,6))
        allocate( this%triad_nn(3,3,this%nND),this%triad_ee(3,3,this%nEL),this%triad_e0(3,3,this%nEL) )
        allocate( this%triad_n1(3,3,this%nEL),this%triad_n2(3,3,this%nEL),this%triad_n3(3,3,this%nEL) )
    
        this%repful(:,:) =0.d0
        this%extful1(:,:)=0.d0
        this%extful2(:,:)=0.d0
    
    !   ===============================================================================================
        call read_structural_datafile(this%jBC(1:this%nND,1:6),this%ele(1:this%nEL,1:5),this%nloc(1:this%nND*6),this%nprof(1:this%nND*6), &
                                      this%nprof2(1:this%nND*6),this%xyzful00(1:this%nND,1:6),this%prop(1:this%nMT,1:10),this%nND, &
                                      this%nEL,this%nEQ,this%nMT,this%nBD,this%nSTF,idat)
        close(idat)

    !   ===============================================================================================
    !   calculate area
        call cptArea(this%areaElem00(1:this%nEL),this%nND,this%nEL,this%ele(1:this%nEL,1:5),this%xyzful00(1:this%nND,1:6))
        nAsfac=sum(this%areaElem00(1:this%nEL))
        this%elmax=maxval(this%areaElem00(1:this%nEL))
        this%elmin=minval(this%areaElem00(1:this%nEL))
        !calculate spanwise length, chord length, aspect ratio
        nLchod  = maxval(this%xyzful00(:,1))-minval(this%xyzful00(:,1))
        lentemp = maxval(this%xyzful00(:,2))-minval(this%xyzful00(:,2))
        if(lentemp .gt. nLchod) nLchod = lentemp
    
    !   loading boundary type*******************************************************************************
        do    iND=1,this%nND
            if(this%jBC(iND,1)==1) this%jBC(iND,1:6)=this%isMotionGiven(1:6)
        enddo
    end subroutine Allocate_solid_
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Allocate memory for solid simulation
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE allocate_solid_memory()
        USE simParam
        implicit none
        integer:: iFish,maxN
        real(8), allocatable:: nAsfac(:),nLchod(:),lentemp(:)
        if(nFish==0) return

        !allocate(BeamInfo(nFish))
        allocate(nAsfac(1:nFish),nLchod(1:nFish),lentemp(1:nFish))
        
        do iFish = 1,nFish
            !call BeamInfo(iFish)%Allocate_solid(FEmeshName(iFish),nAsfac(iFish),nLchod(iFish),lentemp(iFish))
            write(*,*)'read FEMeshFile ',iFish,' end'
        enddo

        ! nND_max=maxval(BeamInfo(:)%nND)
        ! nEL_max=maxval(BeamInfo(:)%nEL)
        ! nMT_max=maxval(BeamInfo(:)%nMT)
        ! nEQ_max=nND_max*6
        !Use the object with the largest area as the reference object
        maxN  = maxloc(nAsfac, dim=1)
        Asfac = nAsfac(maxN)
        Lchod = nLchod(maxN)
    
        if((Lchod-1.0d0)<=1.0d-2)Lchod=1.0d0
        if((maxval(Lspan)-1.0d0)<=1.0d-2)Lspan(maxloc(Lspan))=1.0d0
    
        AR    = maxval(Lspan)**2/Asfac
    END SUBROUTINE allocate_solid_memory

    SUBROUTINE Initialise_(this,pi,time)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer:: iND
        real(8):: pi,time
        this%TTT00(:,:)=0.0d0
        this%TTT00(1,1)=1.0d0
        this%TTT00(2,2)=1.0d0
        this%TTT00(3,3)=1.0d0

        this%XYZ(1:3)=this%XYZo(1:3)+this%XYZAmpl(1:3)*dcos(2.0*pi*this%Freq*time+this%XYZPhi(1:3))
        this%AoA(1:3)=this%AoAo(1:3)+this%AoAAmpl(1:3)*dcos(2.0*pi*this%Freq*time+this%AoAPhi(1:3))

        call AoAtoTTT(this%AoA(1:3),this%TTT0(1:3,1:3))
        call AoAtoTTT(this%AoA(1:3),this%TTTnow(1:3,1:3))
        call get_angle_triad(this%TTT0(1:3,1:3),this%TTTnow(1:3,1:3),this%AoAd(1),this%AoAd(2),this%AoAd(3))

        do iND=1,this%nND
            this%xyzful0(iND,1:3)=matmul(this%TTT0(1:3,1:3),this%xyzful00(iND,1:3))+this%XYZ(1:3)
            this%xyzful0(iND,4:6)=this%AoAd(1:3)
        enddo

        this%xyzful(1:this%nND,1:6)=this%xyzful0(1:this%nND,1:6)
        this%velful(1:this%nND,1:6)=0.0

        this%dspful(1:this%nND,1:6)=0.0
        this%accful(1:this%nND,1:6)=0.0

        if(this%ele(1,4).eq.2)then
            CALL formmass_D(this%ele(1:this%nEL,1:5),this%xyzful0(1:this%nND,1),this%xyzful0(1:this%nND,2),this%xyzful0(1:this%nND,3), &
                            this%prop(1:this%nMT,1:10),this%mss(1:this%nND*6),this%nND,this%nEL,this%nEQ,this%nMT,alphaf)

            do iND = 1, this%nND
                this%mssful(iND,1:6)= this%mss((iND-1)*6+1:(iND-1)*6+6)
                this%grav(iND,1:6)  = this%mssful(iND,1:6)*[g(1),g(2),g(3),0.0d0,0.0d0,0.0d0]
            enddo

            CALL init_triad_D(this%ele(1:this%nEL,1:5),this%xyzful(1:this%nND,1),this%xyzful(1:this%nND,2),this%xyzful(1:this%nND,3), &
                            this%triad_nn(1:3,1:3,1:this%nND),this%triad_n1(1:3,1:3,1:this%nEL),this%triad_n2(1:3,1:3,1:this%nEL), &
                            this%triad_ee(1:3,1:3,1:this%nEL),this%triad_e0(1:3,1:3,1:this%nEL),this%nND,this%nEL)
        else
            do iND = 1, this%nND
                this%mssful(iND,1:6)= 1.0d0 !Culculate accCentM/velCentM/xyzCentM requires mass not zero
                this%grav(iND,1:6)  = this%mssful(iND,1:6)*[g(1),g(2),g(3),0.0d0,0.0d0,0.0d0]
            enddo
        endif

    END SUBROUTINE Initialise_
    !    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    initialize solid field
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE initialize_solid()
        USE simParam
        implicit none
        integer:: iFish
        if(nFish.eq.0) return
    
        do iFish=1,nFish
            call BeamInfo(iFish)%Initialise(pi,time)
        enddo
    END SUBROUTINE initialize_solid

    SUBROUTINE solver()
        USE simParam
        implicit none
        integer:: iFish,isubstep
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iFish)
        do iFish=1,nFish
            call BeamInfo(iFish)%structure(iFish,pi,time,isubstep,subdeltat,iBodyModel(iFish))
        enddo !do iFish=1,nFish
        !$OMP END PARALLEL DO
    END SUBROUTINE

    SUBROUTINE structure_(this,iFish,pi,time,isubstep,subdeltat,iBodyModel)
        implicit none
        class(BeamSolver), intent(inout) :: this
        integer:: iFish,iND,isubstep,iBodyModel
        real(8):: subdeltat,pi,time
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iND)
            if(iBodyModel==1)then     ! rigid body
                !======================================================
                !prescribed motion
                !------------------------------------------------------
                !translational displacement
                this%XYZ(1:3)=this%XYZo(1:3)+this%XYZAmpl(1:3)*dcos(2.0*pi*this%Freq*(time-deltat+isubstep*subdeltat)+this%XYZPhi(1:3))
                !rotational displacement
                this%AoA(1:3)=this%AoAo(1:3)+this%AoAAmpl(1:3)*dcos(2.0*pi*this%Freq*(time-deltat+isubstep*subdeltat)+this%AoAPhi(1:3))
                call AoAtoTTT(this%AoA(1:3),this%TTTnxt(1:3,1:3))
                call get_angle_triad(this%TTT0(1:3,1:3),this%TTTnxt(1:3,1:3),this%AoAd(1),this%AoAd(2),this%AoAd(3))
                !given displacement
                do  iND=1,this%nND
                    this%xyzfulnxt(iND,1:3)=matmul(this%TTTnxt(1:3,1:3),this%xyzful00(iND,1:3))+this%XYZ(1:3)
                    this%xyzfulnxt(iND,4:6)=this%AoAd(1:3)
                enddo
                this%xyzful(1:this%nND,1:6)=this%xyzfulnxt(1:this%nND,1:6)
                !------------------------------------------------------
                !translational velocity
                this%UVW(1:3) =-2.0*pi*this%Freq*this%XYZAmpl(1:3)*dsin(2.0*pi*this%Freq*(time-deltat+isubstep*subdeltat)+this%XYZPhi(1:3))
                !rotational velocity
                this%WWW1(1:3)=-2.0*pi*this%Freq*this%AoAAmpl(1:3)*dsin(2.0*pi*this%Freq*(time-deltat+isubstep*subdeltat)+this%AoAPhi(1:3))
                this%WWW2(1:3)=[this%WWW1(1)*dcos(this%AoA(2))+this%WWW1(3),    &
                                 this%WWW1(1)*dsin(this%AoA(2))*dsin(this%AoA(3))+this%WWW1(2)*dcos(this%AoA(3)),   &
                                 this%WWW1(1)*dsin(this%AoA(2))*dcos(this%AoA(3))-this%WWW1(2)*dsin(this%AoA(3))    ]
                this%WWW3(1:3)=matmul(this%TTTnxt(1:3,1:3),this%WWW2(1:3))
                !given velocity
                do  iND=1,this%nND
                    this%velful(iND,1:3)=[this%WWW3(2)*this%xyzful(iND,3)-this%WWW3(3)*this%xyzful(iND,2),    &
                                           this%WWW3(3)*this%xyzful(iND,1)-this%WWW3(1)*this%xyzful(iND,3),    &
                                           this%WWW3(1)*this%xyzful(iND,2)-this%WWW3(2)*this%xyzful(iND,1)    ]&
                                           + this%UVW(1:3)
                    this%velful(iND,4:6)=this%WWW3(1:3)
                enddo
                !-------------------------------------------------------
            elseif(iBodyModel==2)then !elastic model
                !translational displacement
                this%XYZ(1:3)=this%XYZo(1:3)+this%XYZAmpl(1:3)*dcos(2.0*pi*this%Freq*(time-deltat+isubstep*subdeltat)+this%XYZPhi(1:3))
                !rotational displacement
                this%AoA(1:3)=this%AoAo(1:3)+this%AoAAmpl(1:3)*dcos(2.0*pi*this%Freq*(time-deltat+isubstep*subdeltat)+this%AoAPhi(1:3))
                call AoAtoTTT(this%AoA(1:3),this%TTTnxt(1:3,1:3))
                call get_angle_triad(this%TTT0(1:3,1:3),this%TTTnxt(1:3,1:3),this%AoAd(1),this%AoAd(2),this%AoAd(3))
                !given displacement
                do  iND=1,this%nND
                    this%xyzfulnxt(iND,1:3)=matmul(this%TTTnxt(1:3,1:3),this%xyzful00(iND,1:3))+this%XYZ(1:3)
                    this%xyzfulnxt(iND,4:6)=this%AoAd(1:3)
                enddo
                !-------------------------------------------------------
                !displacement condition
                this%vBC(1:this%nND,1:6) = this%xyzfulnxt(1:this%nND,1:6) - this%xyzful(1:this%nND,1:6)
                !loading vector
                this%lodful(1:this%nND,1:6) = this%extful(1:this%nND,1:6) + this%grav(1:this%nND,1:6) + this%repful(1:this%nND,1:6)
                !-----------------------------------------
                CALL structure_solver(this%jBC(1:this%nND,1:6),this%vBC(1:this%nND,1:6),this%ele(1:this%nEL,1:5), &
                                      this%nloc(1:this%nND*6),this%nprof(1:this%nND*6),this%nprof2(1:this%nND*6), &
                                      this%prop(1:this%nMT,1:10),this%mss(1:this%nND*6), &
                                      this%xyzful0(1:this%nND,1:6),this%xyzful(1:this%nND,1:6),this%dspful(1:this%nND,1:6), &
                                      this%velful(1:this%nND,1:6),this%accful(1:this%nND,1:6),this%lodful(1:this%nND,1:6),  &
                                      subdeltat,dampK,dampM,  &
                                      this%triad_nn(1:3,1:3,1:this%nND),this%triad_ee(1:3,1:3,1:this%nEL), &
                                      this%triad_n1(1:3,1:3,1:this%nEL),this%triad_n2(1:3,1:3,1:this%nEL), &
                                      this%nND,this%nEL,this%nEQ,this%nMT,this%nBD,this%nSTF,NewmarkGamma,NewmarkBeta,dtolFEM,ntolFEM,    &
                                      iFish,this%FishInfo)
            else
                stop 'no define body model'
            endif
            !$OMP END PARALLEL DO
    END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    FEM code for structure
!    input variables (at time t):
!    xyzful0 (const), xyzful, coordinates
!    velful, velocity
!    accful, acceleration
!    lodExteful, externel force
!
!    output variables (at time t+dt):
!    velful
!    accful
!
!    cached variables:
!    dspful displacement at t+dt
!
!    workspace variables:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE structure_solver(jBC,vBC,ele,nloc,nprof,nprof2,prop,mss,xyzful0,xyzful,dspful,velful,accful,lodExteful,deltat,dampK,dampM,     &
                      triad_nn,triad_ee,triad_n1,triad_n2,nND,nEL,nEQ,nMT,nBD,maxstiff,Newmarkdelta,Newmarkalpha, &
                      dtol,iterMax,iFish,FishInfo)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nBD,maxstiff,iFish
    integer:: jBC(nND,6),ele(nEL,5),nloc(nEQ),nprof(nEQ),nprof2(nEQ)
    real(8):: vBC(nND,6),xyzful0(nND,6),xyzful(nND,6),prop(nMT,10)
    real(8):: mss(nEQ),dspful(nND,6),velful(nND,6),accful(nND,6),lodExteful(nND,6),deltat
    real(8):: dampK,dampM

    real(8):: triad_nn(3,3,nND),triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)
!   ----------------------------------------------------------------------------------------------
    real(8):: dspt(nEQ),du(nEQ),ddu(nEQ),lodEffe(nEQ),lodInte(nEQ),lodExte(nEQ),dsp(nEQ),vel(nEQ),acc(nEQ)
    real(8):: stfEffe(maxstiff),stfElas(maxstiff),stfGeom(maxstiff)
    real(8):: wk1(nEQ),wk2(nEQ),velO(nEQ),accO(nEQ)
    real(8):: Newmarkdelta,Newmarkalpha
    real(8):: a0,a1,a2,a3,a4,a5,a6,a7
    real(8):: beta0,beta,gamma,zi,z0
    real(8):: dtol,dnorm,geoFRM(nEL),FishInfo(1:3)
    integer:: i,j,iND,iEQ,iter,iterMax,iloc,ierror,maxramp,iModify,i1,i2
!   -----------------------------------------------------------------------------------------------
    a0 = 1.0d0/(Newmarkalpha*deltat*deltat)
    a1 = Newmarkdelta/(Newmarkalpha*deltat)
    a2 = 1.0d0/(Newmarkalpha*deltat)
    a3 = 1.0d0/(Newmarkalpha*2.0d0) - 1.0d0
    a4 = Newmarkdelta/Newmarkalpha - 1.0d0
    a5 = (Newmarkdelta/Newmarkalpha - 2.0d0)*0.5d0*deltat
    a6 = (1.0d0 - Newmarkdelta)*deltat
    a7 = Newmarkdelta*deltat

    beta0 =1.0d0
    gamma =1.0d0
    maxramp =0

    iModify = 1

    !igeoFRM=27
    !igeoPLT=28
    !open(unit=igeoFRM,file='geoFRM.tmp')
    !open(unit=igeoPLT,file='geoPLT.tmp')

!   ***********************************************************************************************
    do  i=1,nND
        i1 = (i-1)*6+1
        i2 = (i-1)*6+6
        dsp(i1:i2)=dspful(i,1:6)
        vel(i1:i2)=velful(i,1:6)
        acc(i1:i2)=accful(i,1:6)
        lodExte(i1:i2)=lodExteful(i,1:6)
    enddo

    dspt(1:nEQ)=dsp(1:nEQ) ! displacement at time t
!   ------------------------------------------------------------------------------
    iter=0
    dnorm=1.0d0
    do while(dnorm >= dtol .and. iter<= iterMax )
!   ------------------------------------------------------------------------------
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        forces from body stresses
        call body_stress_D( lodInte,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                            ele,prop,triad_n1,triad_n2,triad_ee,    &
                            nND,nEL,nEQ,nMT,geoFRM &
                          )

        do    i= 1, nEQ
            lodEffe(i)=lodExte(i)-lodInte(i)+(a0*(dspt(i)-dsp(i))+a2*vel(i)+a3*acc(i))*mss(i)   &
                                            +(a1*(dspt(i)-dsp(i))+a4*vel(i)+a5*acc(i))*dampM*mss(i)
        enddo


        if    (dampK > 0.0d0) then
            do    i= 1, nEQ
                wk1(i)= a1*(dspt(i)-dsp(i)) +a4*vel(i) +a5*acc(i)
            enddo

            call AxBCOL(stfElas,maxstiff,wk1,wk2,nEQ,nprof,nprof2,nloc)

            do    i= 1, nEQ
                lodEffe(i) = lodEffe(i) + dampK*wk2(i)
            enddo
        endif
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       -------------------------------------------------------------------
        call formstif_s(stfElas,ele,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                        prop,nloc,nND,nEL,nEQ,nMT,maxstiff)

        call formgeom_s(stfGeom,ele,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                        prop,nloc,nND,nEL,nEQ,nMT,maxstiff,geoFRM)
!
        do  i= 1, nEQ
            iloc=nloc(i)
            do    j= 1, nprof(i)
                stfEffe(iloc+j-1)=stfElas(iloc+j-1)+ gamma*stfGeom(iloc+j-1)
            enddo
            stfEffe(iloc) = stfEffe(iloc) + a0*mss(i) + a1*dampM*mss(i)
        enddo

        if    (dampK .gt. 0.0d0) then
            do    i= 1, nEQ
                iloc=nloc(i)
                do    j= 1, nprof(i)
                    stfEffe(iloc+j-1) = stfEffe(iloc+j-1) + a1*dampK*stfElas(iloc+j-1)
                enddo
            enddo
        endif
!       -------------------------------------------------------------------
!       apply boundary conditions
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do  iND=1,nND
        do  j=1  ,6
            if(jBC(iND,j)>0)then
                iEQ=(iND-1)*6+j
                stfEffe(nloc(iEQ))=stfEffe(nloc(iEQ))*1.0d20
                if(iter==0)then
                    lodEffe(iEQ)=stfEffe(nloc(iEQ))*vBC(iND,j)
                else
                    lodEffe(iEQ)=stfEffe(nloc(iEQ))*0.0d0
                endif
            endif
        enddo
        enddo

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       -------------------------------------------------------------------

        call uduCOL_D(stfEffe,maxstiff,nEQ,ierror, nprof,nloc)
        if    (ierror.eq.0) then
            write(*,*)'@@ ERROR: zero diagonal term !!!'
            return
        endif
        call bakCOL_D(stfEffe,maxstiff,lodEffe,nEQ,nBD,du,ierror,nprof,nprof2,nloc)
!       -------------------------------------------------------------------
        if    (iter <= maxramp) then
            zi=2.0d0**(iter)
            z0=2.0d0**(maxramp)
            beta=zi/z0*beta0
        else
            beta=1.0d0*beta0
        endif

        ddu(1:nEQ)= beta*du(1:nEQ)
        dsp(1:nEQ)= dsp(1:nEQ) + ddu(1:nEQ)


        do  i=1,nND
            dspful(i,1:6)=dsp((i-1)*6+1:(i-1)*6+6)
        enddo

        xyzful(1:nND,1:6)=xyzful0(1:nND,1:6)+dspful(1:nND,1:6)

        call update_triad_D(ele,ddu,triad_nn,triad_n1,triad_n2,nND,nEL,nEQ)

        call make_triad_ee(ele,xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3),triad_ee,triad_n1,triad_n2,nND,nEL)

!       -------------------------------------------------------------------
!       test for convergence
        dnorm=dabs(maxval((du(1:nEQ)*beta)**2))
        iter=iter+1
    enddo


    FishInfo(1)=dble(iFish)
    FishInfo(2)=dble(iter)
    FishInfo(3)=dnorm

    velO(1:nEQ) = vel(1:nEQ)
    accO(1:nEQ) = acc(1:nEQ)

    acc(1:nEQ)  = a0*(dsp(1:nEQ) - dspt(1:nEQ)) -a2*velO(1:nEQ) - a3*accO(1:nEQ)
    vel(1:nEQ)  = velO(1:nEQ) + a6*accO(1:nEQ) + a7*acc(1:nEQ)

    do  i=1,nND
        velful(i,1:6)=vel((i-1)*6+1:(i-1)*6+6)
        accful(i,1:6)=acc((i-1)*6+1:(i-1)*6+6)
    enddo

    return
    ENDSUBROUTINE

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ structural DaTafile
    subroutine read_structural_datafile(jBC,ele,nloc,nprof,nprof2,xyzful0,prop,nND,nEL,nEQ,nMT,nBD,maxstiff,idat)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nBD,maxstiff,idat
    integer:: ele(nEL,5),jBC(nND,6),nloc(nND*6),nprof(nND*6),nprof2(nND*6)
    real(8):: xyzful0(nND,6),prop(nMT,10)
!   ---------------------------------------------------------------------------
    integer:: i,j,nbc,node,nmp,ibandh,iend,ibandv,ji1
    character (LEN=50):: endin
!   -----------------------------------------------------------------------------------------------
!   READ  node
    read(idat,*)  nND
    do    i= 1, nND
        read(idat,*) node,xyzful0(node,1),xyzful0(node,2),xyzful0(node,3)
    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------
!   READ elem data
    ! element number, node left, node right, node right, element type, material property index
    read(idat,*) nEL
    do  i= 1, nEL
        read(idat,*) j,ele(j,1:5)
    enddo
    read(idat,'(1a50)') endin

!   -----------------------------------------------------------------------------------------------
!    READ  bcs  default is 0=free
    jBC(1:nND,1:6) = 0
    read(idat,*)  nbc
    do  i=1,nbc
        read(idat,*)node,jBC(node,1),jBC(node,2),jBC(node,3), &
                         jBC(node,4),jBC(node,5),jBC(node,6)

    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------
!   READ element material properties
    ! can have nMT types of material
    ! property data will be overwritten if isKB = 0 or 1
    read(idat,*) nMT
    do    i= 1, nMT
        read(idat,*) nmp,prop(nmp,1:8)
    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------

    nEQ=nND*6

    call maxbnd(ele,nprof,nEL,nEQ,nBD)

!   nprof2
    do i=1,nEQ
       ibandh=1
       iend=i+nBD-1
       if (iend .gt. nEQ) iend=nEQ
       do j=i+1,iend
            ibandv=nprof(j)
            ji1=j-i+1
            if  (ibandv .ge. ji1) then
                ibandh = ji1
            endif
       enddo
       nprof2(i)=ibandh
    enddo

    nloc(1)=1
    do    i=1,nEQ-1
        nloc(i+1) = nloc(i) + nprof(i)
    enddo
    maxstiff=nloc(nEQ)+nprof(nEQ)
    !nprof  :  1  2  3  4  5  6  7  8  9 10 11 12  7  8  9 10 11 12  7  8  9 ...
    !nprof2 : 12 11 10  9  8  7 12 11 10  9  8  7 12 11 10  9  8  7 12 11 10 ...

    !write(*,'(3(A,1x,I8,2x))')'nND=',nND,'nEL=',nEL,'nEQ=',nEQ
    !write(*,'(3(A,1x,I8,2x))')'nMT=',nMT,'nBD=',nBD,'maxstiff=',maxstiff
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   MAX BaNDwidth calculation
!   input: ele elements
!    nEL, number of elements
!    nEQ, number of dof, 6*nEL
!   output:
!    nprof
!    nBD maximum bandwidth
    subroutine maxbnd(ele,nprof,nEL,nEQ,nBD)
    implicit none
    integer:: nBD,nEQ,nEL
    integer:: ele(nEL,5),nprof(nEQ)

    integer:: ipv(18)
    integer:: i,j,n,ihfbnd,idof,jdof,kdof,ieqn1,ieqn2,jband,ieq,ieq2

    nprof(1:nEQ)=0

    ihfbnd=0
    do  n=1,nEL
        idof=(ele(n,1)-1)*6
        jdof=(ele(n,2)-1)*6
        kdof=(ele(n,3)-1)*6
        do  i=1,6
            ipv(i   )=idof+i
            ipv(i+6 )=jdof+i
            ipv(i+12)=kdof+i
        enddo

        do  i=1,18
            ieqn1=ipv(i)
            do  j=i,18
                ieqn2=ipv(j)
                ihfbnd = max0(ihfbnd,iabs(ieqn1-ieqn2))
                jband=abs(ieqn1-ieqn2)+1
                ieq=max(ieqn1,ieqn2)
                if  (jband .gt. nprof(ieq)) then
                    nprof(ieq)=jband
                endif
                ieq2=min(ieqn1,ieqn2)
!               if  (jband .gt. nprof2(ieq2)) then
!                   nprof2(ieq2)=jband
!               endif
            enddo
        enddo

    enddo
    nBD=ihfbnd+1

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   set angle of initial triads
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE init_triad_D(ele,xord0,yord0,zord0,triad_nn,triad_n1,triad_n2,triad_ee,triad_e0,nND,nEL)
    implicit none
    integer:: nND,nEL
    integer:: ele(nEL,5)
    real(8):: xord0(nND),yord0(nND),zord0(nND)
    real(8):: triad_nn(3,3,nND),triad_ee(3,3,nEL),triad_e0(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)
!
    real(8):: dx,dy,dz,xl0
    real(8):: xll1,xmm1,xnn1
    real(8):: dd

    integer:: i,j,n,i1,j1,k1,nELt

!   For each node, set the triad  to global system
    do  n=1,nND

        do i=1,3
        do j=1,3
            triad_nn(i,j,n)=0.0d0
        enddo
        enddo
        triad_nn(1,1,n)=1.0d0
        triad_nn(2,2,n)=1.0d0
        triad_nn(3,3,n)=1.0d0

    enddo
!
!   For each element, calculate the current orientation triad
    do  n=1,nEL

        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
!
        if    (nELt == 2) then
!           frame
            dx = xord0(j1) - xord0(i1)
            dy = yord0(j1) - yord0(i1)
            dz = zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll1=dx/xl0
            xmm1=dy/xl0
            xnn1=dz/xl0
            dd=dsqrt(xll1*xll1+xmm1*xmm1)
            triad_n1(1,1,n)=xll1
            triad_n1(2,1,n)=xmm1
            triad_n1(3,1,n)=xnn1

            if    (dd .lt. 0.001d0) then
                triad_n1(1,2,n)=0.0d0
                triad_n1(2,2,n)=1.0d0
                triad_n1(3,2,n)=0.0d0
                triad_n1(1,3,n)=-xnn1
                triad_n1(2,3,n)=0.00d0
                triad_n1(3,3,n)=0.00d0
            else
                triad_n1(1,2,n)=-xmm1/dd
                triad_n1(2,2,n)=+xll1/dd
                triad_n1(3,2,n)=0.0d0

                triad_n1(1,3,n)=-xll1*xnn1/dd
                triad_n1(2,3,n)=-xmm1*xnn1/dd
                triad_n1(3,3,n)= dd
            endif
!            all element triads have same initial orientation
            triad_n2(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_ee(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_e0(1:3,1:3,n)=triad_n1(1:3,1:3,n)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
    ENDSUBROUTINE init_triad_D


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   update angle of  triads
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE update_triad_D(ele,ddut,triad_nn,triad_n1,triad_n2,nND,nEL,nEQ)
    implicit none
    integer:: nND,nEL,nEQ
    integer:: ele(nEL,5)
    real(8):: ddut(nEQ)
    real(8):: triad_nn(3,3,nND),triad_n1(3,3,nEL),triad_n2(3,3,nEL)

    real(8):: ddutful(nND,6),rr(3,3)
    real(8):: dtx1,dty1,dtz1
    integer:: i,n,i1,j1,k1,node,jdof,ireduc,nELt

    do    i=1,nND*6
        node=(i+5)/6
        jdof=i-(node-1)*6
        ireduc=i
        ddutful(node,jdof)=ddut(ireduc)
    enddo
!
!   For each node, set the triad_nn = [R]triad_nn
    do    n=1,nND
        dtx1=ddutful(n,4)
        dty1=ddutful(n,5)
        dtz1=ddutful(n,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_nn(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_nn(1:3,1:3,n))
    enddo

!   For each element, calculate the current orientation triad
    do  n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
!       n1 node
        dtx1=ddutful(i1,4)
        dty1=ddutful(i1,5)
        dtz1=ddutful(i1,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_n1(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_n1(1:3,1:3,n))
!       n2 node
        dtx1=ddutful(j1,4)
        dty1=ddutful(j1,5)
        dtz1=ddutful(j1,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_n2(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_n2(1:3,1:3,n))
    enddo

    return
    ENDSUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   get orientation of element
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE make_triad_ee(ele,xord,yord,zord,triad_ee,triad_n1,triad_n2,nND,nEL)
    implicit none
    integer:: nND,nEL
    integer:: ele(nEL,5)
    real(8):: xord(nND), yord(nND), zord(nND)
!
    real(8):: triad_aa(3,3)
    real(8):: triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)
    real(8):: rr(3,3),tx,ty,tz
    real(8):: triad_11(3,3),triad_22(3,3)
!
    real(8):: dx,dy,dz,xl0

    real(8):: xll,xmm,xnn,dd,r2e1,r3e1
    integer:: i,j,k,n,i1,j1,k1,nELt
!
!

!   For each element, calculate the current orientation triad
    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
!
        if    (nELt == 2) then
!            frame
            dx = xord(j1) - xord(i1)
            dy = yord(j1) - yord(i1)
            dz = zord(j1) - zord(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl0
            xmm=dy/xl0
            xnn=dz/xl0
            dd =dsqrt(xll*xll+xmm*xmm)
            do    i=1,3
            do    j=1,3
                triad_ee(i,j,n)=0.0d0
            enddo
            enddo
            triad_ee(1,1,n)=xll
            triad_ee(2,1,n)=xmm
            triad_ee(3,1,n)=xnn
!
!            get angle between two triads
            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_n1(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo
            call get_angle_triad( triad_11,triad_22,tx,ty,tz)
!
!           rotate n1 to intermediate
            tx=tx/2.0d0
            ty=ty/2.0d0
            tz=tz/2.0d0
            call finite_rot(tx,ty,tz,rr)
            triad_aa(1:3,1:3)=matmul(rr(1:3,1:3),triad_n1(1:3,1:3,n))
!
!
!           vectors e2 e3
            r2e1 = 0.0d0
            r3e1 = 0.0d0
            do    k=1,3
                r2e1 = r2e1 + triad_aa(k,2)*triad_ee(k,1,n)
                r3e1 = r3e1 + triad_aa(k,3)*triad_ee(k,1,n)
            enddo
            do    j=1,3
                triad_ee(j,2,n)=triad_aa(j,2) - r2e1*(triad_aa(j,1)+triad_ee(j,1,n))/2.0d0
                triad_ee(j,3,n)=triad_aa(j,3) - r3e1*(triad_aa(j,1)+triad_ee(j,1,n))/2.0d0
            enddo
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo
!   end of loop over elements
    return
    ENDSUBROUTINE make_triad_ee

!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE finite_rot(t1,t2,t3,rr)
    implicit none
    real(8):: rr(3,3),rr1(3,3),rr2(3,3),rr3(3,3)
    real(8):: t1,t2,t3,tt,ss,cc,c1,c2

    integer:: i,j
!
    tt=dsqrt( t1**2 + t2**2 + t3**2 )
    ss=dsin(tt)
    cc=dcos(tt)
!
    rr1(1,1)=1.0d0
    rr1(1,2)=0.0d0
    rr1(1,3)=0.0d0
    rr1(2,1)=0.0d0
    rr1(2,2)=1.0d0
    rr1(2,3)=0.0d0
    rr1(3,1)=0.0d0
    rr1(3,2)=0.0d0
    rr1(3,3)=1.0d0
!
    rr2(1,1)=0.0d0
    rr2(1,2)=-t3
    rr2(1,3)= t2
    rr2(2,1)= t3
    rr2(2,2)=0.0d0
    rr2(2,3)=-t1
    rr2(3,1)=-t2
    rr2(3,2)= t1
    rr2(3,3)=0d0
!
    rr3(1,1)=-t3*t3-t2*t2
    rr3(1,2)= t2*t1
    rr3(1,3)= t3*t1
    rr3(2,1)= t1*t2
    rr3(2,2)=-t3*t3-t1*t1
    rr3(2,3)= t3*t2
    rr3(3,1)= t1*t3
    rr3(3,2)= t2*t3
    rr3(3,3)=-t2*t2-t1*t1
!
    do    i=1,3
    do    j=1,3
        if    (tt .lt. 1.0d-10) then
            c1=1.0d0
!!!!        c2=1.0
            c2=0.5d0
        else
            c1 = ss/tt
            c2 = (1.0d0-cc)/tt**2
        endif
        rr(i,j) = rr1(i,j) + rr2(i,j)*c1 + rr3(i,j)*c2
    enddo
    enddo
!
    return
    ENDSUBROUTINE finite_rot


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   GET ANGLE of between TRIADs
!    pass
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE get_angle_triad(triad_n1,triad_n2,tx,ty,tz)
    implicit none
    real(8):: triad_n1(3,3),triad_n2(3,3)
    real(8):: rr(3,3)
    real(8):: tx,ty,tz, dtx,dty,dtz,c1,tt,sint
    integer:: i,j,k
!
!   get angle between two triads
    do    i=1,3
    do    j=1,3
        rr(i,j)=0.0d0
        do    k=1,3
            rr(i,j)=rr(i,j) + triad_n2(i,k)*triad_n1(j,k)
        enddo
    enddo
    enddo

    dtx = (rr(3,2)-rr(2,3))/2.0d0
    dty = (rr(1,3)-rr(3,1))/2.0d0
    dtz = (rr(2,1)-rr(1,2))/2.0d0

    c1=1.0d0
    sint = dsqrt(dtx*dtx+dty*dty+dtz*dtz)

    if (sint .gt. 1.0d0) sint=1.0d0
    tt = dasin(sint)
    if ( sint .lt. 1.0d-6) then
         c1=1.0d0
    else
         c1 = tt/sint
    endif

    tx=c1*dtx
    ty=c1*dty
    tz=c1*dtz

    return
    ENDSUBROUTINE get_angle_triad
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE global_to_local(triad,tx,ty,tz,tx2,ty2,tz2)
    implicit none
    real(8):: triad(3,3)
    real(8):: tx,ty,tz,tx2,ty2,tz2
    tx2 = triad(1,1)*tx+triad(2,1)*ty+triad(3,1)*tz
    ty2 = triad(1,2)*tx+triad(2,2)*ty+triad(3,2)*tz
    tz2 = triad(1,3)*tx+triad(2,3)*ty+triad(3,3)*tz
    return
    ENDSUBROUTINE global_to_local

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Nodal loads due to body stresses
!   input: xord0,yord0,zord0,xord,yord,zord undeformed and deformed coordinates
!   ele elements [nodex, nodey]
!   prop material property
!   nMT number of material types
!   triad_n1,triad_n2,triad_n3,triad_ee,triad_e0
!   output: gforce global force
!    geoFRM
!    geoPLT
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE body_stress_D(    gforce,xord0,yord0,zord0,xord,yord,zord, &
                                ele,prop,triad_n1,triad_n2,triad_ee, &
                                nND,nEL,nEQ,nMT,geoFRM &
                            )
    implicit none
    integer:: nMT,nEL,nND,nEQ
    integer:: ele(nEL,5)
    real(8):: xord0(nND), yord0(nND), zord0(nND)
    real(8):: xord(nND), yord(nND), zord(nND)
!
    real(8):: prop(nMT,10),geoFRM(nEL)

    real(8):: force(18),forceb(18)
    real(8):: gforce(nEQ)
    real(8):: ekb12(12,12)
!
    real(8):: triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)

    real(8):: triad_00(3,3),triad_11(3,3),triad_22(3,3)
    real(8):: rr(3,3)
    real(8):: ub(18),dl
!
    real(8):: dx0,dy0,dz0,du,dv,dw
    real(8):: dx,dy,dz,xl0
    real(8):: tx1,tx2,ty1,ty2,tz1,tz2,tx,ty,tz

    real(8):: e0,g0,a0,b0,r0,zix0,ziy0,ziz0,xl,fxx
    integer:: i,j,k,n,i1,j1,k1,mat,nELt

!    pi=4.0*datan(1.0d0)

    gforce(1:nEQ)=0.0d0

    !rewind(igeoFRM)
    !rewind(igeoPLT)


!   For each element, calculate the nodal forces
    do    n=1,nEL

        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        mat = ele(n,5)
        if    ( nELt == 2) then
!           frame
            e0  =prop(mat,1)
            g0  =prop(mat,2)
            a0  =prop(mat,3)
            r0  =prop(mat,4)
            b0  =prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)

            dx0 = xord0(j1) - xord0(i1)
            dy0 = yord0(j1) - yord0(i1)
            dz0 = zord0(j1) - zord0(i1)
            xl0 = dsqrt(dx0*dx0+dy0*dy0+dz0*dz0) ! initial vector of segment (dx0, dy0, dz0)
!
!           orientation
            du = (xord(j1)-xord0(j1))-(xord(i1)-xord0(i1))
            dv = (yord(j1)-yord0(j1))-(yord(i1)-yord0(i1))
            dw = (zord(j1)-zord0(j1))-(zord(i1)-zord0(i1))

            dx = dx0 + du
            dy = dy0 + dv
            dz = dz0 + dw
            xl =dsqrt(dx*dx+dy*dy+dz*dz) ! deformed vector of segment (dx, dy, dz)
!
            dl = ( (2*dx0+du)*du +(2*dy0+dv)*dv +(2*dz0+dw)*dw )/ (xl+xl0) ! [(dx0 + dx) . du]/(l + l0)
!           get twisting angles
            do    i=1,3
            do    j=1,3
                triad_00(i,j)=triad_ee(i,j,n)
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n1(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)
!
            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

!            non-zero ty1 tz1 u2 tx2 ty2 tz2
            ub(1)=0.0d0
            ub(2)=0.0d0
            ub(3)=0.0d0
            ub(4)=tx1
            ub(5)=ty1
            ub(6)=tz1
!
            ub(7)=dl
            ub(8)=0.0d0
            ub(9)=0.0d0
            ub(10)=tx2
            ub(11)=ty2
            ub(12)=tz2
!
!
!           get current stiffness in local coords. use L0

            call elmstfFRM_D(xl0,zix0,ziy0,ziz0,a0,e0,g0,ekb12,nELt )
!
!           compute axial force
            fxx=dl*e0*a0/xl0
!           save local force for geo stiff
            !write(igeoFRM,'(E25.15)') fxx
            geoFRM(n)=fxx
!           nodal forces in local coords
!           {F}=[k]{u}
            forceb(1:12) =matmul(ekb12(1:12,1:12),ub(1:12))
            forceb(13:18)=0.0d0

!           transform to global
            do  i=1,3
            do  j=1,3
                rr(i,j)=triad_ee(i,j,n)
            enddo
            enddo
            do  i=1,18
                force(i)=0.0d0
            enddo
            do    i=1,3
            do    k=1,3
                force(0+i) = force(0+i) + rr(i,k)*forceb(0+k)
                force(3+i) = force(3+i) + rr(i,k)*forceb(3+k)
                force(6+i) = force(6+i) + rr(i,k)*forceb(6+k)
                force(9+i) = force(9+i) + rr(i,k)*forceb(9+k)
            enddo
            enddo

            call assembFOR(nEQ,gforce,force,i1,j1,k1)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif

    enddo

    return
    ENDSUBROUTINE body_stress_D

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   FORM MASS matrix: only lumped
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine formmass_D(ele,xord,yord,zord,prop,mss,nND,nEL,nEQ,nMT,alphaf)
    implicit none
    integer:: nND,nEL,nEQ,nMT
    integer:: ele(nEL,5)
    real(8):: prop(nMT,10)
!
    real(8):: xord(nND), yord(nND), zord(nND)
    real(8):: mss(nEQ), em(18,18)
    real(8):: em12(12,12)

    integer:: i,i1,j1,k1,mat,nELt
    real(8):: a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl,xll,xmm,xnn

    !ת���������Ӵ�С
    real(8):: alphaf

!   zero array before assembling

    mss(1:nEQ)=0.0d0
!
!   form the element form matrix, and assemble
    do  i=1,nEL
        i1  = ele(i,1)
        j1  = ele(i,2)
        k1  = ele(i,3)
        nELt= ele(i,4) ! element type, 2 line segment
        mat = ele(i,5)
        if (nELt == 2) then
            a0=prop(mat,3)
            r0=prop(mat,4)
            b0=prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl
            call elmmasFRM_D(r0,a0,xl,zix0,em12,alphaf)
            call trans3d_D(xll,xmm,xnn,em12,b0)
            em(1:18,1:18)=0.0d0
            em(1:12,1:12)=em12(1:12,1:12)
            call assembLUM(nEQ,mss,em,i1,j1,k1)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   FORM STIFfness matrix  [K]
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine formstif_s(stf,ele,xord0,yord0,zord0,xord,yord,zord,prop,nloc,nND,nEL,nEQ,nMT,maxstiff)
    implicit none
    integer:: nND,nEL,nEQ,nMT,maxstiff
    integer:: ele(nEL,5), nloc(nEQ)
    real(8):: prop(nMT,10)
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)
    real(8):: stf(maxstiff)
    real(8):: ek(18,18), ek12(12,12)

    integer:: j,k,n,i1,j1,k1,mat,nELt
    real(8):: e0,g0,a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl0,xl,xll,xmm,xnn,xl9

!   initialize [K]  to zero
    stf(1:maxstiff)=0.0d0

!   form each element matrix, and assemble
    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        mat = ele(n,5)
!
        if    (nELt == 2) then
!           frame
!           material props not change
            e0=prop(mat,1)    !1 E modulus
            g0=prop(mat,2)    !2 G shear modulus
            a0=prop(mat,3)    !3 A cross section area
            r0=prop(mat,4)    !4 rho density
            b0=prop(mat,5)    !5 self-rotation angle in degree
            zix0=prop(mat,6)  !6 J torsion rigidity of the cross section
            ziy0=prop(mat,7)  !7 Iy bending rigidity in z direction
            ziz0=prop(mat,8)  !8 Iz bending rigidity in z direction
            dx= xord0(j1) - xord0(i1)
            dy= yord0(j1) - yord0(i1)
            dz= zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
!
!           orientation
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl

!           Calculate the stiffness matrix and assemble
            xl9= xl0
            call elmstfFRM_D(xl9,zix0,ziy0,ziz0,a0,e0,g0,ek12, nELt)
!           use trans3d to transform from local to global stifness
            call trans3d_D(xll,xmm,xnn,ek12,b0)
!
!           expand to [18x18]

            ek(1:18,1:18)=0.0d0

            do    j=1,12
            do    k=1,12
                ek(j,k)=ek12(j,k)
            enddo
            enddo

            call assembCOL(maxstiff,nEQ,stf,ek,i1,j1,k1,nloc)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   FORM GEOMetric stiffness matrices  [KGx],[KGy],[KGxy]
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine formgeom_s(geo,ele,xord0,yord0,zord0,xord,yord,zord,prop,nloc,nND,nEL,nEQ,nMT,maxstiff,geoFRM)
    implicit none
    integer:: nND,nEL,nEQ,nMT,maxstiff
    integer:: ele(nEL,5), nloc(nEQ)
    real(8):: prop(nMT,10),geoFRM(nEL)
    real(8):: eg12(12,12),eg(18,18)
!
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)

    real(8):: geo(maxstiff)

    integer:: n,i1,j1,k1,mat,nELt
    real(8):: e0,g0,a0,b0,dx,dy,dz,xl0,xl,xll,xmm,xnn,sxx,s,xl9

!   initialize [G] to zero
    geo(1:maxstiff)=0.0d0

    !rewind(igeoFRM)
    !rewind(igeoPLT)

    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        mat = ele(n,5)
!
        if (nELt == 2) then
!           constit relation not change
            e0=prop(mat,1)
            g0=prop(mat,2)
            a0=prop(mat,3)
            b0=prop(mat,5)
            dx= xord0(j1) - xord0(i1)
            dy= yord0(j1) - yord0(i1)
            dz= zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
!           orientation
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl
!
!           Calculate the stiffness matrix and assemble
            !read(igeoFRM,*) sxx
            sxx=geoFRM(n)
            s=sxx
            xl9= xl0
            call elmgeomFRM_D(xl9,eg12,s)
!           use trans3d to transform from local to global stifness
            call trans3d_D(xll,xmm,xnn,eg12,b0)

!           expand to [18x18]
            eg(1:18,1:18)=0.0d0
            eg(1:12,1:12)=eg12(1:12,1:12)
            call assembCOL(maxstiff,nEQ,geo,eg ,i1,j1,k1,nloc)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ELeMent MASs matrix for the FRaMe
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine elmmasFRM_D(rho,area,length,zix,em,alphaf)
    implicit none
    real(8):: rho, area, length,zix
    real(8):: em(12,12)
    real(8):: alphaf,roal

    em(1:12,1:12) = 0.0d0
    roal = rho*area*length/2.0d0
    em(1,1)     = roal
    em(2,2)     = roal
    em(3,3)     = roal
    em(4,4)     = roal*zix/area
    em(5,5)     = roal*length*length*alphaf/48d0
    em(6,6)     = roal*length*length*alphaf/48d0
    em(7,7)     = em(1,1)
    em(8,8)     = em(2,2)
    em(9,9)     = em(3,3)
    em(10,10)   = em(4,4)
    em(11,11)   = em(5,5)
    em(12,12)   = em(6,6)
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ELeMent STiFfness for FRaMe
!   calculates the element stiffness matrices. Page 70
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   input: length element length
!   emod E modulus
!   gmod shear modulus
!   area cross section area
!   ix J torsion rigidity of the cross section
!   iz Iz bending rigidity in z direction
!   iy Iy bending rigidity in y direction
    subroutine elmstfFRM_D(length, ix, iy, iz,area, emod, gmod, ek ,nELt)
    implicit none
    real(8):: area, length, ix, iy, iz, emod, gmod
    real(8):: ek(12,12)
    integer:: nELt

    integer:: i,j
    real(8):: emlen,emlen2,emlen3
!
!   initialize all ek elements to zero
    ek(1:12,1:12)=0.0d0
!
!   STIFFNESS matrix in local coordinates
!
    emlen  = emod/length
    emlen2 = emlen/length
    emlen3 = emlen2/length
    if  (nELt == 2) then
        ek(1,1)   =   area*emlen
        ek(2,2)   =   12.0d0*emlen3*iz
        ek(3,3)   =   12.0d0*emlen3*iy
        ek(4,4)   =   gmod*ix/length
        ek(5,5)   =   4.0d0*emlen*iy
        ek(6,6)   =   4.0d0*emlen*iz
!
        ek(2,6)   =   6.0d0*emlen2*iz
        ek(3,5)   =  -6.0d0*emlen2*iy
    else
        write(*,*)'not this nELt:',nELt
        stop
    endif
!
    ek(7,7)   =   ek(1,1)
    ek(8,8)   =   ek(2,2)
    ek(9,9)   =   ek(3,3)
    ek(10,10) =   ek(4,4)
    ek(11,11) =   ek(5,5)
    ek(12,12) =   ek(6,6)
!
    ek(1,7)   =   -ek(1,1)
    ek(2,8)   =   -ek(2,2)
    ek(2,12)  =    ek(2,6)
    ek(3,9)   =   -ek(3,3)
    ek(3,11)  =    ek(3,5)
    ek(4,10)  =   -ek(4,4)
    ek(5,9)   =   -ek(3,5)
    ek(5,11)  =    ek(5,5)/2.0d0
    ek(6,8)   =   -ek(2,6)
    ek(6,12)  =    ek(6,6)/2.0d0
!
    ek(8,12)  =   -ek(2,6)
    ek(9,11)  =   -ek(3,5)
!
!
!   impose the symmetry
    do  i= 1, 12
    do  j= i, 12
        ek(j,i) = ek(i,j)
    enddo
    enddo
!
    return
    end
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ELeMent GEOMetric stiffness matrix for a FRaMe
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! s tension force
    ! element initial length
    subroutine elmgeomFRM_D(length,eg,s)
    implicit none
    real(8):: length,eg(12,12),s
    real(8):: emlenz,alpha
    integer:: i,j
!
!   initialize all eg elements to zero
    eg(1:12,1:12)=0.0d0
!   Stiffness matrix in local coordinates
!** if (s .gt. 200.) s=200.
!   emlenz  =   s/(30.0*length)
!***alpha   =   s*1.0e-0
!   beta    =   1.0
    alpha   =   s*1.0d-3
    emlenz  =   s/length


    eg(1,1)   =   alpha
    eg(2,2)   =   emlenz
    eg(3,3)   =   emlenz
    eg(4,4)   =   alpha
    eg(5,5)   =   0.0d0
    eg(6,6)   =   0.0d0

!
    eg(7,7)   =   eg(1,1)
    eg(8,8)   =   eg(2,2)
    eg(9,9)   =   eg(3,3)
    eg(10,10) =   eg(4,4)
    eg(11,11) =   eg(5,5)
    eg(12,12) =   eg(6,6)
!
    eg(1,7)   =   -eg(1,1)
    eg(2,8)   =   -eg(2,2)
    eg(2,12)  =    eg(2,6)
    eg(3,9)   =   -eg(3,3)
    eg(3,11)  =    eg(3,5)
    eg(4,10)  =   -eg(4,4)
    eg(5,9)   =   -eg(3,5)
    eg(5,11)  =   -eg(5,5)/4.0d0
    eg(6,8)   =   -eg(2,6)
    eg(6,12)  =   -eg(6,6)/4.0d0
!
    eg(8,12)  =   -eg(2,6)
    eg(9,11)  =   -eg(3,5)
!
!   impose the symmetry
    do  i= 1, 12
    do  j= i, 12
        eg(j,i) = eg(i,j)
    enddo
    enddo
!
!   check diagonal terms
    do    i=1,12
        if (dabs(eg(i,i)) .lt. 1.0d-20) eg(i,i)=1.0d-20
    enddo
!
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ASSEMBle element matrices in COLumn form
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine assembCOL(maxstiff,nEQ,aa,a,i1,j1,k1,nloc)
    implicit none
    integer:: maxstiff,nEQ
    real(8):: aa(maxstiff),a(18,18)
    integer:: idof(18)
    integer:: i1,j1,k1,nloc(nEQ)
    integer:: i,j,imax,jmax,ieqn1,ieqn2,jband,iloc
    imax = 0
!
!   Set idof to posible DoF number of each nodes
    if    (j1 .gt. 0  .AND. k1 .gt. 0) then
        imax=18
        jmax=18
    elseif (j1 .eq. 0  .AND. k1 .eq. 0) then
        imax=6
        jmax=6
    endif

    do  i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
!
!   Store the values for individual array in global array
    do  i= 1, imax
        ieqn1 = idof(i)
        do    j= i, jmax
            ieqn2 = idof(j)
            if    (ieqn1 .gt. ieqn2) then
                    jband= (ieqn1-ieqn2)+1
                    iloc = nloc(ieqn1)
                    aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                   aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
            else
                    jband= (ieqn2-ieqn1)+1
                    iloc = nloc(ieqn2)
                    aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                   aa(ieqn1,jband) = aa(ieqn1,jband) + a(i,j)
            endif
        enddo
    enddo
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ASSEMBle LUMped element matrices
    subroutine assembLUM(nEQ,aa,a,i1,j1,k1)
    implicit none
    integer:: nEQ
    real(8):: aa(nEQ),a(18,18)
    integer:: i1,j1,k1
    integer:: idof(18),i,ieqn1
!
!   Set idof to posible DoF number of each nodes
    do  i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
!
!   Store the values for individual array in global array
    do  i= 1, 18
        ieqn1 = idof(i)
        aa(ieqn1) = aa(ieqn1) + a(i,i)
    enddo
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ASSEMBle consistent FORce for pressurized flex plate
!   input: a local value; nEQ length of aa, number of dof (6 * elements);
!   output:: aa global value
    subroutine assembFOR(nEQ, aa, a,i1, j1,k1)
    implicit none
    integer:: nEQ
    real(8):: aa(nEQ  ),a(18)
    integer:: i1,j1,k1
    integer:: idof(18),i,ieqn1
!
!   Set idof to posible DoF number of each nodes
    do    i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
!
!   Store the values for individual array in global array
    do  i= 1, 18
        ieqn1 = idof(i)
        aa(ieqn1) = aa(ieqn1) + a(i)
    enddo
!
    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   TRANSformation of matrix in 3D. L->G
    subroutine trans3d_D(l,m,n,ek,beta)
    implicit none
    real(8):: ek(12,12),rt(3,3),r(3,3),ktemp(12,12)
    real(8):: m,n,l,beta,pi,sb,cb,d
    integer:: i,j,k,j1,j2,ii,jj,in,jn

    pi=4.0d0*datan(1.0d0)
!
    sb=dsin(beta*pi/180.d0)
    cb=dcos(beta*pi/180.d0)
    d=dsqrt(1.0d0-n**2)
!   if (dabs(l).ge. 0.995 .and. dabs(beta).le. 0.01) return
    if (dabs(n).gt.0.995d0) then
           r(1,1)  =  0.0d0
           r(1,2)  =  0.0d0
           r(1,3)  =  n
           r(2,1)  = -n*sb
           r(2,2)  =  cb
           r(2,3)  =  0.0d0
           r(3,1)  = -n*cb
           r(3,2)  = -sb
           r(3,3)  =  0.0d0
    else
           r(1,1)  =  l
           r(1,2)  =  m
           r(1,3)  =  n
           if (dabs(beta) .le. .01d0) then
              r(2,1)  =  -m/d
              r(2,2)  =  l/d
              r(2,3)  =  0.0d0
              r(3,1)  =  -l*n/d
              r(3,2)  =  -m*n/d
              r(3,3)  =  d
           else
              r(2,1)  =  -(m*cb+l*n*sb)/d
              r(2,2)  =  (l*cb-m*n*sb)/d
              r(2,3)  =  d*sb
              r(3,1)  =  (m*sb-l*n*cb)/d
              r(3,2)  =  -(l*sb+m*n*cb)/d
              r(3,3)  =  d*cb
           endif
    endif
!
    do  in=1,3
    do  jn=1,3
        rt(jn,in)=r(in,jn)
    enddo
    enddo
!    take [Rtrans][K][R] using the nature of [R] to speed computation.
!    k is sectioned off into 3x3s then multiplied [rtrans][k][r]
!
    do  i=0,3
    do  j=0,3
        j1=i*3
        j2=j*3
        do    k=1,3
        do    ii=1,3
            ktemp(j1+k,j2+ii)=0.0d0
            do     jj=1,3
                ktemp(j1+k,j2+ii)=ktemp(j1+k,j2+ii)+ek(j1+k,j2+jj)*r(jj,ii)
            enddo
        enddo
        enddo
        do  k=1,3
        do  ii=1,3
            ek(j1+k,j2+ii)=0.0d0
            do  jj=1,3
                ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rt(k,jj)*ktemp(j1+jj,j2+ii)
            enddo
        enddo
        enddo
    enddo
    enddo

    return
    end


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   AxB product for COLumn storage
    subroutine AxBCOL(matrix,maxstiff,vecin,vecout, nEQ,nprof,nprof2,nloc)
    implicit none
    integer:: maxstiff,nEQ
    integer:: nprof(nEQ),nprof2(nEQ),nloc(nEQ)
    real(8):: vecout(nEQ),matrix(maxstiff),vecin(nEQ)
    real(8):: val,valmat
    integer:: i,j,io,is,jlim,iloc
!
    do    i=1,nEQ
        jlim=max(1,(i-nprof(i)+1))
        do  j=jlim,i
            is=i
            io=i-j+1
            if (io .gt. nprof(is)) cycle
            iloc = nloc(is) + io -1
            valmat=matrix(iloc)
            val = vecin(j)
            vecout(i)=vecout(i) + val*valmat
        enddo

        jlim=min(nprof2(i),(nEQ-i+1))
        do  j=2,jlim
            is=i+j-1
            io=j
            if (io .gt. nprof(is)) cycle
            iloc = nloc(is) + io -1
            valmat=matrix(iloc)
            val = vecin(i+j-1)
            vecout(i)=vecout(i) + val*valmat
        enddo
    enddo


    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   UDU decomposition of COLumn profiled system
    subroutine uduCOL_D(a,maxstore,nEQ,imult,nprof,nloc)
    implicit none
    integer:: maxstore,nEQ,imult
    real(8):: a(maxstore)
    integer:: nprof(nEQ), nloc(nEQ)

    real(8):: temp,sum
    integer:: i,j,k,j2,j3,im1,jm1,is,io,iloc,jcol,iloci,ilocj,iloc1

    if (dabs(a(1)) .le. 1.0d-20) then
        imult=0
        return
    endif
!
    if    (nEQ .eq. 1) then
        imult=1
        return
    endif
!

    do  j=2,nEQ
        jm1=j-1
        j2=j-nprof(j)+1
        if    (j2.lt.1) then
            j2=1
        endif
!
!       off-diagonal terms
        if    (jm1.eq.1) then
            is=j
            io=1
            iloc = nloc(is) + io - 1
            sum=a(iloc)
!           sum=a(j,1)
        else
            do  i=j2+1,jm1
                im1=i-1
                is=j
                io=j-i+1
                iloc = nloc(is) + io - 1
                sum=a(iloc   )
!               sum=a(i,j-i+1)
!
                j3=i-nprof(i)+1
                jcol=j3
                if    (j3 .lt. j2) then
                    jcol=j2
                endif
!               do    k=j2,im1
                do  k=jcol,im1
                    is=i
                    io=i-k+1
                    iloci = nloc(is) + io - 1
                    is=j
                    io=j-k+1
                    ilocj = nloc(is) + io - 1
                    sum=sum-a(iloci  )*a(ilocj  )
!                   sum=sum-a(k,i-k+1)*a(k,j-k+1)
                    imult=imult+1
                enddo
                a(iloc   )=sum
!               a(i,j-i+1)=sum
            enddo
            is=j
            io=1
            iloc = nloc(is) + io - 1
            sum=a(iloc   )
!           sum=a(j,1)
        endif
!
!        diagonal terms
        do  k=j2,jm1
            is=j
            io=j-k+1
            ilocj = nloc(is) + io - 1
            is=k
            io=1
            iloc1 = nloc(is) + io - 1
            temp=a(ilocj  )/a(iloc1)
            sum=sum-temp*a(ilocj  )
            a(ilocj  )=temp
!           temp=a(k,j-k+1)/a(k,1)
!           sum=sum-temp*a(k,j-k+1)
!           a(k,j-k+1)=temp
            imult=imult+2
        enddo

        is=j
        io=1
        iloc = nloc(is) + io - 1
        a(iloc   ) = sum
        if (dabs(sum) .le. 1.0d-20) then
            imult=0
            return
        endif
    enddo
!
!
    return
    end

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   BAcK solver of COLumn profiled system
!   input a
!   output wk (du)
    subroutine bakCOL_D(a,maxstore,b,nEQ,nBD,wk,imult,nprof,nprof2,nloc)
    implicit none
    integer:: maxstore,nEQ,nBD,imult
    real(8):: a(maxstore),b(nEQ),wk(nEQ)
    integer:: nprof(nEQ), nprof2(nEQ),nloc(nEQ)

    real(8):: sum
    integer:: i,j,k,i1,k2,jb,jbb,km1,is,io,iloc

!
!   forward substitutions
    do    i=1,nEQ
!        j=i-nBD+1
        j=i-nprof(i)+1
        if    (i.le.nprof(i) ) then
            j=1
        endif
        jb=i-nBD+1
        jbb=jb
        if    (i.le.nBD    ) then
            jbb=1
        endif
        sum=b(i)
        km1=i-1
        if    (j.gt.km1) then
            wk(i)=sum
        else
            do  k=j,km1
                is=i
                io=i-k+1
                iloc = nloc(is) + io - 1
                sum=sum-a(iloc   )*wk(k)
!               sum=sum-a(k,i-k+1)*wk(k)
                imult=imult+1
            enddo
            wk(i)=sum
        endif
    enddo
!
!   middle terms
    do  i=1,nEQ
        is=i
        io=1
        iloc = nloc(is) + io - 1
        wk(i)=wk(i)/a(iloc )
!       wk(i)=wk(i)/a(i,1)
        imult=imult+1
    enddo
!
!   backward substitution
    do  i1=1,nEQ
        i=nEQ-i1+1
        j=i+nprof2(i) -1
        if    (j.gt.nEQ) then
            j=nEQ
        endif
        jb=i+nBD-1
        jbb=jb
        if    (jb.gt.nEQ) then
            jbb=nEQ
        endif
        sum=wk(i)
        k2=i+1
        if    (k2.gt.j) then
            wk(i)=sum
        else
            do  k=k2,j
                is=k
                io=k-i+1
                if (io .gt. nprof(is)) cycle
                iloc = nloc(is) + io - 1
                sum=sum-a(iloc   )*wk(k)
!               sum=sum-a(i,k-i+1)*wk(k)
                imult=imult+1
            enddo

            wk(i)=sum
        endif
    enddo
!
    return
    end
end module SolidSolver