module SolidBody
    use ConstParams
    use SolidSolver
    implicit none
    private
    ! Immersed boundary method parameters
    public :: m_nFish
    integer:: m_nFish,m_ntolLBM
    real(8):: m_dtolLBM, m_IBPenaltyAlpha, m_denIn, m_uvwIn(3), m_Aref, m_Eref, m_Fref, m_Lref, m_Pref, m_Tref, m_Uref
    integer:: m_boundaryConditions(1:6)
    ! nFish     number of bodies
    ! ntolLBM maximum number of iterations for IB force calculation
    ! dtolIBM   tolerance for IB force calculation
    ! Pbeta     coefficient in penalty force calculation
    public :: VBodies,read_solid_files,allocate_solid_memory,Initialise_solid_bodies,Write_solid_v_bodies,FSInteraction_force, &
              Calculate_Solid_params,Solver,Write_solid_cont,Read_solid_cont,write_solid_field,Write_solid_Check, write_solid_Information, &
              calculate_reference_params,set_solidbody_parameters
    type :: VirtualBody
        type(BeamSolver):: rbm
        ! vitural body in which fluid ID
        integer :: v_carrierFluidId
        !!!virtual infomation
        !!!virtual body surface
        integer :: v_nelmts
        integer :: v_type
        ! v_type
        ! -1 body given by surface mesh
        ! 1 plate given by line mesh
        ! 2 rod given by line mesh
        integer :: v_move = 0 ! 0 stationary 1 moving
        integer :: count_Area = 0, count_Interp = 0
        real(8), allocatable :: v_Exyz0(:, :) ! initial element center (x, y, z), required for type 3
        real(8), allocatable :: v_Exyz(:, :) ! element center (x, y, z)
        real(8), allocatable :: v_Ea(:) ! element area
        !area center with equal weight on both sides
        real(8), allocatable :: v_Evel(:, :)
        real(8), allocatable :: v_Eforce(:, :) ! element center (x, y, z)
        integer(2), allocatable :: v_Ei(:, :) ! element stencial integer index [ix-1,ix,ix1,ix2, iy-1,iy,iy1,iy2, iz-1,iz,iz1,iz2]
        real(4), allocatable :: v_Ew(:, :) ! element stential weight [wx-1, wx, wx1, wx2, wy-1, wy, wy1, wy2, wz-1, wz, wz1, wz2]
        !calculated using central linear and angular velocities
        integer,allocatable :: vtor(:)! of size fake_npts
        integer,allocatable :: rtov(:)! of size real_npts+1
        
    contains
        procedure :: Initialise => Initialise_
        procedure :: PlateBuild => PlateBuild_
        procedure :: SurfaceBuild => SurfaceBuild_
        procedure :: UpdatePosVelArea => UpdatePosVelArea_
        procedure :: PlateUpdatePosVelArea => PlateUpdatePosVelArea_
        procedure :: SurfaceBuildPosVelArea => SurfaceBuildPosVelArea_
        procedure :: SurfaceUpdatePosVel => SurfaceUpdatePosVel_
        procedure :: Write_body => Write_body_
        procedure :: PlateWrite_body => PlateWrite_body_
        procedure :: SurfaceWrite_body => SurfaceWrite_body_
        procedure :: UpdateElmtInterp => UpdateElmtInterp_
        procedure :: PenaltyForce => PenaltyForce_
        procedure :: FluidVolumeForce => FluidVolumeForce_
        procedure :: vtor_f => vtor_
    end type VirtualBody
    type(VirtualBody), allocatable :: VBodies(:)

    contains

    SUBROUTINE read_solid_files(filename)
        ! read global body parameters
        implicit none
        character(LEN=40),intent(in):: filename
        character(LEN=256):: buffer
        real(8):: alphaf,NewmarkGamma,NewmarkBeta,dampK,dampM,dtolFEM
        integer:: nfishGroup,isKB,ntolFEM
        integer:: iFish,ifishGroup,numX,numY,numZ
        character(LEN=40) :: t_FEmeshName,keywordstr
        integer:: t_iBodyModel,t_iBodyType,t_isMotionGiven(6)
        real(8):: t_denR,t_psR,t_EmR,t_tcR,t_KB,t_KS,t_St
        real(8):: t_freq,firstXYZ(1:3),deltaXYZ(1:3)
        real(8):: t_initXYZVel(3),t_XYZAmpl(3),t_XYZPhi(3),t_AoAo(3),t_AoAAmpl(3),t_AoAPhi(3)
        integer:: order1=0,order2=0,order3=0,lineX,lineY,lineZ
        character(LEN=40),allocatable:: FEmeshName(:)
        integer,allocatable:: fishNum(:)
        integer,allocatable:: iBodyModel(:),iBodyType(:),isMotionGiven(:,:)
        real(8),allocatable:: denR(:),psR(:),EmR(:),tcR(:),KB(:),KS(:)
        real(8),allocatable:: initXYZVel(:,:),XYZo(:,:),XYZAmpl(:,:),XYZPhi(:,:),freq(:),St(:)
        real(8),allocatable:: AoAo(:,:),AoAAmpl(:,:),AoAPhi(:,:)
        ! read body parameters from inflow file
        open(unit=111, file=filename, status='old', action='read')
        keywordstr = 'SolidBody'
        call found_keyword(111,keywordstr)
        call readNextData(111, buffer)
        read(buffer,*)    m_IBPenaltyAlpha,alphaf
        call readNextData(111, buffer)
        read(buffer,*)    NewmarkGamma,NewmarkBeta
        call readNextData(111, buffer)
        read(buffer,*)    dampK,dampM
        call readNextData(111, buffer)
        read(buffer,*)    dtolFEM,ntolFEM
        call readNextData(111, buffer)
        read(buffer,*)    m_nFish,nfishGroup,isKB
        if(m_IBPenaltyAlpha.le.1.d-6) then
            write(*,*) 'ERROR: IBPenaltyalpha should be positive (default 1)', m_IBPenaltyAlpha
            stop
        endif
        ! set solid solver global parameters
        call Set_SolidSolver_Params(dampK,dampM,NewmarkGamma,NewmarkBeta,alphaf,dtolFEM,ntolFEM,isKB)
        allocate(FEmeshName(m_nFish),fishNum(nfishGroup+1),iBodyModel(m_nFish),iBodyType(m_nFish),isMotionGiven(6,m_nFish), &
                denR(m_nFish),psR(m_nFish),EmR(m_nFish),tcR(m_nFish),KB(m_nFish),KS(m_nFish), &
                initXYZVel(3,m_nFish),XYZo(3,m_nFish),XYZAmpl(3,m_nFish),XYZPhi(3,m_nFish),freq(m_nFish),St(m_nFish), &
                AoAo(3,m_nFish),AoAAmpl(3,m_nFish),AoAPhi(3,m_nFish))
        ! read fish parameters for each type
        fishNum(1)=1
        do ifishGroup = 1,nfishGroup
            call readNextData(111, buffer)
            read(buffer,*)    fishNum(ifishGroup+1),numX,numY,numZ
            call readNextData(111, buffer)
            read(buffer,*)    t_FEmeshName
            call readNextData(111, buffer)
            read(buffer,*)    t_iBodyModel,t_iBodyType
            call readNextData(111, buffer)
            read(buffer,*)    t_isMotionGiven(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    t_isMotionGiven(4:6)
            call readNextData(111, buffer)
            read(buffer,*)    t_denR,t_psR
            call readNextData(111, buffer)
            if(isKB==0) then
                read(buffer,*)    t_EmR,t_tcR
            else
                read(buffer,*)    t_KB,t_KS
            endif
            call readNextData(111, buffer)
            read(buffer,*)    t_freq,t_St
            call readNextData(111, buffer)
            read(buffer,*)    firstXYZ(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    deltaXYZ(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    t_initXYZVel(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    t_XYZAmpl(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    t_XYZPhi(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    t_AoAo(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    t_AoAAmpl(1:3)
            call readNextData(111, buffer)
            read(buffer,*)    t_AoAPhi(1:3)
            if(ifishGroup .lt. nfishGroup) call readequal(111)
            order1 = order1 + fishNum(ifishGroup  );
            order2 = order2 + fishNum(ifishGroup+1);
            ! read parameters for each fish
            do iFish=order1,order2
                FEmeshName(iFish) = t_FEmeshName
                iBodyModel(iFish) = t_iBodyModel
                iBodyType (iFish) = t_iBodyType
                isMotionGiven(1:6,iFish)=t_isMotionGiven(1:6)
                denR(iFish)= t_denR
                psR(iFish) = t_psR
                if(isKB==0) then
                    EmR(iFish) = t_EmR
                    tcR(iFish) = t_tcR
                elseif(isKB==1) then
                    KB(iFish)  = t_KB
                    KS(iFish)  = t_KS
                endif
                freq(iFish) = t_freq
                St(iFish) = t_St
                initXYZVel(1:3,iFish)= t_initXYZVel(1:3)
                XYZAmpl(1:3,iFish) = t_XYZAmpl(1:3)
                XYZPhi(1:3,iFish)  = t_XYZPhi(1:3)
                AoAo(1:3,iFish)    = t_AoAo(1:3)
                AoAAmpl(1:3,iFish) = t_AoAAmpl(1:3)
                AoAPhi(1:3,iFish)  = t_AoAPhi(1:3)
                ! calculate the initial location for each fish
                order3 = iFish - order1
                lineX  = mod(order3,numX)
                lineY  = (order3 - lineX)/numX
                lineZ  = (order3 - mod(order3,numX*numY))/numX*numY
                XYZo(1,iFish) = firstXYZ(1) + deltaXYZ(1) * lineX
                XYZo(2,iFish) = firstXYZ(2) + deltaXYZ(2) * lineY
                XYZo(3,iFish) = firstXYZ(3) + deltaXYZ(3) * lineZ
            enddo
        enddo
        close(111)
        ! allocate bodies memory
        allocate(VBodies(m_nFish))
        do iFish = 1,m_nFish
            VBodies(iFish)%v_type = iBodyType(iFish)
            if (iBodyType(iFish).eq.-1) then
                call SurfacetoBeam_write(FEmeshName(iFish))
            endif
            call VBodies(iFish)%rbm%SetSolver(FEmeshName(iFish),&
                iBodyModel(iFish),isMotionGiven(1:6,iFish), &
                denR(iFish),KB(iFish),KS(iFish),EmR(iFish),psR(iFish),tcR(iFish),St(iFish), &
                freq(iFish),initXYZVel(1:3,iFish),XYZo(1:3,iFish),XYZAmpl(1:3,iFish),XYZPhi(1:3,iFish), &
                AoAo(1:3,iFish),AoAAmpl(1:3,iFish),AoAPhi(1:3,iFish))
        enddo
    end subroutine read_solid_files


    SUBROUTINE calculate_reference_params(flow)
        USE ConstParams
        USE FlowCondition, only: FlowCondType
        implicit none
        type(FlowCondType),intent(inout):: flow
        integer:: iFish
        real(8):: nUref(1:m_nFish)
        ! reference length
        if(m_nFish.eq.0) then
            flow%Lref = 1.d0
        else
            flow%Lref  = flow%Lchod
        endif
        ! reference velocity
        if(flow%UrefType==0) then
            flow%Uref = dabs(flow%uvwIn(1))
        elseif(flow%UrefType==1) then
            flow%Uref = dabs(flow%uvwIn(2))
        elseif(flow%UrefType==2) then
            flow%Uref = dabs(flow%uvwIn(3))
        elseif(flow%UrefType==3) then
            flow%Uref = dsqrt(flow%uvwIn(1)**2 + flow%uvwIn(2)**2 + flow%uvwIn(3)**2)
        elseif(flow%UrefType==4) then
            if (flow%velocityKind .eq. 2) then
                flow%Uref = dabs(flow%shearRateIn(1))  ! flow%shearRateIn(1) is Velocity Amplitude
            else
                write(*,*) 'oscillatory flow must set velocityKind to 2'
                stop
            endif
        elseif(flow%UrefType==5) then
            flow%Uref = flow%Lref * MAXVAL(VBodies(:)%rbm%Freq)
        elseif(flow%UrefType==6) then
            do iFish=1,m_nFish
            nUref(iFish)=2.d0*pi*VBodies(iFish)%rbm%Freq*MAXVAL(dabs(VBodies(iFish)%rbm%xyzAmpl(1:3)))
            enddo
            flow%Uref = MAXVAL(nUref(:))
        elseif(flow%UrefType==7) then
            do iFish=1,m_nFish
            nUref(iFish)=2.d0*pi*VBodies(iFish)%rbm%Freq*MAXVAL(dabs(VBodies(iFish)%rbm%xyzAmpl(1:3)))*2.D0 !Park 2017 pof
            enddo
            flow%Uref = MAXVAL(nUref(1:m_nFish))
        else
            write(*,*) 'use input reference velocity'
        endif
        ! reference time
        if(flow%TrefType==0) then
            flow%Tref = flow%Lref / flow%Uref
        elseif(flow%TrefType==1) then
            flow%Tref = 1 / maxval(VBodies(:)%rbm%Freq)
        else
            write(*,*) 'use input reference time'
        endif
        ! reference acceleration, force, energy, power
        flow%Aref = flow%Uref/flow%Tref
        flow%Fref = 0.5*flow%denIn*flow%Uref**2*flow%Asfac
        flow%Eref = 0.5*flow%denIn*flow%Uref**2*flow%Asfac*flow%Lref
        flow%Pref = 0.5*flow%denIn*flow%Uref**2*flow%Asfac*flow%Uref
        ! fluid viscosity
        flow%nu =  flow%Uref*flow%Lref/flow%Re
        flow%Mu =  flow%nu*flow%denIn
    END SUBROUTINE

    subroutine allocate_solid_memory(Asfac,Lchod,Lspan,AR)
        real(8),intent(out):: Asfac,Lchod,Lspan,AR
        real(8) :: nAsfac(m_nFish),nLchod(m_nFish)
        integer :: iFish,maxN
        write(*,'(A)') '========================================================='
        do iFish = 1,m_nFish
            if (dabs(maxval(VBodies(iFish)%rbm%initXYZVel(1:3))-0.d0) .gt. 1e-5 .or. &
                dabs(maxval(VBodies(iFish)%rbm%XYZAmpl(1:3))-0.d0)    .gt. 1e-5 .or. &
                dabs(maxval(VBodies(iFish)%rbm%AoAAmpl(1:3))-0.d0)    .gt. 1e-5) then
                VBodies(iFish)%v_move = 1
            endif
            call VBodies(iFish)%rbm%Allocate_solid(nAsfac(iFish),nLchod(iFish))
            write(*,*)'read FEMeshFile ',iFish,' end' 
        enddo
        if (m_nFish .gt. 0) then
            !Use the object with the largest area as the reference object
            maxN  = maxloc(nAsfac, dim=1)
            Asfac = nAsfac(maxN)
            Lchod = nLchod(maxN)
            if (VBodies(maxN)%v_type .eq. 1) then
                Lspan = sum(VBodies(maxN)%rbm%r_Lspan(:)+VBodies(maxN)%rbm%r_Rspan(:))/dble(VBodies(maxN)%rbm%nND)
            else
                Lspan = maxval(VBodies(maxN)%rbm%xyzful00(:,3))-minval(VBodies(maxN)%rbm%xyzful00(:,3))
            endif
            if((Lchod-1.0d0)<=1.0d-2)Lchod=1.0d0
            if((Lspan-1.0d0)<=1.0d-2)Lspan=1.0d0
            if (VBodies(maxN)%v_type .eq. 1) then
                AR    = Lspan**2/Asfac
            else
                AR    = 1.d0
            endif
        else
            Asfac = 0.d0
            Lchod = 0.d0
            Lspan = 0.d0
            AR = 0.d0
        endif
    end subroutine allocate_solid_memory

    subroutine set_solidbody_parameters(denIn,uvwIn,BndConds,&
        Aref,Eref,Fref,Lref,Pref,Tref,Uref,ntolLBM,dtolLBM)
        implicit none
        real(8),intent(in):: denIn,uvwIn(1:3),Aref,Eref,Fref,Lref,Pref,Tref,Uref,dtolLBM
        integer,intent(in):: BndConds(1:6),ntolLBM
        real(8):: Lthck,uMax
        ! set gobal parameters
        m_denIn = denIn
        m_uvwIn = uvwIn
        m_boundaryConditions(1:6) = BndConds(1:6)
        m_Aref = Aref
        m_Eref = Eref
        m_Fref = Fref
        m_Lref = Lref
        m_Pref = Pref
        m_Tref = Tref
        m_Uref = Uref
        m_ntolLBM = ntolLBM
        m_dtolLBM = dtolLBM
        call Calculate_Solid_params(uMax,Lthck)
    end subroutine set_solidbody_parameters

    subroutine Initialise_solid_bodies(time,g)
        implicit none
        real(8),intent(in):: time,g(3)
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%Initialise(time,g)
            call VBodies(iFish)%Initialise()
        enddo
    end subroutine Initialise_solid_bodies

    subroutine Initialise_(this)
        ! read beam central line file and allocate memory
        implicit none
        class(VirtualBody), intent(inout) :: this
        if (this%v_type .eq. 1) then
            call this%PlateBuild()
        elseif (this%v_type .eq. -1) then
            call this%SurfaceBuild()
        else
            write(*,*) 'not implemented body type', this%v_type
            stop
        endif
    end subroutine Initialise_

    subroutine Calculate_Solid_params(uMax,Lthck)
        implicit none
        real(8),intent(out):: Lthck,uMax
        integer:: iFish
        real(8):: nLthck(m_nFish)
        uMax = 0.d0
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%calculate_angle_material(m_Lref, m_Uref, m_denIn, uMax, m_uvwIn, nLthck(iFish))
        enddo
        Lthck = maxval(nLthck)
    end subroutine Calculate_Solid_params

    SUBROUTINE Solver(bodies,time,isubstep,deltat,subdeltat)
        implicit none
        integer,intent(in):: bodies(0:m_nFish)
        integer,intent(in):: isubstep
        real(8),intent(in):: time,deltat,subdeltat
        integer:: iFish, i
        !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iFish)
        do i = 1,bodies(0)
            iFish = bodies(i)
            call VBodies(iFish)%rbm%structure(iFish,time,isubstep,deltat,subdeltat)
        enddo !do iFish=1,nFish
        !$OMP END PARALLEL DO
    END SUBROUTINE

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Write_solid_v_bodies(time)
        implicit none
        real(8),intent(in):: time
        integer :: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%Write_body(iFish,time)
        enddo
    end subroutine Write_solid_v_bodies

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write check point file for restarting simulation
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Write_solid_cont(fid)
        implicit none
        integer,intent(in):: fid
        integer:: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%write_solid_temp(fid)
        enddo
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    read check point file for restarting simulation
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Read_solid_cont(fid)
        implicit none
        integer,intent(in):: fid
        integer:: iFish
        do iFish = 1,m_nFish
            call VBodies(iFish)%rbm%read_solid_temp(fid)
        enddo
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write structure field, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_solid_field(time)
        implicit none
        real(8),intent(in):: time
        real(8):: timeTref
    !   -------------------------------------------------------
        integer:: i,iFish
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName
        !==================================================================================================
        timeTref = time/m_Tref
        write(fileName,'(I10)') nint(timeTref*1d5)
        fileName = adjustr(fileName)
        DO  i=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        END DO
    
        do iFish=1,m_nFish
            call VBodies(iFish)%rbm%write_solid(m_Lref,m_Uref,m_Aref,m_Fref,iFish,fileName)
        enddo
    !   =============================================
    end subroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    write solid parameters for checking, tecplot ASCII format
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Write_solid_Check(filename)
        implicit none
        character(len=40):: filename
        character(len=4):: IDstr
        integer:: iFish
        open(111,file=filename,position='append')
        do iFish=1,m_nFish
            write(IDstr,'(I4.4)')iFish
            write(111,'(A,A,A  )')'============================= nFish = ',IDstr,' =============================='
            write(IDstr,'(I4.4)')VBodies(iFish)%v_carrierFluidId
            write(111,'(A,A    )') 'inWhichBlock : ', IDstr
            write(111,'(A,A    )')'---------------------------------------------------------------------------'
            call VBodies(iFish)%rbm%write_solid_params(111)
            call VBodies(iFish)%rbm%write_solid_materials(111)
        enddo
        write(111,'(A      )')'===================================================================='
        close(111)
    end subroutine

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   write solid data
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_solid_Information(time,timeOutInfo,Asfac,solidProbingNum,solidProbingNode)
        implicit none
        integer,intent(in):: solidProbingNum,solidProbingNode(solidProbingNum)
        real(8),intent(in):: time,timeOutInfo,Asfac
        integer:: iFish,fid = 111
        do iFish=1,m_nFish
            call VBodies(iFish)%rbm%write_solid_info(fid,iFish,time,timeOutInfo,m_Tref,m_Lref,m_Uref,m_Aref,m_Fref,m_Pref,m_Eref,Asfac)
            call VBodies(iFish)%rbm%write_solid_probes(fid,iFish,time,solidProbingNum,solidProbingNode(1:solidProbingNum),m_Tref,m_Lref,m_Uref,m_Aref)
        enddo
    end subroutine

    subroutine FSInteraction_force(bodies,dt,dh,xmin,ymin,zmin,xDim,yDim,zDim,uuu,force)
        implicit none
        integer,intent(in):: bodies(0:m_nFish)
        real(8),intent(in):: dt,dh,xmin,ymin,zmin
        integer,intent(in):: xDim,yDim,zDim
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        integer :: i,iFish
        do i = 1,bodies(0)
            iFish = bodies(i)
            call VBodies(iFish)%UpdatePosVelArea()
        enddo
        call calculate_interaction_force(bodies,dt,dh,xmin,ymin,zmin,xDim,yDim,zDim,uuu,force)
    end subroutine

    subroutine PlateUpdatePosVelArea_(this)
        !   compute displacement, velocity, area at surface element center
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        integer :: i,s,cnt,i1,i2
        real(8) :: tmpxyz(3), tmpvel(3), tmpdx(3), dirc(3)
        real(8) :: dh, left, len, dl, ls, area, dirc_norm, IBPenaltyBeta
        IBPenaltyBeta = - m_IBPenaltyalpha* 2.0d0*m_denIn
        do i = 1,this%rbm%nEL
            i1 = this%rbm%ele(i,1)
            i2 = this%rbm%ele(i,2)
            tmpxyz = 0.5d0 * (this%rbm%xyzful(i1,1:3) + this%rbm%xyzful(i2,1:3))
            tmpvel = 0.5d0 * (this%rbm%velful(i1,1:3) + this%rbm%velful(i2,1:3))
            left = 0.5d0 * (this%rbm%r_Lspan(i1) + this%rbm%r_Lspan(i2))
            len = 0.5d0 * (this%rbm%r_Rspan(i1) + this%rbm%r_Rspan(i2)) + left
            dl = len / dble(this%rbm%r_Nspan(i))
            tmpdx = this%rbm%xyzful(i2,1:3) - this%rbm%xyzful(i1,1:3)
            dh = dsqrt(tmpdx(1)*tmpdx(1)+tmpdx(2)*tmpdx(2)+tmpdx(3)*tmpdx(3))
            area = dl * dh * IBPenaltyBeta
            dirc(1:3) = this%rbm%r_dirc(i1,1:3) + this%rbm%r_dirc(i2,1:3)
            dirc_norm = dsqrt(sum(dirc**2))
            if (dirc_norm .gt. 1e-10) then
                dirc = dirc / dirc_norm
            else
                write(*,*) this%rbm%r_dirc(i1,1:3), "and", this%rbm%r_dirc(i2,1:3), "are opposite directions; no unique bisector exists."
                stop
            endif
            cnt = this%rtov(i) - 1
            do s=1,this%rbm%r_Nspan(i)
                ls = dl * (0.5d0 + dble(s-1)) - left
                this%v_Exyz(1:3,cnt+s) = tmpxyz + dirc*ls
                this%v_Evel(1:3,cnt+s) = tmpvel
                this%v_Ea(cnt+s) = area
            enddo
        enddo
    endsubroutine PlateUpdatePosVelArea_

    subroutine SurfaceBuildPosVelArea_(this,Surfacetmpnpts,Surfacetmpnelmts,Surfacetmpxyz,Surfacetmpele)
        !   compute displacement, velocity, area at surface element center
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        integer:: Surfacetmpnpts, Surfacetmpnelmts
        real(8) :: Surfacetmpxyz(3,Surfacetmpnpts)
        integer :: Surfacetmpele(3,Surfacetmpnelmts)
        integer :: i,i1,i2,i3
        real(8) :: A(3),B(3),C(3),tmparea,IBPenaltyBeta
        allocate(this%v_Exyz0(3,this%v_nelmts))
        IBPenaltyBeta = - m_IBPenaltyalpha* 2.0d0*m_denIn
        do i = 1,Surfacetmpnelmts
            i1 = Surfacetmpele(1,i)
            i2 = Surfacetmpele(2,i)
            i3 = Surfacetmpele(3,i)
            A = Surfacetmpxyz(1:3,i1)
            B = Surfacetmpxyz(1:3,i2)
            C = Surfacetmpxyz(1:3,i3)
            call cpt_incenter(this%v_Exyz0(1:3,i))
            call cpt_area(tmparea)
            this%v_Ea(i) = tmparea*IBPenaltyBeta
        enddo
        call this%SurfaceUpdatePosVel()
        contains
        subroutine cpt_incenter(Exyz)
            implicit none
            real(8), intent(out) :: Exyz(3)
            real(8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3
            real(8) :: la, lb, lc, invC
                x1=A(1)
                x2=B(1)
                x3=C(1)
                y1=A(2)
                y2=B(2)
                y3=C(2)
                z1=A(3)
                z2=B(3)
                z3=C(3)
                la = sqrt((x2 - x3)**2 + (y2 - y3)**2 + (z2 - z3)**2)
                lb = sqrt((x3 - x1)**2 + (y3 - y1)**2 + (z3 - z1)**2)
                lc = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
                invC = 1 / (la + lb + lc)
                Exyz(1) = (la*x1 + lb*x2 + lc*x3) * invC
                Exyz(2) = (la*y1 + lb*y2 + lc*y3) * invC
                Exyz(3) = (la*z1 + lb*z2 + lc*z3) * invC
        end subroutine cpt_incenter
        subroutine cpt_area(area)
            implicit none
            real(8), intent(out) :: area
            real(8) :: x1, x2, x3, y1, y2, y3, z1, z2, z3, ax, ay, az
                x1=A(1)
                x2=B(1)
                x3=C(1)
                y1=A(2)
                y2=B(2)
                y3=C(2)
                z1=A(3)
                z2=B(3)
                z3=C(3)
                ax =((z1-z2)*(y3-y2) + (y2-y1)*(z3-z2))
                ay =((x1-x2)*(z3-z2) + (z2-z1)*(x3-x2))
                az =((y1-y2)*(x3-x2) + (x2-x1)*(y3-y2))
                area=dsqrt( ax*ax + ay*ay + az*az)*0.5d0
        end subroutine cpt_area
    end subroutine SurfaceBuildPosVelArea_

    subroutine SurfaceUpdatePosVel_(this)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        integer:: i
        do  i=1,this%v_nelmts
            this%v_Exyz(1:3,i)=matmul(this%rbm%TTTnxt(1:3,1:3),this%v_Exyz0(1:3,i))+this%rbm%XYZ(1:3)
            this%v_Evel(1:3,i)=[this%rbm%WWW3(2)*this%v_Exyz(3,i)-this%rbm%WWW3(3)*this%v_Exyz(2,i),    &
                                this%rbm%WWW3(3)*this%v_Exyz(1,i)-this%rbm%WWW3(1)*this%v_Exyz(3,i),    &
                                this%rbm%WWW3(1)*this%v_Exyz(2,i)-this%rbm%WWW3(2)*this%v_Exyz(1,i)    ]&
                                + this%rbm%UVW(1:3) + this%rbm%initXYZVel(1:3)
        enddo
    end subroutine SurfaceUpdatePosVel_

    subroutine UpdatePosVelArea_(this)
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        if (this%v_type .eq. 1 .and. (this%v_move .eq. 1 .or. this%rbm%iBodyModel .eq. 2 .or. this%count_Area .eq. 0)) then
            call this%PlateUpdatePosVelArea()
            this%count_Area = 1
        elseif (this%v_type .eq. -1 .and. (this%v_move .eq. 1 .or. this%rbm%iBodyModel .eq. 2 .or. this%count_Area .eq. 0)) then
            call this%SurfaceUpdatePosVel()
            this%count_Area = 1
        else
        endif
    end subroutine UpdatePosVelArea_

    FUNCTION vtor_(this,x)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer, intent(in) :: x
        integer :: vtor_
        vtor_ = 0
        if (this%v_type .eq. 1)then
            vtor_ = this%vtor(x)
        elseif (this%v_type .eq. -1)then
            vtor_ = 1
        endif
    ENDFUNCTION vtor_
    FUNCTION rtov_(this,x)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer, intent(in) :: x
        integer :: rtov_
        rtov_ = this%rtov(x)
    ENDFUNCTION rtov_

    subroutine UpdateElmtInterp_(this,dh,xmin,ymin,zmin,xDim,yDim,zDim)
        use omp_lib
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: dh,xmin,ymin,zmin
        integer,intent(in):: xDim,yDim,zDim
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2)
        real(8)::x0,y0,z0,detx,dety,detz,invdh
        integer::i0,j0,k0,iEL,x,y,z,i,j,k
        !==================================================================================================
        invdh = 1.D0/dh
        i0 = floor((this%v_Exyz(1,1) - xmin) * invdh)
        x0 = xmin + dble(i0) * dh
        i0 = i0 + 1
        j0 = floor((this%v_Exyz(2,1) - ymin) * invdh)
        y0 = ymin + dble(j0) * dh
        j0 = j0 + 1
        k0 = floor((this%v_Exyz(3,1) - zmin) * invdh)
        z0 = zmin + dble(k0) * dh
        k0 = k0 + 1
        ! compute the velocity of IB nodes at element center
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,i,j,k,x,y,z,rx,ry,rz,detx,dety,detz,ix,jy,kz)
        do  iEL=1,this%v_nelmts
            call minloc_fast(this%v_Exyz(1,iEL), x0, i0, invdh, i, detx)
            call minloc_fast(this%v_Exyz(2,iEL), y0, j0, invdh, j, dety)
            call minloc_fast(this%v_Exyz(3,iEL), z0, k0, invdh, k, detz)
            call trimedindex(i, xDim, ix, m_boundaryConditions(1:2))
            call trimedindex(j, yDim, jy, m_boundaryConditions(3:4))
            call trimedindex(k, zDim, kz, m_boundaryConditions(5:6))
            this%v_Ei(1:4,iEL) = ix
            this%v_Ei(5:8,iEL) = jy
            this%v_Ei(9:12,iEL) = kz
            do x=-1,2
                rx(x)=Phi(dble(x)-detx)
            enddo
            do y=-1,2
                ry(y)=Phi(dble(y)-dety)
            enddo
            do z=-1,2
                rz(z)=Phi(dble(z)-detz)
            enddo
            this%v_Ew(1:4,iEL) = rx
            this%v_Ew(5:8,iEL) = ry
            this%v_Ew(9:12,iEL) = rz
        enddo
        !$OMP END PARALLEL DO

        contains
        ! return the array(index) <= x < array(index+1)
        ! assum uniform grid around x(x=i0) = x0
        SUBROUTINE minloc_fast(x_, x0_, i0_, invdh_, index_, offset_)
            implicit none
            integer, intent(in):: i0_
            real(8), intent(in):: x_, x0_, invdh_
            integer, intent(out):: index_
            real(8), intent(out):: offset_
            offset_ = (x_ - x0_)*invdh_
            index_ = floor(offset_)
            offset_ = offset_ - dble(index_)
            index_ = index_ + i0_
        END SUBROUTINE
        FUNCTION Phi(x_)
            IMPLICIT NONE
            real(8)::Phi,x_,r
            r=dabs(x_)
            if(r<1.0d0)then
                Phi=(3.d0-2.d0*r+dsqrt( 1.d0+4.d0*r*(1.d0-r)))*0.125d0
            elseif(r<2.0d0)then
                Phi=(5.d0-2.d0*r-dsqrt(-7.d0+4.d0*r*(3.d0-r)))*0.125d0
            else
                Phi=0.0d0
            endif
        ENDFUNCTION Phi
        SUBROUTINE trimedindex(i_, xDim_, ix_, boundaryConditions_)
            implicit none
            integer, intent(in):: boundaryConditions_(1:2)
            integer, intent(in):: i_, xDim_
            integer, intent(out):: ix_(-1:2)
            integer:: k_
            do k_=-1,2
                ix_(k_) = i_ + k_
                if (ix_(k_)<1) then
                    if(boundaryConditions_(1).eq.BCPeriodic) then
                        ix_(k_) = ix_(k_) + xDim_
                    else if((boundaryConditions_(1).eq.BCSymmetric .or. boundaryConditions_(1).eq.BCstationary_Wall) .and. ix_(k_).eq.0) then
                        ix_(k_) = 2
                    else if((boundaryConditions_(1).eq.BCstationary_Wall_halfway) .and. ix_(k_).eq.0) then
                        ix_(k_) = 1
                    else
                        write(*,*) 'index out of xmin bound', ix_(k_)
                        stop
                    endif
                else if(ix_(k_)>xDim_) then
                    if(boundaryConditions_(2).eq.BCPeriodic) then
                        ix_(k_) = ix_(k_) - xDim_
                    else if((boundaryConditions_(2).eq.BCSymmetric .or. boundaryConditions_(2).eq.BCstationary_Wall) .and. ix_(k_).eq.xDim_+1) then
                        ix_(k_) = xDim_ - 1
                    else if((boundaryConditions_(2).eq.BCstationary_Wall_halfway) .and. ix_(k_).eq.xDim_+1) then
                        ix_(k_) = xDim_
                    else
                        write(*,*) 'index out of xmax bound', ix_(k_)
                        stop
                    endif
                endif
            enddo
        END SUBROUTINE trimedindex
    end subroutine UpdateElmtInterp_

    SUBROUTINE calculate_interaction_force(bodies,dt,dh,xmin,ymin,zmin,xDim,yDim,zDim,uuu,force)
        ! calculate elements interaction force using IB method
        IMPLICIT NONE
        integer,intent(in):: bodies(0:m_nFish)
        real(8),intent(in):: dt,dh,xmin,ymin,zmin
        integer,intent(in):: xDim,yDim,zDim
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        !================================
        integer:: iFish, i
        integer:: iterLBM
        real(8):: dmaxLBM,dsum
        real(8)::tol,ntol
        ! update virtual body shape and velocity
        do i = 1, bodies(0)
            iFish = bodies(i)
            if (VBodies(iFish)%v_move .eq. 1 .or. VBodies(iFish)%rbm%iBodyModel .eq. 2 .or. VBodies(iFish)%count_Interp .eq. 0 ) then
                call VBodies(iFish)%UpdateElmtInterp(dh,xmin,ymin,zmin,xDim,yDim,zDim)
                VBodies(iFish)%count_Interp = 1
            endif
            VBodies(iFish)%v_Eforce = 0.0d0
        enddo
        if (bodies(0) .gt. 0) then
            ! calculate interaction force using immersed-boundary method
            iterLBM=0
            dmaxLBM=1d10
            do  while( iterLBM<m_ntolLBM .and. dmaxLBM>m_dtolLBM)
                dmaxLBM = 0.d0
                dsum=0.0d0
                do i = 1, bodies(0)
                    iFish = bodies(i)
                    call VBodies(iFish)%PenaltyForce(dt,dh,xDim,yDim,zDim,tol,ntol,uuu)
                    dmaxLBM = dmaxLBM + tol
                    dsum = dsum + ntol
                enddo
                dmaxLBM=dmaxLBM/(dsum * m_Uref)
                iterLBM=iterLBM+1
            enddo
        endif
        ! update body load and fluid force
        do i = 1, bodies(0)
            iFish = bodies(i)
            VBodies(iFish)%rbm%extful = 0.0d0
            ! to do, consider gravity
        enddo
        do i = 1, bodies(0)
            iFish = bodies(i)
            call VBodies(iFish)%FluidVolumeForce(dh,xDim,yDim,zDim,force)
        enddo
    END SUBROUTINE

    SUBROUTINE FluidVolumeForce_(this,dh,xDim,yDim,zDim,force)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        class(VirtualBody), intent(inout) :: this
        real(8),intent(in):: dh
        integer,intent(in):: xDim,yDim,zDim
        real(8),intent(out)::force(zDim,yDim,xDim,1:3)
        !==================================================================================================
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2),forcetemp(1:3)
        real(8):: forceElemTemp(3),invh3
        !==================================================================================================
        integer::x,y,z,iEL,i1,i2
        !==================================================================================================
        invh3 = (1.d0/dh)**3
        ! compute the velocity of IB nodes at element center
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(1:4,iEL)
            jy = this%v_Ei(5:8,iEL)
            kz = this%v_Ei(9:12,iEL)
            rx = this%v_Ew(1:4,iEL)
            ry = this%v_Ew(5:8,iEL)
            rz = this%v_Ew(9:12,iEL)
            forceElemTemp = this%v_Eforce(1:3,iEL)
            ! update beam load, momentum is not included
            i1=this%rbm%ele(this%vtor_f(iEL),1)
            i2=this%rbm%ele(this%vtor_f(iEL),2)
            this%rbm%extful(i1,1:3) = this%rbm%extful(i1,1:3) + 0.5d0 * forceElemTemp
            this%rbm%extful(i2,1:3) = this%rbm%extful(i2,1:3) + 0.5d0 * forceElemTemp
            forceElemTemp(1:3) = forceElemTemp(1:3) * invh3
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        forceTemp(1:3) = -forceElemTemp(1:3)*rx(x)*ry(y)*rz(z)
                        ! add flow body force
                        force(kz(z),jy(y),ix(x),1:3) = force(kz(z),jy(y),ix(x),1:3) + forceTemp(1:3)
                    enddo
                enddo
            enddo
        enddo
    END SUBROUTINE FluidVolumeForce_

    SUBROUTINE PenaltyForce_(this,dt,dh,xDim,yDim,zDim,tolerance,ntolsum,uuu)
        USE, INTRINSIC :: IEEE_ARITHMETIC
        IMPLICIT NONE
        real(8),intent(in):: dt,dh
        integer,intent(in):: xDim,yDim,zDim
        real(8),intent(inout)::uuu(zDim,yDim,xDim,1:3)
        class(VirtualBody), intent(inout) :: this
        real(8),intent(out)::tolerance, ntolsum
        !======================================================
        integer:: ix(-1:2),jy(-1:2),kz(-1:2)
        real(8):: rx(-1:2),ry(-1:2),rz(-1:2),forcetemp(1:3)
        real(8):: velElemIB(3),forceElemTemp(this%v_nelmts,3),invh3
        !======================================================
        integer:: x,y,z,iEL
        !======================================================
        invh3 = 0.5d0*dt*(1.d0/dh)**3/m_denIn
        tolerance = 0.d0
        ntolsum = dble(this%v_nelmts)
        ! compute the velocity of IB nodes at element center
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iEL,x,y,z,rx,ry,rz,ix,jy,kz,velElemIB,forceTemp) reduction(+:tolerance)
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(1:4,iEL)
            jy = this%v_Ei(5:8,iEL)
            kz = this%v_Ei(9:12,iEL)
            rx = this%v_Ew(1:4,iEL)
            ry = this%v_Ew(5:8,iEL)
            rz = this%v_Ew(9:12,iEL)
            velElemIB(1:3)=0.0d0
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        velElemIB(1:3)=velElemIB(1:3)+uuu(kz(z),jy(y),ix(x),1:3)*rx(x)*ry(y)*rz(z)
                    enddo
                enddo
            enddo
            velElemIB = this%v_Evel(1:3,iEL)-velElemIB(1:3)
            forceTemp(1:3) = velElemIB(1:3)*this%v_Ea(iEL)
            !if ((.not. IEEE_IS_FINITE(forceTemp(1))) .or. (.not. IEEE_IS_FINITE(forceTemp(2))) .or. (.not. IEEE_IS_FINITE(forceTemp(3)))) then
            !    write(*, *) 'Nan found in forceElemTemp', forceTemp
            !    write(*, *) 'Nan found at (ix, jy, kz)', ix(0), jy(0), kz(0)
            !    stop
            !endif
            tolerance = tolerance + dabs(velElemIB(1)) + dabs(velElemIB(2)) + dabs(velElemIB(3))
            this%v_Eforce(1:3,iEL) = this%v_Eforce(1:3,iEL) + forceTemp
            forceElemTemp(iEL,1:3) = forceTemp * invh3
        enddo
        !$OMP END PARALLEL DO
        if (.not. IEEE_IS_FINITE(tolerance)) then
            write(*, *) 'Nan found in PenaltyForce'
            stop
        endif

        ! correct velocity
        do  iEL=1,this%v_nelmts
            ix = this%v_Ei(1:4,iEL)
            jy = this%v_Ei(5:8,iEL)
            kz = this%v_Ei(9:12,iEL)
            rx = this%v_Ew(1:4,iEL)
            ry = this%v_Ew(5:8,iEL)
            rz = this%v_Ew(9:12,iEL)
            do x=-1,2
                do y=-1,2
                    do z=-1,2
                        uuu(kz(z),jy(y),ix(x),1:3) = uuu(kz(z),jy(y),ix(x),1:3) - forceElemTemp(iEL,1:3)*rx(x)*ry(y)*rz(z)
                    enddo
                enddo
            enddo
        enddo
    END SUBROUTINE PenaltyForce_

    subroutine PlateBuild_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer:: i, j
        allocate(this%rtov(this%rbm%nEL+1))
        this%v_nelmts = 0
        do i=1,this%rbm%nEL
            this%rtov(i) = this%v_nelmts + 1
            this%v_nelmts = this%v_nelmts + this%rbm%r_Nspan(i)
        enddo
        this%rtov(this%rbm%nEL+1) = this%v_nelmts + 1
        allocate(this%vtor(this%v_nelmts))
        do i=1,this%rbm%nEL
            do j= this%rtov(i), this%rtov(i+1)-1
                this%vtor(j) = i
            enddo
        enddo
        allocate(this%v_Exyz(3,this%v_nelmts), this%v_Ea(this%v_nelmts), this%v_Eforce(3,this%v_nelmts))
        allocate(this%v_Evel(3,this%v_nelmts), this%v_Ei(12,this%v_nelmts), this%v_Ew(12,this%v_nelmts))
        call this%PlateUpdatePosVelArea()
    end subroutine PlateBuild_

    subroutine SurfaceBuild_(this)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer:: Surfacetmpnpts, Surfacetmpnelmts
        real(8),allocatable :: Surfacetmpxyz(:,:)
        integer,allocatable :: Surfacetmpele(:,:)
        call Read_gmsh(this%rbm%FEmeshName,Surfacetmpnpts,Surfacetmpnelmts,Surfacetmpxyz,Surfacetmpele)
        this%v_nelmts = Surfacetmpnelmts
        allocate(this%v_Exyz(3,this%v_nelmts), this%v_Ea(this%v_nelmts), this%v_Eforce(3,this%v_nelmts))
        allocate(this%v_Evel(3,this%v_nelmts), this%v_Ei(12,this%v_nelmts), this%v_Ew(12,this%v_nelmts))
        call this%SurfaceBuildPosVelArea(Surfacetmpnpts,Surfacetmpnelmts,Surfacetmpxyz,Surfacetmpele)
        deallocate(Surfacetmpxyz,Surfacetmpele)
    end subroutine SurfaceBuild_

    subroutine Read_gmsh(FEmeshName,Surfacetmpnpts,Surfacetmpnelmts,Surfacetmpxyz,Surfacetmpele)
        implicit none
        character (LEN=40),intent(in):: FEmeshName
        character (LEN=40):: gmshName
        integer:: Surfacetmpnpts, Surfacetmpnelmts
        real(8),allocatable :: Surfacetmpxyz(:,:)
        integer,allocatable :: Surfacetmpele(:,:)
        integer :: fileiD = 111, num, i, temp_prop(4)
        i = index(FEmeshName, '.')
        gmshName = FEmeshName(:i) // 'msh'
        open(unit=fileiD, file = trim(adjustl(gmshName)) )! read *.msh file
            ! read nodes
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*) num
            Surfacetmpnpts = num
            allocate(Surfacetmpxyz(3,Surfacetmpnpts))
            do i = 1,Surfacetmpnpts
                read(fileiD,*)num,Surfacetmpxyz(1,i),Surfacetmpxyz(2,i),Surfacetmpxyz(3,i)
            enddo
            ! read element
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*) num
            Surfacetmpnelmts = num
            do i = 1,Surfacetmpnelmts
                read(fileiD,*)num,temp_prop(1:4)
                if((temp_prop(1)==2) .and. (temp_prop(2)==2) .and. (temp_prop(3)==0)) then
                    Surfacetmpnelmts = Surfacetmpnelmts-num+1
                    exit
                endif
            enddo
            backspace(fileiD)
            allocate(Surfacetmpele(3,Surfacetmpnelmts))
            do i = 1,Surfacetmpnelmts
                read(fileiD,*)num,temp_prop(1:4),Surfacetmpele(1,i),Surfacetmpele(2,i),Surfacetmpele(3,i)
            enddo
        close(fileiD)
    endsubroutine Read_gmsh

    subroutine SurfacetoBeam_write(FEmeshName)
        implicit none
        character (LEN=40),intent(inout):: FEmeshName
        integer :: fileiD = 111, i, num
        real(8) :: Surfacetmpxyz(3,3)
        i = index(FEmeshName, '.')
        FEmeshName = FEmeshName(:i) // 'msh'
        open(unit=fileiD, file = trim(adjustl(FEmeshName)) )
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*)
            read(fileiD,*)
            do i = 1,3
                read(fileiD,*)num,Surfacetmpxyz(1,i),Surfacetmpxyz(2,i),Surfacetmpxyz(3,i)
            enddo
        close(fileiD)
        i = index(FEmeshName, '.')
        FEmeshName = FEmeshName(:i) // 'dat'
        open(unit=fileiD, file = trim(adjustl(FEmeshName)))! write *.dat file
            write(fileiD,*) "Frame3D(This is a .dat file converted from .msh file)"
        close(fileiD)
        open(unit=fileiD, file = trim(adjustl(FEmeshName)),position='append')! write *.dat file
            write(fileiD,*) "     3     1     1"
            write(fileiD,*) "END"
            write(fileiD,*) "     3"
            do i = 1,3
                write(fileiD,*) i,Surfacetmpxyz(1,i),Surfacetmpxyz(2,i),Surfacetmpxyz(3,i),"   0.0   0.0   0. 0. 1."
            enddo
            write(fileiD,*) "END"
            write(fileiD,*) "     1     I     J     K  TYPE   MAT   LEN"
            write(fileiD,*) "     1     1     2     3     3     1     0"
            write(fileiD,*) "END"
            write(fileiD,*) "     1  XTRA  YTRA  ZTRA  XROT  YROT  ZROT"
            write(fileiD,*) "     1     1     0     0     0     0     0"
            write(fileiD,*) "END"
            write(fileiD,*) "     1   E           G           A           RHO         GAMMA       IP          IA          IB"
            write(fileiD,*) "     1   0.100D+01   0.100D+01   0.100D+01   0.100D+01   0.000D+00   0.100D+01   0.150D+01   0.500D+00"
            write(fileiD,*) "END"
        close(fileiD)
    endsubroutine SurfacetoBeam_write

    subroutine Write_body_(this,iFish,time)
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer,intent(in) :: iFish
        real(8),intent(in) :: time
        if (this%v_type .eq. 1) then
            call this%PlateWrite_body(iFish,time)
        elseif (this%v_type .eq. -1) then
            call this%SurfaceWrite_body(iFish,time)
        else
            write(*, *) "body type not implemented", this%v_type
        endif
    endsubroutine Write_body_

    subroutine PlateWrite_body_(this,iFish,time)
        ! to do: generate a temporary mesh
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer,intent(in) :: iFish
        real(8),intent(in) :: time
        !   -------------------------------------------------------
        real(8):: timeTref
        integer:: i, i1, i2
        real(8) :: tmpxyz(3)
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName,idstr
        integer,parameter:: idfile=100
        !==========================================================================
        timeTref = time/m_Tref
        write(fileName,'(I10)') nint(timeTref*1d5)
        fileName = adjustr(fileName)
        do  I=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo
        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        !==========================================================================
        open(idfile, FILE='./DatBodySpan/BodyFake'//trim(idstr)//'_'//trim(filename)//'.dat')
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        write(idfile, '(A,I7,A,I7,A)') 'ZONE N=',2*this%rbm%nND,', E=',this%rbm%nEL,', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'

        do i = 1,this%rbm%nNd
            tmpxyz = this%rbm%xyzful(i,1:3)
            write(idfile, *) (tmpxyz - this%rbm%r_Lspan(i) * this%rbm%r_dirc(i,1:3))/m_Lref
            write(idfile, *) (tmpxyz + this%rbm%r_Rspan(i) * this%rbm%r_dirc(i,1:3))/m_Lref
        enddo
        do  i=1,this%rbm%nEL
            i1 = this%rbm%ele(i,1)
            i2 = this%rbm%ele(i,2)
            write(idfile, *) 2*i1 -1, 2*i1,2*i2,2*i2-1
        enddo
        close(idfile)
    end subroutine PlateWrite_body_

    subroutine SurfaceWrite_body_(this,iFish,time)
        ! to do: generate a temporary mesh
        implicit none
        class(VirtualBody), intent(inout) :: this
        integer,intent(in) :: iFish
        real(8),intent(in) :: time
        !   -------------------------------------------------------
        real(8):: timeTref
        integer:: i, i1, i2, i3
        real(8) :: tmpxyz(3)
        integer,parameter::nameLen=10
        character (LEN=nameLen):: fileName,idstr
        integer,parameter:: idfile=100
        integer:: Surfacetmpnpts, Surfacetmpnelmts
        real(8),allocatable :: Surfacetmpxyz(:,:)
        integer,allocatable :: Surfacetmpele(:,:)
        !==========================================================================
        timeTref = time/m_Tref
        write(fileName,'(I10)') nint(timeTref*1d5)
        fileName = adjustr(fileName)
        do  I=1,nameLen
            if(fileName(i:i)==' ')fileName(i:i)='0'
        enddo
        write(idstr, '(I3.3)') iFish ! assume iFish < 1000
        !==========================================================================
        open(idfile, FILE='./DatBodySpan/BodyFake'//trim(idstr)//'_'//trim(filename)//'.dat')
        call Read_gmsh(this%rbm%FEmeshName,Surfacetmpnpts,Surfacetmpnelmts,Surfacetmpxyz,Surfacetmpele)
        do  i=1,Surfacetmpnpts
            Surfacetmpxyz(1:3,i)=matmul(this%rbm%TTTnxt(1:3,1:3),Surfacetmpxyz(1:3,i))+this%rbm%XYZ(1:3)
        enddo
        write(idfile, '(A)') 'variables = "x" "y" "z"'
        write(idfile, '(A,I7,A,I7,A)') 'ZONE N=',Surfacetmpnpts,', E=',Surfacetmpnelmts,', DATAPACKING=POINT, ZONETYPE=FETRIANGLE'
        do i = 1,Surfacetmpnpts
            tmpxyz = Surfacetmpxyz(1:3,i)
            write(idfile, *) tmpxyz/m_Lref
        enddo
        do  i=1,Surfacetmpnelmts
            i1 = Surfacetmpele(1,i)
            i2 = Surfacetmpele(2,i)
            i3 = Surfacetmpele(3,i)
            write(idfile, *) i1, i2, i3
        enddo
        deallocate(Surfacetmpxyz,Surfacetmpele)
        close(idfile)
    end subroutine SurfaceWrite_body_

end module SolidBody
