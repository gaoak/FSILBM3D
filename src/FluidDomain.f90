module FluidDomain
    use ConstParams
    use FlowCondition
    implicit none
    private
    integer :: m_nblock
    type :: LBMBlock
        integer :: npsize
        integer :: ID,iCollidModel,offsetOutput
        integer :: xDim,yDim,zDim
        real(8) :: dh,xmin,ymin,zmin,xmax,ymax,zmax
        integer :: BndConds(1:6)
        real(8) :: params(1:10),Omega,Omega2 ! single time, two time relaxation
        real(8) :: M_COLLID(0:lbmDim,0:lbmDim),M_FORCE(0:lbmDim,0:lbmDim) ! multiple time relaxation
        integer, allocatable :: OMPpartition(:),OMPparindex(:),OMPeid(:)
        real(8), allocatable :: OMPedge(:,:,:)
        real(8), allocatable :: fIn(:,:,:,:),uuu(:,:,:,:),force(:,:,:,:),den(:,:,:)
        real(4), allocatable :: OUTutmp(:,:,:),OUTvtmp(:,:,:),OUTwtmp(:,:,:)
        real(4) :: offsetMoveGrid(1:3),volumeForce(3)
    contains
        procedure :: allocate_fluid => allocate_fluid_
        procedure :: Initialise => Initialise_
        procedure :: calculate_macro_quantities => calculate_macro_quantities_
        procedure :: collision => collision_
        procedure :: streaming => streaming_
        procedure :: ComputeFieldStat => ComputeFieldStat_
        procedure :: write_flow => write_flow_
        procedure :: write_continue => write_continue_
        procedure :: read_continue => read_continue_
        procedure :: update_volumn_force => update_volumn_force_
    end type LBMBlock
    type(LBMBlock), allocatable :: LBMblks(:)

    contains

    SUBROUTINE read_fuild_blocks(filename)
        implicit none
        character(LEN=40),intent(in):: filename
        integer:: iblock
        character(LEN=256):: buffer
        ! creat LBMblks and read parameters
        open(unit=111, file=filename, status='old', action='read')
        call found_keyword(111,'FluidDomain')
        call readNextData(111, buffer)
        read(buffer,*) m_nblock
        allocate(LBMblks(m_nblock))
        do iblock = 1,m_nblock
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%npsize
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%ID,LBMblks(iblock)%iCollidModel,LBMblks(iblock)%offsetOutput
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%xDim,LBMblks(iblock)%yDim,LBMblks(iblock)%zDim
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%dh,LBMblks(iblock)%xmin,LBMblks(iblock)%ymin,LBMblks(iblock)%zmin
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%BndConds(1:6)
            call readNextData(111, buffer)
            read(buffer,*)    LBMblks(iblock)%params(1:10)
            ! calculate the maximum coordinates
            LBMblks(iblock)%xmax = LBMblks(iblock)%xmin + LBMblks(iblock)%dh*(LBMblks(iblock)%xDim-1)
            LBMblks(iblock)%ymax = LBMblks(iblock)%ymin + LBMblks(iblock)%dh*(LBMblks(iblock)%yDim-1)
            LBMblks(iblock)%zmax = LBMblks(iblock)%zmin + LBMblks(iblock)%dh*(LBMblks(iblock)%zDim-1)
        enddo
        close(111)
    END SUBROUTINE

    !=================================================================================================
    ! Functions to iterate over each block
    SUBROUTINE allocate_fuild_memory_blocks()
        implicit none
        integer:: iblock
        do iblock = 1,m_nblock
            call LBMblks(iblock)%allocate_fluid()
        enddo
    END SUBROUTINE

    SUBROUTINE initialise_fuild_blocks(flow)
        implicit none
        type(FlowCondType),intent(inout) :: flow
        integer:: iblock
        do iblock = 1,m_nblock
            call LBMblks(iblock)%initialise(flow)
        enddo
    END SUBROUTINE

    SUBROUTINE read_continue_blocks(filename,step,time)
        implicit none
        character(LEN=40),intent(in):: filename
        integer,intent(out):: step
        real(8),intent(out):: time
        integer:: iblock
        do iblock = 1,m_nblock
            call LBMblks(iblock)%read_continue(filename,step,time)
        enddo
    END SUBROUTINE

    SUBROUTINE update_volumn_force_blocks(time)
        implicit none
        character(LEN=40),intent(in):: filename
        real(8),intent(in):: time
        integer:: iblock
        do iblock = 1,m_nblock
            call LBMblks(iblock)%update_volumn_force(time)
        enddo
    END SUBROUTINE

    SUBROUTINE calculate_macro_quantities_blocks()
        implicit none
        integer:: iblock
        do iblock = 1,m_nblock
            LBMblks(iblock)%calculate_macro_quantities()
        enddo
    END SUBROUTINE

    SUBROUTINE write_flow_blocks(time)
        implicit none
        integer:: iblock
        do iblock = 1,m_nblock
            LBMblks(iblock)%write_flow(time)
        enddo
    END SUBROUTINE

    !=================================================================================================
    ! Functions for single block
    SUBROUTINE allocate_fluid_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: xmin,ymin,zmin,xmax,ymax,zmax
        ! allocate fluid memory
        allocate(this%fIn(this%zDim,this%yDim,this%xDim,0:lbmDim))
        allocate(this%uuu(this%zDim,this%yDim,this%xDim,1:3))
        allocate(this%force(this%zDim,this%yDim,this%xDim,1:3))
        allocate(this%den(this%zDim,this%yDim,this%xDim))
        ! allocate output workspace
        xmin = 1 + this%offsetOutput
        ymin = 1 + this%offsetOutput
        zmin = 1 + this%offsetOutput
        xmax = this%xDim - this%offsetOutput
        ymax = this%yDim - this%offsetOutput
        zmax = this%zDim - this%offsetOutput
        allocate(this%oututmp(zmin:zmax,ymin:ymax,xmin:xmax))
        allocate(this%outvtmp(zmin:zmax,ymin:ymax,xmin:xmax))
        allocate(this%outwtmp(zmin:zmax,ymin:ymax,xmin:xmax))
        ! allocate mesh partition
        allocate(this%OMPpartition(1:this%npsize),this%OMPparindex(1:this%npsize+1),this%OMPeid(1:this%npsize))
        allocate(this%OMPedge(1:this%zDim,1:this%yDim, 1:this%npsize))
        call OMPPrePartition(this%xDim, this%npsize, this%OMPpartition, this%OMPparindex)

        contains

        SUBROUTINE OMPPrePartition(xDim, np, partition, parindex)
            implicit none
            integer:: np, xDim
            integer:: partition(1:np), parindex(1:np+1)
            integer:: psize, p, residual
            psize = xDim/np
            residual = xDim - psize * np
            parindex(1) = 1
            parindex(np+1) = xDim + 1
            do p=1,np
                if (p .gt. np-residual) then
                    partition(p) = psize + 1
                else
                    partition(p) = psize
                endif
                if (p .gt. 1) then
                    parindex(p) = parindex(p-1) + partition(p-1)
                endif
            enddo
        END SUBROUTINE OMPPrePartition
    END SUBROUTINE allocate_fluid_

    SUBROUTINE initialise_(this, flow)
        implicit none
        class(LBMBlock),intent(inout) :: this
        type(FlowCondType),intent(inout) :: flow
        ! select the collision model
        if(this%iCollidModel.eq.1) then
            call calculate_SRT_params()
        elseif(this%iCollidModel.eq.2) then
            call calculate_TRT_params(this%params(1))
        elseif(this%iCollidModel.eq.3) then
            call calculate_MRT_params()
        else
            write(*,*)' collision_step Model is not defined'
        endif
        ! initialize flow information
        call initialize_flow()

        contains

        SUBROUTINE calculate_SRT_params()
            implicit none
            real(8):: tau
            tau   =  flow%nu/(this%dh*Cs2)+0.5d0
            this%Omega =  1.0d0 / tau
        END SUBROUTINE calculate_SRT_params

        SUBROUTINE calculate_TRT_params(lambda)
            implicit none
            real(8):: tau, lambda
            tau   =  flow%nu/(this%dh*Cs2)+0.5d0
            this%Omega =  1.0d0 / tau
            this%Omega2 = 0.d0 ! to be updated
            lambda = 1.0d0/4.0d0;
        END SUBROUTINE calculate_TRT_params

        SUBROUTINE calculate_MRT_params()
            implicit none
            integer:: I
            real(8):: M_MRT(0:lbmDim,0:lbmDim),M_MRTI(0:lbmDim,0:lbmDim),M(0:lbmDim,0:lbmDim)
            real(8):: S_D(0:lbmDim,0:lbmDim),S(0:lbmDim)
            ! calculate MRTM transformation matrix
            DO  I=0,lbmDim
                M_MRT(0,I)=1
                M_MRT(1,I)=19*SUM(ee(I,1:3)**2)-30
                M_MRT(2,I)=(21*SUM(ee(I,1:3)**2)**2-53*SUM(ee(I,1:3)**2)+24)/2.0

                M_MRT(3,I)=ee(I,1)
                M_MRT(5,I)=ee(I,2)
                M_MRT(7,I)=ee(I,3)

                M_MRT(4,I)=(5*SUM(ee(I,1:3)**2)-9)*ee(I,1)
                M_MRT(6,I)=(5*SUM(ee(I,1:3)**2)-9)*ee(I,2)
                M_MRT(8,I)=(5*SUM(ee(I,1:3)**2)-9)*ee(I,3)

                M_MRT(9,I)=3*ee(I,1)**2-SUM(ee(I,1:3)**2)
                M_MRT(11,I)=ee(I,2)**2-ee(I,3)**2

                M_MRT(13,I)=ee(I,1)*ee(I,2)
                M_MRT(14,I)=ee(I,2)*ee(I,3)
                M_MRT(15,I)=ee(I,3)*ee(I,1)

                M_MRT(10,I)=(3*SUM(ee(I,1:3)**2)-5)*(3*ee(I,1)**2-SUM(ee(I,1:3)**2))
                M_MRT(12,I)=(3*SUM(ee(I,1:3)**2)-5)*(ee(I,2)**2-ee(I,3)**2)

                M_MRT(16,I)=(ee(I,2)**2-ee(I,3)**2)*ee(I,1)
                M_MRT(17,I)=(ee(I,3)**2-ee(I,1)**2)*ee(I,2)
                M_MRT(18,I)=(ee(I,1)**2-ee(I,2)**2)*ee(I,3)
            ENDDO
            ! calculate the inverse matrix
            M_MRTI=TRANSPOSE(M_MRT)
            M=MATMUL(M_MRT,M_MRTI)
            DO    I=0,lbmDim
                M_MRTI(0:lbmDim,I)=M_MRTI(0:lbmDim,I)/M(I,I)
            ENDDO
            ! S(0:lbmDim)=Omega ! restore to SRT if S is Omega
            !              0   1  2  3  4  5  6  7  8      9     10      11     12      13        14         15      16  17  18
            S(0:lbmDim)=[  s0,s1,s2,s0,s4,s0,s4,s0,s4,this%Omega,s10,this%Omega,s10,this%Omega,this%Omega,this%Omega,s16,s16,s16]
            ! calculate MRTM collision matrix
            ! IM*S*M
            S_D(0:lbmDim,0:lbmDim)=0.0D0
            DO    i=0,lbmDim
                S_D(i,i)=S(i)
            ENDDO
            this%M_COLLID=MATMUL(MATMUL(M_MRTI,S_D),M_MRT)
            ! calculate MRTM body-force matrix
            ! IM*(I-0.5D0*S)*M=I-0.5*IM*S*M
            S_D(0:lbmDim,0:lbmDim)=0.0D0
            DO    i=0,lbmDim
                S_D(i,i)=1.0d0
            ENDDO
            this%M_FORCE=S_D-0.5*this%M_COLLID
        END SUBROUTINE calculate_MRT_params

        SUBROUTINE initialize_flow()
            implicit none
            real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim)
            real(8):: xCoord,yCoord,zCoord
            integer:: x, y, z
            ! calculating initial flow velocity
            do  x = 1, this%xDim
                xCoord = this%xmin + this%dh * (x - 1);
            do  y = 1, this%yDim
                yCoord = this%ymin + this%dh * (x - 1);
            do  z = 1, this%zDim
                zCoord = this%zmin + this%dh * (x - 1);
                this%uuu(z, y, x, 1) = flow%uvwIn(1) + 0*flow%shearRateIn(1) + yCoord*flow%shearRateIn(2) + zCoord*flow%shearRateIn(3);
                this%uuu(z, y, x, 2) = flow%uvwIn(1) + xCoord*flow%shearRateIn(1) + 0*flow%shearRateIn(2) + zCoord*flow%shearRateIn(3);
                this%uuu(z, y, x, 3) = flow%uvwIn(1) + xCoord*flow%shearRateIn(1) + yCoord*flow%shearRateIn(2) + 0*flow%shearRateIn(3);
            enddo
            enddo
            enddo
            this%den(1:this%zDim,1:this%yDim,1:this%xDim) = flow%denIn
            ! calculating the distribution function
            do  x = 1, this%xDim
            do  y = 1, this%yDim
            do  z = 1, this%zDim
                uSqr           = sum(this%uuu(z,y,x,1:3)**2)
                uxyz(0:lbmDim) = this%uuu(z,y,x,1) * ee(0:lbmDim,1) + this%uuu(z,y,x,2) * ee(0:lbmDim,2)+this%uuu(z,y,x,3) * ee(0:lbmDim,3)
                fEq(0:lbmDim)  = wt(0:lbmDim) * this%den(z,y,x) * (1.0d0 + 3.0d0 * uxyz(0:lbmDim) + 4.5d0 * uxyz(0:lbmDim) * uxyz(0:lbmDim) - 1.5d0 * uSqr)
                
                this%fIn(z,y,x,0:lbmDim)=fEq(0:lbmDim)
            enddo
            enddo
            enddo
        END SUBROUTINE initialize_flow
    END SUBROUTINE initialise_

    !SUBROUTINE evaluateShearVelocity(xCoord,yCoord,zCoord,velocity)
    !    implicit none
    !    real(8):: xCoord,yCoord,zCoord,velocity(1:SpaceDim)
    !    velocity(1) = m_uuuIn(1) + 0 * m_shearRateIn(1) + yCoord * m_shearRateIn(2) + zCoord * m_shearRateIn(3)
    !    velocity(2) = m_uuuIn(2) + xCoord * m_shearRateIn(1) + 0 * m_shearRateIn(2) + zCoord * m_shearRateIn(3)
    !    velocity(3) = m_uuuIn(3) + xCoord * m_shearRateIn(1) + yCoord * m_shearRateIn(2) + 0 * m_shearRateIn(3)
    !END SUBROUTINE

    !SUBROUTINE evaluateOscillatoryVelocity(velocity,time)
    !    implicit none
    !    real(8):: velocity(1:SpaceDim)
    !    real(8):: time
    !    velocity(1) = m_uuuIn(1) + m_VelocityAmp * dcos(2*pi*m_VelocityFreq*time + m_VelocityPhi/180.0*pi)
    !    velocity(2) = m_uuuIn(2)
    !    velocity(3) = m_uuuIn(3)
    !END SUBROUTINE

    !subroutine  initDisturb_(this)
    !    implicit none
    !    class(LBMBlock), intent(inout) :: this
    !    integer:: x, y,z
    !    do  z = 1, this%zDim
    !    do  y = 1, this%yDim
    !    do  x = 1, this%xDim
    !         this%uuu(z,y,x,1)=this%uuu(z,y,x,1)+AmplInitDist(1)*m_Uref*dsin(2.0*pi*waveInitDist*xGrid(x))
    !         this%uuu(z,y,x,2)=this%uuu(z,y,x,2)+AmplInitDist(2)*m_Uref*dsin(2.0*pi*waveInitDist*xGrid(x))
    !         this%uuu(z,y,x,3)=this%uuu(z,y,x,3)+AmplInitDist(3)*m_Uref*dsin(2.0*pi*waveInitDist*xGrid(x))
    !    enddo
    !    enddo
    !    enddo
    !END SUBROUTINE

    SUBROUTINE setBndCond(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        ! set the boundary conditions of the block
        
    END SUBROUTINE setBndCond

    SUBROUTINE calculate_macro_quantities_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer::x,y,z
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
        do  x = 1, this%xDim
        do  y = 1, this%yDim
        do  z = 1, this%zDim
            this%den(z,y,x  )  = (SUM(this%fIn(z,y,x,0:lbmDim)))
            this%uuu(z,y,x,1)  = (SUM(this%fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,1))+0.5d0*this%volumeForce(1)*this%dh)/this%den(z,y,x)
            this%uuu(z,y,x,2)  = (SUM(this%fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,2))+0.5d0*this%volumeForce(2)*this%dh)/this%den(z,y,x)
            this%uuu(z,y,x,3)  = (SUM(this%fIn(z,y,x,0:lbmDim)*ee(0:lbmDim,3))+0.5d0*this%volumeForce(3)*this%dh)/this%den(z,y,x)
            !prs(z,y,x)   = Cs2*(den(z,y,x)-denIn)
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE calculate_macro_quantities_

    SUBROUTINE update_volumn_force_(this,time)
        implicit none
        class(LBMBlock), intent(inout) :: this
        real(8):: time
        this%volumeForce(1) = flow%volumeForceIn(1) + flow%volumeForceAmp * dsin(2.d0 * pi * flow%volumeForceFreq * time + flow%volumeForcePhi/180.0*pi)
        this%volumeForce(2) = flow%volumeForceIn(2)
        this%volumeForce(3) = flow%volumeForceIn(3)
    END SUBROUTINE

    SUBROUTINE addVolumForc(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: x
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x)
        do x=1,this%xDim
            this%force(:,:,x,1) = this%force(:,:,x,1) + this%volumeForce(1)
            this%force(:,:,x,2) = this%force(:,:,x,2) + this%volumeForce(2)
            this%force(:,:,x,3) = this%force(:,:,x,3) + this%volumeForce(3)
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE

    SUBROUTINE collision_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        real(8):: uSqr,uxyz(0:lbmDim),fEq(0:lbmDim),Flb(0:lbmDim),dt3
        integer:: x,y,z
        dt3 = 3.d0*this%dh
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,uSqr,uxyz,fEq,Flb)
        do    x = 1, this%xDim
        do    y = 1, this%yDim
        do    z = 1, this%zDim
            uSqr           = sum(this%uuu(z,y,x,1:3)**2)
            uxyz(0:lbmDim) = this%uuu(z,y,x,1) * ee(0:lbmDim,1) + this%uuu(z,y,x,2) * ee(0:lbmDim,2)+this%uuu(z,y,x,3) * ee(0:lbmDim,3)
            fEq(0:lbmDim)  = wt(0:lbmDim) * this%den(z,y,x) * ( (1.0d0 - 1.5d0 * uSqr) + uxyz(0:lbmDim) * (3.0d0  + 4.5d0 * uxyz(0:lbmDim)) ) - this%fIn(z,y,x,0:lbmDim)
            Flb(0:lbmDim)  = dt3*wt(0:lbmDim)*( &
                              (ee(0:lbmDim,1)-this%uuu(z,y,x,1)+3.d0*uxyz(0:lbmDim)*ee(0:lbmDim,1))*this%force(z,y,x,1) &
                             +(ee(0:lbmDim,2)-this%uuu(z,y,x,2)+3.d0*uxyz(0:lbmDim)*ee(0:lbmDim,2))*this%force(z,y,x,2) &
                             +(ee(0:lbmDim,3)-this%uuu(z,y,x,3)+3.d0*uxyz(0:lbmDim)*ee(0:lbmDim,3))*this%force(z,y,x,3))
            if(this%iCollidModel==1)then
                ! SRT collision
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + this%Omega*fEq(0:lbmDim) + (1.d0-0.5d0*this%Omega)*Flb(0:lbmDim)
            elseif(this%iCollidModel==2)then
                ! TRT collision
                fEq(0) = this%Omega*fEq(0) + (1.d0-0.5d0*this%Omega)*Flb(0)
                uxyz(positivedirs) = 0.5d0*this%Omega *(fEq(positivedirs)+fEq(negativedirs)) + (0.5d0-0.25d0* this%Omega)*(Flb(positivedirs)+Flb(negativedirs))
                uxyz(negativedirs) = 0.5d0*this%Omega2*(fEq(positivedirs)-fEq(negativedirs)) + (0.5d0-0.25d0*this%Omega2)*(Flb(positivedirs)-Flb(negativedirs))
                fEq(positivedirs) = uxyz(positivedirs) + uxyz(negativedirs)
                fEq(negativedirs) = uxyz(positivedirs) - uxyz(negativedirs)
                this%fIn(z,y,x,0:lbmDim) = this%fIn(z,y,x,0:lbmDim) + fEq
            elseif(this%iCollidModel==3)then
                ! MRT collision
                this%fIn(z,y,x,0:lbmDim)=this%fIn(z,y,x,0:lbmDim)+MATMUL( this%M_COLLID(0:lbmDim,0:lbmDim), fEq(0:lbmDim) ) + MATMUL( this%M_FORCE(0:lbmDim,0:lbmDim),Flb(0:lbmDim))
            endif
        enddo
        enddo
        enddo
        !$OMP END PARALLEL DO
    END SUBROUTINE collision_

    !0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    !0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0
    !0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1
    !0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1
    SUBROUTINE streaming_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: i
        do  i=0,lbmDim
            call swapzy(this%fIn, ee(i,3), ee(i,2), i, this%zDim, this%yDim, this%xDim, lbmDim)
            call swapx(this%fIn, ee(i,1), i, this%zDim, this%yDim, this%xDim, lbmDim, this%OMPparindex, this%OMPeid, this%OMPedge)
        enddo

        contains

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

        SUBROUTINE swapx(f, dx, i, zDim, yDim, xDim, lbmDim, parindex, eid, edge)
            implicit none
            integer, intent(in):: dx, i, zDim, yDim, xDim, lbmDim, parindex(:)
            integer, intent(out):: eid(:)
            real(8), intent(out):: edge(:,:,:)
            real(8), intent(inout):: f(1:zDim, 1:yDim,1:xDim,0:lbmDim)
            integer:: p

            if(dx.eq.0) return

            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
            do  p = 1,this%npsize
                if(dx .eq. -1) then
                    call swapxwAtom(f, edge(:,:,p), eid(p), i, zDim, yDim, xDim, lbmDim, parindex(p), parindex(p+1)-1)
                elseif(dx.eq.1) then
                    call swapxeAtom(f, edge(:,:,p), eid(p), i, zDim, yDim, xDim, lbmDim, parindex(p), parindex(p+1)-1)
                endif
            enddo
            !$OMP END PARALLEL DO
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(p)
            do  p = 1,this%npsize
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
        end subroutine
    END SUBROUTINE streaming_


    SUBROUTINE write_flow_(this,time)
        implicit none
        class(LBMBlock), intent(inout) :: this
        real(8), intent(in):: time
        integer:: x,y,z,pid,i
        integer::xmin,xmax,ymin,ymax,zmin,zmax
        integer,parameter::nameLen=10,idfile=100
        character (LEN=nameLen):: fileName
        real(8):: invUref
        real(8):: get_cpu_time, waittime
        integer,dimension(8) :: values0,values1
        call date_and_time(VALUES=values0)
        call mywait()
        call date_and_time(VALUES=values1)
        waittime = get_cpu_time(values1)-get_cpu_time(values0)
        if(waittime.gt.1.d-1) then
            write(*,'(A,F7.2,A)')'Waiting ', waittime, 's for previous outflow finishing.'
        endif
        xmin = 1 + this%offsetOutput
        ymin = 1 + this%offsetOutput
        zmin = 1 + this%offsetOutput
        xmax = this%xDim - this%offsetOutput
        ymax = this%yDim - this%offsetOutput
        zmax = this%zDim - this%offsetOutput
        invUref = 1.d0/flow%Uref
        
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
        do x=xmin, xmax
            do y=ymin, ymax
                do z=zmin, zmax
                    this%oututmp(z,y,x) = this%uuu(z,y,x,1)*invUref
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
        do x=xmin, xmax
            do y=ymin, ymax
                do z=zmin, zmax
                    this%outvtmp(z,y,x) = this%uuu(z,y,x,2)*invUref
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO
        !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z)
        do x=xmin, xmax
            do y=ymin, ymax
                do z=zmin, zmax
                    this%outwtmp(z,y,x) = this%uuu(z,y,x,3)*invUref
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO
        !this%offsetMoveGrid=0.0
        !if(isMoveGrid==1)then
        !    if(isMoveDimX==1) this%offsetMoveGrid(1) = dh*dble(MoveOutputIref(1))
        !    if(isMoveDimY==1) this%offsetMoveGrid(2) = dh*dble(MoveOutputIref(2))
        !    if(isMoveDimZ==1) this%offsetMoveGrid(3) = dh*dble(MoveOutputIref(3))
        !endif
        call myfork(pid)
        if(pid.eq.0) then
            write(fileName,'(I10)') nint(time/flow%Tref*1d5)
            fileName = adjustr(fileName)
            do  i=1,nameLen
                if(fileName(i:i)==' ')fileName(i:i)='0'
            enddo
            open(idfile,file='./DatFlow/Flow'//trim(fileName)//'.plt',form='unformatted',access='stream')
            WRITE(idfile) xmin,xmax,ymin,ymax,zmin,zmax
            WRITE(idfile) this%offsetMoveGrid(1:3)
            write(idfile) this%oututmp,this%outvtmp,this%outwtmp
            close(idfile)
            call myexit(0)
        endif
    END SUBROUTINE write_flow_

    subroutine ComputeFieldStat_(this)
        implicit none
        class(LBMBlock), intent(inout) :: this
        integer:: x,y,z,i
        real(8):: invUref, uLinfty(1:3), uL2(1:3), temp
        invUref = 1.d0/flow%Uref
        uLinfty = -1.d0
        uL2 = 0.d0
        do i=1,3
            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(x,y,z,temp) &
            !$OMP reduction (+: uL2) reduction(max: uLinfty)
            do x=1, this%xDim
                do y=1, this%yDim
                    do z=1, this%zDim
                        temp = dabs(this%uuu(z,y,x,i)*invUref)
                        uL2(i) = uL2(i) + temp * temp
                        if(temp.gt.uLinfty(i)) uLinfty(i) = temp
                    enddo
                enddo
            enddo
            !$OMP END PARALLEL DO
            uL2(i) = dsqrt(uL2(i) / (dble(this%xDim) * dble(this%yDim) * dble(this%zDim)))
        enddo
        write(*,'(A,F18.12)')'FIELDSTAT L2 u ', uL2(1)
        write(*,'(A,F18.12)')'FIELDSTAT L2 v ', uL2(2)
        write(*,'(A,F18.12)')'FIELDSTAT L2 w ', uL2(3)
        write(*,'(A,F18.12)')'FIELDSTAT Linfinity u ', uLinfty(1)
        write(*,'(A,F18.12)')'FIELDSTAT Linfinity v ', uLinfty(2)
        write(*,'(A,F18.12)')'FIELDSTAT Linfinity w ', uLinfty(3)
    endsubroutine ComputeFieldStat_

    SUBROUTINE write_continue_(this,filename,step,time)
        IMPLICIT NONE
        class(LBMBlock), intent(inout) :: this
        character(LEN=40),intent(in):: filename
        integer,intent(in):: step
        real(8),intent(in):: time
        open(unit=13,file=filename,form='unformatted',status='replace')
        write(13) step,time
        write(13) this%fIn
        close(13)
    ENDSUBROUTINE write_continue_

    SUBROUTINE read_continue_(this,filename,step,time)
        IMPLICIT NONE
        class(LBMBlock), intent(inout) :: this
        character(LEN=40),intent(in):: filename
        integer,intent(out):: step
        real(8),intent(out):: time
        open(unit=13,file=filename,form='unformatted',status='old')
        read(13) step,time
        read(13) this%fIn
        close(13)
    ENDSUBROUTINE read_continue_

end module FluidDomain