module SythTurb
    use ConstParams
    use FlowCondition
    implicit none
    private
    public:: InitialiseTurbulentBC, GetDisturbeVelocity
    integer,parameter:: Dim=3,Nk=200
    real(8),parameter:: RL=40000
    ! real(8),parameter:: L=1
    ! real(8),parameter:: nu=1e-6
    real(8):: eta,epsilon,u_eta,k0,kc,K(1:Nk+1),E_k(1:Nk+1)
    real(8):: ukn(1:Nk),kn(1:Nk)
    real(8):: tke,tn,tm
    real(8):: m_nu,m_L
    real(8):: turbIntensity
    real(8),allocatable:: buffer(:,:), rands(:,:)
    integer:: windex, rindex, bsize

    contains

    subroutine InitialiseTurbulentBC(turbIntensity, buffersize)
        implicit none
        if(.not.allocated(buffer)) then
            bsize = buffersize
            allocate(buffer(1:Dim, 1:bsize), rands(0:3,1:bsize))
            rindex = 1
            do windex = 1, bsize
                call UpdateVelocity(windex)
            enddo
            call myfork(pid)
            if(pid.eq.0) then
                do while (.true.)
                    if(windex>bsize) then
                        windex = 1
                        call Rand_Generation()
                    endif
                    call UpdateDisturbeVelocity(buffer(windex))
                    windex = windex + 1
                enddo
            endif
        endif
    end subroutine Initialise

    subroutine GetDisturbeVelocity(vel)
        implicit none
        real(8):: vel(1:Dim), r
        rindex =  round(rands(3,rindex) * bsize)
        if(rindex<1) rindex = 1
        if(rindex>bsize) rindex = bsize
        vel = buffer(:,rindex)
    end subroutine GetTurbVelocity

    subroutine initparams()
        implicit none
        integer:: xDim,yDim,zDim
        integer:: i,j,m
        real(8):: d
        real(8),allocatable:: KN_(:,:,:),Psi(:,:),Sigma(:,:,:),point_uHIT(:,:)
        NX      = xDim
        NY      = yDim
        NZ      = zDim
        Nt      = 1
        m_L     = flow%Lref
        m_nu    = flow%nu
        eta     = m_L/(RL)**(3d0/4d0)
        epsilon = (m_nu**3d0)/(eta**4d0)

        u_eta   = (epsilon*m_nu)**(1d0/4d0)
        k0      = 0.2d0/m_L
        kc      = 2d0/eta

        do i = 1,Nk+1
            d    = (log(kc)-log(k0))*dble(i-1)/dble(Nk)+log(k0)
            K(i) = exp(d)
        enddo
        call Ek()

        allocate(x(NX),y(NY),z(NZ),t(Nt))
        do i = 1,NX
            x(i) = (m_L-0)*dble(i-1)/dble(NX-1)+0d0
        enddo
        do j = 1,NY
            y(j) = (m_L-0)*dble(j-1)/dble(NY-1)+0d0
        enddo
        do m = 1,NZ
            z(m) = (m_L-0)*dble(m-1)/dble(NZ-1)+0d0
        enddo
        if (Nt.eq.1)then
            t(1) = tm
        else
            do i = 1,Nt
                t(i) = (tm-0)*dble(i-1)/dble(Nt-1)+0d0
            enddo
        endif

        ukn = sqrt((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
        kn  = 0.5d0*(K(1:Nk)+K(2:Nk+1))

        tke = sum((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
        tn  = (m_nu/epsilon)**0.5d0
        tm  = 100d0*tn

        allocate(KN_(Nt,Nk,Dim),Psi(Nt,Nk),Sigma(Nt,Nk,Dim),point_uHIT(Nt,Dim))

        if (Dim.eq.3) call Flow_Generation(KN_,Sigma,Psi)

        do i = 1,NX
            do j = 1,NY
                do m = 1,NZ
                    call point_uHIT_Generate(KN_,Sigma,Psi,i,j,m,point_uHIT)
                    uHIT(:,i,j,m,:) = point_uHIT(:,:)
                enddo
            enddo
        enddo
    endsubroutine

    subroutine Ek()
        implicit none
        ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
        ! von Karman spectrum
        real(8):: c,beta,p0,cl,cn
        real(8):: fl(1:Nk+1),fn(1:Nk+1)
        integer:: i

        c    = 1.5d0
        beta = 5.2d0
        p0   = 2.0d0
        cl   = 6.78d0
        cn   = 0.4d0

        do i = 1,Nk+1
            fl(i)  = (K(i)*m_L/sqrt((K(i)*m_L)**2d0+cl))**(5d0/3d0+p0)
            fn(i)  = exp(-beta*(((K(i)*eta)**4d0+cn**4d0)**0.25d0-cn))
            E_k(i) = c*epsilon**(2d0/3d0)*K(i)**(-5d0/3d0)*fl(i)*fn(i)
        enddo
    endsubroutine

    subroutine Rand_Generation(KN_,Sigma,Psi)
        implicit none
        ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
        ! von Karman spectrum
        real(8):: r(Nt,Nk),Theta(Nt,Nk),Phi(Nt,Nk),Alpha(Nt,Nk),k_dir(Nt,Nk,Dim)
        real(8),intent(out):: Psi(Nt,Nk),Sigma(Nt,Nk,Dim)
        real(8):: theta_,psi_,phi_,alpha_,Ry(3,3),Rz(3,3),sigma0(3),sigma_(3)
        real(8),intent(out):: KN_(Nt,Nk,Dim)
        integer:: i,j,m

        ! Dim=3
        call random_number(r)
        Theta(:,:) = asin(r(:,:))
        call random_number(r)
        Psi(:,:) = 2d0*pi*r(:,:)
        call random_number(r)
        Phi(:,:) = 2d0*pi*r(:,:)
        call random_number(r)
        Alpha(:,:) = 2d0*pi*r(:,:)
        Sigma=0d0;
        do i=1,Nt
            do j=1,Nk
                theta_=Theta(i,j);
                psi_=Psi(i,j);
                phi_=Phi(i,j);
                alpha_=Alpha(i,j);

                sigma0=[cos(alpha_),sin(alpha_),0.0d0]
                Ry(1,1) = cos(theta_)
                Ry(1,2) = 0.0d0
                Ry(1,3) = sin(theta_)
                Ry(2,1) = 0.0d0
                Ry(2,2) = 1.0d0
                Ry(2,3) = 0.0d0
                Ry(2,1) =-sin(theta_)
                Ry(3,2) = 0.0d0
                Ry(3,3) = cos(theta_)

                Rz(1,1) = cos(phi_)
                Rz(1,2) =-sin(phi_)
                Rz(1,3) = 0.0d0
                Rz(2,1) = sin(phi_)
                Rz(2,2) = cos(phi_)
                Rz(2,3) = 0.0d0
                Rz(2,1) = 0.0d0
                Rz(3,2) = 0.0d0
                Rz(3,3) = 1.0d0

                sigma_=matmul(matmul(Rz,Ry),sigma0)
                Sigma(i,j,:)=sigma_(:);
            enddo
        enddo
        k_dir=0d0;
        k_dir(:,:,1)=cos(Phi(:,:))*sin(Theta(:,:));
        k_dir(:,:,2)=sin(Phi(:,:))*sin(Theta(:,:));
        k_dir(:,:,3)=cos(Theta(:,:));

        call repmat1dim3(kn,Nk,Nt,1,Dim,KN_)
        KN_ = KN_(:,:,:)*k_dir(:,:,:)

        Sigma(:,:,1)=-sin(Phi(:,:));
        Sigma(:,:,2)= cos(Phi(:,:));

        ! open(111,file='1.dat',position='append')
        ! write(111,'(150E20.10)') uHIT(1,:,:,1,1)
        ! close(111)
    endsubroutine

    subroutine UpdateDisturbeVelocity(uHit)
        !, KN_,Sigma,Psi,i,j,m,point_uHIT)
        implicit none
        !write at buffer(:,windex)
        real(8),intent(in):: KN_(Nt,Nk,Dim),Sigma(Nt,Nk,Dim),Psi(Nt,Nk)
        real(8):: tmp(Nt,Nk),ukn_(Nt,Nk),tmp_(Nt,Nk,Dim),uHIT_(Nt,Nk,Dim)
        real(8),intent(out):: point_uHIT(Nt,Dim)
        integer:: i,j,m

        tmp = cos(KN_(:,:,1)*x(i)+KN_(:,:,2)*y(j)+KN_(:,:,3)*z(m)+Psi(:,:))

        call repmat1dim2(ukn,Nk,Nt,1,ukn_)

        call repmat2dim3(tmp(:,:)*ukn_(:,:),Nt,Nk,1,1,Dim,tmp_)

        uHIT_(:,:,:) = tmp_(:,:,:)*Sigma(:,:,:)
        point_uHIT(:,:) = sum(uHIT_,dim=2)
    endsubroutine

    subroutine repmat1dim3(mat,matDim2,Dim1,Dim2,Dim3,mat_)
        implicit none
        integer:: matDim2,Dim1,Dim2,Dim3
        real(8):: mat(matDim2),mat_(Dim1,matDim2*Dim2,Dim3)
        integer:: a,b,c,p
        do a = 1,Dim1
            do b = 1,Dim2
                do c = 1,Dim3
                    do p = 1,matDim2
                        mat_(a,(b-1)*Dim2+p,c) = mat(p)
                    enddo
                enddo
            enddo
        enddo
    endsubroutine

    subroutine repmat1dim2(mat,matDim2,Dim1,Dim2,mat_)
        implicit none
        integer:: matDim2,Dim1,Dim2
        real(8):: mat(matDim2),mat_(Dim1,matDim2*Dim2)
        integer:: a,b,p
        do a = 1,Dim1
            do b = 1,Dim2
                do p = 1,matDim2
                    mat_(a,(b-1)*Dim2+p) = mat(p)
                enddo
            enddo
        enddo
    endsubroutine

    subroutine repmat2dim3(mat,matDim1,matDim2,Dim1,Dim2,Dim3,mat_)
        implicit none
        integer:: matDim1,matDim2,Dim1,Dim2,Dim3
        real(8):: mat(matDim1,matDim2),mat_(matDim1*Dim1,matDim2*Dim2,Dim3)
        integer:: a,b,c,p,q
        do a = 1,Dim1
            do b = 1,Dim2
                do c = 1,Dim3
                    do p = 1,matDim1
                        do q = 1,matDim2
                            mat_((a-1*Dim1)+p,(b-1)*Dim2+q,c) = mat(p,q)
                        enddo
                    enddo
                enddo
            enddo
        enddo
    endsubroutine

end module SythTurb



