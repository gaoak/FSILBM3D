module SythTurb
    use ConstParams
    implicit none
    private
    public:: Initparams, InitialiseTurbulentBC, SyntheticVelocity
    integer,parameter:: Dim=3,Nk=200
    real(8):: RL
    real(8):: m_nu,m_L,m_U,m_xmin,m_ymin,m_zmin,m_dh
    real(8):: ukn(1:Nk),kn(1:Nk),tke,tn
    integer:: NX,NY,NZ,Nt
    real(8),allocatable:: t(:)
    real(8),allocatable:: vel_buffer(:,:), rands(:,:,:),rand_psi(:)
    integer:: windex, rindex, bsize

    contains

    subroutine InitialiseTurbulentBC(buffersize)
        implicit none
        integer:: pid,buffersize
        if(.not.allocated(vel_buffer)) then
            bsize = buffersize
            allocate(vel_buffer(1:Dim, 1:bsize), rands(0:1,1:Dim,1:Nk),rand_psi(1:Nk))
            rindex = 1
            windex = 1
            call myfork(pid)
            if(pid.eq.0) then
                do while (.true.)
                    if(windex>bsize) then
                        windex = 1
                        call Rand_Generation()
                    endif
                    call UpdateDisturbeVelocity(vel_buffer(:,windex))
                    windex = windex + 1
                enddo
            endif
        endif
    end subroutine InitialiseTurbulentBC

    subroutine SyntheticVelocity(vel)
        implicit none
        real(8):: vel(1:Dim)
        rindex = nint(rand_psi(rindex) * bsize)
        if(rindex<1) rindex = 1
        if(rindex>bsize) rindex = bsize
        vel(:) = vel(:)+vel_buffer(:,rindex)
    end subroutine SyntheticVelocity

    subroutine Initparams(turbIntensity,Uref,Lref,nu,xDim,yDim,zDim,xmin,ymin,zmin,dh)
        implicit none
        real(8):: turbIntensity,Uref,Lref,nu,xmin,ymin,zmin,dh
        integer:: xDim,yDim,zDim
        real(8):: eta,epsilon,u_eta,k0,kc,K(1:Nk+1),E_k(1:Nk+1)
        integer:: i
        real(8):: d
        NX      = xDim
        NY      = yDim
        NZ      = zDim
        Nt      = 1
        m_U     = Uref
        m_L     = Lref
        m_nu    = nu
        m_xmin  = xmin
        m_ymin  = ymin
        m_zmin  = zmin
        m_dh  = dh
        RL      = turbIntensity*m_U*m_L/m_nu
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

        allocate(t(Nt))
        tn  = (m_nu/epsilon)**0.5d0
        t(1) = tn

        ukn = sqrt((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
        kn  = 0.5d0*(K(1:Nk)+K(2:Nk+1))

        tke = sum((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))

        contains

        subroutine Ek()
            implicit none
            ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
            ! von Karman spectrum
            real(8):: c,beta,p0,cl,cn
            real(8):: fl(1:Nk+1),fn(1:Nk+1)

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
    endsubroutine

    subroutine Rand_Generation()
        implicit none
        ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
        ! von Karman spectrum
        real(8):: r(Nt,Nk),Theta(Nt,Nk),Phi(Nt,Nk),Alpha(Nt,Nk),k_dir(Nt,Nk,Dim)
        real(8):: KN_(Nt,Nk,Dim),Psi(Nt,Nk),Sigma(Nt,Nk,Dim)
        real(8):: theta_,psi_,phi_,alpha_,Ry(3,3),Rz(3,3),sigma0(3),sigma_(3)
        integer:: i,j

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

        do i =1,Nk
            rands(0,:,i) = KN_(1,i,:)
            rands(1,:,i) = Sigma(1,i,:)
            rand_psi(i)  = Psi(1,i)
        enddo
    endsubroutine

    subroutine UpdateDisturbeVelocity(uHit)
        implicit none
        real(8):: Sigma(Nt,Nk,Dim)
        real(8):: tmp(Nt,Nk),ukn_(Nt,Nk),tmp_(Nt,Nk,Dim),uHIT_(Nk,Dim)
        real(8),intent(out):: uHit(Dim)
        real(8):: x,y,z
        integer:: i,j,m,dd

        if (NZ.eq.1) then
            i = (windex-1)/NX+1
            j = mod(windex-1,NX)+1
            m = 1
        elseif (NY.eq.1) then
            i = mod(windex-1,NZ)+1
            j = 1
            m = (windex-1)/NZ+1
        elseif (NX.eq.1) then
            i = 1
            j = (windex-1)/NY+1
            m = mod(windex-1,NY)+1
        endif

        x = m_xmin + (i - 1) * m_dh
        y = m_ymin + (j - 1) * m_dh
        z = m_zmin + (m - 1) * m_dh

        tmp(1,:) = cos(rands(0,1,:)*x+rands(0,2,:)*y+rands(0,3,:)*z+rand_psi(:))

        call repmat1dim2(ukn,Nk,Nt,1,ukn_)

        call repmat2dim3(tmp(:,:)*ukn_(:,:),Nt,Nk,1,1,Dim,tmp_)

        do dd = 1,Dim
            Sigma(1,:,dd) = rands(0,dd,:)
        enddo
        uHIT_(:,:) = tmp_(1,:,:)*Sigma(1,:,:)
        uHit(:) = sum(uHIT_(:,:),dim=1)
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



