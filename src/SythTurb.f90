module SythTurb
    use FlowCondition
    implicit none
    private
    public:: SythTurb_Input
    integer,parameter:: Dim=3,Nk=200
    real(8),parameter:: RL=40000
    real(8),parameter:: L=1
    real(8):: nu,eta,epsilon,u_eta,k0,kc,K(1:Nk+1),E_k(1:Nk+1)

    contains

    subroutine SythTurb_Input()
        implicit none
        integer:: i
        real(8):: d
        nu      = 1d-6
        eta     = L/(RL)**(3d0/4d0)
        epsilon = (nu**3d0)/(eta**4d0)

        u_eta   = (epsilon*nu)**(1d0/4d0)
        k0      = 0.2d0/L
        kc      = 2d0/eta

        do i = 1,Nk+1
            d    = (log(kc)-log(k0))*dble(i-1)/dble(Nk)+log(k0)
            K(i) = exp(d)
        enddo
        call Ek()
        if (Dim.eq.2) call Flow_Generation_dim2()
        if (Dim.eq.3) call Flow_Generation_dim3()
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
        p0   = 2d0
        cl   = 6.78d0
        cn   = 0.4d0

        do i = 1,Nk+1
            fl(i)  = (K(i)*L/sqrt((K(i)*L)**2d0+cl))**(5d0/3d0+p0)
            fn(i)  = exp(-beta*(((K(i)*eta)**4d0+cn**4d0)**0.25d0-cn))
            E_k(i) = c*epsilon**(2d0/3d0)*K(i)**(-5d0/3d0)*fl(i)*fn(i)
        enddo
    endsubroutine

    subroutine Flow_Generation_dim2()
        implicit none
        ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
        ! von Karman spectrum
        real(8), parameter:: pi = 3.141592653589793d0
        real(8):: ukn(1:Nk),kn(1:Nk)
        real(8):: tke,tn,tm
        integer:: NX,NY,NZ,Nt
        real(8),allocatable:: x(:),y(:),z(:),t(:)
        real(8),allocatable:: Psi(:,:),Phi(:,:),Sigma(:,:,:),k_dir(:,:,:)
        real(8),allocatable:: KN_(:,:,:),uHIT(:,:,:,:),tmp(:,:),ukn_(:,:),tmp_(:,:,:),uHIT_(:,:,:)
        integer:: i,j,m
        real(8),allocatable:: r(:,:)


        ukn = sqrt((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
        kn  = 0.5d0*(K(1:Nk)+K(2:Nk+1))

        tke = sum((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
        tn  = (nu/epsilon)**0.5d0
        tm  = 100d0*tn
        
        NX = 150
        NY = 150
        NZ = 150
        Nt = 1
        allocate(x(NX),y(NY),z(NZ),t(Nt))
        do i = 1,NX
            x(i) = (L-0)*dble(i-1)/dble(NX-1)+0d0
        enddo
        do j = 1,NY
            y(j) = (L-0)*dble(j-1)/dble(NY-1)+0d0
        enddo
        do m = 1,NZ
            z(m) = (L-0)*dble(m-1)/dble(NZ-1)+0d0
        enddo
        if (Nt.eq.1)then
            t(1) = tm
        else
            do i = 1,Nt
                t(i) = (tm-0)*dble(i-1)/dble(Nt-1)+0d0
            enddo
        endif

        ! Dim=2
        allocate(r(Nt,Nk),Psi(Nt,Nk),Phi(Nt,Nk),Sigma(Nt,Nk,Dim),k_dir(Nt,Nk,Dim))
        call random_number(r)
        Psi(:,:) = 2d0*pi*r(:,:)
        call random_number(r)
        Phi(:,:) = 2d0*pi*r(:,:)
        Sigma=0d0;
        k_dir=0d0;
        k_dir(:,:,1)=cos(Phi(:,:));
        k_dir(:,:,2)=sin(Phi(:,:));

        allocate(KN_(Nt,Nk,Dim),uHIT(Nt,NX,NY,Dim),tmp(Nt,Nk),ukn_(Nt,Nk),tmp_(Nt,Nk,Dim),uHIT_(Nt,Nk,Dim))
        uHIT = 0d0
        call repmat1dim3(kn,Nk,Nt,1,Dim,KN_)
        KN_ = KN_(:,:,:)*k_dir(:,:,:)
    
        Sigma(:,:,1)=-sin(Phi(:,:));
        Sigma(:,:,2)= cos(Phi(:,:));

        do i=1,NX
            do j=1,NY
                tmp = cos(KN_(:,:,1)*x(i)+KN_(:,:,2)*y(j)+Psi(:,:))

                call repmat1dim2(ukn,Nk,Nt,1,ukn_)

                call repmat2dim3(tmp(:,:)*ukn_(:,:),Nt,Nk,1,1,Dim,tmp_)

                uHIT_(:,:,:) = tmp_(:,:,:)*Sigma(:,:,:)
                uHIT(:,i,j,:) = sum(uHIT_,dim=2)
            enddo
        enddo

        ! open(111,file='1.dat',position='append')
        ! write(111,'(150E20.10)') uHIT(1,:,:,1)
        ! close(111)

    endsubroutine

    subroutine Flow_Generation_dim3()
        implicit none
        ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
        ! von Karman spectrum
        real(8), parameter:: pi = 3.141592653589793d0
        real(8):: ukn(1:Nk),kn(1:Nk)
        real(8):: tke,tn,tm
        integer:: NX,NY,NZ,Nt
        real(8),allocatable:: x(:),y(:),z(:),t(:)
        real(8),allocatable:: Theta(:,:),Psi(:,:),Phi(:,:),Alpha(:,:),Sigma(:,:,:),k_dir(:,:,:)
        real(8):: theta_,psi_,phi_,alpha_,Ry(3,3),Rz(3,3),sigma0(3),sigma_(3)
        real(8),allocatable:: KN_(:,:,:),uHIT(:,:,:,:,:),tmp(:,:),ukn_(:,:),tmp_(:,:,:),uHIT_(:,:,:)
        integer:: i,j,m
        real(8),allocatable:: r(:,:)


        ukn = sqrt((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
        kn  = 0.5d0*(K(1:Nk)+K(2:Nk+1))

        tke = sum((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
        tn  = (nu/epsilon)**0.5d0
        tm  = 100d0*tn
        
        NX = 150
        NY = 150
        NZ = 3
        Nt = 1
        allocate(x(NX),y(NY),z(NZ),t(Nt))
        do i = 1,NX
            x(i) = (L-0)*dble(i-1)/dble(NX-1)+0d0
        enddo
        do j = 1,NY
            y(j) = (L-0)*dble(j-1)/dble(NY-1)+0d0
        enddo
        do m = 1,NZ
            z(m) = (L-0)*dble(m-1)/dble(NZ-1)+0d0
        enddo
        if (Nt.eq.1)then
            t(1) = tm
        else
            do i = 1,Nt
                t(i) = (tm-0)*dble(i-1)/dble(Nt-1)+0d0
            enddo
        endif

        ! Dim=3
        allocate(r(Nt,Nk),Theta(Nt,Nk),Psi(Nt,Nk),Phi(Nt,Nk),Alpha(Nt,Nk),Sigma(Nt,Nk,Dim),k_dir(Nt,Nk,Dim))
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

        allocate(KN_(Nt,Nk,Dim),uHIT(Nt,NX,NY,NZ,Dim),tmp(Nt,Nk),ukn_(Nt,Nk),tmp_(Nt,Nk,Dim),uHIT_(Nt,Nk,Dim))
        uHIT = 0d0
        call repmat1dim3(kn,Nk,Nt,1,Dim,KN_)
        KN_ = KN_(:,:,:)*k_dir(:,:,:)

        Sigma(:,:,1)=-sin(Phi(:,:));
        Sigma(:,:,2)= cos(Phi(:,:));

        do i=1,NX
            do j=1,NY
                do m=1,NZ
                    tmp = cos(KN_(:,:,1)*x(i)+KN_(:,:,2)*y(j)+KN_(:,:,3)*z(m)+Psi(:,:))

                    call repmat1dim2(ukn,Nk,Nt,1,ukn_)

                    call repmat2dim3(tmp(:,:)*ukn_(:,:),Nt,Nk,1,1,Dim,tmp_)

                    uHIT_(:,:,:) = tmp_(:,:,:)*Sigma(:,:,:)
                    uHIT(:,i,j,m,:) = sum(uHIT_,dim=2)
                enddo
            enddo
        enddo

        ! open(111,file='1.dat',position='append')
        ! write(111,'(150E20.10)') uHIT(1,:,:,1,1)
        ! close(111)
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



