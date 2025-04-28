program input
    implicit none
    integer,parameter:: Dim=2,Nk=200
    real(8),parameter:: RL=40000
    real(8),parameter:: L=1
    real(8):: nu,n,eps,un,k0,kc,K(1:Nk+1),E_k(1:Nk+1)
    integer:: i
    real(8):: d

    nu  = 1d-6
    n   = L/(RL)**(3d0/4d0)
    eps = (nu**3d0)/(n**4d0)

    un  = (eps*nu)**(1d0/4d0)
    k0  = 0.2d0/L
    kc  = 2d0/n

    do i = 1,Nk+1
        d    = (log(kc)-log(k0))*dble(i-1)/dble(Nk)+log(k0)
        K(i) = exp(d)
    enddo
    call Ek(n,eps,Nk,K,L,E_k)
    call Flow_Generation(Dim,nu,eps,Nk,K,L,E_k)

end program input

subroutine Ek(n,eps,Nk,K,L,E_k)
    implicit none
    ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
    ! von Karman spectrum
    integer:: Nk
    real(8):: L
    real(8):: n,eps,K(1:Nk+1),E_k(1:Nk+1)
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
        fn(i)  = exp(-beta*(((K(i)*n)**4d0+cn**4d0)**0.25d0-cn))
        E_k(i) = c*eps**(2d0/3d0)*K(i)**(-5d0/3d0)*fl(i)*fn(i)
    enddo
endsubroutine


subroutine Flow_Generation(Dim,nu,eps,Nk,K,L,E_k)
    implicit none
    ! Turbulence power spectrum E(K) with -5/3 Kolmogorov spectrum 
    ! von Karman spectrum
    real(8), parameter:: pi = 3.141592653589793d0
    integer:: Dim,Nk
    real(8):: L
    real(8):: nu,eps,K(1:Nk+1),E_k(1:Nk+1)
    real(8):: ukn(1:Nk),kn(1:Nk)
    real(8):: tke,tn,tm
    integer:: N,NX,NY,NZ,Nt
    real(8),allocatable:: x(:),y(:),z(:),t(:)
    real(8),allocatable:: Psi(:,:),Phi(:,:),Sigma(:,:,:),k_dir(:,:,:)
    real(8),allocatable:: KN_(:,:,:),uHIT(:,:,:,:),tmp(:,:),ukn_(:,:),tmp_(:,:,:),uHIT_(:,:,:)
    integer:: i,j
    real(8),allocatable:: r(:,:)


    ukn = sqrt((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
    kn  = 0.5d0*(K(1:Nk)+K(2:Nk+1))

    tke = sum((E_k(1:Nk)+E_k(2:Nk+1))/2d0*(K(2:Nk+1)-K(1:Nk)))
    tn  = (nu/eps)**0.5d0
    tm  = 100d0*tn
    
    N  = 150
    NX = N
    NY = N
    NZ = N
    Nt = 1
    allocate(x(NX),y(NY),z(NZ),t(Nt))
    do i = 1,N
        x(i) = (L-0)*dble(i-1)/dble(N-1)+0d0
        y(i) = (L-0)*dble(i-1)/dble(N-1)+0d0
        z(i) = (L-0)*dble(i-1)/dble(N-1)+0d0
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
    call random_seed()
    call random_number(r)
    Psi(:,:) = 2d0*pi*r(:,:)
    call random_seed()
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

    open(111,file='1.dat',position='append')
    write(111,'(150E20.10)') uHIT(1,:,:,1)
    close(111)

endsubroutine

subroutine repmat1dim3(mat,matDim2,Dim1,Dim2,Dim3,mat_)
    implicit none
    integer:: matDim2,Dim1,Dim2,Dim3
    real(8):: mat(matDim2),mat_(Dim1,matDim2*Dim2,Dim3)
    integer:: i,j,k,m
    do i = 1,Dim1
        do j = 1,Dim2
            do k = 1,Dim3
                do m = 1,matDim2
                    mat_(i,(j-1)*Dim2+m,k) = mat(m)
                enddo
            enddo
        enddo
    enddo
endsubroutine

subroutine repmat1dim2(mat,matDim2,Dim1,Dim2,mat_)
    implicit none
    integer:: matDim2,Dim1,Dim2
    real(8):: mat(matDim2),mat_(Dim1,matDim2*Dim2)
    integer:: i,j,m
    do i = 1,Dim1
        do j = 1,Dim2
            do m = 1,matDim2
                mat_(i,(j-1)*Dim2+m) = mat(m)
            enddo
        enddo
    enddo
endsubroutine

subroutine repmat2dim3(mat,matDim1,matDim2,Dim1,Dim2,Dim3,mat_)
    implicit none
    integer:: matDim1,matDim2,Dim1,Dim2,Dim3
    real(8):: mat(matDim1,matDim2),mat_(matDim1*Dim1,matDim2*Dim2,Dim3)
    integer:: i,j,k,m,n
    do i = 1,Dim1
        do j = 1,Dim2
            do k = 1,Dim3
                do m = 1,matDim1
                    do n = 1,matDim2
                        mat_((i-1*Dim1)+m,(j-1)*Dim2+n,k) = mat(m,n)
                    enddo
                enddo
            enddo
        enddo
    enddo
endsubroutine



