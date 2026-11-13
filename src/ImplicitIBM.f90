! SPDX-License-Identifier: GPL-3.0-or-later
module ImplicitIBM
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ParallelDomain, only: is_root, parallel_abort, sum_matrix, sum_values
    implicit none
    private
    public :: solve_implicit_ibm

contains

    subroutine solve_implicit_ibm(first_x,den,velocity,indices,weights,target,max_iterations,tolerance,uref,correction)
        integer, intent(in) :: first_x,max_iterations
        real(8), intent(in) :: den(:,:,:),target(:,:),tolerance,uref
        real(8), intent(inout) :: velocity(:,:,:,:)
        integer(2), intent(in) :: indices(:,:)
        real(4), intent(in) :: weights(:,:)
        real(8), intent(out) :: correction(:,:)

        integer :: n,nodes,marker,node,x,y,z,ix(4),iy(4),iz(4),i,j,k,s,iter,c
        integer, allocatable :: grid_slot(:,:,:),xyz(:,:),slot(:,:)
        real(8) :: wx(4),wy(4),wz(4),weight,threshold,error,alpha,beta,curvature
        real(8) :: rz(3),next_rz(3)
        logical :: active(3)
        real(8), allocatable :: stencil(:,:),cell(:,:),inverse_den(:),diagonal(:)
        real(8), allocatable :: rhs(:,:),residual(:,:),direction(:,:),product(:,:),preconditioned(:,:)

        if (max_iterations < 1 .or. .not. ieee_is_finite(tolerance) .or. tolerance <= 0.d0 .or. &
            .not. ieee_is_finite(uref) .or. uref <= 0.d0) then
            call parallel_abort('ImplicitDirect IBM requires positive ntolLBM, dtolLBM and Uref')
        endif
        threshold = tolerance * uref
        n = size(target,2)
        allocate(grid_slot(size(den,1),size(den,2),size(den,3)), xyz(3,min(64*n,size(den))))
        allocate(slot(64,n),stencil(64,n),inverse_den(size(xyz,2)),diagonal(n))
        grid_slot = 0
        slot = 0
        stencil = 0.d0
        diagonal = 0.d0
        nodes = 0
        do marker = 1,n
            ix = int(indices(1:4,marker)) - first_x + 1
            iy = int(indices(5:8,marker))
            iz = int(indices(9:12,marker))
            wx = dble(weights(1:4,marker))
            wy = dble(weights(5:8,marker))
            wz = dble(weights(9:12,marker))
            ! Periodic and reflected stencils can visit the same cell more than once.
            call merge_weights(ix,wx)
            call merge_weights(iy,wy)
            call merge_weights(iz,wz)
            s = 0
            do i = 1,4
                x = ix(i)
                if (x < 1 .or. x > size(den,3)) cycle
                do j = 1,4
                    y = iy(j)
                    do k = 1,4
                        z = iz(k)
                        weight = wx(i)*wy(j)*wz(k)
                        if (weight <= 0.d0) cycle
                        node = grid_slot(z,y,x)
                        if (node == 0) then
                            if (.not. ieee_is_finite(den(z,y,x)) .or. den(z,y,x) <= 0.d0) then
                                call parallel_abort('Nonpositive or nonfinite density in ImplicitDirect IBM stencil')
                            endif
                            nodes = nodes + 1
                            node = nodes
                            grid_slot(z,y,x) = node
                            xyz(:,node) = [x,y,z]
                            inverse_den(node) = 1.d0 / den(z,y,x)
                        endif
                        s = s + 1
                        slot(s,marker) = node
                        stencil(s,marker) = weight
                        diagonal(marker) = diagonal(marker) + weight**2 * inverse_den(node)
                    enddo
                enddo
            enddo
        enddo
        deallocate(grid_slot)
        call sum_values(diagonal)
        if (any(diagonal <= 0.d0)) call parallel_abort('ImplicitDirect IBM marker has no Eulerian support')

        allocate(cell(3,nodes),rhs(3,n),residual(3,n),direction(3,n),product(3,n),preconditioned(3,n))
        do node = 1,nodes
            x = xyz(1,node)
            y = xyz(2,node)
            z = xyz(3,node)
            cell(:,node) = velocity(z,y,x,:)
        enddo
        call interpolate(rhs)
        rhs = target - rhs

        ! Solve (W rho^{-1} W^T) correction = target - W velocity.
        ! All ranks hold the marker vectors and own disjoint sets of Eulerian cells.
        correction = 0.d0
        residual = rhs
        do c = 1,3
            preconditioned(c,:) = residual(c,:) / diagonal
        enddo
        direction = preconditioned
        rz = sum(residual * preconditioned,dim=2)
        active = maxval(abs(residual),dim=2) > threshold
        iter = 0
        do while (any(active) .and. iter < max_iterations)
            call spread_correction(direction)
            call interpolate(product)
            do c = 1,3
                if (.not. active(c)) cycle
                curvature = dot_product(direction(c,:),product(c,:))
                if (.not. ieee_is_finite(curvature) .or. curvature <= 0.d0) then
                    call parallel_abort('ImplicitDirect IBM lost positive curvature; check marker spacing and overlap')
                endif
                alpha = rz(c) / curvature
                correction(c,:) = correction(c,:) + alpha * direction(c,:)
                residual(c,:) = residual(c,:) - alpha * product(c,:)
                active(c) = maxval(abs(residual(c,:))) > threshold
                if (active(c)) then
                    preconditioned(c,:) = residual(c,:) / diagonal
                    next_rz(c) = dot_product(residual(c,:),preconditioned(c,:))
                    beta = next_rz(c) / rz(c)
                    direction(c,:) = preconditioned(c,:) + beta * direction(c,:)
                    rz(c) = next_rz(c)
                else
                    direction(c,:) = 0.d0
                endif
            enddo
            iter = iter + 1
        enddo

        ! Re-interpolate the velocity actually passed to collision.
        call spread_correction(correction)
        do node = 1,nodes
            x = xyz(1,node)
            y = xyz(2,node)
            z = xyz(3,node)
            velocity(z,y,x,:) = velocity(z,y,x,:) + cell(:,node)
            cell(:,node) = velocity(z,y,x,:)
        enddo
        call interpolate(product)
        error = maxval(abs(target-product)) / uref
        if (is_root) write(*,'(A,I0,A,ES13.5)') 'ImplicitDirect IBM iterations=',iter,', max slip/Uref=',error
        if (.not. all(ieee_is_finite(product)) .or. .not. ieee_is_finite(error) .or. error > tolerance) then
            call parallel_abort('ImplicitDirect IBM did not reach dtolLBM; check ntolLBM and marker spacing')
        endif

    contains

        subroutine spread_correction(values)
            real(8), intent(in) :: values(:,:)
            integer :: m,entry,p
            cell = 0.d0
            do m = 1,n
                do entry = 1,64
                    p = slot(entry,m)
                    if (p == 0) exit
                    cell(:,p) = cell(:,p) + stencil(entry,m)*values(:,m)
                enddo
            enddo
            do p = 1,nodes
                cell(:,p) = cell(:,p) * inverse_den(p)
            enddo
        end subroutine spread_correction

        subroutine interpolate(values)
            real(8), intent(out) :: values(:,:)
            integer :: m,entry,p
            !$omp parallel do schedule(static) private(m,entry,p)
            do m = 1,n
                values(:,m) = 0.d0
                do entry = 1,64
                    p = slot(entry,m)
                    if (p == 0) exit
                    values(:,m) = values(:,m) + stencil(entry,m)*cell(:,p)
                enddo
            enddo
            !$omp end parallel do
            call sum_matrix(values)
        end subroutine interpolate

    end subroutine solve_implicit_ibm

    subroutine merge_weights(indices,weights)
        integer, intent(in) :: indices(4)
        real(8), intent(inout) :: weights(4)
        integer :: i,j
        do i = 2,4
            do j = 1,i-1
                if (indices(i) /= indices(j)) cycle
                weights(j) = weights(j) + weights(i)
                weights(i) = 0.d0
                exit
            enddo
        enddo
    end subroutine merge_weights

end module ImplicitIBM
