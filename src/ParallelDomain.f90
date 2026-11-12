! SPDX-License-Identifier: GPL-3.0-or-later
module ParallelDomain
#ifdef FSILBM3D_MPI
    use mpi
#endif
    use, intrinsic :: iso_fortran_env, only: error_unit
    implicit none
    private
    public :: rank_id, rank_count, is_root, start_parallel, finish_parallel, parallel_abort, parallel_barrier
    public :: mpi_enabled
    public :: slab_bounds, exchange_plane, sum_values, sum_matrix, sum_face, max_values, gather_field, gather_checkpoint
    integer, protected :: rank_id = 0, rank_count = 1
    logical, protected :: is_root = .true.
#ifdef FSILBM3D_MPI
    logical, parameter :: mpi_enabled = .true.
#else
    logical, parameter :: mpi_enabled = .false.
#endif

contains

    subroutine start_parallel()
#ifdef FSILBM3D_MPI
        integer :: ierr, provided
        call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank_id, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)
        is_root = rank_id == 0
        if (provided < MPI_THREAD_FUNNELED) call parallel_abort('MPI_THREAD_FUNNELED is required')
#endif
    end subroutine

    subroutine sum_matrix(values)
        real(8), intent(inout) :: values(:,:)
#ifdef FSILBM3D_MPI
        integer :: ierr
        if (rank_count > 1) call MPI_Allreduce(MPI_IN_PLACE, values, size(values), &
                                               MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    end subroutine

    subroutine sum_face(values)
        real(8), intent(inout) :: values(:,:,:)
#ifdef FSILBM3D_MPI
        integer :: ierr
        if (rank_count > 1) call MPI_Allreduce(MPI_IN_PLACE, values, size(values), &
                                               MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    end subroutine

    subroutine finish_parallel()
#ifdef FSILBM3D_MPI
        integer :: ierr
        call MPI_Finalize(ierr)
#endif
    end subroutine

    subroutine parallel_abort(message)
        character(*), intent(in) :: message
#ifdef FSILBM3D_MPI
        integer :: ierr
#endif
        write(error_unit,'(A,I0,2A)') 'Rank ', rank_id, ': ', message
#ifdef FSILBM3D_MPI
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
#endif
        stop 1
    end subroutine

    subroutine parallel_barrier()
#ifdef FSILBM3D_MPI
        integer :: ierr
        if (rank_count > 1) call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    end subroutine

    subroutine slab_bounds(nx, rank, first, last)
        integer, intent(in) :: nx, rank
        integer, intent(out) :: first, last
        integer :: width, extra
        width = nx / rank_count
        extra = mod(nx, rank_count)
        first = 1 + rank * width + min(rank, extra)
        last = first + width - 1
        if (rank < extra) last = last + 1
    end subroutine

    subroutine exchange_plane(send, receive, direction)
        real(8), intent(in) :: send(:,:)
        real(8), intent(out) :: receive(:,:)
        integer, intent(in) :: direction
#ifdef FSILBM3D_MPI
        integer :: target, source, ierr
        target = modulo(rank_id + direction, rank_count)
        source = modulo(rank_id - direction, rank_count)
        call MPI_Sendrecv(send, size(send), MPI_DOUBLE_PRECISION, target, 0, &
                          receive, size(receive), MPI_DOUBLE_PRECISION, source, 0, &
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
#else
        receive = send
#endif
    end subroutine

    subroutine sum_values(values)
        real(8), intent(inout) :: values(:)
#ifdef FSILBM3D_MPI
        integer :: ierr
        if (rank_count > 1) call MPI_Allreduce(MPI_IN_PLACE, values, size(values), &
                                               MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
    end subroutine

    subroutine max_values(values)
        real(8), intent(inout) :: values(:)
#ifdef FSILBM3D_MPI
        integer :: ierr
        if (rank_count > 1) call MPI_Allreduce(MPI_IN_PLACE, values, size(values), &
                                               MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
#endif
    end subroutine

    subroutine gather_field(local, global, nx, trim_x)
        real(4), intent(in) :: local(:,:,:)
        real(4), intent(out) :: global(:,:,:)
        integer, intent(in) :: nx, trim_x
#ifdef FSILBM3D_MPI
        integer :: ierr
        integer :: counts(0:rank_count-1), offsets(0:rank_count-1)
        call gather_layout(nx, trim_x, int(size(local,1),8)*size(local,2), counts, offsets)
        call MPI_Gatherv(local, size(local), MPI_REAL, global, counts, offsets, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
#else
        global = local
#endif
    end subroutine

    subroutine gather_checkpoint(local, global, nx)
        real(8), intent(in) :: local(:,:,:)
        real(8), intent(out) :: global(:,:,:)
        integer, intent(in) :: nx
#ifdef FSILBM3D_MPI
        integer :: ierr
        integer :: counts(0:rank_count-1), offsets(0:rank_count-1)
        call gather_layout(nx, 0, int(size(local,1),8)*size(local,2), counts, offsets)
        call MPI_Gatherv(local, size(local), MPI_DOUBLE_PRECISION, global, counts, offsets, &
                         MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#else
        global = local
#endif
    end subroutine

#ifdef FSILBM3D_MPI
    subroutine gather_layout(nx, trim_x, plane_size, counts, offsets)
        integer, intent(in) :: nx, trim_x
        integer(8), intent(in) :: plane_size
        integer, intent(out) :: counts(0:), offsets(0:)
        integer :: r, first, last
        if (plane_size * (nx-2*trim_x) > int(huge(1),8)) then
            call parallel_abort('Output exceeds the MPI_Gatherv count limit')
        endif
        do r=0,rank_count-1
            call slab_bounds(nx, r, first, last)
            first = max(first, 1+trim_x)
            last = min(last, nx-trim_x)
            counts(r) = int(max(0, last-first+1) * plane_size)
            offsets(r) = int(min(first-1-trim_x, nx-2*trim_x) * plane_size)
        enddo
    end subroutine
#endif
end module ParallelDomain
