  module Task
    use mpi
    implicit none
    contains
            
        subroutine GetMaxCoordinates(A, x1, y1, x2, y2)
            implicit none
            real(8), intent(in), dimension(:,:) :: A
            integer(4), intent(out) :: x1, y1, x2, y2
            integer(4) :: n, L, R, Up, Down, m, tmp
            real(8), allocatable :: current_column(:), B(:,:)
            real(8) :: current_sum, max_sum
            logical :: transpos
            integer :: ierr, s1ze, rank, myop, extent, blockcounts(2), offsets(2), types(2), mytype

            call MPI_Comm_rank (MPI_COMM_WORLD, rank, ierr)
            call MPI_Comm_size (MPI_COMM_WORLD, s1ze, ierr)

            m = size(A, dim=1) 
            n = size(A, dim=2) 
            allocate (current_column(m))
            transpos = .FALSE.

            if (m < n) then 
                allocate (B(n, m))
                transpos = .TRUE.   
                B = transpose(A)
                m = size(B, dim=1) 
                n = size(B, dim=2) 
            else
                allocate(B(m, n))
                B = A     
            endif

            max_sum=B(1, 1)
            x1=1
            y1=1
            x2=1
            y2=1

            do L=1+rank, n, s1ze
                current_column = B(:, L)  

                do R=L, size(current_column)
                    if (R > L) then 
                        current_column = current_column + B(:, R)
                    endif

                    call FindMaxInArray (current_column, current_sum, Up, Down)

                    if (current_sum > local_res%local_max_sum) then
                         local_res%local_max_sum = current_sum
                         local_res%bestX1 = Up
                         local_res%bestX2 = Down
                         local_res%bestY1 = L
                         local_res%bestY2 = R
                    endif
                end do
            end do


            if (rank == 0) then                   
                max_sum = res%local_max_sum
                x1 = res%bestX1
                y1 = res%bestY1
                x2 = res%bestX2
                y2 = res%bestY2
            end if


            if (transpos) then  
                tmp = x1
                x1 = y1
                y1 = tmp

                tmp = y2
                y2 = x2
                x2 = tmp
            endif

            deallocate (current_column, B)

        end subroutine GetMaxCoordinates

        subroutine FindMaxInArray(a, Sum, Up, Down)
            real(8), intent(in), dimension(:) :: a
            integer(4), intent(out) :: Up, Down
            real(8), intent(out) :: Sum
            real(8) :: cur_sum
            integer(4) :: minus_pos, i

            Sum = a(1)
            Up = 1
            Down = 1
            cur_sum = 0
            minus_pos = 0

            do i=1, size(a)
                cur_sum = cur_sum + a(i)
                if (cur_sum > Sum) then
                    Sum = cur_sum
                    Up = minus_pos + 1
                    Down = i
                endif

                if (cur_sum < 0) then
                    cur_sum = 0
                    minus_pos = i
                endif

            enddo
        end subroutine FindMaxInArray

  end module Task
