! the code is modified by hongbo base on:
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! sort character(len=leng) array

! Made F conformant by Walt Brainerd

!module qsort_c_module

!implicit none
!public :: QsortC
!private :: Partition

!contains

recursive subroutine QsortC(A,leng)
  character(len=leng), intent(in out), dimension(:) :: A
  integer :: iq
  integer :: leng

  if(size(A) > 1) then
     call Partition(A, iq,leng)
     call QsortC(A(:iq-1),leng)
     call QsortC(A(iq:),leng)
  endif
end subroutine QsortC

subroutine Partition(A, marker,leng)
  character(len=leng), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j,leng
  character(len=leng) :: temp
  character(len=leng) :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

!end module qsort_c_module

!program sortdriver
!  use qsort_c_module
!  implicit none
!  integer, parameter :: r = 10
! ! real, dimension(1:r) :: myarray = &        ! (1:r)
! !    (/0, 50, 20, 25, 90, 10, 5, 1, 99, 75/)
!  character(len=2),dimension(1:r) :: myarray = & 
!  & (/'00', '50', '20', '25', '90', '10', '05', '01', '99', '75'/)
!  print *, "myarray is ", myarray
!  call QsortC(myarray,2)
!  print *, "sorted array is ", myarray
!end program sortdriver

