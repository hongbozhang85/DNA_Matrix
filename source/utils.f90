! the file contains several module : utils, matrix_direct_product, search_substitude




	module utils
	implicit none

	contains
!========================================================================
	subroutine writeinMatrix(dimen,matrix,unitunit) 
	implicit none

! output a Matrix in file or terminal in format of 2 dimension instead of 1 dimension. 

	integer,intent(in) :: dimen,matrix(dimen,dimen)
	integer,optional :: unitunit
	integer :: i

	do i=1,dimen,1
		write(unitunit,*)	 matrix(:,i)
	end do

	end subroutine


	!===========================================================================

	subroutine writeinMatrix_ch(dimen,matrix,unitunit) 
	implicit none

! output a Matrix in file or terminal in format of 2 dimension instead of 1 dimension. 

	integer,intent(in) :: dimen
	character(len=*),intent(in) :: matrix(dimen,dimen)
	integer,optional :: unitunit
	integer :: i

	do i=1,dimen,1
		write(unitunit,*)	 matrix(:,i)
	end do

	end subroutine

	!===========================================================================

	subroutine writeinMatrix_re(dimen,matrix,unitunit) 
	implicit none

! output a Matrix in file or terminal in format of 2 dimension instead of 1 dimension. 

	integer,intent(in) :: dimen
	real,intent(in) :: matrix(dimen,dimen)
	integer,optional :: unitunit
	integer :: i

	do i=1,dimen,1
		write(unitunit,*)	 matrix(:,i)
	end do

	end subroutine
!===============================================================================
	end module utils


!===========================================================================
!===========================================================================
!===========================================================================

	module search_and_substitude
	implicit none

	contains
!===================================================================
	subroutine atcg_to_number(leng,chain_atcg,chain_1234)
	implicit none

	! transfer 'a/t/c/g' to '1/2/3/4'
	integer :: leng !len(chain_atcg)
	character(len=leng),intent(in) :: chain_atcg
	integer,intent(out) :: chain_1234
	character(len=leng) :: buffer
	!character(len=10) :: length_ch
	
	call search_substitute(leng,chain_atcg,buffer,'A','1')
	call search_substitute(leng,buffer,buffer,'T','2')
	call search_substitute(leng,buffer,buffer,'C','3')
	call search_substitute(leng,buffer,buffer,'G','4')

	!write(length_ch,'(I10)') length
	read(buffer,*) chain_1234

	end subroutine


	!====================================================================


	subroutine number_to_atcg(leng,chain_atcg,chain_1234)
	implicit none

	! transfer '1/2/3/4' to 'a/t/c/g'
	integer :: leng !len(chain_atcg)
	character(len=leng),intent(out) :: chain_atcg
	integer(kind=8),intent(in) :: chain_1234
	character(len=leng) :: buffer
	character(len=10) :: leng_ch

	write(leng_ch,'(I10)') leng
	write(buffer,'(I'//trim(adjustl(leng_ch))//')') chain_1234
	
	call search_substitute(leng,buffer,buffer,'1','A')
	call search_substitute(leng,buffer,buffer,'2','T')
	call search_substitute(leng,buffer,buffer,'3','C')
	call search_substitute(leng,buffer,chain_atcg,'4','G')

	!call search_substitute(length,buffer,buffer,'1','A')
	!call search_substitute(length,buffer,buffer,'2','G')
	!call search_substitute(length,buffer,buffer,'3','C')
	!call search_substitute(length,buffer,chain_atcg,'4','T')

	end subroutine


	!====================================================================

	subroutine search_substitute(leng,instring,outstring,search,substitute)
	implicit none

	integer :: leng
	character(len=leng) :: instring
	character(len=leng) :: outstring
	character,intent(in) :: search,substitute
	integer :: i

	! for example 
	! search=a, substitute=1: 'atacg'->'1t1cg'

	outstring = instring
	do i=1,leng,1
		 if ( outstring(i:i)==search ) outstring(i:i)=substitute
	end do
	
	end subroutine 

	!======================================================================
	end module search_and_substitude


!=================================================================================
!=================================================================================
!=================================================================================

	module matrix_direct_product
	use variable
    implicit none
	
	!character(length),allocatable :: atcgmatrix_ch(:,:)
	integer(kind=8),allocatable :: atcgmatrix(:,:)
	
	contains 
!==================================================================

!==================================================================

 	subroutine generate_atcg_matrix_ch(leng)
 	use variable
 	use search_and_substitude
 	implicit none
 !======== generate an ATCG matrix (matrix direct-product)
 	integer :: i,j,leng
 	character(leng) :: atcgmatrix_ch(2**leng,2**leng)
 	!common /n1/ atcgmatrix_ch
 	
    write(*,*) 'generate atcg matrix'
    call generate_atcg_matrix(length,atcgmatrix)
 
 	write(*,*) 'the atcg matrix is (in 1234)'
 	!call writeinMatrix(2**length,atcgmatrix,6)
 
   !======== tranform matrix from 1234 to atcg
 
 	do i=1,2**length,1
 		do j=1,2**length,1
 			call number_to_atcg(length,atcgmatrix_ch(i,j),atcgmatrix(i,j))
 		end do
 	end do
 
 	write(*,*) 'the atcg matrix is (in atcg)'
 	write(*,*) atcgmatrix_ch(2,1)
 	!call writeinMatrix_ch(2**length,atcgmatrix_ch,6)
  end subroutine

!===================================================================
	subroutine generate_atcg_matrix(leng,matrix)							 
	implicit none

	!generate atcg matrix.	such as
  !
	! a | t
	!	-----  1 dim
	! c | g 
	!
	!	aa  at | ta tt
	! ac  ag | tc	tg
	!---------------- 2 dim
	!	ca ct  | ga gt
	!	cc cg	 | gc	gg
	!
	!.... 3dim, 4dim....

	integer(kind=8) :: matrix(2**leng,2**leng),initial(2**leng,2**leng) !length is len(attc...g)
	!integer :: row,column,i !matrix(row,column)
	integer :: leng,i
	initial(1,1)=1
	initial(2,1)=2
	initial(1,2)=3
	initial(2,2)=4
	
	matrix=0

	if (leng==1) then
		matrix=initial
	else
		do i=1,leng-1,1
			call small_to_big(i,leng,initial,matrix)
			initial=matrix
		end do
	end if

	end subroutine


!=================================================================

	subroutine small_to_big(times,leng,small,big)
	implicit none

	! small matrix is n*n, big matrix is 2n*2n
	! map every one element in small matrix to 2 by 2 sub-matrix in big matrix
	!
	! for example
	!
	! a | t
	!	-----  1 dim
	! c | g 
	!
	!	aa  at | ta tt
	! ac  ag | tc	tg
	!---------------- 2 dim
	!	ca ct  | ga gt
	!	cc cg	 | gc	gg
	!
	! map a to 2 by 2 |aa at|
	!									|ac ag|

	integer(kind=8) ::  small(2**leng,2**leng), big(2**leng,2**leng)
	!length is len(at..cg) in small matrix 
	integer :: leng,times,row_s, column_s,row_b,column_b !small(row_s,column_s)

	do row_s=1,2**times,1
		do column_s=1,2**times,1
			row_b = 2*row_s-1
			column_b = 2*column_s-1
			big(row_b,column_b) = 10*small(row_s,column_s)+1
			big(row_b+1,column_b)= 10*small(row_s,column_s)+2
			big(row_b,column_b+1) = 10*small(row_s,column_s)+3
			big(row_b+1,column_b+1) = 10*small(row_s,column_s)+4
		end do
	end do

	end subroutine

	!======================================================================

	subroutine copy_big_to_small(length_s,length_b,small,big) !useless
	implicit none

	!copy left-up part of large matrix(none-zero elements in this question) to small matrix, 
	!such as
	!
	! copy 8 by 8 matrix:
	!
	!	aa  at | ta tt	0000000
	! ac  ag | tc	tg	0000000
	!---------------- 
	!	ca ct  | ga gt	0000000
	!	cc cg	 | gc	gg	0000000
	! 00000000000000000000000
	! 00000000000000000000000
	!
	! to 4 by 4 matrix
	!
	!	aa  at | ta tt
	! ac  ag | tc	tg
	!----------------
	!	ca ct  | ga gt
	!	cc cg	 | gc	gg

	integer :: length_s,length_b, small(2**length_s,2**length_s), big(2**length_b,2**length_b)
	!length is len(at..cg)
	integer :: i,j

	do i=1,2**length_s,1
		do j=1,2**length_s,1
		small(i,j)=big(i,j)
		end do
	end do

	end subroutine

	!===========================================================================
	end module matrix_direct_product

	
