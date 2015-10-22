	subroutine cluster
	use variable
	implicit none
	! cluster 2^n. e.g., n=1, black and white
	
	real :: node(cluster_num)
	real :: temp
	integer :: i,j,k
	
	write(*,*) 'cluster not_log'
	
	max_val = maxval(imagematrix)
	max_loc = maxloc(imagematrix)
	
	do i=1,cluster_num,1
		node(i)= max_val/cluster_num*i
	end do
	
	temp=(max_val+0.0)/cluster_num	
	
	do i=1,2**length,1
		do j=1,2**length,1
			imagematrix_cluster(i,j)=floor(abs(imagematrix(i,j)-0.1)/temp)
		end do
	end do
	
	
		
	end subroutine
!====================================================================
	subroutine cluster_log
	use variable
	implicit none
	! cluster 2^n. e.g., n=1, black and white
	
	!real :: node(cluster_num)
	real :: temp
	integer :: i,j,k
	
	write(*,*) 'cluster log'
	
	max_val = maxval(imagematrix)
	max_loc = maxloc(imagematrix)
	
!	do i=1,cluster_num,1
!		node(i)= max_val/cluster_num*i
!	end do
	
	temp=log(max_val+0.0)/cluster_num	
	
	imagematrix_cluster = 0
	
	do i=1,2**length,1
		do j=1,2**length,1
			if ( imagematrix(i,j) .ne. 0 ) then
			imagematrix_cluster(i,j)=floor(abs(log(imagematrix(i,j)+0.0)-0.0001)/temp)
			endif
		end do
	end do
	
	end subroutine
!=======================================================================
	subroutine cluster_log_minimal
	use variable
	implicit none
	
	real :: temp
	integer :: i,j,k
	
	write(*,*) 'cluster log'
	
	max_val = maxval(imagematrix)
	max_loc = maxloc(imagematrix)
	
	temp=log(max_val+0.0)/cluster_num	
		
	do i=1,2**length,1
		do j=1,2**length,1
			if ( imagematrix(i,j) .ne. 0 ) then
			imagematrix(i,j)=floor(abs(log(imagematrix(i,j)+0.0)-0.0001)/temp)
			endif
		end do
	end do
	
	end subroutine

