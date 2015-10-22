

	module variable
	implicit none

	
	character(len=4),parameter :: order='ATCG'
	integer :: feedback = 1

	!type :: DNA_Data
		integer :: length ! atcg string length , e.g., length=9, AAAAAAAAA
		integer :: len_seq ! length of a given genome sequence
		character(len=256) :: information,datfile,datname,inputpath,outputpath
		integer,allocatable :: imagematrix(:,:) ! pixel level. 2dim image 
		integer,allocatable :: imagematrix_cluster(:,:)
		real,allocatable :: imagematrix_norm(:,:)
		real :: entropy ! pixel level, not clustered. sum_{i} p_i ln p_i, # of {i} is 4**length
		character,allocatable :: seq_origin(:),seq_random(:),seq_mirror(:)
		integer :: cluster_num, max_val,max_loc(2) ! the longest string and its frequency,e.g.,ATATATATA 27813
		! cluster_num : # of julei
		integer :: motif_len ! after clustered, motif_len by motif_len
!		integer,allocatable :: motif_count(:) ! count # of different motifs
		!character,allocatable :: motif_count(:)
		real :: motif_entropy

		
				
	!end type DNA_Data


	end module variable
