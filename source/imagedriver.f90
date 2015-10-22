	program image
	use variable
	use random_version
	use matrix_direct_product
	use DNAmain
	use output
	implicit none

	character(len=256) :: inifile
	integer :: i

	cluster_num = 8
	write(*,*) cluster_num

	call get_command_argument(1,inifile)

	call read_ini(trim(adjustl(inifile))) !read *.ini file
	allocate(seq_origin(len_seq))
	allocate(seq_random(len_seq))
	allocate(seq_mirror(len_seq))

	call readdat !read genome data
	call set_seq
	
	do i = 2,13,1
		length=i
		write(*,*) '**************'
		write(*,*) 'length= ',length
		write(*,*) '**************'
		allocate(imagematrix(2**length,2**length))
		allocate(imagematrix_norm(2**length,2**length))
		allocate(imagematrix_cluster(2**length,2**length))
		allocate(atcgmatrix(2**length,2**length))
		
		call generate_atcg_matrix_ch(length) ! generate atcg matrix
	
		call DNA_2by2_main(seq_origin) ! main part
		call cal_entropy_pixel !cal entropy
		call outputfiles_matrix ! output before cluster
		call cluster ! cluster julei
		call outputfiles_cluster !(length) ! output after cluster, before motif
		call cluster_log
		call outputfiles_cluster_log
		
		call release_buffer_matrix
	end do
	
!	call cal_motif_entropy_method3(motif_len)
!	call outputfiles_motif

	end program
