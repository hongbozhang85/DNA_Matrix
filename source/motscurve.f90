! motif level S curve, S(l_s,l_m,n_c)
	program motscurve
	use variable
	use random_version
	use matrix_direct_product
	use DNAmain
	use output
	implicit none

	character(len=256) :: inifile
	integer :: i

	call get_command_argument(1,inifile)
	call read_ini(trim(adjustl(inifile))) !read *.ini file
	allocate(seq_origin(len_seq))
	allocate(seq_random(len_seq))
	call readdat !read genome data
	call generate_random_seq(len_seq,seq_origin,seq_random)
	
	
	length = 9
	motif_len = 2
!	cluster_num = 4
	
	do cluster_num=2,9,1
	
		write(*,*) 'l_m= ',motif_len,'l_s= ',length,'n_c= ',cluster_num		
		if ( 2**length < motif_len ) stop 'ls is too small for current lm'

		allocate(imagematrix(2**length,2**length))
	
		call DNA_2by2_main(seq_origin) ! main part
	
		call cluster_log_minimal
	
		call cal_motif_entropy_no_output(motif_len)
		
		deallocate(imagematrix)
	end do	


	end program
