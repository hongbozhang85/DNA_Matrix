	program strscurve
	use variable
	use random_version
	use matrix_direct_product
	use DNAmain
	use output
	implicit none
	
	!pixel level s curve

	character(len=256) :: inifile
	
	call get_command_argument(1,inifile)

	call read_ini(trim(adjustl(inifile))) !read *.ini file
	allocate(seq_origin(len_seq))
	allocate(seq_random(len_seq))
!	allocate(seq_mirror(len_seq))
	call readdat !read genome data
	call generate_random_seq(len_seq,seq_origin,seq_random)
	
	do length = 1, 13, 1
		write(*,*) 'l_s= ',length
		allocate(imagematrix(2**length,2**length))
!		call DNA_2by2_main(seq_origin) ! main part
		call DNA_2by2_main(seq_random) ! main part
		call cal_entropy_pixel_minimal
		call outputfiles_matrix ! output before cluster
		deallocate(imagematrix)
	end do
	

	end program
