	program driver
	use variable
	use random_version
	use matrix_direct_product
	use DNAmain
	use output
	implicit none

	character(len=256) :: inifile
	
	cluster_num = 2
	

	call get_command_argument(1,inifile)

	call read_ini(trim(adjustl(inifile))) !read *.ini file
	
	length = 9
	call allocatedat ! allocate imagematrix... seq_xxx... and so on
	call readdat !read genome data
	
	call set_seq !random version
	call generate_atcg_matrix_ch(length) ! generate atcg matrix
	
	call DNA_2by2_main(seq_origin) ! main part
!	call DNA_2by2_main(seq_random)
	call cal_entropy_pixel !cal entropy
	call outputfiles_matrix ! output before cluster
	call cluster ! cluster julei
	call outputfiles_cluster !(length) ! output after cluster, before motif
	call cluster_log
	call outputfiles_cluster_log
	
	do motif_len = 2,30,1
		write(*,*) 'l_m = ', motif_len
		call cal_motif_entropy_method3(motif_len)
		call outputfiles_motif
	end do

	end program
