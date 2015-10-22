
!	module in_out_put
!	
!	contains
!================================================================================================	
	subroutine read_ini(inifile)
	use variable
	implicit none

	character(len=*) :: inifile
	!character(len=256) :: information
	!character(len=256) :: datfile,datname,inputpath,outputpath
	character(len=256) :: temp
	!integer :: len_seq

	open(unit=10,file=inifile)
	read(10,'(a)') information
	read(10,*) len_seq,temp
	read(10,'(a)') datfile
	read(10,'(a)') inputpath
	read(10,'(a)') datname
	read(10,'(a)') outputpath
!	read(10,*) cluster_num
!	read(10,*) motif_len
	close(10)

	if ( feedback==1 ) then
		write(*,*) information
		write(*,*) len_seq
		write(*,*) datfile
		write(*,*) outputpath
		write(*,*) cluster_num
		write(*,*) motif_len
	end if
	
	if ( cluster_num > 9 ) then
    	write(*,*) 'cluster_num is too large, set it smaller than 10, or' 
    	write(*,*) '{modified temptemp in routine motif_to_number_c  '
    	write(*,*) 'and write(45,(I1)) imagematrix_cluster(:,i) in output.f90}'
	end if

	end subroutine
!===============================================================================================
	subroutine allocatedat
	use variable
	use matrix_direct_product
	implicit none
		
	allocate(imagematrix(2**length,2**length))
	allocate(imagematrix_norm(2**length,2**length))
	allocate(imagematrix_cluster(2**length,2**length))
	allocate(seq_origin(len_seq))
	allocate(seq_random(len_seq))
	allocate(seq_mirror(len_seq))
	!allocate(atcgmatrix_ch(2**length,2**length))
	allocate(atcgmatrix(2**length,2**length))
!	allocate(  motif_count( cluster_num**(motif_len**2) )  )	
!	allocate(  motif_count( (2**length-motif_len+1)**2 ) )	

	end subroutine	
!===============================================================================================
	subroutine allocatedat_minimal
	use variable
	use matrix_direct_product
	implicit none
	
	!useful when the ram is not enough
		
	allocate(imagematrix(2**length,2**length))
	!allocate(imagematrix_norm(2**length,2**length))
	!allocate(imagematrix_cluster(2**length,2**length))
	allocate(seq_origin(len_seq))
	allocate(seq_random(len_seq))
	!allocate(seq_mirror(len_seq))
	!allocate(atcgmatrix_ch(2**length,2**length))
	!allocate(atcgmatrix(2**length,2**length))	

	end subroutine	
!===============================================================================================
	subroutine readdat
	use variable
	implicit none

	integer :: i
	
	write(*,*) 'read a DNA sequence in file: '
	write(*,*) trim(adjustl(datfile)) 
			
	open(unit=11,file=trim(adjustl(datfile)))
!	read(11,'(a)') seq_origin
	do i=1,len_seq,1
	 	read(11,'(A1)',advance='NO') seq_origin(i)
	end do
 	close(11)
 	
 	write(*,*) seq_origin(1:100),'......'
 	
 	end subroutine
!===============================================================================================

	subroutine release_buffer_matrix
	use variable
	use matrix_direct_product
	implicit none
	
	!dellocate(imagematrix,imagematrix_cluster,imagematrix_norm)
	deallocate(imagematrix)
	deallocate(imagematrix_norm)
	deallocate(imagematrix_cluster)
	deallocate(atcgmatrix)
	
	return	
	end subroutine
!===============================================================================================
	subroutine release_buffer_all
	implicit none
	
	!dellocate(imagematrix,imagematrix_cluster,imagematrix_norm,seq_origin,seq_random,seq_mirror)
	end subroutine
!===============================================================================================
!	end module in_out_put
