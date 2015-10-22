	module output
	implicit none
	
	character(len=256) :: no_suffix
	
	contains
!===================================================================================================
	subroutine outputfiles_matrix
	use variable
	use random_version
	implicit none
	
	character(len=32) :: date_time
	integer :: i
	character(len=2) :: temp
	
	call name_trunc(datname,no_suffix)
	
	write(temp,'(I2)') length 
	
	open(unit=40,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".image"//trim(adjustl(temp)))
	write(*,*) 'writing to file: ',trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".image"
	do i=1,2**length,1
		write(40,*)	 imagematrix(:,i)
	end do
	close(40)
	
	Call FDate(date_time)
	
	open(unit=41,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".info",position='APPEND')
	write(41,*) 'files of ',trim(adjustl(no_suffix)),' were generated at time: ',date_time	
	write(41,*) information
	write(41,*) 'atcg number count and their prob: ',countnum,'||',prob 
	write(41,*) 'writing to file: ',trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".image",length
	close(41)
	
	open(unit=42,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".entropy",position='APPEND')
	write(42,*) length,entropy
	close(42)
	
	write(*,*) 'output files before cluster: finished'
	
	end subroutine
!=================================================================================================
	subroutine outputfiles_cluster
	use variable
	use matrix_direct_product
	implicit none
	
	integer :: i!,leng
	character(len=2) :: temp,temp1 !
	character(len=256) :: finame
	!character(leng) :: atcgmatrix_ch(2**leng,2**leng)
	!common /n1/ atcgmatrix_ch
	
!	write(*,*) cluster_num
	write(temp,'(I2)') length
	finame=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//'.'//trim(adjustl(temp))//'cluster'	
	write(temp1,'(I2)') cluster_num
!	write(*,*) trim(adjustl(temp1))
	finame=trim(adjustl(finame))//trim(adjustl(temp1))
	
	
	open(unit=43,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".info",position='APPEND')
	write(43,*) 'max value: ',max_val,'string: ',atcgmatrix(max_loc(1),max_loc(2))
	write(43,*) 'writing to file: ',trim(adjustl(finame))
	close(43)
	
	open(unit=44,file=TRIM(adjustl(finame)))
	do i=1,2**length,1
		write(44,*)	 imagematrix_cluster(:,i)
	end do
	close(44)
	
	write(*,*) 'output files after cluster before motif: finished'
	
	end subroutine
!=================================================================================================
	subroutine outputfiles_cluster_log
	use variable
	use matrix_direct_product
	implicit none
	
	integer :: i!,leng
	character(len=2) :: temp,temp1 !
	character(len=10000) :: temp2
	character(len=256) :: finame
	!character(leng) :: atcgmatrix_ch(2**leng,2**leng)
	!common /n1/ atcgmatrix_ch
	
	write(temp,'(I2)') length
	finame=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//'.'//trim(adjustl(temp))//'logcluster'
	write(temp1,'(I2)') cluster_num
	finame=trim(adjustl(finame))//trim(adjustl(temp1))
	
	open(unit=31,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".info",position='APPEND')
	write(31,*) 'max value: ',max_val,'string: ',atcgmatrix(max_loc(1),max_loc(2))
	write(31,*) 'writing to file: ',trim(adjustl(finame))
	close(31)
	
	open(unit=45,file=finame)
	do i=1,2**length,1
		write(45,*)	 imagematrix_cluster(:,i)
	end do
	close(45)
	
	write(*,*) 'output files after cluster before motif log: finished'
	
	end subroutine
!=================================================================================================
	subroutine outputfiles_motif
	use variable
	implicit none
	
	integer :: i
	
	write(*,*) 'this part is done in the module DNAmain, subroutine cal_motif_entropy'	
!	open(unit=46,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".info",position='APPEND')
!	write(46,*) 'motif size: ',motif_len,' by ',motif_len
!	write(46,*) 'cluster into: ',cluster_num,' class'
!	write(46,*) 'the most motif is: ',maxloc(motif_count),' | it appears ',maxval(motif_count),' times'
!	close(46)
!	
!	open(unit=47,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".motif")
!	do i =1,size(motif_count),1
!		write(47,*) i, motif_count(i)
!	end do
!	close(47)
!	
!	open(unit=48,file=trim(adjustl(outputpath))//'/'//trim(adjustl(no_suffix))//".motifentropy")
!	write(48,*) length,cluster_num,motif_len,motif_entropy
!	close(48)	
!	
!	write(*,*) 'output files after motif: finished'
	
	end subroutine
!=================================================================================================
	subroutine name_trunc(string,no_suffix)
	implicit none
	
	character(len=*) :: string
	character(len=*) :: no_suffix
	integer :: pos
	
	pos=index(string,'.')
	!allocate(no_suffix(pos))
	no_suffix=string(1:(pos-1))
	
	end subroutine
!===================================================================================================
	end module
