
	module DNAmain
	use variable
	!use random_version
	
	contains
	!=====================================================================
	!=====================================================================
	subroutine DNA_2by2_main(seq) !(length,len_seq,seq,imagematrix,imagematrix_norm)
	implicit none

	!generate 2-dim image of DNA

	!integer :: imagematrix(2**length,2**length)
	!real :: imagematrix_norm(2**length,2**length)
	character :: seq(len_seq)
	integer :: i,j,k
	integer :: coor_row,coor_col ! coor_col is hang/row; and coor_row is lie/col

	imagematrix=0

	write(*,*) 'main part of the program, generate 2-dim image-matrix of DNA'
	
	coor_row=1
	coor_col=1
	
	do i=1,length,1
		if ( seq(i)== 'A' ) then
			coor_row = coor_row + 0
			coor_col = coor_col + 0
		elseif ( seq(i)== 'T') then
			coor_row = coor_row + 2**(length-i)
			coor_col = coor_col+ 0
		elseif ( seq(i)=='C' ) then
			coor_row = coor_row + 0
			coor_col = coor_col+ 2**(length-i)
		elseif ( seq(i)=='G' ) then
			coor_row = coor_row + 2**(length-i)
			coor_col = coor_col+ 2**(length-i)
		end if
	end do

	imagematrix(coor_row,coor_col) = imagematrix(coor_row,coor_col) + 1


	 
	do k=2,len_seq-length+1,1

		if ( seq(k-1)== 'A' ) then
			coor_row = coor_row - 0
			coor_col = coor_col - 0
		elseif ( seq(k-1)== 'T') then
			coor_row = coor_row - 2**(length-1)
			coor_col = coor_col - 0
		elseif ( seq(k-1)=='C' ) then
			coor_row = coor_row - 0
			coor_col = coor_col - 2**(length-1)
		elseif ( seq(k-1)=='G' ) then
			coor_row = coor_row - 2**(length-1)
			coor_col = coor_col - 2**(length-1)
		end if
		
		coor_row = 2*coor_row - 1
		coor_col = 2*coor_col - 1

		if ( seq(k+length-1)== 'A' ) then
			coor_row = coor_row + 0
			coor_col = coor_col + 0
		elseif ( seq(k+length-1)== 'T') then
			coor_row = coor_row + 1
			coor_col = coor_col+ 0
		elseif ( seq(k+length-1)=='C' ) then
			coor_row = coor_row + 0
			coor_col = coor_col+ 1
		elseif ( seq(k+length-1)=='G' ) then
			coor_row = coor_row + 1
			coor_col = coor_col+ 1
		end if

		imagematrix(coor_row,coor_col) = imagematrix(coor_row,coor_col) + 1
		
	end do

	

	!write(*,*) 'the imagematrix:'
	!call writeinMatrix(2**length,imagematrix,6)
	!write(*,*) 'the normalized imagematrix:'
	!call writeinMatrix_re(2**length,imagematrix_norm,6)
	!!write(11,*) imagematrix_norm
	!!call writeinMatrix_re(2**length,imagematrix_norm,11)
	

	end subroutine
!============================================================================

	subroutine cal_entropy_pixel !(entropy,imagematrix,imagematrix_norm)
	implicit none
	
	! calculate entropy. \sum_{i} p_i ln p_i
	! here {i} are elements of atcgmatrix_ch
	
	!real :: entropy
	!integer :: imagematrix(2**length,2**length)
	!real :: imagematrix_norm(2**length,2**length)
	integer :: i,j
	
	imagematrix_norm = imagematrix/(len_seq-length+1.0)
		
	entropy = 0.0
	do i=1,2**length,1
		do j=1,2**length,1
			if ( imagematrix(i,j)/=0	) then
				entropy = entropy - imagematrix_norm(i,j)*log(imagematrix_norm(i,j))/log(2.0)
			end if
		end do
	end do

	write(*,*) 'the entropy \sum_{i} p_i ln p_i, {i} is # of atcg... # of {i} ~ 4**length'
	write(*,*) entropy

	end subroutine
!============================================================================

	subroutine cal_entropy_pixel_minimal !(entropy,imagematrix,imagematrix_norm)
	implicit none
	
	integer :: i,j
	
	!useful when the ram is not enough
		
	entropy = 0.0
	do i=1,2**length,1
		do j=1,2**length,1
			if ( imagematrix(i,j)/=0	) then
 entropy=entropy-imagematrix(i,j)/(len_seq-length+1.0)*log(imagematrix(i,j)/(len_seq-length+1.0))/log(2.0)
			end if
		end do
	end do

	write(*,*) 'the entropy \sum_{i} p_i ln p_i, {i} is # of atcg... # of {i} ~ 4**length'
	write(*,*) entropy

	end subroutine
!===============================================================================
	subroutine cal_motif_entropy_method3(leng)
	use output
	implicit none

	
	!size of motif_count is (2**length-motif_len+1)**2
	
	integer :: i,j,leng,iprev,numcount
	integer :: motif(leng,leng)
	!!!real :: numb
	character(len=leng*leng) :: numb
!	character(len=leng*leng) :: numbtemp
	real :: temp
	character(len=256) :: nosuf
	!!!real,allocatable :: motif_count_array(:,:)
	character(len=leng*leng),allocatable :: motif_count(:)
	character(len=leng*leng),allocatable :: motif_pattern(:)
	integer,allocatable :: motif_pattern_count(:)
	character(len=7) :: temptemp,temp3
	character(len=128) :: formatstring
!	character(len=leng*leng) :: debugstring(100)
	
	allocate( motif_count( (2**length-motif_len+1)**2 ) )
	call name_trunc(datname,nosuf)
	
	write(temptemp,'(I7)') motif_len*motif_len
	formatstring='(A'//trim(adjustl(temptemp))//',I10)'
	write(temp3,'(I7)') motif_len

		
!	motif_count = 0
	motif_count=''
	do i=1,2**length-motif_len+1,1
		do j=1,2**length-motif_len+1,1
			call cut_motif(i,j,leng,motif)
			call motif_to_number_c(motif,leng,numb)
			motif_count(j+(i-1)*(2**length-motif_len+1))=numb
		end do
	end do
	
!	debugstring = motif_count(1:100)
!	
	open(unit=30,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".templm"//trim(adjustl(temp3)))
!	do i=1,(2**length-motif_len+1)**2,1
	do i=1,(2**length-motif_len+1)**2,1
		write(30,trim(adjustl(formatstring))) motif_count(i),i
	end do
	close(30)
	
	call QsortC(motif_count,leng*leng)
!	call QsortC(debugstring,4)
	
	open(unit=31,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".temp2lm"//trim(adjustl(temp3)))
!	do i=1,(2**length-motif_len+1)**2,1
	do i=1,(2**length-motif_len+1)**2,1
		write(31,trim(adjustl(formatstring))) motif_count(i),i
	end do
	close(31)
	
	open(unit=49,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motif"//trim(adjustl(temp3)))
	iprev=1
	numcount=0
	do i=2,(2**length-motif_len+1)**2,1
		if ( motif_count(i) /= motif_count(i-1) ) then
			write(49,trim(adjustl(formatstring))) motif_count(i-1), i-iprev
			iprev=i
			numcount=numcount+1
		end if
	end do
	write(49,trim(adjustl(formatstring))) motif_count(iprev),(2**length-motif_len+1)**2-iprev+1
	numcount=numcount+1
	close(49)
	
	write(*,*) 'test, motif count: ',motif_count(2)
		
	allocate( motif_pattern(numcount) )
	allocate( motif_pattern_count(numcount) )
	
	motif_pattern_count=0
	motif_pattern = ''
	
	open(unit=50,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motif"//trim(adjustl(temp3)))
	do i=1,numcount,1
		read(50,trim(adjustl(formatstring))) motif_pattern(i),motif_pattern_count(i)
	end do
	close(50)
	
	open(unit=51,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".info",position='APPEND')
	write(51,*) 'motif size: ',motif_len,' by ',motif_len
	write(51,*) 'cluster into: ',cluster_num,' class'
	write(51,*) 'the most motif is: ',motif_pattern(maxloc(motif_pattern_count)), &
		&  ' | it appears ',maxval(motif_pattern_count),' times'
	write(51,*) 'the next most motif is: ',motif_pattern(maxloc(motif_pattern_count(2:))), &
		&  ' | it appears ',maxval(motif_pattern_count(2:)),' times'
	close(51)
	
	motif_entropy=0.0
	
	do i=1,numcount,1
		temp = motif_pattern_count(i)/(2**length-motif_len+1.0)/(2**length-motif_len+1.0)
		motif_entropy = motif_entropy -  temp*log(temp)/log(2.0)
	end do
			
	write(*,*) 'the entropy sum_{i} p_i ln p_i, {i} is motif, # of {i} is cluster_num^(motif_len^2) '
	write(*,*) motif_entropy
	
	open(unit=52,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motifentropy",position='APPEND')
	write(52,*) length,cluster_num,motif_len,motif_entropy
	close(52)	
	
	write(*,*) 'output files after motif: finished'
	
	deallocate(motif_count)
	deallocate( motif_pattern )
	deallocate( motif_pattern_count )
	
	end subroutine
!==============================================================================
	subroutine shanonbox
	use output
	implicit none
	
	character(len=256) :: nosuf
	integer :: t1,t2,t3,t4
	
!	call name_trunc(datname,nosuf)
!	
!	open(unit=53,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motifentropy",position='APPEND')
!	read(53,*) t1,t2,t3,t4
!	close(53)	
	
	return
	
	
	end subroutine
!==============================================================================
	subroutine cal_motif_entropy_no_output(leng)
	use output
	implicit none


	
	!size of motif_count is (2**length-motif_len+1)**2
	
	integer :: i,j,leng,iprev,numcount
	integer :: motif(leng,leng)
	!!!real :: numb
	character(len=leng*leng) :: numb
!	character(len=leng*leng) :: numbtemp
	real :: temp
	character(len=256) :: nosuf
	!!!real,allocatable :: motif_count_array(:,:)
	character(len=leng*leng),allocatable :: motif_count(:)
	character(len=leng*leng),allocatable :: motif_pattern(:)
	integer,allocatable :: motif_pattern_count(:)
	character(len=7) :: temptemp,temp3
	character(len=128) :: formatstring
!	character(len=leng*leng) :: debugstring(100)
	
	allocate( motif_count( (2**length-motif_len+1)**2 ) )
	call name_trunc(datname,nosuf)
	
	write(temptemp,'(I7)') motif_len*motif_len
	formatstring='(A'//trim(adjustl(temptemp))//',I10)'
	write(temp3,'(I7)') cluster_num

		
!	motif_count = 0
	motif_count=''
	do i=1,2**length-motif_len+1,1
		do j=1,2**length-motif_len+1,1
			call cut_motif_minimal(i,j,leng,motif)
			call motif_to_number_c(motif,leng,numb)
			motif_count(j+(i-1)*(2**length-motif_len+1))=numb
		end do
	end do
	
	
	call QsortC(motif_count,leng*leng)
	
	open(unit=10,file=trim(adjustl(outputpath))//'/'//'temp'//trim(adjustl(temp3)))
	iprev=1
	numcount=0
	do i=2,(2**length-motif_len+1)**2,1
		if ( motif_count(i) /= motif_count(i-1) ) then
			write(10,trim(adjustl(formatstring))) motif_count(i-1), i-iprev
			iprev=i
			numcount=numcount+1
		end if
	end do
	write(10,trim(adjustl(formatstring))) motif_count(iprev),(2**length-motif_len+1)**2-iprev+1
	numcount=numcount+1
	close(10)
	
	write(*,*) 'test, motif count: ',motif_count(2)
		
	allocate( motif_pattern(numcount) )
	allocate( motif_pattern_count(numcount) )
	
	motif_pattern_count=0
	motif_pattern = ''
	
	open(unit=11,file=trim(adjustl(outputpath))//'/'//'temp'//trim(adjustl(temp3)))
	do i=1,numcount,1
		read(11,trim(adjustl(formatstring))) motif_pattern(i),motif_pattern_count(i)
	end do
	close(11)
	
!	open(unit=51,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".info",position='APPEND')
!	write(51,*) 'string length',length
!	write(51,*) 'motif size: ',motif_len,' by ',motif_len
!	write(51,*) 'cluster into: ',cluster_num,' class'
!	write(51,*) 'the most motif is: ',motif_pattern(maxloc(motif_pattern_count)), &
!		&  ' | it appears ',maxval(motif_pattern_count),' times'
!	write(51,*) 'the next most motif is: ',motif_pattern(maxloc(motif_pattern_count(2:))), &
!		&  ' | it appears ',maxval(motif_pattern_count(2:)),' times'
!	close(51)
	
	motif_entropy=0.0
	
	do i=1,numcount,1
		temp = motif_pattern_count(i)/(2**length-motif_len+1.0)/(2**length-motif_len+1.0)
		motif_entropy = motif_entropy -  temp*log(temp)/log(2.0)
	end do
			
	write(*,*) 'the entropy sum_{i} p_i ln p_i, {i} is motif, # of {i} is cluster_num^(motif_len^2) '
	write(*,*) motif_entropy
	
	open(unit=52,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motifentropy",position='APPEND')
	write(52,*) length,cluster_num,motif_len,motif_entropy
	close(52)	
	
	write(*,*) 'output files after motif: finished'
	
	deallocate(motif_count)
	deallocate( motif_pattern )
	deallocate( motif_pattern_count )
	
	end subroutine
!==============================================================================
!	subroutine cal_motif_entropy_method2(leng)
!	use output
!	implicit none
!	
!	!size of motif_count is (2**length-motif_len+1)**2
!	
!	integer :: i,j,leng,iprev,numcount
!	integer :: motif(leng,leng)
!	!!!real :: numb
!	integer :: numb
!	real :: temp
!	character(len=256) :: nosuf
!	!!!real,allocatable :: motif_count_array(:,:)
!	integer,allocatable :: motif_count_array(:,:)
!	
!	call name_trunc(datname,nosuf)
!	
!	motif_count = 0
!	
!	do i=1,2**length-motif_len+1,1
!		do j=1,2**length-motif_len+1,1
!			call cut_motif(i,j,leng,motif)
!			call motif_to_number(motif,leng,numb)
!			motif_count(j+i*(2**length-motif_len+1),1)=numb
!			motif_count(j+i*(2**length-motif_len+1),2)=j+i*(2**length-motif_len+1)
!		end do
!	end do
!	
!	call ISORT(motif_count(:,1),motif_count(:,2),(2**length-motif_len+1)**2,2)
!	!!!call SSORT(motif_count(:,1),motif_count(:,2),(2**length-motif_len+1)**2,2)
!	
!	open(unit=49,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motif")
!	iprev=1
!	numcount=0
!	do i=2,(2**length-motif_len+1)**2,1
!		if ( motif_count(i,1) /= motif_count(i-1,1) ) then
!			write(49,*) motif_count(i-1,1), i-iprev
!			iprev=i
!			numcount=numcount+1
!		end if
!	end do
!	write(49,*) motif_count(iprev,1),(2**length-motif_len+1)**2-iprev+1
!	numcount=numcount+1
!	close(49)
!		
!	allocate( motif_count_array(numcount,2) )
!	
!	motif_count_array=0
!	
!	open(unit=50,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motif")
!	do i=1,numcount,1
!		read(50,*) motif_count_array(i,:)
!	end do
!	close(50)
!	
!	open(unit=51,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".info",position='APPEND')
!	write(51,*) 'motif size: ',motif_len,' by ',motif_len
!	write(51,*) 'cluster into: ',cluster_num,' class'
!	write(51,*) 'the most motif is: ',motif_count_array(maxloc(motif_count_array(:,2)),1), &
!		&  ' | it appears ',maxval(motif_count_array(:,2)),' times'
!	write(51,*) 'the next most motif is: ',motif_count_array(maxloc(motif_count_array(2:,2)),1), &
!		&  ' | it appears ',maxval(motif_count_array(2:,2)),' times'
!	close(51)
!	
!	motif_entropy=0.0
!	
!	do i=1,numcount,1
!		temp = motif_count_array(i,2)/(2**length-motif_len+1.0)/(2**length-motif_len+1.0)
!		motif_entropy = motif_entropy -  temp*log(temp)/log(2.0)
!	end do
!			
!	write(*,*) 'the entropy sum_{i} p_i ln p_i, {i} is motif, # of {i} is cluster_num^(motif_len^2) '
!	write(*,*) motif_entropy
!	
!	open(unit=52,file=trim(adjustl(outputpath))//'/'//trim(adjustl(nosuf))//".motifentropy")
!	write(52,*) length,cluster_num,motif_len,motif_entropy
!	close(52)	
!	
!	write(*,*) 'output files after motif: finished'
!	
!	end subroutine
!==============================================================================
!	subroutine cal_motif_entropy_method1(leng)
!	implicit none
!	
!	!size of motif_count is cluster_num^(motif_len^2), cannot use for large motif_len & cluster_num
!	
!	integer :: i,j,leng
!	integer :: motif(leng,leng)
!	integer :: numb
!	real :: temp
!	
!	motif_count = 0
!	
!	do i=1,2**length-motif_len+1,1
!		do j=1,2**length-motif_len+1,1
!			call cut_motif(i,j,leng,motif)
!			call motif_to_number(motif,leng,numb)
!			motif_count(numb+1)=motif_count(numb+1)+1
!		end do
!	end do
!	
!	motif_entropy=0.0
!	
!	do i=1,cluster_num**(motif_len**2),1
!		if ( motif_count(i) /= 0 ) then
!			temp = motif_count(i)/(2**length-motif_len+1.0)/(2**length-motif_len+1.0)
!			motif_entropy = motif_entropy -  temp*log(temp)/log(2.0)
!		end if
!	end do
!			
!	write(*,*) 'the entropy sum_{i} p_i ln p_i, {i} is motif, # of {i} is cluster_num^(motif_len^2) '
!	write(*,*) motif_entropy
!	
!	end subroutine	
!==============================================================================
	subroutine motif_to_number(motif,leng,num)
	implicit none
	
	integer :: motif(leng,leng)
	integer :: i,j,leng ! j is lie/column, i is hang/row
	!!!real :: num
	integer :: num
	
	num=0
	do i=1,leng,1
		do j=1,leng,1
			num=num+motif(j,i)*(cluster_num**(j-1+leng*(i-1)))
		end do
	end do
	
	end subroutine
!===============================================================================
	subroutine motif_to_number_c(motif,leng,num)
	implicit none
	
	integer :: motif(leng,leng)
	integer :: i,j,leng ! j is lie/column, i is hang/row
	!!!real :: num
	character(leng*leng) :: num
	character :: temptemp !
	character(len=2) :: temptemp2
	
	num=''
	do i=1,leng,1
		do j=1,leng,1
			write(temptemp,'(I1)') motif(j,i)
			num=trim(adjustl(num))//temptemp
		end do
	end do
	
!	write(*,*) num
	
	end subroutine

!===============================================================================
	subroutine cut_motif(i,j,leng,motif)
	implicit none
	
	integer :: i,j,k !imagematrix_cluster(i,j)
	integer :: motif(leng,leng),leng
	
	do k=1,leng,1
		motif(:,k)=imagematrix_cluster(i:i+leng-1,j+k-1)
	end do
	
	end subroutine
!===============================================================================
	subroutine cut_motif_minimal(i,j,leng,motif)
	implicit none
	
	integer :: i,j,k !imagematrix_cluster(i,j)
	integer :: motif(leng,leng),leng
	
	do k=1,leng,1
		motif(:,k)=imagematrix(i:i+leng-1,j+k-1)
	end do
	
	end subroutine
!=====================================================================
!order
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
!=====================================================================
	end module
