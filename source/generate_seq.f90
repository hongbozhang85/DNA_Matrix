! prob(4) publicy
module random_version
use variable
implicit none

real :: prob(4)
integer :: countnum(4)

contains
!===========================================================================
subroutine set_seq
implicit none

!type(DNA_Data) :: DNA

 call generate_random_seq(len_seq,seq_origin,seq_random)
 call generate_mirror_seq(len_seq,seq_origin,seq_mirror)

end subroutine
!===========================================================================
subroutine generate_seq(n,option,seq)
implicit none

!random generate a DNA sequence of length n. with prob of atcg (0.25,0.25,0.25,0.25)
!option=T: atcgatcg...
!option=F: 12341234....(with 1=a, 2=t, 3=c, 4=g)

integer,intent(in) :: n
logical,intent(in) :: option
 character :: seq(n)
real :: x
integer :: nucleuo_acid,i

write(*,*) 'generate a sequence. for debug'
 call random_seed()

do i=1,n,1
	call random_number(x)
	nucleuo_acid=int(4*x+1)
	if ( option ) then
		if ( nucleuo_acid==1 ) then
			seq(i) = 'A'
		elseif (nucleuo_acid==2) then
			seq(i) = 'T'
		elseif (nucleuo_acid==3) then
			seq(i) = 'C'
		elseif (nucleuo_acid==4) then
			seq(i) = 'G'
		endif
!		write(*,*) x,nucleuo_acid,seq(i)
	else
		write(seq(i),'(I1)') nucleuo_acid
	endif
end do

!write(*,*) seq

end subroutine generate_seq

!============================================================================

subroutine generate_random_seq(lenseq,seq_in,seq_out)
implicit none

! given seq_in, generate seq_out which has the same probability of a/t/c/g as seq_in
! that is to say there are prob(1) A, prob(2) T, prob(3) C, prob(4) G, in seq_in
! generate seq_out who has the same prob(:) as seq_in

integer,intent(in) :: lenseq
 character,intent(in) :: seq_in(lenseq)
 character :: seq_out(lenseq)
integer :: i     !,countnum(4)
!real :: prob(4)


do i=1,lenseq,1	
	if ( seq_in(i) == 'A' ) then
		countnum(1)=countnum(1)+1
	elseif ( seq_in(i) == 'T' ) then
		countnum(2)=countnum(2)+1
	elseif ( seq_in(i) == 'C' ) then
		countnum(3)=countnum(3)+1
	elseif ( seq_in(i) == 'G' ) then
		countnum(4)=countnum(4)+1
	endif
enddo

prob = (countnum+0.0)/lenseq

write(*,*) 'numbers of atcg in this genome: ',countnum
write(*,*) 'probability of atcg in this genome: ',prob

 call generate_random_seq_given_prob(lenseq,prob,seq_out)

end subroutine

!===========================================================================

subroutine generate_random_seq_given_prob(lenseq,probability,seq_out)
implicit none

! generate seq_out which satisfy distribution of a/t/c/g as prob(4)
! subroutine generate_seq can be viewed as prob(:)=(0.25,0.25,0.25,0.25)

integer,intent(in) :: lenseq
real,intent(in) :: probability(4) ! probability of atcg:  prob(a),prob(t),prob(c),prob(g)
 character :: seq_out(lenseq)
real :: x
integer :: i

 call random_seed()

do i=1,lenseq,1
	call random_number(x)
	if ( x <= probability(1) ) then
		seq_out(i) = 'A'
	elseif ( x <= (probability(1)+probability(2)) ) then
		seq_out(i) = 'T'
	elseif ( x <= (probability(1)+probability(2)+probability(3)) ) then
		seq_out(i) = 'C'
	elseif ( x <= 1.0 ) then
		seq_out(i) = 'G'
	endif
end do
	
	write(*,*) 'complete random version finished'

end subroutine

!========================================================================

subroutine generate_mirror_seq(lenseq,seq_in,seq_out)
implicit none

integer,intent(in) :: lenseq
 character,intent(in) :: seq_in(lenseq)
 character :: seq_out(lenseq)
integer :: i

do i=1,lenseq,1
	seq_out(i) = seq_in(lenseq-i+1)
enddo
	
	write(*,*) 'mirror version finished'

end subroutine

!=======================================================================
end module
