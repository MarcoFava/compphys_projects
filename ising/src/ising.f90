!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c ising.f90
!c
!c Metropolis algorithm to calculate <E>, <M>, in the canonical ensemble
!c (fix T,N,V) with a 2D Ising model
!c
!c Here: K_B = 1
!c       J   = 1
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


subroutine ising(L_in,T_in,nequil_in,nmcs_in,seed_in, eave,e2ave,mave,m2ave,abs_mave,cv,chi,acceptance_prob)
    use common

!  metropolis algorithm for the ising model on a square lattice
    integer, intent(in) :: L_in, nequil_in, nmcs_in
    double precision, intent(in) :: T_in
    integer, dimension(8), intent(in):: seed_in
    double precision, intent(out) :: eave,e2ave,mave,m2ave,abs_mave,cv,chi
    real, intent(out) :: acceptance_prob

  double precision, dimension(5) :: cum
  integer :: imcs,ispin,jspin
  call initial(L_in,T_in,nequil_in,nmcs_in,seed_in,cum)
  !  equilibrate system
  do imcs = 1,nequil
     call metropolis()
  end do
  !  accumulate data while updating spins
  do imcs = 1,nmcs
     call metropolis()
     call data(cum)
  end do
!   call print_output(cum)
  call output(cum,eave,e2ave,mave,m2ave,abs_mave,cv,chi,acceptance_prob)

! write the coordinates of spins up and down on files for plotting
!   open(unit=8,file='ising-up.dat',status='replace')
!   open(unit=9,file='ising-down.dat',status='replace')
!   do jspin = 1,L
!      do ispin = 1,L
!         if(spin(ispin,jspin)==1)write(8,*)ispin,jspin
!         if(spin(ispin,jspin)==-1)write(9,*)ispin,jspin
!      end do
!   end do
!   close(8)
!   close(9)
  
  call finish()
end subroutine ising
