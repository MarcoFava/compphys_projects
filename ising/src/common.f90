module common

  implicit none
  public :: initial,metropolis,DeltaE
  public :: data,output
  double precision, public :: T,E,M
  integer, public, dimension(:,:), allocatable :: spin
  double precision, public, dimension(-8:8) :: w
  integer, public, dimension(8) :: seed
  integer, public :: N,L,nmcs,nequil
  integer, public :: accept
contains

  subroutine initial(L_in,T_in,nequil_in,nmcs_in,seed_in,cum)
    integer, intent(in) :: L_in, nequil_in, nmcs_in
    double precision, intent(in) :: T_in
    integer, dimension(8), intent(in):: seed_in

    double precision, dimension(5), intent (out) :: cum
    integer :: x,y,up,right,sums,i,dE
    real :: rnd
!     print *, "linear dimension of lattice L ="
!     read *, L
!     print *, "reduced temperature T ="
!     read *, T
!     print *, "# MC steps per spin for equilibrium ="
!     read *, nequil
!     print *, "# MC steps per spin for averages ="
!     read *, nmcs
!     print *, "seed (1:8) ="
!     read *, seed

    L = L_in
    T = T_in
    nequil = nequil_in
    nmcs = nmcs_in
    seed = seed_in

    allocate(spin(L,L))     ! PROBLEMA????
    N = L*L


    call random_seed(put=seed)
    M = 0.0d0
    !  random initial configuration
    !  compute initial magnetization
    do y = 1,L
       do x = 1,L
          call random_number(rnd)
          if (rnd < 0.5) then
             spin(x,y) = 1
          else
             spin(x,y) = -1
          end if
          M = M + spin(x,y)
       end do
    end do
    !  compute initial energy
    E = 0.0d0
    do y = 1,L
       !  periodic boundary conditions
       if (y == L) then
          up = 1
       else
          up = y + 1
       end if
       do x = 1,L
          if (x == L) then
             right = 1
          else
             right = x + 1
          end if
          sums = spin(x,up) + spin(right,y)
! calculate the initial energy summing all over pairs
! (for a given spin, consider only the up NN and the right NN
! - NOT the down and the left NN - : each interaction is counted once
          E = E - spin(x,y)*sums
       end do
    end do
!
! calculate the transition probability according
! to the Boltzmann distribution (exp(-deltaE/KT).
! Choosing the interaction parameter  J=1, ***IN CASE OF P.B.C.***
! possible energy variations per spin flip are -8,-4,0,+4,+8:
    do dE = -8,8,4
       w(dE) = exp(-dE/T)
    end do
    accept = 0
    cum = 0.0d0
  end subroutine initial

  subroutine metropolis()
    !  one Monte Carlo step per spin
    integer :: ispin,x,y,dE
    real :: rnd
    do ispin = 1,N
       !     random x and y coordinates for trial spin
       call random_number(rnd)
       x = int(L*rnd) + 1
       call random_number(rnd)
       y = int(L*rnd) + 1
       dE = DeltaE(x,y)
       call random_number(rnd)
       if (rnd <= w(dE)) then
          spin(x,y) = -spin(x,y)
          accept = accept + 1
          M = M + 2*spin(x,y)  ! factor 2 is to account for the variation:
          E = E + dE           ! (-(-)+(+))
       end if
    end do
  end subroutine metropolis

  function DeltaE(x,y) result (DeltaE_result)
    !  periodic boundary conditions
    integer, intent (in) :: x,y
    integer :: DeltaE_result
    integer :: left
    integer :: right
    integer :: up
    integer :: down
    if (x == 1) then
       left = spin(L,y)
       right = spin(2,y)
    else if (x == L) then
       left = spin(L-1,y)
       right = spin(1,y)
    else
       left = spin(x-1,y)
       right = spin(x+1,y)
    end if
    if (y == 1) then
       up = spin(x,2)
       down = spin(x,L)
    else if (y == L) then
       up = spin(x,1)
       down = spin(x,L-1)
    else
       up = spin(x,y+1)
       down = spin(x,y-1)
    end if
    DeltaE_result = 2*spin(x,y)*(left + right + up + down)
! also here the factor 2 is to account for the variation
  end function DeltaE

  subroutine data(cum)
    !  accumulate data after every Monte Carlo step per spin
    double precision, dimension(5), intent (inout) :: cum
    cum(1) = cum(1) + E
    cum(2) = cum(2) + E*E
    cum(3) = cum(3) + M
    cum(4) = cum(4) + M*M
    cum(5) = cum(5) + abs(M)
  end subroutine data

  subroutine print_output(cum)
    double precision, dimension(5), intent (inout) :: cum
    double precision :: eave,e2ave,mave,m2ave,abs_mave,cv,chi
    real :: acceptance_prob
    acceptance_prob = real(accept)/real(N)/real(nmcs+nequil)
    eave     = cum(1)/real(N)/real(nmcs)    ! to avoid interger overflow
    e2ave    = cum(2)/real(N*N)/real(nmcs)
    mave     = cum(3)/real(N)/real(nmcs)
    m2ave    = cum(4)/real(N*N)/real(nmcs)
    abs_mave = cum(5)/real(N)/real(nmcs)
    cv       = (e2ave - eave*eave)/T/T*real(N)
    chi      = (m2ave - mave*mave)/T*real(N)
    print *, "temperature                =", T
    print *, "acceptance probability     =", acceptance_prob
    print *, "mean energy per spin                =", eave
    print *, "mean squared energy per spin        =", e2ave
    print *, "mean magnetization  per spin        =", mave
    print *, "mean squared magnetization per spin =", m2ave
    print *, "mean |magnetization| per spin       =", abs_mave
    print *, "specific heat per spin              =", cv
    print *, "susceptibility per spin             =", chi
  end subroutine print_output

  subroutine output(cum,eave,e2ave,mave,m2ave,abs_mave,cv,chi,acceptance_prob)
    double precision, dimension(5), intent (in) :: cum
    double precision, intent(out) :: eave,e2ave,mave,m2ave,abs_mave,cv,chi
    real, intent(out) :: acceptance_prob
    acceptance_prob = real(accept)/real(N)/real(nmcs+nequil)
    eave            = cum(1)/real(N)/real(nmcs)    ! to avoid interger overflow
    e2ave           = cum(2)/real(N*N)/real(nmcs)
    mave            = cum(3)/real(N)/real(nmcs)
    m2ave           = cum(4)/real(N*N)/real(nmcs)
    abs_mave        = cum(5)/real(N)/real(nmcs)
    cv              = (e2ave - eave*eave)/T/T*real(N)
    chi             = (m2ave - mave*mave)/T*real(N)
  end subroutine output

  subroutine finish()
!     print *, "acceptance probability =", real(accept)/real(N)/real(nmcs+nequil)
    DEALLOCATE(spin)
  end subroutine finish

end module common
