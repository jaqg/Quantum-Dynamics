!-----------------!
module IO
!-----------------!
   ! Module for Input/Output procedures
   contains

   ! -------------------------------------------!
   function file_exists(filename) result(res)
   ! -------------------------------------------!
      !
      ! Function to check if a file exists
      !
      implicit none
      character(len=*),intent(in) :: filename
      logical :: res

      ! Check if the file exists
      inquire(file=trim(filename), exist=res)

   end function
   ! -------------------------------------------!

   ! -------------------------------------------!
   subroutine read_coefficients(filename, n_HA, coeff)
   ! -------------------------------------------!
      !
      ! Subroutine to read the coefficients, c_n, from a file
      !
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: n_HA 
      real(kind=8), dimension(:), allocatable, intent(out) :: coeff
      integer :: i, n, uf, stat
      real(kind=8) :: cn 

      ! Chech if the file exists
      if(.not. file_exists(filename)) then
         write(*,*) 'File ', trim(filename), ' does not exist'
         stop
      endif
      
      ! Open input file
      open(newunit=uf, file=trim(filename), status='old', action='read')

      ! Read number of eigenfunctions to include in the WF
      read(uf,*)
      read(uf,*) n_HA

      ! Get dimensions of coeff
      read(uf,*)
      n = -1
      stat = 0
      do while(stat == 0)
         n = n + 1
         read(uf,*,iostat=stat) cn
      enddo
                          
      ! Allocate coeff
      allocate(coeff(n))

      ! Rewind input file unitfile
      rewind(uf)

      ! Skip rows
      do i=1, 3
         read(uf,*)
      end do

      ! Read the coefficients into the array coeff
      do i=1, n
         read(uf,*) coeff(i)
      end do
      
      !
      return
   end subroutine read_coefficients 
   ! -------------------------------------------!
end module IO
!-----------------!

!-----------------!
module parameters
!-----------------!

   real(8), parameter :: length=5.12d0   ! length of the box (in Bohr)
   real(8), parameter :: mass = 1822.88839d0 != 1 atomic mass unit (=1 g/mol)
   real(8), parameter :: pi=3.141592653589793d0
   real(8), parameter :: au2kcalmol=627.509d0, fs2au=41.341373336561d0, au2invcm=219474.63d0
   real(8), parameter :: hbar=1.d0  ! atomic units
   real(8) :: angfreq, barrier
   character(len=80) :: wavepacket_type, potentialtype 

end module parameters
!---------------------!

!-----------------!
module modulo_OA1D
!-----------------!
   use parameters
   !
   ! Modulo para el calculo de las funciones propias del oscilador armonico
   ! unidimensional
   !
   implicit none
   contains
   ! function E_OA1D(n, hbar, k, m) result(energia)
   !    !
   !    ! Funcion para el calculo de las energias propias del oscilador
   !    ! armonico unidimensional
   !    !
   !    implicit none
   !    integer, intent(in) :: n
   !    real*8, intent(in) :: hbar, k, m
   !    real*8 :: omega, energia
   !
   !    omega = dsqrt(k/m)
   !    energia = (n*1.d0 + 0.5d0)*hbar*omega
   !
   !    return
   ! end function E_OA1D
   !--------------------------------------------!
   recursive function hermp(n, x) result(res_hp)
   !--------------------------------------------!
      !
      ! Function to compute the Hermite polynomials using the recurrence relation
      ! H_n(x) = 2xH_n-1(x) - 2(n-1)H_n-2(x)
      ! taking into account that
      ! H_0(x) = 1
      ! H_1(x) = 2x
      !
      ! n: n-th Hermite polynomial
      ! x: independent variable
      !
      ! n -> integer
      ! x -> real
      !
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: x
      real(kind=8) :: hp1, hp2, res_hp

      if (n==0) then
         res_hp = 1.d0
      elseif (n==1) then
         res_hp = 2.d0 * x
      else
         hp1 = 2.d0 * x * hermp(n-1, x)
         hp2 = 2.d0 * dble(n-1) * hermp(n-2, x)
         res_hp = hp1 - hp2
      end if

      return
   end function hermp
   !--------------------------------------------!
    
   !--------------------------------------!
   function phi(n, x, m, omega) result(res_phi)
   !--------------------------------------!
      !
      ! Function to compute the eigenstates, phi_n, of the 1D harmonic oscillator
      !
      ! n: n-th eigenstate
      ! x: position (independent variable)
      !
      ! n -> integer
      ! x -> real
      !
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: x, m, omega
      real(kind=8) :: alpha, Nv, res_phi

      ! For the 1D harmonic oscillator
      alpha = m * omega/hbar

      ! Note: factorial(n) = Gamma(n+1), where n is an integer number
      Nv = (alpha/pi)**0.25d0 * 1.d0/dsqrt(2.**n * Gamma(dble(n+1)))

      ! Eigenfunction
      res_phi = Nv * exp(-alpha * x**2/2.d0) * hermp(n, dsqrt(alpha) * x)

      return
   end function phi
   !--------------------------------------!

   !--------------------------------------------------------------------!
   subroutine initpsi_harmonic(n_HA, npoints, dx, coeff, m, omega, psi0)
   !--------------------------------------------------------------------!
      implicit none
      integer, intent(in) :: n_HA, npoints 
      real(8), intent(in) :: dx, coeff(:), m, omega
      complex(8), dimension(:) :: psi0
      integer :: i, j, k, ierr 
      real(kind=8) :: x, suma, psisq, normalization_cons

      ! Check that the specified number of eigenfunctions is less or equal
      ! than the total number of coefficients
      if (n_HA > size(coeff)) then
         write(6,'(a)') 'ERROR initpsi_harmonic: n_HA > size(coeff)'
         stop
      end if

      ! Main loop for the wavefunction over position
      do i=-npoints/2+1,npoints/2
         ! Compute the position
         x=dble(i)*dx

         ! Store first the possitive positions
         if (i>0) then
            j=i
         else     
            j=i+npoints
         endif

         ! Compute the wavefunction as a LC of HA eigenfunctions
         suma = 0.d0
         do k=1, n_HA
            suma = suma + coeff(k) * phi(k, x, m, omega)
         end do
         psi0(j) = suma
      end do

      ! Normalize the wavefunction:
      ! Compute the sum |Psi|^2
      psisq = 0.d0
      do j = 1, npoints
         psisq = psisq + abs(psi0(j))**2
      end do
      psisq = psisq * dx

      ! Compute the normalization constant
      normalization_cons = dsqrt(1.d0/(psisq))

      ! Normalize the wavefunction
      psi0 = normalization_cons * psi0

   end subroutine initpsi_harmonic
   !--------------------------------------------------------------!

end module
!-----------------!

!-----------------!
program propagate
!-----------------!
   use IO
   use parameters
   use modulo_OA1D
   implicit none

   integer :: npoints,ntime,snapshot,i, n_HA
   real(8) :: alpha,dt,dx,x0, energy_t0, m, omega
   real(8), dimension(:), allocatable :: pot, kin, t, norm, energy, &
   & survival_probability, eigenval_spectrum, coeff(:)
   complex(8), dimension(:), allocatable :: psi, psi0, exppot, expkin, correlation

   open(unit=10,file='wavepacket_propagate')
   read(10,*) wavepacket_type       !Type of wavepacket: Gaussian or harmonic
   close(10)

   open(unit=10,file='wavepacket_analysis')
   read(10,*) npoints               !Number of lattice points
   read(10,*) x0                    !Initial position
   read(10,*) alpha                 !Governs the initial width of the wave packet
   read(10,*) dt                    !Propagation time step
   read(10,*) ntime                 !Number of propagation steps
   read(10,*) snapshot              !snapshot frequency 
   close(10)

   open(unit=11,file='potential')
   read(11,*) potentialtype         !harmonic or double well potential
   read(11,*) angfreq               !Angular frequency for harmonic potential (in fs^(-1))
   read(11,*) barrier               !Height of barrier in double well potential (in kcal/mol)
   close(11)

   dt=dt*fs2au                        !convert femtoseconds to atomic units
   angfreq=angfreq/fs2au              !convert femtoseconds to atomic units

   allocate(psi(npoints),psi0(npoints))
   allocate(pot(npoints),exppot(npoints))
   allocate(kin(npoints),expkin(npoints))
   allocate(t(ntime),norm(ntime),energy(ntime),survival_probability(ntime))
   allocate(correlation(ntime), eigenval_spectrum(ntime))

   dx=length/dble(npoints)

      ! Compute initial wavefunction, psi0
   m = 1.d0
   omega = 1.d0
   if (wavepacket_type=='gaussian' .or. wavepacket_type=='Gaussian') then
      call initpsi_gaussian(npoints,dx,alpha,x0,psi0)
   elseif (wavepacket_type=='harmonic' .or. wavepacket_type=='Harmonic') then
      call read_coefficients('LC-coefficients', n_HA, coeff)
      call initpsi_harmonic(n_HA, npoints, dx, coeff, m, omega, psi0)
   endif

   !Initialize the FFT
   call fourier(0,npoints,psi0)

   !Calculate the kinetic and potential operators
   call operators(npoints,dx,dt,pot,kin,exppot,expkin)

   open(unit=12,file='norm')
   open(unit=13,file='energy')
   open(unit=14,file='survival')

   write(12,'(a,12x,a)') 'Time', 'Norm'
   write(13,'(a,12x,a,2x,a,2x,a)') 'Time', 'Energy (Ha)', 'Energy (kcal/mol)',&
   & 'E(t) - E(t=0) (kcal/mol)'
   ! write(14,'(a,12x,a,8x,a)') 'Time', 'C(t)', 'S(t)'
   write(14,'(a,12x,a)') 'Time', 'S(t)'

   !Set the wavepacket psi at t=0 equal psi0
   psi=psi0
   ! Compute the energy at t=0
   call calcenergy(npoints, psi, dx, pot, kin, energy_t0)
   ! Start the propagation
   do i=1,ntime
      t(i)=dble(i)*dt
      psi=psi*exppot                                   !Multiply psi with exp(-i*dt*potential)
      call fourier(1,npoints,psi)                      !Forward FFT to momentum space
      psi=psi*expkin                                   !Multiply psi with the exp(-i*dt*kinetic operator)
      call fourier(-1,npoints,psi)                     !Backward FFT to position space
      call calcnorm(npoints, psi, dx, t(i), norm(i))
      write(12,'(f10.2,2x,f10.4)') t(i), norm(i)
      call calcenergy(npoints, psi, dx, pot, kin, energy(i))
      write(13,'(f10.2,2x,f10.4,3x,f10.4,10x,f10.4)') t(i), energy(i), &
      & energy(i) * au2kcalmol, (energy(i) - energy_t0) * au2kcalmol
      call autocorrelation(npoints, psi0, psi, dx, correlation(i), &
      & survival_probability(i))
      write(14,'(f10.2,2x,f10.4)') t(i), survival_probability(i)
   end do                                            !End propagation
   close(12)
   close(13)
   close(14)

   call spectrum(ntime, correlation, energy, t, dt, eigenval_spectrum)

   open(unit=15,file='eigenvalue-spectra')
   write(15,'(a,4x,a,4x,a,4x,a)') 'Energy (kcal/mol)', 'Energy (cm-1)', &
   'C(E) (kcal/mol)', 'C(E) (cm-1)'
   do i=1,ntime
      write(15,'(f10.4,10x,f10.4,4x,f10.4,13x,f15.4)') energy(i), &
      & energy(i) * au2invcm, eigenval_spectrum(i), &
      & survival_probability(i) * au2invcm
   end do
   close(15)

   deallocate(psi,psi0)
   deallocate(pot,exppot)
   deallocate(kin,expkin)

end program propagate
!---------------------!

!---------------------------------------------------!
subroutine initpsi_gaussian(npoints,dx,alpha,x0,psi0)
!---------------------------------------------------!
   implicit none

   integer :: i,j,npoints
   real(8) :: alpha,x,x0,dx,norm
   complex(8) :: psi0(npoints)

   norm=0.d0
   do i=-npoints/2+1,npoints/2
      x=dble(i)*dx
      if (i>0) then
         j=i
      else     
         j=i+npoints
      endif
      psi0(j)=exp(-alpha*(x-x0)**2)
      norm=norm+abs(psi0(j))**2*dx
   end do
   norm=1.d0/dsqrt(norm)
   do i=1,npoints
      psi0(i)=psi0(i)*norm
   end do

end subroutine initpsi_gaussian
!-------------------------------!

!---------------------------------------------------------!
subroutine operators(npoints,dx,dt,pot,kin,exppot,expkin)
!---------------------------------------------------------!
   use parameters
   implicit none

   integer :: i,j,npoints
   real(8) :: x,p,b,dt,dx,dp,pot(npoints),kin(npoints)
   complex(8) :: exppot(npoints),expkin(npoints)

   dp=2.d0*pi/length
   do i=-npoints/2+1,npoints/2
      x=dble(i)*dx
      p=dble(i-1)*dp
      if (i>0) then
         j=i
      else
         j=i+npoints
      endif
      if (potentialtype=='harmonic') then
         pot(j)=0.5d0*mass*angfreq**2*x**2
      elseif (potentialtype=='doublewell') then
         pot(j)=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)/au2kcalmol
      endif
      kin(j)=0.5d0*p**2/mass
      exppot(j)=exp(-dt*(0,1)*pot(j))
      expkin(j)=exp(-dt*(0,1)*kin(j))
   end do

end subroutine operators
!------------------------!

!-----------------------------------!
subroutine fourier(dir,npoints,psi)
!-----------------------------------!
   implicit none

   integer :: i,npoints,dir
   real(8) :: nr
   complex(8) :: psi(npoints)
   real(8), allocatable, save :: wsave(:)       

   if (dir==1) then
      call dcfftf(npoints,psi,wsave)
      nr=1.d0/sqrt(dble(npoints))
      do i=1,npoints
         psi(i)=psi(i)*nr
      end do
   elseif (dir==-1) then
      call dcfftb(npoints,psi,wsave)
      nr=1.d0/sqrt(dble(npoints))
      do i=1,npoints
         psi(i)=psi(i)*nr
      end do
   elseif (dir==0) then
      if (allocated(wsave)) deallocate(wsave)
      allocate(wsave(4*npoints+20))
      call dcffti(npoints,wsave)
   endif

end subroutine fourier
!----------------------!

!---------------------!
subroutine calcnorm(n, psi, dx, t, norm)
!---------------------!
   use parameters
   implicit none
   integer, intent(in) :: n 
   complex(kind=8), dimension(n), intent(in) :: psi 
   real(kind=8), intent(in) :: dx, t 
   real(kind=8), intent(out) :: norm 
   integer :: i 

   norm = 0.d0
   do i = 1, n
      norm = norm + abs(psi(i))**2
   end do
   norm = norm * dx

end subroutine calcnorm
!-----------------------!

!-----------------------!
subroutine calcenergy(n, psi, dx, pot, kin, E)
!-----------------------!
   use parameters
   implicit none
   integer, intent(in) :: n 
   complex(kind=8), dimension(n), intent(in) :: psi
   real(kind=8), dimension(n), intent(in) :: pot, kin
   real(kind=8), intent(in) :: dx
   real(kind=8), intent(out) :: E 
   integer :: i 

   E = 0.d0

   ! First compute the kinetic energy contribution in momentum space
   call fourier(1,n,psi)  !Forward FFT to momentum space
   do i = 1, n
      E = E + conjg(psi(i)) * kin(i) * psi(i)
   end do
   ! Backward psi to position space to compute the potential contribution
   call fourier(-1,n,psi)  !Backward FFT to position space
   do i = 1, n
      E = E + conjg(psi(i)) * pot(i) * psi(i)
   end do
   E = E * dx

end subroutine calcenergy
!-------------------------!

!----------------------------!
subroutine autocorrelation(n, psi0, psi, dx, C, S)
!----------------------------!
   use parameters
   implicit none
   integer, intent(in) :: n 
   complex(kind=8), dimension(n), intent(in) :: psi0, psi
   real(kind=8), intent(in) :: dx
   complex(kind=8), intent(out) :: C
   real(kind=8), intent(out) :: S
   integer :: i 
   
   ! Compute the autocorrelation function
   C = 0.d0
   do i = 1, n
      C = C + conjg(psi0(i)) * psi(i)
   end do
   C = C * dx

   ! Compute the survival probability
   S = abs(C)**2

end subroutine autocorrelation
!------------------------------!

!---------------------!
subroutine spectrum(n, C, E, t, dt, spect)
!---------------------!
   use parameters
   implicit none
   integer, intent(in) :: n 
   complex(kind=8), dimension(n), intent(in) :: C
   real(kind=8), dimension(n), intent(in) :: E, t
   real(kind=8), intent(in) :: dt
   real(kind=8), dimension(n), intent(out) :: spect 
   integer :: i 

   spect = 0.d0
   do i = 1, n
      ! Note: (0,1) is the imaginary unit, i
      spect(i) = C(i) * exp((0,1) * E(i) * t(i))
   end do
   spect = 1.d0/(2.d0 * pi) * spect * dt

end subroutine spectrum
!-----------------------!
