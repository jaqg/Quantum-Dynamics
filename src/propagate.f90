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
   real(8), parameter :: au2kcalmol=627.509d0  ! = kcalmol-1/au
   real(8), parameter :: fs2au=41.341373336561d0
   real(8), parameter :: c=2.99792458D8  ! m/s
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

   integer :: i, j 
   integer :: npoints,ntime,snapshot, n_HA
   real(kind=8) :: m, omega 
   real(8) :: alpha,dt,t,dx,x0, psisq, normalization_cons, T_period, E_diff
   real(8), allocatable :: pot(:),kin(:),psisquare(:), coeff(:)
   complex(8), allocatable :: psi(:),psi0(:),exppot(:),expkin(:)
   real(kind=8) :: hermite, phi1

   open(unit=10,file='wavepacket_propagate')
   read(10,*) wavepacket_type       !Type of wavepacket: Gaussian or harmonic
   read(10,*) npoints               !Number of lattice points
   read(10,*) x0                    !Initial position
   read(10,*) alpha                 !Governs the initial width of the wave packet
   read(10,*) dt                    !Propagation time step
   read(10,*) ntime                 !Number of propagation steps
   read(10,*) snapshot              !snapshot frequency 
   close(10)

   open(unit=11,file='potential')
   read(11,*) potentialtype         !harmonic or double well potential
   read(11,*) angfreq               !Angular frequency for harmonic potential
   read(11,*) barrier               !Height of barrier in double well potential (in kcal/mol)
   close(11)

   dt=dt*fs2au                        !convert femtoseconds to atomic units
   angfreq=angfreq/fs2au              !convert femtoseconds to atomic units

   allocate(psi(npoints),psi0(npoints))
   allocate(pot(npoints),exppot(npoints))
   allocate(kin(npoints),expkin(npoints))
   allocate(psisquare(npoints))

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

   ! Questions
   write(6,'(a,f10.4)') 'Distance between adjancent grid points: ', dx
   write(6,'(a,f10.4)') 'Initial position of the centre of the WP: ', x0
   write(6,'(a,f10.4,a)') 'Total propagation time: ', &
   & (dble(ntime) * dt) / fs2au, ' fs'
   write(6,'(a,f10.4,a)') 'Interval of the snapshots: ', &
   & snapshot * dt / fs2au, ' fs'
   T_period = (2.d0 * pi / angfreq)
   write(6,'(a,f10.4,a)') 'Period T of the prob density: ', &
   & T_period / fs2au, ' fs'
   write(6,'(a,f10.4,a)') 'Energy of the barrier in the double-well: ', &
   & barrier, ' kcal/mol'

   psi=psi0                                            !Set the wavepacket psi at t=0 equal psi0
   ! Loop over time
   do i=0,ntime                                        !Start propagation
      t=i*dt
      if (i>0) then
         psi=psi*exppot                                !Multiply psi with exp(-i*dt*potential)
         call fourier(1,npoints,psi)                   !Forward FFT to momentum space
         psi=psi*expkin                                !Multiply psi with the exp(-i*dt*kinetic operator)
         call fourier(-1,npoints,psi)                  !Backward FFT to position space
      endif
      if (mod(i,snapshot)==0) then                     !Take a snapshot if the remainder of i/snapshot equals 0
         call initgraph(i/snapshot,t)                  !Initialize graph
         psisquare=(abs(psi))**2
         call graphpot(dx,npoints)                     !Plot the potential
         call graphpsi(dx,npoints,psisquare)           !Plot |psi|^2
      endif
   end do                                              !End propagation

   ! Calculate the energy difference between two stationary states
   E_diff = 1.d0/((T_period/fs2au) * 1.D-15 * c) * 1.D-2  ! cm-1
   write(6,'(a,f10.4,a)') 'Energy difference between two stationary states: ',&
   & E_diff, ' cm-1'

   deallocate(psi,psi0)
   deallocate(pot,exppot)
   deallocate(kin,expkin)
   deallocate(psisquare)

end program propagate
!---------------------!

!--------------------------------------------!
subroutine initpsi_gaussian(npoints,dx,alpha,x0,psi0)
!--------------------------------------------!
   implicit none

   integer :: i,j,npoints
   real(8) :: alpha,x,x0,dx,psisq,normalization_cons
   complex(8) :: psi0(npoints)

   do i=-npoints/2+1,npoints/2
      x=dble(i)*dx
      if (i>0) then
         j=i
      else     
         j=i+npoints
      endif
      psi0(j)=exp(-alpha*(x-x0)**2)
   end do

   ! Normalize the wavefunction
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

end subroutine initpsi_gaussian
!----------------------!

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
      ! Note: (0,1) is the imaginary unit, i
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
