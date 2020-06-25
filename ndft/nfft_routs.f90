module nfft_routs

use precision

integer, private, save	:: n		! number of input samples
integer, private, save	:: mtotal	! number of kx values
integer, private, save	:: mfirst	! first kx index
integer, private, save	:: mlast	! last  kx index
integer, private, save	:: q		! # points Gauss filter x-domain
integer, private, save	:: mc		! # Fourier comp. in comp. with FFT

integer     , private, allocatable  :: istart(:)	! start values
real(single), private, allocatable  :: gaussx(:,:)	! Gauss in x
real(single), private, allocatable  :: gaussk(:)	! Gauss in kx


interface nfft_init
!   module procedure nfft_init_si, nfft_init_dp
    module procedure nfft_init_si
end interface

interface nfft
!   module procedure nfft_si, nfft_dp
    module procedure nfft_si
end interface

contains

!***********************************************************************
!***	routine nfft for computation of ndft via the FFT
!***	in single and double precision versions
!***	and corresponding initialization routines
!***********************************************************************


subroutine nfft_init_si(x,deltakx,mtotal_in,qin,fin)
!***********************************************************
!
!  TITLE      : nfft
!
!  KEYWORDS   : non-equisistant fast Fourier transform
!
!  DESCRIPTION: 
!		initialization for the computation of
!		pout(m) = sum(n=1,N) pin(n) exp{jm*deltakx*x(n))
!		m = mfirst,...,mlast , where
!		mfirst = -mtotal/2
!		mlast  =  mtotal/2 - 1
!		and mtotal is restricted to be even
!
!		precision is single
!
!  LANGUAGE   : Fortran 90
!
!  CALLING SEQUENCE:    call nfft_init(x,deltakx,mtotal_in,qin,fin)
!
!  NOTES      :	
!
!----------------------------------------------------------------------
!
!  EXTNL. REFERENCES:	none
!
!  SYSTEM UTILITY
!  REFERENCES       :   none
!
!  NOTES            :
!
!
!  REVISION HISTORY:
!  VERSION      COMMENT                      AUTHOR          DATE
!    1.0        creation                     Adri Duijndam   25 March 1997
!
!*******************************************************************

implicit none

!--------- interface ----------------------------------

real(single)   , intent(in)   :: x(:)	   ! positions of input samples
real(single)   , intent(in)   :: deltakx   ! delta k_x 
integer, intent(in)           :: mtotal_in ! total number of kx values
integer, intent(in), optional :: qin	   ! lenght Gauss pulse (input->q)
real   , intent(in), optional :: fin	   ! factor for FFT     (input->f)

!--------- local declarations ----------------------------------

real(single), parameter :: pi=3.141592654_single
integer, parameter      :: qdefault = 10	! default for q
real   , parameter      :: fdefault = 2.	! default for f

real         :: f			 ! factor for FFT
real(single) :: xrange, xcomp(size(x))	 ! modified x values
real(single) :: dx			 ! delta x of comp. grid
real(single) :: kx(-mtotal_in/2:mtotal_in/2-1) ! k_x values
integer      :: i,j,ic,m		 ! loop indices
integer      :: iend			 ! end of a loop indicator
real(single) :: km,b			 ! max freq. + Gauss par

real(single), allocatable :: xdif(:)	 ! work array for comp of Gaussx

!-----------------------------------------------------------------------
!***	get sizes and check whether n>0 and whether M>=0 and m is even
!-----------------------------------------------------------------------

n= size(x   )

if (n < 0) then
    print *,'nfft_init_si: # input samples ndft less than 0'
    return
endif

if (mtotal_in < 0) then
    print *,'nfft_init_si: # Fourier components less than 0'
    return
endif

if (mod(mtotal_in,2) /= 0) then
    print *,'nfft_init_si: # of Fourier components is not even'
    return
else
    mtotal = mtotal_in
endif

mfirst = -mtotal/2
mlast  =  mtotal/2 - 1

!-----------------------------------------------------------------------
!***	get optional parameters if present, otherwise use default values
!-----------------------------------------------------------------------

if (present(qin)) then
    q = qin
else
    q = qdefault
endif

if (present(fin)) then
    f = fin
else
    f = fdefault
endif

!-----------------------------------------------------------------------
!***	fill kx array
!-----------------------------------------------------------------------

kx = (/ (m, m=mfirst,mlast) /) * deltakx

!-----------------------------------------------------------------------
!***	define periodicity and map x positions into this range
!-----------------------------------------------------------------------

xrange = 2.*pi/deltakx
xcomp = modulo(x,xrange)

!-----------------------------------------------------------------------
!***	define computational grid
!-----------------------------------------------------------------------

mc = nint( f*mtotal )			! total # Fourier components
dx = 2.*pi / (real(mc)*deltakx)		! delta x computational grid

!-----------------------------------------------------------------------
!***	define Gauss function
!-----------------------------------------------------------------------

km = real(mtotal/2) * deltakx			! max. Fourier variable
b  = (2*(f**2) - f) * km**2 / (pi*real(q))	! char. var. Gauss pulse

!-----------------------------------------------------------------------
!***	distr. input samples over computational grid, with Gauss pulse
!-----------------------------------------------------------------------

if ( allocated(istart) ) deallocate(istart)
if ( allocated(gaussx) ) deallocate(gaussx)
if ( allocated(gaussk) ) deallocate(gaussk)

allocate( istart(1:n) )
allocate( gaussx(q,n) )
allocate( gaussk(mfirst:mlast) )

allocate( xdif(1:q) )

do i = 1,n
    istart(i)   = floor(xcomp(i)/dx) - q/2 + 1
    iend        = istart(i) + q - 1
    xdif        = (/ (ic, ic = istart(i),iend) /) * dx - xcomp(i)
    gaussx(:,i) = dx * exp(-b*xdif**2)
enddo

deallocate( xdif )

!-----------------------------------------------------------------------
!***	array for eventual correction of Fourier domain
!-----------------------------------------------------------------------

gaussk = sqrt(b/pi) * exp( kx**2 / (4*b) )

end subroutine nfft_init_si


!***********************************************************
!***********************************************************
!	ROUTINE NFFT_SI
!***********************************************************
!***********************************************************

subroutine nfft_si(pin,pout)
!***********************************************************
!
!  TITLE      : nfft
!
!  KEYWORDS   : non-equisistant fast Fourier transform
!
!  DESCRIPTION: 
!		computes
!		pout(m) = sum(n=1,N) pin(n) exp{jm*deltakx*x(n))
!		m = mfirst,...,mlast , where
!		mfirst = -mtotal/2
!		mlast  =  mtotal/2 - 1
!		where mtotal has been defined by the initialization routine
!		and is restricted to be even
!
!		precision is single
!
!  LANGUAGE   : Fortran 90
!
!  CALLING SEQUENCE:    call nfft_si(pin,pout)
!
!  NOTES      :	 initialization by nfft_init_si
!
!----------------------------------------------------------------------
!
!  EXTNL. REFERENCES:	none
!
!  SYSTEM UTILITY
!  REFERENCES       :   none
!
!  NOTES            :
!
!
!  REVISION HISTORY:
!  VERSION      COMMENT                      AUTHOR          DATE
!    1.0        creation                     Adri Duijndam   27 March 1997
!
!*******************************************************************

use fft_gen

implicit none

!--------- interface ----------------------------------

complex(single), intent(in) :: pin(:)		! input data
complex(single), intent(out):: pout(mfirst:)	! output data

!--------- local declarations ----------------------------------

integer		:: i				! loop index
complex(single) :: pgrid(-q/2-1:mc+q/2)		! computational grid

!-----------------------------------------------------------------------
!***	get sizes and check whether n>0 and whether m>=0 and m is even
!-----------------------------------------------------------------------

if (size(pin) /= n) then
    print *,'nfft_si: size of input data incorrect!!'
    return
endif

if (size(pout) < mtotal) then
    print *,'nfft_si: size of output data incorrect!!'
    return
endif

!-----------------------------------------------------------------------
!***	distr. input samples over computational grid, with Gauss pulse
!-----------------------------------------------------------------------

pgrid = 0.			! intialize to 0.

do i = 1,n
    pgrid(istart(i):istart(i)+q-1) = pgrid(istart(i):istart(i)+q-1) + &
				     pin(i) * gaussx(:,i)
enddo

!-----------------------------------------------------------------------
!***	map ends onto actual computational grid
!-----------------------------------------------------------------------

pgrid(0:q/2-1)     = pgrid(0:q/2-1) + pgrid(mc:mc+q/2-1)  ! right -> start
pgrid(mc-q/2:mc-1) = pgrid(mc-q/2:mc-1) + pgrid(-q/2:-1)  ! left  -> end

!-----------------------------------------------------------------------
!***	(Fast) Fourier transform over comp. grid
!-----------------------------------------------------------------------

call fftb_gen( pgrid(0:mc-1) )		! FFT to go to Fourier domain

!-----------------------------------------------------------------------
!***	correct Fourier domain
!-----------------------------------------------------------------------

pout(mfirst:-1   ) = pgrid(mc-mtotal/2 :mc-1)
pout(0     :mlast) = pgrid(0:mlast)

pout = gaussk * pout

end subroutine nfft_si

end module nfft_routs

