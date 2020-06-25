module fft_gen

use precision

interface fftb_gen
    module procedure fftb_1d_single, fftb_1d_double
    module procedure fftb_2d_single, fftb_2d_double
endinterface

contains

subroutine fftb_1d_single(data)
!***********************************************************
!
!  TITLE      : fftb_gen
!
!  KEYWORDS   : FFT
!
!  DESCRIPTION: general input backward fft
!		at this point in time, only complex data implemented
!		numerical precision: kind = single
!
!  LANGUAGE   : Fortran90
!
!  CALLING SEQUENCE:    call fftb_gen(data)
!
!  ARGUMENTS: see interface
!
!  NOTES      :	
!
!
!----------------------------------------------------------------------
!
!  EXTNL. REFERENCES:	
!
!  SYSTEM UTILITY
!  REFERENCES       :   
!
!  NOTES            :
!
!
!  REVISION HISTORY:
!  VERSION      COMMENT                         AUTHOR          DATE
!    1.0        creation                        Adri Duijndam   feb 1997
!
!*******************************************************************

!------- interface ------------------------------

complex(single), intent(inout) :: data(:)

!------- local declarations ------------------------------

logical					:: initialized = .false.
real(single) ,	allocatable,save	:: work(:)
integer					:: n
integer,	save			:: nprevious

!-----------------------------------------------------------------------
!***	get size of data set
!-----------------------------------------------------------------------

n = size(data)

!-----------------------------------------------------------------------
!***	initialize work array if necessary
!-----------------------------------------------------------------------

if (.not. initialized) then
    allocate(work(1:4*n+15))
    call cffti(n,work)
    initialized = .true.
    nprevious = n
elseif (n /= nprevious) then
    deallocate(work)
    allocate(work(1:4*n+15))
    call cffti(n,work)
    nprevious = n
endif

!-----------------------------------------------------------------------
!***	the actual fft
!-----------------------------------------------------------------------

call cfftb(n,data,work)

end subroutine fftb_1d_single

!=======================================================================
!=======================================================================

subroutine fftb_1d_double(data)
!***********************************************************
!
!  TITLE      : fftb_gen
!
!  KEYWORDS   : FFT
!
!  DESCRIPTION: general input backward fft
!		at this point in time, only complex data implemented
!		numerical precision: kind = double
!
!  LANGUAGE   : Fortran90
!
!  CALLING SEQUENCE:    call fftb_gen(data)
!
!  ARGUMENTS: see interface
!
!  NOTES      :	
!
!
!----------------------------------------------------------------------
!
!  EXTNL. REFERENCES:	
!
!  SYSTEM UTILITY
!  REFERENCES       :   
!
!  NOTES            :
!
!
!  REVISION HISTORY:
!  VERSION      COMMENT                         AUTHOR          DATE
!    1.0        creation                        Adri Duijndam   feb 1997
!
!*******************************************************************

!------- interface ------------------------------

complex(double), intent(inout) :: data(:)

!------- local declarations ------------------------------

logical					:: initialized = .false.
real(double) ,	allocatable,save	:: work(:)
integer					:: n
integer,	save			:: nprevious

!-----------------------------------------------------------------------
!***	get size of data set
!-----------------------------------------------------------------------

n = size(data)

!-----------------------------------------------------------------------
!***	initialize work array if necessary
!-----------------------------------------------------------------------

if (.not. initialized) then
    allocate(work(1:4*n+15))
    call zffti(n,work)
    initialized = .true.
    nprevious = n
elseif (n /= nprevious) then
    deallocate(work)
    allocate(work(1:4*n+15))
    call zffti(n,work)
    nprevious = n
endif

!-----------------------------------------------------------------------
!***	the actual fft
!-----------------------------------------------------------------------

call zfftb(n,data,work)

end subroutine fftb_1d_double

!=======================================================================
!=======================================================================

subroutine fftb_2d_single(data)
!***********************************************************
!
!  TITLE      : fftb_gen
!
!  KEYWORDS   : FFT
!
!  DESCRIPTION: general input backward 2D fft
!		at this point in time, only complex data implemented
!		numerical precision: kind = single
!
!  LANGUAGE   : Fortran90
!
!  CALLING SEQUENCE:    call fftb_gen(data)
!
!  ARGUMENTS: see interface
!
!  NOTES      :	
!
!
!----------------------------------------------------------------------
!
!  EXTNL. REFERENCES:	
!
!  SYSTEM UTILITY
!  REFERENCES       :   
!
!  NOTES            :
!
!
!  REVISION HISTORY:
!  VERSION      COMMENT                         AUTHOR          DATE
!    1.0        creation                        Adri Duijndam   April 1997
!
!*******************************************************************

!------- interface ------------------------------

complex(single), intent(inout) :: data(:,:)

!------- local declarations ------------------------------

logical					:: initialized = .false.
real(single) ,	allocatable,save	:: work1(:),work2(:)
integer					:: n1,n2
integer,	save			:: nprevious1, nprevious2
integer					:: i1, i2

!-----------------------------------------------------------------------
!***	get size of data set
!-----------------------------------------------------------------------

n1 = size(data,dim=1)
n2 = size(data,dim=2)

!-----------------------------------------------------------------------
!***	initialize work array if necessary
!-----------------------------------------------------------------------

if (.not. initialized) then

    allocate(work1(1:4*n1+15))
    allocate(work2(1:4*n2+15))
    call cffti(n1,work1)
    call cffti(n2,work2)
    initialized = .true.
    nprevious1 = n1
    nprevious2 = n2

elseif (n1 /= nprevious1 .or. &
        n2 /= nprevious2      ) then

    deallocate(work1)
    deallocate(work2)

    allocate(work1(1:4*n1+15))
    allocate(work2(1:4*n2+15))
    call cffti(n1,work1)
    call cffti(n2,work2)
    nprevious1 = n1
    nprevious2 = n2

endif

!-----------------------------------------------------------------------
!***	the actual fft
!-----------------------------------------------------------------------

do i2 = 1,n2				! for all indices, dimension 2:
    call cfftb(n1,data(:,i2),work1)	! FFT along dimension 1
enddo					! 

do i1 = 1,n1				! for all indices, dimension 1:
    call cfftb(n2,data(i1,:),work2)	! FFT along dimension 2
enddo

end subroutine fftb_2d_single

!=======================================================================
!=======================================================================

subroutine fftb_2d_double(data)
!***********************************************************
!
!  TITLE      : fftb_gen
!
!  KEYWORDS   : FFT
!
!  DESCRIPTION: general input backward 2D fft
!		at this point in time, only complex data implemented
!		numerical precision: kind = double
!
!  LANGUAGE   : Fortran90
!
!  CALLING SEQUENCE:    call fftb_gen(data)
!
!  ARGUMENTS: see interface
!
!  NOTES      :	
!
!
!----------------------------------------------------------------------
!
!  EXTNL. REFERENCES:	
!
!  SYSTEM UTILITY
!  REFERENCES       :   
!
!  NOTES            :
!
!
!  REVISION HISTORY:
!  VERSION      COMMENT                         AUTHOR          DATE
!    1.0        creation                        Adri Duijndam   April 1997
!
!*******************************************************************

!------- interface ------------------------------

complex(double), intent(inout) :: data(:,:)

!------- local declarations ------------------------------

logical					:: initialized = .false.
real(double) ,	allocatable,save	:: work1(:),work2(:)
integer					:: n1,n2
integer,	save			:: nprevious1, nprevious2
integer					:: i1, i2

!-----------------------------------------------------------------------
!***	get size of data set
!-----------------------------------------------------------------------

n1 = size(data,dim=1)
n2 = size(data,dim=2)

!-----------------------------------------------------------------------
!***	initialize work array if necessary
!-----------------------------------------------------------------------

if (.not. initialized) then

    allocate(work1(1:4*n1+15))
    allocate(work2(1:4*n2+15))
    call zffti(n1,work1)
    call zffti(n2,work2)
    initialized = .true.
    nprevious1 = n1
    nprevious2 = n2

elseif (n1 /= nprevious1 .or. &
        n2 /= nprevious2      ) then

    deallocate(work1)
    deallocate(work2)

    allocate(work1(1:4*n1+15))
    allocate(work2(1:4*n2+15))
    call zffti(n1,work1)
    call zffti(n2,work2)
    nprevious1 = n1
    nprevious2 = n2

endif

!-----------------------------------------------------------------------
!***	the actual fft
!-----------------------------------------------------------------------

do i2 = 1,n2				! for all indices, dimension 2:
    call zfftb(n1,data(:,i2),work1)	! FFT along dimension 1
enddo					! 

do i1 = 1,n1				! for all indices, dimension 1:
    call zfftb(n2,data(i1,:),work2)	! FFT along dimension 2
enddo

end subroutine fftb_2d_double


end module fft_gen

