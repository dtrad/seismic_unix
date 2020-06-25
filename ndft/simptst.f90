program simptst

use nfft_routs

implicit none

real, dimension(8) :: xco = (/ 1.1, 1.9, 2.9, 4., 5.1, 6.3, 6.9, 8.1 /)
complex, dimension(8) :: indat = (/(1,1), (-1,1), (.8,.3), (10,-1), &
                                   (2,2), (-2,2), (.4,.6), (20,-2) /)
complex, dimension(8) :: outdat 

call nfft_init(xco, 0.78, 8, 10, 2.0)
call nfft(indat,outdat)

print *,outdat

endprogram simptst
