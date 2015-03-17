MODULE lfsr_mod
! L'Ecuyer's 1999 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! This version requires that the default integer type is of 32-bits
! http://www.ozemail.com.au/~milleraj
! http://users.bigpond.net.au/amiller/
! Latest revision - 12 January 2001

IMPLICIT NONE
!INTEGER, PARAMETER :: dp = 8 !SELECTED_REAL_KIND(14, 60)
! These are unsigned integers in the C version
INTEGER, SAVE      :: s1 = 153587801, s2 = -759022222, s3 = 1288503317, &
                      s4 = -1718083407

CONTAINS

subroutine init_lfsr_time()
  integer seed(8), i
  integer,dimension(8) :: values

  ! using keyword arguments          
  call date_and_time(VALUES=values)  

  !CALL RANDOM_SEED(SIZE=K)
  seed(:) = 0D0            
  do i = 1, 8              
    seed(i) = 100000+120010*values(9-i)
  enddo                                
  call init_seeds(seed(1), seed(2), seed(3), seed(4))

end subroutine init_lfsr_time

SUBROUTINE init_seeds(i1, i2, i3, i4)
IMPLICIT NONE

INTEGER, INTENT(IN) :: i1, i2, i3, i4

s1 = i1
s2 = i2
s3 = i3
s4 = i4
IF (IAND(s1,  -2) == 0) s1 = i1 - 1023
IF (IAND(s2,  -8) == 0) s2 = i2 - 1023
IF (IAND(s3, -16) == 0) s3 = i3 - 1023
IF (IAND(s4,-128) == 0) s4 = i4 - 1023

RETURN
END SUBROUTINE init_seeds



FUNCTION lfsr113() RESULT(random_numb)
! Generates a random number between 0 and 1.  Translated from C function in:
! Reference:
! L'Ecuyer, P. (1999) `Tables of maximally equidistributed combined LFSR
! generators', Math. of Comput., 68, 261-269.

! The cycle length is claimed to be about 2^(113) or about 10^(34).
! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).(2^25 - 1)

IMPLICIT NONE
REAL *8 :: random_numb

INTEGER   :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(s1,6), s1), -13)
s1 = IEOR( ISHFT( IAND(s1,-2), 18), b)
b  = ISHFT( IEOR( ISHFT(s2,2), s2), -27)
s2 = IEOR( ISHFT( IAND(s2,-8), 2), b)
b  = ISHFT( IEOR( ISHFT(s3,13), s3), -21)
s3 = IEOR( ISHFT( IAND(s3,-16), 7), b)
b  = ISHFT( IEOR( ISHFT(s4,3), s4), -12)
s4 = IEOR( ISHFT( IAND(s4,-128), 13), b)

! The constant below is the reciprocal of (2^32 - 1)
random_numb = IEOR( IEOR( IEOR(s1,s2), s3), s4) * 2.3283064365E-10 + 0.5

RETURN
END FUNCTION lfsr113

END MODULE lfsr_mod
