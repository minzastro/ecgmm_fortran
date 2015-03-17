module ecgmm
use lfsr_mod
use sorting
!This is Fortran-90 implementation of ECGMM algorithm presented in:
!Hao, J. et al. (2009).
!Precision Measurements of the Cluster Red Sequence Using an Error-Corrected Gaussian Mixture Model.
!The Astrophysical Journal, 702(1), 745–758.
!
! Implemented by A.A. Mints
!
! Designed to be used with Python via f2py
!
! See formulas in the paper.

implicit none

  real*8, parameter :: PI = 3.1415926535897932d0
  real*8, parameter :: SQRT_TWOPI = 2.506628274631d0
  real*8, parameter :: PARAM_PRECISION = 1d-6
  integer, parameter :: MAX_BOOTSTRAP = 50
  integer, parameter :: MAX_ITERATIONS = 1000 ! Unused
  real, parameter :: MIN_SIGMA = 1e-10 ! Minimum value of sigma

contains

real*8 function get_p(y, delta2, mu)
  !Formula (A4)
  real*8, intent(in) :: y, delta2, mu
    get_p = exp( - 0.5*(y-mu)*(y-mu)/delta2) / (SQRT_TWOPI * dsqrt(delta2))
end function get_p

subroutine p_zyt(m, y, delta, n, w, mu, sigma, p_zyt_out)
  !Formula (A5)
  integer, intent(in) :: m ! Number of datapoints
  integer, intent(in) :: n ! Number of gaussians
  real*8, intent(in) :: y(m), delta(m) ! Values and errors
  real*8, intent(in) :: w(n), mu(n), sigma(n) ! Weight, offset and width
  real*8, intent(out) :: p_zyt_out(m, n)
  integer i, j
  real*8 p(n), sump
    do j = 1, m
      do i = 1, n
        p(i) = w(i)*get_p(y(j), (delta(j)**2 + sigma(i)**2), mu(i))
      enddo
      sump = 1./sum(p(:))
      if (sump.gt.HUGE(1.0)) then ! For very small p's
        sump = 1.
      endif
      p_zyt_out(j, :) = p(:)*sump
    enddo
end subroutine p_zyt

subroutine one_step(m, y, delta, n, w, mu, sigma, t_w, t_mu, t_sigma)
  integer, intent(in) :: m ! Number of datapoints
  integer, intent(in) :: n ! Number of gaussians
  real*8, intent(in) :: y(m), delta(m) ! Values and errors
  real*8, intent(in) :: w(n), mu(n), sigma(n) ! Weight, offset and width
  real*8, intent(out) :: t_w(n), t_mu(n), t_sigma(n) ! Weight, offset and width
  integer i
  real*8 p_zyt_out(m, n)
  real*8 p_zyt_2(m, n), sum_p_zyt_2
  real*8 p_zyt_3(m, n)
    call p_zyt(m, y, delta, n, w, mu, sigma, p_zyt_out)
    do i = 1, n
      ! Update to speed-up (A9) and (A11)
      p_zyt_2(1:m, i) = p_zyt_out(1:m, i)/( 1.0 + (delta(1:m)/sigma(i))**2)
      p_zyt_3(1:m, i) = p_zyt_out(1:m, i)/( 1.0 + (delta(1:m)/sigma(i))**2)**2
    enddo
    do i = 1, n
      sum_p_zyt_2 = 1.0 / sum(p_zyt_2(:, i))
      t_mu(i) = sum(y(:)*p_zyt_2(:, i)) * sum_p_zyt_2 ! (A9)
      t_sigma(i) = dsqrt(sum( (y(:) - mu(i))**2 * p_zyt_3(:, i)) * sum_p_zyt_2) ! (A11)
      if (t_sigma(i).lt.MIN_SIGMA) then ! To prevent overfitting...
        t_sigma(i) = MIN_SIGMA
      endif
      t_w(i) = sum(p_zyt_out(1:m, i))/dble(m) ! (A14)
    enddo
end subroutine one_step

subroutine iterate(m, y, delta, n, w, mu, sigma, iter)
  integer, intent(in) :: m ! Number of datapoints
  integer, intent(in) :: n ! Number of gaussians
  integer, intent(in) :: iter
  real*8, intent(in) :: y(m), delta(m) ! Values and errors
  real*8, intent(inout) :: w(n), mu(n), sigma(n) ! Weight, offset and width
  integer j
  real*8 t_w(n), t_mu(n), t_sigma(n)
    do j = 1, iter
        call one_step(m, y, delta, n, w, mu, sigma, t_w, t_mu, t_sigma)
        !write(33, *) w, mu, sigma, t_w, t_mu, t_sigma
        w = t_w
        mu = t_mu
        sigma = t_sigma
    enddo
end subroutine iterate

logical function is_precise(n, w_old, w_new)
  !Check if all values are changed by relative amount
  !smaller than PARAM_PRECISION.
  integer, intent(in) :: n
  real*8, intent(in) :: w_old(n), w_new(n)
  if ( all( w_old(:).eq.0d0) ) then
    if ( all( w_new(:).eq.0d0) ) then
      is_precise = .True.
    else
      is_precise = .False.
    endif
  elseif (all( abs((w_old(:) - w_new(:))/w_old(:)) .lt. PARAM_PRECISION)) then
    is_precise = .True.
  else
    is_precise = .False.
  endif
end function is_precise

subroutine iterate_safe(m, y, delta, n, w, mu, sigma, iter)
  integer, intent(in) :: m ! Number of datapoints
  integer, intent(in) :: n ! Number of gaussians
  integer, intent(in) :: iter ! Maximum number of iterations allowed
  real*8, intent(in) :: y(m), delta(m) ! Values and errors
  real*8, intent(inout) :: w(n), mu(n), sigma(n) ! Weight, offset and width

  integer j
  real*8 t_w(n), t_mu(n), t_sigma(n)
    j = 0
    do while(j.le.iter)
      call one_step(m, y, delta, n, w, mu, sigma, t_w, t_mu, t_sigma)
      ! Check if desired precision is reached (all |dA/A| < \delta)
      if ((is_precise(n, w, t_w)).and.(is_precise(n, mu, t_mu)).and.(is_precise(n, sigma, t_sigma))) then
        exit
      endif
      w = t_w
      mu = t_mu
      sigma = t_sigma
      j = j + 1
    enddo
end subroutine iterate_safe

subroutine iterate_bootstrap(m, y, delta, n, w, mu, sigma, iter)
  integer, intent(in) :: m ! Number of datapoints
  integer, intent(in) :: n ! Number of gaussians
  integer, intent(in) :: iter ! Maximum number of iterations allowed
  real*8, intent(in) :: y(m), delta(m) ! Values and errors
  real*8, intent(inout) :: w(n), mu(n), sigma(n) ! Weight, offset and width

  integer i, j, k
  integer mask(m, MAX_BOOTSTRAP)
  integer isortedw(n)
  real*8 x_w(n, MAX_BOOTSTRAP), x_mu(n, MAX_BOOTSTRAP), x_sigma(n, MAX_BOOTSTRAP)
  real*8 t_w(n, MAX_BOOTSTRAP), t_mu(n, MAX_BOOTSTRAP), t_sigma(n, MAX_BOOTSTRAP)

    call init_lfsr_time()
    do i = 1, m
      do j = 1, MAX_BOOTSTRAP
        mask(i, j) = int(m*lfsr113() + 1)
      enddo
    enddo
    do i = 1, MAX_BOOTSTRAP
      x_w(:, i) = w(:)
      x_mu(:, i) = mu(:)
      x_sigma(:, i) = sigma(:)
    enddo
    j = 0
    do while(j.le.iter)
      do k = 1, MAX_BOOTSTRAP
        call one_step(m, y, delta, n, x_w(:, k), x_mu(:, k), x_sigma(:, k), &
                      t_w(:, k), t_mu(:, k), t_sigma(:, k))
      enddo
      ! Check if desired precision is reached (all |dA/A| < \delta)
      if ((is_precise(n*MAX_BOOTSTRAP, x_w, t_w)).and. &
          (is_precise(n*MAX_BOOTSTRAP, x_mu, t_mu)).and. &
          (is_precise(n*MAX_BOOTSTRAP, x_sigma, t_sigma))) then
        exit
      endif
      x_w = t_w
      x_mu = t_mu
      x_sigma = t_sigma
      j = j + 1
    enddo
    do k = 1, MAX_BOOTSTRAP
      call ssort_index(x_w(:, k), x_w(:, k), isortedw(:), n)
      x_mu(:, k) = x_mu(isortedw(:), k)
      x_sigma(:, k) = x_sigma(isortedw(:), k)
    enddo
    do i = 1, n
      w(i) = sum(x_w(i, :))/dble(MAX_BOOTSTRAP)
      mu(i) = sum(x_mu(i, :))/dble(MAX_BOOTSTRAP)
      sigma(i) = sum(x_sigma(i, :))/dble(MAX_BOOTSTRAP)
    enddo
end subroutine iterate_bootstrap


subroutine iterate_fixed(m, y, delta, n, w, mu, sigma, fixed, iter)
  integer, intent(in) :: m ! Number of datapoints
  integer, intent(in) :: n ! Number of gaussians
  integer, intent(in) :: iter
  real*8, intent(in) :: y(m), delta(m) ! Values and errors
  logical, intent(in) :: fixed(3) ! Flags for varying w, mu and sigma
  real*8, intent(inout) :: w(n), mu(n), sigma(n) ! Weight, offset and width
  integer i, j
  real*8 t_w(n), t_mu(n), t_sigma(n)
        j = 0
    do while(j.le.iter)
        call one_step(m, y, delta, n, w, mu, sigma, t_w, t_mu, t_sigma)
        !write(88, *) sigma, t_sigma
      ! Check if desired precision is reached (all |dA/A| < \delta)
        if ((is_precise(n, w, t_w).or.(fixed(1))).and. &
            (is_precise(n, mu, t_mu).or.(fixed(2))).and. &
            (is_precise(n, sigma, t_sigma).or.(fixed(3)))) then
          exit
        endif
        if (.not.fixed(1)) then
          w = t_w
        endif
        if (.not.fixed(2)) then
          mu = t_mu
        endif
        if (.not.fixed(3)) then
          sigma = t_sigma
        endif
        j = j + 1
    enddo
end subroutine iterate_fixed

end module ecgmm
