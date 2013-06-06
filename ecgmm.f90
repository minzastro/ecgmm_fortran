module ecgmm
implicit none

  real*8, parameter :: PI = 3.1415926535897932d0
  real*8, parameter :: SQRT_TWOPI = 2.506628274631d0

contains

real*8 function get_p(y, delta2, mu)
  real*8, intent(in) :: y, delta2, mu
    get_p = exp( - 0.5*(y-mu)*(y-mu)/delta2) / (SQRT_TWOPI * dsqrt(delta2))
end function get_p

subroutine p_zyt(m, y, delta, n, w, mu, sigma, p_zyt_out)
  integer, intent(in) :: m
  integer, intent(in) :: n
  real*8, intent(in) :: y(m), delta(m)
  real*8, intent(in) :: w(n), mu(n), sigma(n)
  real*8, intent(out) :: p_zyt_out(m, n)
  integer i, j
  real*8 p(n), sump
    do j = 1, m
      do i = 1, n
        p(i) = w(i)*get_p(y(j), (delta(j)**2 + sigma(i)**2), mu(i))
      enddo
      sump = 1./sum(p(:))
      p_zyt_out(j, :) = p(:)*sump
    enddo
end subroutine p_zyt

subroutine iterate(m, y, delta, n, w, mu, sigma, iter)
  integer, intent(in) :: m
  integer, intent(in) :: n
  integer, intent(in) :: iter
  real*8, intent(in) :: y(m), delta(m)
  real*8, intent(inout) :: w(n), mu(n), sigma(n)
  integer i, j, k
  real*8 t_w(n), t_mu(n), t_sigma(n)
  real*8 p_zyt_out(m, n)
  real*8 p_zyt_2(m, n), sum_p_zyt_2
    do j = 1, iter
        call p_zyt(m, y, delta, n, w, mu, sigma, p_zyt_out)
        do i = 1, n
          p_zyt_2(1:m, i) = p_zyt_out(1:m, i)/( 1.0 + (delta(1:m)/sigma(i))**2)
        enddo
        do i = 1, n
          sum_p_zyt_2 = 1.0 / sum(p_zyt_2(:, i))
          t_mu(i) = sum(y(:)*p_zyt_2(:, i)) * sum_p_zyt_2
          t_sigma(i) = dsqrt(sum( (y(:) - mu(i))**2 * p_zyt_2(:, i)) * sum_p_zyt_2)
          t_w(i) = sum(p_zyt_out(1:m, i))/dble(m)
        enddo
        w = t_w
        mu = t_mu
        sigma = t_sigma
    enddo
end subroutine iterate

end module
