program ExponentialGrowth

implicit none

integer, parameter ::                 &
    sp = selected_real_kind(6, 37),   &
    dp = selected_real_kind(15, 307), &
    ep = selected_real_kind(18,4931), &
    qp = selected_real_kind(33,4931), &
    p = sp ! Q1.1 --> here you can choose the desired precision (single, double, extended, quadruple)

integer :: i0
real(p) :: d, r, id

print *,'number of initial infections ?'
read *, i0

print *,'number of days ?'
read *, d

print *,'daily infection rate ?'
read *, r

id = i0*(1 + r)**d
write(*, '(a, f7.1, a, f30.15)') 'Number of infections after ', d, ' days :', id

end program