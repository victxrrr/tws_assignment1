program Siqrd
    use Solvers

    integer, parameter ::             &
    sp = selected_real_kind(6, 37),   &
    dp = selected_real_kind(15, 307), &
    ep = selected_real_kind(18,4931), &
    qp = selected_real_kind(33,4931), &
    p  = sp

    integer , parameter :: N = 150
    real(p)            :: T, x(N+1,5), tk(N+1)
    T = 30._p
    

    call EulerForward(T, N, tk, x)
    call Heun(T, N, tk, x)
    call EulerBackward(T, N, tk, x, 1.0e-5)

end program