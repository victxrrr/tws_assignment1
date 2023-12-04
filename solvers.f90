module Solvers
    use Solver_gfortran, only: solve

    implicit none   
    save
    private
    public EulerForward, Heun, EulerBackward

    integer, parameter ::             &
    sp = selected_real_kind(6, 37),   &
    dp = selected_real_kind(15, 307), &
    ep = selected_real_kind(18,4931), &
    qp = selected_real_kind(33,4931), &
    p = sp

    real(p) :: beta, mu, gamma, alpha, delta

contains

    !> Auxiliary function for initializing the state vector of the SIQRD model
    !! @return x0 = [S_0, I_0, Q_0, R_0, D_0]
    function initialize() result(x0)
        real(p) :: x0(5)
        x0 = 0._p

        open(7,file='parameters.in')
        read(7,*) beta,mu,gamma,alpha,delta,x0(1),x0(2)
        close(7)
    end function

    !> Auxiliary function for extracting the system parameters from parameters.in  describing the SIQRD model 
    !! and for returning the corresponding evalutation of xk as a function of these parameters
    !! @param xk state vector at time tk
    !! @return eval = f([params], xk)
    function f(xk) result(eval)
        real(p), intent(in)          :: xk(:)
        real(p), dimension(size(xk)) :: eval

        eval(1) = -beta*xk(1)*xk(2)/(xk(1) + xk(2) + xk(4)) + mu*xk(4)
        eval(2) = (beta*xk(1)/(xk(1) + xk(2) + xk(4)) - gamma - delta - alpha)*xk(2)
        eval(3) = delta*xk(2) - (gamma + alpha)*xk(3)
        eval(4) = gamma*(xk(2) + xk(3)) - mu*xk(4)
        eval(5) = alpha*(xk(2) + xk(3))

    end function

    !> Auxiliary function for extracting the system paramteres from parameters.in describing the SIQRD model
    !! and for returning the jacobian matrix of f evaluated at xk
    !! @param xk
    !! @return eval = df/dx([params], xk)
    function jac_f(xk) result(eval)
        real(p), intent(in)                    :: xk(:)
        real(p), dimension(size(xk), size(xk)) :: eval
        real(p)                                :: tmp

        tmp = (xk(1) + xk(2) + xk(4))**2

        eval(1,1)   = -beta * (xk(2)*xk(2) + xk(2)*xk(4))/tmp
        eval(2,1)   = -eval(1,1)
        eval(3:,1)  = 0._p
        eval(1,2)   = -beta * (xk(1)*xk(1) + xk(1)*xk(4))/tmp
        eval(2,2)   = -eval(2,1)
        eval(3,2)   = delta
        eval(4,2)   = gamma
        eval(5,2)   = alpha
        eval(1:2,3) = 0._p
        eval(3,3)   = -(gamma + alpha)
        eval(4,3)   = gamma
        eval(5,3)   = alpha
        eval(1,4)   = beta*xk(1)*xk(2)/tmp + mu
        eval(2,4)   = -(eval(4,1) - mu)
        eval(3,4)   = 0._p
        eval(4,4)   = -mu
        eval(5,4)   = 0._p
        eval(:,5)   = 0._p
   
    end function

    !> Simple vector's euclidian norm computation
    !! @param x input vector
    !! @return n norm of x
    function norm(x) result(n)
        real(p), intent(in) :: x(:)
        real(p)             :: n
        integer             :: k

        n = 0._p
        do k = 1,size(x)
            n =  n + x(k)*x(k)
        enddo
        n = sqrt(n)
        
    end function
            

    !> Euler's forward method tailored to the SIQRD model
    !! @param T simulation horizon
    !! @param N number of discretization steps
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    subroutine EulerForward(T, N, tk, x)
        real(p), intent(in)  :: T
        integer, intent(in)  :: N
        real(p), intent(out) :: tk(:), x(:,:)
        integer              :: k

        tk(1) = 0._p
        x(1,:) = initialize()
        do k = 1,N
            tk(k+1) = (T/N)*k
            x(k+1,:) = x(k,:) + (T/N)*f(x(k,:))
        enddo

        call PrintAndSave(tk, x)
        
    end subroutine

    !> Heun's method tailored to the SIQRD model
    !! @param T simulation horizon
    !! @param N number of discretization steps
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    subroutine Heun(T, N, tk, x)
        real(p), intent(in)  :: T
        integer, intent(in)  :: N
        real(p), intent(out) :: tk(:), x(:,:)
        integer              :: k

        tk(1) = 0._p
        x(1,:) = initialize()
        do k = 1,N
            tk(k+1) = (T/N)*k
            x(k+1,:) = x(k,:) + (T/N)*(0.5*f(x(k,:)) + 0.5*f(x(k,:) + (T/N)*f(x(k,:))))
        enddo

        call PrintAndSave(tk, x)

    end subroutine

    !> Euler's backward method tailored to the SIQRD model
    !! @param T simulation horizon
    !! @param N number of discretization steps
    !! @param tk array containing the N + 1 equidistant grid points
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    !! @param tol tolerance for the stopping criterion of Newton's method
    subroutine EulerBackward(T, N, tk, x, tol)
        real(p), intent(in)                     :: T, tol
        integer, intent(in)                     :: N
        real(p), intent(out)                    :: tk(:), x(:,:)
        real(p), dimension(size(x,2),size(x,2)) :: A
        real(p)                                 :: bx(size(x,2))
        integer                                 :: k, i, j

        j = 0

        tk(1) = 0._p
        x(1,:) = initialize()
        do k = 1,N
            tk(k+1) = (T/N)*k
            
            ! Euler's method with Backward Error stopping criterion
            x(k+1,:) = x(k,:)
            do while (.true.)
                bx = x(k, :) + (T/N)*f(x(k+1,:)) - x(k+1,:)
                if (norm(bx) < tol) exit
                A = (T/N)*jac_f(x(k+1,:))
                do i = 1, size(A,1)
                    A(i,i) = A(i,i) - 1._p
                enddo
    
                call solve(A,bx)
                x(k+1,:) = x(k+1,:) - bx
            enddo        
        enddo

        call PrintAndSave(tk, x)

    end subroutine

    !> Auxiliary function for displaying the calculated solution on the flow output and 
    !! saving it in a sol.dat file
    !! @param x 2D-array containing the solution approximation at time 0, T/N, 2T/N, ... , T
    subroutine PrintAndSave(tk, x)
        real(p), intent(in) :: tk(:), x(:,:)
        integer             :: k

        open(7, file='sol.dat')
        do k = 1, size(x,1)
            write(*, '(f12.2, 5f16.5)') tk(k), x(k, 1), x(k, 2), x(k, 3), x(k, 4), x(k, 5)
            write(7, '(f12.2, 5f16.5)') tk(k), x(k, 1), x(k, 2), x(k, 3), x(k, 4), x(k, 5)
        enddo
        close(7)
    end subroutine

end module