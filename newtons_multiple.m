function [x0, exit_flag, num_calls] = newton_solver_mult(fun,x0,dxtol,ftol,max_iter, dxmax, j_type)
    iter = 1;
    num_calls = 0;
    exit_flag = 0;

    if j_type == 0

        [f_x, j] = fun(x0);
    else
      
        f_x = fun(x0);
        [j, j_num_calls] = j_approx(fun, x0);
        num_calls = num_calls + j_num_calls + 1;
    end

    x1 = x0 - j\f_x;


    while sum(abs(f_x)) > ftol && sum(abs(x1-x0)) > dxtol

        if sum(abs(x1-x0)) > dxmax
            exit_flag = 1;
            break
        end

        x0 = x1;
        if j_type == 0

            [f_x, j] = fun(x0);
        else
            f_x = fun(x0);

            [j, j_num_calls] = j_approx(fun, x0);
            num_calls = num_calls + j_num_calls + 1;
        end
        
       if det(j*j') == 0
            exit_flag = 1;
            break
       end

       x1 = x0 - j\f_x;

       iter = iter + 1;

       if iter >= max_iter
            exit_flag = 2;
            break
       end
    end
end

% x_guess = [1;2;3];
% [root, exit] = newton_solver_mult(@test_function02, x_guess, 1e-14, 1e-14, 200, 1e4, 1);
% [test, ~] = test_function02(root)