function [j, num_calls] = j_approximate(fun, x)
    
    num_calls = 0;
    alpha = 1e-6;

    f0 = fun(x);
    num_calls = num_calls + 1;

    in_len = length(x);
    out_len = length(f0);

    dx = zeros(size(x));

    j = zeros([out_len, in_len]);

    for i = 1:in_len

        dx(i) = alpha;

        f_left = fun(x-dx);
        f_right = fun(x+dx);
        num_calls = num_calls + 2;

        j(:, i) = (f_right - f_left) / (2*alpha);

        dx = zeros(size(x));

    end

end

%test_numerical_jacobian(@j_approximate)