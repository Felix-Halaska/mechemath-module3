function G = backward_euler(rate_func_in, t, XA, h, XB)

    G = XA + h*rate_func_in(t+h, XB) - XB;
    
end