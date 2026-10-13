function G = implicit_midpoint(rate_func_in, t, XA, h, XB)

    G = XA + h*rate_func_in(t+h/2, (XA+XB)/2) - XB;
    
end