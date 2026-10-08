%Runs numerical integration using explicit midpoint approximation
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0';X1';X2';...;(X_end)'] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_midpoint_fixed_step_integration_test(rate_func_in,tspan,X0,h_ref)
    n_steps = (tspan(2)-tspan(1))/h_ref + 1;

    h_avg = (tspan(2) - tspan(1))/n_steps;
    num_evals = 0;
    XA = X0;
    t_list = linspace(tspan(1), tspan(2), n_steps);
    X_list = [];
    
    for i = 1:n_steps
        [XB, evals] = explicit_midpoint_step(rate_func_in, t_list(i), XA, h_avg);
        num_evals = num_evals + evals;
        X_list = [X_list, XB];
        XA = XB;
    end

end

% 
% [t_list_small, X_list_small, h_avg, num_evals] = explicit_midpoint_fixed_step_integration_test(@rate_func01, [0,4*pi], 1, 0.001);
% [t_list_med, X_list_med, h_avg, num_evals] = explicit_midpoint_fixed_step_integration_test(@rate_func01, [0,4*pi], 1, 0.01);
% [t_list_large, X_list_large, h_avg, num_evals] = explicit_midpoint_fixed_step_integration_test(@rate_func01, [0,4*pi], 1, 0.1);
% 
% num_evals
% 
% figure();
% t_space = linspace(0,4*pi,2*pi/0.01);
% true_sol = cos(t_space);
% plot(t_space, true_sol, LineWidth=3, Color="k")
% hold on
% plot(t_list_small, X_list_small, '--', LineWidth=2, Color="r")
% plot(t_list_med, X_list_med, '--', LineWidth=2, Color="g")
% plot(t_list_large, X_list_large, '--', LineWidth=2, Color="b")
% hold off
% legend(["True","h = 0.01", "h = 0.05", "h = 0.1"], Interpreter="latex", Location="southeast", FontSize=14)
% xlabel("x (-)", Interpreter="latex", FontSize=14)
% ylabel("f(x) (-)", Interpreter="latex", FontSize=14)
% axis([0 4*pi -1.5 1.5])
% title("Explicit Midpoint Method Timestep", Interpreter='latex', FontSize=16)
% ax = gca;
% ax.FontSize = 14;
% ax.XAxis.TickLabelInterpreter = 'Latex';
% ax.YAxis.TickLabelInterpreter = 'Latex';
% exportgraphics(gca, 'plots/explicit_midpoint.png', 'Resolution',300)