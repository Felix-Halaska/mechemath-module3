t_ref = [0, .492];
h_step = logspace(-5, 0, 100);

func1_XA = solution01(t_ref(1));
func2_XA = solution02(t_ref(1));
func1_gt = solution01(t_ref(end));
func2_gt = solution02(t_ref(end));

func1_fe_err_list = [];
func1_mp_err_list = [];

func2_fe_err_list = [];
func2_mp_err_list = [];

fe_evals_list = [];
mp_evals_list = [];


for i = 1:length(h_step)
    [t_list_fe, x_list_fe, ~, num_evals_fe] = forward_euler_fixed_step_integration(@rate_func01, t_ref, func1_XA, h_step(i));
    [t_list_mp, x_list_mp, ~, num_evals_mp] = explicit_midpoint_fixed_step_integration(@rate_func01, t_ref, func1_XA, h_step(i));

    XB_func1_fe = x_list_fe(end);
    XB_func1_mp = x_list_mp(end);

    func1_fe_err_list = [func1_fe_err_list, norm(XB_func1_fe - func1_gt)];
    func1_mp_err_list = [func1_mp_err_list, norm(XB_func1_mp - func1_gt)];

    fe_evals_list = [fe_evals_list, num_evals_fe];
    mp_evals_list = [mp_evals_list, num_evals_mp];

    [t_list_fe2, x_list_fe2, ~, num_evals_fe2] = forward_euler_fixed_step_integration(@rate_func02, t_ref, func2_XA, h_step(i));
    [t_list_mp2, x_list_mp2, ~, num_evals_mp2] = explicit_midpoint_fixed_step_integration(@rate_func02, t_ref, func2_XA, h_step(i));

    XB_func2_fe = x_list_fe2(end);
    XB_func2_mp = x_list_mp2(end);

    func2_fe_err_list = [func2_fe_err_list, norm(XB_func2_fe - func2_gt)];
    func2_mp_err_list = [func2_mp_err_list, norm(XB_func2_mp - func2_gt)];

end
filtered_h_step = h_step(h_step <= 1e-2);
filt_length = length(filtered_h_step);

filtered_fe_err = func1_fe_err_list(1:filt_length);
filtered_mp_err = func1_mp_err_list(1:filt_length);


mp_p = polyfit(log(filtered_h_step), log(filtered_mp_err), 1)
fe_p = polyfit(log(filtered_h_step), log(filtered_fe_err), 1)


f_func1_mp = exp(polyval(mp_p, log(filtered_h_step)));
f_func1_fe = exp(polyval(fe_p, log(filtered_h_step)));


% 
% func1_mp_p = polyfit(log(h_step), log(func1_mp_err_list), 1)
% f_func1_mp = exp(polyval(func1_mp_p, log(h_step)));
% 
% func1_fe_p = polyfit(log(h_step), log(func1_fe_err_list), 1)
% f_func1_fe = exp(polyval(func1_fe_p, log(h_step)));
% 
% mp_eval_p = polyfit(log(mp_evals_list), log(mp_err_list), 1)
% f_mp = exp(polyval(mp_eval_p, log(mp_evals_list)));
% 
% fe_eval_p = polyfit(log(fe_evals_list), log(fe_err_list), 1)
% f_fe = exp(polyval(fe_eval_p, log(fe_evals_list)));

% func2_mp_p = polyfit(log(h_step), log(func2_mp_err_list), 1)
% f_func2_mp = exp(polyval(func2_mp_p, log(h_step)));
% 
% func2_fe_p = polyfit(log(h_step), log(func2_fe_err_list), 1)
% f_func2_fe = exp(polyval(func2_fe_p, log(h_step)));



figure();
scatter(filtered_h_step, filtered_mp_err, 60,'filled',"r")
hold on
scatter(filtered_h_step, filtered_fe_err,60,'filled',"b")
plot(filtered_h_step, f_func1_mp, '-', LineWidth=3, Color="k")
plot(filtered_h_step, f_func1_fe, '-', LineWidth=3, Color="k")
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend(["Explicit Midpoint Measured","Forward Euler Measured", "Explicit Midpoint Fit", "Forward Euler Fit"], Location = "northwest", Interpreter="latex", FontSize=14)
xlabel("Time Step - $h$ (-)", Interpreter="latex", FontSize=14)
ylabel("Global Truncation Error (-)", Interpreter="latex", FontSize=14)
% axis([0 4*pi -1.5 1.5])
title("Global Truncation Error vs Step Size - Explicit Methods", Interpreter='latex', FontSize=18)
ax = gca;
ax.FontSize = 16;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
exportgraphics(gca, 'plots/explicit_methods_global_error_testfunc1.png', 'Resolution',300)

% figure();
% scatter(mp_evals_list, func1_mp_err_list, "r")
% hold on
% scatter(fe_evals_list, func1_fe_err_list,"b")
% plot(mp_evals_list, f_mp, '--', LineWidth=3, Color="r")
% plot(fe_evals_list, f_fe, '--', LineWidth=3, Color="b")
% hold off
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% legend(["Explicit Midpoint Raw","Forward Euler Raw", "Explicit Midpoint Fit", "Forward Euler Fit"], Location = "northeast", Interpreter="latex", FontSize=14)
% xlabel("Number of Rate Function Calls (-)", Interpreter="latex", FontSize=14)
% ylabel("Error (-)", Interpreter="latex", FontSize=14)
% % axis([0 4*pi -1.5 1.5])
% title("Explicit Methods Global Error", Interpreter='latex', FontSize=16)
% ax = gca;
% ax.FontSize = 14;
% ax.XAxis.TickLabelInterpreter = 'Latex';
% ax.YAxis.TickLabelInterpreter = 'Latex';
% % exportgraphics(gca, 'plots/explicit_methods_global_error__func_calls_testfunc1.png', 'Resolution',300)



