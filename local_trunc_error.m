t_ref = 0.492;
h_step = logspace(-5, 1, 100);

func1_XA = solution01(t_ref);
func2_XA = solution02(t_ref);

func1_fe_err_list = [];
func1_mp_err_list = [];
func2_fe_err_list = [];
func2_mp_err_list = [];
func1_gt_err_list = [];
func2_gt_err_list = [];


for i = 1:length(h_step)
    func1_gt = solution01(t_ref + h_step(i));
    func2_gt = solution02(t_ref + h_step(i));

    [XB_func1_fe, ~] = forward_euler_step(@rate_func01, t_ref, func1_XA, h_step(i));
    [XB_func1_mp, ~] = explicit_midpoint_step(@rate_func01, t_ref, func1_XA, h_step(i));

    func1_fe_err_list = [func1_fe_err_list, norm(XB_func1_fe - func1_gt)];
    func1_mp_err_list = [func1_mp_err_list, norm(XB_func1_mp - func1_gt)];

    [XB_func2_fe, ~] = forward_euler_step(@rate_func02, t_ref, func2_XA, h_step(i));
    [XB_func2_mp, ~] = explicit_midpoint_step(@rate_func02, t_ref, func2_XA, h_step(i));

    func2_fe_err_list = [func2_fe_err_list, norm(XB_func2_fe - func2_gt)];
    func2_mp_err_list = [func2_mp_err_list, norm(XB_func2_mp - func2_gt)];

    func1_gt_err_list = [func1_gt_err_list, norm(func1_gt - func1_XA)];
    func2_gt_err_list = [func2_gt_err_list, norm(func2_gt - func2_XA)];

end

filtered_h_step = h_step(h_step <= 1e0);
filt_length = length(filtered_h_step);

filtered_fe_err = func1_fe_err_list(1:filt_length);
filtered_mp_err = func1_mp_err_list(1:filt_length);

filtered_gt_err = func1_gt_err_list(1:filt_length);

mp_p = polyfit(log(filtered_h_step), log(filtered_mp_err), 1)
fe_p = polyfit(log(filtered_h_step), log(filtered_fe_err), 1)

gt_p = polyfit(log(filtered_h_step), log(filtered_gt_err), 1)

f_func1_mp = exp(polyval(mp_p, log(filtered_h_step)));
f_func1_fe = exp(polyval(fe_p, log(filtered_h_step)));
f_func1_gt = exp(polyval(gt_p, log(filtered_h_step)));


figure();
scatter(filtered_h_step, filtered_mp_err, "filled", "r")
hold on
scatter(filtered_h_step, filtered_fe_err,"filled", "g")
scatter(filtered_h_step, filtered_gt_err, "filled", "b");
plot(filtered_h_step, f_func1_mp, '--', LineWidth=2, Color="k")
plot(filtered_h_step, f_func1_fe, '--', LineWidth=2, Color="k")
plot(filtered_h_step, f_func1_gt, '--', LineWidth=2, Color="k")
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend(["Explicit Midpoint Raw","Forward Euler Raw", "Ground Truth Raw", "Explicit Midpoint Fit", "Forward Euler Fit", "Ground Truth Fit"], Location = "northwest", Interpreter="latex", FontSize=14)
xlabel("h step (-)", Interpreter="latex", FontSize=14)
ylabel("error (-)", Interpreter="latex", FontSize=14)
% axis([0 4*pi -1.5 1.5])
title("Test Func 1: Explicit Methods Local Error", Interpreter='latex', FontSize=16)
ax = gca;
ax.FontSize = 14;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
exportgraphics(gca, 'plots/explicit_methods_local_error_testfunc1.png', 'Resolution',300)



