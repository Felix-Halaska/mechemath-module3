t_ref = 0.492;
h_step = logspace(-5, 1, 100);

XA = solution02(t_ref);

test_func = @rate_func02;

fe_err_list = [];
mp_err_list = [];
gt_err_list = [];
be_err_list = [];
imp_err_list = [];

for i = 1:length(h_step)
    gt = solution02(t_ref + h_step(i));

    [XB_fe, ~] = forward_euler_step(test_func, t_ref, XA, h_step(i));
    [XB_mp, ~] = explicit_midpoint_step(test_func, t_ref, XA, h_step(i));
    [XB_be, ~] = backward_euler_step(test_func, t_ref, XA, h_step(i));
    [XB_imp, ~] = implicit_midpoint_step(test_func, t_ref, XA, h_step(i));


    fe_err_list = [fe_err_list, norm(XB_fe - gt)];
    mp_err_list = [mp_err_list, norm(XB_mp - gt)];
    be_err_list = [be_err_list, norm(XB_be - gt)];
    imp_err_list = [imp_err_list, norm(XB_imp - gt)];

    gt_err_list = [gt_err_list, norm(gt - XA)];
end

filtered_h_step = h_step(h_step <= 1e-1);
filt_length = length(filtered_h_step);

filtered_fe_err = fe_err_list(1:filt_length);
filtered_mp_err = mp_err_list(1:filt_length);
filtered_be_err = be_err_list(1:filt_length);
filtered_imp_err = imp_err_list(1:filt_length);
filtered_gt_err = gt_err_list(1:filt_length);


mp_p = polyfit(log(filtered_h_step), log(filtered_mp_err), 1)
fe_p = polyfit(log(filtered_h_step), log(filtered_fe_err), 1)
be_p = polyfit(log(filtered_h_step), log(filtered_be_err), 1)
imp_p = polyfit(log(filtered_h_step), log(filtered_imp_err), 1)

gt_p = polyfit(log(filtered_h_step), log(filtered_gt_err), 1)


figure();
scatter(filtered_h_step, filtered_mp_err, "r")
hold on
scatter(filtered_h_step, filtered_fe_err, "g")
scatter(filtered_h_step, filtered_be_err, "b");
scatter(filtered_h_step, filtered_imp_err, "k");
% scatter(h_step, gt_err_list);
hold off
set(gca,'xscale','log')
set(gca,'yscale','log')
legend(["Explicit Midpoint","Forward Euler", "Backward Euler", "Implicit Midpoint", "Ground Truth"], Location = "northwest", Interpreter="latex", FontSize=14)
xlabel("h step (-)", Interpreter="latex", FontSize=14)
ylabel("error (-)", Interpreter="latex", FontSize=14)
axis([1e-5 1e-1 1e-15 1e0])
title("All Methods Local Error", Interpreter='latex', FontSize=16)
ax = gca;
ax.FontSize = 14;
ax.XAxis.TickLabelInterpreter = 'Latex';
ax.YAxis.TickLabelInterpreter = 'Latex';
% exportgraphics(gca, 'plots/all_methods_local_error_testfunc1.png', 'Resolution',300)



