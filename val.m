parameters = {'kp';'Tl';'Ti';'Tp';'zeta';'omega'};
A = [Kp_A;0;0;tau_A;zeta_A;omega_A];
B = [x(1);x(2);0;x(3);x(4);x(5)];
C = [Kp_C;TL_C;TI_C;tau_C;zeta_C;omega_C];
T = table(A,B,C,'RowNames',parameters);
disp(T)