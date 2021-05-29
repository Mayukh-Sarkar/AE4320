parameters = {'kp';'Tl';'Ti';'Tp';'zeta';'omega'};
A = [y(1);0;0;y(2);y(3);y(4)];
B = [x(1);x(2);0;x(3);x(4);x(5)];
C = [z(1);z(2);z(3);z(4);z(5);z(6)];
T = table(A,B,C,'RowNames',parameters);
disp(T)