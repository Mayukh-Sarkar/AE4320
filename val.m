parameters = {'kp';'Tl';'Ti';'Tp';'zeta';'omega';'cost';'validation'};
A = [y(1);0;0;y(2);y(3);y(4);modelA;modelAv];
B = [x(1);x(2);0;x(3);x(4);x(5);modelB;modelBv];
C = [z(1);z(2);z(3);z(4);z(5);z(6);modelC;modelCv];
%D = [w(1);w(2);w(3);w(4);w(5);w(6);modelD];
D = [w(1);0;w(2);w(3);w(4);w(5);modelD;modelDv];

T = table(A,B,C,D,'RowNames',parameters);
disp(T)