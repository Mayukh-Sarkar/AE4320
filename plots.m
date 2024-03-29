% f_t peaks plot
figure(1)

loglog(omega,abs(ft_dft(1:N/2)),'k',omega(pklocs),ft_pks,'ok')
xlabel('\omega[rad/s]')
ylabel('f_t')
legend('f_t','excitation frequecy')
grid on

%f_d peaks plot
figure(2)

loglog(omega,abs(fd_dft(1:N/2)),'k',omega(pklocsf),fd_pks,'ok')
xlabel('\omega[rad/s]')
ylabel('f_d')
legend('f_d','excitation frequecy')
grid on

%cpsd of ft and error at prominent peaks-excitation freq
figure(3)

loglog(omega,Sef(1:N/2),'k',omega(pklocs),Sef(pklocs),'ok')

legend('cpsd of e and ft','prominent peaks')
xlabel('\omega[rad/s]')
ylabel('|S_e_f_t(j \omega)|')
grid on

%cpsd of ft and u at excitation freq, prominet peaks
figure(4)
loglog(omega,Suf(1:N/2),'k',omega(pklocs),Suf(pklocs),'ok')
legend('cpsd of u and ft','prominent peaks')
xlabel('\omega[rad/s]')
ylabel('|S_u_f_t(j \omega)|')
grid on

% pilot dynamics FRF
figure(5)
subplot(2,1,1)
loglog(omega(pklocs),magh,'ok')
xlabel(' \omega[rad/s]')
ylabel('|H_p(j \omega)|')
axis (10.^[-0.5 1.4 0.2 1.4]) 
grid on
subplot(2,1,2)
semilogx(omega(pklocs),phaseh,'*k')
xlabel(' \omega[rad/s]')
ylabel('\angle H_p(j \omega)')
grid on

%paramter estimation and fitting of models A,B and C

figure(6)
set(gcf, 'Position', [100 100 700 650])
subplot(2,1,1)
loglog(omega(pklocs), magh,'ok')
hold on
loglog(omega_test, data_magA(1,:),'-k')
 hold on
plot(omega_test, data_magB(1,:),'--k')
 hold on
loglog(omega_test, data_magC(1,:),'-.k')
hold on
%loglog(omega_test, data_magD(1,:),'-r')
hold off
axis(10.^[-.5 1.5 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
legend('Pilot FRF','Model A','Model B','Model C','Location','northwest')
ah=gca; 
set(ah,'Fontsize',12)

grid on

subplot(2,1,2)
semilogx(omega(pklocs), phaseh,'*k')
hold on
semilogx(omega_test, data_phaseA(1,:),'-k')
 hold on
semilogx(omega_test, data_phaseB(1,:),'--k')
 hold on
semilogx(omega_test, data_phaseC(1,:),'-.k')
hold on
%semilogx(omega_test, data_phaseD(1,:),'-r')
hold off
axis([10.^-.5 10.^1.5 -300 200])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
legend('Pilot FRF','Model A','Model B','Model C','Location','northeast')
ah=gca; 
set(ah,'Fontsize',12)
grid on

%cpsd of e,u and ft using validation data set and comparing it to the cpsd
%from identification data set
figure(7)
subplot(2,1,1)
loglog(omega,Sef_v(1:N/2),'r',omega(pklocs),Sef_v(pklocs),'or')
hold on
loglog(omega,Sef(1:N/2),'k',omega(pklocs),Sef(pklocs),'ok')
hold off
legend('cpsd of e and ft of validation data','prominent peaks validation data','cpsd of e and ft of identification data','prominent peaks identification data')

xlabel('\omega[rad/s]')
ylabel('|S_e_f_t(j \omega)|')
grid on
subplot(2,1,2)
loglog(omega,Suf_v(1:N/2),'r',omega(pklocs),Suf_v(pklocs),'or')
hold on
loglog(omega,Suf(1:N/2),'k',omega(pklocs),Suf(pklocs),'ok')
legend('cpsd of u and ft of validation data','prominent peaks validation data','cpsd of u and ft of identification data','prominent peaks identification data')
xlabel('\omega[rad/s]')
ylabel('|S_u_f_t(j \omega)|')
grid on

% calucation pilot dynamcis FRF using validation data set and comparing it
% to the Hp of identification data set
figure(8)
subplot(2,1,1)
loglog(omega(pklocs),magh,'ok')
hold on
loglog(omega(pklocs),magh_v,'sk')
hold off
xlabel(' \omega[rad/s]')
ylabel('|H_p(j \omega)|[-]')
legend('Magnitude from identification data','Magnitude from validation data')
axis (10.^[-0.5 1.4 0.2 1.4]) 
ah=gca; 
set(ah,'Fontsize',12)
grid on
subplot(2,1,2)
semilogx(omega(pklocs),phaseh,'*k')
hold on
semilogx(omega(pklocs),phaseh_v,'dk')
legend('Phase from identification data','Phase from validation data')
xlabel(' \omega[rad/s]')
ylabel('\angle H_p(j \omega)[deg]')
ah=gca; 
set(ah,'Fontsize',12)
grid on

% fitting of model D and comparing it to B,C
figure(9)
set(gcf, 'Position', [100 100 700 650])
subplot(2,1,1)
loglog(omega(pklocs), magh,'ok')
hold on
%loglog(omega_test, data_magA(1,:),'-k')
 hold on
plot(omega_test, data_magB(1,:),'--k')
 hold on
loglog(omega_test, data_magC(1,:),'-.k')
hold on
loglog(omega_test, data_magD(1,:),'-k')
hold off
axis(10.^[-.5 1.5 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
legend('Pilot FRF','Model B','Model C','Model D','Location','northwest')
ah=gca; 
set(ah,'Fontsize',12)

grid on

subplot(2,1,2)
semilogx(omega(pklocs), phaseh,'*k')
hold on
%semilogx(omega_test, data_phaseA(1,:),'-k')
 hold on
semilogx(omega_test, data_phaseB(1,:),'--k')
 hold on
semilogx(omega_test, data_phaseC(1,:),'-.k')
hold on
semilogx(omega_test, data_phaseD(1,:),'-k')
hold off
axis([10.^-.5 10.^1.5 -300 200])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
legend('Pilot FRF','Model B','Model C','Model D','Location','northeast')
ah=gca; 
set(ah,'Fontsize',12)
grid on
