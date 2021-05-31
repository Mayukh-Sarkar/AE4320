% 
% figure(1)
% 
% loglog(omega,abs(ft_dft(1:N/2)),'k',omega(pklocs),ft_pks,'ok')
% xlabel('\omega[rad/s]')
% ylabel('f_t')
% legend('f_t','excitation frequecy')
% grid on
% 
% figure(3)
% 
% loglog(omega,Sef(1:N/2),'k',omega(pklocs),Sef(pklocs),'ok')
% 
% legend('cpsd of e and ft','prominent peaks')
% xlabel('\omega[rad/s]')
% ylabel('|S_e_f_t(j \omega)|')
% grid on
% figure(4)
% loglog(omega,Suf(1:N/2),'k',omega(pklocs),Suf(pklocs),'ok')
% legend('cpsd of u and ft','prominent peaks')
% xlabel('\omega[rad/s]')
% ylabel('|S_u_f_t(j \omega)|')
% grid on
% 
% figure(5)
% subplot(2,1,1)
% loglog(omega(pklocs),magh,'ok')
% xlabel(' \omega[rad/s]')
% ylabel('|H_p(j \omega)|')
% axis (10.^[-0.5 1.4 0.2 1.4]) 
% grid on
% subplot(2,1,2)
% semilogx(omega(pklocs),phaseh,'*k')
% xlabel(' \omega[rad/s]')
% ylabel('\angle H_p(j \omega)')
% grid on

figure(6)
set(gcf, 'Position', [100 100 700 650])
subplot(2,1,1)
loglog(omega(pklocs), magh,'ok')
hold on
%loglog(omega_test, data_magA(1,:),'-k')
% hold on
plot(omega_test, data_magB(1,:),'--k')
% hold on
%loglog(omega_test, data_magC(1,:),'-.k')
hold off
axis(10.^[-.5 1.5 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
%legend('Pilot FRF','Model A','Model B','Model C','Location','northwest')
ah=gca; 
set(ah,'Fontsize',12)

grid on

subplot(2,1,2)
semilogx(omega(pklocs), phaseh,'*k')
hold on
 %semilogx(omega_test, data_phaseA(1,:),'-k')
% hold on
semilogx(omega_test, data_phaseB(1,:),'--k')
% hold on
%semilogx(omega_test, data_phaseC(1,:),'-.k')
hold off
axis([10.^-.5 10.^1.5 -300 200])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
%legend('Pilot FRF','Model A','Model B','Model C','Location','northeast')
ah=gca; 
set(ah,'Fontsize',12)
grid on
