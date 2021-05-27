
figure(1)

loglog(omega,abs(ft_dft(1:N/2)),'k',omega(pklocs),ft_pks,'ok')
xlabel('\omega[rad/s]')
ylabel('f_t')
legend('f_t','excitation frequecy')

figure(3)

loglog(omega,Sef(1:N/2),'k',omega(pklocs),Sef(pklocs),'ok')

legend('cpsd of e and ft','prominent peaks')
xlabel('\omega[rad/s]')
ylabel('|S_e_f_t(j \omega)|')
grid on
figure(4)
loglog(omega,Suf(1:N/2),'k',omega(pklocs),Suf(pklocs),'ok')
legend('cpsd of u and ft','prominent peaks')
xlabel('\omega[rad/s]')
ylabel('|S_u_f_t(j \omega)|')
grid on

figure(5)
loglog(omega(pklocs),magh,'ok')

xlabel(' \omega[rad/s]')
ylabel('|H_p(j \omega)|')
axis (10.^[-0.5 1.4 0.2 1.4]) 
figure(6)
semilogx(omega(pklocs),phaseh,'ok')
xlabel(' \omega[rad/s]')
ylabel('\angle H_p(j \omega)')
