load('ae4320_dataset1.mat')
%% excitation frequecny for FRF
dt = t(2) -t(1);
N1 = length(t);
Fs = 1/dt;
freq = (Fs/N1)*(1:N1);
omega = 2*pi*freq;
[ft_pks, ft_pklocs] = findpeaks(ft);
figure(1)
subplot(2,1,1)
plot(freq,ft,freq(ft_pklocs),ft_pks,'or')
xlabel('time')
ylabel('forcing signal')
subplot(2,1,2)
plot(t(ft_pklocs),ft_pks)
xlabel('time')
ylabel('forcing signal')
%axis tight

[fd_pks, fd_pklocs] = findpeaks(fd);
figure(2)
subplot(2,1,1)
plot(t,fd,t(fd_pklocs),fd_pks,'or')
xlabel('time')
ylabel('forcing signal')
subplot(2,1,2)
plot(t(fd_pklocs),fd_pks)
xlabel('time')
ylabel('forcing signal')
%axis tight

%% Estimate of pilot control dynamics
error = identification.e;
input = identification.u;
N = length(ft_pks);
e_dft = fft(error(ft_pklocs),N);
u_dft = fft(input(ft_pklocs),N);
ft_dft = fft(ft_pks,N);
Seft = conj(e_dft(1:N)).*ft_dft(1:N)/N;
Suft = conj(u_dft(1:N)).*ft_dft(1:N)/N;


Hp = Suft./Seft;
magh = abs(Hp);
phaseh = 180*angle(Hp)/pi;
figure(3)
loglog( freq(ft_pklocs), magh) ;
figure(4)
semilogx( freq(ft_pklocs), phaseh) ;



 