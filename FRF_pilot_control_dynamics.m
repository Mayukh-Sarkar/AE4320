load('ae4320_dataset1.mat')
%% excitation frequecny for FRF
dt = t(2) -t(1);
N1 = length(t);
Fs = 1/dt;
freq = (Fs/N1)*(1:N1);% freq in HZ
omega = 2*pi*freq;
[ft_pks, ft_pklocs] = findpeaks(ft); % identifying the excitation freq. Peaks of the forcing function
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

%% Estimate of cross psd
error = identification.e;
input = identification.u;
N = length(ft_pks);
Freq = (Fs/N)*(1:N);

e_dft = fft(error(ft_pklocs),N); % Fourier transform of the error signal at excitation freq
u_dft = fft(input(ft_pklocs),N); % Fourier transform of the input signal at excitation freq
ft_dft = fft(ft_pks,N);  % Fourier transform of the forcing function
Seft = conj(e_dft(1:N)).*ft_dft(1:N)/N; % cross power spectral density of error and ft
Suft = conj(u_dft(1:N)).*ft_dft(1:N)/N; %cross power spectral density of input and ft

figure(3)
subplot(2,1,1)
loglog(Freq,abs(Seft))
subplot(2,1,2)
loglog(Freq,abs(Suft))


%% Pilot dynamics
Hp = Suft./Seft; %pilot dynamics at excitation freq
magh = abs(Hp); % magnitude of Hp
phaseh = 180*angle(Hp)/pi; % phase of hp
figure(4)
loglog( Freq, magh) ;
figure(5)
semilogx( Freq, phaseh) ;



 