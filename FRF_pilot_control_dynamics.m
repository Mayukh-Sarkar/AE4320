clear;
clc;
load('ae4320_dataset1.mat')
N = 2^13;
t = t(1:N);
ft = ft(end-N+1:end);
fd = fd(end-N+1:end);
error = identification.e;
input = identification.u;
error = error(end-N+1:end);
input = input(end-N+1:end);
%% excitation frequecny for FRF
f_t = ft;
dt = t(2) -t(1);
N = length(t);
Fs = 1/dt;
%freq = (Fs/N)*(1:N);% freq in HZ

% Display results


% [fd_pks, fd_pklocs] = findpeaks(fd);
% figure(2)
% subplot(2,1,1)
% plot(t,fd,t(fd_pklocs),fd_pks,'or')
% xlabel('time')
% ylabel('forcing signal')
% subplot(2,1,2)
% plot(t(fd_pklocs),fd_pks)
% xlabel('time')
% ylabel('forcing signal')
% axis tight

%% Estimate of cross psd

Freq = (Fs/N)*(1:N/2);
omega = 2*pi*Freq;

e_dft = fft(error,N); % Fourier transform of the error signal at excitation freq
u_dft = fft(input,N); % Fourier transform of the input signal at excitation freq
ft_dft = fft(ft,N);  % Fourier transform of the forcing function
fd_dft = fft(fd,N);
[ft_pks , pklocs] = findpeaks(abs(ft_dft(1:N/2)),"MinPeakHeight", 0.1);
[fd_pks , pklocsf] = findpeaks(abs(fd_dft(1:N/2)),"MinPeakHeight", 0.1);

Seft = (e_dft(1:N)).*ft_dft(1:N)/N; % cross power spectral density of error and ft
Suft = (u_dft(1:N)).*ft_dft(1:N)/N; %cross power spectral density of input and ft
Sef =  abs(Seft);
Suf =  abs(Suft);



%% Pilot dynamics
omegaf = omega(pklocs)';
Hp = Suft(pklocs)./Seft(pklocs); %pilot dynamics at excitation freq
magh = abs(Hp); % magnitude of Hp
phaseh = rad2deg(unwrap(angle(Hp))); % phase of hp

CF = max(abs(ft))/rms(ft);
sprintf('The crest factor is %d',CF)


