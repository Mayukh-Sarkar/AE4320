load('ae4320_dataset1.mat')
%% excitation frequecny for FRF
dt = t(2) -t(1);
N = length(t);
Fs = 1/dt;
freq = (Fs/N)*(1:N);% freq in HZ
omega = 2*pi*freq;
[ft_pks , pklocs] = findpeaks(ft,"MinPeakProminence",0.3,"MinPeakHeight", 0.8);
pklocs(end) = N(end);
pklocs(3) =[];
ft_pks(end) = ft(end);
ft_pks(3) =[];
% Display results
figure(1)
plot(freq,ft)
hold on
scatter(freq(pklocs),ft_pks,'or','filled')
hold off
grid on
xlabel('frequency(Hz)')
ylabel('f_t')
legend('f_t','excitation frequecy')

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
error = identification.e;
input = identification.u;

Freq = (Fs/N)*(1:N/2);

e_dft = fft(error,N); % Fourier transform of the error signal at excitation freq
u_dft = fft(input,N); % Fourier transform of the input signal at excitation freq
ft_dft = fft(ft,N);  % Fourier transform of the forcing function
Seft = conj(e_dft(1:N)).*ft_dft(1:N)/N; % cross power spectral density of error and ft
Suft = conj(u_dft(1:N)).*ft_dft(1:N)/N; %cross power spectral density of input and ft
Sef =  abs(Seft);
Suf =  abs(Suft);
%[pks,pklocs] = findpeaks(Sef(1:N/2),'MinPeakHeight',8);
Sef_pklocs = [9,18,37,56,72,99,139,188,261,308];
Suf_pklocs = [9,18,37,56,72,99,139,188,261,308];
figure(2)
loglog(Freq,Sef(1:N/2),Freq(Sef_pklocs),Sef(Sef_pklocs),'or')
legend('cpsd of e and ft','prominent peaks')
xlabel('frequecy[Hz]')
ylabel('|S_e_f_t(jw)|')
grid on
figure(3)
loglog(Freq,Suf(1:N/2),Freq(Sef_pklocs),Suf(Sef_pklocs),'or')
legend('cpsd of u and ft','prominent peaks')
xlabel('frequecy[Hz]')
ylabel('|S_e_f_t(jw)|')
grid on


%% Pilot dynamics
 Hp = Suft./Seft; %pilot dynamics at excitation freq
 magh = abs(Hp(Sef_pklocs)); % magnitude of Hp
 phaseh = 180*angle(Hp(Sef_pklocs))/pi; % phase of hp
 figure(4)
 loglog( Freq(Sef_pklocs), magh,'or') ;
 figure(5)
 semilogx( Freq(Sef_pklocs), phaseh,'or') ;
% 
% 
% 
%  