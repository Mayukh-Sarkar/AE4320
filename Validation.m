load('ae4320_dataset1.mat')
N = 2^13;
t = t(1:N);
ft = ft(end-N+1:end);
error_v = validation.e;
input_v = validation.u;
error_v = error_v(end-N+1:end);
input_v = input_v(end-N+1:end);
f_t = ft;
dt = t(2) -t(1);
N = length(t);
Fs = 1/dt;

Freq = (Fs/N)*(1:N/2);
omega = 2*pi*Freq;

e_v_dft = fft(error_v,N); % Fourier transform of the error signal at excitation freq
u_v_dft = fft(input_v,N); % Fourier transform of the input signal at excitation freq
ft_dft = fft(ft,N);  % Fourier transform of the forcing function
[ft_pks , pklocs] = findpeaks(abs(ft_dft(1:N/2)),"MinPeakHeight", 0.1);

Seft_v = conj(e_v_dft(1:N)).*ft_dft(1:N)/N; % cross power spectral density of error and ft
Suft_v = conj(u_v_dft(1:N)).*ft_dft(1:N)/N; %cross power spectral density of input and ft
Sef_v =  abs(Seft_v);
Suf_v =  abs(Suft_v);



%% Pilot dynamics
omegaf = omega(pklocs)';
Hp_v = Suft_v(pklocs)./Seft_v(pklocs); %pilot dynamics at excitation freq
magh_v = abs(Hp_v); % magnitude of Hp
phaseh_v = 180*angle(Hp_v)/pi; % phase of hp


 



