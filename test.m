
% N_runs = 5 ;
N_participants = 1 ;
N_parameters = 5 ;
N_runs = 1 ;
human_par_M_Kp = zeros(N_runs,N_participants) ;
human_par_M_TL = zeros(N_runs,N_participants) ;
human_par_M_tau_p = zeros(N_runs,N_participants) ;
human_par_M_zeta_nm = zeros(N_runs,N_participants) ;
human_par_M_omega_nm = zeros(N_runs,N_participants) ;
omega_test = logspace(-1, 1.5, 200) ;
options = struct('MaxFunEvals', 15000,'MaxIter', 15000);
omega = omegaf;
data_mag1 = zeros(N_runs, length(omega)) ;
data_phase1 = zeros(N_runs, length(omega)) ;
data_mag2 = zeros(N_runs, length(omega_test)) ;
data_phase2 = zeros(N_runs, length(omega_test)) ;

for i = 1 : N_runs
    phase_M_1 = phaseh ;
    mag_M_1 = magh ;


    for m = N_participants
        H_pe = zeros(length(omega),1) ; 
        fun = @(x) 0 ;

        Phase = phaseh;
        Mag = magh;
        for k = 1 : length(omega)
            H_pe(k) = Mag(k) * (cos(Phase(k)*pi/180) + 1j * sin(Phase(k)*pi/180)) ;
            g = @(x) (abs(H_pe(k) - x(1)*(1 + x(2)*(1j*omega(k))) * exp(-1j*omega(k)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omega(k) + (1j*omega(k))^2))))^2;
            fun = @(x) fun(x) + g(x) ;
            data_mag1(i,k) = Mag(k) ;
            data_phase1(i,k) = Phase(k) ;
        end

        x0 = [3, 1, 0.35, 0.5,0] ;
        % [x,fval,exitflag,output] = fminsearch(f, x0, options);
%         x = fminsearch(f,x0,options);
        lb = [0, 0, 0, 0, 5];
        ub = [100, 10, 10, 1, 30];
        x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

        human_par_M_Kp(i,m) = x(1) ;
        human_par_M_TL(i,m) = x(2) ;
        human_par_M_tau_p(i,m) = x(3) ;
        human_par_M_zeta_nm(i,m) = x(4) ;
        human_par_M_omega_nm(i,m) = x(5) ;

        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;
        
        get_i = 100000 ;
        get_ii = 100000 ;
        get_iii = 100000 ;
        count = 0 ;
        for l = 1:length(omega)
            Hpe_model = x(1)*(1 + x(2)*(1j*omega_test(l))) * exp(-1j*omega_test(l)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
            mag_out(l) = abs(Hpe_model) ;
            phase_out(l) = angle(Hpe_model)*180/pi ;
            if (phase_out(l) < -150) && (count == 0)
                get_i = l ;
                count = 1 ;
            end
            if (l >= get_i) && (phase_out(l) > -150)
                phase_out(l) = phase_out(l) - 360 ;
            end
            if (phase_out(l) < -480) && (count == 1)
                get_ii = l ;
                count = 2 ;
            end
            if (l >= get_ii) && (phase_out(l) > -480)
                phase_out(l) = phase_out(l) - 360 ;
            end
%             if (phase_out(l) < -880) && (count == 2)
%                 get_iii = l ;
%                 count = 3 ;
%             end
%             if (l >= get_iii) && (phase_out(l) > -880)
%                 phase_out(l) = phase_out(l) - 360 ;
%             end
            data_mag2(i,l) = mag_out(l) ;
            data_phase2(i,l) = phase_out(l) ;
        end
        
%         title = "Bode plot for run #" + num2str(i) + " from participant #" + num2str(m) ;
% %         figure(N_participants*(i-1)+m)
%         figure(4+i)
%         subplot(2,1,1)
%         loglog(omega, mag_M(:,m), 'o')
%         hold on
%         loglog(omega_test, mag_out)
%         hold off
%         axis(10.^[-1 1.5 -1 2])
%         legend('Experiment','Parameter estimation','Location','southwest')
% 
%         subplot(2,1,2)
%         semilogx(omega, phase_M(:,m), 'o')
%         hold on
%         semilogx(omega_test, phase_out)
%         hold off
%         axis([10.^-1 10.^1.5 -360 180])
%         legend('Experiment','Parameter estimation','Location','southwest')
%         sgtitle(title)
    end
end

%% Bode plots M
% title = "Bode plot for run 1 and 60 from participant 1 in Moving simulator" ;
figure(1)
set(gcf, 'Position', [100 100 700 650])
subplot(2,1,1)
loglog(omega, data_mag1(1,:),'ok')
hold off
axis(10.^[-1 1.5 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
ah=gca; 
set(ah,'Fontsize',12)
% legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southwest')

subplot(2,1,2)
semilogx(omega, data_phase1(1,:),'ok')
hold on
semilogx(omega_test, data_phase2(1,:),'-k')
hold off
axis([10.^-1 10.^1.5 -360 180])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southwest')
ah=gca; 
set(ah,'Fontsize',12)
% sgtitle(title)