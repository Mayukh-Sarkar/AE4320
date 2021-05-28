for i = 1 : 1
    phase_M_1 = phaseh ;


    mag_M_1 = magh ;


    for m = 1 : 1
        H_pe = zeros(length(omegaf),1) ; 
        fun = @(y) 0 ;

        Phase = phase_M_1;
        Mag = mag_M_1;
        for k = 1 : length(omegaf)
            H_pe(k) = Mag(k) ;
            g = @(y) (abs(H_pe(k) - y(1) * exp(-1j*omegaf(k)*y(2)) * (y(4)^2/(y(4)^2 + 2*y(4)*y(3)*1j*omegaf(k) + (1j*omegaf(k))^2))))^2;
            data_mag1(i,k) = Mag(k) ;
            data_phase1(i,k) = Phase(k) ;
        end

        %x0 = [3,0.35,0.5,1,15] ;
        x0 = [2.5,0.15,0.25,15];
        %x0=[1.5,0.3,0.1,1,10];
        %[x,fval,exitflag,output] = fminsearch(g, x0, options);
        %x = fminsearch(g,x0,options);
        lb = [0, 0, 0, 5];
        ub = [100, 10, 1, 30];
        y = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

        Kp_B(i,m) = y(1) ;
        
        Tp_B(i,m) = y(2) ;
        zeta_B(i,m) = y(3) ;
        omega_B(i,m) = y(4) ;

        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;
        
        get_i = 100000 ;
        get_ii = 100000 ;
        get_iii = 100000 ;
        count = 0 ;
        for l = 1:length(omega_test)
            Hpe_model = y(1) * exp(-1j*omega_test(l)*y(2)) * (y(4)^2/(y(4)^2 + 2*y(4)*y(3)*1j*omega_test(l) + (1j*omega_test(l))^2));
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
            data_magA(i,l) = mag_out(l) ;
            data_phaseA(i,l) = phase_out(l) ;
        end
        
    end
end

%% Bode plots M
% title = "Bode plot for run 1 and 60 from participant 1 in Moving simulator" ;
figure(2)
set(gcf, 'Position', [100 100 700 650])
subplot(2,1,1)
loglog(omegaf,magh,'ok')
hold on
loglog(omega_test, data_magA(1,:),'-k')
% hold on
% plot(omega, data_mag1(60,:),'*k')
% hold on
% loglog(omega_test, data_mag2(60,:),'--k')
hold off
axis(10.^[-0.5 1.3 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
ah=gca; 
set(ah,'Fontsize',12)
% legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southomegafest')

subplot(2,1,2)
semilogx(omegaf,phaseh,'ok')
hold on
semilogx(omega_test, data_phaseA(1,:),'-k')
% hold on
% semilogx(omega, data_phase1(60,:),'*k')
% hold on
% semilogx(omega_test, data_phase2(60,:),'--k')
hold off
axis([10.^-0.5 10.^1.3 -180 180])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
%legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southomegafest')
ah=gca; 
set(ah,'Fontsize',12)
% sgtitle(title)
