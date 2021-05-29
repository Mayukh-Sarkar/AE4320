%%%% Specify function for optimization

omega_test = logspace(-0.5, 1.5, 100);

options = struct('MaxFunEvals', 15000,'MaxIter', 15000);

data_mag1 = zeros(1, length(omega)) ;
data_phase1 = zeros(1, length(omega)) ;

phase_M_1 = phaseh ;
mag_M_1 = magh ;
H_pe = zeros(length(omegaf),1) ; 
Phase = phase_M_1;
Mag = mag_M_1;


for i = 1 : 1
    phase_M_1 = phaseh ;


    mag_M_1 = magh ;


    for m = 1 : 1
        H_pe = zeros(length(omegaf),1) ; 
        fun = @(x) 0 ;

        Phase = phase_M_1;
        Mag = mag_M_1;
        for k = 1 : length(omegaf)
            H_pe(k) = Mag(k) ;
            g = @(x) (abs(H_pe(k) - x(1)*(1 + x(2)*(1j*omegaf(k))) * exp(-1j*omegaf(k)*x(3)) * (x(5)^2/(x(5)^2 + 2*x(4)*x(5)*1j*omegaf(k) + (1j*omegaf(k))^2))))^2;
            data_mag1(i,k) = Mag(k) ;
            data_phase1(i,k) = Phase(k) ;
        end

        %x0 = [3,0.35,0.5,1,15] ;
        x0 = [2.5,0.4,0.1,1,15];
        %x0=[1.5,0.3,0.1,1,10];
        %[x,fval,exitflag,output] = fminsearch(g, x0, options);
        %x = fminsearch(g,x0,options);
        lb = [0, 0, 0, 0, 5];
        ub = [100, 10, 10, 1, 30];
        x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

        Kp_B(i,m) = x(1) ;
        TL_A(i,m) = x(2) ;
        Tp_A(i,m) = x(3) ;
        zeta_A(i,m) = x(4) ;
        omega_A(i,m) = x(5) ;

        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;
        
        get_i = 100000 ;
        get_ii = 100000 ;
        get_iii = 100000 ;
        count = 0 ;
        for l = 1:length(omega_test)
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
            data_magB(i,l) = mag_out(l) ;
            data_phaseB(i,l) = phase_out(l) ;
        end
        
    end
end


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
        y0 = [2.5,0.15,0.25,15];
        %x0=[1.5,0.3,0.1,1,10];
        %[x,fval,exitflag,output] = fminsearch(g, x0, options);
        %x = fminsearch(g,x0,options);
        lb = [0, 0, 0, 5];
        ub = [100, 10, 1, 30];
        y = fmincon(fun,y0,[],[],[],[],lb,ub,[],options);

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

for i = 1 : 1
    phase_M_1 = phaseh ;


    mag_M_1 = magh ;


    for m = 1 : 1
        H_pe = zeros(length(omegaf),1) ; 
        fun = @(z) 0 ;

        Phase = phase_M_1;
        Mag = mag_M_1;
        for k = 1 : length(omegaf)
            H_pe(k) = Mag(k) ;
            g = @(z) (abs(H_pe(k) - z(1)*(1 + z(2)*(1j*omegaf(k))) /(1+z(3)*(1j*omegaf(k)))* exp(-1j*omegaf(k)*z(4)) * (z(6)^2/(z(6)^2 + 2*z(5)*z(6)*1j*omegaf(k) + (1j*omegaf(k))^2))))^2;
            data_magC(i,k) = Mag(k) ;
            data_phaseC(i,k) = Phase(k) ;
        end

        z0 = [2.5,0.7,0.4,0.2,1,15] ;
        %[x,fval,exitflag,output] = fminsearch(g, x0, options);
        %x = fminsearch(g,x0,options);
        lb = [0, 0, 0, 0,0, 5];
        ub = [100, 10, 10,10, 1, 30];
        z = fmincon(fun,z0,[],[],[],[],lb,ub,[],options);

        Kp(i,m) = z(1) ;
        TL(i,m) = z(2) ;
       Ti(i,m) = z(3) ;
        tau_p(i,m) = z(4) ;
        zeta(i,m) = z(5) ;
        omega(i,m) = z(6) ;

        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;
        
        get_i = 100000 ;
        get_ii = 100000 ;
        get_iii = 100000 ;
        count = 0 ;
        for l = 1:length(omega_test)
            Hpe_model = z(1)*(1 + z(2)*(1j*omega_test(l)))/(1+z(3)*(1j*omega_test(k))) * exp(-1j*omega_test(l)*z(4)) * (z(6)^2/(z(6)^2 + 2*z(6)*z(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
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
            data_magC(i,l) = mag_out(l) ;
            data_phaseC(i,l) = phase_out(l) ;
        end
        
    end
end


val

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
hold off
axis(10.^[-1 1.5 -1 2])
xlabel("\omega [rad/s]")
ylabel("|H_{p} (j \omega)| [-]")
ah=gca; 
set(ah,'Fontsize',12)


subplot(2,1,2)
semilogx(omega(pklocs), phaseh,'ok')
hold on
semilogx(omega_test, data_phaseA(1,:),'-k')
hold on
semilogx(omega_test, data_phaseB(1,:),'--k')
hold on
semilogx(omega_test, data_phaseC(1,:),'-.k')
hold off
axis([10.^-1 10.^1.5 -360 180])
xlabel("\omega [rad/s]")
ylabel("\angle H_{p} (j \omega) [deg]")
%legend('Test data run 1','Test data run 60','Model run 1','Model run 60','Location','southomegafest')
ah=gca; 
set(ah,'Fontsize',12)