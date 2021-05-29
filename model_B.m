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
        TL_B(i,m) = x(2) ;
        Tp_B(i,m) = x(3) ;
        zeta_B(i,m) = x(4) ;
        omega_B(i,m) = x(5) ;

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


