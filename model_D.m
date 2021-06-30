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
        %fun = @(x) 0 ;

        Phase = phase_M_1;
        Mag = mag_M_1;
        for k = 1 : length(omegaf)
            H_pe(k) = Mag(k) ;
            g = @(w) sum(abs(H_pe(k) - w(1)*1/(1 + w(2)*(1j*omega_test(l))) * exp(-1j*omegaf(k)*w(3)) * (w(5)^2/(w(5)^2 + 2*w(4)*w(5)*1j*omegaf(k) + (1j*omegaf(k))^2)))^2);
            data_phase1(i,k) = Phase(k) ;
        end
        %w0 =  [2,0.7,0.4,0.2,0.5,10]   ;  
        w0 =  [1.5,0.4,0.2,0.2,10]   ; 
        %[x,fval,exitflag,output] = fminsearch(g, x0, options);
        [w,fvalD] = fminsearch(g,w0,options);
        %lb = [0, 0, 0, 0, 5];
        %ub = [100, 10, 10, 1, 30];
        %x = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

        Kp_B(i,m) = w(1) ;
        %TL_B(i,m) = w(2) ;
        Ti_D(i,m) = w(2);
        Tp_B(i,m) = w(3) ;
        zeta_B(i,m) = w(4) ;
        omega_B(i,m) = w(5) ;

        phase_out = zeros(length(omega_test),1) ;
        mag_out = zeros(length(omega_test),1) ;
        
        get_i = 100000 ;
        get_ii = 100000 ;
        get_iii = 100000 ;
        count = 0 ;
        for l = 1:length(omega_test)
            Hpe_model = w(1)*1/(1 + w(2)*(1j*omega_test(l))) * exp(-1j*omega_test(l)*w(3)) * (w(5)^2/(w(5)^2 + 2*w(4)*w(5)*1j*omega_test(l) + (1j*omega_test(l))^2));
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
            data_magD(i,l) = mag_out(l) ;
            data_phaseD(i,l) = phase_out(l) ;
        end
        
    end
end

%modelD = sum(abs(H_pe(k) - w(1)*1/(1 + w(2)*(1j*omega_test(l))) * exp(-1j*omegaf(k)*w(3)) * (w(5)^2/(w(5)^2 + 2*w(4)*w(5)*1j*omegaf(k) + (1j*omegaf(k))^2))));
modelD =fvalD;
