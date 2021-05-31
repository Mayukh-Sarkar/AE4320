for i = 1 : 1
    phase_M_1 = phaseh ;


    mag_M_1 = magh ;


    for m = 1 : 1
        H_pe = zeros(length(omegaf),1) ; 
        %fun = @(y) 0 ;

        Phase = phase_M_1;
        Mag = mag_M_1;
        for k = 1 : length(omegaf)
            H_pe(k) = Mag(k) ;
            h = @(y) (abs(H_pe(k) - y(1) * exp(-1j*omegaf(k)*y(2)) * (y(4)^2/(y(4)^2 + 2*y(4)*y(3)*1j*omegaf(k) + (1j*omegaf(k))^2))))^2;
            data_mag1(i,k) = Mag(k) ;
            data_phase1(i,k) = Phase(k) ;
        end

        %y0 = [3,0.35,0.5,1,15] ;
        y0 = [2.5,0.15,0.25,15]; %cost    2.0166e-09
        %y0 = [0,0,0,0];
        %y0=[1.5,0.3,0.1,1,10];
        %[x,fval,exitflag,output] = fminsearch(g, x0, options);
        y = fminsearch(h,y0,options);
        %lb = [0, 0, 0, 5];
        %ub = [100, 10, 1, 30];
        %y = fmincon(h,y0,[],[],[],[],lb,ub,[],options);

        Kp_A(i,m) = y(1) ;
        
        Tp_A(i,m) = y(2) ;
        zeta_A(i,m) = y(3) ;
        omega_A(i,m) = y(4) ;

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

