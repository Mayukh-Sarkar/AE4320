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
            f = @(z) (abs(H_pe(k) - z(1)*(1 + z(2)*(1j*omegaf(k))) /(1+z(3)*(1j*omegaf(k)))* exp(-1j*omegaf(k)*z(4)) * (z(6)^2/(z(6)^2 + 2*z(5)*z(6)*1j*omegaf(k) + (1j*omegaf(k))^2))))^2;
            data_magC(i,k) = Mag(k) ;
            data_phaseC(i,k) = Phase(k) ;
        end
        
        %z0 = [2.5,0.35,0.1,0.2,0.5,15] ; % cost 3.0141e-10
        %z0 = [3,0.35,0.1,0.2,1,15] ; %cost power -9
        z0 = [2,0.7,0.4,0.2,0.5,10] ;  %cost 3.9903e-10
        %z0 = [2,0.7,0.3,0.2,0.5,10] ; %power in power -8
        %[x,fval,exitflag,output] = fminsearch(g, x0, options);
        z = fminsearch(f,z0,options);
        %z0 = [0,0,0,0,0,0];
        lb = [0, 0, 0, 0,0, 5];
        ub = [100, 10, 10,10, 1, 30];
        %z = fmincon(g,z0,[],[],[],[],lb,ub,[],options);

        Kp_C(i,m) = z(1) ;
        TL_C(i,m) = z(2) ;
       Ti_C(i,m) = z(3) ;
        tau_C(i,m) = z(4) ;
        zeta_C(i,m) = z(5) ;
        omega_C(i,m) = z(6) ;

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


