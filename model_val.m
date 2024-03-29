% cost for validation data set
for i = 1 : 1
    phase_M_1 = phaseh ;


    mag_M_1 = magh_v ;
    
   


    for m = 1 : 1
        H_pe = zeros(length(omegaf),1) ; 
        
        %fun = @(y) 0 ;

        Phase = phase_M_1;
        Mag = mag_M_1;
    %% cost function for the validation data set
        for k = 1 : length(omegaf)
            H_pe(k) = Mag(k) ;
            
            
            Av = @(j)sum(abs(H_pe(k) - j(1) * exp(-1j*omegaf(k)*j(2)) * (j(4)^2/(j(4)^2 + 2*j(4)*j(3)*1j*omegaf(k) + (1j*omegaf(k))^2)))^2);
            Bv = @(xv) sum(abs(H_pe(k) - xv(1)*(1 + xv(2)*(1j*omegaf(k))) * exp(-1j*omegaf(k)*xv(3)) * (xv(5)^2/(xv(5)^2 + 2*xv(4)*xv(5)*1j*omegaf(k) + (1j*omegaf(k))^2)))^2);
            Cv = @(zv)  sum(abs(H_pe(k) - zv(1)*(1 + zv(2)*(1j*omegaf(k))) /(1+zv(3)*(1j*omegaf(k)))* exp(-1j*omegaf(k)*zv(4)) * (zv(6)^2/(zv(6)^2 + 2*zv(5)*zv(6)*1j*omegaf(k) + (1j*omegaf(k))^2)))^2);
            Dv = @(wv) sum(abs(H_pe(k) - wv(1)*1/(1 + wv(2)*(1j*omega_test(l))) * exp(-1j*omegaf(k)*wv(3)) * (wv(5)^2/(wv(5)^2 + 2*wv(4)*wv(5)*1j*omegaf(k) + (1j*omegaf(k))^2)))^2);



            %             data_mag1(i,k) = Mag(k) ;
%             data_phase1(i,k) = Phase(k) ;
        end
       j0 = [2.5,0.25,0.2,15];
       [j,fvalAv] = fminsearch(Av,j0,options);
       xv0 = [2.7,0.3,0.2,0.2,10] ;  
       [xv,fvalBv] = fminsearch(Bv,xv0,options);
       zv0 = [2,0.7,0.4,0.2,0.5,10] ;  
       [zv,fvalCv] = fminsearch(Cv,zv0,options);
       wv0 =  [1.5,0.4,0.2,0.2,10]   ; 
       [wv,fvalDv] = fminsearch(Dv,wv0,options);
  




       
        
    end
end
%% cost function values for all the models using validation data set
modelAv = fvalAv;
modelBv = fvalBv;
modelCv = fvalCv;
modelDv = fvalDv;
