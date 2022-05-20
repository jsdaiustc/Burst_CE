function [H]=Burst_CE(Y,X,N,N_point,etc)

% N: number of antennas
% N_point: number of grid points
% etc:  number of active grid points

% Note that the ccurrently used grid  uniformly covers the range [-1,1], which is different from the definition ([-90,90]) in Section III-A of TSP2018_Dai . 
% The two forms are equvalent. But the first one can connect with the DFT basis directly (with N=N_point).

[T,M]=size(Y);
search_area=[-1:2/N_point:1];
reslu=2/N_point;
F=exp(-1i*pi*(0:N-1)'*search_area)/sqrt(N);
a=1e-10;b=1e-10; 
maxiter=500;
tol=1e-5;

%initialization
converged = false;
iter = 0;
Aw=X*F;
alpha0=1;
delta_last=100*ones(length(search_area),1);
delta=ones(length(search_area),1);
delta_inv=1./delta;
Z=ones(length(delta),3);
Z= diag(  1./sum(Z,2) ) *   Z;
deltal=[delta(2:length(delta));delta(1)];
deltar=[delta(end);delta(1:length(delta)-1)];
diag_inv=diag( 1./(  Z(:,1).*deltar + Z(:,2).*delta +  Z(:,3).*deltal   ) );
Phi_delta = Aw *  diag_inv;
V_temp= 1/alpha0*eye(T) + Phi_delta * Aw';
Sigma = diag_inv  -Phi_delta' * (V_temp \Phi_delta);
mu = alpha0 * Sigma * Aw' * Y;
Phi=Aw;  

while ~converged
  
   switch iter-floor(iter/5)*5
     
       case 0
               %update alpha
               resid=Y-Phi*mu;
               term2=sum(diag(Phi* Sigma*Phi'));
               alpha0=( M*T + a )/( b +  norm(resid(:), 'fro')^2  +   real(term2) );
        case 1
               %calculate mu and Sigma
               deltal=[delta(2:length(delta));delta(1)];
               deltar=[delta(end);delta(1:length(delta)-1)];
               diag_inv=diag( 1./(  Z(:,1).*deltar + Z(:,2).*delta +  Z(:,3).*deltal   ) );
               Phi_delta = Phi *  diag_inv;
               V_temp= 1/alpha0*eye(T) + Phi_delta * Phi';
               Sigma = diag_inv  -Phi_delta' * (V_temp \Phi_delta);
               mu = alpha0 * Sigma * Phi' * Y;
       case 2
               %update delta
               delta_last = delta_inv;
               sum_mu=sum( mu.*conj(mu), 2);
               temp=sum_mu + M*real(diag(Sigma));    
               tc=real(temp);
               tcl=[tc(2:length(tc));tc(1)];
               tcr=[tc(end);tc(1:length(tc)-1)];
               Z_2=Z(:,2);
               Z_1l=[Z(2:length(delta),1);Z(1,1)];
               Z_3r=[Z(end,3);Z(1:length(delta)-1,3)];
               c_k = a + Z_1l + Z_2 + Z_3r;
               d_k = b + Z_1l.*tcl + Z_2.*tc + Z_3r.*tcr;
               delta=c_k./d_k;
               ln_delta=   psi( c_k  ) -  log(  d_k  );           
               delta_inv=1./delta;

       case 3
               %update Z
                ln_deltal=[ln_delta(2:length(ln_delta));ln_delta(1)];
                ln_deltar=[ln_delta(end);ln_delta(1:length(ln_delta)-1)];
                deltal=[delta(2:length(delta));delta(1)];
                deltar=[delta(end);delta(1:length(delta)-1)];
                t1 = ln_deltar - tc.*deltar ;
                t2 = ln_delta  - tc.*delta ;
                t3 = ln_deltal - tc.*deltal ;
                et=[t1,t2,t3];
                for kk=1:length(delta)
                     temp_p(kk,:)= exp(et(kk,:)-max(max(et(kk,:)))); 
                end  
                Z= diag(  1./sum(temp_p,2) ) *   temp_p;
       case 4
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% grid refine
                sum_mu=sum( mu.*conj(mu), 2);
                resid=Y-Phi*mu;
                Pm=sum_mu;
                [~,sort_ind]=sort(Pm, 'descend');    
                index_amp = sort_ind(1:etc);
                tempPS=Phi*  Sigma(:,index_amp) ;       
                for j=1:length(index_amp)
                    ii=index_amp(j);
                    ai=F(:,ii);
                    mut=mu(ii,:);
                    Sigmat=Sigma(:,ii);
                    c1=mut*mut' +  M* Sigmat(ii);
                    c1=abs(c1)*(-alpha0);
                    Yti=resid +  Phi(:, ii)*mu(ii,:);
                    c2=  M*(  tempPS(:,j) - Phi(:,ii)*Sigmat(ii) )  -Yti*(mut');
                    c2= c2*(-alpha0);
                     phii=Phi(:,ii);
                    sinta= search_area(ii); costa=cos(asin(sinta));
                    c3=(-1i*pi* costa)*[0:N-1]';    %c3=(-1i*2*pi*gap_d/sqrt(N))*[0:N-1]';
                    tt1=  X*(c3.*ai); 
                    f1= tt1'*phii*c1  +   tt1'*c2;
                    f1= 2*real(f1);
                    angle_cand= sign(f1)*reslu/100 ;       
                    sin_add = search_area(ii) + angle_cand;
                    search_area(ii)=sin_add;
                    ro1=exp(-sin_add*pi*1i);            
                    F(:,ii)=ro1.^((0:N-1)')/sqrt(N);
                    Aw(:,ii)=X* F(:,ii);
                end    
                Phi=Aw;  
   end
        % stopping criteria
        erro=norm(delta_inv - delta_last)/norm(delta_last);
        if erro < tol || iter >= maxiter
            converged = true;
        end
         iter = iter + 1;
end

H=F*mu;










