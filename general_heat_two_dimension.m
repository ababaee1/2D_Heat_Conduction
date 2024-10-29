function [TT,tt,nt,r,z,kappa,cp,ro, E, alpha_thermal, nu,index_heated_top,all_index_top]=general_heat_two_dimension_FGM_GDQ_TID(nr,nz,dt,a,b,hh,kisi,right_edge,top_edge,bottom_edge,left_edge,partial_ratio_top,partial_ratio_bottom)
%% material propeties 
format long
clc;
t=293;

% 'insul' --> zero heat flux 
% 'convection' --> heat convection with air
% 'prescribed' --> prescribed temprature
% right_edge='prescribed';
% top_edge='prescribed';
% bottom_edge='prescribed';
% left_edge='insul';

[c0z,c1z,c2z,z]= GDQ_heat_1D(nz,hh);
[c0r,c1r,c2r,r]= GDQ_1D(nr,a,b);

% kappa SUS 304
 s1=-1.4087; % zaribe sabet
 s2=1.3982;
 s3=0.2543;
 s4=-0.6260;
 s5=0.2334;
 s6=0.4256;
 s7=-0.4658;
 s8=0.1650;
 s9=-0.0199;

% kappa AISI 1020
 p1 =181.41306;  %zarib sabet
 p2 =- 544.33993; 
 p3 =  675.78083;
 p4 = - 441.94194 ;
 p5 = 161.06483; 
 p6 = - 31.07912;
 p7 = 2.48398; 
 p8 = 0; 
 p9 = 0; 


%cp SUS 304 

aa=	22.006100;
bb=	-127.552800;
cc=	303.647000;
dd=	-381.009800;
e=	274.032800;
f=	-112.921200;
g=	24.759300;
h=	-2.239153;
i=	0.000000;


%cp AISI 1020

a11=	- 220.09324;
b11=	+ 628.28238;
c11=	- 740.82380;
d11=	+ 466.14910;
e11=	- 164.95367;
f11=	+31.12036;
g11=	-2.44526;
h11=	0;
i11=	0;


%CTE for SUS 304
A=	-2.95540000E+02;
B=	-3.98110000E-01;
C=	9.26830000E-03;
D=	-2.02610000E-05;
EE=	1.71270000E-08;

%CTE for AISI 1020

pp1 = +8.034e-09   ;	
pp2 =-1.108e-05 ;	
pp3 =+0.006069 ; 	
pp4 =-0.36 ; 	
pp5 =-196.1;  	
       
% E for SUS 304
aa_1=	2.1005930E+02;
bb_1=	1.5348830E-01;
cc_1=	-1.6173900E-03;
dd_1=	5.1170600E-06;
ee_1=	-6.1546000E-09;

% E for AISI 1020 

pp11 =1.4e-09	;
pp22 =-1.36e-06	;
pp33 =0.0004651 ;	
pp44 =-0.1125	;
pp55 =221.3;	

% Nu for SUS 304 and AISI 1020

PP1 =1.664e-12;  	
PP2 =-1.462e-09; 	
PP3 =4.721e-07 	;
PP4 = -1.043e-05 ; 	
PP5 = 0.2766 	;


%50<T(K)<300 (for CTE)
 alpha_t=(10^-5)*(B+ 2*C*t +3*D*(t.^2) +4*EE*(t.^3));%1/k
 alpha_b=(10^-5)*(pp4 +2*pp3*t +3*pp2*(t.^2) +4*pp1*(t.^3)); %1/k			

 %57<T(K)<300 (for E)
 
 Et=10^9*(aa_1 +(bb_1*t) +cc_1*(t.^2) +dd_1*(t.^3) +ee_1*(t.^4)); %pa
 Eb=10^9*(pp11*t.^4 + pp22*t.^3 + pp33*t.^2 + pp44*t +pp55); %Pa 
 
 %57<T(K)<300 (for E)
 nu_t= PP1*(t.^4) + PP2*(t.^3) + PP3*(t.^2) + PP4*t + PP5 ;				
 nu_b=PP1*(t.^4) + PP2*(t.^3) + PP3*(t.^2) + PP4*t + PP5 ;


% Temperature dependency for properties 

 %50<T(K)<300 (for conductivity)
 kappa_t= 10.^(s9*((log10(t)).^8)+ s8*((log10(t)).^7)+ s7*((log10(t)).^6) +s6*((log10(t)).^5) + s5*((log10(t)).^4) +s4*((log10(t)).^3) + s3*((log10(t)).^2) +s2*((log10(t))) + s1);
 kappa_b=10.^(p9*((log10(t)).^8)+ p8*((log10(t)).^7)+ p7*((log10(t)).^6) +p6*((log10(t)).^5) + p5*((log10(t)).^4) +p4*((log10(t)).^3) + p3*((log10(t)).^2) +p2*((log10(t))) + p1);
 
 %50<T(K)<300 (for specific heat) 
 cp_t=10.^(i*((log10(t)).^8)+ h*((log10(t)).^7)+ g*((log10(t)).^6) +f*((log10(t)).^5) + e*((log10(t)).^4) +dd*((log10(t)).^3) + cc*((log10(t)).^2) +bb*((log10(t))) + aa);
 cp_b=10.^(i11*((log10(t)).^8)+ h11*((log10(t)).^7)+ g11*((log10(t)).^6) +f11*((log10(t)).^5) + e11*((log10(t)).^4) +d11*((log10(t)).^3) + c11*((log10(t)).^2) +b11*((log10(t))) + a11);
 
 ro_t=7850; %kg/m^3
 ro_b=7850; %kg/m^3
    
 cp_cm=cp_t-cp_b;
 kappa_cm=kappa_t-kappa_b;
 ro_cm=ro_t-ro_b;
 
 cp_old=cp_b +cp_cm.*(0.5 +z'./hh).^kisi ;
 kappa_old=kappa_b +kappa_cm.*(0.5 +z'./hh).^kisi; 
 ro_old=ro_b +ro_cm*(0.5 +z'./hh).^kisi ;
 

 Ecm=Et-Eb;
 alpha_cm=alpha_t-alpha_b;
 nu_cm=nu_t-nu_b;  

 
 alpha_old=alpha_b +alpha_cm.*(0.5 +z'./hh).^kisi ;
 E_old=Eb +Ecm.*(0.5 +z'./hh).^kisi ;
 nu_old=nu_b +nu_cm.*(0.5 +z'./hh).^kisi ;
 
 
 %%
tf=1; %final time
t0=0; %initial time
nt=[(tf-t0)/dt]+1; % to compute 10 time steps k=[1,11]

%dz=(z/(n-1)); %dz is space domain descritization step 

k=zeros(nr*nz,nr*nz,nt);
k_bar=zeros(nr*nz,nr*nz,nt);
k_hat=zeros(nr*nz,nr*nz,nt);
k_r=zeros(nr*nz,nr*nz,nt);  %stiffness for radial derivatives 
k_z=zeros(nr*nz,nr*nz,nt);  % stiffness for z derivatives 
c=zeros(nr*nz,nr*nz,nt);
UUU=zeros(nr*nz,1,nt);
F=zeros(nr*nz,1,nt);  % F matrix RHS
F_hat=zeros(nr*nz,1,nt);  % F_hat  matrix RHS

conv_coef=5;  % W/m^2. K   Heat transfer coeficient of air in convection  
T_inf= 293; %%  temperature of far field in air 


alpha=0.5; % crank-Nicolson scheme parameter (always stable) 
a1=alpha*dt; 
a2=(1-alpha)*dt;

UUU(:,1,1)=293; %initial condition value for temprature T=300K 


% changing the properties from vectors (1D) to 2D as eveything is 2D now ! 

% here the FGM is only 1D through the thickness (z)
% here the FGM is only 1D through the thickness (z)
for j=1:nz
    ro(j,1 : nr)=ro_old(j);
    cp(j,1 : nr)=cp_old(j);
    kappa(j,1 : nr)=kappa_old(j);
    alpha_thermal(j,1:nr)=alpha_old(j);
    E(j,1:nr)=E_old(j);
    nu(j,1:nr)=nu_old(j);
end


%form C matrix which is -Ro*Cp
for i=1:nr*nz

    c(i,i,:)=-ro(i)*cp(i);
  
end


% heated_top_index
index_top=find(partial_ratio_top*a>=r);
index_heated_top= 1:nz:(index_top(end)-1)*nz +1;

% non-heated_top_index
all_index_top= 1:nz:(nr-1)*nz +1 ; 
index_non_heated_top=setdiff(all_index_top,index_heated_top);



% solution of Parabolic Eq with GDQ and Crank Nicolson 
for i=1:nt
    
  
    % k([1+(j-1)*nr : j*nr], [1+(j-1)*nr : j*nr],i) =  (1/r(j)).*(kappa(:,j).*(c2r(j,:) + (r(j)).*(c2r(j,:)*c0z)))  +(c0r(j,:)*c1z*kappa(:,j)).*(c0r(j,:)*c1z) + kappa(:,j).*(c0r(j,:)*c2z);
    
     C_0 = repmat({c0z},nr,1);
    k_z_0 = blkdiag(C_0{:});
    
    C_1 = repmat({c1z},nr,1);
    k_z_1 = blkdiag(C_1{:});
    
    C_2 = repmat({c2z},nr,1);
    k_z_2 = blkdiag(C_2{:});
% %     
   RR=kron(r', ones(nr,1));
 %   RR=repelem(r',nr,1);
    EYE=eye(nr);
% %     
     k_r(:,:,i)=(kappa(:)./RR).*(kron(c1r,EYE)) + kappa(:).*(kron(c2r,EYE));   
     k_z(:,:,i)= (k_z_1*kappa(:)).*(k_z_1) + kappa(:).*(k_z_2);
% %     
     k(:,:,i)=k_z(:,:,i)+k_r(:,:,i) ;
%k_r(:,:,i) +
%
% kk_new(:,:,i)=double(k_new);
    % Initial Boundary condition
    %%
          
    if i==1  % Initial Boundary condition for whole nodes 
     c(:,:,i)=0;
     for j=1:nr*nz
      k(j,:,i)=0; 
      k(j,j,i)=1 ;
%       kk_new(j,:,i)=0; 
%       kk_new(j,j,i)=1 ;
      F(j,1,i)=293;  % Initial Boundary condition fow whole nodes 
     end 
    end
    
    if i>=2
        % Boundary conditions for 2D heat conduction in r and z directions 
     
 
        k_hat(:,:,i)=c(:,:,i) +a1*k(:,:,i);
        
%----------------------------------------------------------------------------        
        % MAIN BC for k_hat
 %----------------------------------------------------------------------------

        Conv_r= kappa(:).*(kron(c1r,EYE)) +conv_coef.*kron(c0r,EYE);
        inusl_r=kron(c1r,EYE); 
        inusl_z=k_z_1;
        % @ top face for convection with air 
        % given that the z is positive upward --> -kdT/dz -h(T-T_inf)=0
        Conv_z_top= kappa(:).*(k_z_1) +conv_coef.*k_z_0;
        
        % @ bottom face for convection with air 
        % given that the z is positive upward --> +kdT/dz -h(T-T_inf)=0
        Conv_z_bottom= kappa(:).*(k_z_1) -conv_coef.*k_z_0;
        
  
        % RIGHT edge 
        for jj= (nr-1)*nz +1:1:nr*nz
            if right_edge=="prescribed"
        k_hat(jj,:,i)=0; 
        k_hat(jj,jj,i)=1;
            elseif right_edge=="conv"
        k_hat(jj,:,i)=Conv_r(jj,:);
            elseif right_edge=="insul"
        k_hat(jj,:,i)=inusl_r(jj,:);
            end
        end
        
        % LEFT edge 
        for jj= 1:1:nz
            if left_edge=="prescribed"
        k_hat(jj,:,i)=0; 
        k_hat(jj,jj,i)=1;
            elseif left_edge=="insul"
        k_hat(jj,:,i)=inusl_r(jj,:); 
            end
        end

  
        % Bottom edge  
        index_bottom=find(partial_ratio_bottom*a>=r); 
        for jj= nz:nz:(index_bottom(end))*nz 
            if bottom_edge=="prescribed"
        k_hat(jj,:,i)=0; 
        k_hat(jj,jj,i)=1;
            elseif bottom_edge=="conv"
        k_hat(jj,:,i)=Conv_z_bottom(jj,:);
            elseif bottom_edge=="insul"
        k_hat(jj,:,i)=inusl_z(jj,:);
            end
        end

        
        % TOP edge 
        % Heated part 
        for jj= index_heated_top
            k_hat(jj,:,i)=0; 
            k_hat(jj,jj,i)=1;  %top face prescribed value (heated part)
        end
        
        % Non-heated part
        for jj=index_non_heated_top
            if top_edge=="prescribed" 
            k_hat(jj,:,i)=0; 
            k_hat(jj,jj,i)=1;  %top face prescribed value (non-heated part)
            elseif top_edge=="conv"
            k_hat(jj,:,i)=Conv_z_top(jj,:);    %top face heat convection with air  (non- heated part)     
            elseif top_edge=="insul"
            k_hat(jj,:,i)=inusl_z(jj,:);
            end
        end
       
%%  -----------------------------------------------------------------------
%        k_bar(:,:,i-1)=c(:,:,i) -a2*c(:,:,i)*inv(c(:,:,i-1))*k(:,:,i-1);  % changes in the second term 
%         F_hat(:,1,i)=k_bar(:,:,i-1)*UUU(:,1,i-1) +a1*F(:,1,i) +a2*c(:,:,i)*inv(c(:,:,i-1))*F(:,1,i-1);  % changes in the last term
       
       k_bar(:,:,i-1)=c(:,:,i) -a2*k(:,:,i-1);          
       F_hat(:,1,i)=k_bar(:,:,i-1)*UUU(:,1,i-1) +a1*F(:,1,i) +a2*F(:,1,i-1);
%%-------------------------------------------------------------------------        
        %MAIN BC for F_hat
%----------------------------------------------------------------------------

         % RIGHT edge         
        for jj= (nr-1)*nz +1:1:nr*nz
            if right_edge=="prescribed"
        F_hat(jj,1,i)=293;  %right edge prescribed value 
            elseif right_edge=="conv"
        F_hat(jj,1,i)=conv_coef*T_inf;
            elseif right_edge=="insul"
        F_hat(jj,1,i)=0;
            end
        end

        % LEFT edge 
        
        for jj= 1:1:nz
            if left_edge=="prescribed"
        F_hat(jj,1,i)=293;
            elseif left_edge=="insul"
        F_hat(jj,1,i)=0;  %left edge zero heat flux in radial dT/dr=0 @r=0  
            end
        end


        % Bottom edge  
        index_bottom=find(partial_ratio_bottom*a>=r); 
        for jj= nz:nz:(index_bottom(end))*nz  
            if bottom_edge=="prescribed"
        F_hat(jj,1,i)=293;  %lower face prescribed value  
        elseif bottom_edge=="conv"
        F_hat(jj,1,i)=-conv_coef*T_inf;     % check the sign !! 
            elseif bottom_edge=="insul"
        F_hat(jj,1,i)=0;        
            end
        end 
    

         % TOP edge  heated part 
        for jj= index_heated_top
            F_hat(jj,1,i)=93;  %top face prescribed value (heated part)
        end
        
        
        % Top face non-heated part 
         for jj=index_non_heated_top
             if top_edge=="prescribed" 
           F_hat(jj,1,i)=293;           %top face prescribed value (non-heated part)
              elseif top_edge=="conv"
           F_hat(jj,1,i)=conv_coef*T_inf;  %top face air heat convection (non- heated part)    
               elseif top_edge=="insul"
            F_hat(jj,1,i)=0;        
             end
         end
       
        %%
        UUU(:,1,i)=inv(k_hat(:,:,i))*F_hat(:,1,i); % Temprature distrubution   
        
        tt(i)=t0+(i-1)*dt;
    end
    
    
end

%

% figure(5)
% hold on
% for i=1:length(tt)
% y(1,i)=UUU(nr*ceil(nr/2) - (ceil(nr/2)-1),1,i);
% end
% 
% plot(tt,y)
% 
figure(3)
% hold off
TT = reshape(UUU(:,1,:),[nz,nr,nt]);
[R,Z]=meshgrid(r,z);
surf(R,Z,TT(:,:,end))
view(2)
end
