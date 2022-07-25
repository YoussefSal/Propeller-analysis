function CP_1()

clearvars
clc
close all

fprintf('Welcome, Dr.Hesham ! \n')

%% variables declaration
s = input('Enter area in ft^2: ');  %ft^2 
b = input('Enter span in ft: ');    %ft

AR=b^2/s;
fprintf('\n Aspect Ratio= %d \n \n',AR)

C_Do = input('Enter C_Do: ');
e = input('Enter Oswald efficiency factor (e): ');
C_L_Max = input('Enter Maximum lift coefficient from (C_L,alpha) graph: ');

Wo = input('Enter gross weight in Ib: ');       %Ib
n = input('Enter Propeller efficiency: ');
f = input('Enter fuel capacity in gal: ');      %gallon
Wf = f*5.64; %fuel weight in Ib
Wi= Wo - Wf; 
c = input('Enter fuel consumption Ib/(hp)(h): ');
Piston = input('Enter power of one piston (hp): ');

V = input('Enter Aircraft Velocity in ft/s: '); 
V = round(V);
h = input('Enter gliding height in ft = ');
rho_0 = 0.002378;           % density @ sea-level
rho_20 = 12.67e-4 ;         % density at any height to calc slope t_min   

%% Powers & Thrusts

P_a=Piston*n;    % power available Ib.ft/s
                        
v_points=0:400;    % velocity vector intialization

C_L=Wo/0.5/rho_0./v_points.^2/s; % Vec

k = 1/e/AR/pi;
C_D=C_Do+(k*C_L.^2);             % Vec

T_R=Wo*C_D./C_L;      % Thrust required Vector 

P_R=T_R.*v_points/550;    % Power required Vector

T_a=P_a*550./v_points;    % Thrust required Vector

%% Velocities

V_stall = sqrt(2*Wo/(rho_0*s*C_L_Max));
% 
 [~ , idx] = min(abs(P_R - P_a)); %gets nearest value in P_R matrix equal to P_a
 V_max = v_points(idx);


function [f] = solve_function(v_point)
f = (0.5*rho_0*v_point(1).^2*s*C_Do + (2*k*s./(rho_0*v_point(1)))*(Wo/s)^2)-P_a*550;
end

beta = 20;  %intial value for intersection
warning('off')
intersection = fsolve(@solve_function, beta,optimoptions('fsolve','Display','off'));


if (intersection <= V_stall)
         
    V_min = V_stall;            

elseif (intersection > V_stall)
    
    V_min = intersection;
end

[~,idx2] = min(abs(v_points-intersection));

figure
plot(v_points,T_R) 
hold on
plot(v_points,T_a)
hold on

plot(V_max,T_R(V_max),'ko',intersection,T_a(idx2),'ko')
text(V_max,T_R(V_max)+1200,'\downarrow V_m_a_x')
text(intersection,T_a(idx2)+1400,'\downarrow V_m_i_n')

xlim([10,400])
ylim([-0.5e4,2.5e4])
xlabel('V ft/s')
ylabel('Thrust Ib')
grid on 
legend('T_R','T_A')
title('Thrust & Velocity')

v_L_D_max_0=v_points(P_R == min(P_R));

figure
plot(v_points,P_R);
hold on
plot(v_points,P_a*ones(length(v_points)), intersection,P_a,'ko', V_min,P_a,'ko', V_max,P_a,'ko', v_L_D_max_0,min(P_R),'ko');

text(intersection+2,P_a+25,'V_m_i_n')
text( V_min+2,P_a+25,'V_s_t_a_l_l')
text(V_max+2,P_a+25,'V_m_a_x')
text(v_L_D_max_0+2,min(P_R)-25,'V_P_m_i_n')


xlim([10,400])
xlabel('V ft/s')
ylabel('Power (hp)')
grid on 
legend('P_R','P_A','location','southeast')
title('Power & Velocity')

%% Maximum theta
V_thetamax = 4*(Wo/s)*k/(rho_0*n*(P_a*550/Wo));
V_thetamax = round(V_thetamax);
theta_climb_max=asind(T_a(V_thetamax)/Wo-0.5*rho_0*V_thetamax^2*(s/Wo)*C_Do-(Wo/s)*2*k/(rho_0*V_thetamax^2));

%% Rate of climb

excesspower_0=(P_a-P_R)*550;         
R_C_0=excesspower_0/Wo;         % Rate of climb at sea-level
R_C_Max_0=max(R_C_0);           % Maximum '' '' ''

v_R_C_Max = sqrt(2*Wo/(rho_0*s)*sqrt(k/3/C_Do)); %velocity at maximum R/C

%R_C_0 = tan(theta_climb_max)*v_points;
%x = find (R_C_0 == R_C_Max_0); 
%plot(v_points,R_C_0);

figure
plot(v_points,R_C_0,v_R_C_Max,R_C_Max_0,'ko');
text(v_R_C_Max,R_C_Max_0-5,'\uparrow (R/C)_m_a_x')

xlim([10 400])
xlabel('V ft/s')
ylabel('(R/C) ft/s')
grid on 
title('Rate of climb')

%% Maximum ranges and endurance
R_gliding = h/sqrt(4*k*C_Do);

c = c*1/550/3660; %fuel consumption unit convertion  
R_max = n/c*C_L(V)/C_D(V)*log(Wo/Wi); %Maximum Range

C_D_Max = (C_Do + k*C_L_Max^2); 
E = n/c*C_L_Max^(3/2)/C_D_Max*sqrt(2*rho_0*s)*(Wi^(-0.5)-Wo^(-0.5)); %maximum Endurance in seconds

%% Absolute and service ceiling

P_a_alt=P_a*(rho_20/rho_0);         %Power available at any altitude
P_R_alt=P_R*(rho_0/rho_20)^0.5;     %Power required at any altitude
excesspower_20=(P_a_alt-P_R_alt)*550;
R_C_20=excesspower_20/Wo;           % R/C at 20e3 ft
R_C_20_Max = max(R_C_20);

b_slope = 20e3/(R_C_20_Max - R_C_Max_0); 
h_abs = -R_C_Max_0*b_slope;

height = h_abs + b_slope*R_C_0;        

[~, idx1] = min( abs(R_C_0 - (5/3)) ); % nearest value to 5/3 which is the service rate of climb 
R_C_abs = R_C_0(idx1);
h_service = height(idx1);

figure
plot(R_C_0,height, 0,h_abs,'ko', R_C_Max_0,0,'ko', R_C_abs,h_service,'ko');
xlim([0 30])
text(0.4,h_abs,'\leftarrow h_a_b_s_o_l_u_t_e')
text(R_C_abs+0.4,h_service,'\leftarrow h_s_e_r_v_i_c_e')
text(R_C_Max_0,1e3,'\downarrow (R/C)_m_a_x')
xlabel('(R/C) ft/s')
ylabel('H ft')
grid on 
title('Absolute & service ceiling')

%% minimum time

slope2 = 1/b_slope;
h_vec =0:1e3:27e3;      %height intialization 
R_C = R_C_Max_0 + slope2*h_vec; 

t_min =1/slope2*(log(R_C_Max_0+slope2*20e3)-log(R_C_Max_0));

figure
plot(h_vec,R_C, h_abs,0,'ko', 0,R_C_Max_0,'ko');
ylim([0 30])
text(h_abs,1.5,'h_a_b_s_o_l_u_t_e')
text(1e3,R_C_Max_0,'\leftarrow (R/C)_m_a_x')
xlabel('H ft')
ylabel('(R/C) ft/s ')
grid on
title('Minimum time derivation')

%% Table 
  figure
col = {'performance variable','calculated value','Unit'};
row = {'Aspect ratio','Lift coeffecient','Drag coeffecient', 'Thrust available',...;
        'Thrust required','Power available','Power required','Maximum speed','Minimum speed','Stall speed',...;
        'Rate of climb at sea-level','Maximum Rate of climb','Velocity at maximum (R/C)',...;
        'Gliding Range','Maximum Range','Absolute ceiling','Service ceiling','Maximum climb angle'...;
        'Minimum time of climb','Maximum endurance'};
dat={'AR',AR,'~'; 'C_L',C_L(V),'~'; 'C_D',C_D(V),'~';
    'T_a',T_a(V),'Ib'; 'T_r',T_R(V),'Ib';
    'P_a',P_a,'hp'; 'P_r',P_R(V),'hp';
    'V_max',V_max,'ft/s'; 'V_min',V_min,'ft/s'; 'V_stall',V_stall,'ft/s'; 
    '(R/C)o',R_C_0(V),'ft/s'; '(R/C)max',R_C_Max_0,'ft/s'; 
    'V_(R/C)max',v_R_C_Max,'ft/s'; 'R_gliding',R_gliding,'ft'; 'R_max',R_max,'ft';
    'h_absolute',h_abs,'ft'; 'h_service',h_service,'ft'; 'Theta_max',theta_climb_max,'degree';
    't_min',t_min/60,'min';'E_max',E/3600,'hr'};
 uitable ('columnname',col,'rowname',row,'data',dat,'Position',[0 15 568 382],'Fontsize',8);
 
 end