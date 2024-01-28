function [charge,P0] = ChargeDensity(E1,E2,E3,E4,P1,P2,P3,P4)
%%% Input:
%%% Electric fields measured by each of the four MMS/EDP satellites  
%%% E1,E2,E3,E4  
%%% the unit is mV/m 
%%% Satellite location measured by each of the four MMS satellites  
%%% P1,P2,P3,P4  
%%% the unit is km 
%%%%%%%%%%%%%%%%%%%
%%% Output
%%% Charge density and position at the centre of the satellite constellation
%%% Charge density e m^-3
%%% Position km

P0=(P1+P2+P3+P4)/4;%% position at the centre of the satellite constellation

r1=P1(1:3)-P0(1:3);
r2=P2(1:3)-P0(1:3);
r3=P3(1:3)-P0(1:3);
r4=P4(1:3)-P0(1:3);
clear P1 P2 P3 P4

%%%The electric field gradient tensor is calculated with reference to 
%%%"Shen, C., Li, X., Dunlop, M., Liu, Z. X., Balogh, A., Baker, D. N., . . . Wang, X. (2003). 
%%% Analyses on the geometrical structure of magnetic field in the current sheet based on cluster measurements. 
%%% doi:https://doi.org/10.1029/2002JA009612".

R(:,:)=(r1(1:3)'*r1(1:3)+r2(1:3)'*r2(1:3)+r3(1:3)'*r3(1:3)+r4(1:3)'*r4(1:3))/4;%volume tensor 1/4sigma(rak*raj)
R0(:,:)=R(:,:)\eye(3);%********The inverse of a volume tensor

ER(:,:)=(E1(1:3)'*r1(1:3)+E2(1:3)'*r2(1:3)+E3(1:3)'*r3(1:3)+E4(1:3)'*r4(1:3))/4;%tensor 1/4sigma(Bai*rak)
nabla(:,:)=ER(:,:)*R0(:,:);
charge=squeeze(nabla(1,1)+nabla(2,2)+nabla(3,3))*8.854/0.1602;%%% Output units are e m^-3

end

