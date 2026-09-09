
function [V,p,Cp,Profile] = physical_quantities(M,rho,u)

%Profil kridla - body
Bmark = BoundaryPts(M);
Profile.pts = find(Bmark == 1);
Profile.coord = M.POS(Profile.pts,:);

Profile.upper = find( Profile.coord(:,2) >= 0 );
Profile.lower = find( Profile.coord(:,2) <= 0 );
Profile.length = unique( Profile.coord(:,1) );

%Speed
V.upper = sqrt( u.uux(Profile.upper).^2 + u.uuy(Profile.upper).^2 );
V.lower = sqrt( u.uux(Profile.lower).^2 + u.uuy(Profile.lower).^2 );
V.wing = sqrt( u.uux(Profile.pts).^2 + u.uuy(Profile.pts).^2 );
V.Omega = sqrt( u.uux.^2 + u.uuy.^2 );

%Pressure
p.upper = 1/2 * rho * ( u.U.^2 -  abs(V.upper).^2);
p.lower = 1/2 * rho * ( u.U.^2 -  abs(V.lower).^2);
p.wing =  1/2 * rho * ( u.U.^2 -  abs(V.wing).^2);
p.Omega = 1/2 * rho * ( u.U.^2 - abs(V.Omega).^2);
p.rozdil = p.lower - p.upper;

%Pressure coefficient
Cp.upper = (p.upper) / ( 0.5 * rho * u.U.^2);
Cp.lower = (p.lower) / ( 0.5 * rho * u.U.^2);
Cp.ratio = Cp.upper ./ Cp.lower;
Cp.diff = Cp.upper - Cp.lower;

end
