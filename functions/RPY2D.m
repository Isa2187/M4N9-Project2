function [M] = RPY2D(Xi,Xj,a)

% This function evaluates the RPY tensor assuming the particles i and j are
% in the xy-plane.
%
% Inputs: Xi: Position of particle i.
%         Xj: Position of particle j.
%          a: radius of the particles.
%
% Output:  M: 2x2 RPY tensor.
%

DX = Xi - Xj;
r = norm(DX);
Dxh = DX./(r+10^-16);

fac1 = 1/(8*pi*r);
fac2 = 1/(6*pi*a);

if(r > 2*a)
    M(1,1) = fac1*(1 + 2*a^2/(3*r^2)) + fac1*(1 - 2*a^2/(r^2))*Dxh(1)*Dxh(1);
    M(1,2) = fac1*(1 - 2*a^2/(r^2))*Dxh(1)*Dxh(2);
    M(2,2) = fac1*(1 + 2*a^2/(3*r^2)) + fac1*(1 - 2*a^2/(r^2))*Dxh(2)*Dxh(2);
    M(2,1) = M(1,2);
else
    M(1,1) = fac2*(1 - 9*r/(32*a)) + fac2*3*r/(32*a)*Dxh(1)*Dxh(1);
    M(1,2) = fac2*3*r/(32*a)*Dxh(1)*Dxh(2);
    M(2,2) = fac2*(1 - 9*r/(32*a)) + fac2*3*r/(32*a)*Dxh(2)*Dxh(2);
    M(2,1) = M(1,2);
end

end

