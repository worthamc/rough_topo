function BVPtest

% testing bvp4c

% constants and vectors

% initial guess

gamma = 15;
solinit = bvpinit(linspace(0,pi,10),@mat4init,gamma);

% solver
sol = bvpsolver(@VertStrct,@VertStructBC,solinit);

fprintf('The ___ eigenvalue is approximately %7.3f.\n',sol.parameters)

zint = linspace(0,zmin);
Szint = deval(sol,zint);
figure;
plot(Szint(1,:),zint);
title('Eigenfunction'); 
xlabel('amplitude');
ylabel('depth');

% initial guess

% equation definition
    function dPhidz = VertStruct(z,Phi,gamma)
        dPhidz = [ Phi(2)
         -(lambda - 2*q*cos(2*x))*y(1) ];
    end

% boundary condition definition

end












