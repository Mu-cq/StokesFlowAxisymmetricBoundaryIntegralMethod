%draw shape like Campas

function R = residualsShapeCampas(alpha,rhoQ,beta,V0,PARAM)

%dependent vaibales
rho = rhoQ(1:end-1);
Q = rhoQ(end);

%compute first derivative of rho
rhop = PARAM.D1*rho;

%compute residuals
R = beta*cos(alpha).^2+Q-1./rhop.*cos(alpha)-sin(alpha)./rho;
R = [R; VolumePolarAxisSpectral(alpha,rho,PARAM)-V0];