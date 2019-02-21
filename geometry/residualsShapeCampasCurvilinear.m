%draw shape like Campas

function R = residualsShapeCampasCurvilinear(t,alphaRhoQ,beta,V0,PARAM)

%dependent vaibales
alpha = alphaRhoQ(1:numel(t));
rho = alphaRhoQ(numel(t)+1:end-1);
Q = alphaRhoQ(end);

%compute first derivative of rho
rhop = PARAM.D1*rho;
alphap = PARAM.D1*alpha;

%compute residuals
R = beta*cos(alpha).^2+Q-alphap./rhop.*cos(alpha)-sin(alpha)./rho;
R = [R; VolumePolarAxisSpectral(alpha,rho,PARAM)-V0];