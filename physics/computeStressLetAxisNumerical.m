%compute force on panel

function S = computeStressLetAxisNumerical(x,y,fxfy,PARAM)

xDecay = logspace(3,5,100);
yDecay = zeros(1,numel(xDecay));

%compute velocity field
[~,~,uField,vField] = computeVelPressField(xDecay,yDecay,x,y,fxfy,zeros(2*numel(xDecay),1),PARAM,0,0);
Uabs = sqrt(uField.^2+vField.^2);

%fit
ft = fittype({'x^(-2)'});
f1 = fit(xDecay',uField',ft);
S = f1.a;

% figure
% loglog(xDecay,Uabs,'o')
% hold on
% loglog(xDecay,abs(S)*xDecay.^(-2))
% legend('Numerical','Fitting')
% grid on
% xlabel('d')
% xlabel('|U|')
% title(['S=' num2str(S)])

%S = -S*8*pi/3/A;