%remesh drop spectrally

function [y2,halt] = remeshDropSpectral(y,PARAM)

xMode = y(1:2:end-1);
yMode = y(2:2:end);
y2 = zeros(numel(y),1);

if PARAM.remesh==1
        
    %remesh mantainig spectral accuracy
    if PARAM.legendre==2

          x = @(t) interpLegendreZeroOne(t,xMode);
          y = @(t) interpLegendreZeroOne(t,yMode);

          %go on GLL grid for remesshing
          if PARAM.legendre==1
              PARAM.D1 = PARAM.D1gll;
              PARAM.D2 = PARAM.D2gll;
              PARAM.WG = PARAM.WGgll;
              PARAM.PPP = PARAM.PPPgll;
              PARAM.legendre = 2;
              PARAM.t = PARAM.tgll;
          end
          PPP = PARAM.PPP;

          %compute arclenght
          [x,y,halt] = remeshSpectralFun(PARAM.t,x,y,PARAM);

          %compute coeffs
          xMode = LegendreSerieSpectralXY(x,PPP,PARAM);
          yMode = LegendreSerieSpectralXY(y,PPP,PARAM);

    elseif PARAM.legendre==0

          x = chebcoeffs2chebvals(xMode);
          y = chebcoeffs2chebvals(yMode);

          x = chebfun(x,[0 1]);
          y = chebfun(y,[0 1]);

          %compute arclenght
          [x,y,halt] = remeshSpectralFun(PARAM.t,x,y,PARAM);
          x = chebfun(x,[0 1]);     y = chebfun(y,[0 1]);

          xMode = chebcoeffs(x);
          yMode = chebcoeffs(y);

    elseif PARAM.legendre==1

          x = interpLegendreZeroOne(PARAM.tcheb,xMode);
          y = interpLegendreZeroOne(PARAM.tcheb,yMode);

          x = chebfun(x,[0 1]);
          y = chebfun(y,[0 1]);

          %go on GLL grid for remesshing
          if PARAM.legendre==1
              PARAM.D1 = PARAM.D1cheb;
              PARAM.D2 = PARAM.D2cheb;
              PARAM.WG = PARAM.WGcheb;
              PARAM.legendre = 0;
              PARAM.t = PARAM.tcheb;
          end

          %compute arclenght
          [x,y,halt] = remeshSpectralFun(PARAM.t,x,y,PARAM);
          x = chebfun(x,[0 1]);     y = chebfun(y,[0 1]);

          xMode = chebcoeffs(x);
          yMode = chebcoeffs(y);

          xMode = cheb2leg(xMode);
          yMode = cheb2leg(yMode);

    end

elseif PARAM.remesh==2
    
    %remesh using SPlines
    if PARAM.legendre==0
        
          %this never fails but is less accurate
          halt = 0;

          %t equispaced
          tHere = linspace(0,1,numel(xMode));
          x = chebcoeffs2chebvals(xMode);
          y = chebcoeffs2chebvals(yMode);
          x = chebfun(x,[0 1]);
          y = chebfun(y,[0 1]);
          x = x(tHere);
          y = y(tHere);
          
          %compute splines coefficients
          l = int_spline_symmetric_lenght(x,y);
          
          %fit spline to arc lenght
          x = spline(l,x,PARAM.t*l(end));
          y = spline(l,y,PARAM.t*l(end));
          
          %first and last point are on the axis
          y([1 end]) = [0 0];
          
          %chebfun
          x = chebfun(x,[0 1]);     y = chebfun(y,[0 1]);

          %modes coeff
          xMode = chebcoeffs(x);
          yMode = chebcoeffs(y);
          
    elseif PARAM.legendre==1
        
          warning('the selected remesh is not efficient with non Lobatto points')
        
          %this never fails but is less accurate
          halt = 0;

          %t equispaced
          tHere = linspace(0,1,numel(xMode));
          x = interpLegendreZeroOne(tHere,xMode);
          y = interpLegendreZeroOne(tHere,yMode);
          
          %compute splines coefficients
          l = int_spline_symmetric_lenght(x',y');
          
          %fit spline to arc lenght
          x = spline(l,x,PARAM.t*l(end));
          y = spline(l,y,PARAM.t*l(end));

          %modes coefficients
          xMode = LegendreSerieSpectralXY(x,PARAM.PPP,PARAM);
          yMode = LegendreSerieSpectralXY(y,PARAM.PPP,PARAM);

    else
        
       error('no remeshing implemented') 
        
    end
    
end

%de-aliasing
xMode(PARAM.dealiasing+1:end) = 0;
yMode(PARAM.dealiasing+1:end) = 0;

y2(1:2:end-1) = xMode;
y2(2:2:end) = yMode;
        
