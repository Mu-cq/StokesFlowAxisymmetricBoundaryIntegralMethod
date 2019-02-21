function figuresize( w , h , u )
p = inputParser;
p.addRequired('width', @(x) isnumeric(x) && all(size(x)==1) );
p.addRequired('height',@(x) isnumeric(x) && all(size(x)==1) );
p.addOptional('units','centimeters',...
  @(x) any(strcmpi(x,{'normalized','centimeters','inches','points'})) );

p.parse( w, h, u );
w = p.Results.width;
h = p.Results.height;
u = p.Results.units;

p = 0.01;

set(gcf,'Units',u);
screenpos = get(gcf,'Position');

set(gcf,...
  'Position',[screenpos(1:2) w h],...
  'PaperUnits',u,...
  'PaperPosition',[p*w p*h w h],...
  'PaperSize',[w*(1+2*p) h*(1+2*p)]);

end