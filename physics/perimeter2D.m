%compute perimeter of a 2D figure using splines approximation

function P = perimeter2D(x,y)

  %interface normal and curvature
  [ax, bx, cx, dx, ay, by, cy, dy] = my_spline_periodic (x, y);
  %compute splines coordinates
  t = 0:0.1:0.9;
  ttt = repmat(t,1,numel(ax));
  axxx = reshape(repmat(ax,numel(t),1),1,numel(ax)*numel(t));
  bxxx = reshape(repmat(bx,numel(t),1),1,numel(bx)*numel(t));
  cxxx = reshape(repmat(cx,numel(t),1),1,numel(cx)*numel(t));
  dxxx = reshape(repmat(dx,numel(t),1),1,numel(dx)*numel(t));
  ayyy = reshape(repmat(ay,numel(t),1),1,numel(ay)*numel(t));
  byyy = reshape(repmat(by,numel(t),1),1,numel(by)*numel(t));
  cyyy = reshape(repmat(cy,numel(t),1),1,numel(cy)*numel(t));
  dyyy = reshape(repmat(dy,numel(t),1),1,numel(dy)*numel(t));
  
  %splines coordinates
  xxx = [axxx+bxxx.*ttt+cxxx.*ttt.^2+dxxx.*ttt.^3 x(end)];
  yyy = [ayyy+byyy.*ttt+cyyy.*ttt.^2+dyyy.*ttt.^3 y(end)];
  
  %subdivision length
  dl = sqrt(diff(xxx).^2+diff(yyy).^2);
  P = sum(dl);

end