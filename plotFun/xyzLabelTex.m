%plot with Latex symbols

function xyzLabelTex(xlab,ylab,zlab)

xlabel(['$' xlab '$'],'interpreter','latex')
ylabel(['$' ylab '$'],'interpreter','latex')
zlabel(['$' zlab '$'],'interpreter','latex')