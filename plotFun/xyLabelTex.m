%plot with Latex symbols

function xyLabelTex(xlab,ylab)

xlabel(['$' xlab '$'],'interpreter','latex')
ylabel(['$' ylab '$'],'interpreter','latex')