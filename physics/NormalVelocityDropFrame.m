%normal velocity in drop frame

function [UnFrame,Vdrop] = NormalVelocityDropFrame(x,y,u,v)

    Un = DropNormalVelocity(x,y,u,v);
    Vdrop = DropVelocityAxis(x,y,Un);
    UnFrame = DropNormalVelocity(x,y,u-Vdrop,v);

end