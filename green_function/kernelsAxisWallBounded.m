function [MXX,MXY,MYX,MYY,QXXX,QXXY,QXYX,QXYY,QYXX,QYXY,QYYX,QYYY] ...
...
  = kernelsAxisWallBounded(X,Y,X0,Y0,wallPos)

        %check that the geometry id on the right side of the wall
        if min(min(X0))<wallPos
            error('Geometry has to be on the right side of the wall')
        end
        
        h = X0-wallPos;
        X0im = 2*wallPos-X0;
        
        [MstXX,MstXY,MstYX,MstYY,QstXXX,QstXXY,QstXYX,QstXYY,QstYXX,QstYXY,QstYYX,QstYYY] =...
            sgf_ax_fs_vect3 (X,Y,X0,Y0);
        
        [MstXXim,MstXYim,MstYXim,MstYYim,QstXXXim,QstXXYim,QstXYXim,QstXYYim,QstYXXim,QstYXYim,QstYYXim,QstYYYim] =...
            sgf_ax_fs_vect3 (X,Y,X0im,Y0);
        
        [MdXX,MdXY,MdYX,MdYY,QdXXX,QdXXY,QdXYX,QdXYY,QdYXX,QdYXY,QdYYX,QdYYY,...
            MsdXX,MsdXY,MsdYX,MsdYY,QsdXXX,QsdXXY,QsdXYX,QsdXYY,QsdYXX,QsdYXY,QsdYYX,QsdYYY] =...
            kernelsAxisDipolesAndDoublet (X,Y,X0im,Y0);
        
        %Assembly single layer
        MXX = MstXX - MstXXim + 2*h.^2.*MdXX - 2*h.*MsdXX;
        MXY = MstXY - MstXYim + 2*h.^2.*MdXY - 2*h.*MsdXY;
        MYX = MstYX - MstYXim + 2*h.^2.*MdYX - 2*h.*MsdYX;
        MYY = MstYY - MstYYim + 2*h.^2.*MdYY - 2*h.*MsdYY;
        
        %Assembly doubel layer
        QXXX = QstXXX - QstXXXim + 2*h.^2.*QdXXX - 2*h.*QsdXXX;
        QXXY = QstXXY - QstXXYim + 2*h.^2.*QdXXY - 2*h.*QsdXXY;
        QXYX = QstXYX - QstXYXim + 2*h.^2.*QdXYX - 2*h.*QsdXYX;
        QXYY = QstXYY - QstXYYim + 2*h.^2.*QdXYY - 2*h.*QsdXYY;
        QYXX = QstYXX - QstYXXim + 2*h.^2.*QdYXX - 2*h.*QsdYXX;
        QYXY = QstYXY - QstYXYim + 2*h.^2.*QdYXY - 2*h.*QsdYXY;
        QYYX = QstYYX - QstYYXim + 2*h.^2.*QdYYX - 2*h.*QsdYYX;
        QYYY = QstYYY - QstYYYim + 2*h.^2.*QdYYY - 2*h.*QsdYYY;
        
return
