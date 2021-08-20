%%
%  Transform a set of point given in lon lat to XY
%%
function [xi,yi]=transformLonLatToXY(boxSquare,lengthX,lengthY,loni,lati)

 csii =  [ -1 -1  1  1 ];
 ethai = [ -1  1  1 -1 ];
 
 boxlon=boxSquare(:,1)';
 boxlat=boxSquare(:,2)';
 
 for i=1:length(loni)
     aux(i,1:2)= ethaandcsi(csii, ethai, loni(i),lati(i), boxlon, boxlat)';
     aux(i,1)  = (aux(i,1)+1)*lengthX*.5;
     aux(i,2)  = (aux(i,2)+1)*lengthY*.5;
 end
 
 xi=aux(:,1);
 yi=aux(:,2);
 
 return