
%**************************************************************************
%
%   load_database_displaytopography:
%
%   reset = 0 no
%           1 yes
%
%**************************************************************************
function handles=load_database_displaytopography(handles,reset)

databasePath=handles.databasepath;

%Load the topography
topo=loadsurface([databasePath '/surf_air.in']);
topo.limYMin=min(topo.loni);
topo.limXMin=min(topo.lati);
topo.limYMax=max(topo.loni);
topo.limXMax=max(topo.lati);
handles.topo=topo;
resampleF=handles.resamplefactortopo;

% plot a resampled version of the topography
a=length(handles.topo.lati);
b=length(handles.topo.loni);

set(gcf,'renderer','zbuffer');
if(reset==0)
    hold on;
else
    hold off;
end

contourf(handles.topo.loni(1:resampleF:b), handles.topo.lati(1:resampleF:a),...
    handles.topo.raster(1:resampleF:a,1:resampleF:b),20,'LineStyle','none');

if(reset==1)
    hold on;
end
load colormapregion.txt
xlabel('Y or Lon')
ylabel('X or Lat')
axis([topo.limYMin topo.limYMax topo.limXMin topo.limXMax 0 10000])
view(2);
colormap(colormapregion)
colorbar
grid on
box on
hold off
