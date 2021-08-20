
%**************************************************************************
%
%   add_references_geographic_political: 
%    
%**************************************************************************
function handles=add_references_geographic_political(handlesinput,color)

path=handlesinput.geographicreferencespath;

set(gcf,'renderer','zbuffer');
hold on

% add references
gatrib(1).array=load([path '/stateboundariesnan.txt']);
gatrib(2).array=load([path '/embaymentlonlatborder.txt']);
if(color==0)
    gatrib(1).color=[.1 0 .9];
    gatrib(2).color=[0 0 0];
    linewidth=1;
else
    gatrib(1).color=[1 1 1];
    gatrib(2).color=[1 1 1];
    linewidth=2;
end

for i=1:length(gatrib)
    plot3(gatrib(i).array(:,1),gatrib(i).array(:,2),gatrib(i).array(:,2)*0+1000,...
        'color',gatrib(i).color,'LineWidth',linewidth)
end

if(color==0)
    colorfonts=[.25 .25 .25];
    text(-92.5,35.,10000,'Mississippi Embayment', 'Rotation',55,'FontSize',10,...
        'Color',[0 0 0])
    text(-88, 35.7,1000,'Tennessee','Color',colorfonts, 'FontSize',10);
    text(-87, 37.3,1000,'Kentucky','Color',colorfonts, 'FontSize',10);
    text(-86.8, 39,1000,'Indiana','Color',colorfonts, 'FontSize',10);
    text(-89.3, 39.3,1000,'Illinois','Color',colorfonts, 'FontSize',10);
    text(-92, 38,1000,'Missouri','Color',colorfonts, 'FontSize',10);
    text(-92.5, 34.,1000,'Arkansas','Color',colorfonts, 'FontSize',10);
    text(-90., 34.,1000,'Mississippi','Color',colorfonts, 'FontSize',10);
    text(-87.25, 34.,1000,'Alabama','Color',colorfonts, 'FontSize',10);
end
hold off
