%**************************************************************************
%
%   load_handles_default: input parameters to initiallize the database
%    
%**************************************************************************
function handles = load_handles_default(handles)
handles.update = 0;
handles.databasepath = ' ';
handles.resamplefactortopo = 10;

handles.geographicreferencespath = ' ';

handles.doiplotelevationordepth = 1;
handles.doiplotvp = 0;
handles.doiplotvs = 0;
handles.doiplotrho = 0;
handles.doiplotunits = 0;

%% Cross-sections
handles.crosssectionquery.points.ya = 0;
handles.crosssectionquery.points.xa = 0;
handles.crosssectionquery.points.zMax = 0;
handles.crosssectionquery.points.yb = 0;
handles.crosssectionquery.points.xb = 0;
handles.crosssectionquery.points.zMin = 0;
handles.crosssectionquery.points.samplesalong = 100;
handles.crosssectionquery.points.samplesZorDepth = 100;
handles.crosssectionquery.points.depthMin = 0;
handles.crosssectionquery.points.depthMax = 0;


%% Point operations
handles.pointquery.doiplotvp = 0;
handles.pointquery.doiplotvs = 0;
handles.pointquery.doiplotrho = 0;
handles.pointquery.doiplotunit = 0;
handles.pointquery.x = 0;
handles.pointquery.y = 0;
handles.pointquery.zMin = 0;
handles.pointquery.zMax = 1000;

handles.pointquery.update = 0;
handles.pointquery.numSamples = 100;

% Map views
handles.mapviewquery.pivot.y = 0;
handles.mapviewquery.pivot.x = 0;
handles.mapviewquery.lengthx = 1;
handles.mapviewquery.lengthy = 1;
handles.mapviewquery.samplesx = 100;
handles.mapviewquery.samplesy = 100;
handles.mapviewquery.pivot.z = 0;

handles = handles;

return





