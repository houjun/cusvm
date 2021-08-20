% GEODATAQGUI M-file for geodataqgui.fig
%      GEODATAQGUI, by itself, creates a new GEODATAQGUI or raises the existing
%      singleton*.
%
%      H = GEODATAQGUI returns the handle to a new GEODATAQGUI or the handle to
%      the existing singleton*.
%
%      GEODATAQGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GEODATAQGUI.M with the given input arguments.
%
%      GEODATAQGUI('Property','Value',...) creates a new GEODATAQGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before geodataqgui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop. All inputs are passed to geodataqgui_OpeningFcn via varargin.
%
%      *See GUI Options on MATLAB's GUIDE's Tools menu. x Choose "GUI allows only i
%      one instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
function varargout = geodataqgui(varargin)

% Edit the above text to modify the response to help geodataqgui

% Last Modified by GUIDE v2.5 27-Oct-2011 17:48:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @geodataqgui_OpeningFcn, ...
                   'gui_OutputFcn',  @geodataqgui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before geodataqgui is made visible.
function geodataqgui_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles=load_handles_default(handles);
guidata(hObject, handles);

function varargout = geodataqgui_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% *************************************************************************
%                          INPUT EDITTEXT BOXES 
% *************************************************************************
function y=pointquery_textinordisplay_xorlat_Callback(hObject,eventdata,handles)
handles.pointquery.x=str2double(get(hObject,'String'));
if(handles.pointquery.update==1)
   handles.pointquery.x=handles.pointquery.ginput(2);
   set(handles.pointquery_textinordisplay_xorlat,'String',num2str(handles.pointquery.x));
end

set(handles.pointquery_textinordisplay_xorlat,'String',num2str(handles.pointquery.x));
handles.pointquery.points.update=0;
guidata(hObject, handles);
y=handles;

function pointquery_textinordisplay_xorlat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function y=pointquery_textinordisplay_yorlon_Callback(hObject, eventdata, handles)
handles.pointquery.y=str2double(get(hObject,'String'));
if(handles.pointquery.update==1)
   handles.pointquery.y=handles.pointquery.ginput(1);
   set(handles.pointquery_textinordisplay_yorlon,'String',num2str(handles.pointquery.y));
end

set(handles.pointquery_textinordisplay_yorlon,'String',num2str(handles.pointquery.y));
handles.pointquery.update=0;
guidata(hObject, handles);
y=handles;

function pointquery_textinordisplay_yorlon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pointquery_textinordisplay_Zmin_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.pointquery.zMin=str2double(get(hObject,'String'));
guidata(hObject, handles);

function pointquery_textinordisplay_Zmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pointquery_textinordisplay_Zmax_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.pointquery.zMax=str2double(get(hObject,'String'));
guidata(hObject, handles);

function pointquery_textinordisplay_Zmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pointquery_textinordisplay_numSamples_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.pointquery.numSamples=str2double(get(hObject,'String'));
guidata(hObject, handles);

function pointquery_textinordisplay_numSamples_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit12_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.mapviewquery.pivot.z=str2double(get(hObject,'String'));
guidata(hObject, handles);

function edit12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit13_Callback(hObject, eventdata, handles)

function edit13_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%**************************************************************************
%                           CROSS-SECTION  
%**************************************************************************

%************************** POINT A X OR LAT ******************************
function y=textcrosssectionpoinaxorlat_Callback(hObject, eventdata, handles)
handles.crosssectionquery.points.xa=str2double(get(hObject,'String'));
if(handles.crosssectionquery.points.update==1)
   handles.crosssectionquery.points.xa=handles.crosssectionquery.points.ginput(1,2);
   set(handles.textcrosssectionpoinaxorlat,'String',num2str(handles.crosssectionquery.points.xa));
end

set(handles.textcrosssectionpoinaxorlat,'String',handles.crosssectionquery.points.xa);
handles.crosssectionquery.points.update=0;
guidata(hObject, handles);
y=handles;

function textcrosssectionpoinaxorlat_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** POINT A Y OR LON ******************************
function y=textcrosssectionpoinayorlon_Callback(hObject, eventdata, handles)
handles.crosssectionquery.points.ya=str2double(get(hObject,'String'));
if(handles.crosssectionquery.points.update==1)
   handles.crosssectionquery.points.ya=handles.crosssectionquery.points.ginput(1,1);
   set(handles.textcrosssectionpoinayorlon,'String',num2str(handles.crosssectionquery.points.ya));
end
set(handles.textcrosssectionpoinayorlon,'String',handles.crosssectionquery.points.ya);
handles.crosssectionquery.points.update=0;
guidata(hObject, handles);
y=handles;

function textcrosssectionpoinayorlon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** POINT B X OR LAT ******************************
function y=textcrosssectionpoinbxorlat_Callback(hObject, eventdata, handles)
handles.crosssectionquery.points.xb=str2double(get(hObject,'String'));
if(handles.crosssectionquery.points.update==1)
   handles.crosssectionquery.points.xb=handles.crosssectionquery.points.ginput(2,2);
   set(handles.textcrosssectionpoinbxorlat,'String',num2str(handles.crosssectionquery.points.xb));
end
set(handles.textcrosssectionpoinbxorlat,'String',handles.crosssectionquery.points.xb);
handles.crosssectionquery.points.update=0;
guidata(hObject, handles);
y=handles;

function textcrosssectionpoinbxorlat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** POINT B Y OR LON ******************************
function y=textcrosssectionpoinbyorlon_Callback(hObject, eventdata, handles)
handles.crosssectionquery.points.yb=str2double(get(hObject,'String'));
if(handles.crosssectionquery.points.update==1)
   handles.crosssectionquery.points.yb=handles.crosssectionquery.points.ginput(2,1);
   set(handles.textcrosssectionpoinbyorlon,'String',num2str(handles.crosssectionquery.points.yb));
end
set(handles.textcrosssectionpoinbyorlon,'String',handles.crosssectionquery.points.yb);
handles.crosssectionquery.points.update=0;
guidata(hObject, handles);
y=handles;

function textcrosssectionpoinbyorlon_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** SAMPLES HORIZONTAL*****************************
function y=samplescrosssectionhorizontal_Callback(hObject, eventdata, handles)
handles.crosssectionquery.points.samplesalong=str2double(get(hObject,'String'));
guidata(hObject, handles);
y=handles;

function samplescrosssectionhorizontal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** TEXT RADIOBUTTONELEVATION MIN*****************************
function textelevationmin_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.crosssectionquery.points.zMin=str2double(get(hObject,'String'));
guidata(hObject, handles);

function textelevationmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** TEXT RADIOBUTTONELEVATION MAX*****************************
function textelevationmax_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.crosssectionquery.points.zMax=str2double(get(hObject,'String'));
guidata(hObject, handles);

function textelevationmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** TEXT RADIOBUTTONDEPTH MIN ********************************
function textdepthmin_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.crosssectionquery.points.depthMin=str2double(get(hObject,'String'));
guidata(hObject, handles);

function textdepthmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** TEXT RADIOBUTTONDEPTH MAX ********************************
function textdepthmax_Callback(hObject, eventdata, handles)
get(hObject,'String');
handles.crosssectionquery.points.depthMax=str2double(get(hObject,'String'));
guidata(hObject, handles);

function textdepthmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************** SAMPLES Z *************************************
function samplesz_Callback(hObject, eventdata, handles)
handles.crosssectionquery.points.samplesZorDepth=str2double(get(hObject,'String'));
guidata(hObject, handles);
y=handles;

function samplesz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% *************************************************************************
% *************************************************************************
%
%                  INITIALIZATION AND RETRIEVING PARAMETERS
%
% *************************************************************************
% *************************************************************************

%**************************************************************************
%                   LOAD DATABASE PATH
%**************************************************************************
function loaddatabasepath_Callback(hObject, eventdata, handles)
if(handles.update==1)   
   set(handles.loaddatabasepath,'String',handles.databasepath);
end
set(handles.loaddatabasepath,'String',handles.databasepath);
handles.update=0;

guidata(hObject, handles);

function loaddatabasepath_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%**************************************************************************
%                  UPDATE RESPAMPLE TOPO FACTOR 
%**************************************************************************
function loadsampletopofactor_Callback(hObject, eventdata, handles)
handles.resamplefactortopo=str2double(get(hObject,'String'));
guidata(hObject, handles);

function loadsampletopofactor_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%**************************************************************************
%                  UPDATE REFERENCES PATHS
%**************************************************************************
function loadgeographicreferencespath_Callback(hObject, eventdata, handles)
if(handles.update==1)   
   set(handles.loadgeographicreferencespath,'String',handles.geographicreferencespath);
end
guidata(hObject, handles);

function loadgeographicreferencespath_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%**************************************************************************
%               INITIALIZE DATABASE AND DISPLAY TOPOGRAPHY                 
%**************************************************************************
function initializedisplaydatabase_Callback(hObject, eventdata, handles)
reset=0;
handles.databasepath=uigetdir();
handles.update=1;
loaddatabasepath_Callback(hObject, eventdata, handles);
load_database_displaytopography(handles,reset);
guidata(hObject, handles);

%**************************************************************************
%                            ZOOM  
%**************************************************************************
function zoom_Callback(hObject, eventdata, handles)
zoom on

%**************************************************************************
%                          DISPLAY REFERENCES                 
%**************************************************************************
function addgeographicpoliticalreferences_Callback(hObject, eventdata, handles)
handles.geographicreferencespath=uigetdir();
handles.update=1;
loadgeographicreferencespath_Callback(hObject, eventdata, handles);
color=0; % gui topography 
add_references_geographic_political(handles,color);

%**************************************************************************
%                          RESET DISPLAY
%**************************************************************************
function resetdisplay_Callback(hObject, eventdata, handles)
reset=1;
load_database_displaytopography(handles,reset);

%**************************************************************************
%                          PICK CROSS-SECTION
%**************************************************************************
function pickcrossection_Callback(hObject, eventdata, handles)
handles.crosssectionquery.points.ginput=ginput(2);
handles.crosssectionquery.points.update=1;
guidata(hObject, handles);
handles=textcrosssectionpoinaxorlat_Callback(hObject, eventdata, handles);
handles.crosssectionquery.points.update=1;
handles=textcrosssectionpoinayorlon_Callback(hObject, eventdata, handles);
handles.crosssectionquery.points.update=1;
handles=textcrosssectionpoinbxorlat_Callback(hObject, eventdata, handles);
handles.crosssectionquery.points.update=1;
handles=textcrosssectionpoinbyorlon_Callback(hObject, eventdata, handles);
handles.crosssectionquery.points.update=1;

hold on;
scatter3(handles.crosssectionquery.points.ya,handles.crosssectionquery.points.xa,10000,'ko','Filled');
scatter3(handles.crosssectionquery.points.yb,handles.crosssectionquery.points.xb,10000,'ko','Filled');

plot3(handles.crosssectionquery.points.ginput(:,1),handles.crosssectionquery.points.ginput(:,2),...
    handles.crosssectionquery.points.ginput(:,1)*0+10000,'k');

hold off;

function pushbutton4_Callback(hObject, eventdata, handles)
function pushbutton5_Callback(hObject, eventdata, handles)
function pushbutton6_Callback(hObject, eventdata, handles)
function pushbutton7_Callback(hObject, eventdata, handles)
function pushbutton8_Callback(hObject, eventdata, handles)
function pushbutton9_Callback(hObject, eventdata, handles)
function pushbutton10_Callback(hObject, eventdata, handles)

function pushbutton2_Callback(hObject, eventdata, handles)
function pushbutton3_Callback(hObject, eventdata, handles)

function pushbutton15_Callback(hObject, eventdata, handles)
function pushbutton14_Callback(hObject, eventdata, handles)
function pushbutton13_Callback(hObject, eventdata, handles)
function pushbutton12_Callback(hObject, eventdata, handles)
function pushbutton11_Callback(hObject, eventdata, handles)
function pushbutton19_Callback(hObject, eventdata, handles)
function pushbutton20_Callback(hObject, eventdata, handles)

function pushbutton23_Callback(hObject, eventdata, handles)

% ***************PLOT PROPERTIES CROSS-SECTIONS ***************************

function plotpropertycrosssection_Callback(hObject, eventdata, handles)

% Write the input file
fileName='input.in';
fp=fopen('input.in','w');
fprintf(fp,'\n hororvert   = 1');
fprintf(fp,'\n originLon   = %f', handles.crosssectionquery.points.ya);
fprintf(fp,'\n originLat   = %f', handles.crosssectionquery.points.xa);
fprintf(fp,'\n originDepth = %f', handles.crosssectionquery.points.zMax);
fprintf(fp,'\n finalLon    = %f', handles.crosssectionquery.points.yb);
fprintf(fp,'\n finalLat    = %f', handles.crosssectionquery.points.xb);
fprintf(fp,'\n finalDepth  = %f', handles.crosssectionquery.points.zMin);
fprintf(fp,'\n numAlong    = %d', handles.crosssectionquery.points.samplesalong);
fprintf(fp,'\n numD        = %d', handles.crosssectionquery.points.samplesZorDepth);
fclose(fp);
program= './geodataq';
arg0=' 1 ';
if(handles.doiplotelevationordepth==0)  % elevation
    arg1.string=' 0 ';  % velorveldepthorbdrortopo
    arg1.val=0;
else
    arg1.string=' 1 ';                              % depth
    arg1.val=1;
end
arg2= handles.databasepath;
arg3='./input.in';
arg4='./';

space=' ';
if(system( [program arg0 space arg1.string space arg2 space arg3 space arg4 ])~=0)
    error('Error');
end
   
% **************************DISPLAY ***************************************
toPlot=[handles.doiplotvp handles.doiplotvs handles.doiplotrho handles.doiplotunits];
propertyName(1).string='vp';
propertyName(2).string='vs';
propertyName(3).string='rho';
propertyName(4).string='containingunit';

totalPlots=sum(toPlot);

if(totalPlots > 0)         
    colormapD=jet(100);
    iPlot=0;
    for i=1:4
        
        if(toPlot(i)==1)
            figureR;
            fpProperty= fopen([propertyName(i).string '.out'],'r');
            hov       = fscanf(fpProperty,'%d',1);
            origin    = fscanf(fpProperty,'%e',3);
            final     = fscanf(fpProperty,'%e',3);
            nGrid     = fscanf(fpProperty,'%d',2);
            if(i==4)
                numberOfUnits=fscanf(fpProperty,'%d',1);
                for iUnits=1:numberOfUnits;
                    units(iUnits).name=fscanf(fpProperty,'%s',1);
                end
                units(numberOfUnits+1).name=' ';
            end
            property  = fscanf(fpProperty,'%e',[nGrid(2),nGrid(1)]);
            
            z=linspace(origin(3),final(3),nGrid(2));
            y=linspace(origin(1),final(1),nGrid(1));
            x=linspace(origin(2),final(2),nGrid(1));

            if(i~=4)
                k=find(property < 0);
                property(k)=NaN;
            end
            if abs(final(1)-origin(1)) > abs(final(2)-origin(2))
            	contourf(y,z,property,40,'LineStyle','none');
		xlab = 'Lon';
            else
            	contourf(x,z,property,40,'LineStyle','none');
		xlab = 'Lat';
            end
            view(2);
            if(i==4) % quantize
                maximum=max(max(property));
                minimum=min(min(property));
                divisionsCM=maximum-minimum+1;
                colormapQ=jet(numberOfUnits+1);
                colormapQ(1,:)=[1 1 1];
                colormap(colormapQ(minimum+1:maximum+1,:));
                caxis([minimum-.5 maximum+.5]);
                h=colorbar('YTickLabel',{units(minimum+1:maximum+1).name});
                                                
            else
                maximum=max(max(property));
                minimum=min(min(property)); 
                if(maximum ~= minimum)                    
                    caxis([minimum maximum]);                
                end
                colormap(colormapD);
                colorbar                
            end
                        
            box on
            title(propertyName(i).string)
            ylabel('Elevation (m asl)')
            xlabel(xlab)
            min(min(property));
            if(handles.doiplotelevationordepth==1)
                axis ij;
                ylabel('Depth (m)')
            end
            
        end
    end
end


function plotpropertycrosssectiondepth_Callback(hObject, eventdata, handles)


% Write the input file
fileName='input.in';
fp=fopen('input.in','w');
fprintf(fp,'\n hororvert   = 1');
fprintf(fp,'\n originLon   = %f', handles.crosssectionquery.points.ya);
fprintf(fp,'\n originLat   = %f', handles.crosssectionquery.points.xa);
fprintf(fp,'\n originDepth = %f', handles.crosssectionquery.points.depthMin);
fprintf(fp,'\n finalLon    = %f', handles.crosssectionquery.points.yb);
fprintf(fp,'\n finalLat    = %f', handles.crosssectionquery.points.xb);
fprintf(fp,'\n finalDepth  = %f', handles.crosssectionquery.points.depthMax);
fprintf(fp,'\n numAlong    = %d', handles.crosssectionquery.points.samplesalong);
fprintf(fp,'\n numD        = %d', handles.crosssectionquery.points.samplesZorDepth);
fclose(fp);
program= './geodataq';
arg0=' 1 ';
arg1=' 1 '; % different options to output data (1) is depth
arg2= handles.databasepath;
arg3='./input.in';
arg4='./';

space=' ';
if(system( [program arg0 space arg1 space arg2 space arg3 space arg4 ])~=0)
    error('Error');
end
   
% **************************DISPLAY ***************************************
toPlot=[handles.doiplotvp handles.doiplotvs handles.doiplotrho handles.doiplotunits];
propertyName(1).string='vp.out';
propertyName(2).string='vs.out';
propertyName(3).string='rho.out';
propertyName(4).string='containingunit.out';
for i=1:4
    if(toPlot(i)==1)
        figureR;
        fpProperty= fopen(propertyName(i).string,'r');
        hov       = fscanf(fpProperty,'%d',1);
        origin    = fscanf(fpProperty,'%e',3);
        final     = fscanf(fpProperty,'%e',3);
        nGrid     = fscanf(fpProperty,'%d',2);
        property  = fscanf(fpProperty,'%e',[nGrid(2),nGrid(1)]);
        
        z=linspace(origin(3),final(3),nGrid(2));
        y=linspace(origin(1),final(1),nGrid(1));
        x=linspace(origin(2),final(2),nGrid(1));
        
        surf(y,z,property)
        view(2)
        shading interp
        colorbar
        box on
        title('Vs')
        ylabel('Elevation (m asl)')
        xlabel('Lon')
        axis ij
        min(min(property));
    end
end

function doIPlotVp_Callback(hObject, eventdata, handles)
handles.doiplotvp=get(hObject,'Value');
guidata(hObject, handles);

function doIPlotVs_Callback(hObject, eventdata, handles)
handles.doiplotvs=get(hObject,'Value');
guidata(hObject, handles);
 
function doIPlotDensity_Callback(hObject, eventdata, handles)
handles.doiplotrho=get(hObject,'Value');
guidata(hObject, handles);

function doIPlotUnit_Callback(hObject, eventdata, handles)
handles.doiplotunits=get(hObject,'Value');
guidata(hObject, handles);

% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5

% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6

% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7

% --- Executes on button press in checkbox8.
function checkbox8_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox8

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9

% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10

% ******   POINT RELATED FUNCTIONS ****************************************

function pointquery_chooseproperty_doiplotelevation_Callback(hObject, eventdata, handles)
handles.pointquery.doiplotelevation=get(hObject,'Value');
guidata(hObject, handles);

function pointquery_chooseproperty_doiplotdepth_Callback(hObject, eventdata, handles)
handles.pointquery.doiplotdepth=get(hObject,'Value');
guidata(hObject, handles);

function general_doiplotvp_Callback(hObject, eventdata, handles)
handles.doiplotvp=get(hObject,'Value');
guidata(hObject, handles);

function point_doiplotvs_Callback(hObject, eventdata, handles)
handles.pointquery.doiplotvs=get(hObject,'Value');
guidata(hObject, handles);

function point_doiplotrho_Callback(hObject, eventdata, handles)
handles.pointquery.doiplotrho=get(hObject,'Value');
guidata(hObject, handles);

function point_doiplotunit_Callback(hObject, eventdata, handles)
handles.pointquery.doiplotunit=get(hObject,'Value');
guidata(hObject, handles);

function pointquery_displayproperty_Callback(hObject, eventdata, handles)
% Write the input file
fileName='input.in';
fp=fopen('input.in','w');
fprintf(fp,'\n hororvert   = 1');
fprintf(fp,'\n originLon   = %f', handles.pointquery.y);
fprintf(fp,'\n originLat   = %f', handles.pointquery.x);
fprintf(fp,'\n originDepth = %f', handles.pointquery.zMin);
fprintf(fp,'\n finalLon    = %f', handles.pointquery.y);
fprintf(fp,'\n finalLat    = %f', handles.pointquery.x);
fprintf(fp,'\n finalDepth  = %f', handles.pointquery.zMax);
fprintf(fp,'\n numAlong    = %d', 1);
fprintf(fp,'\n numD        = %d', handles.pointquery.numSamples);
fclose(fp);
program= './geodataq';
arg0=' 1 ';  % type of query
if(handles.doiplotelevationordepth==0);  % elevation
    arg1.string=' 0 ';  % velorveldepthorbdrortopo
    arg1.val=0;
else
    arg1.string=' 1 ';                              % depth
    arg1.val=1;
end

arg2= handles.databasepath;
arg3='./input.in';
arg4='./';

space=' ';

if(system( [program arg0 space arg1.string space arg2 space arg3 space arg4 ])~=0)
    error('Error')
end

% **************************DISPLAY ***************************************
toPlot=[handles.doiplotvp handles.doiplotvs ...
    handles.doiplotrho handles.doiplotunits];
propertyName(1).string='vp';
propertyName(2).string='vs';
propertyName(3).string='rho';
propertyName(4).string='containingunit';

totalPlots=sum(toPlot);

if(totalPlots > 0)
    figureR;
    iPlot=0;
    for i=1:4
        if(toPlot(i)==1)
            iPlot=iPlot+1;
            subplot(1,totalPlots,iPlot)
            fpProperty= fopen([propertyName(i).string '.out'],'r');
            hov       = fscanf(fpProperty,'%d',1);
            origin    = fscanf(fpProperty,'%e',3);
            final     = fscanf(fpProperty,'%e',3);
            nGrid     = fscanf(fpProperty,'%d',2);
            if(i==4)
                numberOfUnits=fscanf(fpProperty,'%d',1);
                for iUnits=1:numberOfUnits;
                    units(iUnits).name=fscanf(fpProperty,'%s',1);
                end
            end
            property  = fscanf(fpProperty,'%e',[nGrid(2),nGrid(1)]);
            
            z=linspace(origin(3),final(3),nGrid(2));
            y=linspace(origin(1),final(1),nGrid(1));
            x=linspace(origin(2),final(2),nGrid(1));
            
            if(i~=4)
                k=find(property < 0);
                property(k)=NaN;
            end
            plot(property,z)
            
            box on
            title(propertyName(i).string)
            
            if(arg1.val==0)
                ylabel('Elevation (m asl)')
            else
                ylabel('Depth (m)')
                axis ij
            end
        end
    end
end

function pointquery_choosepoints_Callback(hObject, eventdata, handles)
handles.pointquery.ginput=ginput(1);
handles.pointquery.update=1;
guidata(hObject, handles);
handles=pointquery_textinordisplay_xorlat_Callback(hObject, eventdata, handles);
handles.pointquery.update=1;
handles=pointquery_textinordisplay_yorlon_Callback(hObject, eventdata, handles);

hold on;
plot3(handles.pointquery.y,handles.pointquery.x,10000,'ko','LineStyle','none','MarkerFaceColor',[0 0 0]);

hold off;

% --- Executes during object creation, after setting all properties.
function uipanel10_CreateFcn(hObject, eventdata, handles)

function uipanel10_SelectionChangeFcn(hObject, eventdata, handles)

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobuttonelevation'
        handles.doiplotelevationordepth=0;
        guidata(hObject, handles);        
    case 'radiobuttondepth'
        handles.doiplotelevationordepth=1;
        guidata(hObject, handles);
    otherwise
        
end

function pointquery_chooseproperty_doiplotdepth_DeleteFcn(hObject, eventdata, handles)

% *************FUNCTIONS TO QUERY MAP VIEW ********************************

function pick_pivot_mapview_Callback(hObject, eventdata, handles)
handles.mapviewquery.ginput=ginput(1);
handles.mapviewquery.update=1;
guidata(hObject, handles);
handles=mapviewquery_textinordisplay_xorlatmin_Callback(hObject, eventdata, handles);
handles.mapviewquery.update=1;
handles=mapviewquery_textinordisplay_yorlonmin_Callback(hObject, eventdata, handles);

hold on;
plot3(handles.mapviewquery.pivot.y,handles.mapviewquery.pivot.x,10000,'ko','LineStyle','none','MarkerFaceColor',[0 0 0]);

hold off;

function y=mapviewquery_textinordisplay_xorlatmin_Callback(hObject, eventdata, handles)
handles.mapviewquery.pivot.x=str2double(get(hObject,'String'));
if(handles.mapviewquery.update==1)
   handles.mapviewquery.pivot.x=handles.mapviewquery.ginput(2);
   set(handles.mapviewquery_textinordisplay_xorlatmin,'String',num2str(handles.mapviewquery.pivot.x));
end

set(handles.mapviewquery_textinordisplay_xorlatmin,'String',num2str(handles.mapviewquery.pivot.x));
handles.mapviewquery.points.update=0;
guidata(hObject, handles);
y=handles;

function mapviewquery_textinordisplay_xorlatmin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%function mapviewquery_textinordisplay_xorlatmin_Callback(hObject, eventdata, handles)

function y=mapviewquery_textinordisplay_yorlonmin_Callback(hObject, eventdata, handles)
handles.mapviewquery.pivot.y=str2double(get(hObject,'String'));
if(handles.mapviewquery.update==1)
   handles.mapviewquery.pivot.y=handles.mapviewquery.ginput(1);
   set(handles.mapviewquery_textinordisplay_yorlonmin,'String',num2str(handles.mapviewquery.pivot.y));
end

set(handles.mapviewquery_textinordisplay_yorlonmin,'String',num2str(handles.mapviewquery.pivot.y));
% draw a square

x0=handles.mapviewquery.pivot.x;
y0=handles.mapviewquery.pivot.y;
lx=handles.mapviewquery.lengthx;
ly=handles.mapviewquery.lengthy;

x=[x0 x0+lx x0+lx x0 x0];
y=[y0 y0 y0+ly y0+ly y0];
hold on
plot3(y,x,x*0+1000,'k')
hold off
handles.mapviewquery.points.update=0;
guidata(hObject, handles);
y=handles;

function mapviewquery_textinordisplay_yorlonmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mapviewquery_textinordisplay_lengthx_Callback(hObject, eventdata, handles)
get(hObject,'String')
handles.mapviewquery.lengthx=str2double(get(hObject,'String'));
guidata(hObject, handles);

function mapviewquery_textinordisplay_lengthx_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mapviewquery_textinordisplay_lengthy_Callback(hObject, eventdata, handles)
get(hObject,'String')
handles.mapviewquery.lengthy=str2double(get(hObject,'String'));
guidata(hObject, handles);

function mapviewquery_textinordisplay_lengthy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mapviewquery_textinordisplay_samplesx_Callback(hObject, eventdata, handles)
get(hObject,'String')
handles.mapviewquery.samplesx=str2double(get(hObject,'String'));
guidata(hObject, handles);

function mapviewquery_textinordisplay_samplesx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mapviewquery_textinordisplay_samplesy_Callback(hObject, eventdata, handles)
get(hObject,'String')
handles.mapviewquery.samplesy=str2double(get(hObject,'String'));
guidata(hObject, handles);

function mapviewquery_textinordisplay_samplesy_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton16_Callback(hObject, eventdata, handles)

% Write the input file
fileName='input.in';
fp=fopen('input.in','w');
fprintf(fp,'\n hororvert   = 0');
fprintf(fp,'\n originLon   = %f', handles.mapviewquery.pivot.y);
fprintf(fp,'\n originLat   = %f', handles.mapviewquery.pivot.x);
fprintf(fp,'\n originDepth = %f', handles.mapviewquery.pivot.z);
fprintf(fp,'\n lengthLat   = %f', handles.mapviewquery.lengthx);
fprintf(fp,'\n lengthLon   = %f', handles.mapviewquery.lengthy);
fprintf(fp,'\n numLon      = %f', handles.mapviewquery.samplesy);
fprintf(fp,'\n numLat      = %f', handles.mapviewquery.samplesx);
fclose(fp);
program= './geodataq';
arg0=' 1 ';
if(handles.doiplotelevationordepth==0);  % elevation
    arg1.string=' 0 ';
    arg1.val=0;
else
    arg1.string=' 1 ';                              % depth
    arg1.val=1;
end

arg2= handles.databasepath;
arg3='./input.in';
arg4='./';

space=' ';
if(system( [program arg0 space arg1.string space arg2 space arg3 space arg4 ])~=0)
    error('Error')
end

% **************************DISPLAY ***************************************
toPlot=[handles.doiplotvp handles.doiplotvs handles.doiplotrho handles.doiplotunits];
propertyName(1).string='vp';
propertyName(2).string='vs';
propertyName(3).string='rho';
propertyName(4).string='containingunit';

totalPlots=sum(toPlot);

if(totalPlots > 0)        
    colormapD=jet(100);    
    iPlot=0;
    for i=1:4
        
        if(toPlot(i)==1)
            figureR;
            fpProperty   = fopen([propertyName(i).string '.out'],'r');
            hov          = fscanf(fpProperty,'%d',1);
            origin       = fscanf(fpProperty,'%e',3);
            nGrid        = fscanf(fpProperty,'%d',3);
            lengthLonLat = fscanf(fpProperty,'%e',2);
            depth        = fscanf(fpProperty,'%e',1);
            if(i==4)
                numberOfUnits=fscanf(fpProperty,'%d',1);
                for iUnits=1:numberOfUnits;
                    units(iUnits).name=fscanf(fpProperty,'%s',1);
                end
            end
            
            property  = fscanf(fpProperty,'%e',[nGrid(2),nGrid(1)]);
            
            final=origin+[lengthLonLat; 0];
            
            z=linspace(origin(3),final(3),nGrid(2));
            y=linspace(origin(1),final(1),nGrid(1));
            x=linspace(origin(2),final(2),nGrid(1));
            
            if(i~=4)
                k=find(property < 0);
                property(k)=NaN;
            end            
            contourf(y,x,property,10,'LineColor','none');
            view(2)
            if(i==4) % quantize
                maximum=max(max(property));
                minimum=min(min(property));
                divisionsCM=maximum-minimum+1;
                colormapQ=jet(numberOfUnits+1);
                colormapQ(1,:)=[1 1 1];
                colormap(colormapQ(minimum+1:maximum+1,:));
                caxis([minimum-.5 maximum+.5]);
                colorbar('YTickLabel',{units(minimum+1:maximum+1).name})                
            else
                maximum=max(max(property));
                minimum=min(min(property));
                caxis([minimum maximum]);
                colormap(colormapD);
                colorbar
            end
                        
            box on
            title([propertyName(i).string])
            ylabel('Lat')
            xlabel('Lon')
            min(min(property));
            color=1; % property
            add_references_geographic_political(handles,color);
        end
    end
end

% --- Executes when selected object is changed in uipanel12.
function uipanel12_SelectionChangeFcn(hObject, eventdata, handles)
if hObject== handles.radiobuttonelevation
        handles.doiplotelevationordepth=0;
        guidata(hObject, handles);
elseif hObject == handles.radiobuttondepth
        handles.doiplotelevationordepth=1;
        guidata(hObject, handles);
        
guidata(hObject, handles);       

end
