clear
close all

% Creates test VAWT geometry file

% Add geom creation scripts to path
path(path,'../../CreateGeom');

% Params
R=0.5;            % Center radius (ft)
HR=2;           % Height to radius ratio 
CRr=0.28;        % Root chord to radius
eta=.5;             % Blade mount point ratio (mount point behind leading edge as a fraction of chord)
NBlade=3;
NBElem=6;
NStrut=3;       % number of struts
NSElem=5;
CRs=CRr;        % strut chord to radius
TCs=0.20;        % strut thickness to chord

% Output filename
FN='unh-rvat.geom';

% Plot data?
PlotTurbine=0;

% Convert
dToR=pi/180;

% Create basic parabolic blade VAWT
Type='VAWT';
BShape=0;
T=CreateTurbine(NBlade,NBElem,NStrut,NSElem,R,[],[],[],Type,1,CRr,HR,eta,BShape,CRs,TCs);

% Write geom file
WriteTurbineGeom(FN,T);

% Plot if desired
if PlotTurbine
    
    % Plot animated turbine rotation
    XLim=[-4,4];
    YLim=[-2,4];
    ZLim=[-4,4];
    
    % Plot controls
    PlotVec=1;
    SFVec=.5;
    Trans=.5;
    
    hf=figure(1);
    set(hf,'Position',[303   124   956   610]) 
    set(gca,'Position',[5.2743e-002  5.1245e-002  8.9979e-001  8.8141e-001])
    set(gca,'CameraPosition',[-52.1999   30.4749   62.2119])
    set(gca,'CameraUpVector',[1.8643e-001  9.7433e-001 -1.2615e-001])
    set(gca,'CameraViewAngle',6.3060e+000)
    grid on
    set(gcf,'Color','white');
    hl=light('Position',[-1,0,0]);
    set(gca,'Color','white');
    set(gca,'DataAspectRatio',[1,1,1])
    set(gca,'XLim',XLim,'YLim',YLim,'ZLim',ZLim)
    
    HIn=[];
    PhasePlot=linspace(0,2*pi,150);
    for i=1:length(PhasePlot)
       H=PlotTurbineGeom(T,hf,PhasePlot(i),HIn,Trans,PlotVec,SFVec);
       HIn=H;
       pause(.01);
    end
    
end


