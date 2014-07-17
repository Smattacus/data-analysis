function varargout = GUI_001(varargin)
% GUI_001 MATLAB code for GUI_001.fig
%      GUI_001, by itself, creates a new GUI_001 or raises the existing
%      singleton*.
%
%      H = GUI_001 returns the handle to a new GUI_001 or the handle to
%      the existing singleton*.
%
%      GUI_001('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_001.M with the given input arguments.
%
%      GUI_001('Property','Value',...) creates a new GUI_001 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_001_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_001_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_001

% Last Modified by GUIDE v2.5 16-Jul-2013 12:30:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_001_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_001_OutputFcn, ...
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


% --- Executes just before GUI_001 is made visible.
function GUI_001_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_001 (see VARARGIN)
% Choose default command line output for GUI_001
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
dm = createDM('Dye_611.66005-means.mat', true, 0.01);
handles.testdata = dm;
guidata(hObject, handles);


% UIWAIT makes GUI_001 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_001_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%a1 = subplot(2, 1, 1, 'Parent', handles.uipanel1);
%semilogy(handles.testdata(1,:), handles.testdata(2,:));
%a2 = subplot(2,1,2, 'Parent', handles.uipanel1);
%semilogy(handles.testdata(1,:), handles.testdata(2,:) * 0);
ipf(handles.testdata);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes3);
semilogy(handles.testdata(1,:), handles.testdata(2,:));


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
% Interprets key presses from the Figure window.
global X Y xx yy xo dx start FitResults NumPeaks NumTrials Shape residual FIXEDPARAMETERS
global delta AA xxx PEAKHEIGHTS AUTOZERO extra MeanFitError SavedSignal logplot et
% When a key is pressed, interprets key and calls corresponding function.
% Note: If you don't like my key assignments, you can change the numbers
% in the case statements here to re-assign that function to any other key.
NumTrialsBoot=1;
key= eventdata.Character;
if isscalar(key),
    ly=length(Y);
    switch double(key),
        case 28
            % Pans down when right arrow pressed.
            xo=xo-dx/20;
            if xo<1,xo=1;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 29
            % Pans up when left arrow pressed.
            xo=xo+dx/20;
            if xo>ly,xo=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 91
            % Nudge down 1 point when [ pressed.
            xo=xo-1;
            if xo<1,xo=1;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 93
            % Nudge up 1 point when ] pressed.
            xo=xo+1;
            if xo>ly,xo=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 44
            % Pans down when > key pressed.
            xo=xo-dx;
            if xo<1,xo=1;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 46
            % Pans up when < key pressed.
            xo=xo+dx;
            if xo>ly,xo=ly;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 31
            % Zooms up when up arrow pressed.
            dx=dx+dx/10;
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 30
            % Zooms down when down arrow pressed.
            dx=dx-dx/10;
            if dx<2,dx=2;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 47
            % Zooms 2% up when / pressed.
            dx=dx+ly/50;
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 39
            % Zooms 2% down when ' pressed.
            dx=dx-ly/50;
            if dx<2,dx=2;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 109
            % When 'M' is pressed, toggles on/off log plot mode
            if logplot==0,
                logplot=1;
            else
                logplot=0;
            end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 27 % When 'ESC' key is pressed, resets pan and zoom
            xo=length(Y)/2; % Initial Pan setting
            dx=length(Y)/4; % Initial Zoom setting
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 102
            % When 'f' key is pressed, tweaks start values, computes and plots fit.
            startnow=start;
            delta=(max(xx)-min(xx))/100;
            for k=1:2:2*NumPeaks,
                startnow(k)=start(k)+randn*delta;
            end
            [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,startnow,extra);
        case 99
            % When 'c' key is pressed, user clicks graph to enter start positons,
            % then fit is computed and graph re-drawn.
            % Acquire first-guess peak positions from user mouse pointer
            figure(1);subplot(2,1,1);xlabel('Click on the estimated positions of each proposed component peak.')
            [clickX,clickY] = ginput(NumPeaks);
            % Create a "start" vector using these peak positions, with peak
            % widths
            n=max(xx)-min(xx);
            width=n/(5*NumPeaks);
            start=[];
            for k=1:NumPeaks,
                start=[start clickX(k) width];
            end
            if Shape==11||12,
                fixedstart=[];
                for pk=1:NumPeaks,
                    fixedstart(pk)=start(2*pk-1);
                end
                peakfit([xx;yy],0,0,NumPeaks,Shape,extra,1,start,AUTOZERO,FIXEDPARAMETERS);
            end
            % [FitResults,MeanFitError]=peakfit([xx;yy],0,0,NumPeaks,Shape,extra,1,start);
            [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
        case 98
            % When 'b' key is pressed, user clicks graph before and after peaks
            % to enter background points, then fit is computed and graph re-drawn.
            figure(1);subplot(2,1,1);xlabel('Click on the baseline to the LEFT the peak(s).')
            [X1,Y1] = ginput(1);
            figure(1);subplot(2,1,1);xlabel('Now click on the baseline to the RIGHT the peak(s).')
            [X2,Y2] = ginput(1);
            n=length(xx);
            %  Create "start" for this number of peaks
            yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
            [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
        case 8
            % When 'Backspace' key is pressed, user clicks the graph
            % along the presumed background points, then the program
            % subtracts the interploated background between the points.
            SavedSignal=Y;
            BaselinePoints=input('Number of baseline points to click): ');
            if isempty(BaselinePoints),BaselinePoints=8;end
            % Acquire background points from user mouse clicks
            subplot(2,1,2)
            title(['Click on ' num2str(BaselinePoints) ' points on the baseline between the peaks.'])
            bX=[];bY=[];
            for g=1:BaselinePoints;
                [clickX,clickY] = ginput(1);
                bX(g)=clickX;
                bY(g)=clickY;
                xlabel(['Baseline point '  num2str(g) ' / ' num2str(BaselinePoints) ])
            end
            yy=Y;
            for k=1:length(bX)-1,
                fp=val2ind(X,bX(k)); % First point in segment
                lp=val2ind(X,bX(k+1));  % Last point in segment
                % Subtract piecewise linear background from Y
                yy(fp:lp)=Y(fp:lp)-((bY(k+1)-bY(k))/(bX(k+1)-bX(k))*(X(fp:lp)-bX(k))+bY(k));
            end
            Y=yy;
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 92
            % When '\' key is pressed, restoreed original signal
            Y=SavedSignal;
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case {49,50,51,52,53,54,55,56,57}
            % When a number key is pressed, sets the number of peaks to that number.
            n=key-48;
            ShapeString='';
            if round(n)~=NumPeaks,
                NumPeaks=round(n);
                switch Shape
                    case 1
                        ShapeString='Gaussian';
                    case 2
                        ShapeString='Lorentzian';
                    case 3
                        ShapeString='logistic';
                    case 4
                        ShapeString='Pearson7';
                    case 5
                        ShapeString='ExpGaussian';
                    case 6
                        ShapeString='Equal width Gaussians';
                    case 7
                        ShapeString='Equal width Lorentzians';
                    case 8
                        ShapeString='Equal-width ExpGauss.';
                    case 9
                        ShapeString='Exponental pulse';
                    case 10
                        ShapeString='Sigmoid';
                    case 11
                        ShapeString='Fixed-width Gaussian';
                    case 12
                        ShapeString='Fixed-width Lorentzian';
                    case 13
                        ShapeString='Gauss/Lorentz blend';
                    case 14
                        ShapeString='bifurcated Gaussian';
                    case 15
                        ShapeString='bifurcated Lorentzian';
                    case 16
                        ShapeString='Fixed-position Gaussians';
                    case 17
                        ShapeString='Fixed-position Lorentzians';
                    otherwise
                        ShapeString='';
                end % switch
                subplot(2,1,2)
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
            end % if
            n=max(xx)-min(xx);
            % Create a start value for this number of peaks
            start=[];
            startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker);
                start=[start markx n/(5*NumPeaks)];
            end
            PEAKHEIGHTS=zeros(NumPeaks,1);
            AA=zeros(NumPeaks,200);
            RedrawSignal(X,Y,xo,dx,NumPeaks);
            xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
        case 48  % Zero key is pressed to enter number of peaks
            NumPeaks=input('Number of peaks: ');
            if isempty(NumPeaks),NumPeaks=1;end
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='logistic';
                case 4
                    ShapeString='Pearson7';
                case 5
                    ShapeString='ExpGaussian';
                case 6
                    ShapeString='Equal width Gaussians';
                case 7
                    ShapeString='Equal width Lorentzians';
                case 8
                    ShapeString='Equal-width ExpGauss.';
                case 9
                    ShapeString='Exponental pulse';
                case 10
                    ShapeString='Sigmoid';
                case 11
                    ShapeString='Fixed-width Gaussian';
                case 12
                    ShapeString='Fixed-width Lorentzian';
                case 13
                    ShapeString='Gauss/Lorentz blend';
                case 14
                    ShapeString='bifurcated Gaussian';
                case 15
                    ShapeString='bifurcated Lorentzian';
                case 16
                    ShapeString='Fixed-position Gaussians';
                case 17
                    ShapeString='Fixed-position Lorentzians';
                otherwise
                    ShapeString='';
            end % switch
            % Create a start value for this number of peaks
            n=max(xx)-min(xx);
            start=[];
            startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+min(xx);
            for marker=1:NumPeaks,
                markx=startpos(marker);
                start=[start markx n/(5*NumPeaks)];
            end
            PEAKHEIGHTS=zeros(NumPeaks,1);
            AA=zeros(NumPeaks,200);
            RedrawSignal(X,Y,xo,dx,NumPeaks);
            xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
        case {103,108,111,112,101,104,59,106,117,115,71,76,96,72,58,80,123}
            % Selects peak shape when G,L,O,P,E,H,:,J,U,S, `, Shift-H or : key is pressed
            switch key
                case 103 % When 'g' key is pressed, peak shape is set to Gaussian.
                    n=1;
                case 108 % When 'l' key is pressed, peak shape is set to Lorentzian.
                    n=2;
                case 111 % When 'o' key is pressed, peak shape is set to Logistic.
                    n=3;
                case 112 % When 'p' key is pressed, peak shape is set to Pearson.
                    n=4;
                case 101 % When 'e' key is pressed, peak shape is set to exponentally-broadened gaussian.
                    n=5;
                case 104 % When 'h' key is pressed, peak shape is set to equal width gaussians.
                    n=6;
                case 59 % When ';' key is pressed, peak shape is set to equal width gaussians.
                    n=7;
                case 106 % When 'J' key is pressed, peak shape is set to exponentally-broadened equal width gaussians.
                    n=8;
                case 117 % When 'U' key is pressed, peak shape is set to exponental pulse.
                    n=9;
                case 115 % When 'S' key is pressed, peak shape is set to Sigmoid.
                    n=10;
                case 71 % When 'Shift-G' key is pressed, peak shape is set to FWGaussian.
                    n=11;
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 76 % When ';' key is pressed, peak shape is set to FWLorentzian.
                    n=12;
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 96 % When '`' key is pressed, peak shape is set to Gauss/Lorentz blend
                    n=13;
                case 72 % When 'Shift-H' key is pressed, peak shape is set to bifurcated Gaussian
                    n=14;
                case 58 % When ':' key is pressed, peak shape is set to bifurcated Lorentzian
                    n=15;
                case 80 % When 'Shift-P' key is pressed, peak shape is set to
                    n=16;
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                    end
                case 123 % When 'Shift-[' key is pressed, peak shape is set to
                    n=17;
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                    end
                otherwise
            end
            % switch
            if round(n)~=Shape,
                Shape=round(n);
                switch Shape
                    case 1
                        ShapeString='Gaussian';
                    case 2
                        ShapeString='Lorentzian';
                    case 3
                        ShapeString='logistic';
                    case 4
                        ShapeString='Pearson';
                    case 5
                        ShapeString='ExpGaussian';
                    case 6
                        ShapeString='Equal-width Gaussian';
                    case 7
                        ShapeString='Equal-width Lorentzian';
                    case 8
                        ShapeString='Equal-width ExpGauss.';
                    case 9
                        ShapeString='Exponental pulse';
                    case 10
                        ShapeString='Sigmoid';
                    case 11
                        ShapeString='Fixed-width Gaussian';
                    case 12
                        ShapeString='Fixed-width Lorentzian';
                    case 13
                        ShapeString='Gauss/Lorentz blend';
                    case 14
                        ShapeString='bifurcated Gaussian';
                    case 15
                        ShapeString='bifurcated Lorentzian';
                    case 16
                        ShapeString='Fixed-position Gaussians';
                    case 17
                        ShapeString='Fixed-position Lorentzians';
                    otherwise
                end % switch
                figure(1);subplot(2,1,2)
                xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
            end % if
        case 45 % if '-' key pressed
            disp(' ')
            disp('Select the peak shape of the model (type 1-17 and press Enter key):')
            disp('Gaussian: y=exp(-((x-pos)./(0.6005612.*width)) .^2)')
            disp('  Gaussians with independent positions and widths : 1 (default)')
            disp('  Exponentional-broadened Gaussian : 5 ')
            disp('  Exponentional-broadened equal-width Gaussian : 8')
            disp('  Gaussians with the same widths : 6')
            disp('  Gaussians with preset fixed widths : 11')
            disp('  Fixed-position Gaussians : 16 ')
            disp('  Asymmetrical Gaussians with unequal half-widths on both sides : 14')      
            disp('Lorentzian: y=ones(size(x))./(1+((x-pos)./(0.5.*width)).^2)')
            disp('  Lorentzians with independent positions and widths : 2')
            disp('  Equal-width Lorentzians : 7  ')
            disp('  Fixed-width Lorentzian : 12')
            disp('  Fixed-position Lorentzian : 17')
            disp('  Asymmetrical Lorentzians with unequal half-widths on both sides : 15')
            disp('Blended sum of Gaussian and Lorentzian functions : 13')
            disp('Logistic: n=exp(-((x-pos)/(.477.*wid)).^2); y=(2.*n)./(1+n) : 3  ')
            disp('Pearson: y=ones(size(x))./(1+((x-pos)./((0.5.^(2/m)).*wid)).^2).^m : 4')
            disp('Exponential pulse: y=exp(-tau1.*x).*(1-exp(-tau2.*x)) : 9')
            disp('Sigmoid: y=1/2 + 1/2* erf((x-t1)/sqrt(2*t2)) : 10')
            disp(' ')
            Shapeinput=input('Peak shape number: ');
            if isempty(Shape),
            else
                Shape=Shapeinput;
            end
            if Shape>17, Shape=17;end
            if Shape<1, Shape=1;end
            if isempty(Shape),Shape=1;end
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='logistic';
                case 4
                    ShapeString='Pearson';
                case 5
                    ShapeString='ExpGaussian';
                case 6
                    ShapeString='Equal-width Gaussian';
                case 7
                    ShapeString='Equal-width Lorentzian';
                case 8
                    ShapeString='Equal-width ExpGauss.';
                case 9
                    ShapeString='Exponental pulse';
                case 10
                    ShapeString='Sigmoid';
                case 11
                    ShapeString='Fixed-width Gaussian';
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 12
                    ShapeString='Fixed-width Lorentzian';
                    inputwidth=input('Peak width: ');
                    if isempty(inputwidth),
                    else
                        FIXEDPARAMETERS=inputwidth;
                    end
                case 13
                    ShapeString='Gauss/Lorentz blend';
                case 14
                    ShapeString='bifurcated Gaussian';
                case 15
                    ShapeString='bifurcated Lorentzian';
                case 16
                    ShapeString='Fixed-position Gaussians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                    end
                case 17
                    ShapeString='Fixed-position Lorentzians';
                    inputpositions=input('Peak positions as a vector, e.g. [200 400 600]: ');
                    if isempty(inputpositions),
                    else
                        FIXEDPARAMETERS=inputpositions;
                    end
                otherwise
            end % switch
            figure(1);subplot(2,1,2)
            xlabel(['Number of peaks = ' num2str(NumPeaks) '    Shape = ' ShapeString  ] )
        case 88 % Shift X
            extra=input('Value of extra parameter:');
            if Shape==4||5||13||14||15, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 97
            % When 'a' key is pressed, increases "extra" by 5%
            extra=extra+.05*extra;
            if extra==0, extra=.01;end
            if Shape==4||5||13||14||15, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 122
            % When 'z' key is pressed, decreases "extra" by 5%
            extra=extra-.05*extra;
            if extra==0, extra=.01;end
            if Shape==4||5||13||14||15, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 65
            % When 'A' key is pressed, increases "extra" by 0.5%
            extra=extra+.005*extra;
            if extra==0, extra=.001;end
            if Shape==4||5||13||14||15, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end
        case 90
            % When 'Z' key is pressed, decreases "extra" by 0.5%
            extra=extra-.005*extra;
            if extra==0, extra=.001;end
            if Shape==4||5||13||14||15, % Only do a re-fit if shape uses Extra parameter.
                [FitResults,MeanFitError]=FitAndPlot(xx,yy,NumPeaks,Shape,delta,start,extra);
            else
                [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
            end  
        case 114
            % When 'r' key is pressed, prints out table of fit results
            disp(['Fitting Error = ' num2str(MeanFitError) '%     Elapsed time = '  num2str(et) ' sec.' ])
            if Shape==9||Shape==10, % Pulse and Sigmoid only
                disp('         Peak#     Tau1         Height       Tau2          Area');
            else
                if Shape==0, % For future use
                    disp('         Peak#      Position       Height      Area');
                else
                    disp('         Peak#   Position       Height       Width         Area');
                end
            end
            disp(FitResults)
        case 100
            % When 'd' key is pressed, prints out table of model peak data
            disp('     x          y')
            disp([xx' yy'])
            switch NumPeaks,
                case 1
                    disp('     x           peak 1')
                    disp([xxx' PEAKHEIGHTS(1)*AA(1,:)' ])
                case 2
                    disp('     x           peak 1        peak 2')
                    disp([xxx' PEAKHEIGHTS(1)*AA(1,:)' PEAKHEIGHTS(2)*AA(2,:)'])
                case 3
                    disp('     x           peak 1        peak 2       peak 3')
                    disp([xxx' PEAKHEIGHTS(1)*AA(1,:)' PEAKHEIGHTS(2)*AA(2,:)' PEAKHEIGHTS(3)*AA(3,:)'])
                case 4
                    disp('     x           peak 1        peak 2       peak 3       peak 4')
                    disp([xxx' PEAKHEIGHTS(1)*AA(1,:)' PEAKHEIGHTS(2)*AA(2,:)' PEAKHEIGHTS(3)*AA(3,:)'  PEAKHEIGHTS(4)*AA(4,:)'])
                case 5
                    disp('     x           peak 1        peak 2       peak 3       peak 4       peak 5')
                    disp([xxx' PEAKHEIGHTS(1)*AA(1,:)' PEAKHEIGHTS(2)*AA(2,:)' PEAKHEIGHTS(3)*AA(3,:)'  PEAKHEIGHTS(4)*AA(4,:)'  PEAKHEIGHTS(5)*AA(5,:)'])
                case 6
                    disp('     x           peak 1        peak 2       peak 3       peak 4       peak 5       peak 6')
                    disp([xxx' PEAKHEIGHTS(1)*AA(1,:)' PEAKHEIGHTS(2)*AA(2,:)' PEAKHEIGHTS(3)*AA(3,:)'  PEAKHEIGHTS(4)*AA(4,:)'  PEAKHEIGHTS(5)*AA(5,:)'  PEAKHEIGHTS(6)*AA(6,:)'])
            end
        case 61 % When '=' key is pressed
            PEAKHEIGHKTS
        case 116
            % When 't' key is pressed, steps through AUTOZERO modes
            AUTOZERO=AUTOZERO+1;
            if AUTOZERO==3,AUTOZERO=0;end
            [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks);
        case 107
            % When 'k' key is pressed, prints out table of keyboar commands
            disp('KEYBOARD CONTROLS (Version 9.2):')
            disp(' Pan signal left and right...Coarse: < and >')
            disp('                             Fine: left and right cursor arrow keys')
            disp('                             Nudge: [ ] ')
            disp(' Zoom in and out.............Coarse zoom: / and "  ')
            disp('                             Fine zoom: up and down cursor arrow keys')
            disp(' Resets pan and zoom.........ESC')
            disp(' Select # of peaks...........Number keys 1-9, or press 0 key to enter manually')
            disp(' Select peak shape from menu - (minus or hyphen), then type number and Enter')
            disp(' Select peak shape by key....g Gaussian')
            disp('                             h equal-width Gaussians')
            disp('                             G (Shift-G) fixed-width Gaussians')
            disp('                             Shift P fixed-position Gaussians')
            disp('                             H (Shift-H) bifurcated Gaussians (a,z keys adjust shape)')
            disp('                             e Exponential-broadened Gaussian')
            disp('                             j exponential-broadened equal-width Gaussians')
            disp('                                 (a,z keys adjust broadening)')
            disp('                             l Lorentzian')
            disp('                             ; equal-width Lorentzians')
            disp('                             Shift [ fixed-position Lorentzians')
            disp('                             : bifurcated Lorentzians (a,z keys adjust shape)')
            disp('                             L (Shift-L) Fixed-width Lorentzians')
            disp('                             o LOgistic')
            disp('                             p Pearson (a,z keys adjust shape)')
            disp('                             u exponential pUlse  y=exp(-tau1.*x).*(1-exp(-tau2.*x))')
            disp('                             s Sigmoid: y=1/2 + 1/2*erf((x-t1)/sqrt(2*t2))')
            disp('                             ~ Gauss/Lorentz blend (a/z keys adjust shape)')
            disp(' Fit.........................f')
            disp(' Select autozero mode........t  selects none, linear, or quadratic autozero')
            disp(' Toggle log y mode OFF/ON....m  Plot linear or log Y axis on lower graph')
            disp(' 2-point Baseline............b, then click left and right baseline')
            disp(' Background subtraction......Backspace, then click baseline at multiple points')
            disp(' Restore original baseline...\  to cancel previous background subtraction')
            disp(' Cursor start positions......c, then click on each peak position')
            disp(' Enter value of ''Extra''......Shift-x, type value, press Enter.')
            disp(' Adjust ''Extra'' up/down......a,z: 5% change; upper case A,Z: 0.5% change.')
            disp(' Print parameters & results..q')
            disp(' Print fit results only......r')
            disp(' Compute bootstrap stats.....v  Estimates standard deViations of parameters.')
            disp(' Fit single bootstrap........n  Extracts and Fits siNgle bootstrap sub-sample.')
            disp(' Plot signal in figure 2.....y')
            disp(' Print model data table......d')
            disp(' Refine fit..................x Takes best of 10 trial fits (change in line 177)')
            disp(' Print peakfit function......w  Print peakfit function with all input arguments')
        case 120 % When 'x' key is pressed, calls peakfit to take best of 'NumTrials' trial fits
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            disp(['Best of ' num2str(NumTrials) ' trial fits.' ])
             if Shape==16||Shape==17, % Fixed-position shapes
                    fixedparameters=FIXEDPARAMETERS;
                else
                     fixedparameters=FIXEDPARAMETERS;
             end 
            tic
            [FitResults,MeanFitError]=peakfit([X',Y'],center,window,NumPeaks,Shape,extra,NumTrials,start,AUTOZERO,fixedparameters);
            et=toc;
            disp(['Percent Fitting Error ' num2str(MeanFitError) '%     Elapsed time = '  num2str(et) ' sec.' ])
            if Shape==9||Shape==10, % Pulse and Sigmoid only
                disp('         Peak#     Tau1         Height       Tau2          Area');
            else
                if Shape==0, % For future use
                    disp('         Peak#      Position       Height      Area');
                else
                    disp('         Peak#   Position       Height       Width         Area');
                end
            end
            disp(FitResults)
            figure(1)
        case 119
            % When 'W' is pressed, prints out peakfit functions with arguments in command window
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            FirstGuess=[];
            for peaknumber=1:NumPeaks,
                FirstGuess=[FirstGuess FitResults(peaknumber,2) FitResults(peaknumber,4)];
            end
             disp(['ipf(datamatrix,' num2str(center) ',' num2str(window) ')']);
            if Shape==11||12,
               disp(['[FitResults,MeanFitError]=peakfit(datamatrix,' num2str(center) ',' num2str(window) ',' num2str(NumPeaks) ',' num2str(Shape) ',' num2str(extra) ',' num2str(NumTrials) ', [' num2str(FirstGuess) '], ' num2str(AUTOZERO) ', ' num2str(FIXEDPARAMETERS) ' ) ' ]);
            else
               disp(['[FitResults,MeanFitError]=peakfit(datamatrix,' num2str(center) ',' num2str(window) ',' num2str(NumPeaks) ',' num2str(Shape) ',' num2str(extra) ',' num2str(NumTrials) ', [' num2str(FirstGuess) '], ' num2str(AUTOZERO) ' ) ']);
            end
         case 113
            % When 'q' key is pressed, prints out fitting parameters
            switch Shape
                case 1
                    ShapeString='Gaussian';
                case 2
                    ShapeString='Lorentzian';
                case 3
                    ShapeString='Logistic';
                case 4
                    ShapeString='Pearson';
                case 5
                    ShapeString='Exponentially-broadened Gaussian';
                case 6
                    ShapeString='Equal width Gaussians';
                case 7
                    ShapeString='Equal width Lorentzians';
                case 8
                    ShapeString='Equal-width ExpGauss.';
                case 9
                    ShapeString='Exponental pulse';
                case 10
                    ShapeString='Sigmoid';
                case 11
                    ShapeString='Fixed-width Gaussian';
                case 12
                    ShapeString='Fixed-width Lorentzian';
                case 13
                    ShapeString='Gauss/Lorentz blend';
                case 14
                    ShapeString='bifurcated Gaussian';
                case 15
                    ShapeString='bifurcated Lorentzian';
                case 16
                    ShapeString='Fixed-position Gaussians';
                case 17
                    ShapeString='Fixed-position Lorentzians';
                otherwise
            end % switch
            disp('------------------------------------------------------------------')
            disp(['Peak Shape = ' ShapeString])
            switch AUTOZERO,
                case 0
                    disp('Autozero OFF')
                case 1
                    disp('Linear autozero')
                case 2
                    disp('Quadratic autozero')
            end
            disp(['Number of peaks = ' num2str(NumPeaks)])
            if Shape==4||Shape==5||Shape==8||Shape==13||Shape==14||Shape==15, disp(['Extra = ' num2str(extra)]), end
            if Shape==11||Shape==12, disp(['Fixed Peak Width = ' num2str(FIXEDPARAMETERS)]), end
            disp(['Fitted x range = ' num2str(min(xx)) ' - ' num2str(max(xx)) ' (dx=' num2str(max(xx)-min(xx)) ')  (Center=' num2str((max(xx)+min(xx))/2) ')  ' ])
            disp(['Percent Fitting Error = ' num2str(MeanFitError) '%     Elapsed time = '  num2str(et) ' sec.' ])
            if Shape==9||Shape==10,
                disp('         Peak#     Tau1         Height       Tau2          Area');
            else
                disp('         Peak#   Position       Height       Width         Area');
            end
            disp(FitResults)
        case 121    % When 'Y' key is pressed (Added on version 5)
            figure(2) % Plot the entire signal cleanly in Figure window 2
            plot(X,Y)
            axis([X(1) X(length(X)) min(residual) max(Y)]);
            hold on
            for m=1:NumPeaks,
                % Add the individual component peaks in green lines
                plot(xxx,PEAKHEIGHTS(m)*AA(m,:),'g')
            end
            % Show residual in red
            plot(xx,residual,'r')
            hold off
            title('Blue=Original Data; Green=Computed Fit; Red=Residual (Fit-Data)')
            % figure(1)
        case 118 % When 'v' key is pressed (Added on version 8)
            NumTrialsBoot=input('Number of fit trials per bootstrap sample (0 to cancel): ');
            if isempty(NumTrialsBoot),NumTrialsBoot=1;end
            % EstTime=2.4+NumTrialsBoot.*(length(xx)./1436).*NumPeaks;
            EstTime=round((0.71581.*NumTrialsBoot).*(0.00069642.*length(xx)).*(5.4659.*NumPeaks)+2.5);
            % lengthxx=length(xx)
            disp(['Estimated time: ',num2str(EstTime) ' seconds']);
            if NumTrialsBoot,
                disp('Computing bootstrap sampling statistics....May take several minutes.')
                tic;
                BootstrapResultsMatrix=zeros(5,100,NumPeaks);
                BootstrapErrorMatrix=zeros(1,100,NumPeaks);
                center=(max(xx)+min(xx))/2;
                window=max(xx)-min(xx);
                clear bx by
                cutoff=0.5;
                for trial=1:100,
                    n=1;
                    bx=X';
                    by=Y';
                    while n<length(X)-1,
                        if rand>cutoff,
                            bx(n)=X(n+1);
                            by(n)=Y(n+1);
                        end
                        n=n+1;
                    end
                    [FitResults,BootFitError]=peakfit([bx,by],center,window,NumPeaks,Shape,extra,NumTrialsBoot,start,AUTOZERO,FIXEDPARAMETERS);
                    for peak=1:NumPeaks,
                        BootstrapResultsMatrix(:,trial,peak)=FitResults(peak,:);
                        BootstrapErrorMatrix(:,trial,peak)=BootFitError;
                    end
                end
                for peak=1:NumPeaks,
                    disp(' ')
                    if Shape==9||Shape==10, % Pulse and Sigmoid only
                        disp('         Peak#     Tau1         Height       Tau2          Area');
                    else
                        if Shape==0, % For future use
                            disp('         Peak#      Position       Height      Area');
                        else
                            disp('         Peak#   Position       Height       Width         Area');
                        end
                    end
                    BootstrapMean=mean(real(BootstrapResultsMatrix(:,:,peak)'));
                    BootstrapSTD=std(BootstrapResultsMatrix(:,:,peak)');
                    BootstrapIQR=iqr(BootstrapResultsMatrix(:,:,peak)');
                    PercentRSD=100.*BootstrapSTD./BootstrapMean;
                    PercentIQR=100.*BootstrapIQR./BootstrapMean;
                    MaxError=max(real(BootstrapErrorMatrix(:,:,peak)'));
                    MinError=min(real(BootstrapErrorMatrix(:,:,peak)'));
                    disp(['Bootstrap Mean: ', num2str(BootstrapMean(2:5))])
                    disp(['Bootstrap STD:  ', num2str(BootstrapSTD(2:5))])
                    disp(['Bootstrap IQR:  ', num2str(BootstrapIQR(2:5))])
                    disp(['Percent RSD:    ', num2str(PercentRSD(2:5))])
                    disp(['Percent IQR:    ', num2str(PercentIQR(2:5))])
                    % mean(PercentIQR(2:5)./PercentRSD(2:5))
                end
                toc;
                disp(['Min/Max Fitting Error of the 100 bootstrap samples: ', num2str(MaxError) '/' num2str(MinError) ])
                disp('-------------------------------------------------------------------')
            end
            figure(1)
            title('ipf 9.2  Typical Bootstrap sample fit')
        case 110 % When 'n' key is pressed (Added on version 8)
            disp('Fit to single bootstrap sample')
            tic;
            center=(max(xx)+min(xx))/2;
            window=max(xx)-min(xx);
            clear bx by
            cutoff=0.5;
            n=1;
            bx=X';
            by=Y';
            while n<length(X)-1,
                if rand>cutoff,
                    bx(n)=X(n+1);
                    by(n)=Y(n+1);
                end
                n=n+1;
            end
            [FitResults,MeanFitError]=peakfit([bx,by],center,window,NumPeaks,Shape,extra,NumTrialsBoot,start,AUTOZERO,FIXEDPARAMETERS);
            disp(['Fitting Error = ' num2str(MeanFitError) '%'])
            if Shape==9||Shape==10,
                disp('         Peak#     Tau1         Height       Tau2          Area');
            else
                disp('         Peak#   Position       Height       Width         Area');
            end
            disp(FitResults)
            figure(1)
            title('ipf 9.2  Single Bootstrap sample fit')
        otherwise
            UnassignedKey=double(key)
            disp('Press k to print out list of keyboard commands')
    end % switch
end % if

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--------------------------------------------------------------------------
function [xx,yy,start]=RedrawSignal(X,Y,xo,dx,NumPeaks)
% Plots the entire signal (X,Y) in the lower half of the plot window and an
% isolated segment (xx,yy) in the upper half, controlled by Pan and Zoom
% keys. Top half of the figure shows original signal
global AUTOZERO AA xxx PEAKHEIGHTS logplot
minX=min(X);maxX=max(X);minY=min(Y);maxY=max(Y);
Startx=round(xo-(dx/2));
Endx=abs(round(xo+(dx/2)-1));
if Endx>length(Y),Endx=length(Y);end
if Startx<1,Startx=1;end
PlotRange=Startx:Endx;
if (Endx-Startx)<2, PlotRange=xo:xo+2;end
xx=X(PlotRange);
yy=Y(PlotRange);
X1=min(xx);
X2=max(xx);
Y1=min(Y);
Y2=max(Y);

% Remove linear baseline from data segment if AUTOZERO==1
X1=min(xx);
X2=max(xx);
bkgsize=round(length(xx)/10);
if bkgsize<2,bkgsize=2;end
if AUTOZERO==1, % linear autozero operation
    Y1=mean(yy(1:bkgsize));
    Y2=mean(yy((length(xx)-bkgsize):length(xx)));
    yy=yy-((Y2-Y1)/(X2-X1)*(xx-X1)+Y1);
end % if

lxx=length(xx);
bkgsize=5;
if AUTOZERO==2, % Quadratic autozero operation
    XX1=xx(1:round(lxx/bkgsize));
    XX2=xx((lxx-round(lxx/bkgsize)):lxx);
    Y1=yy(1:round(length(xx)/bkgsize));
    Y2=yy((lxx-round(lxx/bkgsize)):lxx);
    bkgcoef=polyfit([XX1,XX2],[Y1,Y2],2);  % Fit parabola to sub-group of points
    bkg=polyval(bkgcoef,xx);
    yy=yy-bkg;
end % if autozero
hold off
clf
handles.uipanel1;subplot(2,1,1); % Select upper window
plot(xx,yy,'b.'); % Plot the original signal in blue in upper window
xlabel('Line up the estimated peak positions roughly with the vertical lines')
lyy=min(yy);
uyy=max(yy)+(max(yy)-min(yy))/10;
switch AUTOZERO,
    case 0
        title('ipf 9.2 Autozero OFF. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 1
        title('ipf 9.2 Linear autozero. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
    case 2
        title('ipf 9.2 Quadratic autozero. Pan and Zoom to isolate peaks to be fit in upper window.')
        if lyy<uyy;axis([X(Startx) X(Endx) lyy uyy ]);end
end
hold on
% Mark starting peak positions with vertical dashed lines in upper window
% Determine locations for peak (vertical line) markers
n=X2-X1;
width=n/(5*NumPeaks);
start=[];
startpos=[n/(NumPeaks+1):n/(NumPeaks+1):n-(n/(NumPeaks+1))]+X1;
for marker=1:NumPeaks,
    markx=startpos(marker);
    start=[start markx width];
    plot([markx markx],[lyy uyy],'m--')
end % for marker
hold off
%
% Bottom half of the figure shows full signal in either linear or log mode
% as set by the M key.
subplot(2,1,2);cla
if logplot,
    semilogy(X,abs(Y))  % Graph the signal with linear Y axis
    ylabel('Log y mode')
    axis([X(1) X(length(X)) min(abs(Y)) max(Y)]); % Update plot
else
    plot(X,Y)  % Graph the signal with linear Y axis
    ylabel('Linear y mode')
    axis([X(1) X(length(X)) min(Y) max(Y)]); % Update plot
end
title('# peaks 1-9,0   Shapes: g h G H e j l ; L : o p u s `   Fit=f   t=autozero   m=linear/log')
hold on
for marker=1:NumPeaks,
    markx=startpos(marker);
    plot([markx markx],[minY maxY],'m--')
end % for marker
% Added in version 5
for m=1:NumPeaks,
    lengthxxx=length(xxx);
    lengthmodel=length(PEAKHEIGHTS(m)*AA(m,:));
    if lengthxxx==lengthmodel,
        plot(xxx,PEAKHEIGHTS(m)*AA(m,:),'g')  % Plot the individual component peaks in green lines
    end
end
% Mark the limits of the upper windows on the lower whole-signal plot
plot([X1 X1],[minY maxY],'g--')
plot([X2 X2],[minY maxY],'g--')
hold off
xlabel(['Press K to print out all keyboard commands'])
% ----------------------------------------------------------------------
