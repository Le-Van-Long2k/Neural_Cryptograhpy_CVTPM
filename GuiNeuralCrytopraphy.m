function varargout = GuiNeuralCrytopraphy(varargin)
% GUINEURALCRYTOPRAPHY MATLAB code for GuiNeuralCrytopraphy.fig
%      GUINEURALCRYTOPRAPHY, by itself, creates a new GUINEURALCRYTOPRAPHY or raises the existing
%      singleton*.
%
%      H = GUINEURALCRYTOPRAPHY returns the handle to a new GUINEURALCRYTOPRAPHY or the handle to
%      the existing singleton*.
%
%      GUINEURALCRYTOPRAPHY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUINEURALCRYTOPRAPHY.M with the given input arguments.
%
%      GUINEURALCRYTOPRAPHY('Property','Value',...) creates a new GUINEURALCRYTOPRAPHY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GuiNeuralCrytopraphy_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GuiNeuralCrytopraphy_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GuiNeuralCrytopraphy

% Last Modified by GUIDE v2.5 20-May-2021 15:35:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GuiNeuralCrytopraphy_OpeningFcn, ...
                   'gui_OutputFcn',  @GuiNeuralCrytopraphy_OutputFcn, ...
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


% --- Executes just before GuiNeuralCrytopraphy is made visible.
function GuiNeuralCrytopraphy_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GuiNeuralCrytopraphy (see VARARGIN)

% Choose default command line output for GuiNeuralCrytopraphy
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GuiNeuralCrytopraphy wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GuiNeuralCrytopraphy_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function img_neural_CreateFcn(hObject, eventdata, handles)
% hObject    handle to img_neural (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate img_neural
img_neural = imread('E:\20202\Cryptography\BTL\CodeMatlabNeuralCryptography\Captur2525252e.PNG');
imshow(img_neural);
