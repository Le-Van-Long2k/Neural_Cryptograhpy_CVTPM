function varargout = Part1(varargin)
% PART1 MATLAB code for Part1.fig
%      PART1, by itself, creates a new PART1 or raises the existing
%      singleton*.
%
%      H = PART1 returns the handle to a new PART1 or the handle to
%      the existing singleton*.
%
%      PART1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PART1.M with the given input arguments.
%
%      PART1('Property','Value',...) creates a new PART1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Part1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Part1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Part1

% Last Modified by GUIDE v2.5 27-May-2021 20:38:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Part1_OpeningFcn, ...
                   'gui_OutputFcn',  @Part1_OutputFcn, ...
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


% --- Executes just before Part1 is made visible.
function Part1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Part1 (see VARARGIN)

% Choose default command line output for Part1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Part1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%axes_imgNeural = imread('E:\20202\Cryptography\BTL\CodeMatlabNeuralCryptography\Captur2525252e.PNG');
%imshow(axes_imgNeural);
% --- Outputs from this function are returned to the command line.
function varargout = Part1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global K N KxN Ta Tb Wa Wb S Key L X xma xmb ED RTa RTb RWa RWb RKey RX Rxma Rxmb RED;
KxN = 0;
Ta = 0;
Tb = 0;
Wa = 0;
Wb = 0;
S = 0;
Key = 0;
L = 0;
N = 0;
K = 0;
xma = 0;
xmb = 0;

RTa = 0;
RTb = 0;
RWa = 0;
RWb = 0;
RKey = 0;
Rxma = 0;
Rxmb = 0;
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_Chart.
function pushbutton_Chart_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Chart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global S ED xx RED;
title('\bfSynchronization CVTPM Plot');
xlabel('Synchronization time');
ylabel('ED');
xx = 1:1:S;
plot(xx,ED,'b',xx,RED,'r');
legend('CVTPM', 'TPM');

% --- Executes on button press in pushbutton_Init.
function pushbutton_Init_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Init (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global K N KxN Ta Tb Wa Wb S Key L X xma xmb ED RTa RTb RWa RWb RKey RX Rxma Rxmb RED;
ED(1,1) = 0;
RED(1,1) = 0;
N = get(handles.edit_N, 'string');
N = str2num(N);
K = get(handles.edit_K, 'string');
K = str2num(K);
KxN = K*N;
S = get(handles.edit_S, 'string');
S = str2num(S);
set(handles.edit_KxN, 'string', KxN);

L = get(handles.edit_L, 'string');
L = str2num(L);

Wa = weight(L,K,N);
RWa = Rweight(L,K,N);
wa = mat2str(Wa);
set(handles.edit_Wa, 'string', wa);

Wb = weight(L,K,N);
RWb = Rweight(L,K,N);
wb = mat2str(Wb);
set(handles.edit_Wb, 'string', wb);

    for i = 1:K
        for j = 1:N
            ED(1,1) = ED(1,1) + abs((Wa(i,j))-(Wb(i,j)));
            RED(1,1) = RED(1,1) + abs((RWa(i,j))-(RWb(i,j)));
        end
    end

set(handles.edit_ED, 'string', ED(1,1));
set(handles.edit_KeyA, 'string', mat2str(Wa));
set(handles.edit_KeyB, 'string', mat2str(Wb));
set(handles.edit_xma, 'string', mat2str(xma));
set(handles.edit_xmb, 'string', mat2str(xmb));
set(handles.edit_Ta, 'string', mat2str(Ta));
set(handles.edit_Tb, 'string', mat2str(Tb));

function edit_N_Callback(hObject, eventdata, handles)
% hObject    handle to edit_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_N as text
%        str2double(get(hObject,'String')) returns contents of edit_N as a double


% --- Executes during object creation, after setting all properties.
function edit_N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_K_Callback(hObject, eventdata, handles)
% hObject    handle to edit_K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_K as text
%        str2double(get(hObject,'String')) returns contents of edit_K as a double


% --- Executes during object creation, after setting all properties.
function edit_K_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_K (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_L_Callback(hObject, eventdata, handles)
% hObject    handle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_L as text
%        str2double(get(hObject,'String')) returns contents of edit_L as a double


% --- Executes during object creation, after setting all properties.
function edit_L_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_S_Callback(hObject, eventdata, handles)
% hObject    handle to edit_S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_S as text
%        str2double(get(hObject,'String')) returns contents of edit_S as a double


% --- Executes during object creation, after setting all properties.
function edit_S_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ED_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ED as text
%        str2double(get(hObject,'String')) returns contents of edit_ED as a double


% --- Executes during object creation, after setting all properties.
function edit_ED_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Wa_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Wa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Wa as text
%        str2double(get(hObject,'String')) returns contents of edit_Wa as a double


% --- Executes during object creation, after setting all properties.
function edit_Wa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Wa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Wb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Wb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Wb as text
%        str2double(get(hObject,'String')) returns contents of edit_Wb as a double


% --- Executes during object creation, after setting all properties.
function edit_Wb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Wb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_KeyA_Callback(hObject, eventdata, handles)
% hObject    handle to edit_KeyA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_KeyA as text
%        str2double(get(hObject,'String')) returns contents of edit_KeyA as a double


% --- Executes during object creation, after setting all properties.
function edit_KeyA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_KeyA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_KxN_Callback(hObject, eventdata, handles)
% hObject    handle to edit_KxN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_KxN as text
%        str2double(get(hObject,'String')) returns contents of edit_KxN as a double


% --- Executes during object creation, after setting all properties.
function edit_KxN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_KxN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_KeyB_Callback(hObject, eventdata, handles)
% hObject    handle to edit_KeyB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_KeyB as text
%        str2double(get(hObject,'String')) returns contents of edit_KeyB as a double


% --- Executes during object creation, after setting all properties.
function edit_KeyB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_KeyB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%cac ham cua CVTPM
function [W] = weight(L1,K1,N1)
% ham tao weight
K1 = int32(K1);
N1 = int32(N1);
L1 = int32(L1); 
for i=1:K1
   for j=1:N1
      W(i,j)= randi([-L1,L1]) + 1i*randi([-L1,L1]); 
   end
end

function [X] = xInput(K,N)
% ham tao matran X input
 for i=1:K
    for j=1:N
       X(i,j)= reIm() + 1i * reIm() ; 
    end
 end

function ab = reIm()
% ham tao -1 or 1  cho phan thuc, ao
    o = randi([0,1]);
    if o == 0
        ab = -1;
    else 
        ab = 1;
    end


function [xm] = Xm(N,K,W,X)
%Tinh h cua mot hidden unit thu kth


for i=1:K
    %Tinh phan thuc, ao
Re = 0;
Im = 0;
 for j=1:N  
    Re = Re + real(W(i,j)) * real(X(i,j));    
    Im = Im + imag(W(i,j)) * imag(X(i,j));    
 end
hk = 0;
%Tinh Hk
hk(i,1) = (Re + 1i*Im)/sqrt(N);
%   disp(hk);
thuc = sign(real(hk(i,1)));
ao = sign(imag(hk(i,1)));
if thuc == 0
    thuc = -1;
end
if ao ==0
    ao = -1;
end
xm(i,1) = thuc + 1i*ao;
end

function [T] = TinhTo(xm,K)
    Re = real(xm(1,1)) ;    
    Im = imag(xm(1,1)) ;  
 for i=2:K  
    Re = Re*real(xm(i,1)) ;    
    Im = Im*imag(xm(i,1)) ;    
 end

 T = Re + Im*1i;


 function [G] = g(omega,L)
%Guarantee that each elements of weight is in the range [-L,L]
%   g(omega,L) = sign(omega)*L if abs(omega)>L else g = omega
if abs(omega) > L
    G = sign(omega)*L;
else
    G = omega;
end

function [G] = Heaviside(a,b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if (a == b)
    G = 1;
else
    G = 0;
end

function [W] = update_weight(w,K,N,L,xm,to_1,to_2,X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
W = w;
if (real(to_1) == real(to_2) && imag(to_1) == imag(to_2))
    for i=1:K
        if (real(to_1) == real(xm(i,1))) && (imag(to_1) == imag(xm(i,1)))
            for j=1:N
                real_w = real(w(i,j));
                imag_w = imag(w(i,j));

                real_w = g(real_w + real(X(i,j))*real(to_1)*Heaviside(real(xm(i,1)), real(to_1))*Heaviside(real(to_1), real(to_2)),L);
                imag_w = g(imag_w + imag(X(i,j))*imag(to_1)*Heaviside(imag(xm(i,1)), imag(to_1))*Heaviside(imag(to_1), imag(to_2)),L);
                w(i,j) = real_w + 1i*imag_w;
            end
        elseif (real(to_1) == real(xm(i,1))) && (imag(to_1) ~= imag(xm(i,1)))
            for j=1:N
                real_w = real(w(i,j));
                imag_w = imag(w(i,j));
                real_w = g(real_w + real(X(i,j))*real(to_1)*Heaviside(real(xm(i,1)), real(to_1))*Heaviside(real(to_1), real(to_2)),L);
                w(i,j) = real_w + 1i*imag_w;
            end
        elseif (real(to_1) ~= real(xm(i,1))) && (imag(to_1) == imag(xm(i,1)))
            for j=1:N
                real_w = real(w(i,j));
                imag_w = imag(w(i,j));
              imag_w = g(imag_w + imag(X(i,j))*imag(to_1)*Heaviside(imag(xm(i,1)), imag(to_1))*Heaviside(imag(to_1), imag(to_2)),L);
                w(i,j) = real_w + 1i*imag_w;
            end
        end
    end
end
W = w;

%cac ham cua TPM
function [W] = Rweight(L1,K1,N1)
% ham tao weight
K1 = int32(K1);
N1 = int32(N1);
L1 = int32(L1); 
for i=1:K1
   for j=1:N1
      W(i,j)= randi([-L1,L1]); 
   end
end

function [X] = RxInput(K,N)
% ham tao matran X input
 for i = 1:K
    for j = 1:N
       X(i,j)= RreIm() ; 
    end
 end

function ab = RreIm()
% ham tao -1 or 1  cho phan thuc, ao
    o = randi([0,1]);
    if o == 0
        ab = -1;
    else 
        ab = 1;
    end


function [xm] = RXm(N,K,W,X)
%Tinh h cua mot hidden unit thu kth
% Kth = input('Nhap thu tu hidden unit: ');
xm = 0;
for i=1:K
    %Tinh phan thuc
Re = 0;

 for j=1:N  
    Re = Re + W(i,j) * X(i,j);        
 end
 
if Re == 0
    Re = -1;
end

xm(i) = (sign(Re)) ;
end

function [T] = RTinhTo(xm,K)
Re = real(xm(1));

 for i=2:K  
    Re = Re*real(xm(i)) ;      
 end

 T = Re ;
 
 function [G] = Rg(omega,L)
%Guarantee that each elements of weight is in the range [-L,L]
%   g(omega,L) = sign(omega)*L if abs(omega)>L else g = omega
if abs(omega) > L
    G = sign(omega)*L;
else
    G = omega;
end

function [G] = RHeaviside(a,b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if (a == b)
    G = 1;
else
    G = 0;
end

function [W] = Rupdate_weight(w,K,N,L,xm,to_1,to_2,X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
W = w;
if (to_1 == to_2)
    for i=1:K
        if (to_1 == xm(i))
            for j=1:N
                W(i,j) = Rg((W(i,j) + to_1*X(i,j)*RHeaviside(xm(i),to_1)*RHeaviside(to_1,to_2)),L) ;
            end
        end
    end
end


% --- Executes on button press in pushbutton_ED.
function pushbutton_ED_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global K N KxN Ta Tb Wa Wb S Key L X xma xmb ED EDx RTa RTb RWa RWb RKey RX Rxma Rxmb RED REDx;

for m = 1:S
    X = xInput(K,N);
    RX = RxInput(K,N);
    
    xma = Xm(N,K,Wa,X);
    xmb = Xm(N,K,Wb,X);
    Rxma = RXm(N,K,RWa,RX);
    Rxmb = RXm(N,K,RWb,RX);
    
    Ta = TinhTo(xma,K);
    Tb = TinhTo(xmb,K);
    RTa = RTinhTo(Rxma,K);
    RTb = RTinhTo(Rxmb,K);

   Wa = update_weight(Wa,K,N,L,xma,Ta,Tb,X);
   Wb = update_weight(Wb,K,N,L,xmb,Ta,Tb,X);
   RWa = Rupdate_weight(RWa,K,N,L,Rxma,RTa,RTb,RX);
   RWb = Rupdate_weight(RWb,K,N,L,Rxmb,RTa,RTb,RX);
   EDx = 0;
   REDx = 0;
   
    for i = 1:K
        for j = 1:N
            EDx = EDx + abs((Wa(i,j))-(Wb(i,j)));
            REDx = REDx + abs((RWa(i,j))-(RWb(i,j)));
        end
    end
    ED(m) = EDx;
    RED(m) = REDx;
    
    
set(handles.edit_ED, 'string', ED(m));  
set(handles.edit_KeyA, 'string', mat2str(Wa));
set(handles.edit_KeyB, 'string', mat2str(Wb));
set(handles.edit_X, 'string', mat2str(X));
set(handles.edit_xma, 'string', mat2str(xma));
set(handles.edit_xmb, 'string', mat2str(xmb));
set(handles.edit_Ta, 'string', mat2str(Ta));
set(handles.edit_Tb, 'string', mat2str(Tb));
pause(0.0001);
end


function edit_X_Callback(hObject, eventdata, handles)
% hObject    handle to edit_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_X as text
%        str2double(get(hObject,'String')) returns contents of edit_X as a double


% --- Executes during object creation, after setting all properties.
function edit_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xma as text
%        str2double(get(hObject,'String')) returns contents of edit_xma as a double


% --- Executes during object creation, after setting all properties.
function edit_xma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xmb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xmb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xmb as text
%        str2double(get(hObject,'String')) returns contents of edit_xmb as a double


% --- Executes during object creation, after setting all properties.
function edit_xmb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xmb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Ta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Ta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Ta as text
%        str2double(get(hObject,'String')) returns contents of edit_Ta as a double


% --- Executes during object creation, after setting all properties.
function edit_Ta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Ta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Tb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Tb as text
%        str2double(get(hObject,'String')) returns contents of edit_Tb as a double


% --- Executes during object creation, after setting all properties.
function edit_Tb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Tb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
