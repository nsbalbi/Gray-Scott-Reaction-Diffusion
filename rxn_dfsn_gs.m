function rxn_dfsn_gs(f,k,Du,Dv,options)
% RXN_DFSN_GS Generates and records Gray-Scott reaction diffusion
% models using Euler's method.
%
%    RXN_DFSN_GS() generates a Gray-Scott reaction-diffusion
%    model with sample coefficients. The result is saved as an MP4 named 
%    'rxn_dfsn_gs.mp4'.
%
%    RXN_DFSN_GS(f,k) generates a Gray-Scott reaction-diffusion
%    model with standard diffusion coefficients, a rate of conversion k,
%    and a feed rate f. The result is saved as an MP4 named 
%    'rxn_dfsn_gs.mp4'. See "Inputs" for details.
%
%    RXN_DFSN_GS(f,k,Du,Dv) generates a Gray-Scott reaction-diffusion
%    model with the diffusion coefficients Du and Dv, rate of conversion k,
%    and feed rate f. The result is saved as an MP4 named 
%    'rxn_dfsn_gs.mp4'. See "Inputs" for details. 
% 
%    RXN_DFSN_GS(...,PARAM1,VAL1,PARAM2,VAL2,...) are additional
%    name-value pairs that can be used to change default function
%    parameters. See "Parameters" for details.
%
% Inputs
%    Du,Dv: double. Diffusion coefficients for species u and v. All 
%       equation constants are based off of the standard set of equations 
%       for Gray-Scott systems. Default values are 1 and 0.5, respectively.
%    f: Feed rate. Recommended value between 0 and 0.1. Default is 0.04 
%    k: Rate of conversion. Recommended value between 0 and 0.1. Default is
%       0.0636
%
% Parameters
%    'GridSize': int. Edge length for the square grid used for the
%       model. Note generation time increases exponentially with grid size.
%       All initial states work with a grid size >= 70. (Default = 200)
%    'InitState': char. Intial model state. Possible inputs listed below.
%          'square'- Centered square concentration spot of species v. (Default)
%       'wavefront'- Simulated wave front of species v.
%             'bar'- Centered rectangular concentration spot of species v.
%          'vSpots'- 10 random square concentration spots of species v.
%          'uSpots'- 10 random square concentration spots of species u.
%    'NumIterations': int. Total number of model iterations. Each iteration
%       represents a step of Euler's method. (Default = 10000)
%    'FrameSpacing': int. Number of iterations between captured frames 
%       (Default = 20)
%    'dt': double. Time increment between iterations of Euler's method. 
%       (Default = 1)
%    'Colormap': char. Colormap used for visualization. Supports all 
%       default MATLAB colormaps (Default = 'jet')
%    'FileName': char. Name and file format of created video. See the
%       MATLAB help site for a list of supported video export formats.
%       (Default = 'rxn_dfsn_gs.mp4')
%
%
% Examples
%    Example 1: Generate and record default model.
%       rxn_dfsn_gs();
%    Example 2: Generate and record model with custom f and k constants.
%       rxn_dfsn_gs(0.022,0.051);
%    Example 3: Generate and record model with custom constants.
%       rxn_dfsn_gs(0.022,0.055,1,0.3);
%    Example 4: Generate and record model with custom f,k and initial state.
%       rxn_dfsn_gs(0.01,0.041,'InitState','wavefront');
%    Example 5: Generate and record model with custom constants and
%    multiple custom parameters.
%       rxn_dfsn_gs(0.09,0.059,1,0.45,...
%                   'InitState','uSpots',...
%                   'FrameSpacing',60,...
%                   'NumIterations',20000);       

arguments
   f double = 0.04
   k double = 0.0636
   Du double = 1
   Dv double = 0.5
   options.GridSize int64 = 200
   options.InitState char = 'square'
   options.NumIterations int64 = 10000
   options.FrameSpacing int64 = 20
   options.dt double = 1
   options.Colormap char = 'jet'
   options.FileName char = 'rxn_dfsn_gs.mp4'
end

%% Initialization

% Default cell values
u = ones(options.GridSize);
v = zeros(options.GridSize);
if (strcmp(options.InitState,'square')) % Square initial state
    u(options.GridSize/2-5:options.GridSize/2+5,options.GridSize/2-5:options.GridSize/2+5) = 0;
    v(options.GridSize/2-5:options.GridSize/2+5,options.GridSize/2-5:options.GridSize/2+5) = 1;
elseif (strcmp(options.InitState,'wavefront')) % Wavefront initial state
    u(options.GridSize/2-30:options.GridSize/2+30,options.GridSize/2-6:options.GridSize/2-4) = 0.75;
    u(options.GridSize/2-30:options.GridSize/2+30,options.GridSize/2-3:options.GridSize/2-1) = 0.5;
    u(options.GridSize/2-30:options.GridSize/2+30,options.GridSize/2:options.GridSize/2+2) = 0.25;
    u(options.GridSize/2-30:options.GridSize/2+30,options.GridSize/2+3:options.GridSize/2+5) = 0;
    v(options.GridSize/2-30:options.GridSize/2+30,options.GridSize/2+3:options.GridSize/2+5) = 1;
elseif (strcmp(options.InitState,'bar')) % Bar initial state
    u(options.GridSize/2-10:options.GridSize/2+10,options.GridSize/2-2:options.GridSize/2+2) = 0;
    v(options.GridSize/2-10:options.GridSize/2+10,options.GridSize/2-2:options.GridSize/2+2) = 1;
elseif (strcmp(options.InitState,'vSpots')) % Random spots of v initial state
    for i = 1:10
        xrand = randi([5,options.GridSize-5]);
        yrand = randi([5,options.GridSize-5]);
        v(xrand-4:xrand+4,yrand-4:yrand+4) = 1;
        u(xrand-4:xrand+4,yrand-4:yrand+4) = 0;
    end
elseif (strcmp(options.InitState,'uSpots')) % Random spots of u initial state
    u = zeros(options.GridSize);
    v = ones(options.GridSize);
    for i = 1:10
        xrand = randi([5,options.GridSize-5]);
        yrand = randi([5,options.GridSize-5]);
        u(xrand-4:xrand+4,yrand-4:yrand+4) = 1;
        v(xrand-4:xrand+4,yrand-4:yrand+4) = 0;
    end
end

writer = VideoWriter(options.FileName); % Define VideoWriter object
open(writer);
figure(1) % Open figure for image plotting
colormap(options.Colormap) % Set colormap

%% Iterative Modeling

% Model iteration loop
for i = 1:options.NumIterations
   % Calculates instantaneous rate of change
   duT = Du.*(laplacian(u))-u.*v.^2+f.*(1.-u);
   dvT = Dv.*(laplacian(v))+u.*v.^2-(f+k).*v;
   % Increments values using Euler's method
   u = u + options.dt*duT;
   v = v + options.dt*dvT;
   % Records frame every 'FrameSpacing' iterations
   if (mod(i,options.FrameSpacing) == 0)
       image(u,'CDataMapping','scaled'); % Create image based from u conc
       frame = getframe(1); % Grab frame
       writeVideo(writer,frame); % Write frame to video
   end
end

close(writer); % Close videoWriter
fprintf('Done!\n'); % Notify user when recording is complete

end
