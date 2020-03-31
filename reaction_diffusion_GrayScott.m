function reaction_diffusion_GrayScott(gridSize,Du,Dv,f,k,dt,initState,numIterations,frameRate)
% Generates Gray-Scott Reaction Diffusion Model, Exports Video
% gridSize: grid size. Standard values are 100, 200, 400
% Du, Dv: diffusion coefficients. Standard values are 1, 0.5
% f, k: feed and kill rates. Standard values are between 0.01, 0.1
% dt: time increment between iterations. Standard value of 1
% initState: init reaction state. "square", "bar", "wavefront" accepted
% numIterations: number of iterations. Standard values between 5000, 25000
% frameRate: number of iterations between frames. Standard value of 25

% Default cell values
u = ones(gridSize);
v = zeros(gridSize);

if (initState == "wavefront")
    % Wavefront initial state
    u(gridSize/2-30:gridSize/2+30,gridSize/2-6:gridSize/2-4) = 0.75;
    u(gridSize/2-30:gridSize/2+30,gridSize/2-3:gridSize/2-1) = 0.5;
    u(gridSize/2-30:gridSize/2+30,gridSize/2:gridSize/2+2) = 0.25;
    u(gridSize/2-30:gridSize/2+30,gridSize/2+3:gridSize/2+5) = 0;
    v(gridSize/2-30:gridSize/2+30,gridSize/2+3:gridSize/2+5) = 1;
elseif (initState == "bar") 
    % Bar initial state
    u(gridSize/2-10:gridSize/2+10,gridSize/2-2:gridSize/2+2) = 0;
    v(gridSize/2-10:gridSize/2+10,gridSize/2-2:gridSize/2+2) = 1;
else 
    % Square initial state
    u(gridSize/2-5:gridSize/2+5,gridSize/2-5:gridSize/2+5) = 0;
    v(gridSize/2-5:gridSize/2+5,gridSize/2-5:gridSize/2+5) = 1;
end

% Sets colormap
colormap jet

% Defines VideoWriter object
% writer = VideoWriter('rxn_dfsn_gs.avi');
writer = VideoWriter(strcat('rxn_dfsn_gs_',num2str(f),'_',num2str(k),'.avi'));
open(writer);

% Opens figure for image plotting
figure(1)
% Creates image based off of u concentrations
image(u,'CDataMapping','scaled');
% Grabs frame
frame = getframe(1);
% Writes frame to video
writeVideo(writer,frame);

% Iterates model
for i = 1:numIterations
   % Calculates laplacian
   uLaplace = laplacian(u);
   vLaplace = laplacian(v);
   % Calculates instantaneous rate of change
   duT = Du.*(uLaplace)-u.*v.^2+f.*(1.-u);
   dvT = Dv.*(vLaplace)+u.*v.^2-(f+k).*v;
   % Increments values using Euler's method
   u = u + dt*duT;
   v = v + dt*dvT;
   % Records frame every "frameRate" iterations
   if (mod(i,frameRate) == 0)
       image(u,'CDataMapping','scaled');
       frame = getframe(1);
       writeVideo(writer,frame);
   end
end

% Closes videoWriter and downloads video
close(writer);

fprintf('Done!\n');

end