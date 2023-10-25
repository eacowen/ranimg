clear all
close all

% Adapted from ranimgcl.for which was used for Cowen & Monismith 1997
% Copyright Edwin A. Cowen October 24, 2023
% Note setting SVX=-SUY yields solid body rotation about the center point 
% (SUYC, SVXC)

% Set input parameters 
IFILESTART = 1;    % Starting FILE number
PREFIX = 'im_';    % 3 character prefix for image file name
NUMSTP = 10;       % Number of steps (image pairs)
NPERWN = 12;       % Number of particles per 32 x 32 subwindow
NLENX = 512;       % Number of columns in image
NLENY = 512;       % Number of rows in image
SIGPAR = 0.5;      % STDEV of particle image diameter (particle diameter is 4*SIGPAR
MEANNS = 50;       % Mean image background noise level in counts
SDNS = 2.;         % STDEV of image background noise level
BITDEPTH = 12;     % Bitdepth of modeled images (typically, 8, 10, 12, 14, 16)
LENPAR = 7;        % Length of particle array (square array)
INTPAR = 200;      % Mean maximum intensity of particles
ZMEAN = .05;       % Mean fraction of out of plane motion  (0. - 1.)
XTRANS = 10;       % Min translation in x direction
XTRSTP = 0.1;      % Translation in x direction step
YTRANS = 0.;       % Translation in y direction
SX = 0.0;          % dudx (displacement gradient)
SY = 0.0;          % dvdy (displacement gradient)
SVX = -0.0;        % dvdx (displacement gradient) 
SUY = 0.0;         % dudy (displacement gradient)
SXC = 256.;        % Center of dudx gradient in pixels
SYC = 256.;        % Center of dvdy gradient in pixels
SVXC = 256.;       % Center of dvdx gradient in pixels
SUYC = 256.;       % Center of dudy gradient in pixels

% Adjust the values for the physical location of the center based on my
% calculation method
SVXC = -SVX * SVXC;
SUYC = -SUY * SUYC;
SXC = -SX * SXC;
SYC = -SY * SYC;

% Estimate maximum displacement of a partcile
dispMaxX=ceil(abs(XTRANS)+max(abs(SUYC),abs(SXC)));
dispMaxY=ceil(abs(YTRANS)+max(abs(SVXC),abs(SYC)));

% Calculate total number of particles
NUMPAR = round((NLENX+2*dispMaxX) * (NLENY+2*dispMaxY) * NPERWN / 1024);
% Set figure colorbar limits
cmax=INTPAR+MEANNS+4*SDNS;
cmin=-cmax;

% Change BITDEPTH to maximum pixel intensity
BITDEPTH=2^BITDEPTH-1;

% Generate image pairs
for ILOOP = 1:NUMSTP
    % Initialize array of arguments for particle intensities
    % Generate random particle positions within the image boundaries
    X = rand(NUMPAR,1)*(NLENX+2*dispMaxX)-dispMaxX;  % Initialize X to invalid values
    Y = rand(NUMPAR,1)*(NLENY+2*dispMaxY)-dispMaxY;  % Initialize Y to invalid values

    X0=X; %Image 1 particle locations in x (j)
    Y0=Y; %Image 1 particle locations in y (i)
    X=X0-SUY*Y0+XTRANS-SUYC+SX*X0+SXC; %Image 2 particle locations in x (j)
    Y=Y0-SVX*X0-YTRANS-SVXC+SY*Y0+SYC; %Image 2 particle locations in y (i)

    % Remove particles impacted by out-of-plane motion
    IZOUT=round(ZMEAN*NUMPAR); % Number of particles to remove from image 2
    disp([num2str(IZOUT),' out of plane of ', num2str(NUMPAR),...
        ' particles in image greater region'])
    X(1:IZOUT)= rand(IZOUT,1)*NLENX;  % replace displaced partciles with
    Y(1:IZOUT) = rand(IZOUT,1)*NLENY; % random X and Y locations

    % Generate Images
    [I1,numParticles1] = GENIMG(NLENX,NLENY,SIGPAR,NUMPAR,X0,Y0,INTPAR,LENPAR,SDNS,MEANNS,BITDEPTH);
    [I2,numParticles2] = GENIMG(NLENX,NLENY,SIGPAR,NUMPAR,X,Y,INTPAR,LENPAR,SDNS,MEANNS,BITDEPTH);

    % Displace next image XTRANS by XTRSTP
    XTRANS=XTRANS+XTRSTP; 

    % Plot difference of particle intensities between Image 2 and Image 1    
    figure(1)
    D=I2-I1;
    set(gcf,"Position",[500 100 1000 1000])
    axis image
    imagesc(D,[cmin cmax])
    colormap("gray")
    colorbar
    set(gca,'fontsize',18)
    xlabel('$x$ or $j$ coordinate (pixels)','Interpreter','latex',FontSize=28)
    ylabel('$y$ or $i$ coordinate (pixels)','Interpreter','latex',FontSize=28)
    title('Image 2 - Image 1 (Black particles are image 1, white image 2')
    
    % Write image files 
    s1=sprintf( '%04d', IFILESTART); %image 1 number
    IFILESTART=IFILESTART+1;
    s2=sprintf( '%04d', IFILESTART); %image 2 number
    disp(['iteration ',num2str(ILOOP),'  x-displacement = ',num2str(XTRANS)...
        ,'  y-displacement = ', num2str(YTRANS)])
    disp(['Image 1 number of particles ' num2str(numParticles1)])
    disp(['Image 2 number of particles ' num2str(numParticles2)])
    disp(' ')
    fname1=[PREFIX s1]; %image 1 name
    fname2=[PREFIX s2]; %image 2 name
    save(fname1,"I1") % save image 1 to .mat
    save(fname2,"I2") % save image 2 to .mat
end 

function [image,numParticles]=GENIMG(NLENX,NLENY,SIGPAR,NUMPAR,X,Y,INTPAR,LENPAR,SDNS,MEANNS,BITDEPTH)
% Function to generate images of particles with integrated Gaussian pixel
% intensities
%
% image - 2D array to store the image
% NLENX, NLENY - Size of the image in x and y directions
% SIGPAR - Standard deviation of particle diameter
% NUMPAR - Number of particles
% X, Y - Particle positions
% INTPAR - Mean maximum intensity of particles
% LENPAR - Length of particle array
% SIG2 - Variance of background noise
% MEANNS - Mean value of background noise

% Initialize image with noise characteristics as specified
image=round(MEANNS+SDNS*randn([NLENY,NLENX]));

% Ensure the intensity is within a valid range
ind = find(image<0);
image(ind)=0;
ind = find(image>BITDEPTH);
image(ind)=BITDEPTH;

% Set up constants for integration of assumed Gaussian light scattering profile
FAC = 1 / (sqrt(2) * SIGPAR);
ARG = (erf(0.5 * FAC) - erf(-0.5 * FAC))^2;
FACINT = INTPAR / ARG;
% Define integration limits for particle region
LIMUP = (LENPAR - 1) / 2;
LIMLOW = -LIMUP;
% Determine particle iamge location releative to pixel array grid
IXBAR=round(X); % Neartest integer of pixel center in X
IYBAR=round(Y); % Neartest integer of pixel center in Y
XDIF=X-IXBAR; % dX of pixel center in X from IXBAR
YDIF=Y-IYBAR; % dY of pixel center in Y from IYBAR
ARGX1 = 0.5 - XDIF; 
ARGX2 = -0.5 - XDIF;
ARGY1 = 0.5 - YDIF;
ARGY2 = -0.5 - YDIF;

% Integrate Gaussian image partcile of appropriate diameter and replace
% background noise with it if particle intensity + noise is greater than
% noise
for I = 1:NUMPAR
    for J = LIMLOW:LIMUP
        ERFCTJ = erf((J + ARGX1(I)) * FAC) - erf((J + ARGX2(I)) * FAC);
        for K = LIMLOW:LIMUP
            ERFCTK = erf((K + ARGY1(I)) * FAC) - erf((K + ARGY2(I)) * FAC);
            BKNOIS = MEANNS+SDNS*randn();
            INTENS = FACINT * ERFCTJ * ERFCTK + BKNOIS;
            if INTENS < 0
                INTENS = 0;
            elseif INTENS > BITDEPTH
                INTENS = BITDEPTH;
            end
            IX = IXBAR(I) + J;
            IY = IYBAR(I) + K;
            % Verify particle coordinates are within image and update
            % intensity
            if IX >= 1 && IX <= NLENX && IY >= 1 && IY <= NLENY
                if image(IY, IX) < INTENS
                    image(IY, IX) = INTENS;
                end
            end
        end
    end
end
isInImage=(X >= 0.5) .* (X <= NLENX+0.5) .* (Y >= 0.5) .* (Y <= NLENY+0.5);
ind=find(isInImage==1);
x=X(ind);
y=Y(ind);
% save the exact coordinates of images with their centers considered within
% the image based on 0.5 <-> NLENX/y+0.5.
% Note this file overwrites each image generation so it is only the last
% image written (final I2) unless this is modified.
save exactXY x y
numParticles=sum(isInImage);
end



