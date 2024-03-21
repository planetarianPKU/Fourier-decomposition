%clear all wins and vars,
close all; clear; clc;

%read a image in the directionary
I = imread("H:\Crips\Mona_Lisa,_by_Leonardo_da_Vinci,_from_C2RMF_retouched.jpg");

I = rgb2gray(I);
%Display the original image
imshow(I);

% Compute the Fourier transform of the image
FT_I = fftn(I);

% Shift the zero-frequency component to the center of the spectrum
FT_I_shifted = fftshift(FT_I);

% Compute the log magnitude for better visualization
magnitude = log( abs(FT_I_shifted));

% Display the adjusted magnitude imag
imshow(magnitude/max(max(magnitude))*1);
%%
% Set the scale factor
scale = 0.25;   % an be {1, 0.5, or 0.25} to ensure the size of image is integer

% Resize the image
I = imresize(double(I), scale);

% Calculate the mean of the image
I_mean = mean(I(:));

% Center the image by subtracting its mean
%Because we focus on the change of image, not it average value
I = I - I_mean;

% Get the size of the image
N = size(I);

% Find the minimum and maximum values in the image
minI = min(I(:));
maxI = max(I(:));

% Set the viewing range for the image display
viewRange = [(1.25*minI - 0.25*maxI), (1.25*maxI - 0.25*minI)];

% Set the speed for the demonstration
demoSpeed = 20;
%%
% Configure the figure window size and background color
figure
set(gcf,'unit','normalized','position',[0,0,1,1]);
set(gcf,'color','white');

% Recompute the Fourier transform for the centered image
FT_I = fft2(I);
FT_I_center = fftshift(FT_I);

% Initialize arrays for the Fourier transform and its centered version
FT_new = zeros(size(FT_I));
FT_new_center = zeros(size(FT_I_center));

%Magnitude of 2-D spectrum
ABS_FT_I = abs(FT_I);
ABS_FT_I_center = abs(FT_I_center);

%spectrum of the extracted Nth group of wavenumbers
FT_cur = zeros(size(FT_new));
%accumulated image
FT_cur_center = zeros(size(FT_new_center));
%extracted image
FT_cur_center_now = zeros(size(FT_new_center));

%image of the extracted Nth group of wavenumbers
newIm = zeros(size(FT_new));
newFT_center = zeros(size(FT_new_center));

%Set the viewing range for the spectrum image display
maxFT = max(ABS_FT_I_center(:));
viewRange_FT = [0 log(maxFT+1)];

%%
n = 1; nWaves = 1;
%Begin processing each frequency component
while n < numel(I)
    % Logic to increment the wave number in each loop
    FT_cur = 0*FT_cur;
    %Extract_this_step
    FT_cur_center_now = 0*FT_cur_center_now;
    if n > 20;  nWaves = 10 ; end
    if n > 200; nWaves = 100; end
    if n > 2000; nWaves = 1000; end
    
    % Select the strongest frequency components for processing,
    % in other words, low-frequency energy of the current image
    for p = 1:nWaves*2
        [a,b] = find(ABS_FT_I == max(ABS_FT_I(:)), 1, 'first');
        ABS_FT_I(a,b) = 0;
        FT_cur(a,b) = FT_I(a,b);
    end
    % Same as above, but for the centered spectrum to display the spectrum
    for p = 1:nWaves*2
        [a,b] = find(ABS_FT_I_center == max(ABS_FT_I_center(:)), 1, 'first');
        ABS_FT_I_center(a,b) = 0;
        FT_cur_center_now(a,b) = FT_I_center(a,b);
        FT_cur_center(a,b) = FT_I_center(a,b);
    end
    
    
    % Perform the inverse Fourier transform on the selected frequency components
    I_cur = ifft2(FT_cur);
    
    % Prepare the canvas for the spectrum and image reconstruction display
    canvas_FT = (cat(2,(ABS_FT_I_center),log(maxFT+1)*2/3*ones(size(I)),(abs(FT_cur_center))));
    canvasShow_FT = canvas_FT;
    canvasShow_FT(:,1:N(2)) = canvasShow_FT(:,1:N(2)) + ((FT_cur_center_now));
    
    % Note that the spectrum canvasShow_FT is transferred to
    subplot(2,1,1); imagesc(log(abs(canvasShow_FT)+1), viewRange_FT); axis equal tight; colormap(gray);
    
    % Prepare the canvas for image reconstruction display
    canvas = cat(2, real(I - newIm - I_cur), zeros(size(I)), newIm);
    canvasShow = canvas;
    canvasShow(:,1:N(2)) = canvasShow(:,1:N(2)) + I_cur;
    subplot(2,1,2); imagesc(real(canvasShow), viewRange); axis equal tight; colormap gray;
    
    
    % Remaining loop and plotting logic continues...
    
    % Determine the title based on the number of waves processed
    if n < 20
        title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
    else
        title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);
    end
    % Set the x-axis label for the subplot showing the images
    
    xlabel('residual image (left)                  extracted wave (middle)                 accumulated image (right)', 'FontSize',12);
    % Pause for a certain time, which decreases as n increases, to control the demonstration speed
    
    pause((1/demoSpeed)*exp(-n/2));
    % Iterate through a set of positions for displaying the extracted wave in the middle section
    for L = [(min(max(1, round((1 + sin(linspace(-pi/2, pi/2, 40/min(20, n))))/2*N(2))), N(2)+1))   (N(2)+1)]
        % Update the Fourier transform canvas for the current wave
        
        canvasShow_FT = canvas_FT;
        canvasShow_FT(:,L:(L+(N(2)-1))) = canvasShow_FT(:,L:(L+(N(2)-1))) + ((FT_cur_center_now));
        subplot(2,1,1); imagesc(log(abs(canvasShow_FT)+1), viewRange_FT); axis equal tight; colormap gray;
        
        % Display the updated Fourier transform canvas
        canvasShow = canvas;
        canvasShow(:,L:(L+(N(2)-1))) = canvasShow(:,L:(L+(N(2)-1))) + I_cur;
        subplot(2,1,2); imagesc(real(canvasShow), viewRange); axis equal tight; colormap gray;
        % Update the title for the image display
        if n < 20
            title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
        else
            title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);
        end
        xlabel('residual image (left)                  extracted wave (middle)                 accumulated image (right)', 'FontSize',12);
        pause((1/demoSpeed)*exp(-(n+4)));
        drawnow
    end
    if n < 20
        title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
    else
        title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);
    end
    xlabel('residual image (left)             extracted wave (middle)             accumulated image (right)', 'FontSize',12);
    pause((1/demoSpeed)*exp(-n/2));
    
    % Repeat similar steps as above for a different set of positions in the image
    % This block is similar to the previous one, updating and displaying the changes
    % in the Fourier transform and image canvases at different positions
    % It ensures that the visualization smoothly transitions across the image
    for L = N(2) + [(min(max(1, round((1 + sin(linspace(-pi/2, pi/2, 40/min(20, n))))/2*N(2))), N(2)+1))   (N(2)+1)]
        
        canvasShow_FT = canvas_FT;
        canvasShow_FT(:,L:(L+(N(2)-1))) = canvasShow_FT(:,L:(L+(N(2)-1))) + ((FT_cur_center_now));
        subplot(2,1,1); imagesc(log(abs(canvasShow_FT)+1), viewRange_FT); axis equal tight; colormap gray;
        
        canvasShow = canvas;
        canvasShow(:,L:(L+(N(2)-1))) = canvasShow(:,L:(L+(N(2)-1))) + I_cur;
        subplot(2,1,2); imagesc(real(canvasShow), viewRange); axis equal tight; colormap gray;
        if n < 20
            title(strcat('Wave number: ', num2str(n)), 'FontSize',14);
        else
            title(strcat('Wave numbers: ', num2str(n), ' - ', num2str(min(numel(I),n+nWaves-1))), 'FontSize',14);
        end
        xlabel('residual image (left)                  extracted wave (middle)                 accumulated image (right)', 'FontSize',12);
        pause((1/demoSpeed)*exp(-(n+4)));
        drawnow
    end
    % Add the current image part to the accumulated image
    newIm = newIm + I_cur;
    % Add the current Fourier part to the accumulated Fourier transform
    newFT_center = newFT_center + FT_cur_center;
    % Increment the wave counter by the number of waves processed
    n = n + nWaves;
end

