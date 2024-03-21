# Fourier-decomposition
% Demonstration of fourier decomposition. 
% Goal: to visualize that (2D) images can be decomposed 
% into sinusoidal waves. 
% ----------
% Implementation: At each 'iteration' the most prominent 
% wave* is taken away from the input image, and set aside.
% Where the residual image gradually loses contrast, the removed 
% waves together start forming a new reconstructed 
% image (initially low pass filtered). After 20 
% waves the rest of the waves is removed in groups of 
% increasing size.
% (*together with its 'conjugate' wave). 
% ----------
% Author: Job Bouwman
% date: 20 july 2013 
% ----------

% Demonstration of fourier decomposition, showing both the image domain
% and spectrum domain
% modified from https://www.mathworks.com/matlabcentral/fileexchange/42776-fourier-decomposition-demo
% add the spectrum domain and the code comments

% ----------
% Contributor: Sun
% date: 21 March 2024 
% ----------
