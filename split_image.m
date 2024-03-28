function [mask, vis] = split_image(imo)
% splits a segmented image (for example an image that has been processed by
% ilastik already) into two approximately symmetric regions
% by returning a mask.
% Assumes background pixels to be <=1, and all other labels being
% higher in value than 1.
%
% CAUTION: This function makes assumptions about the content of the image
% being symmetric and that only one large region gets split.
% This function will fail to produce anything useful if images do not have
% symmetric content.
%
% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vis = [];
mask = [];
im = imo;
%binarize
im(im<=1) = 0;
im(im>0) = 255;

%aggessive imclose and keep largest component
im = imclose(im, strel('disk', 10));
BW = bwareafilt(logical(im), 1);  % get the largest area

%produce convex hull and calculate resulting object properties --- we
%should also experiment with just calculating properties directly
%BW = bwconvhull(BW);
s  = regionprops(BW, 'Orientation', 'MajorAxisLength', ...
    'MinorAxisLength', 'Eccentricity', 'Centroid');


%% 
centroid = s.Centroid;
% ellipse outline
phi = linspace(0,2*pi,500);
cosphi = cos(phi);
sinphi = sin(phi);
xbar = s.Centroid(1);
ybar = s.Centroid(2);
a = s.MajorAxisLength/2;
b = s.MinorAxisLength/2;
theta = pi*s.Orientation/180;
R = [ cos(theta)   sin(theta)
    -sin(theta)   cos(theta)];
xy = [a*cosphi; b*sinphi];
xy = R*xy;
x = xy(1,:) + xbar;
y = xy(2,:) + ybar;

% generate the dividing line equation
d = sqrt((xbar-x(:)).^2 + (ybar-y(:)).^2);
dind = find(d==min(d), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Use this if not registering
mind = dind;
X = [xbar x(mind)];
Y = [ybar y(mind)];
coefficients = polyfit(X, Y, 1);
a = coefficients (1);
b = coefficients (2);
X = [1 size(imo,2)];
Y = a.*X + b;
% generate the dividing mask
Xm = [X size(imo,1) 1];
Ym = [Y 1 1];
mask = poly2mask(Xm, Ym, size(imo,1), size(imo,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % %%%% Experiment: Keep this commented out
% % % We could register images over a range of angles and use the score to 
% % % assess asymmetry.
% % vec = -20:5:20;
% % dind_vec = zeros(length(vec), 2);    
% % for mix = -20:5:20
% %     mind = dind+mix;
% %     X = [xbar x(mind)];
% %     Y = [ybar y(mind)];
% %     coefficients = polyfit(X, Y, 1);
% %     a = coefficients (1);
% %     b = coefficients (2);
% %     X = [1 size(imo,2)];
% %     Y = a.*X + b;
% %     % generate the dividing mask
% %     Xm = [X size(imo,1) 1];
% %     Ym = [Y 1 1];
% %     mask = poly2mask(Xm, Ym, size(imo,1), size(imo,2));
% %     
% %     im1 = double(imo).*double(mask);
% %     im2 = double(imo).*double(~mask);
% %     score = register_pair(im1, im2);
% %     dind_vec(mix,:) = [score mind];
% % end

if nargout>1 % then we need to return a quality control figure
    h = figure('visible', 'off');
    hold on;
    imshow(mat2gray(imo+uint8(mask)));
    hold on;
    plot(x,y,'r','LineWidth',2);
    line(X,Y); % plot the line
    plot(centroid(:,1), centroid(:,2), 'b*'); % plot centroid
    hold off
    vis = getframe(gca);
    vis = vis.cdata;
    close;
end
close;


%% below function not currently used. ----- experimental
function score = register_pair(fixed, moving)

imshowpair(fixed, moving,'Scaling','joint');

tformEstimate = imregcorr(moving,fixed);
Rfixed = imref2d(size(fixed));
movingReg = imwarp(moving,tformEstimate,'OutputView',Rfixed);

imshowpair(fixed,movingReg,'montage')



%% Below didn't work on first try --- would need to be optimized
[optimizer, metric] = imregconfig('multimodal');

optimizer.InitialRadius = 0.009;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 300;

movingRegistered = imregister(moving, fixed, 'affine', optimizer, metric);


imshowpair(fixed, movingRegistered,'Scaling','joint');
imshowpair(moving, movingRegistered,'Scaling','joint');







