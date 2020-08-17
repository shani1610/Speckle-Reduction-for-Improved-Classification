function Features = FeaturesExtraction(I)
%% Features Extraction
% calculating numeric values of features of Ultrasound breast lesions.
% Inbal Aharoni and Shani Israelov, Technion 2020.

show_images = false; 
%figure(4);

%% Preprocessing
% Read the desired image  

%I_crop = rgb2gray(I);
%I2 = imread('breast1_norm_croped.png');

%I2 = imread('FI1.png');
%I_crop = rgb2gray(I2);
I_crop=I;

if show_images
    figure(1);
    imshow(I_crop);
    F(1) = getframe; % For animated image
end
[row_size, col_size] = size(I_crop);

% Segmantation
bw = I_crop > 100;
bw = bwareaopen(bw,1000);
bw=~bw;
bw = bwareaopen(bw,800); % Remove all object containing fewer than 800 pixels
bw2 = imfill(bw, 'holes');
bw2 = imclearborder(bw2); % Clear objects that are close to the bonderies 

% Creating one object from two 
SE1 = strel('disk',6);
mask1 = imdilate(bw2,SE1);
if show_images 
    figure(2);
    imshow(mask1);
    title('After imdilate');
end 
mask2 = imfill(mask1, 'holes');
SE1 = strel('disk',6);
mask3 = imerode(mask2, SE1);
% if show_images 
%     figure(3);
%     imshow(mask3);
%     title('After imerode');
% end 
bw = bwareaopen(mask3,1000); % Remove all object containing fewer than 1000 pixels.
[B,L,objectsNum,A] = bwboundaries(bw); % A is the adjacency matrix.
% Select only the region and not the boundarys inside the region:
B = B(1:objectsNum);
L(L>objectsNum) = 0;
% Display the label matrix and draw each boundary:
%figure(4);
%imshow(label2rgb(L,@jet,[.5 .5 .5])); 
%F(2) = getframe;
%hold on
%title('Label');
for k = 1:length(B)
  boundary = B{k};
  if show_images 
      plot(boundary(:,2),boundary(:,1),'w','LineWidth',1)
  end
end

stats = regionprops(L,'all');

% Write the coreesponding number inside each object

for k = 1:numel(stats)
    center= stats(k).Centroid;
    if show_images 
        text(center(1), center(2), sprintf('%d', k), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end
end

%% ===================Shape =========================
% 1.Roundness - Elliptic-normalized circumference 
% 2.ENC
% 3.LW
% 4.LN - Nnumber of lobes corresponds to the concave points of the lesion
% 5.LI - The distribution size of the lobes in a lesion
%% Roundness
roundness_array = zeros(1,(length(B)));
perimeter_array = zeros(1,(length(B)));
area_array = zeros(1,(length(B)));
for k = 1:length(B)
  boundary = B{k};
  delta_sq = diff(boundary).^2;    
  perimeter = sum(sqrt(sum(delta_sq,2)));
  perimeter_array(k) = perimeter;
  area_array(k) = stats(k).Area; % Number of pixels in the region
  Dmax = stats(k).MajorAxisLength;
  metric = (4*area_array(k))/(pi*(Dmax^2));
  metric_string = sprintf('%0.2f',metric);
  roundness_array(k) = metric;
end

%% ENC
ENC_array = zeros(1,length(B));
for k = 1:length(B)
    major = stats(k).MajorAxisLength;
    minor = stats(k).MinorAxisLength;
    primeter_ellipse = 2*pi*(sqrt(major^2 + minor^2)/2);
    ENC_array(k) = perimeter_array(k)/primeter_ellipse;
end   


%% LW
LW_array = zeros(1,length(B));
for k = 1:length(B)
    LW_array(k) = stats(k).MajorAxisLength/stats(k).MinorAxisLength;
end


%% LN
% Multiply the convex hull of the region by the not of the region.
if length(B)>0
    CH1 = stats(1).ConvexImage;
    im = stats(1).Image;
    im = ~im;
    side_blobs = im.*CH1; 
    side_blobs = bwareaopen(side_blobs,350); % Remove all regoins that are smaller than 100 pixel.
    if show_images
        figure(5);
        imshow(side_blobs);
        F(3) = getframe;
    end
    % Count how many regions surround the blobs 
    [B_LI,L_LI,num_blobs,A_LI] = bwboundaries(side_blobs);
end

%% LI
%% =============== texture ===============
% 1.Lesion Boundary
% 2.Average
% 3.Variance
% 4.Contrast
% 5.PAF - Posterior Acoustic Feature
%% Lesion Boundary
% Calculate the area around the region
SE1 = strel('square',11);
mask4_1 = imdilate(bw2,SE1);
if show_images  
    figure(6)
    imshow(mask4_1); 
    F(4) = getframe;
end 
mask4_2 = ~mask3;
mask4_3 = mask4_2.*mask4_1;
if show_images
    figure(7)
    imshow(mask4_3);
    F(5) = getframe;
end 
[B_around,L_around,n_around,A_around] = bwboundaries(mask4_3); %n is the number of objects %A is the  adjacency matrix
B_around = B_around(1:n_around);
L_around(L_around>n_around) = 0;

sum_values_array_around = zeros(1,(length(B_around)));
avg_tissue_array_around = zeros(1,(length(B_around)));
area_array_around = zeros(1,(length(B_around)));
stats_around = regionprops(L_around,'all');
for k = 1:length(B_around)
  area_array_around(k) = stats_around(k).Area;
end

% Calculating the average intensity around the tissue 
for k = 1:length(B_around)
    for row=1:size(L_around,1)
        for col=1:size(L_around,2)
         if (L_around(row,col)==k)
             sum_values_array_around(k)=  sum_values_array_around(k) + double(I_crop(row,col));
         end
        end
    end
    avg_array_around(k) = sum_values_array_around(k)/area_array_around(k);
end

% Calculating the average intensity of a strip inside the tissue
% Making a mask to show the strip inside the tissue
mask4_1 = imerode(mask3,SE1);
mask4_1 = ~mask4_1;
mask4_2 = mask3.*mask4_1;
if show_images
    figure(8)
    imshow(mask4_2);
    F(6) = getframe;
end
[B_strip_inside,L_strip_inside,n_strip_inside,A_strip_inside] = bwboundaries(mask4_2); %n is the number of objects %A is the  adjacency matrix
B_strip_inside = B_strip_inside(1:n_strip_inside);
L_strip_inside(L_strip_inside>n_strip_inside) = 0;
sum_values_array_strip_inside = zeros(1,(length(B_strip_inside)));
avg_tissue_array_strip_inside = zeros(1,(length(B_strip_inside)));
area_array_strip_inside = zeros(1,(length(B_strip_inside)));
stats_strip_inside = regionprops(L_strip_inside,'all');
for k = 1:length(B_strip_inside)
  area_array_strip_inside(k) = stats_strip_inside(k).Area;
end

% Calculating the average intensity of a strip inside the tissue 
for k = 1:length(B_strip_inside)
    for row=1:size(L_strip_inside,1)
        for col=1:size(L_strip_inside,2)
         if (L_strip_inside(row,col)==k)
             sum_values_array_strip_inside(k)=  sum_values_array_strip_inside(k) + double(I_crop(row,col));
         end
        end
    end
    avg_array_strip_inside(k) = sum_values_array_strip_inside(k)/area_array_strip_inside(k);
end

%% Important:
% Creating image that containing only the biggest object
biggest_object=bwpropfilt(bw,'Area',1); % If we want to analyze the second big we change to 2
biggest_obj_image = I_crop;
biggest_obj_image(biggest_object == 0) = biggest_object(biggest_object == 0);
if show_images
    figure(9);
    imshow(biggest_obj_image);
    F(7) = getframe;
    impixelinfo;
    title('Image conatins the biggest object');
end
%% Echo pattern
% Creating a gray-level co-occurrence matrix (GLCM) from the image that
% contains only the main (biggest) object
[glcms_biggest_obj,SI_biggest] = graycomatrix(biggest_obj_image(:,:,1));
stats_glcms_biggest_obj = graycoprops(glcms_biggest_obj);
if show_images
    figure(10);
    imshow(rescale(SI_biggest));
    title('rescale SI');
    F(8) = getframe;
end
%% Average 
glcms_biggest_obj(1,1)=0;
sum_rows=sum(glcms_biggest_obj');
sum_cols=sum(sum_rows);
average= sum_cols;
%% Variance 
variance=0;
for row=1:1:8
   for col=1:1:8
       variance=variance+((row-average)^2)*glcms_biggest_obj(row,col);
   end
end
%% Contrast 
contrast=stats_glcms_biggest_obj.Contrast;
%% PAF
% Average intensity array
sum_values_array = zeros(1,(length(B)));
for k = 1:length(B)
    for row=1:size(L,1)
        for col=1:size(L,2)
         if (L(row,col)==k)
             sum_values_array(k)=  sum_values_array(k) + double(biggest_obj_image(row,col));
         end
        end
    end
    avg_array(k) = sum_values_array(k)/area_array(k); %area_array contains number of pixels in the region
end
% First calculating avg_PA
avg_PA = zeros(1,(length(B)));
sum_pa=0;
N_pa=0;
for k = 1:length(B)
    major= stats(k).MajorAxisLength;
    minor= stats(k).MinorAxisLength;
    if minor<100 
        height=minor;
    else
        height=100; %i think 100 is too small
    end
    x_coord=center(1);
    y_coord=round(center(2));
    possible_left_point=round(x_coord-(1/3)*major);
    if  possible_left_point<1 ,left_point=1; else ,left_point=possible_left_point; end
    possible_right_point=round(x_coord+(1/3)*major);
    if possible_right_point>col_size ,right_point=col_size; else ,right_point=possible_right_point; end
    possible_height_PA=round(y_coord+height);
    if possible_height_PA>row_size ,height_PA=row_size; else ,height_PA=possible_height_PA; end
    for y=y_coord:height_PA
        for x=left_point:right_point
%             sum_pa=sum_pa+double(I_crop(x,y));
            N_pa=N_pa+1;
        end
    end
%     avg_PA(k)=sum_pa/N_pa; 
end
PAF_array = zeros(1,(length(B)));
for k = 1:length(B)
   PAF_array(k)=abs(avg_PA(k)-avg_array(k));
end
%% Calculate how many 'holes' are inside the reigon
% mask
%bw_repmat=double(repmat(bw,[1 1 1]));
% bw_repmat=repmat(bw,[1 1 1]);
bw_repmat=uint8(repmat(bw,[1 1 3]));
mask=I_crop.*bw_repmat; % not necessary one object
bw_repmat=uint8(repmat(bw,[1 1 3]));

if show_images 
    figure(11);
    imshow(mask,[]);
    F(9) = getframe;
    impixelinfo;
    title('mask');
end
mask_bright = rgb2gray(mask);
bw_bright = mask_bright > 140;
bw_bright = bwareaopen(bw_bright,100); %Remove all object containing fewer than 100 pixels.
bw_bright=~bw_bright;
bw_bright = bwareaopen(bw_bright,50); %Remove all object containing fewer than 50 pixels.
[B_bright,L_bright,n_bright,A_bright] = bwboundaries(bw_bright);
% Select only the region and not the boundarys inside the region
B_bright = B_bright(1:n_bright);
L_bright(L_bright>n_bright) = 0;
if show_images 
    figure(12);
    imshow(bw_bright);
    F(10) = getframe;
    title('bw_bright');
end
num_holes = length(B_bright);
%% Animated Image 
% axis tight manual % This ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
% for objectsNum = 1:1:10
%       % Capture the plot as an image 
%       im = frame2im(F(objectsNum)); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if objectsNum == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
%  end
%% Exporting data to excel
% if length(B)>0
%     dataMatrix=[ roundness_array(1); ENC_array(1); LW_array(1) ; num_blobs; 0; ... 
%         avg_array_strip_inside(1) ; average ; variance ; ...  % Not sure about avg_array_strip_inside
%         contrast ; PAF_array(1) ; num_holes ];
%     rowDataNames= [ "Roundness" ; "ENC" ; "LW" ; "LN" ; "LI" ;...
%         "Lesion Boundary" ; "Average" ; "Variance" ;...
%         "Contrast" ; "PAF" ; "Holes Number"];
%     colDataNames="object1";
%     xlswrite('Breast2_denoised_croped.xls',rowDataNames,1,'A2:A12');
%     xlswrite('Breast2_denoised_croped.xls',colDataNames,1,'B1:B1');
%     xlswrite('Breast2_denoised_croped.xls',dataMatrix,1,'B2:B12');
% end
%% Creating a sturct
Features.Roundness = roundness_array(1)
Features.ENC = ENC_array(1);
Features.LW = LW_array(1);
Features.NumBlobs = num_blobs;
Features.AvgAround = avg_array_around(1);
Features.AvgInside = avg_array_strip_inside(1);
Features.Average = average;
Features.Variance = variance;
Features.Contrast = contrast;
Features.PAF = PAF_array(1);
Features.NumHoles = num_holes;
end

