function [aligned] = align_stacks_simple(data)
%Author: Nuria Vendrell Llopis, Carmen Diaz Verdugo and Emre Yaksi
%Reference: Reiten et al, Current Biology, 2016, DOI 10.1016/j.cub.2016.11.036   
%Modified by Christa Ringers and Nathalie Jurisch-Yaksi
%Reference: Ringers et al, BioRXiv, 2021, https://doi.org/10.1101/2021.11.23.469646 

% im_vector_orig=im_vector_orig %% Only to play around EY 01.07.12
warning off all


first_lsm = single(data(:,:,1));   %the first raw image will be the control one
LSMs = horzcat(data);
LSMs= single(LSMs(:,:,1:size(data,3)));
control_map = first_lsm;

%init_variable
im_raw_vector_corr(size(LSMs,3)).im = first_lsm; %to initialize the struct
im_raw_vector_corr(1,1).im = first_lsm;
min_y = 1;
max_y = size(data,1);
min_x = 1;
max_x = size(data,2);
im_aux = first_lsm;

try
    progress()
    for i = 2:size(LSMs,3)
        %disp (['Images left ', num2str(size(LSMs,3) - i)])
        progress(i,size(LSMs,3),10)
        stack = LSMs(:,:,i);
        [~, imOut] = alignImagesCPU(control_map,stack,16,2,'translation',1); %This parameters seems to work; don't give the mask
        x = find(isnan(imOut(ceil(max_y/2),:)));
        y = find(isnan(imOut(:,ceil(max_x/2))));
        min_x = max([x(x<ceil(max_x/2))+1, min_x]);
        max_x = min([x(x>floor(max_x/2))-1, max_x]);
        min_y = max([y(y<ceil(max_y/2))+1; min_y]);
        max_y = min([y(y>floor(max_y/2))-1; max_y]);
        imOut(isnan(imOut)) = 0;
        imOut(isinf(imOut)) = 0;
        im_aux(:,:,i) = imOut;
    end
catch
    sprintf('Aligned is too large maybe')
    
    
end

    
%After obtaining the min and max of x and y we can reshape the matrix
for i = 1:size(LSMs,3)
    im_raw_vector_corr(1,i).im = im_aux(min_y:max_y,min_x:max_x,i);
end


% Put the matrix in the required format
out=im_raw_vector_corr;
aligned=zeros (size(out(1).im,1),size(out(1).im,2),length(out)); % allocation for speed

for i=1:length(out)
   aligned(:,:,i)=horzcat(out(i).im);
end  