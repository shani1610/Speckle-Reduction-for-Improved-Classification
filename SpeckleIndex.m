function SI = SpeckleIndex(Im)
    % Convert to gray and double
    %Im = rgb2gray(Im);
    Im = im2double(Im);
    
    % Calculate Speckle Index
    [m, n] = size(Im);
    w = 5; % window size For example 3x3 or 5x5
    sum = 0;
    for i=ceil(w/2):n-floor(w/2)
      for j=ceil(w/2):m-floor(w/2)
          neighbours=Im(j-floor(w/2):j+floor(w/2), i-floor(w/2):i+floor(w/2));
          if size(neighbours)~=[w,w]
              error('size of neighbours matrix is not as expected');
          end
          local_std=std2(neighbours);
          neighbours(w+ceil(w/2))=[];
          local_mean=mean(neighbours);
          local_val = local_std/local_mean;
          if isnan(local_val) || local_mean==0
              continue;
          end
          sum = sum + local_val;
      end
    end
    SI = (1/(m*n))*sum;
end

