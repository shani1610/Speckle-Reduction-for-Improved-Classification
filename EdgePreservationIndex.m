function EPI = EdgePreservationIndex(originalIm,despeckledIm)
    % Convert to gray and double 
    originalIm = rgb2gray(originalIm);
    originalIm = im2double(originalIm);
%    despeckledIm = rgb2gray(despeckledIm);
    despeckledIm = im2double(despeckledIm);
    
    % Calculate Edge Preservation Index
    [m, n] = size(originalIm);
    [m_d, n_d] = size(despeckledIm);
    if m~=m_d || n~=n_d 
        error('The images are not from the same size');
    end
    p_despeckledIm = sum(sum(abs(diff(despeckledIm,1,2)),2));
    p_originalIm = sum(sum(abs(diff(originalIm,1,2)),2));
    EPI =  p_despeckledIm/p_originalIm;
end

