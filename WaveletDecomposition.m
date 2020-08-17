function FI= WaveletDecomposition(Im,level,wave_name)
    Im = im2double(rgb2gray(Im));
    size_Im = size(Im);
%--------------Single level of decomposition------------------------------
    [cA1,cH1,cV1,cD1] = dwt2(Im,wave_name);
    [cV1,cH1, cD1] = activate_threshold(cA1,cV1,cH1,cD1);
    if (level == 1)
        FI = idwt2(cA1,cH1,cV1,cD1,wave_name); %Inverse DWT     
    end
%--------------Second level of decomposition------------------------------
    if ((level == 2)||(level == 3))
        [cA2,cH2,cV2,cD2] = dwt2(cA1,wave_name); %Discrete Wavelet Transform
        [cV2,cH2, cD2] = activate_threshold(cA2,cV2,cH2,cD2);

        if(level == 2)
            FI = idwt2(cA2,cH2,cV2,cD2,wave_name); %Inverse DWT
            FI = imresize(FI,size(cA1));
            FI = idwt2(FI,cH1,cV1,cD1,wave_name); %Inverse DWT
        end
    end
%--------------third level of decomposition------------------------------
    if (level == 3)
        [cA3,cH3,cV3,cD3] = dwt2(cA2,wave_name); %Discrete Wavelet Transform
        [cV3,cH3, cD3] = activate_threshold(cA3,cV3,cH3,cD3);
         FI = idwt2(cA3,cH3,cV3,cD3,wave_name); %Inverse DWT
         FI = imresize(FI,size(cA2));
         FI = idwt2(FI,cH2,cV2,cD2,wave_name); %Inverse DWT
         FI = imresize(FI,size(cA1));
         FI = idwt2(FI,cH1,cV1,cD1,wave_name); %Inverse DWT
    end
    FI = imresize(FI,size_Im); % Resizing filtered image to preprocessed image
end
  
%% 
function [cV1,cH1, cD1] = activate_threshold(cA1,cV1,cH1,cD1)
    
    v_cH1 = var(cH1, 0,[1 2]);
    v_cV1 = var(cV1, 0,[1 2]);
    v_cD1 = var(cD1, 0,[1 2]);   
    
    meancH1 = mean(cH1,'all'); 
    meancV1 = mean(cV1,'all'); 
    meancD1 = mean(cD1,'all'); 
    
    [r c] = size(cA1);
    
    tv1 = sqrt((meancH1+v_cH1)*2*log(r*c));%TODO - check the variance
    tv2 = sqrt((meancV1+v_cV1)*2*log(r*c));
    tv3 = sqrt((meancD1+v_cD1)*2*log(r*c));

    threshold_vector1 = [ tv1 tv2 tv3 ]; %cH, cV, vD 
    
    cV1 = wthresh(cV1,'s',threshold_vector1(2)); %Eliminating HL
    cH1 = wthresh(cH1,'s',threshold_vector1(1)); %Eliminating LH
    cD1 = wthresh(cD1,'s',threshold_vector1(3)); %Eliminating HH
end

