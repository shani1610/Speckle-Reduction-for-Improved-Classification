%description- for the image input, 
%Returns parameters that describe the distance from the original image and
%the noise reduction.
%prints 4 different transform of the image. 
function [] = ImageDescription(OriginalIm)

    figure(20);
    FI1 = WaveletDecomposition(OriginalIm,2,'sym4');
    [Roundness1, ENC1, LW1, NumBlobs1, Average1, Variance1, Contrast1, PAF1, NumHoles1, EPI1, SI1, MSE1, PSNR1] = QuantitativeValues(OriginalIm, FI1);
    subplot(2,2,1);
    imshow(FI1);
    title(['sym4 - level 2']);
    
    FI2= WaveletDecomposition(OriginalIm,2,'dmey');
    [Roundness2, ENC2, LW2, NumBlobs2, Average2, Variance2, Contrast2, PAF2, NumHoles2, EPI2, SI2, MSE2, PSNR2] = QuantitativeValues(OriginalIm, FI2);
    subplot(2,2,2);
    imshow(FI2);
    title(['dmey - level 2']);
    
    FI3 = WaveletDecomposition(OriginalIm,2,'db4');
    [Roundness3, ENC3, LW3, NumBlobs3, Average3, Variance3, Contrast3, PAF3, NumHoles3, EPI3, SI3, MSE3, PSNR3] = QuantitativeValues(OriginalIm, FI3);
    subplot(2,2,3);
    imshow(FI3);    
    title(['db4 - level 2']);

    FI4 = WaveletDecomposition(OriginalIm,3,'db4');
    [Roundness4, ENC4, LW4, NumBlobs4, Average4, Variance4, Contrast4, PAF4, NumHoles4, EPI4, SI4, MSE4, PSNR4] = QuantitativeValues(OriginalIm, FI3);
    subplot(2,2,4);
    imshow(FI4);
    title(['db4 - level 3']);
    
    %OBNLM
    fimgd = SpeckleRemovalOBNLM(OriginalIm);
    [Roundness5, ENC5, LW5, NumBlobs5, Average5, Variance5, Contrast5, PAF5, NumHoles5, EPI5, SI5, MSE5, PSNR5] = QuantitativeValues(OriginalIm, fimgd);

    %Original
    OriginalIm_rgb2gray = rgb2gray(OriginalIm);
    OriginalIm_rgb2gray = im2double(OriginalIm_rgb2gray);
    [Roundness6, ENC6, LW6, NumBlobs6, Average6, Variance6, Contrast6, PAF6, NumHoles6, EPI6, SI6, MSE6, PSNR6] = QuantitativeValues(OriginalIm, OriginalIm_rgb2gray);
    
    figure(22)
    subplot(3,1,1);
    imshow(FI1);
    title(['Wavelet (sym4,level 2)']);
    subplot(3,1,2);
    imshow(fimgd);
    title(['OBNLM']);
    subplot(3,1,3);
    imshow(OriginalIm_rgb2gray);
    title(['Original']);

    %% Exporting data to excel
    dataMatrix=[ Roundness1,Roundness2, Roundness3, Roundness4, Roundness5, Roundness6; ...
                  ENC1, ENC2, ENC3, ENC4, ENC5, ENC6;...
                  LW1, LW2, LW3, LW4, LW5, LW6;...
                  NumBlobs1, NumBlobs2, NumBlobs3, NumBlobs4, NumBlobs5, NumBlobs6;...
                  Average1, Average2, Average3, Average4, Average5, Average6;...
                  Variance1, Variance2, Variance3, Variance4, Variance5, Variance6;...
                  Contrast1, Contrast2, Contrast3, Contrast4, Contrast5, Contrast6;...
                  PAF1, PAF2, PAF3, PAF4, PAF5, PAF6;...
                  NumHoles1, NumHoles2, NumHoles3, NumHoles4, NumHoles5, NumHoles6;...
                  SI1, SI2, SI3, SI4, SI5, SI6;...
                  EPI1, EPI2, EPI3, EPI4, EPI5, EPI6 ;...
                  MSE1, MSE2, MSE3, MSE4, MSE5, MSE6 ;...
                  PSNR1, PSNR2, PSNR3, PSNR4, PSNR5, PSNR6 ];
    rowDataNames = [ "Roundness" ; "ENC" ; "LW" ; "Number of Blobs"; "Average" ; "Variance" ; "Contrast" ; "PAF" ; "Holes Number"; ...
                     "Speckle Index"; "Edge Preservation Index"; "MSE"; "PSNR" ];
    colDataNames= [ "sym4" , "dmey" , "db4" , "db4 Level 3", "OBNLM" , 'original' ];
    xlswrite('method_comparison.xls',rowDataNames,1,'A2:A14');
    xlswrite('method_comparison.xls',colDataNames,1,'B1:G1');
    xlswrite('method_comparison.xls',dataMatrix,1,'B2:G14');
end


    
