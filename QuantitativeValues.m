function [Roundness, ENC, LW, NumBlobs, Average, Variance, Contrast, PAF, NumHoles, EPI, SI, MSE, PSNR] = QuantitativeValues(OriginalIm,FI)
    imwrite(FI, 'FI.png');
    I = imread('FI.png');
    Features = FeaturesExtraction(I);
    Roundness = Features.Roundness;
    ENC = Features.ENC;
    LW = Features.LW;
    NumBlobs = Features.NumBlobs;
    Average = Features.Average;
    Variance = Features.Variance;
    Contrast = Features.Contrast;
    PAF = Features.PAF;
    NumHoles = Features.NumHoles;
    
    EPI = EdgePreservationIndex(OriginalIm,FI); % first original, second despeckled 
    SI = SpeckleIndex(FI);
    
    diffFI = abs(double(OriginalIm)-double(FI));
    MSE =sum(sum(diffFI))/numel(diffFI); % Mean squared error
    MSE=MSE(1);
    diffFI = rgb2gray(diffFI);
    diffFI = imadjust(diffFI);
    metrics = my_MetricsMeasurement(OriginalIm,FI);
    PSNR = metrics.PSNR;

end

