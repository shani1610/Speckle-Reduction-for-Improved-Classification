# Effect of despeckle filtering on classification of breast ultrasound images

Students: [Inbal Aharoni](mailto:AHARONINBAL@gmail.com),  [Shani Israelov](mailto:shani1610@gmail.com), Supervised by: [Shira Nemirovsky-Rotman](mailto:)
August 2020, VISL lab, Electrical Engineering department, Technion

## Indroduction
In this project, different speckle reduction techniques have been studied and implemented in software,
and then tested on breast images in order to examine the impact on lesions classification as malignant or benign.
Speckle pattern tends to obscure edges and reduce the image contrast and therefore is interrupting the radiologists in the medical diagnosis.
The first method is OBLMN Algorithm which we didn’t implemented ourselves
and the second is Adaptive Wavelet Thresholding which we implemented in MATLAB. 
We examine the impact on lesions classification in two ways. 
The first one is quantitative Features' Extraction. The features divided into morphological and textural features that characterize the mass.
The second way to examine the impact is by radiologist professional opinion. Finally, we made conclusions. 
the code was made to make a comparison between methods to speckle denoising in breast uldrasound images. 
the images are given in a folder. 

## how to run? 
to run the code you need to run the Test.m file.  make sure the that specified in the Test file is included in the folder. 

## whats included? 
1.	Test
2.	ImageDescription
3.	WaveletDecomposition - our implementation to decomposition, up to 3 levels. 
4. SpeckleRemovalOBNLM - Bayesian NLM filter made by:
  %     %  *  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
  %     %  *  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
  %     %  *  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, *
  %     %  *  Avril 2008        
  %  
5. 	QuantitativeValues:
  % *	FeatureExtraction - returns a sturct that includes morphogical and texturical features
  % *	EdgePreservationIndex
  % *	SpeckleIndex
  % *	My_MerticsMeasurment - MSE, PSNR calculation made by:     
  %     %  *  ASHISH MESHRAM (meetashish85@gmail.com; www.facebook.com/ashishmeet)       

