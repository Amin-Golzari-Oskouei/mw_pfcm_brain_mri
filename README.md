# mw_pfcm_brain_mri
An Adaptive Multi-Weighted Possibilistic Fuzzy C-Means Clustering Approach for Brain Tumor Segmentation in Magnetic Resonance Imaging

The possibilistic fuzzy c-means (PFCM) method is a well-known approach to practical engineering for image segmentation. However, it still has some limitations, including 1) considering equal importance to all image features, 2) sensitivity to initialization, and 3) the neighborhood information of the image pixels is not considered in the objective function. To alleviate the mentioned challenges of possibilistic fuzzy c-means, this paper proposes a novel adaptive multi-weighted possibilistic fuzzy c-means image segmentation approach. The proposed possibilistic fuzzy objective function is composed of 1) a feature weighting factor to distinguish between the importance of different features, 2) a cluster weighting schema to mitigate the initialization sensitivity, 3) a spatial constraint term containing local information of pixels to reduce the algorithm's sensitivity to noise. In addition, the curvelet transform is employed in this research to enhance the edges in brain Magnetic Resonance Imaging (MRI) images, their contrast, and brightness intensity. Moreover, an efficient combination of image features increases the segmentation accuracy. Finally, mathematical analyses are provided to obtain the required updating functions. In summary, the artificial intelligence-powered segmentation method integrates adaptive weighting and spatial constraints to better process and interpret complex brain images and, as a result, develop efficient computer-aided diagnosis (CAD) systems. The proposed method's performance is evaluated and compared with state-of-the-art methods. The results show that the proposed method outperforms the competitors.

# Used Dataset:
To evaluate the proposed method, the benchmark dataset has been used. There are some images of this dataset in the uploaded file (Dataset folder). Whole dataset is available on: https://www.kaggle.com/datasets/mateuszbuda/lgg-mri-segmentation.

# Overview of the mw_pfcm_brain_mri:
<img width="630" height="215" alt="Untitled" src="https://github.com/user-attachments/assets/23ed51b3-a276-43e1-8ee9-48011f6ea65b" />


# Comment:
The repository file includes the MATLAB implementation of the mw_pfcm_brain_mri algorithm.

Comments are written for all steps of the algorithm for better understanding the code. Also, a demo is implemented for ease of running, which runs by importing the data and other necessary algorithm parameters.
