The alignment of carbon nanotube (CNT) in films or fibers is a very significant parameter that demonstrates the quality and determines the application. In this paper, a simple, faster, precise, and economical approach based on image processing method is proposed to provide an end-to-end solution for determining the orientation of CNTs in films within a short amount of time with no expertise required. As preprocessing steps, Laplacian edge enhancement filter is applied to the image for improving the appearance of edge regions and thereafter the image is partitioned into multiple image blocks to capture the nanoscale orientation characteristics. The extraction of textural component from this edge enhanced image blocks further alleviates the noisy smoothing components and improves the fine scale details present in the nanostructures. Afterwards, 2D-fast Fourier transform (2D-FFT) is applied on the textural component of edge enhanced image blocks rather than the original image to determine the orientation distribution of CNT films which is later utilized to estimate 2D orientation parameter. Additionally, a Fourier subtraction technique is introduced to discard the holes, cracks, or contaminants in the films which may create unwanted alignment information. Extensive analysis on not only scanning electron microscopy but also optical images of CNT films with different order of alignment was performed to quantify the local and global alignment in those samples. To show the effectiveness of our method, the results were corroborated with state-of-the-art experimental techniques. 

To simulate the program run "main.m" ( can be found in Matlab code" folder.

The results found from each of the steps are depicted below:  

![image_for_paper_4_updated_new_reduced_1](https://user-images.githubusercontent.com/28448474/133827070-32faf541-d40d-49f5-9dd2-a0fa5ae8428a.jpg)

Different analyses of our method and the results found from absorption spectroscopy and THz-TDS measurements are presented below:

![new_image_1](https://user-images.githubusercontent.com/28448474/133827364-242c9af1-b031-456a-b545-1d1b7c456f3f.jpg)
![THz_Data_1](https://user-images.githubusercontent.com/28448474/133827373-da7ee95f-22b9-44e7-ba0e-2dd2404d3e3d.jpg)


