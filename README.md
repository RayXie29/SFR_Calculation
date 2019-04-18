**#ISO 12233:1999 SFR  calculation**
<br /> 
SFR is the short for "Spatial Frequency Response". It represents the response of an imaging system in different spatial frequency. Imagining we use a almost perfect imaging system to take a picture of an pattern which has strong pixel value chages(from 0 to 255 or 255 to 0). The pixel value change of edge will vary very fast and strong. But if we use a bad imaging system to take the same picture, the pixel value change of the edge would be slow and small. Because a bad imaging system will have serious diffraction while the light incident through lens of imaging system lens. This will scatter the light and let the light received by several pixels. So the resoultion of the edge part will decrease. And the SFR is using this value change of edge to estimate the quality of imaging system. <br />
<br />
<br />
Normally, the SFR calculation will use image of ISO12233 test chart which took by the imaging system we want to test. And the roi of slanted-edge pattern in ISO12233 test chart will be used for SFR calculation.
<br />
<br />
Here is the original image which took by my web-cam
![alt text](https://raw.githubusercontent.com/RayXie29/SFR_Calculation/master/imgs/original_img.bmp)
<br />
<br />
I using opencv mouse event to crop the ROI<br />
![alt text](https://raw.githubusercontent.com/RayXie29/SFR_Calculation/master/imgs/cropping.PNG)
<br />
<br />
Here is the ROI I cropped<br />
![alt text](https://raw.githubusercontent.com/RayXie29/SFR_Calculation/master/imgs/roi.PNG)
<br />
After cropping the ROI, SFR calculation will begin and it will auto save a csv file of mtf values in different spatial frequency.
<br />
<br />
Reference: ISO12233:1999(E)<br />
