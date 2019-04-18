#include "ISOsfr.h"
#include <opencv2/opencv.hpp>
#include <complex>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>



void de_Gamma(cv::Mat &Src, double gamma)
{
	if (Src.channels() != 1) { return; }

	for (int i = 0; i < Src.rows; ++i)
	{
		uchar *SrcP = Src.ptr(i);
		for (int j = 0; j < Src.cols; ++j)
		{
			SrcP[j] = 255 * (pow((double)SrcP[j] / 255, 1 / gamma));
		}
	}
}

void SLR(std::vector<double> &Cen_Shifts, std::vector<double> &y_shifts, double *a, double *b)
{
	//a -> intercept, b->slope of slanted edge

	int y_size = y_shifts.size();
	//using simple linear regression to solve equation and get slope and intercept
	double xsquare = 0, xavg = 0, yavg = 0;
	int i;
	for (i = 0; i < y_size; ++i)
	{
		yavg += y_shifts[i];
		xavg += Cen_Shifts[i];
	}
	xavg /= (double)y_size;
	yavg /= (double)y_size;

	//simple linear regession
	for (i = 0; i < y_size; ++i)
	{
		double temp = (Cen_Shifts[i] - xavg);
		*b += temp * (y_shifts[i] - yavg);
		xsquare += temp * temp;
	}

	*b /= xsquare;
	*a = yavg - (*b)*xavg;
}

std::vector<double> CentroidFind(cv::Mat &Src,std::vector<double> &y_shifts, double *CCoffset)
{
	std::vector<double> Cen_Shifts(Src.rows);
	int i, j, height = Src.rows, width = Src.cols;

	cv::Mat tempSrc(Src.size(), CV_8UC1);

	//Do the bilaterFilter on template ROI image to make sure we can find the slanted edge more acurrate
	cv::bilateralFilter(Src, tempSrc, 3, 200, 200);

	// calculate the centroid of every row in ROI
	// centroid formuld: Total moments/Total amount
	for (i = 0; i < height; ++i)
	{
		double molecule = 0, denominator = 0, temp = 0;
		uchar *tempSrcPtr = tempSrc.ptr<uchar>(i);
		for (j = 0; j < width - 1; ++j)
		{
			temp = (double)tempSrcPtr[j + 1] - (double)tempSrcPtr[j];
			molecule += temp * j;
			denominator += temp;
		}
		Cen_Shifts[i] = molecule / denominator;
	}

	//Eliminate the noise far from slant edge(+/- 10 pixels)
	tempSrc = Src.clone();
	cv::GaussianBlur(Src, Src, cv::Size(3, 3), 0);
	for (i = 0; i < height; ++i)
	{
		uchar *SrcPtr = Src.ptr<uchar>(i);
		uchar *tempPtr = tempSrc.ptr<uchar>(i);
		for (j = int(Cen_Shifts[i]) - 10; j < int(Cen_Shifts[i]) + 10; ++j)
		{
			SrcPtr[j] = tempPtr[j];
		}
	}

	// check whether the edge in image is too close to the image corners
	if (Cen_Shifts[0] < 2.0 || width - Cen_Shifts[0] < 2.0)
	{
		std::cerr << "The edge in ROI is too close to the image corners" << std::endl;
		return {};
	}

	if (Cen_Shifts[height - 1] < 2 || width - Cen_Shifts[height - 1] < 2)
	{
		std::cerr << "The edge in ROI is too close to the image corners" << std::endl;
		return {};
	}

	int half_ySize = height / 2;
	int CC = Cen_Shifts[half_ySize]; //Centroid of image center
	*CCoffset = CC;
	for (i = 0; i < height; ++i)
	{
		Cen_Shifts[i] -= CC; //Calculate the Shifts between Centroids and Centroid of image center		
		y_shifts[i] = i - half_ySize; //Calculate the shifts between height of image center and each row
	}

	return Cen_Shifts;
}


std::vector<double> OverSampling(cv::Mat &Src, double slope, double CCoffset, int height, int width, int *SamplingLen)
{
	std::vector<double> RowShifts(height, 0);

	int i, j, k;
	int halfY = height >> 1;

	//calculate the pixel shift of each row to align the centroid of each row as close as possible
	for (i = 0; i < height; ++i) { RowShifts[i] = (double)(i - halfY) / slope + CCoffset;  }

	//DataMap -> to record the index of pixel after shift
	//Datas -> to record the original pixel value 
	std::vector<double> DataMap(height*width, 0);
	std::vector<double> Datas(height*width, 0);

	for (i = 0, k = 0; i < height; ++i)
	{
		int baseindex = width * i;
		for (j = 0; j < width; ++j)
		{
			DataMap[baseindex+j] = j - RowShifts[i];
			Datas[baseindex + j] = int(Src.at<uchar>(i, j));
		}
	}

	std::vector<double> SamplingBar(*SamplingLen, 0);
	std::vector<int> MappingCount(*SamplingLen, 0);

	//Start to mapping the original data to 4x sampling data and record the count of each pixel in 4x sampling data
	for (i = 0; i < height*width; ++i)
	{
		int MappingIndex = static_cast<int>(4 * DataMap[i]);
		if (MappingIndex >= 0 && MappingIndex < *SamplingLen)
		{
			SamplingBar[MappingIndex] = SamplingBar[MappingIndex] + Datas[i];
			MappingCount[MappingIndex]++;
		}
	}

	//average the pixel value in 4x sampling data, if the pixel value in pixel is zero, copy the value of close pixel
	for (i = 0; i < *SamplingLen; ++i) {
		j = 0;
		k = 1;
		if (MappingCount[i] == 0) {

			if (i == 0) {
				while (!j)
				{
					if (MappingCount[i + k] != 0)
					{
						SamplingBar[i] = SamplingBar[i + k] / ((double)MappingCount[i + k]);
						j = 1;
					}
					else ++k;
				}
			}
			else {
				while (!j && ((i - k) >= 0))
				{
					if (MappingCount[i - k] != 0)
					{
						SamplingBar[i] = SamplingBar[i - k];   /* Don't divide by counts since it already happened in previous iteration */
						j = 1;
					}
					else ++k;
				}
				if ((i - k) < 0)
				{
					k = 1;
					while (!j)
					{
						if (MappingCount[i + k] != 0)
						{
							SamplingBar[i] = SamplingBar[i + k] / ((double)MappingCount[i + k]);
							j = 1;
						}
						else ++k;
					}
				}
			}
		}
		else
			SamplingBar[i] = (SamplingBar[i]) / ((double)MappingCount[i]);
	}

	// reduce the length of sampling data 
	// because the datas at the edge are only matters, we truncate the data close to the two side of sampling data
	// which has very small contribution to the result

	//truncating the data smaller than 10% and bigger than 90% of original length
	int originalSamplingLen = *SamplingLen;
	*SamplingLen = *SamplingLen * 0.8;

	std::vector<double> deSampling(*SamplingLen,0);
	//derivative sampling data(which is ESF) to get the line spread function
	for (i = originalSamplingLen * 0.1, j = 1; i < originalSamplingLen, j < *SamplingLen; ++i,++j)
	{
		deSampling[j] = SamplingBar[i + 1] - SamplingBar[i];
	}
	
	return deSampling;
}


std::vector<double> HammingWindows(std::vector<double> &deSampling, int SamplingLen)
{
	int i, j;
	std::vector<double> tempData(SamplingLen);

	//We want to shift the peak data to the center of line spread function data
	//Because we will do the hamming window later, this will keep the important data away from filtering
	//In case there are two peaks, we use two variable to record the peak data position
	int L_location = -1, R_location = -1;

	double SamplingMax = 0;
	for (i = 0; i < SamplingLen; ++i)
	{
		if (fabs(deSampling[i]) > fabs(SamplingMax))
			SamplingMax = deSampling[i];
	}

	for (i = 0; i < SamplingLen; ++i)
	{
		if (deSampling[i] == SamplingMax)
		{
			if (L_location < 0) { L_location = i; }
			R_location = i;
		}
	}

	//the shift amount 
	int PeakOffset = (R_location + L_location) / 2 - SamplingLen / 2;

	if (PeakOffset)
	{
		for (i = 0; i < SamplingLen; ++i)
		{
			int newIndex = i - PeakOffset;
			if (newIndex >= 0 && newIndex < SamplingLen) { tempData[newIndex] = deSampling[i]; }
		}
	}
	else
	{
		for (i = 0; i < SamplingLen; ++i) { tempData[i] = deSampling[i]; }
	}

	//do the hamming window filtering
	for (int i = 0; i < SamplingLen; ++i)
	{
		tempData[i] = tempData[i] * (0.54 - 0.46*cos(2 * M_PI*i / (SamplingLen - 1)));
	}

	return tempData;
}

void DFT(std::vector<double> &data, int size)
{
	int i, j;
	std::complex<double> *arr = new std::complex<double>[size];
	for (i = 0; i < size; ++i)
	{
		arr[i] = std::complex<double>(data[i], 0);
	}

	for (i = 0; i < size / 2.0; ++i)
	{
		std::complex<double> temp = 0;
		for (j = 0; j < size; ++j)
		{
			double w = 2 * 3.1415*i*j / size;
			std::complex<double> deg(cos(w), -sin(w));
			temp += arr[j] * deg;
		}
		data[i] = sqrt(temp.real()*temp.real() + temp.imag()*temp.imag());
	}
}

void ReduceRows(double slope, int *ImgHeight)
{
	double tempSlope = fabs(slope);
	int cycs = (*ImgHeight) / tempSlope;
	if (tempSlope*cycs < *ImgHeight) { *ImgHeight = tempSlope * cycs; }
}

double SFRCalculation(cv::Mat &ROI, double gamma)
{
	if (ROI.empty())
	{
		std::cerr << "Open the ROI image error" << std::endl;
		return 0;
	}

	int height = ROI.rows, width = ROI.cols;
	//Do the gamma decoding to eliminate the gamma encoded by camera device 
	de_Gamma(ROI, gamma);
	int i, j;

	double slope = 0, intercept = 0;
	
	//Center centroid offset
	double CCoffset = 0;
	std::vector<double> y_shifts(height);
	
	////Calculate the shifts between Centroids and Centroid of image center	
	std::vector<double> Cen_Shifts = CentroidFind(ROI, y_shifts, &CCoffset);
	if (Cen_Shifts.empty()) { return 0; }

	//simple linear regression for slanted edge fitting
	SLR(Cen_Shifts, y_shifts, &intercept, &slope);

	//Truncate the number of rows of data to largest slope cycle which will have an integer number of full phase rotations
	ReduceRows(slope, &height);

	//update the CCoffset to the offset between original mid point of image and reference mid point we calculated
	CCoffset = CCoffset + 0.5 + intercept - width / 2;

	//Mapping the pixel value of original image into a sampling data which the length is 4 times of original image width
	//This step is for concentrating the amount of change of original pixel values
	int SamplingLen = width * 4;
	std::vector<double> OverSamplingData = OverSampling(ROI, slope, CCoffset, height, width, &SamplingLen);

	//Using hamming window to filter the ripple signal of two side of data
	OverSamplingData = HammingWindows(OverSamplingData, SamplingLen);

	//decrete four transform
	DFT(OverSamplingData, SamplingLen);
	width = int(SamplingLen / 4);
	double maxData = 0;
	for (i = 0; i < SamplingLen; ++i)
	{
		if (OverSamplingData[i] != 0)
		{
			maxData = OverSamplingData[i];
			break;
		}
	}

	for (int i = 0; i < SamplingLen; ++i) { OverSamplingData[i] /= maxData; }

	std::fstream mtf_file("./ref/mtf.csv", std::ios::out);
	
	for (i = 0; i <= width; ++i)
	{
		double frequency = (double)i / width;
		mtf_file << frequency << "," << OverSamplingData[i] << "\n";
	}

	mtf_file.close();
}