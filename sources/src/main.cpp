//Copyright (C) 2016, Eloïse Grossiord <eloise.grossiord@gmail.com>
//This program is free software: you can use, modify and/or
//redistribute it under the terms of the GNU General Public
//License as published by the Free Software Foundation, either
//version 3 of the License, or (at your option) any later
//version. You should have received a copy of this license along
//this program. If not, see <http://www.gnu.org/licenses/>.

#include <ctime>
#include <nifti1_io.h>
#include "cgraph.h"
#include "cgraphwatcher.h"
#include "colorordering.h"
#include "shaping.h"
#include "Types.h"
#include <typeinfo>
#include "ConnectedComponents.h"


using namespace std;

// for LibTIM classes
using namespace LibTIM;

class graphWatcher : public CGraphWatcher{
public :
    void progressUpdate(int ncur, int nfinal) {
    }

    void progressUpdate() {
        curProgress++;
    }
    graphWatcher(int finalProgress) : CGraphWatcher(finalProgress) {}
};



template <class T1, class T2>
Image <T2> adjustContrast(Image <T1> &im, T2 A, T2  B)
{
    Image<T2> imRes(im.getSize());

    T1 b=im.getMax();
    T1 a=im.getMin();

    if( (b-a)!=0)
    {
        double ratio=(double)(B-A)/(b-a);
        for(int i=0; i<im.getBufSize(); i++) {
            imRes(i)=(T2)(A+(im(i)-a)*ratio);
        }
    }
    return imRes;
}

template <typename T>
Image<T> readRestrainedNifti(char *filename, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, nifti_1_header* nfHeader=NULL)
{

    nifti_image *nim=NULL;
    nim=nifti_image_read(filename,1);
    if (nfHeader != NULL)
        *nfHeader = nifti_convert_nim2nhdr(nim);
    if(nim!=0) 
    {
        Image<T> inputImage(nim->nx,nim->ny,nim->nz);

	inputImage.setSpacingX(nim->dx);
	inputImage.setSpacingY(nim->dy);
	inputImage.setSpacingZ(nim->dz);
	
	int iBeg = xmin + (ymin * inputImage.getSizeX() ) + (zmin*inputImage.getSizeX()*inputImage.getSizeY());
	int iEnd = xmax + (ymax * inputImage.getSizeX() ) + (zmax*inputImage.getSizeX()*inputImage.getSizeY());


	for(int x = 0; x<inputImage.getSizeX(); x++)
	{
	    for(int y = 0; y<inputImage.getSizeY(); y++)
	    {
		for(int z = 0; z<inputImage.getSizeZ(); z++)
		{
		    int ind = x + (y*inputImage.getSizeX() ) + (z*inputImage.getSizeX()*inputImage.getSizeY());
		    if((ind >= iBeg) && (ind <= iEnd))
		    	inputImage(ind) = ((T *)nim->data)[ind];
		    else
		    	inputImage(ind) = 0;
		}
	    }
	}
        free(nim);
        return inputImage;
    }
    else {
        std::cout << "read error \n";
        exit(1);
    }
}

template<typename T>
void getSpacingAndOrigin(char *filename,float *quatern_b, float *quatern_c, float *quatern_d, float *qoffset_x, float *qoffset_y, float *qoffset_z, float *qfac, nifti_1_header* nfHeader=NULL)
{

    nifti_image *nim=NULL;
    nim=nifti_image_read(filename,1);
    if (nfHeader != NULL)
 	*nfHeader = nifti_convert_nim2nhdr(nim);
    if(nim!=0) 
    {   	
	*quatern_b = nim->quatern_b;
	*quatern_c = nim->quatern_c;    
	*quatern_d = nim->quatern_d;
	*qoffset_x = nim->qoffset_x;
	*qoffset_y = nim->qoffset_y;
	*qoffset_z = nim->qoffset_z;

	free(nim);
    	return;
    }
    else
    {
	std::cout << "read error \n";
        exit(1);
    }
}

template <typename T>
Image<T> readNifti(char *filename,nifti_1_header* nfHeader=NULL)
{
    nifti_image *nim=NULL;
    nim=nifti_image_read(filename,1);
    if (nfHeader != NULL)
        *nfHeader = nifti_convert_nim2nhdr(nim);
    if(nim!=0) {
        Image<T> inputImage(nim->nx,nim->ny,nim->nz);

	inputImage.setSpacingX(nim->dx);
	inputImage.setSpacingY(nim->dy);
	inputImage.setSpacingZ(nim->dz);
	 
 	for(int i=0; i<inputImage.getBufSize(); i++) {
             inputImage(i)=((T *)nim->data)[i];
         }
        free(nim);
        return inputImage;
    }
    else {
        std::cout << "read error \n";
        exit(1);
    }

}
 
template<typename PixelDataType>
void writeNifti(Image<PixelDataType> &im, const char *filename,float *quatern_b, float *quatern_c, float *quatern_d, float *qoffset_x, float *qoffset_y, float *qoffset_z, float *qfac, nifti_1_header* nfHeader=NULL)
{
    nifti_image *nim = NULL;
    if(nfHeader != NULL)
    {
	nim = nifti_convert_nhdr2nim(*nfHeader, filename); //cree image nifti avec le header de l'image d'entree
	int test = nifti_set_filenames(nim,filename,0,nim->byteorder);
    }
    else
    {
	int t[8] = {3,im.getSizeX(), im.getSizeY(), im.getSizeZ(), 0,0,0,0};
	if(typeid(PixelDataType) == typeid(U16))
         {       
                 nim = nifti_make_new_nim(t, 512, 1);
         }
         
         else if(typeid(PixelDataType) == typeid(U8))
         {
                 nim = nifti_make_new_nim(t, 2, 1);
         }
         else if(typeid(PixelDataType) == typeid(S8))
         {
                 nim = nifti_make_new_nim(t, 256, 1);
         }
         else if(typeid(PixelDataType) == typeid(S32))
         {
                 nim = nifti_make_new_nim(t, 8, 1);
         }
         else if(typeid(PixelDataType) == typeid(U32))
         {
                 nim = nifti_make_new_nim(t, 16, 1);
         }
         else
         {
                 assert(false);
         }
    }

//  Set spacing
    nim->dx = im.getSpacingX();
    nim->dy = im.getSpacingY();
    nim->dz = im.getSpacingZ();
    // Fill
    (nim->data) = (PixelDataType*)im.getData();
    int error = nifti_set_filenames(nim, filename, 0, nim->byteorder);
    if(0 == error)
    {
	
    	// set spacing and origin
    	nim->qform_code = 1;
        nim->sform_code = 1;
 
	nim->quatern_b = *quatern_b;
	nim->quatern_c = *quatern_c;
	nim->quatern_d = *quatern_d;
	nim->qoffset_x = *qoffset_x;
	nim->qoffset_y = *qoffset_y;
 	nim->qoffset_z = *qoffset_z;
	nim->qfac = *qfac;

	nifti_image_write(nim);
    }
    return;
}

// Keep band=0,1,2 and return associate unsigned char image
Image<unsigned char> keepBand(Image<RGB> &im, int band)
{
    assert(band==0 || band==1 || band==2);
    Image<unsigned char> imRes(im.getSize());

    for(int i=0; i<imRes.getBufSize(); i++) {
        imRes(i)=im(i)[band];
    }

    return imRes;
}

Image<unsigned char> reconstructFullImage(Image<unsigned char> & pet,Image<unsigned char> &mask, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    Image<unsigned char> imRes(pet.getSize());
    imRes.setSpacingX(pet.getSpacingX());
    imRes.setSpacingY(pet.getSpacingY());
    imRes.setSpacingZ(pet.getSpacingZ());
    imRes.fill(0);

    int indBig = 0;
    int indCrop = 0;
    int cropSizeX  = xmax-xmin+1;
    int cropSizeY = ymax-ymin+1;
    int cropSizeZ = zmax-zmin+1;
    // reconstruct Full Image with PET dimensions
    for(int x = xmin; x<xmax+1; x++)
    {
	for(int y = ymin; y<ymax; y++)
	{
	    for(int z = zmin; z<zmax+1; z++)
	    {
		indBig = x + (y*imRes.getSizeX() ) + (z*imRes.getSizeX()*imRes.getSizeY());
		indCrop = (x-xmin) + ((y-ymin)*cropSizeX) + ((z-zmin)*cropSizeX*cropSizeY);
	   	imRes(indBig) = mask(indCrop);
	    }
	}
    }

    return imRes;
}

Image<unsigned char> getLabeledMask(Image<unsigned char> &im, FlatSE connexity)
{
      std::cout<<"get Labeled Mask" << std::endl;
      Image<unsigned char> imRes(im.getSize());
      imRes.fill(0);
      // Construction binary mask
      unsigned char value = 0;
      for(int i=0; i<im.getBufSize(); i++)
      {
         if(im(i) > value)
         {
             //imRes(i)=(unsigned char)(im(i));
            imRes(i) = (unsigned char)255;
         }
      }
 
      // Label into CC
      Image<TLabel> imLabel(im.getSize());
      imLabel.fill(0);
      imLabel = labelConnectedComponents(imRes, connexity);
 
      // Convert into unsigned char
      Image<unsigned char> imLabeledMask(im.getSize());
      imLabeledMask.fill(0);
      unsigned char value_ = 0;
      for(int j= 0; j<imLabeledMask.getBufSize(); j++)
      {
         value_ = static_cast<unsigned char>(imLabel(j));
         imLabeledMask(j) = value_;
      }
 
     //return imRes;
     return imLabeledMask;
}

// Associate im1 and im2 to red band and green band of result image
Image<RGB> mergeBands(Image<unsigned char> &im1, Image <unsigned char> &im2)
{

    Image<RGB> imRes(im1.getSize());

    for(int i=0; i<imRes.getBufSize(); i++) {
        // red
        imRes(i)[0]=im1(i);
        // green
        imRes(i)[1]=im2(i);
        //blue
        imRes(i)[2]=0;
    }

    return imRes;
}

Image<RGB> createRestrainedRGBImage(Image<unsigned char> &im1, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax)
{
    TSize sizeX = xmax - xmin + 1;
    TSize sizeY  = ymax - ymin + 1;
    TSize sizeZ  = zmax - zmin + 1;
    std::cout << "im1.getSize() " << im1.getSize() << std::endl;
    
    Image<RGB> imRes(sizeX,sizeY,sizeZ);

    for(int x = xmin; x<xmax+1; x++)
    {	for(int y = ymin; y<ymax+1; y++)
	{
	    for(int z = zmin; z<zmax+1; z++)
	    {	// red : PET
		imRes(x,y,z)[0] = im1(x,y,z);
		// blue
		imRes(x,y,z)[1] = 0;
		// green
		imRes(x,y,z)[2] = 0;		
       	    }
  	}
    }
    return imRes;
}

Image<RGB> createRGBImage(Image<unsigned char> &im1)
{
    Image<RGB> imRes(im1.getSize());
    for(int i=0; i<imRes.getBufSize(); i++) {
        // red : pet
        imRes(i)[0]=im1(i);
        // green
        imRes(i)[1]=0;
        //blue
        imRes(i)[2]=0;
    }
    
    return imRes;
}

/**
* Image attribute filtering based on component-tree using the shaping framework **/

int main(int argc, char *argv[])
{
    if(argc!= 16)
    {
        cout<<"Usage: " << argv[0] << " <pet_image_discretized.nii> <intensity_mean> <area_min> <area_max> <height_min> <elong_min> <shaping area> <shaping contrast> <xmin> <xmax> <ymin> <ymax> <zmin> <zmax> <output_name> \n";
        exit(1);
    }
    
    clock_t t1,t2;
    t1=clock();

    /* Declaration of: **/
    int xMin = atoi(argv[9]);
    int xMax = atoi(argv[10]);
    int yMin = atoi(argv[11]);
    int yMax = atoi(argv[12]);
    int zMin = atoi(argv[13]);
    int zMax = atoi(argv[14]);


    /* - ctImage - petImage **/
    nifti_1_header petHeader;// = NULL;

    float quatern_b, quatern_c, quatern_d, qoffset_x, qoffset_y, qoffset_z, qfac = 0;
    
    /* Read input :  discretized PET SUV image (unsigned char) **/
    Image <unsigned char> petImageOri = readNifti<unsigned char>(argv[1], &petHeader);

    getSpacingAndOrigin<unsigned char>(argv[1], &quatern_b, &quatern_c, &quatern_d, &qoffset_x, &qoffset_y, &qoffset_z,&qfac, &petHeader);
    
    /* Crop image **/
    Image<unsigned char> croppedImage = petImageOri.crop(xMin,xMax+1,yMin,yMax+1, zMin, zMax+1);
    croppedImage.setSpacingX(petImageOri.getSpacingX());
    croppedImage.setSpacingY(petImageOri.getSpacingY());
    croppedImage.setSpacingZ(petImageOri.getSpacingZ());
    std::cout <<"croppedImage : Image dimensions : " << static_cast<int>(croppedImage.getSizeX()) << " " << static_cast<int>(croppedImage.getSizeY()) << " " << static_cast<int>(croppedImage.getSizeZ()) << "\n";

    /* Create multiband (RGB) image **/
    Image<RGB> PETImage = createRGBImage(croppedImage);

    /* Declaration of : **/
    // -connexity (26-adjacency)
    FlatSE connexity;
    connexity.make3DN26();
    // threshold on mean region intensity
    int intensityMean = atoi(argv[2]);
    // threshold on min and max region volume
    int areaMin=atoi(argv[3]);
    int areaMax = atoi(argv[4]);
    // threshold on height (i.e. contrast = difference between node attribute and leaf attribute)
    int contrastMin=atoi(argv[5]);
    // threshold on elongation i.e. compacity
    int elongMin=atoi(argv[6]);
    // threshold on volume in 2nd tree (i.e. number of nodes in a given node of the 2nd tree)
    int shapingArea=atoi(argv[7]);
    // threshold on contrast in 2nd tree
    int shapingContrast=atoi(argv[8]);
    const char *output_name = argv[15];
    
    std::cout<<"cgraph construction " << std::endl;
    CGraph *cgraph = new CGraph(PETImage, connexity);
    //std::cout << "number of cgraph nodes " << cgraph->graph.size() << std::endl;

    /* Track computation progress **/
    graphWatcher *myWatcher = new graphWatcher(PETImage.getBufSize());
    /* Set marginal ordering on RGB colour space **/
    ColorMarginalOrdering  *order = new ColorMarginalOrdering();
    /* Compute \ddot component-graph **/
    cgraph->computeGraph(order,myWatcher);


    /* Compute attributes (elongation, here) (the computation for other attributes is commented) **/
    std::cout<<"computeMU"<<std::endl;
    cgraph->computeMU(croppedImage.getSpacingX(), croppedImage.getSpacingY(), croppedImage.getSpacingZ());
    std::cout<<"computeElongation"<<std::endl;
    cgraph->computeElongation();
    // scgraph->writeDot("cgraph.dot");
    
     /* First tree filtering on nodes intensity (SUV) and volume for clinical use  (other possibilities are commented) **/
     cgraph->intensityFiltering(intensityMean); // (SUV < 3 ) in Fiji
     //cgraph->elongFiltering(elongMin);
     //cgraph->areaFiltering(areaMin,areaMax); // area: number of voxels
     cgraph->volumeFiltering(areaMin,areaMax); // volume in mL
     cgraph->contrastFiltering(contrastMin);

     /* Second tree construction for shaping **/
     std::cout << "computeShaping\n";
     Shaping *shaping = new Shaping(cgraph);
     shaping->computeShaping();
    
     /* Second tree attributes computation **/
     shaping->computeArea(); // number of nodes of the first tree in a node of the second tree
     shaping->computeContrast();
     //shaping->writeDot("shaping_beforeFiltering.dot");

     /* Second tree filtering :  compacity, Area and contrast filtering **/
     shaping->elongFiltering(elongMin);
     shaping->areaFiltering(shapingArea);
     shaping->contrastFilteringMax(shapingContrast); // height filtering -> if not, set to 0 in the function call
    
     /*  Write graph in dot format **/
     //shaping->writeDot("shapingFiltered.dot");

     /*   Reconstruct graph after filtering **/
     shaping->constructGraph();

     /* Write graph in dot format : Filtered nodes are shown in gray **/
     //cgraph->writeDot("cgraph.dot");
    
    /* Compute resulting image from filtered cgraph **/
    Image<RGB> imResult=shaping->cgraph->constructImage(order);

    //Image<unsigned char> petImageOriginal = keepBand(PETImage,0);
    //Image<unsigned char> petMask = getLabeledMask(petImageOriginal,connexity);
     
    Image<unsigned char> imResultSecondBand = keepBand(imResult,0);
    Image<unsigned char> imLabMask = getLabeledMask(imResultSecondBand, connexity);

  	imLabMask.setSpacingX(petImageOri.getSpacingX());
	imLabMask.setSpacingY(petImageOri.getSpacingY());
	imLabMask.setSpacingZ(petImageOri.getSpacingZ());   

    std::cout << "Mask dimensions : " << static_cast<int>(imLabMask.getSizeX()) << " " << static_cast<int>(imLabMask.getSizeY()) << " " << static_cast<int>(imLabMask.getSizeZ()) << "\n";
    Image<unsigned char> imLabMaskFinal = reconstructFullImage(petImageOri,imLabMask, xMin, xMax, yMin, yMax, zMin, zMax);

     /* Write Output **/
     writeNifti(imLabMaskFinal, output_name,&quatern_b, &quatern_c, &quatern_d, &qoffset_x, &qoffset_y, &qoffset_z, &qfac);


    t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<"Run : "<< diff / CLOCKS_PER_SEC << " seconds " << endl;    

    delete cgraph;
    delete shaping;
    return 1;
}

