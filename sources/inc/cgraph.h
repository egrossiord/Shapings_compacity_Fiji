//Copyright (C) 2013, Benoît Naegel <b.naegel@unistra.fr> and modified by Eloïse Grossiord <eloise.grossiord@gmail.com>
//This program is free software: you can use, modify and/or
//redistribute it under the terms of the GNU General Public
//License as published by the Free Software Foundation, either
//version 3 of the License, or (at your option) any later
//version. You should have received a copy of this license along
//this program. If not, see <http://www.gnu.org/licenses/>.


#ifndef CGRAPH_H
#define CGRAPH_H

#include <vector>
#include <map>

#include "Image.h"
#include "FlatSE.h"
#include "ragraph.h"
#include "cgraphwatcher.h"
#include "colorordering.h"
#include "utils.h"

using namespace LibTIM;
using namespace std;

/** Component-graph storage and computation **/

class CGraph
{
public:
    /** A node of the directed graph */
    struct Node {
        int index;
        RGB color;
        RGB dispColor;
        Node(int index, RGB color, int area) {this->index=index; this->color=color; this->area=area; contrast=0; active=true; active2=true;
        m100=m010=m001=m110=m011=m101=m200=m020=m002=mu200=mu020=mu002=mu110=mu011=mu101=sumf=sumff=SUVmax=SUVmean=SUVmin=sumfPET=sumffPET=sumfCT=sumffCT=0; elong=elong2=flatn=flatn2=0;fVolumemL=0; }
        std::vector<Node *> childs;
        std::vector<Node *> fathers;
        // list of flat-zones belonging to node and having same value
        std::vector<int > regions;
        std::vector<int> pixels;
        int area;
        float fVolumemL;
        int contrast;

        // add
        int elong;
        long double elong2;
        long double flatn;
        long double flatn2;
        long double noncompacity;

        // moments (dans chaque direction)
        long double m100;
        long double m010;
		long double m001;
		long double m110;
		long double m011;
		long double m101;
		long double m200;
		long double m020;
		long double m002;

        float sumf;
        float sumff;
        double SUVmax;
        double SUVmean;
        double SUVmin;
		double sumfPET;
		double sumffPET;
        float sumfCT;
        double sumffCT;
		long double variancePET;
        long double varianceCT;
		long double meanPET;
        long double meanCT;
        long double skewness;
        long double kurtosis;
        long double eccentricity;
        long double solidity;
        long double elongation;
        long double aspectRatio;
        long double flatness;
        
        long double sphericity;
        long double sparsness; // wilkinson
        //Image<double> m_glcm;
		// ITK PET texture features
        long double energy;
        long double entropy;
        long double inverseDifferenceMoment;
        long double inertia;
        long double clusterShade; // skewness of GLCM
        long double clusterProminence; // kurtosis of GLCM
        // ITK CT Texture Features
        long double energyCT;
        long double entropyCT;
        long double inverseDifferenceMomentCT;
        long double inertiaCT;
        long double clusterShadeCT;
        long double clusterProminenceCT;
        bool bVTTumourTag;
        bool bVTRecoverAlreadyTreated;
        U16 nTumourLabel;
        U16 ntmpTumourLabel;
        float fVTRecovery;
        int nPixelsListSize;
        // attributes for nodes selection (see paper BN/NP interactive segmentation for nodes selection in respect with GT)
        int n;
        int ps;
        double calpha;
    
		// mu
		long double mu200;
		long double mu020;
		long double mu002;
		long double mu110;
		long double mu011;
		long double mu101;

        //moment d'inertie (first Hu's invariant moment)
        double I;

        bool active;
        // if active2==false, the node will not be reconstructed in the result image
        // (used for graph filtering before shaping)
        bool active2;

        void addChild(Node *child) {
            this->childs.push_back(child);
            child->fathers.push_back(this);
        }
    };

public :
    RAGraph *rag; /*!< Region adjacency graph on which is computed the component-graph */

    Image<RGB> imSource; /*! PET Source image */
    //Image<float> imSUVPET;
    //Image<S16> imCT; /*! CT Source image */
    //Image<U16> imVTLabels; // Images de labels
    FlatSE connexity;   /*!< Underlying connexity (4- or 8-) */

    /** Container of the nodes of the (directed) graph **/
    std::vector<Node *> graph;
    /** Root of the graph **/
    Node *root; 

public:

    CGraph(Image <RGB> &imSource, FlatSE &connexity){
        this->imSource = imSource;
        this->connexity=connexity;
    
        this->rag=new RAGraph(imSource,connexity);
        }
    CGraph(RAGraph *rag) : rag(rag) {}
    ~CGraph() { delete rag;}

    Node *componentGraphNaive(FlatSE &connexity);

    /** Component-graph \ddot G **/
    int computeGraph(ColorOrdering *order, CGraphWatcher *watcher);
    /** Component-graph G **/
    int computeGraphFull(ColorOrdering *order, CGraphWatcher *watcher);

    /** Return synthetic images */
    static Image<RGB> syntheticImage();
    static Image<RGB> syntheticImage2();

    /** Write graph in dot format **/
    int writeDot(const char *filename);
    /** Write graph attributes in csv file **/
    int writeAttributes(const char *filename, char* patientName);

	/** Methods : computation of attributes */
	//int computeAttributes(Node *n);
    //int computeArea(Node *n);
	//int computeContrast(Node *n);
	//int computeVolume(Node *n);
	//int computeSubNodes(Node *n);
	//long double computeM010(Node *n);
	//long double computeM100(Node *n);
	//long double computeM001(Node *n);
	//long double computeM020(Node *n);
	//long double computeM200(Node *n);
	//void computeInertiaMoment(Node *n);
	void update_attributes(Node *n, vector< Point<int> > &pixels); // Tcoord=int
    //void update_pixel(Node *n, int x1, int y1, int z1, int x2, int y2, int z2);
	void computeMU(float spacingX, float spacingY, float spacingZ);
	void computeElongation();

	//void computeGLCM(Node *n, vector< Point<int> > &pixels, int nSize);
	//void compute_attributes_from_glcm(Node *n);

    void intensityFiltering(float intensityMin);
    // Fonction de filtrage sur le nombre de voxels minimal dans une composante
    void areaFiltering(int areaMin);
    void elongFiltering(int elongMin);

    Image<RGB> constructImage(ColorOrdering *order);
    /** Paint node index into Image**/
    Image <RGB> constructNode(int index);
    Image <RGB> constructNode2(ColorOrdering *order,int index);

    /** Helper functions **/
    vector<CGraph::Node *> computeComponents(ColorOrdering *order, OrderedQueue<RGB> &pq);
    static bool notComparable(Image<RGB> &im, Point<TCoord> &p, Point<TCoord> &q);

    bool isEqual(Node *n, Node *m);
    bool isIncluded(Node *n, Node *m);
    bool isIncludedFast(Node *n, Node *m, vector<bool> &tmp);
    void computeLinks(ColorOrdering *order, vector<Node *> &nodes);
    Node *addRoot(vector<Node *> &nodes);
    vector<Node *> minimalElements(ColorOrdering *order, vector<Node *> &nodes, vector<bool> &tmp);
    void computeTransitiveReduction(ColorOrdering *order, vector<Node *> &nodes);

    int computeGraphInverse(CGraphWatcher *watcher);
    Image<RGB> constructImageInf();
    // Fonction de filtrage sur le nombre min et max de voxels au sein d'une composante
    void areaFiltering(int areaMin, int areaMax);
    // Fonction de filtrage sur le volume (mL) min et max d'une composante
    void volumeFiltering(float volumeMin, float volumeMax);
    bool isLTE(RGB &v, RGB &w);
    Image<RGB> constructImageInverse();

    void contrastFiltering(int contrastMin);
    void contrastFiltering(int contrastMin, int contrastMax);
    void resetFiltering();

    /** Adaptive filtering **/
    void adaptiveAreaFiltering(int p);
    void adaptiveContrastFiltering(int p);

    bool isIncluded(Node *n, Node *m, vector<bool> &tmp);
    Image<RGB> constructImageTest(ColorOrdering *order);
private:
    void paintNode(Image<RGB> &im, Node *n, RGB &value);
    void paintNodeSup(Image<RGB> &imRes, Node *n, RGB &value);

    vector<Node *> computeComponents(Image<RGB> &im, FlatSE &connexity);

};


#endif // CGRAPH_H
