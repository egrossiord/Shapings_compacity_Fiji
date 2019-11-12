//Copyright (C) 2013, Benoît Naegel <b.naegel@unistra.fr>, modified by Eloïse Grossiord <eloise.grossiord@gmail.com>
//This program is free software: you can use, modify and/or
//redistribute it under the terms of the GNU General Public
//License as published by the Free Software Foundation, either
//version 3 of the License, or (at your option) any later
//version. You should have received a copy of this license along
//this program. If not, see <http://www.gnu.org/licenses/>.


#include "ragraph.h"
#include "cgraph.h"
#include "utils.h"
#include <fstream>
#include <glob.h>
//#include <stdint.h>
#include <limits>

using namespace std;
#define INT8_MAX 127



inline std::vector<std::string> glob(const std::string& path){
    using namespace std;
    glob_t glob_result;
    glob(path.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}


Image<RGB> CGraph::syntheticImage()
{
    int a[]={0,0,0};
    int b[]={1,0,0};
    int c[]={0,1,0};
    int d[]={1,1,0};
    int e[]={1,1,1};
    int f[]={0,0,1};

    RGB A(a);
    RGB B(b);
    RGB C(c);
    RGB D(d);
    RGB E(e);
    RGB F(f);

    TSize size[]={7,1,1};
    TSpacing spacing[]={1.0,1.0,1.0};
    const RGB data[]={ A,E,C,F,B,E,A};
    Image<RGB> imTmp(size,spacing,data);
    return imTmp;
}


Image<RGB> CGraph::syntheticImage2()
{
    int a[]={0,0,0};
    int b[]={1,0,0};
    int c[]={0,1,0};
    int d[]={1,1,0};
    int e[]={1,1,1};
    int f[]={0,0,1};

    RGB A(a);
    RGB B(b);
    RGB C(c);
    RGB D(d);
    RGB E(e);
    RGB F(f);

    TSize size[]={21,1,1};
    TSpacing spacing[]={1.0,1.0,1.0};
    const RGB data[]={  A,C,C,C,B,E,D,D,C,E,E,D,B,E,B,C,C,C,B,B,A};

    Image<RGB> imTmp(size,spacing,data);
    return imTmp;
}


/**
*   Computation of component-graph \ddot G
**/
int CGraph::computeGraph(ColorOrdering *order, CGraphWatcher *watcher)
{

    OrderedQueue <int> pq;

    vector<bool> processed(rag->nodes.size());

    /** Warning : nodes[0] does not exist !!!! */
    for(int i=1; i<rag->nodes.size(); i++) {

        RGB value=rag->nodes[i]->color;

        int priority = order->getPriority(value);

        // put in the prior queue the flat zone number
        pq.put(priority,i);
        processed[i]=false;
    }
    //std::cout<<"end sorting" << std::endl;

    // index from flat-zone number to corresponding node
    vector <Node *> regionToNode(rag->nodes.size());
    regionToNode.assign(regionToNode.size(),0);

    std::queue<int> fifo;
    std::vector<bool> infifo(regionToNode.size());

    while(!pq.empty() ) {
        int i=pq.get();

        if(watcher!=0)
            {
            watcher->progressUpdate();
            }

        if(regionToNode[i]==0) {
            /** i is a canonic region (i==regionToNode[i]->index- */
            RAGraph::Vertex *curVertex=rag->nodes[i];
            Node *curNode = new Node(curVertex->index,curVertex->color,curVertex->pixels.size());
	    update_attributes(curNode,curVertex->pixels);
			//computeGLCM(curNode, curVertex->pixels, 64);
           

            curNode->regions.push_back(i);
            regionToNode[i]=curNode;

            int valueMax = static_cast<int>(curNode->color[0]+curNode->color[1]+curNode->color[2]);

            /** on lance la propagation */
            infifo.assign(infifo.size(), false);
            fifo.push(i);
            infifo[i]=true;
            int M=0;
            while(!fifo.empty()) {
                int curpt=fifo.front();
                M++;
                fifo.pop();
                RAGraph::Vertex *tmpVertex=rag->nodes[curpt];
                for(int v=0; v<tmpVertex->allNb.size(); v++ ) {
                    int nb=tmpVertex->allNb[v];
                    if(!infifo[nb] && order->islessequal(rag->nodes[i]->color , rag->nodes[nb]->color)) {
                         curNode->regions.push_back(nb);
                        if(order->isequal(rag->nodes[i]->color , rag->nodes[nb]->color)) {
                            /** ZP qui appartient au même noeud : on fusionne*/

                            regionToNode[nb]=curNode;

                            //std::cout << "Fusion de " << i << " et " << nb << "\n";
                        }
                        else{
                            /** zp potentiellement descendant du noeud courant : on teste si un père du noeud est comparable avec curNode
                        *  Si c'est le cas, a n'est pas un fils direct de curNode
                        */
                            Node *tmp=regionToNode[nb];

                            valueMax = std::max(valueMax,static_cast<int>(tmp->color[0]+tmp->color[1]+tmp->color[2]));

                            assert(tmp!=0);

                            bool isChild=true;
                            for(int a=0; a<tmp->fathers.size(); a++) {
                                if(order->islessequal( curNode->color , tmp->fathers[a]->color) ) {
                                    isChild=false;
                                    break;
                                }
                            }
                            if(isChild) {
                                curNode->childs.push_back(tmp);
                                tmp->fathers.push_back(curNode);

                            }

                        }

                        // Calcul des attributs necessaires au cours de la construction
                        curNode->area+=rag->nodes[nb]->pixels.size();
                        //std::cout<<"curNode->area "<<curNode->area<<std::endl;
                        curNode->contrast=valueMax-(curNode->color[0]+curNode->color[1]+curNode->color[2]);
                        assert(curNode->contrast >= 0);
                        update_attributes(curNode,rag->nodes[nb]->pixels);
                        //computeGLCM(curNode, rag->nodes[nb]->pixels,64);

                        fifo.push(nb);
                        infifo[nb]=true;
                    }
                }


            }
            processed[i]=true;
        }
    }

    /** Assign a new index to each node */
    //std::cout << "assign new index to each node" << std::endl;
    for(int i=0; i<regionToNode.size(); i++) {
        if(regionToNode[i]!=0 && regionToNode[i]->index==i) {
            graph.push_back(regionToNode[i]);
            graph[graph.size()-1]->index=graph.size()-1;
        }

    }

    /** Add fictitious root*/
    //std::cout << "Add fictitious root" << std::endl;
    RGB value(0,0,0);
    //BN FIX : on ajoute le noeud fictif à l'index graph.size() (dernier noeud de la liste)
    //root=new Node(-1,value,imSource.getBufSize());

    root=new Node(graph.size(),value,imSource.getBufSize());
    graph.push_back(root);

        for (int i=0; i<graph.size(); i++) {

            if(graph[i]->fathers.size()==0) {
                root->addChild(graph[i]);
            }
        }

    return 1;
}



bool CGraph::isLTE(RGB &v, RGB &w)
{
    if(v[0]>=w[0] && v[1]>=w[1] && v[2]>=w[2] )
        return true;
    else return false;
}

/**
*   Computation of component-graph \ddot G using inverse order
**/
int CGraph::computeGraphInverse(CGraphWatcher *watcher) {

    OrderedQueue <int> pq;

    vector<bool> processed(rag->nodes.size());

    //std::cout << rag->nodes.size() << "\n";


    /** Warning : nodes[0] does not exist !!!! */
    for(int i=1; i<rag->nodes.size(); i++) {

        RGB value=rag->nodes[i]->color;

        int R=value[0];
        int G=value[1];
        int B=value[2];
        // put in the prior. queue the flat zone number
        pq.put((R+G+B),i);
        processed[i]=false;
    }

    // index from flat-zone number to corresponding node
    vector <Node *> regionToNode(rag->nodes.size());
    regionToNode.assign(regionToNode.size(),0);

    std::queue<int> fifo;
    std::vector<bool> infifo(regionToNode.size());

    std::vector<Node *> potentialChilds;

    while(!pq.empty() ) {
        int i=pq.get();

        if(watcher!=0)
            {
            watcher->progressUpdate();
            }
        //std::cout << "Région visitée : " << i << "\n";

        if(regionToNode[i]==0) {
            /** i est une région canonique (i==regionToNode[i]->index- */
            RAGraph::Vertex *curVertex=rag->nodes[i];
            Node *curNode=new Node(curVertex->index,curVertex->color,curVertex->pixels.size());

            curNode->regions.push_back(i);
            regionToNode[i]=curNode;

            potentialChilds.clear();

            /** on lance la propagation */
            infifo.assign(infifo.size(), false);
            fifo.push(i);
            infifo[i]=true;
            int M=0;
            while(!fifo.empty()) {
                int curpt=fifo.front();
                M++;
                fifo.pop();
                RAGraph::Vertex *tmpVertex=rag->nodes[curpt];
                for(int v=0; v<tmpVertex->allNb.size(); v++ ) {
                    int nb=tmpVertex->allNb[v];
                    if(!infifo[nb] && isLTE(rag->nodes[i]->color, rag->nodes[nb]->color))
                    {
                        if(rag->nodes[i]->color == rag->nodes[nb]->color) {
                        /** ZP qui appartient au même noeud : on fusionne*/
                        curNode->regions.push_back(nb);
                        regionToNode[nb]=curNode;

                        /** Mise à jour de l'aire */

                        //std::cout << "Fusion de " << i << " et " << nb << "\n";
                        }
                        else {
                        /** zp potentiellement descendant du noeud courant : on teste si un père du noeud est comparable avec curNode
                        *  Si c'est le cas, a n'est pas un fils direct de curNode
                        */
                            Node *tmp=regionToNode[nb];
                            bool isChild=true;
                            for(int a=0; a<tmp->fathers.size(); a++) {
                                if(isLTE(curNode->color , tmp->fathers[a]->color) ) {
                                    isChild=false;
                                    break;
                                }
                            }
                            if(isChild) {
                                curNode->childs.push_back(tmp);
                                tmp->fathers.push_back(curNode);

                            }
                        //potentialChilds.push_back(tmp);
                        //std::cout << "\t Ajoute " << nb << "\n";
                        }
                        curNode->area+=rag->nodes[nb]->pixels.size();
                        fifo.push(nb);
                        infifo[nb]=true;
                    }
                }


            }
        processed[i]=true;
        }
    }

    /** Assign a new index to each node */
    for(int i=0; i<regionToNode.size(); i++) {
        if(regionToNode[i]!=0 && regionToNode[i]->index==i) {
            graph.push_back(regionToNode[i]);
            graph[graph.size()-1]->index=graph.size()-1;
        }

    }

    /** Add fictitious root*/

    RGB value(255,255,255);

    root=new Node(-1,value,imSource.getBufSize());

        for (int i=0; i<graph.size(); i++) {

            if(graph[i]->fathers.size()==0) {
                root->addChild(graph[i]);
            }
        }

    return 1;
}



void CGraph::intensityFiltering(float intensityMin)
{
    int nEliminatedNodes = 0;
    for(int i= 0; i<graph.size()-1; i++)
    {
	Node *n = graph[i];
 	float mean = n->sumf/n->area;
 	if(mean < intensityMin)
	{
	    n->active = false;
	    n->active2 = false;
	    nEliminatedNodes++;
	}	
	
    }
    return;
}

void CGraph::elongFiltering(int elongMin)
{
    for(int i=0; i<graph.size()-1; i++)
    {
        Node *n=graph[i];

        if(n->elong<elongMin) {n->active=false; n->active2=false;}
    }
    return;
}


void CGraph::areaFiltering(int areaMin)
{

    for(int i=0; i<graph.size()-1; i++) {
        Node *n=graph[i];
        if(n->area<areaMin) {
            n->active=false;
            n->active2=false;
            }
    }
    
    return;

}

void CGraph::volumeFiltering(float volumeMin, float volumeMax)
{
	/* Filtering on nodes from the graph based on a volume attribute : */
	/* nodes with volume (mL) < volumeMin or volume (mL) > volumeMax are filtered out*/
 	std::cout << " Filtering of nodes with volume lower than " << volumeMin << " and volume greater than " << volumeMax << std::endl;
	for(int i=0; i<graph.size()-1; i++)
	{
		Node *n = graph[i];
		if(n->fVolumemL<volumeMin || n->fVolumemL>volumeMax)
		{
			n->active = false; // supprime du premier arbre
			n->active2 = false; // supprime du deuxieme arbre
		}	
	}
	return;
}



void CGraph::areaFiltering(int areaMin, int areaMax)
{


    for(int i=0; i<graph.size()-1; i++) {
        Node *n=graph[i];
        if(n->area<areaMin || n->area>areaMax) {
            n->active=false;
            n->active2=false;
        }
    }

   // if(root->area<areaMin || root->area>areaMax) root->active=false;
    return;
}


void CGraph::contrastFiltering(int contrastMin)
{

    for(int i=0; i<graph.size()-1; i++) {
        Node *n=graph[i];
        if(n->contrast<contrastMin) {
            n->active=false; n->active2=false;}
    }

    //if(root->contrast<contrastMin) {root->active=false; root->active2=false;}
    return;
}


void CGraph::resetFiltering()
{
/** Initialisation*/
    for(int i=0; i<graph.size(); i++) {graph[i]->active=true;}
    
    return;
}

/**
    Keep the p % nodes having the biggest area
*/
void CGraph::adaptiveAreaFiltering(int p)
{


    std::map<int, int> histo;
    for(int i=0; i<graph.size(); i++) {
        Node *n=graph[i];
        histo[-n->area]++;
    }

    histo[-root->area]++;

    int totalNodes=graph.size()+1;
    int bestNodes=(p*totalNodes)/100;

    int count=0;

    int threshold;

    std::map<int,int>::iterator it;
    for(it=histo.begin(); it!=histo.end(); ++it) {
        count+=it->second;
        if(count> bestNodes )
            {
            threshold=-it->first;
            break;
            }
    }
    areaFiltering(threshold);
    
    return;

}


/**
    Keep the p % nodes having the biggest area
*/
void CGraph::adaptiveContrastFiltering(int p)
{
    std::map<int, int> histo;
    for(int i=0; i<graph.size(); i++) {
        Node *n=graph[i];
        histo[-n->contrast]++;
    }

    histo[-root->contrast]++;

    int totalNodes=graph.size()+1;
    int bestNodes=(p*totalNodes)/100;

    int count=0;

    int threshold;

    std::map<int,int>::iterator it;
    for(it=histo.begin(); it!=histo.end(); ++it) {
        count+=it->second;
        if(count> bestNodes )
            {
            threshold=-it->first;
            break;
            }
    }
    contrastFiltering(threshold);
    
    return;

}


void CGraph::contrastFiltering(int contrastMin, int contrastMax)
{


    for(int i=0; i<graph.size(); i++) {
        Node *n=graph[i];
        if(n->contrast<contrastMin || n->contrast>contrastMax) n->active=false;
    }

    if(root->contrast<contrastMin || root->contrast>contrastMax) root->active=false;
    
    return;
}

Image<RGB> CGraph::constructImageInf()
{
    
    std::queue<Node *> fifo;

    vector <bool> active(graph.size());
    active.assign(active.size(),true);

    OrderedQueue <Node *> pq;

    vector<bool> processed(graph.size());

    RGB rgbmin(0,0,0);

    /** Warning : nodes[0] does not exist !!!! */
    for(int i=0; i<graph.size(); i++) {

        graph[i]->dispColor=rgbmin;

        RGB value=graph[i]->color;

        int R=value[0];
        int G=value[1];
        int B=value[2];
        pq.put(-(R+G+B),graph[i]);
        processed[i]=false;
    }

    /** Special case for the root*/
    root->dispColor=rgbmin;

    while(!pq.empty()) {
        Node *n=pq.get();

        if(n->active==false && n->fathers.size()>1) {

            int nactive=0;
            RGB value(255,255,255);
            for(int j=0; j<n->fathers.size(); j++) {
                if(n->fathers[j]->active==true)
                {
                    nactive++;
                    value[0]=std::min(value[0],n->fathers[j]->color[0]);
                    value[1]=std::min(value[1],n->fathers[j]->color[1]);
                    value[2]=std::min(value[2],n->fathers[j]->color[2]);
                }
            }
            //qDebug() << "Node: " << n->label << " "<<  nactive << "\n";
            if(nactive>1) {
                for(int j=0; j<n->fathers.size(); j++) {

                    if(n->fathers[j]->active==true) {
                        n->fathers[j]->dispColor=value;
                    }
                }
            }
        }
    }

    Image<RGB> imRes=imSource;
    imRes.fill(0);

    for(int i=0; i<graph.size(); i++) {


        RGB value=graph[i]->color;

        int R=value[0];
        int G=value[1];
        int B=value[2];
        pq.put((R+G+B),graph[i]);
        processed[i]=false;
    }

    /** Special case for the root*/
    pq.put(-1,root);

    while(!pq.empty()) {
        Node *curNode=pq.get();
        if(curNode->active==true && curNode->dispColor!=rgbmin) {
            
            curNode->dispColor=curNode->color;
        }

        paintNode(imRes,curNode, curNode->dispColor);

        for(int c=0; c<curNode->childs.size(); c++) {
            Node *curChild=curNode->childs[c];
            if(curChild->active==false) {
//            if(curNode->dispColor>curChild->dispColor)
               curChild->dispColor=curNode->dispColor;
            }

        }
    }

    return imRes;

}

void CGraph::paintNode(Image <RGB> &imRes, Node *n, RGB &value)
{
    for(int r=0; r<n->regions.size(); r++) 
    {
        vector<Point<TCoord> > points=rag->nodes[n->regions[r]]->pixels;
        for(int p=0; p<points.size(); p++)
        {
                imRes(points[p])=value;
        }
    }
    return;
}


void CGraph::paintNodeSup(Image <RGB> &imRes, Node *n, RGB &value)
{
    for(int r=0; r<n->regions.size(); r++) 
    {
        assert(n->regions[r]!=0);
        vector<Point<TCoord> > points=rag->nodes[n->regions[r]]->pixels;
        for(int p=0; p<points.size(); p++) 
	{
            RGB imValue = imRes(points[p]);
            imValue[0]=std::max(imValue[0],value[0]);
            imValue[1]=std::max(imValue[1],value[1]);
            imValue[2]=std::max(imValue[2],value[2]);
	    
            imRes(points[p])=imValue;
        }
    }
    
    return;
}

 
Image<RGB> CGraph::constructImage(ColorOrdering *order)
{
    Image<RGB> imRes=imSource;
    imRes.fill(0);
    std::queue<Node *> fifo;

    OrderedQueue <Node *> pq;

    vector<bool> processed(graph.size());

    RGB rgbmin(0,0,0);

    for(int i=0; i<graph.size()-1; i++) {

        RGB value=graph[i]->color;

        pq.put(-order->getPriority(value),graph[i]); // rangement dans pq des noeuds de color par ordre croissant de la somme des 3 bandes (i.e. des min aux max d'intensite)
        processed[i]=false;
    }

    /** Special case for the root*/
    pq.put(-1,root);

    assert(root->active==true);
    while(!pq.empty()) 
    {
        Node *curNode=pq.get();
        if(curNode->active==true ) {
            
            curNode->dispColor=curNode->color;
        }
        
        paintNode(imRes,curNode, curNode->dispColor);

	// on parcours les fils de curNode pour les remonter a la valeur de curNode
	if(curNode->childs.size() != 0)
	{
	    for(int c = 0; c<curNode->childs.size(); c++)
	    {
		Node *curChild=curNode->childs[c];
		if(curChild->active==false)
		{
		    curChild->dispColor = curNode->dispColor;
		}
	    }
	}

    }
    return imRes;
}


Image<RGB> CGraph::constructImageTest(ColorOrdering *order)
{
    Image<RGB> imRes=imSource;
    imRes.fill(0);
    std::queue<Node *> fifo;

    OrderedQueue <Node *> pq;

    vector<bool> processed(graph.size());

    RGB rgbmin(0,0,0);

    for(int i=0; i<graph.size()-1; i++) {

        RGB value=graph[i]->color;

        pq.put(-order->getPriority(value),graph[i]); // rangement dans pq des noeuds de color par ordre croissant de la somme des 3 bandes (i.e. des min aux max d'intensite)
        processed[i]=false;
    }

    /** Special case for the root*/
    pq.put(-1,root);


    while(!pq.empty()) {
        Node *curNode=pq.get();
        if(curNode->active==true ) {
           curNode->dispColor=curNode->color;
        }
        // BN FIX
        std::cout <<"Paint node " << curNode->index <<  "\n";
        paintNode(imRes,curNode, curNode->color);

        // on parcours les fils de curNode pour les remonter a la valeur de curNode
        for(int c=0; c<curNode->childs.size(); c++) {
            Node *curChild=curNode->childs[c];
            if(curChild->active==false) {
            if(curNode->dispColor>curChild->dispColor)
               curChild->dispColor=curNode->dispColor;
            }

        }
    }

    return imRes;
}

Image<RGB> CGraph::constructNode(int index)
{
    Image<RGB> imRes(this->imSource.getSize());
    imRes.fill(0);
    Node *n=graph[index];

    RGB color(n->color[0],n->color[1],0);
    paintNode(imRes,graph[index],color);
    return imRes;
}

Image<RGB> CGraph::constructNode2(ColorOrdering *order,int index)
{
    Image<RGB> imRes(this->imSource.getSize());
    imRes.fill(0);
    for(int i=0; i<graph.size(); i++)
        if(i==index)
            graph[i]->active=true;
        else
            graph[i]->active=false;
    graph[graph.size()-1]->active=true;

    imRes=this->constructImage(order);
    return imRes;
}


Image<RGB> CGraph::constructImageInverse()
{
    Image<RGB> imRes=imSource;
    imRes.fill(0);
    std::queue<Node *> fifo;

    OrderedQueue <Node *> pq;

    vector<bool> processed(graph.size());

    RGB rgbmin(0,0,0);

    /** Warning : nodes[0] does not exist !!!! */
    for(int i=0; i<graph.size(); i++) {

        graph[i]->dispColor=rgbmin;

        RGB value=graph[i]->color;

        int R=value[0];
        int G=value[1];
        int B=value[2];
        pq.put(-(R+G+B),graph[i]);
        processed[i]=false;
    }

    /** Special case for the root*/
    pq.put(-1,root);


    while(!pq.empty()) {
        Node *curNode=pq.get();

        if(curNode->active==true ) {
            curNode->dispColor=curNode->color;
        }

        paintNode(imRes,curNode, curNode->dispColor);

        for(int c=0; c<curNode->childs.size(); c++) {
            Node *curChild=curNode->childs[c];
            if(curChild->active==false) {
//            if(curNode->dispColor>curChild->dispColor)
               curChild->dispColor=curNode->dispColor;
            }

        }
    }

    return imRes;
}



bool CGraph::notComparable(Image<RGB> &im, Point<TCoord> &p, Point<TCoord> &q)
{
    if(!(im(p)<im(q) ) && !(im(p)>im(q)) && !(im(p)==im(q)))
        return true;
    else return false;

}

//int CGraph::writeAttributes(const char *filename, char* patientName)
//{
//    
//    if(0 != root)
//    {
//        std::ofstream outputFile(filename,std::ios::out | std::ios::app );
//        if(!outputFile)
//        {
//            std::cerr << "File writeAttributes error \n";
//            return 0;
//        }
//    
//        
//        int labelMax = graph.size()+1;
//        bool isActive[labelMax];
//        for(int i=0; i<labelMax; i++) isActive[i]=true;
//        bool isWritten[labelMax];
//        for(int j = 0; j<labelMax; j++) isWritten[j] = false;
//        
//        std::cout << string(patientName) << std::endl;
//
//        std::queue<Node *> qNodes;
//        qNodes.push(root);
//        
//        while(!qNodes.empty())
//        {
//            
//            Node *n = qNodes.front();
//            qNodes.pop();
//            std::cout << "n->index " << n->index << std::endl;
//            std::stringstream nodeAttribute;
//            nodeAttribute.str(" "); // on vide nodeAttribute avt de stocker les nouvelles valeurs
//            nodeAttribute << string(patientName) << ";" << n->index << ";" << n->SUVmax << ";" << n->SUVmin << ";" << n->SUVmean << ";" << n->area << ";" << n->contrast << ";" << n->elong << ";" << n->elong2 << ";" << n->flatn << ";" << n->flatn2 << ";" << n->noncompacity << ";" << n->eccentricity << ";" << n->sphericity << ";" << n->elongation << ";" << n->aspectRatio << ";" << n->flatness << ";" << n->sparsness << ";"<< n->meanPET << ";" << n->variancePET << ";" << n->meanCT << ";" << n->varianceCT << ";" << n->skewness << ";" << n->kurtosis << ";" << n->energy << ";" << n->entropy << ";" << n->inverseDifferenceMoment << ";" << n->inertia << ";" << n->clusterShade << ";" << n->clusterProminence << ";" << n->energyCT << ";" << n->entropyCT << ";" << n->inverseDifferenceMomentCT << ";" << n->inertiaCT << ";" << n->clusterShadeCT << ";" << n->clusterProminenceCT << ";" << n->bVTTumourTag << std::endl;
//            
//            if(!isWritten[n->index])
//            {
//
//                outputFile << std::string(nodeAttribute.str()) << ";" << std::endl;
//                std::cout << std::string(nodeAttribute.str()) << ";" << std::endl;
//                isWritten[n->index] = true;
//            }
//                
//            //Pareil pour enfant
//                for(int c = 0; c < n->childs.size(); c++)
//                {
//                    
//                    if(n->childs[c]->index != n->index)
//                    {
//                        std::stringstream nodeAttributeChild;
//                        nodeAttributeChild.str(""); // on vide
//                    
//                        nodeAttributeChild << string(patientName) << ";" << n->childs[c]->index << ";" << n->childs[c]->SUVmax << ";" << n->childs[c]->SUVmin << ";" << n->childs[c]->SUVmean << ";" << n->childs[c]->area << ";" << n->childs[c]->contrast << ";" << n->childs[c]->elong << ";" << n->childs[c]->elong2 << ";" << n->childs[c]->flatn << ";" << n->childs[c]->flatn2 << ";" << n->childs[c]->noncompacity << ";" << n->childs[c]->eccentricity << ";" << n->childs[c]->sphericity << ";" << n->childs[c]->elongation << ";" << n->childs[c]->aspectRatio << ";" << n->childs[c]->flatness << ";" << n->childs[c]->sparsness << ";"<< n->childs[c]->meanPET << ";" << n->childs[c]->variancePET << ";" << n->childs[c]->meanCT << ";" << n->childs[c]->varianceCT << ";" << n->childs[c]->skewness << ";" << n->childs[c]->kurtosis << ";" << n->childs[c]->energy << ";" << n->childs[c]->entropy << ";" << n->childs[c]->inverseDifferenceMoment << ";" << n->childs[c]->inertia << ";" << n->childs[c]->clusterShade << ";" << n->childs[c]->clusterProminence << ";" << n->childs[c]->energyCT << ";" << n->childs[c]->entropyCT << ";" << n->childs[c]->inverseDifferenceMomentCT << ";" << n->childs[c]->inertiaCT << ";" << n->childs[c]->clusterShadeCT << ";" << n->childs[c]->clusterProminenceCT << ";" << n->childs[c]->bVTTumourTag << std::endl;
//                        
//                        if(!isWritten[n->childs[c]->index])
//                        {
//                            outputFile << std::string(nodeAttributeChild.str()) <<  ";" << std::endl;
//                            std::cout << std::string(nodeAttributeChild.str()) << ";" << std::endl;
//                            isWritten[n->childs[c]->index] = true;
//                        }
//                        
//                        
//                        if(isActive[n->childs[c]->index]==true)
//                        {
//                            qNodes.push(n->childs[c]);
//                            isActive[n->childs[c]->index]=false;
//                        }
//                    }
//                }
//        }
//        
//        outputFile.close();
//        return 1;
//    }
//    
//    else
//    {
//        return 0;
//    }
//   
//}

int CGraph::writeDot(const char *filename)
{
    if(root!=0)
    {
        std::ofstream file(filename,std::ios_base::trunc  | std::ios_base::binary);
        if(!file)
        {
            std::cerr << "File I/O error\n";
            return 0;
        }

        /** Maximum label of the graph **/
        int labelMax=graph.size()+1;
        bool isActive[labelMax];
        for(int i=0; i<labelMax; i++) isActive[i]=true;

        file << "digraph G {\n";

        std::queue <Node *> fifo;
        fifo.push(root);
        while(!fifo.empty() )
        {
            Node *tmp=fifo.front();
            fifo.pop();

            std::stringstream nodeName;
            nodeName.str("");
            nodeName << "\"" << tmp->index << ":(" <<(int)tmp->color[0] << "," <<
                        (int)tmp->color[1]<<","<<
                        (int)tmp->color[2]<< " e=" << tmp->elong << ")\"";
            if(!tmp->active)
                file << "\t" << nodeName.str() << "[style=filled, fillcolor=grey];\n";

            for(int i=0; i<tmp->childs.size(); i++)
            {
                std::stringstream nodeNameChild;
                nodeNameChild << "\"" << tmp->childs[i]->index << ":(" <<(int)tmp->childs[i]->color[0] << "," <<
                                 (int)tmp->childs[i]->color[1]<<","<<
                                 (int)tmp->childs[i]->color[2]<< " e=" << tmp->childs[i]->elong  << ")\"";

                file << "\t" <<
                        nodeName.str() << "->" << nodeNameChild.str() << ";\n";

                if(isActive[tmp->childs[i]->index]==true)
                {
                    fifo.push(tmp->childs[i]);
                    isActive[tmp->childs[i]->index]=false;
                }
            }

        }

        file << "}\n";

        file.close();
        return 1;
    }
    else
        return 0;
}



/** check if set n is equal to m */
bool CGraph::isEqual(Node *n, Node *m) {

    if(n->area!=m->area) return false;

    for(int i=0; i<n->pixels.size(); i++) {
        TOffset curPixel=n->pixels[i];
        bool curPixelIncluded=false;
        for(int j=0; j<m->pixels.size(); j++) {
            if(curPixel==m->pixels[j])
                {
                curPixelIncluded=true;
                }
        }
        if(curPixelIncluded==false) return false;
    }

    return true;
}



/** Compute Attributes
*/

//int CGraph::computeArea(Node *n)
//{
	//if(n!=0)
	//{
		//Node::ContainerChilds::iterator it;
		//for(it=n->childs.begin(); it!=n->childs.end(); ++it)
			//{
				//n->area+=computeArea(*it);
			//}
		//return n->area;
	//}
	//// error
	//else return -1;
//}


//long double CGraph::computeM100(Node *n)
//{
	//if(n!=0)
	//{
		//Node::ContainerChilds::iterator it;
		//for(it=n->childs.begin(); it!=n->childs.end(); ++it)
		//{
			//n->m100+=computeM100(*it);
		//}
		//return n->m100;
	//}
	//// error
	//else return -1;
//}

//long double CGraph::computeM010(Node *n)
//{
	//if(n!=0)
	//{
		//Node::ContainerChilds::iterator it;
		//for(it=n->childs.begin(); it!=n->childs.end(); ++it)
		//{
			//n->m010+=computeM01(*it);
		//}
		//return n->m010;
	//}
	//// error
	//else return -1;
//}

//long double CGraph::computeM001(Node *n)
//{
	//if(n!=0)
	//{
		//Node::ContainerChilds::iterator it;
		//for(it=n->childs.begin(); it!=n->childs.end(); ++it)
		//{
			//n->m001+=computeM001(*it);
		//}
		//return n->m001;
	//}
	//// error
	//else return -1;
//}

//long double CGraph::computeM200(Node *n)
//{
	//if(n!=0)
	//{
		//Node::ContainerChilds::iterator it;
		//for(it=n->childs.begin(); it!=n->childs.end(); ++it)
		//{
			//n->m200+=computeM200(*it);
		//}
		//return n->m200;
	//}
	//// error
	//else return -1;
//}

//long double CGraph::computeM020(Node *n)
//{
	//if(n!=0)
	//{
		//Node::ContainerChilds::iterator it;
		//for(it=n->childs.begin(); it!=n->childs.end(); ++it)
		//{
			//n->m020+=computeM020(*it);
		//}
		//return n->m020;
	//}
	//// error
	//else return -1;
//}

//long double CGraph::computeM002(Node *n)
//{
	//if(n!=0)
	//{
		//Node::ContainerChilds::iterator it;
		//for(it=n->childs.begin(); it!=n->childs.end(); ++it)
		//{
			//n->m002+=computeM002(*it);
		//}
		//return n->m002;
	//}
	//// error
	//else return -1;
//}


//void CGraph::computeInertiaMoment(Node *n)
//{
	//std::queue<Node *> fifo;
	//fifo.push(n);

	//while(!fifo.empty())
	//{
		//Node *tmp=fifo.front();
		//fifo.pop();

		//long double xmoy=tmp->m100/tmp->area;
		//long double ymoy=tmp->m010/tmp->area;
		//long double zmoy=tmp->m001/tmp->area;

		//long double eta200=(tmp->m200-xmoy*tmp->m100)/(tmp->area*tmp->area);
		//long double eta020=(tmp->m020-ymoy*tmp->m010)/(tmp->area*tmp->area);
		//long double eta002=(tmp->m002-zmoy*tmp->m001)/(tmp->area*tmp->area);
		//tmp->I=100*(eta200+eta020+eta002);

		//std::vector<Node *>::iterator it;
		//for(it=tmp->childs.begin(); it!=tmp->childs.end(); ++it)
				//fifo.push(*it);
	//}
//}



//void CGraph::computeAttributes(Node *n)
//{
	//if(n!=0)
	//{
		////n->area=computeArea(n);
        ////n->contrast=computeContrast(n);
        ////n->volume=computeVolume(n);

        ////computeComplexityAndCompacity(n);
        ////computeBoundingBox(n);
        ////n->subNodes=computeSubNodes(n);
        ////n->m010=computeM010(n);
        ////n->m100=computeM100(n);
        ////n->m200=computeM200(n);
        ////n->m020=computeM020(n);
        ////n->m002=computeM002(n);
        ////computeInertiaMoment(n);
        //computeElongation(n);
	//}
//}

//void CGraph::update_pixel(Node *n, int x1, int y1, int z1, int x2, int y2, int z2)
//{
//    // Check si neighbour (x2,y2,z2) existe
//    if((x2<0) or (x2>imSource.getSizeX()) or (y2<0) or (y2<imSource.getSizeY()) or (z2<0) or (z2>imSource.getSizeZ()))
//        return;
//
//    unsigned char pixelValue;
//    unsigned char neighbourValue;
//
//    // valeurs du point et de son voisin dans l'image TEP
//    pixelValue = imSource(x1,y1,z1)[1];
//    neighbourValue = imSoruce(x2,y2,z2)[1];
//
//    // on se place dans la glcm a l'indice (pixelValue,neighbourValue)
//    n->m_glcm(pixelValue, neighbourValue)++;
//}



void CGraph::update_attributes(Node *n, vector< Point<int> > &pixels )
{
	//std::cout<<"nb de pixels du noeud : "<<pixels.size()<<std::endl;
    //n->sumfPET = 0;
    //n->sumffPET = 0;
    //n->sumfCT = 0;
    //n->sumffCT = 0;
   
    //n->sumf = 0;
    //n->sumff = 0; 

	for(int i=0; i<pixels.size(); ++i)
	{
		Point<int> imCoord = pixels[i];

		n->m100 += imCoord.x; //sumX
		n->m010 += imCoord.y; //sumY
		n->m001 += imCoord.z; //sumZ

		n->m200 += imCoord.x*imCoord.x; //sumXX
		n->m020 += imCoord.y*imCoord.y; //sumYY
		n->m002 += imCoord.z*imCoord.z; //sumZZ

		n->m110 += imCoord.x*imCoord.y; //sumXY
		n->m101 += imCoord.x*imCoord.z; //sumXZ
		n->m011 += imCoord.y*imCoord.z; //sumYZ

		n->sumf += (imSource(imCoord)[0]);
	  	n->sumff += imSource(imCoord)[0]*imSource(imCoord)[0];
        //n->sumfPET += float(imSUVPET(imCoord));
		//n->sumffPET += float(imSUVPET(imCoord))*float(imSUVPET(imCoord)); // somme des intensites
        
        //n->sumfCT += imCT(imCoord);
        //n->sumffCT += imCT(imCoord)*imCT(imCoord);

    }
    
    return;
}



//void CGraph::compute_attributes_from_GLCM(Node *n)
//{
//
//    n->entropy = compute_entropy(n);
//
//}

///** Entropy : mesure le desordre dans l'image et augmente lorsque la texture est aleatoire*/
//double compute_entropy(Node *n)
//{
//
//    double entropy = 0;
//    if(0 !=n)
//    {
//        for(int x=0; x<nNg; x++)
//        {
//            for(int y =0; y<nNg; y++)
//            {
//                entropy += (n->m_glcm(y,x)*std::log(n->m_glcm(y,x));
//            }
//        }
//    }
//
//    return -entropy;
//}

///** Dissimilarity : mesure la variation des niveaux de gris dans l'imagee*/
//double compute_dissimilarity(Node *n)
//{
//    double dissimilarity = 0;
//    if(0 !=n)
//    {
//        for(int x=0; x<nNg; x++)
//        {
//            for(int y =0; y<nNg; y++)
//            {
//                dissimilarity += std::abs(y-x)*n->m_glcm(y,x);
//            }
//        }
//
//    }
//    return dissimilarity;
//}



//// compute Grey-Level Co-Occurrence Matrix
//void CGraph::computeGLCM(Node *n, vector< Point<int> > &pixels, int nNg)
//{
//    int d = 1; // distance choisie pour parcours de l'image
//    // Parcours des pixels du noeuds
//    for(int i = 0; i<pixels.size(); i++)
//    {
//        Point<int> pCoord = pixels[i];
//        int x_p = pCoord.x;
//        int y_p = pCoord.y;s
//        int z_p = pCoord.z;
//
//        // calcul de la matrice de co-occurrence ds les 13 directions
//        update_pixel(n,x_p,y_p,z_p, x,y-1,z);
//        update_pixel(n,x_p,y_p,z_p, x-1,y-1,z);
//        update_pixel(n,x_p,y_p,z_p, x+1,y,z);
//        update_pixel(n,x_p,y_p,z_p, x+1,y+1,z);
//        update_pixel(n,x_p,y_p,z_p, x,y+1,z);
//        update_pixel(n,x_p,y_p,z_p, x+1,y-1,z+1);
//        update_pixel(n,x_p,y_p,z_p, x+1, y, z+1);
//        update_pixel(n,x_p,y_p,z_p, x+1, y+1, z+1);
//        update_pixel(n,x_p,y_p,z_p, x,y+1,z+1);
//        update_pixel(n,x_p,y_p,z_p, x,y-1,z-1);
//        update_pixel(n,x_p,y_p,z_p, x+1,y-1,z-1);
//        update_pixel(n,x_p,y_p,z_p, x+1,y,z-1);
//        update_pixel(n,x_p,y_p,z_p, x+1,y+1,z-1);
//    }
//
//    // Normalisation de m_glcm
//    double dNormalizationFactor = nNg*(nNg-1);
//    for(int x=0; x<nNg; x++)
//    {
//        for(int y =0; y<nNg; y++)
//        {
//            n->m_glcm(y,x) = n->m_glcm(y,x)/dNormalizationFactor;
//        }
//    }
//
//}

void CGraph::computeMU(float spacingX, float spacingY, float spacingZ)
{
    for(int i=0; i<graph.size()-1; ++i) // parcours de tous les noeuds du graph
	{
		Node *n = graph[i];
		////std::cout<<"n->area"<<n->area<<std::endl;s
		// update des coefficients
        if (n->area > 0.0)
        {
	    n->fVolumemL = n->area*(spacingX*spacingY*spacingZ)*0.001; // facteur 0.001 pour passer de mm3 en mL
	    n->mu200 = n->m200/n->area + 1/12.0 - n->m100*n->m100/(n->area*n->area); //(sumXX-(sumX*sumX/V(C) + V(C)/12)
            n->mu020 = n->m020/n->area + 1/12.0 - n->m010*n->m010/(n->area*n->area);
            n->mu002 = n->m002/n->area + 1/12.0 - n->m001*n->m001/(n->area*n->area);

            n->mu110 = n->m110/n->area - (n->m100 * n->m010)/(n->area*n->area);
            n->mu011 = n->m011/n->area - (n->m010 * n->m001)/(n->area*n->area);
            n->mu101 = n->m101/n->area - (n->m100 * n->m001)/(n->area*n->area);
        }
	}
    
    return;
}

/** Compute elongation attribute on the component-graph
*/
void CGraph::computeElongation()
{

	double sumXX,sumYY,sumZZ,sumXY,sumXZ,sumYZ;
    double A,B,C,D,lambda[3],temp;
	int solutions;
    const double g_Pi = 3.14159265358979323846;
    for(int i=0; i<graph.size()-1; ++i) // parcours de tous les noeuds du graph (sauf noeud racine à la fin)
	{
		Node *tmp=graph[i];
        int emax = 0;
        int emin = 256;

        sumXX=tmp->mu200;
        sumYY=tmp->mu020;
        sumZZ=tmp->mu002;
        sumXY=tmp->mu110;
        sumXZ=tmp->mu101;
        sumYZ=tmp->mu011;

        /** Inertia tensor parameters : coefficients of characteristic polynomial of inertia matrix */
        A = -1.0;
        B = sumXX + sumYY + sumZZ;

        C = -sumXX*sumYY - sumXX*sumZZ - sumYY*sumZZ +
            sumXY*sumXY + sumXZ*sumXZ + sumYZ*sumYZ;

        D = sumXX*sumYY*sumZZ - sumXX*sumYZ*sumYZ
            - sumZZ*sumXY*sumXY + 2*sumXY*sumYZ*sumXZ -
            sumYY*sumXZ*sumXZ;


        // Solve inertia tensor
        if(A == A && B == B && C == C && D == D) //test nan values : check attribution
        {
            SolveCubic(A,B,C,D,&solutions,lambda);
        }

        // Solutions : number of solutions
        if (solutions == 3)
        {
            // sort the solutions : sort Eigen values (lambda) such that (lambda[0] > lambda[1] > lamba[2])
            //-------------------------------
            if (fabs(lambda[1]) > fabs(lambda[0]))
            {
                temp = lambda[1];
                lambda[1] = lambda[0];
                lambda[0] = temp;
            }

            if (fabs(lambda[2]) > fabs(lambda[0]))
            {
                temp = lambda[2];
                lambda[2] = lambda[0];
                lambda[0] = temp;
            }

            else if (fabs(lambda[2]) > fabs(lambda[1]))
            {
				temp = lambda[2];
                lambda[2] = lambda[1];
                lambda[1] = temp;
            }

            //-------------------------------
            // Attribute computation based on eigen values
            //-------------------------------
            // On multiplie par 100 pour passer de float en int
            // elongation along the major axis : computed as the ratio of the first two eigenvalues of the covariance matrix (e1/e2)
            //tmp->elong = (fabs(lambda[0]/lambda[1])-1)*100000;
            tmp->elong = (fabs(lambda[2]/lambda[0]))*1000;

            /// HACK !!!!
            if(tmp->elong>255) tmp->elong=255;
            if(tmp->elong<0) tmp->elong=0;

            if (tmp->elong > emax)
            {
                emax = tmp->elong;
            }
            if (tmp->elong < emin)
            {
                emin = tmp->elong;
            }

            // les noeuds les plus compacts seront mis dans les feuilles du max-tree
            //tmp->elong=255-tmp->elong;
            
            // elongation along  the 2nd axis : computed as e1/sqrt(e2*e3).
            tmp->elong2 = (fabs(lambda[0]/sqrt(lambda[1]*lambda[2]))-1); //*100000

            // flatness along the major axis
            tmp->flatn = (fabs(lambda[1]/lambda[2])-1);

            // flatness along the 2nd axis
            tmp->flatn2 = (fabs(sqrt(lambda[0]*lambda[1])/lambda[2])-1);
        
            tmp->noncompacity=((sumXX+sumYY+sumZZ+tmp->area/4.0f)/pow(float(tmp->area),1.666f)-1.0);
            
            // Add 02/03
            tmp->eccentricity = std::sqrt((lambda[0]-lambda[2])/lambda[0]);
            tmp->sphericity = std::pow((lambda[1]*lambda[2]/sqrt(lambda[0])), 1/3.);
            tmp->elongation = fabs(lambda[0]/lambda[2]);
            tmp->aspectRatio = fabs(lambda[0]/lambda[1]);
            tmp->flatness = fabs(lambda[1]/lambda[2]);
            // computation of sparsness (Wilkinson) : ratio of the expected volume (from the ellipse) and the actual volume
            long double d0 = std::sqrt(20*abs(lambda[0])/tmp->area);
            long double d1 = std::sqrt(20*abs(lambda[1])/tmp->area);
            long double d2 = std::sqrt(20*abs(lambda[2])/tmp->area);
            
            tmp->sparsness = g_Pi*d0*d1*d2/(6*tmp->area);

            if(tmp->elong != tmp->elong) tmp->elong=0;
            if(tmp->elong2 != tmp->elong2) tmp->elong2=0;
            if(tmp->flatn != tmp->flatn) tmp->flatn=0;
            if(tmp->flatn2 != tmp->flatn2) tmp->flatn2=0;
            
            if(tmp->eccentricity != tmp->eccentricity) tmp->eccentricity=0;
            if(tmp->sphericity != tmp->sphericity) tmp->sphericity=0;
            if(tmp->elongation != tmp->elongation) tmp->elongation=0;
            if(tmp->aspectRatio != tmp->aspectRatio) tmp->aspectRatio=0;
            if(tmp->flatness != tmp->flatness) tmp->flatness=0;
            if(tmp->noncompacity != tmp->noncompacity) tmp->noncompacity=0;
            
		}

		else // solutions = 0 or solutions = 1
        {
            // elong and flatn = 1 (instead of 0) so that it doesn't interfer in other attributes computation
            tmp->elong = 0;
            tmp->elong2 = 0;
            tmp->flatn = 0;
            tmp->flatn2 = 0;
            
            tmp->eccentricity = 0;
            tmp->sphericity = 0;
            tmp->elongation = 0;
            tmp->aspectRatio = 0;
            tmp->flatness = 0;
            tmp->noncompacity= 0;
            
        }

	}
    
    return;
}

/** Test if set n is included in m
*/
bool CGraph::isIncluded(Node *n, Node *m, vector<bool> &tmp)
{
    tmp.assign(tmp.size(), false);
    for(int i=0; i<m->pixels.size(); i++)
        tmp[m->pixels[i]]=true;

    for(int i=0; i<n->pixels.size(); i++)
        if(tmp[n->pixels[i]]==false) return false;
    return true;
}

void CGraph::computeLinks(ColorOrdering *order, vector <Node *> &nodes)
{
    vector<bool> tmp(imSource.getBufSize(),false);
    for(int i=0; i<nodes.size(); i++) {
        Node *n=nodes[i];
        for(int j=0; j<nodes.size(); j++) {
            if(j!=i && nodes[j]->area<=nodes[i]->area && order->islessequal(nodes[i]->color,nodes[j]->color)) {
                Node *m=nodes[j];
                if(!order->isequal(m->color,n->color) && isIncluded(m,n,tmp)) {
                    n->addChild(m);
                }
            }
        }
    }
    
    return;
}



CGraph::Node * CGraph::addRoot(vector <Node *> &nodes)
{
    Node *root=new Node(-1,0 ,0);

    for (int i=0; i<nodes.size(); i++) {
        if(nodes[i]->fathers.size()==0) {
            root->addChild(nodes[i]);
        }
    }
    return root;
}


vector <CGraph::Node *> CGraph::minimalElements(ColorOrdering *order, vector <Node *> &nodes, vector <bool> &tmp) {
    vector <Node *> res;

    vector <bool> active(nodes.size(), true);
    for(int j=0; j<nodes.size(); j++) {
        RGB i=nodes[j]->color;
        for(int k=0; k<nodes.size(); k++) {
            if(k!=j) {
                RGB value2=nodes[k]->color;
                if(order->islessequal(i,value2) && isIncluded(nodes[k],nodes[j], tmp) ) {
                    active[k]=false;
                }
            }
        }
    }

    for(int j=0; j<nodes.size(); j++) {
        if(active[j]) res.push_back(nodes[j]);
    }

    return res;
}

/** Compute transitive reduction of graph from its root
* For each node:
*  - keep as childs only the minimal elements of the childs
*/
void CGraph::computeTransitiveReduction(ColorOrdering *order, vector<Node *> &nodes)
{
    vector<bool> tmp(imSource.getBufSize(),false);
    for(int i=0; i<nodes.size(); i++) {
        Node *curNode=nodes[i];
        curNode->childs=minimalElements(order,curNode->childs, tmp);
    }
    
    return;
}

/** Compute the nodes for G and \dot G component-graph
**/
vector<CGraph::Node *> CGraph::computeComponents(ColorOrdering *order, OrderedQueue <RGB> &pq)
{
    int dx=imSource.getSizeX();
    int dy=imSource.getSizeY();

    vector <Node *> nodes;

    Image<bool> active(imSource.getSize());
    Image<bool> inQueue(imSource.getSize());

    inQueue.fill(false);
    std::queue<Point<TCoord > > fifo;

    int index=1;

    while(!pq.empty())
    {
        RGB value=pq.get();

        active.fill(true);

        int ncomposantes=0;

        for(int y=0; y<dy; y++)
            for(int x=0;x<dx; x++) {

                if(active(x,y) && order->islessequal(value,imSource(x,y))) {
                /** Construct a new node and begin a propagation
                **/
                    Point <TCoord> p(x,y);
                    inQueue.fill(false);

                    fifo.push(p);
                    ncomposantes++;
                    inQueue(p)=true;

                    Node *n=new Node(index,value,0);
                    index++;

                    while(!fifo.empty()) {
                        Point<TCoord> curPt=fifo.front();
                        fifo.pop();

                        active(curPt)=false;

                        n->area++;
                        n->pixels.push_back(imSource.getOffset(curPt));
                        for(int i=0; i<connexity.getNbPoints(); i++) {
                            Point<TCoord> q=curPt+connexity.getPoint(i);

                            if(imSource.isPosValid(q)) {

                                if(inQueue(q)==false && order->islessequal(n->color,imSource(q)) ) {
                                    fifo.push(q);
                                    inQueue(q)=true;
                                }
                            }
                        }
                    }
                    nodes.push_back(n);
                }

            }
        
    }

    return nodes;
}



/**
* Compute component-graph G
* -compute all components
* -compute inclusion relations/ construct links
* -compute the transitive reduction
**/
int CGraph::computeGraphFull(ColorOrdering *order, CGraphWatcher *watcher)
{
    Image <RGB> im=this->imSource;
    OrderedQueue <RGB> pq;

    vector<Node *> nodes;

    vector <RGB> colorProcessed;

    /** Put all colors in priority queue
    **/
    for(int r=0; r<256; r++)
        for(int g=0; g<256; g++)
            for(int b=0; b<256; b++)
            {
                RGB value(r,g,b);
                pq.put(order->getPriority(value),value);

            }

    clock_t c1=clock();

    nodes=computeComponents(order,pq);


    clock_t c2=clock();

    c1=clock();
    computeLinks(order,nodes);
    c2=clock();

    for(int i=0; i<nodes.size(); i++)
        graph.push_back(nodes[i]);
    root=addRoot(nodes);

    writeDot("fulllinks.dot");

   computeTransitiveReduction(order, nodes);

    writeDot("fullgraph.dot");
    return 0;
}
