//Copyright (C) 2014, Benoît Naegel <b.naegel@unistra.fr> and modified by Eloïse Grossiord <eloise.grossiord@gmail.com>
//This program is free software: you can use, modify and/or
//redistribute it under the terms of the GNU General Public
//License as published by the Free Software Foundation, either
//version 3 of the License, or (at your option) any later
//version. You should have received a copy of this license along
//this program. If not, see <http://www.gnu.org/licenses/>.

#include "shaping.h"
#include "cgraph.h"

Shaping::Shaping()
{
}

/** Compute the max-tree of the underlying weighted graph */
void Shaping::computeShaping()
{
    totalNodes=0;
    init();
    computeTree();
    delete[] hq;
}



void Shaping::link_node(ShapingNode *tree, ShapingNode *child)
{
    child->father=tree;
    tree->childs.push_back(child);
}

Shaping::ShapingNode *Shaping::new_node(int weight, int index)
{

    ShapingNode *res=new ShapingNode(index,weight,0,totalNodes);
    totalNodes++;

    return res;
}

// index=index of the node in the graph
void Shaping::update_attributes(ShapingNode *n, int index)
{
   //special case for the fictitious root ...
   //if(index!=cgraph->graph.size()-1)
        n->nodes.push_back(index);
   n->area++;
}

int Shaping::flood(int h)
{
    int m;

    while(!hq[h].empty())
    {

        //std::cout << "Processing level : " << h << "\n";
        int p=hq[h].front();
        Node *curNode=cgraph->graph[p];
        hq[h].pop();

        assert(p>=0 && p<STATUS.getBufSize());
        STATUS(p)=number_nodes[h];

        if(index[h][STATUS(p)]==0)
        {
            index[h][STATUS(p)]=this->new_node(indexToH(h),STATUS(p));
        }

        update_attributes(index[h][STATUS(p)],p);

        // visit all neighbors i.e. fathers...

        for(int i=0; i<curNode->fathers.size(); i++) {
            {
                Node *n=curNode->fathers[i];

                int q=n->index;

                if(q<0 || q>=STATUS.getBufSize()) {
                    std::cout << " WARNING !! q = " << q << " GRAPH SIZE=" << cgraph->graph.size() <<" STATUS SIZE : " << STATUS.getBufSize()<<"\n";
                }
                assert(q>=0 && q<STATUS.getBufSize());
                if(STATUS(q)==ACTIVE)
                {

                    hq[hToIndex(n->elong)].push(q);
                    STATUS(q)=NOT_ACTIVE;

                    node_at_level[hToIndex(n->elong)]=true;

                    if(n->elong>curNode->elong)
                    {

                        m=hToIndex(n->elong);

                        do 	{

                            m=this->flood(m);

                        } while(m!=h);
                    }
                }
            }
        }

        // ... and childs
        for(int i=0; i<curNode->childs.size(); i++) {
            Node *n=curNode->childs[i];

            int q=n->index;

            if(q<0 || q>=STATUS.getBufSize()) {
                std::cout << " WARNING !! q = " << q << " GRAPH SIZE=" << cgraph->graph.size() <<" STATUS SIZE : " << STATUS.getBufSize()<<"\n";
            }
            assert(q>=0 && q<STATUS.getBufSize());
            if(STATUS(q)==ACTIVE)
            {

                hq[hToIndex(n->elong)].push(q);
                STATUS(q)=NOT_ACTIVE;

                node_at_level[hToIndex(n->elong)]=true;

                if(n->elong>curNode->elong)
                {

                    m=hToIndex(n->elong);

                    do 	{

                        m=this->flood(m);

                    } while(m!=h);
                }
            }
        }
    }
    //End of recursion: we have reached a regional maximum
    number_nodes[h]=number_nodes[h]+1;

    m=h-1;
    while(m>=hToIndex(hMin) && node_at_level[m]==false) m--;

    if(m>=hToIndex(hMin) )
    {
        int i=number_nodes[h]-1;
        int j=number_nodes[m];
        if(index[m][j]==0)
        {
            index[m][j]=new_node(indexToH(m),j);
        }

        this->link_node(index[m][j],index[h][i]);
    }
    else
    {
        //The father of root is itself
        index[hToIndex(hMin)][0]->father=index[hToIndex(hMin)][0];
    }
    node_at_level[h]=false;
    return m;
}


void Shaping::computeTree()
{
    //Put the first pixel with value hMin in the queue
    for(int i=0; i<cgraph->graph.size(); i++) {
        Node *n=cgraph->graph[i];
            if(n->elong==hMin && STATUS(i)==ACTIVE)
                {
                hq[hToIndex(hMin)].push(i);
                break;
                }
    }

    node_at_level[hToIndex(hMin)]=true;

    this->flood(hToIndex(hMin));

    root=index[hToIndex(hMin)][0];

    //for(int i=0; i<STATUS.getBufSize();i++) assert(STATUS(i)==NOT_ACTIVE);

}


void Shaping::init()
{
    // work with this->graph in input

    STATUS.setSize(cgraph->graph.size(),1,1 );
    STATUS.fill(ACTIVE);

    int min,max;
    min=std::numeric_limits<int>::max();
    max=std::numeric_limits<int>::min();
    // Compute min and max weight

    for(int i=0; i<cgraph->graph.size() ; i++) 
    {
        Node *n=cgraph->graph[i];
  	//std::cout << "n->elong " << n->elong << std::endl;
        if(n->elong > max) max = n->elong;
        if(n->elong < min) min = n->elong;
    }

    this->hMin=min;
    this->hMax=max;
    this->numberOfLevels=hMax-hMin+1;

    index.resize(numberOfLevels);

    hq = new std::queue<int> [numberOfLevels];

    //we take a (max-min+1) * (number of grey-levels at level h)
    // so we compute histogram

    int *histo;
    histo =new int[numberOfLevels];
    for(int i=0; i<numberOfLevels; i++)
        histo[i]=0;


    for(int i=0; i<cgraph->graph.size() ; i++)
        {
            Node *n = cgraph->graph[i];
            //int v = n->elong-hMin);
	    // transformer la difference de double vers int
	    int v = (int) ((n->elong-hMin) + 0.5);
            histo[v]++;
        }

    for(int i=0; i<numberOfLevels; i++)
        {
            int sizeOfLevelH=histo[i];

            index[i].resize(sizeOfLevelH);

            for(int j=0; j<sizeOfLevelH; j++)
                 index[i][j]=0;
        }

    this->number_nodes.resize(numberOfLevels);
    this->node_at_level.resize(numberOfLevels);

    for(int i=0; i<numberOfLevels; i++)
        {
        this->number_nodes[i]=0;
        this->node_at_level[i]=false;
        }

    delete[] histo;
}


int Shaping::computeAreaHelper(ShapingNode *tree)
{
    if(tree!=0)
        {

        for(int i=0; i<tree->childs.size(); i++)
            {
            tree->area+=computeAreaHelper(tree->childs[i]);
            }
        return tree->area;
        }
    // error
    else return -1;
}

int Shaping::computeArea()
{
    root->area=computeAreaHelper(root);
}

int Shaping::computeContrastHelper(ShapingNode *node)
{
    if(node!=0)
        {
        int current_level=node->contrast;
        int current_max=0;
        int current_contrast=0;

        for(int i=0; i< node->childs.size(); i++)
            {
            current_contrast=(node->childs[i]->contrast-current_level)+computeContrastHelper(node->childs[i]);
            if(current_contrast>current_max)
                current_max=current_contrast;
            }
        node->contrast=current_max;
        return node->contrast;
        }
    // error
    else return -1;
}


int Shaping::computeContrast()
{
    root->contrast=computeContrastHelper(root);
}


void Shaping::areaFiltering(int areaMin, int areaMax)
{
    if(root!=0)
    {
        std::queue <ShapingNode *> fifo;
        fifo.push(root);

        while(!fifo.empty() )
        {
            ShapingNode *curNode=fifo.front();
            fifo.pop();
            if(curNode->area < areaMin || curNode->area > areaMax)
            {
                curNode->active=false;
            }

            for(int i=0; i<curNode->childs.size(); i++)
            {
                fifo.push(curNode->childs[i]);
            }
        }
    }
}


void Shaping::areaFiltering(int areaMin)
{
    int nbDeletedNode = 0;
    if(root!=0)
    {
        std::queue <ShapingNode *> fifo;
        fifo.push(root);

        while(!fifo.empty() )
        {
            ShapingNode *curNode=fifo.front();
            fifo.pop();
            if(curNode->area < areaMin ) // on elimine noeuds dont area est inf au seuil cad contenant un nombre de noeuds inf au seuil
            {
                nbDeletedNode++;
                curNode->active=false;
                std::cout << "Area FILTERING : area=" <<curNode->area << "active = false\n";
            }

            for(int i=0; i<curNode->childs.size(); i++)
            {
                fifo.push(curNode->childs[i]);
            }
        }
    }
    return;
}


// height filtering
void Shaping::contrastFiltering(int tMin, int tMax)
{
    if(root!=0)
    {
        std::queue <ShapingNode *> fifo;
        fifo.push(root);

        while(!fifo.empty() )
        {
            ShapingNode *curNode=fifo.front();
            fifo.pop();
            if(curNode->contrast < tMin || curNode->contrast > tMax)
            {
                curNode->active=false;
                std::cout << "Contrast FILTERING : contrast=" <<curNode->contrast << "active = false\n";
            }

            for(int i=0; i<curNode->childs.size(); i++)
            {
                fifo.push(curNode->childs[i]);
            }
        }
    }
    return;
}


// height filtering
void Shaping::contrastFiltering(int tMin)
{
    if(root!=0)
    {
        std::queue <ShapingNode *> fifo;
        fifo.push(root);

        while(!fifo.empty() )
        {
            ShapingNode *curNode=fifo.front();
            fifo.pop();
            if(curNode->contrast < tMin )
            {
                curNode->active=false;
            }

            for(int i=0; i<curNode->childs.size(); i++)
            {
                fifo.push(curNode->childs[i]);
            }
        }
    }
    return;
}

// height max filtering (non-increasing)
void Shaping::contrastFilteringMax(int tMax)
{
    std::cout<<" CONTRAST FILTERING SUR ARBRE DU SHAPING" << std::endl;
    int nbDeletedNode = 0;
    if(root!=0)
    {
        std::queue <ShapingNode *> fifo;
        fifo.push(root);

        while(!fifo.empty() )
        {
            ShapingNode *curNode=fifo.front();
            fifo.pop();
            if(curNode->contrast > tMax ) // on elimine noeuds dont la hauteur est superieure au seuil. Hauteur de 0 aux feuilles puis croissante qd on va vers racine
            {

                curNode->active=false;
                nbDeletedNode++;
            }

            for(int i=0; i<curNode->childs.size(); i++)
            {
                fifo.push(curNode->childs[i]);
            }

        }
    }
    return;
}

// elong filtering (non-increasing)
void Shaping::elongFiltering(int tMin)
{
	std::cout << "ENTER ELONGFILTERING WITH SHAPING" << std::endl;
    if(root!=0)
    {
        std::queue <ShapingNode *> fifo;
        fifo.push(root);

        while(!fifo.empty())
        {
            ShapingNode *curNode=fifo.front();
            fifo.pop();
            if(curNode->weight < tMin ) // on elimine noeud dont elong est inf au seuil
            {
                curNode->active=false; // desactive
            }

            for(int i=0; i<curNode->childs.size(); i++)
            {
                fifo.push(curNode->childs[i]);
            }

        }
    }
}


// shaping to graph
// modify cgraph in place (activate set of nodes activated in the shaping)
void Shaping::constructGraph()
{
    std::cout << "RECONSTRUCT GRAPH\n";
    // FIX
    // !! le dernier noeud est la racine fictive du graphe, il ne faut donc pas le parcourir ni le changer !!!
    //
    // deactivate all nodes

    for(int i=0; i< cgraph->graph.size()-1; i++) {
        cgraph->graph[i]->active=false;
    }
    assert(cgraph->graph[cgraph->graph.size()-1]->active==true);
    assert(root!=0);
    // tree breadth first
    int nb_active_nodes = 0;
    int number_active_nodes = 0;
    if(root!=0)
    {
        std::cout << "root != 0 " << std::endl;
        std::queue <ShapingNode *> fifo;
        fifo.push(root);

        while(!fifo.empty() )
        {
            ShapingNode *curNode=fifo.front();
            fifo.pop();

            if(curNode->active==true)
            {
                //std::cout<< "enter dans le if de constructGraph" << std::endl;
                nb_active_nodes +=1;
                for(int i=0; i<curNode->nodes.size(); i++)
                {
                    cgraph->graph[curNode->nodes[i]]->active=cgraph->graph[curNode->nodes[i]]->active2;

                }
            }
            for(int i=0; i<curNode->childs.size(); i++)
            {
                fifo.push(curNode->childs[i]);
            }
        }
    }
    assert(cgraph->graph[cgraph->graph.size()-1]->active==true);

    for(int i=0; i< cgraph->graph.size(); i++) {

        if(cgraph->graph[i]->active==true)
        {number_active_nodes++;}
    }

    std::cout<< "the number of nodes in the initial graph is " << cgraph->graph.size() << std::endl;
    std::cout << "pass exterieur : the number of active nodes is " << number_active_nodes << std::endl;

}

int Shaping::writeDot(const char *filename)
{
    if(root!=0)
        {
        std::ofstream file(filename,std::ios_base::trunc  | std::ios_base::binary);
        if(!file)
            {
            std::cerr << "File I/O error\n";
            return 0;
            }

        file << "digraph G {\n";

        std::queue <ShapingNode *> fifo;
        fifo.push(root);
        while(!fifo.empty() )
            {
            ShapingNode *tmp=fifo.front();
            fifo.pop();

            if(tmp->active==true)
                {
                // write father->son relation if the node is not the root
                if(tmp->father!=tmp)
                {
                file << "\t" << " \"" << tmp->father->uniqueId << ",h=" << tmp->father->weight << "," <<  " a= " << tmp->father->area << "," << "c= " << tmp->father->contrast << "\" "
                   << "->" << " \"" << tmp->uniqueId << ",h=" << tmp->weight << ","   << " a= " << tmp->area << "," << "c= " << tmp->contrast <<
                    "\" " <<
                    ";\n";
                }
                }
                // push the childs
                for(int i=0; i<tmp->childs.size(); i++)
                    {
                    fifo.push(tmp->childs[i]);
                    }


            }

        file << "}\n";

        file.close();
        return 1;
        }
    else
        return 0;
}
