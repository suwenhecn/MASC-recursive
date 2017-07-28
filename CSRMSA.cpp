#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>  
#include <vector>
#include <string> 
//#include <omp.h>
#define MAX_CHAR 100
using namespace std; 
  
/*#define GET_ARRAY_LEN(char **array,int len){
int len = (sizeof(array)/sizeof(array[0]));
}*/
struct SuffixTreeNode {
    struct SuffixTreeNode *children[MAX_CHAR];
    //pointer to other node via suffix link
    struct SuffixTreeNode *suffixLink;
    /*(start, end) interval specifies the edge, by which the
     node is connected to its parent node. Each edge will
     connect two nodes,  one parent and one child, and
     (start, end) interval of a given edge  will be stored
     in the child node. Lets say there are two nods A and B
     connected by an edge with indices (5, 8) then this
     indices (5, 8) will be stored in node B. */
    int start;
    int *end;
    /*for leaf nodes, it stores the index of suffix for
      the path  from root to leaf*/
    int suffixIndex;
};
  
typedef struct SuffixTreeNode Node;
  
class SuffixTree
{
public:
    char text[40000]; //Input string
    Node *root; //Pointer to root node
  
/*lastNewNode will point to newly created internal node,
  waiting for it's suffix link to be set, which might get
  a new suffix link (other than root) in next extension of
  same phase. lastNewNode will be set to NULL when last
  newly created internal node (if there is any) got it's
  suffix link reset to new internal node created in next
  extension of same phase. */
  Node *lastNewNode ;
  Node *activeNode ;
  
/*activeEdge is represeted as input string character
  index (not the character itself)*/
  int activeEdge;
  int activeLength ;
  
// remainingSuffixCount tells how many suffixes yet to
// be added in tree
 int remainingSuffixCount;
 int leafEnd;
 int *rootEnd;
 int *splitEnd ;
 int size ; //Length of input string
 int size1; //Size of 1st string
 SuffixTree(){
     root = NULL; 
     lastNewNode = NULL;
      activeNode = NULL;
     activeEdge = -1;
     activeLength = 0;
     remainingSuffixCount = 0;
     leafEnd = -1;
     rootEnd = NULL;
    splitEnd = NULL;
    size = -1; 
    size1 = 0; 
 }
 SuffixTree(int s,int s1,string t){
    root = NULL; 
     lastNewNode = NULL;
      activeNode = NULL;
     activeEdge = -1;
     activeLength = 0;
     remainingSuffixCount = 0;
     leafEnd = -1;
     rootEnd = NULL;
    splitEnd = NULL;
    size=s;
    size1=s1;
    strcpy(text,t.c_str());
 }
 Node *newNode(int start, int *end)
{
    Node *node =(Node*) malloc(sizeof(Node));
    int i;
    for (i = 0; i < MAX_CHAR; i++)
          node->children[i] = NULL;
  
    /*For root node, suffixLink will be set to NULL
    For internal nodes, suffixLink will be set to root
    by default in  current extension and may change in
    next extension*/
    node->suffixLink = root;
    node->start = start;
    node->end = end;
  
    /*suffixIndex will be set to -1 by default and
      actual suffix index will be set later for leaves
      at the end of all phases*/
    node->suffixIndex = -1;
    return node;
}
  
int edgeLength(Node *n) {
    if(n == root)
        return 0;
    return *(n->end) - (n->start) + 1;
}
  
int walkDown(Node *currNode)
{
    /*activePoint change for walk down (APCFWD) using
     Skip/Count Trick  (Trick 1). If activeLength is greater
     than current edge length, set next  internal node as
     activeNode and adjust activeEdge and activeLength
     accordingly to represent same activePoint*/
    if (activeLength >= edgeLength(currNode))
    {
        activeEdge += edgeLength(currNode);
        activeLength -= edgeLength(currNode);
        activeNode = currNode;
        return 1;
    }
    return 0;
}
  
void extendSuffixTree(int pos)
{
    /*Extension Rule 1, this takes care of extending all
    leaves created so far in tree*/
    leafEnd = pos;
  
    /*Increment remainingSuffixCount indicating that a
    new suffix added to the list of suffixes yet to be
    added in tree*/
    remainingSuffixCount++;
  
    /*set lastNewNode to NULL while starting a new phase,
     indicating there is no internal node waiting for
     it's suffix link reset in current phase*/
    lastNewNode = NULL;
  
    //Add all suffixes (yet to be added) one by one in tree
    while(remainingSuffixCount > 0) {
        if (activeLength == 0)
            activeEdge = pos; //APCFALZ
  
        // There is no outgoing edge starting with
        // activeEdge from activeNode
        if (activeNode->children[text[activeEdge]] == NULL)
        {
            //Extension Rule 2 (A new leaf edge gets created)
            activeNode->children[text[activeEdge]] =
                                          newNode(pos, &leafEnd);
  
            /*A new leaf edge is created in above line starting
             from  an existng node (the current activeNode), and
             if there is any internal node waiting for it's suffix
             link get reset, point the suffix link from that last
             internal node to current activeNode. Then set lastNewNode
             to NULL indicating no more node waiting for suffix link
             reset.*/
            if (lastNewNode != NULL)
            {
                lastNewNode->suffixLink = activeNode;
                lastNewNode = NULL;
            }
        }
        // There is an outgoing edge starting with activeEdge
        // from activeNode
        else
        {
            // Get the next node at the end of edge starting
            // with activeEdge
            Node *next = activeNode->children[text[activeEdge]];
            if (walkDown(next))//Do walkdown
            {
                //Start from next node (the new activeNode)
                continue;
            }
            /*Extension Rule 3 (current character being processed
              is already on the edge)*/
            if (text[next->start + activeLength] == text[pos])
            {
                //If a newly created node waiting for it's 
                //suffix link to be set, then set suffix link 
                //of that waiting node to curent active node
                if(lastNewNode != NULL && activeNode != root)
                {
                    lastNewNode->suffixLink = activeNode;
                    lastNewNode = NULL;
                }
 
                //APCFER3
                activeLength++;
                /*STOP all further processing in this phase
                and move on to next phase*/
                break;
            }
  
            /*We will be here when activePoint is in middle of
              the edge being traversed and current character
              being processed is not  on the edge (we fall off
              the tree). In this case, we add a new internal node
              and a new leaf edge going out of that new node. This
              is Extension Rule 2, where a new leaf edge and a new
            internal node get created*/
            splitEnd = (int*) malloc(sizeof(int));
            *splitEnd = next->start + activeLength - 1;
  
            //New internal node
            Node *split = newNode(next->start, splitEnd);
            activeNode->children[text[activeEdge]] = split;
  
            //New leaf coming out of new internal node
            split->children[text[pos]] = newNode(pos, &leafEnd);
            next->start += activeLength;
            split->children[text[next->start]] = next;
  
            /*We got a new internal node here. If there is any
              internal node created in last extensions of same
              phase which is still waiting for it's suffix link
              reset, do it now.*/
            if (lastNewNode != NULL)
            {
            /*suffixLink of lastNewNode points to current newly
              created internal node*/
                lastNewNode->suffixLink = split;
            }
  
            /*Make the current newly created internal node waiting
              for it's suffix link reset (which is pointing to root
              at present). If we come across any other internal node
              (existing or newly created) in next extension of same
              phase, when a new leaf edge gets added (i.e. when
              Extension Rule 2 applies is any of the next extension
              of same phase) at that point, suffixLink of this node
              will point to that internal node.*/
            lastNewNode = split;
        }
  
        /* One suffix got added in tree, decrement the count of
          suffixes yet to be added.*/
        remainingSuffixCount--;
        if (activeNode == root && activeLength > 0) //APCFER2C1
        {
            activeLength--;
            activeEdge = pos - remainingSuffixCount + 1;
        }
        else if (activeNode != root) //APCFER2C2
        {
            activeNode = activeNode->suffixLink;
        }
    }
}
  
void print(int i, int j)
{
    int k;
    for (k=i; k<=j && text[k] != '#'; k++)
        printf("%c", text[k]);
    if(k<=j)
        printf("#");
}
  
//Print the suffix tree as well along with setting suffix index
//So tree will be printed in DFS manner
//Each edge along with it's suffix index will be printed
void setSuffixIndexByDFS(Node *n, int labelHeight)
{
    if (n == NULL)  return;
  
    if (n->start != -1) //A non-root node
    {
        //Print the label on edge from parent to current node
        //Uncomment below line to print suffix tree
        //print(n->start, *(n->end));
    }
    int leaf = 1;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != NULL)
        {
            //Uncomment below two lines to print suffix index
         //   if (leaf == 1 && n->start != -1)
           //     printf(" [%d]\n", n->suffixIndex);
  
            //Current node is not a leaf as it has outgoing
            //edges from it.
            leaf = 0;
            setSuffixIndexByDFS(n->children[i], labelHeight +
                                  edgeLength(n->children[i]));
        }
    }
    if (leaf == 1)
    {
        for(i= n->start; i<= *(n->end); i++)
        {
            if(text[i] == '#')
            {
                n->end = (int*) malloc(sizeof(int));
                *(n->end) = i;
            }
        }
        n->suffixIndex = size - labelHeight;
        //Uncomment below line to print suffix index
       // printf(" [%d]\n", n->suffixIndex);
    }
}
  
void freeSuffixTreeByPostOrder(Node *n)
{
    if (n == NULL)
        return;
    int i;
    for (i = 0; i < MAX_CHAR; i++)
    {
        if (n->children[i] != NULL)
        {
            freeSuffixTreeByPostOrder(n->children[i]);
        }
    }
    if (n->suffixIndex == -1)
        free(n->end);
    free(n);
}
  
/*Build the suffix tree and print the edge labels along with
suffixIndex. suffixIndex for leaf edges will be >= 0 and
for non-leaf edges will be -1*/
void buildSuffixTree()
{
    int i;
    rootEnd = (int*) malloc(sizeof(int));
    *rootEnd = - 1;
    
    /*Root is a special node with start and end indices as -1,
    as it has no parent from where an edge comes to root*/
    root = newNode(-1, rootEnd);
    activeNode = root; //First activeNode will be root
    for (i=0; i<size; i++){
        extendSuffixTree(i);
    }
    int labelHeight = 0;
    setSuffixIndexByDFS(root, labelHeight);
}
 
int doTraversal(Node *n, int labelHeight, int* maxHeight, int* substringStartIndex)
{
    if(n == NULL)
    {
        return NULL;
    }
    int i=0;
    int ret = -1;
    if(n->suffixIndex < 0) //If it is internal node
    {
        for (i = 0; i < MAX_CHAR; i++)
        {
            if(n->children[i] != NULL)
            {
                ret = doTraversal(n->children[i], labelHeight + 
                    edgeLength(n->children[i]), 
                    maxHeight, substringStartIndex);
                 
                if(n->suffixIndex == -1)
                    n->suffixIndex = ret;
                else if((n->suffixIndex == -2 && ret == -3) ||
                    (n->suffixIndex == -3 && ret == -2) || 
                    n->suffixIndex == -4)
                {
                    n->suffixIndex = -4;//Mark node as XY
                    //Keep track of deepest node
                    if(*maxHeight < labelHeight)
                    {
                        *maxHeight = labelHeight;
                        *substringStartIndex = *(n->end) - 
                            labelHeight + 1;
                    }
                }
            }
        }
    }
    else if(n->suffixIndex > -1 && n->suffixIndex < size1)//suffix of X
        return -2;//Mark node as X
    else if(n->suffixIndex >= size1)//suffix of Y
        return -3;//Mark node as Y
    return n->suffixIndex;
}
 
char* getLongestCommonSubstring()
{
    int maxHeight = 0;
    int substringStartIndex = 0;
    doTraversal(root, 0, &maxHeight, &substringStartIndex);
    char LCSubstring[20000];
    strcpy(LCSubstring, ""); 
   // printf("LCSubstring=%d\n",strlen(LCSubstring ));
   // printf("LCS size is ==%d\n",strlen(LCSubstring ));
    int k=0;
   // printf("maxHeight=%d\n", maxHeight);
   // printf("maxHeight=%d\n", maxHeight);
    for (k=0; k<maxHeight; k++){
         sprintf(LCSubstring+k,"%c",text[k+substringStartIndex]);
        // printf("k==%d\n", k);
         //  printf("LCS=%d\n",strlen(LCSubstring ));
    }
    return LCSubstring;
}

string LCS(){
    int maxHeight = 0;
    int substringStartIndex = 0;
    doTraversal(root, 0, &maxHeight, &substringStartIndex);
     
    int k;
    char ret[20000];
    //sprintf(ret,"%s",)
    for (k=0; k<maxHeight; k++)
        sprintf(ret+k,"%c",text[k+substringStartIndex]);
//        printf("%c", text[k + substringStartIndex]);
    
    string ret_str=string(ret);
    //cout<<"ret_str is "<<ret_str<<endl;
    if(k==0)
        return "";
    else 
        return ret_str;
}  
    
};


int theshold=20;
int total_build_tree_time=0;
int total_time_lcs=0;


class aligning
{
public:
    string s1;
    string s2;
    const int match=10;
    const int dismatch=-5;
    const int gap=-10;
    int match_score(char a, char b){
        if(a==b) 
         return match;
        else if(a=='-'||b=='-')
         return gap;
        else
          return dismatch;
    }
    void needle(){
        int m=s1.length();
        int n=s2.length();
        vector<vector<int> > score;
        int i,j;
        for (i=0;i<m+1;i++){
            vector<int > temp;
            for (j=0;j<n+1;j++){
                temp.push_back(0);
            }
            score.push_back(temp);
        }

        int s;
        int x;
        int l;
        for(i=0;i<m+1;i++){
            score[i][0]=gap*i;
            //gap*i;
        }
        for(j=0;j<n+1;j++){
           score[0][j] =gap*j;
           //gap*j;
        }
        for(i=1;i<m+1;i++){
            for(j=1;j<n+1;j++){
                s= score[i - 1][j - 1] + match_score(s1[i-1], s2[j-1]);
                x= score[i - 1][j] + gap;
                l= score[i][j - 1] + gap;
                score[i][j] =s;
                if(score[i][j]<x) 
                    score[i][j]=x;
                if(score[i][j]<l)
                    score[i][j]=l;
            }
        }
        string align1, align2;
        char in='-';
        int score_current,score_diagonal,score_up,score_left;
        i=m;
        j=n;
        while(i > 0 && j > 0){
            score_current = score[i][j];
            score_diagonal = score[i-1][j-1];
            score_up = score[i][j-1];
            score_left = score[i-1][j];
            if(score_current == score_diagonal+match_score(s1[i-1], s2[j-1])){
                align1=align1+s1[i-1];
                align2=align2+s2[j-1];
                i--;
                j--;
            }
            else if(score_current == score_left + gap){
                align1=align1+s1[i-1];
                align2=align2+in;
                i--;
            }
            else if(score_current == score_up + gap){
                align1=align1+in;
                align2=align2+s2[j-1];
                j--;
            }
        }
        while(i > 0){
            align1=align1+s1[i-1];
            align2=align2+in;
            i--;
        }
        while(j > 0){
            align1=align1+in;
            align2=align2+s2[j-1];
            j--;
        }
        for(i=0;i<align1.length();i++){
            if(i<s1.length())
            s1[i]=align1[align1.length()-1-i];
            else
             s1=s1+align1[align1.length()-1-i];
        }
        for(i=0;i<align2.length();i++){
            if(i<s2.length())
            s2[i]=align2[align2.length()-1-i];
            else
             s2=s2+align2[align2.length()-1-i];
        }
        return;
    }
};

vector<vector<string> > readfile(vector<vector<string> > &whole){
    ifstream infile;
    string line="1x.fasta";  
    infile.open(line.c_str());
    if(!infile){
        cerr<<"fail to open file\n";
    }
    vector<string> name;
    vector<string> alignment;
    while(getline(infile,line)) {
        if(line[0]=='>'){
            name.push_back(line);
        }
        else{
            alignment.push_back(line);
        }

    }  
    whole.push_back(name);
    whole.push_back(alignment);
    infile.close(); 
    return whole;
}

void Writefile(vector<string> name,vector<string> aligned){
    ofstream in;
    in.open("putout.txt",ios::trunc); 
    int i;
    string center=aligned[aligned.size()-1];
    in<<name[0]<<"\n";
    in<<center<<"\n";
    for(int i=0;i<aligned.size()-1;i++){
       in<<name[i+1]<<"\n";
       in<<aligned[i]<<"\n";
    }
    in.close();//关闭文件
}

double Score(vector<string> aligned,double score){
    int len=aligned.size();
    printf("len==%d\n",len );
    int leng=aligned[0].length();
    printf("length==%d\n",leng);
    int num=len*(len-1)/2;
    int i,j,n;
    score=0.0;
    for(i=0;i<len;i++){
        for(j=i+1;j<len;j++){
            double tmp_score=0.0;
            for(n=0;n<leng;n++){
                if(aligned[i][n]!=aligned[j][n]){
                    if(aligned[i][n]=='-'|| aligned[j][n]=='-'){
                       // printf("i %c j %c\n", aligned[i][n],aligned[j][n]);
                        tmp_score=tmp_score+2;
                    }
                    else{
                       // printf("i %c j %c\n", aligned[i][n],aligned[j][n]);
                        tmp_score++;
                    }
                }
            }
            tmp_score=tmp_score/num;
            score+=tmp_score;
        }
    }
    //score=score/len;
    return score;

}

vector<string>  pair_wise_alignment(string center,string p,vector<string> vtwo){
    vector<string> vvtwo1;
    vector<string> vvtwo2;
    string prefix_1;
    string prefix_2;
    string post_1;
    string post_2;
    string aligned_center;
    string aligned_p;
    string temp=center+"#"+p+"$";
    int size1=center.length()+1;
    int size=center.length()+1+p.length()+1;
    SuffixTree tree=SuffixTree(size,size1,temp);
    tree.buildSuffixTree();
    /*string s1="abcab";
    string s2="abc";
    Tree.building_tree(s1,s2);*/
    char *LCS=tree.getLongestCommonSubstring();
    string lcs=string(LCS);
   // printf("lcs=%d\n",lcs.length() );
    tree.freeSuffixTreeByPostOrder(tree.root);

    int start_pos_1=center.find(lcs);
    int start_pos_2=p.find(lcs);
    prefix_1=center.substr(0,start_pos_1);
    prefix_2=p.substr(0,start_pos_2);
    if(start_pos_1+lcs.length()<=center.length()-1)
       post_1=center.substr(start_pos_1+lcs.length());
    if(start_pos_2+lcs.length()<=p.length()-1)
       post_2=p.substr(start_pos_2+lcs.length());
   
    if (prefix_1.length()>theshold && prefix_2.length()>theshold){
        vvtwo1=pair_wise_alignment(prefix_1,prefix_2,vvtwo1);
        prefix_1=vvtwo1[0];
        prefix_2=vvtwo1[1];
    }
    else{
        aligning align;  
        align.s1=prefix_1;
        align.s2=prefix_2;
        align.needle();
        prefix_1=align.s1;
        prefix_2=align.s2;

    }
    if(post_1.length()>theshold && post_2.length()>theshold){
        vvtwo2=pair_wise_alignment(post_1,post_2,vvtwo2);
        post_1=vvtwo2[0];
        post_2=vvtwo2[1];
    }
    else{
        aligning align;
        align.s1=post_1;
        align.s2=post_2;
        align.needle();
        post_1=align.s1;
        post_2=align.s2;
    }
    if(prefix_1.length()>0)
        aligned_center=prefix_1+lcs;
    else
        aligned_center=lcs;
    if(post_1.length()>0)
        aligned_center=aligned_center+post_1;
    if(prefix_2.length()>0)
        aligned_p=prefix_2+lcs;
    else
        aligned_p=lcs;
    if(post_2.length()>0)
        aligned_p=aligned_p+post_2;

    vtwo.push_back(aligned_center);
   // printf("center=%d\n",aligned_center.length() );
    vtwo.push_back(aligned_p);
    return vtwo;
}

vector<string> merge(vector<string> center_set,vector<string> other_set){
    string center;
    vector<string> merged;
    int iter1;int iter2;
    int flag=0;
    int come=1;
    char *in="-";
    int i=0;
    int j=0;
   // printf("begin while\n");
    while(come==1){
        iter1=i;
        iter2=i;
        for(j=0;j<center_set.size();j++){
           // printf("i==%dj==%d\n",i, j);
          if(center_set[j].length()>i){
            if(center_set[j][i]=='-'){
                flag++;
                break;
            }
          }
        }
        if(flag!=0){
         for(j=0;j<center_set.size();j++){
            if (center_set[j].length()<=i){
                center_set[j]=center_set[j]+'-';
                other_set[j]=other_set[j]+'-';
                //printf("other insert =%c\n",other_set[j][i] );
            }
            else if(center_set[j][i]!='-'){
              //  printf("pre len is %d\n",center_set[j].length() );
            center_set[j].insert(iter1,in);
            other_set[j].insert(iter2,in);
           // printf("%d other insert =%c\n",j,other_set[j][i] );
           // printf("now len is %d\n",center_set[j].length() );
            }
         }
         center=center+in;
        }
     
        else{
            for(j=0;j<center_set.size();j++){
                if (center_set[j].length()<=i){
                    center_set[j]=center_set[j]+'-';
                    other_set[j]=other_set[j]+'-';
                }
            }
            for(j=0;j<center_set.size();j++){
                if(center_set[j].length()>i)
                    break;
            }
            center=center+center_set[j][i];
        }

        flag=0;
        i++;
        come=0;

        for(j=0;j<center_set.size();j++){
            if(center_set[j].length()>i){
                come=1;
                break;
            }
        }
    }
   // printf("while is over\n");
    other_set.push_back(center);
    return other_set;
}

vector<string> mutiple_alignment(vector<string> aligned,vector<string> alignment){
    string center=alignment[0];
    vector<string> center_set;//存放比对结束后的中心序列
    vector<string> other_set;//存放比对后的非中心序列
    double score=0.0;
   for(int i=1;i<alignment.size();i++){ 
       string p=alignment[i];
       string aligned_center;
       string aligned_p;
       vector<string> vtwo;
       long start,end;
       start=clock();
       vtwo=pair_wise_alignment(center,p,vtwo);
       end=clock();
       total_build_tree_time+=(end-start);
       aligned_center=vtwo[0];
       aligned_p=vtwo[1];
       center_set.push_back(aligned_center);
       other_set.push_back(aligned_p);
       //printf("time of %d pair_alignment is %fseconds\n",i,(float)(end-start)/CLOCKS_PER_SEC);
   }

    long start,end;
    start=clock();
    aligned=merge(center_set,other_set);
    end=clock();
    printf("time of merge is %fseconds\n",(float)(end-start)/CLOCKS_PER_SEC);
    score=Score(aligned,score);
    printf("score==%lf\n", score);
    return aligned;
}


int main(){
   vector<vector<string> > whole;
	long start,end;
    start=clock();
    whole=readfile(whole);
    end=clock();
    printf("==========\n") ;
    printf("time of reading file is %fseconds\n",(float)(end-start)/CLOCKS_PER_SEC);
    printf("==========\n") ;

    vector<string> name=whole[0];
    vector<string> alignment=whole[1];
    vector<string> aligned;

    start=clock();
    aligned=mutiple_alignment(aligned,alignment);
    end=clock();
    printf("time of mutiple_alignment is %fseconds\n",(float)(end-start)/CLOCKS_PER_SEC);

    printf("==========\n") ;
    start=clock();
    Writefile(name,aligned);
    end=clock();
    printf("time of Writing file is %fseconds\n",(float)(end-start)/CLOCKS_PER_SEC);
    printf("==========\n") ;
    printf("time of building_tree is %fseconds\n",(float)(total_build_tree_time)/CLOCKS_PER_SEC);
    printf("==========\n") ;
    return 0;
}