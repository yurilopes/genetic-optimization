#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct _BinaryNode{
    int             id; //This should not be used by the programmer
    double          value;
    struct          _BinaryNode* left;
    struct          _BinaryNode* right;
} BinaryNode;

typedef struct _BinaryTree{
    int lastId;
    BinaryNode *firstNode;
} BinaryTree;

BinaryTree *newTree(){
    BinaryTree *result=malloc(sizeof(BinaryTree));
    result->lastId=-1;
    result->firstNode=NULL;
    return result;
}

BinaryNode *newNode(double value, int id){
    BinaryNode *result=malloc(sizeof(BinaryNode));
    memset(result, 0, sizeof(BinaryNode));
    result->value=value;
    result->id=id;
    return result;
}

BinaryNode *GetleftmostLeaf(BinaryNode *Node){
    if(Node==NULL)
        return NULL;
    if(Node->left!=NULL)
        return(GetleftmostLeaf(Node->left));
    else
        return(Node);
}

void RemoveNode(BinaryNode **ParentNode, int id){
    BinaryNode **ChildNode, *tNode, Node;
    int ChildNodes=0;
    if(*ParentNode==NULL) //Nesse caso a árvore está vazia
        return;
    while(1){
        if((*ParentNode)->id==id){ //Encontrei o nó alvo
            ChildNodes=0;

            if((*ParentNode)->left!=NULL){
                ChildNodes++;
                ChildNode=&(*ParentNode)->left;
            }
            if((*ParentNode)->right!=NULL){
                ChildNodes++;
                ChildNode=&(*ParentNode)->right;
            }

            if(ChildNodes==0){ //Sou uma folha, devo virar NULL
                    free(*ParentNode);
                    *ParentNode=NULL;
                    return;
            }

            if(ChildNodes==1){ //Tenho só um filho, guardado em ChildNode
                    tNode=*ParentNode;
                    /*
                     Preciso disso, pois ChildNode representa uma aresta de ParentNode.
                     Como dou free em ParentNode, esse dado se perde, logo *ChildNode se altera
                    */
                    *ParentNode=*ChildNode;
                    free(tNode);
                    return;
            }

            if(ChildNodes==2){ //Tenho dois filhos
                memcpy(&Node, GetleftmostLeaf(*ChildNode), sizeof(BinaryNode));
                RemoveNode(ParentNode, Node.id);
                (*ParentNode)->id=Node.id;
                return;
            }

        }

        if(id>(*ParentNode)->id){ //Se a chave é maior que a minha, então varro para a direita
            if((*ParentNode)->right!=NULL){ //Se há direita, então vou percorrê-la
                ParentNode=&(*ParentNode)->right;
                continue;
            }
            else //Nesse caso a chave não existe na árvore
                return;
        }
        else{ //Se a chave é menor que a minha, então varro para a esquerda
            if((*ParentNode)->left!=NULL){ //Se há esquerda, então vou percorrê-la
                ParentNode=&(*ParentNode)->left;
                continue;
            }
            else //Nesse caso a chave não existe na árvore
                return;
        }
    }
}

BinaryNode *PutNode(BinaryNode **ParentNode, double value, int id){
    if(*ParentNode==NULL){
        *ParentNode=newNode(value, id); //Adiciono a nova folha
        return *ParentNode;
    }
    else{
        if((*ParentNode)->id==id) //A chave já existe na árvore
            return *ParentNode;
        if(id<(*ParentNode)->id)
            PutNode(&(*ParentNode)->left, value, id);
        else
            PutNode(&(*ParentNode)->right, value, id);
    }
}


BinaryNode *putInBinaryTree(BinaryTree *tree, double value){
    if(tree==NULL)
        return;
    tree->lastId++;
    return PutNode(&tree->firstNode, value, tree->lastId);
}

void removeFromBinaryTree(BinaryTree *tree, BinaryNode *node){
    if(tree==NULL)
        return;
    RemoveNode(&tree->firstNode, node->id);
    return;
}
// ----------------------------------------------------------------

void printTreeNodes(BinaryNode* Node){
    if(Node==NULL)
        return;
    printf("ID: %.02d - V: %f\n", Node->id, Node->value); //Pré-ordem
    if(Node->left!=NULL)
        printTreeNodes(Node->left);
    if(Node->right!=NULL)
        printTreeNodes(Node->right);
}

void printTree(BinaryTree* tree){
    if(tree!=NULL)
        printTreeNodes(tree->firstNode);
}
