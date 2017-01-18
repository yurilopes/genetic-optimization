using namespace std;

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <list>

#include "genetic_operations.h"
#include "population.h"

//#include "binarytree.h"
/*
    BinaryTree *tree;
    BinaryNode *node;

    tree=newTree();

    putInBinaryTree(tree, 0);
    putInBinaryTree(tree, 1);
    node=putInBinaryTree(tree, 2);
    putInBinaryTree(tree, 3);
    putInBinaryTree(tree, 2);

    removeFromBinaryTree(tree, node);

    printTree(tree);
*/

/*
    Tenho uma população definida por Linked Lists de variáveis do problema.
*/

int main(){
    /*
    A population is list of a lisft of ints
    Each individual is a list of ints
    Each problem variable is a int
    */
    Population population(2, 2);

    population.printPopulation();

    //initializePopulation(population, 2, 2);

    //printPopulation(population);

    //crossOver(population.front(), population.back());


    return 0;
}
