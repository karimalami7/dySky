#include "common.h"


void listeAttributsPresents(Space subspace, Space d, vector<Space> &result){
/*  retourne la liste des attributs d'un sous espace
    le résultat se présente sous la forme d'un vector
    où chaque cellule correspond au numéro de l'attribut 1, 2, 3, ...,d
    la première valeur [0] est le nombre d'attributs */
    vector<Space> sortie;
    bitset<SPACESIZE> subspace_aux=subspace;
    Space j;
    for (j=0;j<d;j++) if (subspace_aux[j]) sortie.push_back(j+1);
    result.swap(sortie);
}