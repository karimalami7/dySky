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

long fact(const int n)
{
    long res = 1;
    for(int i = 2; i <= n; i++)
        res*=i;
    return res;
}

void generate_all_orders(Config *cfg, vector<vector<Order>> &all_orders){
    for (int i=0; i<cfg->dyDim_size; i++){
        for (int best_value=0; best_value<cfg->dyDim_val;best_value++){
            for (int worst_value = 0; worst_value <cfg->dyDim_val; ++worst_value){
                if(best_value!=worst_value){
                    all_orders[i].push_back(Order(best_value,worst_value));
                }
            }
        }
    }
}