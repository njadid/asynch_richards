#ifndef MISC_H
#define MISC_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <system.h>
#include <stdlib.h>
#include <stdio.h>

//void CalculateWidth(Link** root,unsigned int N);
void CreateGraph(Link** sys,unsigned int N);
//void CreateGraphRain(Link** sys,unsigned int N,GlobalVars* GlobalVars);
void CalcHortonOrder(Link** sys,unsigned int N,unsigned int* order,unsigned short int* complete);
//int* Partition_System_File(Link** sys,unsigned int N,Link** leaves,unsigned int numleaves,unsigned int** my_sys,unsigned int* my_N,unsigned int* my_max_nodes,TransData* my_data,short int *getting);
void CreateStrComplete(Link** sys,unsigned int N);

#endif //MISC_H