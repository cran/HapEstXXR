/*******************************************************************************
i_quiksort

to sort integer list

February 12, 2010

*******************************************************************************/
#include <stdio.h>
#include "R.h"

void i_swap(int l[], int i, int j)
{  int dummy;
   dummy=l[j];
   l[j]=l[i];
   l[i]=dummy;
 }

int i_partions(int l[],int low,int high)
{
    int prvotkey=l[low];
    while (low<high)
    {
        while (low<high && l[high]>=prvotkey)
            --high;
        i_swap(l,high,low);
        while (low<high && l[low]<=prvotkey)
            ++low;
        i_swap(l,high,low);
    }

    return low;
}

void i_qsort(int l[],int low,int high)
{
    int prvotloc;
    if(low<high)
    {
        prvotloc=i_partions(l,low,high);
        i_qsort(l,low,prvotloc);
        i_qsort(l,prvotloc+1,high);
    }
}

void i_quicksort(int l[],int n)
{
    i_qsort(l,0,n);
}
