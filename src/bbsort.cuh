/*
 * Authored by: Chen, Shifu
 * 
 * Email: chen@gmtk.org / sf.chen@ieee.org
 *
 * The code is distributed under BSD license, you are allowed to use, modify or sell this code, but a statement is required if you used this code any where.
 * 
 */
#ifndef _BBSORT_H_
#define _BBSORT_H_

#define BLOCK_SIZE 512

#define DISORDERLY 0
#define NEARLY_SORTED 1
#define AUTO_EVALUATE 2

void bbSort(float4* dData,int number,int listOrder=AUTO_EVALUATE);

#endif _BBSORT_H_
