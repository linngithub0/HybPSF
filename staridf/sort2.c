#include "nrutil.h"
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50
void sort2(unsigned long n, float arr[], float brr[])
//Sorts an array arr[1..n] into ascending order using Quicksort, while making the corresponding
//rearrangement of the array brr[1..n] .
{
	unsigned long i,ir=n,j,k,l=1,*istack;
	int jstack=0;
	float a,b,temp;
	istack=lvector(1,NSTACK);
	for (;;) {//Insertion sort when subarray small enough.
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				b=brr[j];
				for (i=j-1;i>=l;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
					brr[i+1]=brr[i];
				}
	            arr[i+1]=a;
	            brr[i+1]=b;
	        }
	        if (!jstack) {
	        	free_lvector(istack,1,NSTACK);
	        	return;
	        }
	        ir=istack[jstack];//Pop stack and begin a new round of parti-tioning.
	        l=istack[jstack-1];
	        jstack -= 2;
	    }
	    else {
	        k=(l+ir) >> 1;//Choose median of left, center and right el-rearrange so that a[l] ≤ a[l+1] ≤ a[ir].\
	        ements as partitioning element a. Also
	        SWAP(arr[k],arr[l+1])
	        
	        SWAP(brr[k],brr[l+1])
	        
	        if (arr[l] > arr[ir]) {
	    	    SWAP(arr[l],arr[ir])
	    	    SWAP(brr[l],brr[ir])
	        }
	    		if (arr[l+1] > arr[ir]) {
	    			SWAP(arr[l+1],arr[ir])
	    			SWAP(brr[l+1],brr[ir])
	    		}
	    		if (arr[l] > arr[l+1]) {
	    			SWAP(arr[l],arr[l+1])
	    			SWAP(brr[l],brr[l+1])
	    		}
	    		i=l+1;
	    		//Initialize pointers for partitioning.
	    		j=ir;
	    		a=arr[l+1];
	    		//Partitioning element.
	    		b=brr[l+1];
	    		for (;;) {
	        		//Beginning of innermost loop.
	        		do i++; while (arr[i] < a);//Scan up to find element > a.
	        		do j--; while (arr[j] > a);//Scan down to find element < a.
	        		if (j < i) break;//Pointers crossed. Partitioning complete.
	        		SWAP(arr[i],arr[j])//Exchange elements of both arrays.
	        		SWAP(brr[i],brr[j])
	    		}//End of innermost loop.
	    		arr[l+1]=arr[j];//Insert partitioning element in both arrays.
	    		arr[j]=a;
	    		brr[l+1]=brr[j];
	    		brr[j]=b;
	    		jstack += 2;//Push pointers to larger subarray on stack, process smaller subarray immediately.
	    		if (jstack > NSTACK) nrerror("NSTACK too small in sort2.");
	    		if (ir-i+1 >= j-l) {
	    			istack[jstack]=ir;
	    			istack[jstack-1]=i;
	    			ir=j-1;
	    		} 
	    		else {
	    			istack[jstack]=j-1;
	    			istack[jstack-1]=l;
	    			l=i;
	    		}
			}
		}
}