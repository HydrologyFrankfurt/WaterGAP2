#include <stdio.h>
#define M 7
#define NSTACK 50

void index(int nn, short int arr[], unsigned short int indx[])
	 /* indexes an array 'arr[0..nn-1]', i.e. outputs the array  */
	 /* 'indx[0..nn-1]' such that 'arr[indx[j]]' is in ascending */
	 /* order for 'j=0,..,N'.                                    */
	 /* The input quantities nn and arr are not changed.         */
	 /* Numerical Recipes in C, 2nd edition, chapter 8.4         */
{
	int i, indxt, ir = nn - 1, itemp, j, k, l = 0;
	int jstack = 0;
	float a;
	int istack[NSTACK];


	for (j = 0; j <= nn - 1; j++)
		indx[j] = j;
	for (;;) {
		if (ir - l < M) {
			for (j = l + 1; j <= ir; j++) {
				indxt = indx[j];
				a = arr[indxt];
				for (i = j - 1; i >= 1; i--) {
					if (arr[indx[i]] <= a)
						break;
					indx[i + 1] = indx[i];
				}
				indx[i + 1] = indxt;
			}
			if (jstack == 0)
				break;
			ir = istack[jstack--];
			l = istack[jstack--];
		} else {
			k = (l + ir) >> 1;

			itemp = indx[k];
			indx[k] = indx[l + 1];
			indx[l + 1] = itemp;

			if (arr[indx[l + 1]] > arr[indx[ir]]) {

				itemp = indx[l + 1];
				indx[l + 1] = indx[ir];
				indx[ir] = itemp;
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				itemp = indx[l];
				indx[l] = indx[ir];
				indx[ir] = itemp;
			}
			if (arr[indx[l + 1]] > arr[indx[l]]) {
				itemp = indx[l + 1];
				indx[l + 1] = indx[l];
				indx[l] = itemp;
			}
			i = l + 1;
			j = ir;
			indxt = indx[l];
			a = arr[indxt];
			for (;;) {
				do
					i++;
				while (arr[indx[i]] < a);
				do
					j--;
				while (arr[indx[j]] > a);
				if (j < i)
					break;
				itemp = indx[i];
				indx[i] = indx[j];
				indx[j] = itemp;
			}
			indx[l] = indx[j];
			indx[j] = indxt;
			jstack += 2;
			if (jstack > NSTACK)
				printf("NSTACK too small in indexx.");
			if (ir - i + 1 >= j - l) {
				istack[jstack] = ir;
				istack[jstack - 1] = i;
				ir = j - 1;

			} else {
				istack[jstack] = j - 1;
				istack[jstack - 1] = l;
				l = i;
			}
		}
	}
}
