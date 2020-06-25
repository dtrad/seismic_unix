/* Reading formatted file data with fscanf(). */
#include <stdlib.h>
#include <stdio.h>

int main(void)
{
     float f1, f2, f3, f4, f5;
     FILE *fp;

     if ( (fp = fopen("input.txt", "r")) == NULL)
     {
         fprintf(stderr, "Error opening file.\n");
         exit(1);
     }

     fscanf(fp, "%f %f %f %f %f", &f1, &f2, &f3, &f4, &f5);
     printf("The values are %f, %f, %f, %f, and %f.\n",
             f1, f2, f3, f4, f5);

     fclose(fp);
     return 0;
}
