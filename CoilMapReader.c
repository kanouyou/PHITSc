extern struct{
	double xcc, ycc, zcc, bxc, byc, bzc;
} rzforc_;

#define nSegments 5
enum SEGMENTS {cap1, ms1, ts1, coli1, det1};

char* SegNames[nSegments] = { "cap_1.map", "ts_1.map", "ts_2.map", "bs_1.map", "det_1.map"};
//#define nSegments 1
//enum SEGMENTS {cap1};

//char* SegNames[nSegments] = { "cap_1.map"};

struct FileHeader {
	int i_max;
	long nR, nZ, nSheets, magic;
	double x_min, x_max, y_min, y_max, z_min, z_max, innerRadius, outerRadius, length, dR, dZ, error;
} SegmentFileHeaders[nSegments];

float *SegMapX[nSegments];
float *SegMapY[nSegments];
float *SegMapZ[nSegments];
float *SegMapBx[nSegments];
float *SegMapBy[nSegments];
float *SegMapBz[nSegments];
struct FileHeader *SegFileHeaderPtr[nSegments];

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void readFile(const char * filename, struct FileHeader* fh, float* mapX, float* mapY, float* mapZ, float* mapBx, float* mapBy, float* mapBz) {
#ifdef WIN32
	FILE *in = fopen(filename,"rb");
#else
	FILE *in = fopen(filename,"r");
#endif

	int i = 0;
	char str_temp[140];

	// Read in first 4 lines (this is a header that just gives definitions of the variables)
	for (i = 0; i <=3; i++)
	{
		fgets(str_temp, 140, in);
	}

	// Read in each element to the necessary array
	fh->x_min = 999999; fh->y_min = 999999;	fh->z_min = 999999; // something large
	fh->x_max = -99999; fh->y_max = -99999;	fh->z_max = -99999; // something small
	i=0;
	while(feof(in) != 1)
	{
		fscanf(in, "%f", mapX + i);
		fscanf(in, "%f", mapY + i);
		fscanf(in, "%f", mapZ + i);

		// Update min/max coordinates
		//		if ( i > fh->i_max -5) printf("%f %f %f\n", *(mapX+i), *(mapY+i), *(mapZ+i));
		if (*(mapX+i) < fh->x_min) fh->x_min = *(mapX+i);
		if (*(mapY+i) < fh->y_min) fh->y_min = *(mapY+i);
		if (*(mapZ+i) < fh->z_min) fh->z_min = *(mapZ+i);
		if (*(mapX+i) > fh->x_max) fh->x_max = *(mapX+i);
		if (*(mapY+i) > fh->y_max) fh->y_max = *(mapY+i);
		if (*(mapZ+i) > fh->z_max) fh->z_max = *(mapZ+i);

		fscanf(in, "%f", mapBx + i);
		fscanf(in, "%f", mapBy + i);
		fscanf(in, "%f", mapBz + i);
		*(mapBz + i) = *(mapBz + i) * -1; // need to flip Bz to get correct direction
		//printf("Bx = %f \n", *(mapBx + i));
		i++;
	}

	// Save the minimum x, y, z values (the first element in each array)
	/*	fh->x_min = *(mapX + fh->i_max -2); // directions in reverse order - min value is the last element (last line is blank and want the (n-1)th element 
	fh->y_min = *(mapY + 0);
	fh->z_min = *(mapZ + 0);

	// Save the maximum x, y, z values (the final element in each array)
	fh->x_max = *(mapX + 0); // directions in reverse order - max value is the first element
	fh->y_max = *(mapY + fh->i_max -2);
	fh->z_max = *(mapZ + fh->i_max -2);
	*/
	fclose(in);
	in = 0;

	//	printf("x_min=%f, y_min=%f, z_min=%f\n", *(mapX + 0), *(mapY + 0), *(mapZ + fh->i_max -2));
	//	printf("x_min=%f, y_min=%f, z_min=%f\n", fh->x_min, fh->y_min, fh->z_min);
	
	// Rpint out first 10 lines to make sure I;ve done it right
	/*	for(i = 0; i < 10; i++){
		printf("x = %f", *(mapX + i));
		printf(" y = %f", *(mapY + i));
		printf(" z = %f", *(mapZ + i));
		printf(" Bx = %f", *(mapBx + i));
		printf(" By = %f", *(mapBy + i));
		printf(" Bz = %f\n", *(mapBz + i));
	}*/
	
	//printf("file was opened \n");
}

// Read the number of entries in the file (for malloc)
void ReadFileSize(const char * filename, struct FileHeader * fh)
{
#ifdef WIN32
	FILE *in = fopen(filename,"rb");
#else
	FILE *in = fopen(filename,"r");
#endif
	int i = 0;
	char str_temp[200];
	
	// Read in first 4 lines (this is a header that just gives definitions of the variables)
	for (i = 0; i <=3; i++)
	{
		fgets(str_temp, 200, in);
	}
	fh->i_max = 0;
	i = 0;

	while(feof(in) != 1)
	{
		fgets(str_temp, 200, in);
		i++;
	}
	fh->i_max = i;
	//	printf("i_max = %d", i);

	fclose(in);
	in = 0;
}


void findFieldAtPoint(float x, float y, float z, float *fieldPtr, float *mapX, float *mapY, float *mapZ, float *mapBx, float *mapBy, float *mapBz, struct FileHeader* fh){
	
	double localX = x;
	double localY = y;
	double localZ = z;
	
	// Deal with -ve y positions since we only get given magnetic field values for +ve y
	double sign = 1.0;
	if(localY < 0.0) {
		sign = -sign;
		localY = -localY;
	}

	double dX = 10; // 10mm (from file - each point is separated by 10mm in each direction)
	double dY = 10;
	double dZ = 10;
	long nX = (fh->x_max - fh->x_min)/dX +1; // +1 so that we include the minimum value
	long nY = (fh->y_max - fh->y_min)/dY +1;
	long nZ = (fh->z_max - fh->z_min)/dZ +1;

	// Find the right elements for the given position
	// Find integer part of steps (want to round down) before multiplying by nX, nY etc (otherwise fractional parts might multiple up to integer parts and get the wrong element)
	int x_steps = (abs(localX - *(mapX+0))/dX);
	int y_steps = (abs(localY - *(mapY+0))/dY);
	int z_steps = (abs(localZ - *(mapZ+0))/dZ);
	int i = x_steps + y_steps*nX + z_steps*nY*nX;
	//	printf("x_min=%f, y_min=%f, z_min=%f\n", fh->x_min, fh->y_min, fh->z_min);
	//	printf("x_steps=%d, y_steps=%d, z_steps=%d\n", x_steps, y_steps, z_steps);
/*	for (i = 0; i < fh->i_max; i++)
	{
		//printf("localX=%f ", localX);
		if ( *(mapX + i) <= localX && (*(mapX+i)+dX) > localX && *(mapY + i) <= localY && (*(mapY+i)+dY) > localY && *(mapZ + i) <= localZ && (*(mapZ+i)+dZ) > localZ)
			break;
	}*/
		//printf("i_max = %d, i = %d\n", fh->i_max, i);
		//printf("localX=%f localY=%f localZ=%f dX=%lf dY=%lf dZ=%lf, nX=%d, nY=%d, nZ=%d\n", localX, localY, localZ, dX, dY, dZ, nX, nY, nZ);
	
	// if we have found the correct position (i.e. the position is in this segment)
	if(i < fh->i_max -1){
		float fx = 1.0 - abs(localX - *(mapX+i)) / dX;
		float fy = 1.0 - abs(localY - *(mapY+i)) / dY;
		float fz = 1.0 - abs(localZ - *(mapZ+i)) / dZ;
				//printf("localX=%f localY=%f localZ=%f\nmapX=%f mapY=%f mapZ=%f\ndX=%lf dY=%lf dZ=%lf\nnX=%d, nY=%d, nZ=%d\nfx=%f fy=%f fz=%f \n", localX, localY, localZ, *(mapX+i), *(mapY+i), *(mapZ+i), dX, dY, dZ, nX, nY, nZ, fx, fy, fz);
				//printf("(Before): Bx = %f, By = %f, Bz = %f\n", *(mapBx + i), *(mapBy + i), *(mapBz + i));

		float Bx = *(mapBx + i)*fx*fy*fz; // point (i, j, k)
		float By = *(mapBy + i)*fx*fy*fz; // point (i, j, k)
		float Bz = *(mapBz + i)*fx*fy*fz; // point (i, j, k)

		// point (i, j, k+1)
		if ( (i + nX*nY) < fh->i_max -1) {
		  Bx += *(mapBx + i + nX*nY)*fx*fy*(1.0-fz);
		  By += *(mapBy + i + nX*nY)*fx*fy*(1.0-fz);
		  Bz += *(mapBz + i + nX*nY)*fx*fy*(1.0-fz);
		}

		// point (i, j+1, k)
		if ( (i + nX) < fh->i_max -1) {
		  Bx += *(mapBx + i + nX)*fx*(1.0-fy)*fz;
		  By += *(mapBy + i + nX)*fx*(1.0-fy)*fz;
		  Bz += *(mapBz + i + nX)*fx*(1.0-fy)*fz;
		}

		// point (i+1, j, k)
		if ( (i + 1) < fh->i_max -1) {
		  Bx += *(mapBx + i + 1)*(1.0-fx)*fy*fz;
		  By += *(mapBy + i + 1)*(1.0-fx)*fy*fz;
		  Bz += *(mapBz + i + 1)*(1.0-fx)*fy*fz;
		}

		// point (i+1, j, k+1)
		if ( (i + 1 + nX*nY) < fh->i_max -1) {
		  Bx += *(mapBx + i + 1 + nX*nY)*(1.0-fx)*fy*(1.0-fz);
		  By += *(mapBy + i + 1 + nX*nY)*(1.0-fx)*fy*(1.0-fz);
		  Bz += *(mapBz + i + 1 + nX*nY)*(1.0-fx)*fy*(1.0-fz);
		}

		// point (i, j+1, k+1)
		if ( (i + nX + nX*nY) < fh->i_max -1) {
		  Bx += *(mapBx + i + nX + nX*nY)*fx*(1.0-fy)*(1.0-fz);
		  By += *(mapBy + i + nX + nX*nY)*fx*(1.0-fy)*(1.0-fz);
		  Bz += *(mapBz + i + nX + nX*nY)*fx*(1.0-fy)*(1.0-fz);
		}

		// point (i+1, j+1, k)
		if ( (i + 1 + nX) < fh->i_max -1) {
		  Bx += *(mapBx + i + 1 + nX)*(1.0-fx)*(1.0-fy)*fz;
		  By += *(mapBy + i + 1 + nX)*(1.0-fx)*(1.0-fy)*fz;
		  Bz += *(mapBz + i + 1 + nX)*(1.0-fx)*(1.0-fy)*fz;
		}

		// point (i+1, j+1, k+1)
		if ( (i + 1 + nX + nY*nX) < fh->i_max -1) {
		  Bx += *(mapBx + i + 1 + nX + nY*nX)*(1.0-fx)*(1.0-fy)*(1.0-fz);
		  By += *(mapBy + i + 1 + nX + nY*nX)*(1.0-fx)*(1.0-fy)*(1.0-fz);
		  Bz += *(mapBz + i + 1 + nX + nY*nX)*(1.0-fx)*(1.0-fy)*(1.0-fz);
		}
				
				//printf("(After): Bx = %f, By = %f, Bz = %f\n", Bx, By, Bz);

		// Assign the magnetic field values at this point to the array
		*(fieldPtr) = Bx; // there will be a sign change
		*(fieldPtr + 1) = sign * By;
		*(fieldPtr + 2) = Bz;
		//		printf("(After): Bx=%lf By=%lf Bz= %lf \n", *(fieldPtr), *(fieldPtr+1), *(fieldPtr+2));
	}
	else { printf("Position not in this segment\n"); }
		//printf("\n");
}

void readfieldfrommap_() {
	float field[3] = {0., 0., 0.};
	float *fieldPtr = field;
	//int printFieldToFile = 1; // if = 0 don't print to .txt file, if = 1 then print field to .txt file

	// Get the directory where the fieldmap is stored from the environment variable
	char* fieldmap_location = "/Users/YeYang/Documents/COMET/information/fieldMaps/140624/";
	int iSeg = 0;
	for (iSeg = 0; iSeg < nSegments; ++iSeg) {
	  char seg_map_location[120];
	  strcpy(seg_map_location, fieldmap_location);
	  strcat(seg_map_location, SegNames[iSeg]);
	  // Read in the file if it hasn't been read before
	  if (SegMapX[iSeg] == 0 || SegMapY[iSeg] == 0 || SegMapZ[iSeg] == 0 || SegMapBx[iSeg] == 0 || SegMapBy[iSeg] == 0 || SegMapBz[iSeg] == 0)
	    {

	      SegFileHeaderPtr[iSeg] = &SegmentFileHeaders[iSeg];

	      // Find out how many lines of data there are
	      ReadFileSize(seg_map_location, SegFileHeaderPtr[iSeg]);

	      // Allocate enough memory for all the files
	      SegMapX[iSeg] = malloc(SegFileHeaderPtr[iSeg]->i_max * sizeof (float));
	      SegMapY[iSeg] = malloc(SegFileHeaderPtr[iSeg]->i_max * sizeof (float));
	      SegMapZ[iSeg] = malloc(SegFileHeaderPtr[iSeg]->i_max * sizeof (float));
	      SegMapBx[iSeg] = malloc(SegFileHeaderPtr[iSeg]->i_max * sizeof (float));
	      SegMapBy[iSeg] = malloc(SegFileHeaderPtr[iSeg]->i_max * sizeof (float));
	      SegMapBz[iSeg] = malloc(SegFileHeaderPtr[iSeg]->i_max * sizeof (float));
	      
	      // Read the data in
	      readFile(seg_map_location, SegFileHeaderPtr[iSeg], SegMapX[iSeg], SegMapY[iSeg], SegMapZ[iSeg], SegMapBx[iSeg], SegMapBy[iSeg], SegMapBz[iSeg]);
	    }
	}

		//printf("x=%f, y=%f, z=%f\n", rzforc_.xcc, rzforc_.ycc, rzforc_.zcc);

	for (iSeg = 0; iSeg < nSegments; ++iSeg) {
		
	  // Only try and find field at that point if it is in the segment
	  if (rzforc_.xcc>=SegFileHeaderPtr[iSeg]->x_min && rzforc_.xcc<=SegFileHeaderPtr[iSeg]->x_max && rzforc_.ycc>=(-SegFileHeaderPtr[iSeg]->y_max) && rzforc_.ycc<=SegFileHeaderPtr[iSeg]->y_max && rzforc_.zcc>=SegFileHeaderPtr[iSeg]->z_min && rzforc_.zcc<=SegFileHeaderPtr[iSeg]->z_max)
	    {
	      findFieldAtPoint(rzforc_.xcc, rzforc_.ycc, rzforc_.zcc, fieldPtr, SegMapX[iSeg], SegMapY[iSeg], SegMapZ[iSeg], SegMapBx[iSeg], SegMapBy[iSeg], SegMapBz[iSeg], SegFileHeaderPtr[iSeg]);
	    
	      //printf("readfieldfrommap\n pos = %f %f %f\n field = %f %f %f\n", rzforc_.xcc, rzforc_.ycc, rzforc_.zcc, field[0], field[1], field[2]);
	    }

	  if (iSeg == ts1) {
	    //	    printf("Before: %f\n", field[1]);
	    field[1] += 0.0561;
	    //	    printf("After: %f\n", field[1]);
	  }
	}
	
	//printf("R=%lf Z=%lf Br=%lf Bz= %lf \n", rzforc_.rc, rzforc_.zc, field[0], field[1]);	
	rzforc_.bxc = field[0];
	rzforc_.byc = field[1];
	rzforc_.bzc = field[2];

}
