#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "params.h"


/* Calculate distance and azimuth/back-azimuth between two (lon,lat) coordinates 
 * using delaz.f subroutine
 * Shoeball and Ellsworth locations
 */

// function declaration
void delaz_(float *lat1, float *lon1, float *lat2, float *lon2, float *dist, float *az, float *baz);
int getcols( const char * const line, const char * const delim, char ***out_storage);
//void strip(char *s);
void assign_cols_flatfile(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec);
void assign_cols_catalog(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec, float *evMagSource, char *magString);
char * replace(char const * const original, char const * const pattern, char const * const replacement );
int compute_epochTime(int yearIn, int monthIn, int dayIn, int hourIn, int minIn, int secIn);


/*--------------------------------------------------------------------------*/
void assign_cols_flatfile(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec)
/*--------------------------------------------------------------------------*/
{
//
  *evMag=atof(columns[6]);
  *evYear=atoi(columns[0]);
  *evMon=atoi(columns[1]);
  *evDay=atoi(columns[2]);
  *evHour=atoi(columns[3]);
  *evMin=atoi(columns[4]);
  *evSec=atof(columns[5]);
//  fprintf(stderr,"assign_cols_flatfile, year/month/day/hour/min/sec/Mag: %d %d %d %d %d %.2f %f\n", *evYear, *evMon, *evDay, *evHour, *evMin, *evSec, *evMag);

}

/*--------------------------------------------------------------------------*/
void assign_cols_catalog(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec, float *evMagSource, char *magString)
/*--------------------------------------------------------------------------*/
{
//
  *evMag=atof(columns[0]);
  *evYear=atoi(columns[4]);
  *evMon=atoi(columns[5]);
  *evDay=atoi(columns[6]);
  *evHour=atoi(columns[7]);
  *evMin=atoi(columns[8]);
  *evSec=atof(columns[9]);
  *evMagSource=atof(columns[15]);
  magString=strcpy(magString,columns[13]);
//  fprintf(stderr,"assign_cols_catalog, year/month/day/hour/min/sec/Mag/MagSource/MagString: %d %d %d %d %d %.2f %f %f %s\n", *evYear, *evMon, *evDay, *evHour, *evMin, *evSec, *evMag, *evMagSource, magString);

}


/*--------------------------------------------------------------------------*/
int main (int argc, char *argv[])
/*--------------------------------------------------------------------------*/
{
  FILE *fp_eventFlatFile, *fp_catalogFile, *fp_outputFile;
  int hlines, cnt1;
  int cols_found;
  int evYear, evMon, evDay, evHour, evMin;
  int evYear2, evMon2, evDay2, evHour2, evMin2;
  int epochTimeFlatFile, epochTimeCatalog; 
  int diffSec;
  float evMag, evMag2, evSec, evSec2, evMagSource;
  float diffMag;
  char eventFlatFile[200], catalogFile[200], outputFile[200];
  char magString[20];
  char buff[BUFFLEN], buff2[BUFFLEN];
  char **columns;
  char **columns2;
  char delim[] = ",";


/* CHECK INPUT ARGUMENTS */
  if ( argc != 4 ) {
    fprintf(stderr,"USAGE: %s [events from flatfile] [catalog file] [output file]\n", argv[0]);
    fprintf(stderr,"(e.g., %s Events_complete_20170915.csv emm_c2_OK_KS_201702_wyeck_8_29_2017.csv Event_complete_associate_mag_ES.csv\n", argv[0]);
    exit(1);
  }
  sscanf(argv[1],"%s", eventFlatFile);
  sscanf(argv[2],"%s", catalogFile);
  sscanf(argv[3],"%s", outputFile);

// open files
  if ((fp_eventFlatFile = fopen(eventFlatFile, "r")) == NULL) {
    fprintf(stderr,"Could not open event list from flatfile, %s\n", eventFlatFile);
    exit(0);
  }
  if ((fp_catalogFile = fopen(catalogFile, "r")) == NULL) {
    fprintf(stderr,"Could not open catalog, %s\n", catalogFile);
    exit(0);
  }
  fp_outputFile = fopen(outputFile, "w");

// READ/APPEND EVENT-FLATFILE
// header lines
  hlines=2;
  for (cnt1=0; cnt1<hlines; cnt1++) {
    fgets(buff,BUFFLEN,fp_eventFlatFile);
    fprintf(fp_outputFile,"%s",buff);
  }
  cnt1=0;
// loop over events from flatfile
  while( fgets(buff,BUFFLEN,fp_eventFlatFile) ) {
    if ( strlen(buff) > BUFFLEN ) {
      fprintf(stderr,"Increase BUFFLEN from %d.\n", (int)BUFFLEN);
      exit(1);
    }
    buff[strcspn(buff, "\r\n")] = 0;
//     
    columns = NULL;
    cols_found = getcols(buff, delim, &columns);
    assign_cols_flatfile(columns, &evMag, &evYear, &evMon, &evDay, &evHour, &evMin, &evSec);
    free(columns);
    epochTimeFlatFile=compute_epochTime(evYear,evMon,evDay,evHour,evMin,(int)evSec);
// read through hazard catalog until match input event from flatfile
    while( fgets(buff2,BUFFLEN,fp_catalogFile) ) {
      buff2[strcspn(buff2, "\r\n")] = 0;
      columns2 = NULL;
      cols_found = getcols(buff2, delim, &columns2);
      assign_cols_catalog(columns2, &evMag2, &evYear2, &evMon2, &evDay2, &evHour2, &evMin2, &evSec2, &evMagSource, magString);
      free(columns2);
      epochTimeCatalog=compute_epochTime(evYear2,evMon2,evDay2,evHour2,evMin2,(int)evSec2);
      diffSec=abs(epochTimeFlatFile-epochTimeCatalog);
      diffMag=fabsf(evMag-evMag2);
      if ( (diffMag<0.1) && (diffSec<1)) {
        fprintf(stderr,"MATCH: ");
        fprintf(stderr,"dM=%.1f, dTime=%d\n",diffMag, diffSec);
        fprintf(stderr,"%s\n", buff); 
        fprintf(stderr,"%s\n\n", buff2); 
        fprintf(fp_outputFile,"%s,%.1f,%s\n",buff,evMagSource,magString);
        break;
      }
    }
    rewind(fp_catalogFile);
  }


// close file
  fclose(fp_eventFlatFile); 
  fclose(fp_catalogFile);
  fclose(fp_outputFile);


  return 0;
}
