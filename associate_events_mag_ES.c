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
void assign_cols_flatfile(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec, float *staLat, float *staLon);
void assign_cols_catalog(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec, float *evMagSource, char *magString, float *catLon, float *catLat, float *catDep, float *hdLon, float *hdLat, float *hdDep, float *hypoddLon, float *hypoddLat, float *hypoddDep, float *hypoESLon, float *hypoESLat, float *hypoESDep);
char * replace(char const * const original, char const * const pattern, char const * const replacement );
int compute_epochTime(int yearIn, int monthIn, int dayIn, int hourIn, int minIn, int secIn);
void calc_distances(float *stLon, float *stLat, float *evLon, float *evLat, float *evDep, float *epiDist, float *hypoDist);


/*--------------------------------------------------------------------------*/
void assign_cols_flatfile(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec, float *staLat, float *staLon)
/*--------------------------------------------------------------------------*/
{
//
  *evMag=atof(columns[20]);
  *evYear=atoi(columns[14]);
  *evMon=atoi(columns[15]);
  *evDay=atoi(columns[16]);
  *evHour=atoi(columns[17]);
  *evMin=atoi(columns[18]);
  *evSec=atof(columns[19]);
  *staLat=atof(columns[5]);
  *staLon=atof(columns[6]);
//  fprintf(stderr,"assign_cols_flatfile, year/month/day/hour/min/sec/Mag: %d %d %d %d %d %.2f %f\n", *evYear, *evMon, *evDay, *evHour, *evMin, *evSec, *evMag);

}

/*--------------------------------------------------------------------------*/
void calc_distances(float *stLon, float *stLat, float *evLon, float *evLat, float *evDep, float *epiDist, float *hypoDist)
/*--------------------------------------------------------------------------*/
{
  float lon1, lat1, lon2, lat2, dist, az, baz, dep;

//
  lon1=*stLon;
  lat1=*stLat;
  lon2=*evLon;
  lat2=*evLat;
  dep=*evDep;

//
  delaz_(&lat2,&lon2,&lat1,&lon1,&dist,&az,&baz);
  *epiDist=dist;
  *hypoDist=sqrt(dist*dist+dep*dep);

}

/*--------------------------------------------------------------------------*/
void assign_cols_catalog(char **columns, float *evMag, int *evYear, int *evMon, int *evDay, int *evHour, int *evMin, float *evSec, float *evMagSource, char *magString, float *catLon, float *catLat, float *catDep, float *hdLon, float *hdLat, float *hdDep, float *hypoddLon, float *hypoddLat, float *hypoddDep, float *hypoESLon, float *hypoESLat, float *hypoESDep)
/*--------------------------------------------------------------------------*/
{
  int len_hd, len_hypodd, len_hypoES;

//
  *evMag=atof(columns[0]);
  *evYear=atoi(columns[4]);
  *evMon=atoi(columns[5]);
  *evDay=atoi(columns[6]);
  *evHour=atoi(columns[7]);
  *evMin=atoi(columns[8]);
  *evSec=atof(columns[9]);
  *evMagSource=atof(columns[13]);
  magString=strcpy(magString,columns[13]);
// locations
  *catLon=atof(columns[1]);
  *catLat=atof(columns[2]);
  *catDep=atof(columns[3]);
//
  len_hd=strlen(columns[24]);
  len_hypodd=strlen(columns[30]);
  len_hypoES=strlen(columns[35]);
//  len_hypoES=strlen("test");

//  fprintf(stderr,"assign_cols_catalog, year/month/day/hour/min/sec/Mag/MagSource/MagString: %d %d %d %d %d %.2f %f %f %s\n", *evYear, *evMon, *evDay, *evHour, *evMin, *evSec, *evMag, *evMagSource, magString);
//fprintf(stderr,"assign_cols_catalog, cat:\n%f %f %f\n", *catLon, *catLat, *catDep);
//fprintf(stderr,"string:%s-%s-%s-\n", columns[24], columns[30], columns[35]);
//fprintf(stderr,"string lengths: %d %d %d\n", len_hd, len_hypodd, len_hypoES);

  if (len_hd>2) {
    *hdLon=atof(columns[24]);
    *hdLat=atof(columns[25]);
    *hdDep=atof(columns[26]);
  }
  else {
    *hdLon=-999;
    *hdLat=-999;
    *hdDep=-999;
  }
  if (len_hypodd>2) {
    *hypoddLon=atof(columns[30]);
    *hypoddLat=atof(columns[31]);
    *hypoddDep=atof(columns[32]);
  }
  else {
    *hypoddLon=-999;
    *hypoddLat=-999;
    *hypoddDep=-999;
  }
  if (len_hypoES>2) {
    *hypoESLon=atof(columns[35]);
    *hypoESLat=atof(columns[36]);
    *hypoESDep=atof(columns[37]);
  }
  else {
    *hypoESLon=-999;
    *hypoESLat=-999;
    *hypoESDep=-999;
  }
//  fprintf(stderr,"assign_cols_catalog, year/month/day/hour/min/sec/Mag/MagSource/MagString: %d %d %d %d %d %.2f %f %f %s\n", *evYear, *evMon, *evDay, *evHour, *evMin, *evSec, *evMag, *evMagSource, magString);
//  if (*hypoddLon > -999 ) {
//    fprintf(stderr,"assign_cols_catalog, cat:\n%f %f %f\n", *catLon, *catLat, *catDep);
//    fprintf(stderr,"assign_cols_catalog, locations cat/HD/hypoDD/hypoES:\n%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n", *catLon, *catLat, *catDep, *hdLon, *hdLat, *hdDep, *hypoddLon, *hypoddLat, *hypoddDep, *hypoESLon, *hypoESLat, *hypoESDep);
//    exit(1);
//  }

}


/*--------------------------------------------------------------------------*/
int main (int argc, char *argv[])
/*--------------------------------------------------------------------------*/
{
  FILE *fp_eventFlatFile, *fp_catalogFile, *fp_outputFile;
  int hlines, cnt1;
  int evYear, evMon, evDay, evHour, evMin;
  int evYear2, evMon2, evDay2, evHour2, evMin2;
  int epochTimeFlatFile, epochTimeCatalog; 
  int diffSec, cols_found;
  float stLon, stLat;
  float evMag, evMag2, evSec, evSec2, evMagSource;
  float catLon, catLat, catDep, hdLon, hdLat, hdDep, hypoddLon, hypoddLat, hypoddDep, hypoESLon, hypoESLat, hypoESDep;
  float catEpi, catHypo, hdEpi, hdHypo, hypoddEpi, hypoddHypo, hypoESEpi, hypoESHypo;
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
    buff[strcspn(buff, "\r\n")] = 0;
    if ( cnt1 == 0) {
      fprintf(fp_outputFile,"%s\n",buff);
    }
    else {
      fprintf(fp_outputFile,"%s,",buff);
    }
  }
// header from catalog file, ES
  hlines=1;
  for (cnt1=0; cnt1<hlines; cnt1++) {
    fgets(buff,BUFFLEN,fp_catalogFile);
    buff[strcspn(buff, "\r\n")] = 0;
    fprintf(fp_outputFile,"%s,R_epi-cat,R_hypo-cat,R_epi-HD,R_hypo-HD,R_epi-hypoDD,R_hypo-hypoDD,R_epi-hypoES,R_hypo-hypoES\n",buff);
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
    if ( cols_found>0 ) {
      assign_cols_flatfile(columns, &evMag, &evYear, &evMon, &evDay, &evHour, &evMin, &evSec, &stLat, &stLon);
      free(columns);
    }
    else {
      fprintf(stderr,"Problem reading csv file!!\n");
      exit(1);
    }
    epochTimeFlatFile=compute_epochTime(evYear,evMon,evDay,evHour,evMin,(int)evSec);
// read through hazard catalog until match input event from flatfile
    while( fgets(buff2,BUFFLEN,fp_catalogFile) ) {
      buff2[strcspn(buff2, "\r\n")] = 0;
      columns2 = NULL;
      cols_found = getcols(buff2, delim, &columns2);
      assign_cols_catalog(columns2, &evMag2, &evYear2, &evMon2, &evDay2, &evHour2, &evMin2, &evSec2, &evMagSource, magString, &catLon, &catLat, &catDep, &hdLon, &hdLat, &hdDep, &hypoddLon, &hypoddLat, &hypoddDep, &hypoESLon, &hypoESLat, &hypoESDep);
      free(columns2);
      epochTimeCatalog=compute_epochTime(evYear2,evMon2,evDay2,evHour2,evMin2,(int)evSec2);
      diffSec=abs(epochTimeFlatFile-epochTimeCatalog);
      diffMag=fabsf(evMag-evMag2);
      if ( (diffMag<0.1) && (diffSec<1)) {
//        fprintf(stderr,"MATCH: ");
//        fprintf(stderr,"dM=%.1f, dTime=%d\n",diffMag, diffSec);
//        fprintf(stderr,"%s\n", buff); 
//        fprintf(stderr,"%s\n\n", buff2); 
//        fprintf(fp_outputFile,"%s,%.1f,%s\n",buff,evMagSource,magString);
//
        fprintf(fp_outputFile,"%s,%s",buff,buff2);
// add epicentral and hypocentral distances to flatfile
//  float catEpi, catHypo, hdEpi, hdHypo, hypoddEpi, hypoddHypo, hypoESEpi, hypoESHypo;
        calc_distances(&stLon, &stLat, &catLon, &catLat, &catDep, &catEpi, &catHypo);
        fprintf(fp_outputFile,",%.2f,%.2f",catEpi,catHypo);
        if ( hdLon > -999 ) { 
          calc_distances(&stLon, &stLat, &hdLon, &hdLat, &hdDep, &hdEpi, &hdHypo);
          fprintf(fp_outputFile,",%.2f,%.2f",hdEpi,hdHypo);
        }
        else {
          fprintf(fp_outputFile,",,");
        }
        if ( hypoddLon > -999 ) { 
          calc_distances(&stLon, &stLat, &hypoddLon, &hypoddLat, &hypoddDep, &hypoddEpi, &hypoddHypo);
          fprintf(fp_outputFile,",%.2f,%.2f",hypoddEpi,hypoddHypo);
        }
        else {
          fprintf(fp_outputFile,",,");
        }
        if ( hypoESLon > -999 ) { 
          calc_distances(&stLon, &stLat, &hypoESLon, &hypoESLat, &hypoESDep, &hypoESEpi, &hypoESHypo);
          fprintf(fp_outputFile,",%.2f,%.2f,\n",hypoESEpi,hypoESHypo);
        }
        else {
          fprintf(fp_outputFile,",,,\n");
        }
//        fprintf(stderr,"%s,%s\n",buff,buff2);
//        fprintf(stderr,"%s,%.1f,%s\n",buff,evMagSource,magString);
        break;
      }
    }
    rewind(fp_catalogFile);
// read header line
    fgets(buff2,BUFFLEN,fp_catalogFile);
    cnt1++;
//if ( cnt1 > 1000 ) exit(1);
    if ( cnt1%10000 == 0 ) fprintf(stderr,"...%d\n", cnt1);
  }


// close file
  fclose(fp_eventFlatFile); 
  fclose(fp_catalogFile);
  fclose(fp_outputFile);


  return 0;
}
