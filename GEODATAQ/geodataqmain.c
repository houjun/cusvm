/******************************************************************************
 * Code: geodataq_version2
 *
 * Spatial Velocity Model Database Manager. 
 *
 * Contact:
 * Leonardo Ramirez-Guzman
 * lramirezg@iingen.unam.mx or leoramirezg@gmail.com
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "geodataq.h"

void Init_Database(database_t *current, char *databasePath, int visualcheckingoffon);
int extract_print_point_query(char **argv, database_t *database, int squeeze);
int extract_print_bedrockelevation(char **argv, database_t *database);
int extract_print_topoelevation(char **argv, database_t *database);
int extract_print_bedrockdepth(char **argv, database_t *database);
int extract_print_properties_grid(char **argv, database_t *database, int elevordepthorsqueeze);

/*
 * print_usage()
 */
void  print_usage(){
  fputs ( "Usage: geodataq typeofquery velorveldepthorbdrortopo db in rout \n "
	  "              typeofquery:   (0) points or (1) lines or planes\n"
          "                                 Points is only applicable to velorveldepthorbdrortopo 0, 1,\n"
          "                                 and 2\n"
	  " velorveldepthorbdrortopo:   (0) velocity & density (elevation above sea level) \n"
	  "                             (1) velocity & density (depth from surface) \n"
	  "                             (2) velocity & density (depth from surface) \n"
          "                                   topography squeezed to mean sea level) \n"
          "                                   Squeezing scales the depth by: \n"
          "                                   (50e3-topo[lon_point][lat_point])/(50e3-elevation_max) \n"
	  "                             (3) elevation topo (relative to datum, positive down)\n"
	  "                             (4) elevation bedrock (relative to datum, positive down)\n"
	  "                             (5) depth to bedrock\n"
	  "                       db:   path to a database\n"
	  "                       in:   input file, query specifications\n"
	  "                     rout:   output path (directory), \n"
	  "                             output file named output.out\n",stderr);
  exit(1);
}


/***************************************************************************** 
 *                             Main                                          
 *                                                                           
 *   Notes:                                                                  
 *         Boundaries of units have to be given in geodetic coordinates,     
 *         i.e. lon,lat, elevation (msl), they will be shifted to datum      
 *         to use depth from theDatum value.                                 
 *                                                                           
 *****************************************************************************/
int main (int argc, char **argv){

  int typeofquery, velorveldepthorbdrortopo, squeeze;    
  FILE *fpInputQuery;
  static char inputPlanesQueryPath[256], outputPath[256],databasePath[256];
  database_t cuscvm; 

  /* check input args */
  if (argc != 6) print_usage();
   
  typeofquery             =atoi(argv[1]);
  velorveldepthorbdrortopo=atoi(argv[2]);
   
  if(typeofquery == 0 && velorveldepthorbdrortopo > 2){ /* point by point */
      fprintf(stdout,"** Invalid typeofquery (%d) for velorveldepthorbdrortopo (%d)\n\n",typeofquery,velorveldepthorbdrortopo);
      print_usage(); 
  }
  strcpy(databasePath          ,argv[3]);
  strcpy(inputPlanesQueryPath  ,argv[4]);
  strcpy(outputPath            ,argv[5]);
    
  fpInputQuery = fopen(inputPlanesQueryPath,"r");
    
  Init_Database(&cuscvm,databasePath,1);    
  fprintf(stdout,"\n"); fflush(stdout);

  if(typeofquery == 0){ /* point by point */
    switch(velorveldepthorbdrortopo) 
      { /*  properties can include velocity, density or containing unit */
      case 0: extract_print_point_query(argv,&cuscvm,0);break; 
      case 1: extract_print_point_query(argv,&cuscvm,1);break; 
      case 2: extract_print_point_query(argv,&cuscvm,2);break; 
      }
  }
  else{
    switch(velorveldepthorbdrortopo)
      {	/*  properties can include velocity, density or containing unit */
      case 0: extract_print_properties_grid (argv,&cuscvm,0);break;
      case 1: extract_print_properties_grid (argv,&cuscvm,1);break;
      case 2: extract_print_properties_grid (argv,&cuscvm,2);break;
      case 3: extract_print_topoelevation   (argv,&cuscvm)  ;break;
      case 4: extract_print_bedrockelevation(argv,&cuscvm)  ;break;
      case 5: extract_print_bedrockdepth    (argv,&cuscvm)  ;break;
      }
  }    
  return 0;
} 

