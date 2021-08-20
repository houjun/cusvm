#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Stdio.h"
#include "geodataq.h"
#include "geominterp.h"

/***************************************************************************** 
 *                                                                           *
 *                            Auxiliary functions                            *
 *                                                                           *
 ****************************************************************************/ 
void linspace(double xa,double xb,int numpoints,double *result){
    int i;
    double dx;
    
    if(numpoints==1)
	dx=0;
    else
	dx=(xb-xa)/(numpoints-1);      

    for (i=0; i<numpoints; i++)	
	result[i]=xa+i*dx;     
    return;   
}

/***************************************************************************** 
 *                                                                           *
 *                           Query functions                                 *
 *                                                                           *
 ****************************************************************************/ 

/***************************************************************************** 
 *                                                                           *
 *  velocity : *  the function assumed                                       *
 *               v=vmin*pow(vmax/vmin,pow(z/h,1))                            *
 *                                                                           *
 *****************************************************************************/
double velocity(double vmin,double vmax,double z, double exponent){
    double val;    
    if(vmin == 0){
	fprintf(stdout, " Err vmin = 0 ");
	exit(1);
    }      
    val = vmin*pow(vmax/vmin,pow(z,exponent));    
    return val;    
}

/***************************************************************************** 
 *                                                                           *
 *  property_function: *  the function assumed                               *
 *               v=vmin*pow(vmax/vmin,pow(z/h,1))                            *
 *                                                                           *
 *****************************************************************************/
double property_function(int functiontype,polynomial_t *poly, double pmin,double pmax,
			 double z, double exponent){
    int i;
    double val;    

    /* Polynomial */
    switch(functiontype){
    case 0:
	val=0;
	for(i=0;i<poly->order+1;i++)
	    val=val+poly->coefficient[i]*pow(z,poly->exponent[i]);
	break;
    case 1:
       if(pmin == 0){
	fprintf(stdout, " Err vmin = 0 ");
	exit(1);
       }      
       val = pmin*pow(pmax/pmin,pow(z,exponent));    
       break;
    }

    return val;    
}

/***************************************************************************** 
 *                                                                           *
 *  density_fvp : rho as a function of vp in m/s (Brocher, 2008)             *
 *                                                                           *
 *****************************************************************************/
double density_fvp(double vp){
    double vpInKm,density;

    vp  = vp/1000;
    
    density = vp*(1.6612+vp*(-.4721+vp*(.0671+vp*(-.0043+.000106*vp))));

    return density*1000;
}    

/***************************************************************************** 
 *                                                                           *
 *  vp_fvs : in m/s (Brocher, 2008)                                          *
 *                                                                           *
 *****************************************************************************/
double vp_fvs(double vs){
    double vp;

    vs  = vs/1000;
    
    vp = 0.9409+2.0947*vs-0.8206*pow(vs,2)+0.2683*pow(vs,3)-.0251*pow(vs,4); 
   

    return vp*1000;
}    

/***************************************************************************** 
 *                                                                           *
 *  weathering function: *  the function assumed                             *
 *               v=vmin*pow(vmax/vmin,pow(z/h,1))                            *
 *                                                                           *
 *****************************************************************************/
double weathering(double min,double max,double z, double exponent){
    double val;    
    if(min == 0){
	fprintf(stdout, " Err vmin = 0 ");
	exit(1);
    }      
    val = min*pow(max/min,pow(z,exponent));    
    return val;    
}

/***************************************************************************** 
 *                                                                           *
 * Query_Topo_Bedrock_DepthToBedrock :                                       *
 *                                                                           *
 *****************************************************************************/
/*
 *  casetbdb = 0 topo
 *             1 bedrock
 *             2 depthtobedrock
 */
void Query_Topo_Bedrock_DepthToBedrock (double x,double y,double *z , database_t *database,
					int casetbdb){   

    int iPolyg,iDb;
    double geologicunitbottom, geologicunittop=0, region;
    double depth,topo,bedrockelevation,depthtobedrock;
    meshdb_t *stencilplane;
    surface_t *surface, *velocitySurface;   
    geologicunit_t *geologicunit;
    
    /* if the query is out of the database return the value of the last point in the direction
       that has not been passed if none it will return the value of the upper rigth corner*/
    
    if(database->xmin     > x  )x=database->xmin;
    if(database->xmax     < x  )x=database->xmax;
    if(database->ymin     > y  )y=database->ymin;
    if(database->ymax     < y  )y=database->ymax;
    
    /* GO THROUGH THE GEOLOGIC UNITS */
    for (iDb=0; iDb<database->howmanygeologicunits; iDb++){     	
	
	geologicunit = &(database->geologicunits[iDb]);
	surface = &(geologicunit->reflectorsurfaces);
	
	if(geologicunit->isbottomfreesurface ==1){
	    topo=
		Interp_In_Ele_Global( x,y,surface->zgrid,surface->xgrid,surface->ygrid,
				     surface->nxgrid,surface->nygrid );
	    if(casetbdb==0){
		*z=topo;
		return;
	    }
	}
	
	if(geologicunit->isbottombedrock ==1){
	    depth=
		Interp_In_Ele_Global( x,y,surface->zgrid,surface->xgrid,surface->ygrid,
				      surface->nxgrid,surface->nygrid );
	    bedrockelevation= depth;	
	    depthtobedrock  =-topo+bedrockelevation; /* axis are flipped */
	    if(casetbdb==1){
		*z= bedrockelevation;
		return;
	    }
	    if(casetbdb==2){
		*z= depthtobedrock ;
		return;
	    }
	}	
    }
    return;   
}

/***************************************************************************** 
 *                                                                           *
 *  Single_Search :                                                          *
 *                                                                           *
 *****************************************************************************/
int Single_Search(double x,  double z, double y,
		   double *vp, double *vs, double *rho , int *containingunit,
		   double *qp, double *qs, database_t *database, int elevordepthorsqueeze){   
    int iPolyg,iDb;
    double geologicunitbottom, geologicunittop=0, region, zused ;
    meshdb_t *stencilplane;
    surface_t *surface, *velocitySurface;
    geologicunit_t *geologicunit;
    double topo,bedrock,tht50,thmt50,maxTopo, depthInput;
    double vs30, surficialVsMin,surficialVsUsed,vpvssurface,vpvsdeep,vpvsfactor,slope,thresholdVs;
    double bedrockthreshold,factorVs30, depth;
    double weatheringfactor;
    polynomial_t *polynomial;
    int typeoffunction;

    /* Identify the topography and bedrock  */
    Query_Topo_Bedrock_DepthToBedrock( x,y,&topo   ,database,0);    
    Query_Topo_Bedrock_DepthToBedrock( x,y,&bedrock,database,1); 

    switch(elevordepthorsqueeze){
    case 0:
	depth=database->datum-z;
	break;
    case 1:	
	depth=topo+z;
	if(z < database->depth0)
	    depth=depth+database->depth0;
	break;
    case 2: /* testing */
	/* Compute topo to 50km thickness*/
	tht50=50e3-topo;	
	/* Compute maxtopo to 50km thickness*/
	maxTopo=(database->geologicunits[0]).reflectorsurfaces.maxsurfaceval;
	thmt50=50e3-maxTopo;       
	/* Compute the true value of depth in the database */ 
	depth=(tht50/thmt50)*z+topo;
	break;
    }    

    /* if the query is out of the database return the value of the last point 
       in the direction that has not been passed if none it will return the 
       value of the upper right corner*/    
    if(database->xmin     > x  )x=database->xmin;
    if(database->xmax     < x  )x=database->xmax;
    if(database->ymin     > y  )y=database->ymin;
    if(database->ymax     < y  )y=database->ymax;
    if(database->depthmin > depth)depth=database->depthmin;
    if(database->depthmax < depth)depth=database->depthmax;
    
    /* ****************** GO THROUGH THE GEOLOGIC UNITS  *************************/
    geologicunittop=0;
    for (iDb=0; iDb<database->howmanygeologicunits; iDb++){       
 	geologicunit = &(database->geologicunits[iDb]);
	surface = &(geologicunit->reflectorsurfaces);
	geologicunitbottom= 
	    Interp_In_Ele_Global( x,y,surface->zgrid,surface->xgrid,surface->ygrid,
				  surface->nxgrid,surface->nygrid );

	if(geologicunitbottom < topo)
	    geologicunitbottom=topo;

	if(depth < geologicunitbottom && ((database->geologicunits[iDb]).isexcluded ==0 ) ){

	    *containingunit=iDb;

	    switch(geologicunit->znormalizationtype){
		    
	    case 0: /* Z IS REFERED TO THE TOP OF THE UNIT */
		zused = depth-geologicunittop;
		break;		    
	    case 1: /* NORMALIZED TO THE TOP OF THE GEOLOGIC UNIT */
		zused = (depth-geologicunittop)/geologicunitbottom;
		break;		    
	    case 2: /*  ABSOLUTE FROM THE SURFACE OF THE DATABASE*/
		zused = depth;
		break;		    
	    case 3: /* Z IS REFERED TO THE TOP OF THE UNIT */
		zused = depth/ geologicunitbottom;
		break;		    
	    case 4: /* Z IS normalized within the unit*/
		zused = (depth-geologicunittop)/(geologicunitbottom-geologicunittop);
		break;
	    case 5: /* Z IS REFERED TO THE TOPOGRAPHY*/
		zused = (depth-topo);
		break;	   
	    }
	    if( zused < 0){
		fprintf(stdout,"Error");
		exit(1);
	    }

	    if(geologicunit->profiletype==0){
		
		velocitySurface= &((database->geologicunits[iDb]).veldistribsurface);		
		region=Interp_In_Ele_Global_Vel(x,y,iDb,velocitySurface->zgrid,
						velocitySurface->xgrid,
						velocitySurface->ygrid,
						velocitySurface->nxgrid,
						velocitySurface->nygrid,
						geologicunit->numprofiles,geologicunit->velprofiles,
						vp,vs,rho,zused );
		
		if(*vp<=0 || *vs<=0 || *rho<=0){
		    fprintf(stdout,"Error: vp=%f vs=%f rho=%f zused=%f surface=%d"
			    ,*vp,*vs,*rho,zused,iDb);
		    exit(1);
		}
	    }

	    if(geologicunit->profiletype==1){
		typeoffunction = (database->geologicunits[iDb]).velprofilefunction;
		switch (typeoffunction) 
		    {
		    case 0:
			polynomial=&(database->geologicunits[iDb]).poly;    				    
			*vs=property_function(typeoffunction, polynomial,0,0,zused,0);			
			/*vpvssurface=(database->geologicunits[iDb]).vpvsa;
			  vpvsdeep   =(database->geologicunits[iDb]).vpvsb;
			  if(elevordepthorsqueeze == 2)
			  zused = (depthInput)/fabs(geologicunitbottom-topo);
			  else
			  zused = (depth-topo)/fabs(geologicunitbottom-topo);
			  
			  if(zused < 0)zused=0;
			  if(zused >1 )zused=1;
			  
			  vpvsfactor=(vpvsdeep-vpvssurface)*zused+vpvssurface;*/
			
			*vp=vp_fvs(*vs);			
			*rho=density_fvp(*vp);
			break;
		    case 1:			
			*vs =velocity(geologicunit->minvs,geologicunit->maxvs  ,
				      zused, geologicunit->exponentprofile );
			*vp =velocity(geologicunit->minvp,geologicunit->maxvp  ,
				      zused, geologicunit->exponentprofile );
			*rho=velocity(geologicunit->minrho,geologicunit->maxrho,
				      zused, geologicunit->exponentprofile );
			
			break;
		    }
	    }
	    
	    if(geologicunit->profiletype==2){
		
		velocitySurface= &((database->geologicunits[iDb]).vs30surface);
		surficialVsMin= (database->geologicunits[iDb]).minvs;
		thresholdVs= (database->geologicunits[iDb]).maxvs;
		vpvssurface=(database->geologicunits[iDb]).vpvsa;
		vpvsdeep   =(database->geologicunits[iDb]).vpvsb;
		factorVs30=(database->geologicunits[iDb]).factorvs30;
		typeoffunction = (database->geologicunits[iDb]).velprofilefunction;
		
		vs30=Interp_In_Ele_Global( x,y,velocitySurface->zgrid,velocitySurface->xgrid,
					   velocitySurface->ygrid,velocitySurface->nxgrid,
					   velocitySurface->nygrid );
		
		if(typeoffunction == 0){
		    polynomial=&(database->geologicunits[iDb]).poly;    				    
		    *vs=vs30*property_function(typeoffunction, polynomial,0,0,zused,0);
		}

		/* use the same variable zused to get the vp to vs ratio note it is normalized*/
		if(elevordepthorsqueeze == 2)
		    zused = (depthInput)/fabs(geologicunitbottom-topo);
		else
		    zused = (depth-topo)/fabs(geologicunitbottom-topo);
		
		if(zused < 0)zused=0;
		if(zused >1 )zused=1;
		
		vpvsfactor=(vpvsdeep-vpvssurface)*zused+vpvssurface;

		*vp=(*vs)*vpvsfactor;
                *rho=density_fvp(*vp);	

	    }
	    
	    /* WEATHERING  AND THRESHOLDS*/
	    if(geologicunit->profiletype==1 || geologicunit->profiletype==2 ){
		
		if(*vs < (database->geologicunits[iDb]).minvs)
		    *vs=(database->geologicunits[iDb]).minvs;
		if(*vs > (database->geologicunits[iDb]).maxvs)
		    *vs=(database->geologicunits[iDb]).maxvs;
		
		if(*vp < (database->geologicunits[iDb]).minvp)
		    *vp=(database->geologicunits[iDb]).minvp;
		if(*vp > (database->geologicunits[iDb]).maxvp)
		    *vp=(database->geologicunits[iDb]).maxvp;

		if(*rho < (database->geologicunits[iDb]).minrho)
		    *rho=(database->geologicunits[iDb]).minrho;
		if(*rho > (database->geologicunits[iDb]).maxrho)
		    *rho=(database->geologicunits[iDb]).maxrho;

		zused= (depth-topo)/fabs(geologicunitbottom-topo);
		weatheringfactor=weathering(geologicunit->weatheringfactor,1,zused,geologicunit->weatheringexp);
		*vs=*vs*weatheringfactor;
		*vp=*vp*weatheringfactor;
		*rho=*rho*weatheringfactor;

		if(*vs < (database->geologicunits[iDb]).minvs)
		    *vs=(database->geologicunits[iDb]).minvs;
		if(*vs > (database->geologicunits[iDb]).maxvs)
		    *vs=(database->geologicunits[iDb]).maxvs;
		
		if(*vp < (database->geologicunits[iDb]).minvp)
		    *vp=(database->geologicunits[iDb]).minvp;
		if(*vp > (database->geologicunits[iDb]).maxvp)
		    *vp=(database->geologicunits[iDb]).maxvp;
		
		if(*rho < (database->geologicunits[iDb]).minrho)
		    *rho=(database->geologicunits[iDb]).minrho;
		if(*rho > (database->geologicunits[iDb]).maxrho)
		    *rho=(database->geologicunits[iDb]).maxrho;
	    }	    

	    return 0;
	}
	
	if(geologicunitbottom>geologicunittop && ((database->geologicunits[iDb]).isexcluded ==0))
	    geologicunittop=geologicunitbottom;
    }
    
    /* if the query is out of the database return the value of the last point in the direction
       that has not been passed if none it will return the value of the upper rigth corner, if
       the query goes upto here something is wrong and will print the default values*/      

    *vp=-99999;
    *vs=-99999;
    *rho=-99999;
    return 0;   
}

/***************************************************************************** 
 *                                                                           *
 *   extract_print_point_query:  Function to query properties at a point     *
 *                                                                           *
 *****************************************************************************/
int extract_print_point_query(char **argv, database_t *database, int squeeze){
    static char auxString[256];
    double Lon,Lat,depth,vp,vs,rho,qp,qs;
    int containingunit;
    int iStation, numberofstations;
    FILE *fpOutProp[3],*fpOutQpQs, *fpInputQuery;
    
    strcpy(auxString,argv[4]);
    fpInputQuery = fopen(auxString,"r");    
    sprintf(auxString, "%s/output.out", argv[5]);
    fpOutProp[0] = fopen(auxString,"w");
    fscanf( fpInputQuery,"%d",&numberofstations);
    for (iStation=0; iStation<numberofstations;iStation++){
	fscanf(fpInputQuery,"%lf %lf %lf", &Lon,&Lat,&depth);
	Single_Search( Lat,depth,Lon,
		       &vp, &vs, &rho , &containingunit, &qp, &qs, database,squeeze); 
	/* THE DAMPING IS NOT IMPLEMENTED IT IS DECIDED BY THE USER */
	
	fprintf(fpOutProp[0],"\n%f %f %f %f %f %f %d", Lon,Lat,depth,vp,vs,rho,
                                                       containingunit);	
    }        
    
    fclose(fpOutProp[0]);
    
    return 1;
}

/***************************************************************************** 
 *                                                                           *
 *         extract_print_bedrockelevation                                    *
 *                                                                           *
 *****************************************************************************/
int extract_print_bedrockelevation(char **argv, database_t *database){
    
    double Lon,Lat,depth,vp,vs,rho,qp,qs;
    double stepLon,stepLat,stepD,lengthLon,lengthLat, totDepth;
    double originLon, originLat, originDepth;
    double X,Y,Z,xRotated[3];
    double finalLon, finalLat, finalDepth;

    int numAlong;
    int i, j,k, iStation, numberofstations,numLon,numLat,numD,iLon,iLat,iD;
    int typeofquery;
    int hov; /*HORIZONTAL (0) OR VERTICAL (1) PLANE TO QUERY */
    int velorveldepthorbdrortopo;

    FILE *fpOutProp[3],*fpOutQpQs, *fpInputQuery;
    FILE *fpBdr;
    static char inputPlanesQueryPath[256], outputPath[256],auxString[256];
    
    double *longitudes, *latitudes, *depths;
    
    sprintf(auxString, "%s/bedrock.out", argv[5]);
    fpBdr=fopen(auxString,"w");
    
    strcpy(inputPlanesQueryPath,argv[4]);
    fpInputQuery = fopen(inputPlanesQueryPath,"r");    
    
    /* DEFINITION OF THE BOX WHERE WE ARE QUERYING THE DATA */
    /* IT CAN BE A HORIZONTAL PLANE ( 0 ) OR A VERTICAL PLANE (1)*/	
    Parse_Text( fpInputQuery,"hororvert"  , 'i',&hov);
    Parse_Text( fpInputQuery,"originLon"  , 'd',&originLon);	 
    Parse_Text( fpInputQuery,"originLat"  , 'd',&originLat);
    
    if(hov == 0 ){
	Parse_Text( fpInputQuery,"lengthLon", 'd',&lengthLon);
	Parse_Text( fpInputQuery,"lengthLat", 'd',&lengthLat);
	Parse_Text( fpInputQuery,"numLon", 'i',&numLon);
	Parse_Text( fpInputQuery,"numLat", 'i',&numLat);
	if(numLon > 1) stepLon = lengthLon /(numLon-1);
	else           stepLon = lengthLon;
	if(numLat > 1) stepLat = lengthLat /(numLat-1);
	else           stepLat = lengthLat;
	numD = 1;
	stepD   = 1;
	
	/* OUTPUT PLANES */
	/* WRITE HEADER  */
	fprintf(fpBdr,"%d %e %e 0 %d %d 1 %e %e %e",hov,originLon, originLat,
		numLon,    numLat,lengthLon, lengthLat,0.0);	    	
	
	for (iLon=0;iLon < numLon; iLon++){	
	    Lon=originLon+stepLon*iLon;		
	    for (iLat=0;iLat < numLat; iLat++){	    
		Lat= originLat+stepLat*iLat;		
		Query_Topo_Bedrock_DepthToBedrock( Lat,Lon,&depth, database,1);		
		fprintf(fpBdr,"%e ",depth);		
	    }		
	}
    }
    if(hov == 1 ){
	
	Parse_Text( fpInputQuery,"finalLon"  , 'd',&finalLon);
	Parse_Text( fpInputQuery,"finalLat"  , 'd',&finalLat);
	Parse_Text( fpInputQuery,"numAlong", 'i',&numAlong);
	longitudes=malloc(sizeof(double)*numAlong);
	latitudes =malloc(sizeof(double)*numAlong);
	
	linspace(originLon  ,finalLon,numAlong,longitudes);
	linspace(originLat  ,finalLat,numAlong,latitudes );
	
	/* PRINT HEADER */
	fprintf(fpBdr,"%d %e %e %e %e %d \n",hov,originLon,originLat,
		finalLon,finalLat, numAlong);
	
	for (i=0;i<numAlong;i++){
	    Query_Topo_Bedrock_DepthToBedrock(latitudes[i],longitudes[i],&depth,
					      database,1);		
	    fprintf(fpBdr,"%e ",depth);	    
	}		
    }   
    return 0;
}

/***************************************************************************** 
 *                                                                           *
 *       extract_print_topoelevation                                         *
 *                                                                           *
 *****************************************************************************/
int extract_print_topoelevation(char **argv, database_t *database){
    
    double Lon,Lat,depth,vp,vs,rho,qp,qs;
    double stepLon,stepLat,stepD,lengthLon,lengthLat,totDepth;
    double originLon,originLat,originDepth;
    double X,Y,Z,xRotated[3];
    double finalLon,finalLat,finalDepth;

    int numAlong;
    int i,j,k,iStation,numberofstations,numLon,numLat,numD,iLon,iLat,iD;
    int typeofquery;
    int hov; /*HORIZONTAL (0) OR VERTICAL (1) PLANE TO QUERY */
    int velorveldepthorbdrortopo;

    FILE *fpOutProp[3],*fpOutQpQs, *fpInputQuery;
    FILE *fpBdr;
    static char inputPlanesQueryPath[256], outputPath[256],auxString[256];

    double *longitudes, *latitudes, *depths;

    sprintf(auxString, "%s/topo.out", argv[5]);
    fprintf(stdout,"%s\n",auxString);
    fpBdr=fopen(auxString,"w");

    strcpy(inputPlanesQueryPath,argv[4]);
    fpInputQuery = fopen(inputPlanesQueryPath,"r");
        
    /* DEFINITION OF THE BOX WHERE WE ARE QUERYING THE DATA */
    /* IT CAN BE A HORIZONTAL PLANE ( 0 ) OR A VERTICAL PLANE (1)*/	
    Parse_Text( fpInputQuery,"hororvert"  , 'i',&hov);
    Parse_Text( fpInputQuery,"originLon"  , 'd',&originLon);	 
    Parse_Text( fpInputQuery,"originLat"  , 'd',&originLat);

    if(hov == 0 ){
	Parse_Text( fpInputQuery,"lengthLon", 'd',&lengthLon);
	Parse_Text( fpInputQuery,"lengthLat", 'd',&lengthLat);
	Parse_Text( fpInputQuery,"numLon", 'i',&numLon);
	Parse_Text( fpInputQuery,"numLat", 'i',&numLat);
	if(numLon > 1) stepLon = lengthLon /(numLon-1);
	else           stepLon = lengthLon;
	if(numLat > 1) stepLat = lengthLat /(numLat-1);
	else           stepLat = lengthLat;
	numD = 1;
	stepD   = 1;
	
	/* OUTPUT PLANES */
	/* WRITE HEADER  */
	fprintf(fpBdr,"%d %e %e 0 %d %d 1 %e %e %e",hov,originLon, originLat,
		numLon,    numLat,lengthLon, lengthLat,0.0);	    	
	
	for (iLon=0;iLon < numLon; iLon++){	
	    Lon=originLon+stepLon*iLon;		
	    for (iLat=0;iLat < numLat; iLat++){	    
		Lat= originLat+stepLat*iLat;
		Query_Topo_Bedrock_DepthToBedrock(Lat,Lon,&depth, database,0);		
		fprintf(fpBdr,"%e ",depth);		
	    }		
	}
    }
    if(hov == 1 ){	
	Parse_Text( fpInputQuery,"finalLon"  , 'd',&finalLon);
	Parse_Text( fpInputQuery,"finalLat"  , 'd',&finalLat);
	Parse_Text( fpInputQuery,"numAlong", 'i',&numAlong);
	longitudes=malloc(sizeof(double)*numAlong);
	latitudes =malloc(sizeof(double)*numAlong);
	
	linspace(originLon  ,finalLon,numAlong,longitudes);
	linspace(originLat  ,finalLat,numAlong,latitudes );
	
	/* PRINT HEADER */
	fprintf(fpBdr,"%d %e %e %e %e %d \n",hov,originLon,originLat,
		finalLon,finalLat, numAlong);

	for (i=0;i<numAlong;i++){
	    Query_Topo_Bedrock_DepthToBedrock(latitudes[i],longitudes[i],&depth,
					      database,0);		
	    fprintf(fpBdr,"%e ",depth);	    
	}   
    }   
    return 0;
}

/***************************************************************************** 
 *                                                                           *
 *       extract_print_bedrockdepth                                          *
 *                                                                           *
 *****************************************************************************/
int extract_print_bedrockdepth(char **argv, database_t *database){
    
    double Lon,Lat,depth,vp,vs,rho,qp,qs;
    double stepLon,stepLat,stepD,lengthLon,lengthLat, totDepth;
    double originLon, originLat, originDepth;
    double X,Y,Z,xRotated[3];
    double finalLon, finalLat, finalDepth;

    int numAlong;
    int i, j,k, iStation, numberofstations,numLon,numLat,numD,iLon,iLat,iD;
    int typeofquery;
    int hov; /*HORIZONTAL (0) OR VERTICAL (1) PLANE TO QUERY */
    int velorveldepthorbdrortopo;

    FILE *fpOutProp[3],*fpOutQpQs, *fpInputQuery;
    FILE *fpBdr;
    static char inputPlanesQueryPath[256], outputPath[256],auxString[256];

    double *longitudes, *latitudes, *depths;

    sprintf(auxString, "%s/bedrockdepth.out", argv[5]);
    fpBdr=fopen(auxString,"w");

    strcpy(inputPlanesQueryPath,argv[4]);
    fpInputQuery = fopen(inputPlanesQueryPath,"r");
    
    /* DEFINITION OF THE BOX WHERE WE ARE QUERYING THE DATA */
    /* IT CAN BE A HORIZONTAL PLANE ( 0 ) OR A VERTICAL PLANE (1)*/	
    Parse_Text( fpInputQuery,"hororvert"  , 'i',&hov);
    Parse_Text( fpInputQuery,"originLon"  , 'd',&originLon);	 
    Parse_Text( fpInputQuery,"originLat"  , 'd',&originLat);

    if(hov == 0 ){
	Parse_Text( fpInputQuery,"lengthLon", 'd',&lengthLon);
	Parse_Text( fpInputQuery,"lengthLat", 'd',&lengthLat);
	Parse_Text( fpInputQuery,"numLon", 'i',&numLon);
	Parse_Text( fpInputQuery,"numLat", 'i',&numLat);
	if(numLon > 1) stepLon = lengthLon /(numLon-1);
	else           stepLon = lengthLon;
	if(numLat > 1) stepLat = lengthLat /(numLat-1);
	else           stepLat = lengthLat;
	numD = 1;
	stepD   = 1;
	
	/* OUTPUT PLANES */
	/* WRITE HEADER  */
	fprintf(fpBdr,"%d %e %e 0 %d %d 1 %e %e %e",hov,originLon, originLat,
		numLon,    numLat,lengthLon, lengthLat,0.0);	    	
	
	for (iLon=0;iLon < numLon; iLon++){	
	    Lon=originLon+stepLon*iLon;		
	    for (iLat=0;iLat < numLat; iLat++){	    
		Lat= originLat+stepLat*iLat;
		Query_Topo_Bedrock_DepthToBedrock(Lat,Lon,&depth, database,2);		
		fprintf(fpBdr,"%e ",depth);		
	    }		
	}
    }
    if(hov == 1 ){
	
	Parse_Text( fpInputQuery,"finalLon"  , 'd',&finalLon);
	Parse_Text( fpInputQuery,"finalLat"  , 'd',&finalLat);
	Parse_Text( fpInputQuery,"numAlong", 'i',&numAlong);
	longitudes=malloc(sizeof(double)*numAlong);
	latitudes =malloc(sizeof(double)*numAlong);
	
	linspace(originLon  ,finalLon,numAlong,longitudes);
	linspace(originLat  ,finalLat,numAlong,latitudes );
	
	/* PRINT HEADER */
	fprintf(fpBdr,"%d %e %e %e %e %d \n",hov,originLon,originLat,
		finalLon,finalLat, numAlong);

	for (i=0;i<numAlong;i++){
	    Query_Topo_Bedrock_DepthToBedrock(latitudes[i],longitudes[i],&depth,database,2);		
	    fprintf(fpBdr,"%e ",depth);	    
	}    
    }   
    return 0;
}

/***************************************************************************** 
 *                                                                           *
 *        Function to query properties grids                                 *
 *                                                                           *
 *****************************************************************************/
int extract_print_properties_grid(char **argv, database_t *database,int elevordepthorsqueeze){
    
    double Lon,Lat,depth,vp,vs,rho,qp,qs;
    double stepLon,stepLat,stepD,lengthLon,lengthLat, totDepth;
    double originLon, originLat, originDepth;
    double finalLon, finalLat, finalDepth;
    double *longitudes, *latitudes, *depths;
    int containingunit;
    int numAlong;
    int i, j, k,iUnit, iStation, numberofstations,numLon,numLat,numD,iLon,iLat,iD;
    int typeofquery;
    int hov; /*HORIZONTAL (0) OR VERTICAL (1) PLANE TO QUERY */

    FILE *fpOutProp[4],*fpOutQpQs, *fpInputQuery;
    FILE *fpBdr;
    static char inputPlanesQueryPath[256], outputPath[256],auxString[256];

    strcpy(inputPlanesQueryPath,argv[4]);
    fpInputQuery = fopen(inputPlanesQueryPath,"r");

    /* AUXILIARY FILES TO OUTPUT PROPERTIES */
    sprintf(auxString, "%s/vp.out", argv[5]);
    fpOutProp[0] = fopen(auxString,"w");
    sprintf(auxString, "%s/vs.out", argv[5]);
    fpOutProp[1] = fopen(auxString,"w");
    sprintf(auxString, "%s/rho.out", argv[5]);
    fpOutProp[2] = fopen(auxString,"w");
    sprintf(auxString, "%s/containingunit.out", argv[5]);
    fpOutProp[3] = fopen(auxString,"w");

    sprintf(auxString, "%s/qpqs.out", argv[5]);
    fpOutQpQs  = fopen(auxString,"w");
    
    /* DEFINITION OF THE BOX WHERE WE ARE QUERYING THE DATA */
    /* IT CAN BE A HORIZONTAL PLANE ( 0 ) OR A VERTICAL PLANE (1)*/	
    Parse_Text( fpInputQuery,"hororvert"  , 'i',&hov);
    Parse_Text( fpInputQuery,"originLon"  , 'd',&originLon);	 
    Parse_Text( fpInputQuery,"originLat"  , 'd',&originLat);
    Parse_Text( fpInputQuery,"originDepth", 'd',&originDepth);
        
    if(hov == 0 ){
	Parse_Text( fpInputQuery,"lengthLon", 'd',&lengthLon);
	Parse_Text( fpInputQuery,"lengthLat", 'd',&lengthLat);
	Parse_Text( fpInputQuery,"numLon", 'i',&numLon);
	Parse_Text( fpInputQuery,"numLat", 'i',&numLat);
	if(numLon > 1) stepLon = lengthLon /(numLon-1);
	else           stepLon = lengthLon;
	if(numLat > 1) stepLat = lengthLat /(numLat-1);
	else           stepLat = lengthLat;
	numD = 1;
	stepD   = 1;
	
	/* OUTPUT PLANES */
	/* WRITE HEADER  */
	for (i=0;i<4;i++){
	    fprintf(fpOutProp[i],"%d %e %e %e %d %d %d %e %e %e\n",hov,
		    originLon, originLat,originDepth,
		    numLon,    numLat,       numD,
		    lengthLon, lengthLat,   totDepth);
	    if(i==3){
		fprintf(fpOutProp[i],"%d\n",database->numberofobjects);
		for (iUnit =0;iUnit<database->numberofobjects; iUnit++) {
		    fprintf(fpOutProp[i],"%s \n",database->geologicunits[iUnit].geolunitname);
		}
	    }	    
	}
	for (iLon=0;iLon < numLon; iLon++){	    	
	    Lon=originLon+stepLon*iLon;		
	    for (iLat=0;iLat < numLat; iLat++){	    
		Lat= originLat+stepLat*iLat;
		
		for (iD=0; iD < numD ; iD++){	      
		    depth=originDepth+stepD*iD;
		    Single_Search( Lat,depth,Lon,
				   &vp, &vs, &rho , &containingunit, 
                                   &qp, &qs, database,elevordepthorsqueeze);		
		    fprintf(fpOutProp[0]," %e ", vp );
		    fprintf(fpOutProp[1]," %e ", vs );
		    fprintf(fpOutProp[2]," %e ", rho);
		    fprintf(fpOutProp[3]," %d ", containingunit);
		}	
	    }
	}	
    }

    if(hov == 1 ){
	Parse_Text( fpInputQuery,"finalLon"  , 'd',&finalLon);
	Parse_Text( fpInputQuery,"finalLat"  , 'd',&finalLat);
	Parse_Text( fpInputQuery,"finalDepth", 'd',&finalDepth);
	Parse_Text( fpInputQuery,"numAlong", 'i',&numAlong);
	Parse_Text( fpInputQuery,"numD"    , 'i',&numD);	
	longitudes=malloc(sizeof(double)*numAlong);
	latitudes =malloc(sizeof(double)*numAlong);
	depths    =malloc(sizeof(double)*numD);
	linspace(originLon  ,finalLon,numAlong,longitudes);
	linspace(originLat  ,finalLat,numAlong,latitudes );
	linspace(originDepth,finalDepth,  numD,depths    );
	/* PRINT HEADER */
	for (i=0;i<4;i++){
	    fprintf(fpOutProp[i],"%d %e %e %e %e %e %e %d %d \n",hov,
		    originLon, originLat,originDepth,finalLon,finalLat,finalDepth,
		    numAlong,       numD);
	    
	    if(i==3){
		fprintf(fpOutProp[i],"%d\n",database->numberofobjects);	
		for (iUnit =0;iUnit<database->numberofobjects; iUnit++) {
		    fprintf(fpOutProp[i],"%s \n",database->geologicunits[iUnit].geolunitname);
		}
	    
	    }

	}
	for (i=0;i<numAlong;i++)
	    
	    for (j=0;j<numD;j++){				   
		Single_Search( latitudes[i],depths[j],longitudes[i],
			       &vp, &vs, &rho ,&containingunit, &qp, &qs,
			       database,elevordepthorsqueeze); 
		fprintf(fpOutProp[0]," %e ", vp );
		fprintf(fpOutProp[1]," %e ", vs );
		fprintf(fpOutProp[2]," %e ", rho);
		fprintf(fpOutProp[3]," %d ", containingunit);
	    }	
    }    
    return 0;
}
