/****************************************************************************
 *
 * MODULE:       r.openess 
 * AUTHOR(S):    Francisco Alonso-Sarría
 * PURPOSE:      A modification of the Yokohama et al (2002) Openess parameter.
 * COPYRIGHT:    
 *
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 *****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <math.h>
#include <stdbool.h>

#define PI 3.141592

//#define MAIN
//#include "method.h"

double min(double a,double b){if (a<b){return(a);}else{return(b);}}
double max(double a,double b){if (a>b){return(a);}else{return(b);}}

//Calculate difference among angles
float difang(float a1,float a2){ 
     float d;
     if(fabs(a2-a1)<180) d=fabs(a2-a1); else d=360-fabs(a2-a1);
     return(d);
}

// Obtain dx and dy equivalent from a r.watershed flow direction
void asp2dxdy(int asp,int *dx, int *dy){ 
          if(asp<4){*dy=-1;*dx=(asp-2)*-1;}
          if(asp>4 && asp<8){*dy=1;*dx=asp-6;}
          if(asp==4){*dy=0;*dx=-1;}
          if(asp==8){*dy=0;*dx=1;}
}      
 
// Obtain the flow direction that a given neighbour cell (defined by dx and dy) should have to drain into the analyzed cell
int dxdy2antiasp(int dx,int dy){
    int r;
    if(dy==-1)r=dx*(-1)+6;
    if(dy==0)r=dx*(-2)+6;
    if(dy==1)r=dx+2;
    return(r);
}

// Transform dx y dy in angle (0-360)
float dxdy2ang(int dx, int dy){
   float ddx,ddy, ang;
   ddx=(float)dx;
   ddy=(float)dy;
   if (dx==0) ddx=0.000001;
   if (dy==0) ddy=0.000001;
   if ( (ddy>0) && (ddx>0) ) ang=90 - 180*atan2(ddy,ddx)/PI;
   if ( (ddy<0) && (ddx>0) ) ang=90+fabs(180*atan2(ddy,ddx)/PI);
   if ( (ddy<0) && (ddx<0) ) ang=90+fabs(180*atan2(ddy,ddx)/PI);
   if ( (ddy>0) && (ddx<0) ) ang=450 - 180*atan2(ddy,ddx)/PI;
   return(ang);
}

// Calculate weights for angles
void ponderAngulos(float angup, float angdown, int tam, float *angulos, float *wdown, float *wup, float *w){
     int i;
     float an1, an2, da1, da2, da3, da4, ext1, ext2, swdown, swup, sw, *w1, *w2;     

     w1 = (float *) G_malloc(tam * sizeof(float)); 
     w2 = (float *) G_malloc(tam * sizeof(float)); 

     if (angdown>angup){
          an1=(angdown+angup)/2;an2=angdown+((360-angdown)+angup)/2;if(an2>360) an2=an2-360;
          ext1=angdown-angup-90;ext2=180-ext1;
     }
     if (angdown<angup){
          an1=(angup+angdown)/2;an2=angup+((360-angup)+angdown)/2;if(an2>360) an2=an2-360;
          ext1=angup-angdown-90;ext2=180-ext1;
     }

     //printf("\n\nponderando... up=%f down=%f an1=%f an2=%f ext1=%f ext2=%f\n",angup, angdown, an1, an2, ext1, ext2);
     swdown=0;swup=0;sw=0;
     for (i=0; i< tam; i=i+1){
         *(wdown+i) = 0;*(wup+i) = 0;*(w1+i) = 0;*(w2+i) = 0;
         da1=difang(*(angulos+i),angdown);
         if (da1<90) *(wdown+i) = max(0,cos(2*PI*da1/180));
         da2=difang(*(angulos+i),angup);
         if (da2<90) *(wup+i) = max(0,cos(2*PI*da2/180));
         da3=difang(*(angulos+i),an1);
         if (da3<(ext1/2)) *(w1+i) = max(0,cos(2*PI*da3/180)); 
         da4=difang(*(angulos+i),an2);
         if (da4<(ext2/2)) *(w2+i) = max(0,cos(2*PI*da4/180)); 
         *(w+i)=*(w1+i) + *(w2+i);
         swdown=swdown + *(wdown+i); swup=swup+*(wup+i); sw=sw + *(w+i);
    
        // printf("    %d   %f    %f %f %f %f       %f %f %f %f \n",i,*(angulos+i),da1,da2,da3,da4, *(wup+i),*(wdown+i),*(w1+i),*(w2+i));
     }
     for (i=0; i< tam; i=i+1){*(wdown+i)=*(wdown+i)/swdown;*(wup+i)=*(wup+i)/swup;*(w+i)=*(w+i)/sw;}

     G_free(w1);
     G_free(w2);
}

int main(int argc, char *argv[]){

    struct GModule *module;
    struct Option *mapab, *mapad;
    struct Option *mapaomx, *mapa8mx, *mapaupmx, *mapadnmx, *mapabnmx;
    struct Option *mapaomn, *mapa8mn, *mapaupmn, *mapadnmn, *mapabnmn;
    struct Option *ventana, *mode1;

    // Output maps are FCELL, input map may be any type
    CELL  *datosd;
    CELL  *datosb_c;
    FCELL *datosb_f;
    DCELL *datosb_d; 
    CELL  *entrada_dir;
    CELL  *entrada_c;
    FCELL *entrada_f;
    DCELL *entrada_d; 
    FCELL *datosomx, *datos8mx, *datosupmx, *datosdnmx, *datosbnmx, *datosomn, *datos8mn, *datosupmn, *datosdnmn, *datosbnmn;
    FCELL *ww, *w,*wup,*wdown, *ll, *slo, *dist, *mxsl, *mnsl, *ang, *angr;       // weights (1/0 in circle, all 1 in square)

    FCELL zz, z0, c, tang,dcc,drr, angup, angdown, mediamx, mediamn, mediamx8, mediamn8, resol, sl, MDE, minslo,upmx,upmn,dnmx,dnmn,bnmn,bnmx;
    CELL i, j, h, hh, b,*tams;

    bool mdr;

    int pt1, pt2, *dxx, *dyy, *En8, *borde;

    struct Cell_head region;

    char *mapsetb, *mapsetd, *name_b, *name_d;
    char *name_8mx, *name_omx, *name_upmx, *name_dnmx, *name_bnmx;
    char *name_8mn, *name_omn, *name_upmn, *name_dnmn, *name_bnmn;

    int tipob, tipod,row,col,nn,vv,dc,dr,dx,dy,rc,rr,posts,tam, modo;
    int fdd,fdb,fdomx,fd8mx,fdupmx,fddnmx,fdbnmx,fdomn,fd8mn,fdupmn,fddnmn,fdbnmn;  //File identifier

    //double pi2=2*3.141592,a;

    G_gisinit(argv[0]);

    // Parameter specifications
    module = G_define_module();
    //module->keywords = _("raster, statistics");
    module->description =_("Calculates openness parameter in a DEM");
    
    G_add_keyword(_("raster"));
    G_add_keyword(_("focal"));
    G_add_keyword(_("geomorphometry"));

    mapab=G_define_standard_option(G_OPT_R_INPUT);
    mapab -> key = "elev";
    mapab -> description = "Elevation map"; 

    mapad=G_define_standard_option(G_OPT_R_INPUT);
    mapad -> key = "flowdir";
    mapad -> required = NO;
    mapad -> description = "Flow direction map";

    mapa8mx=G_define_standard_option(G_OPT_R_OUTPUT);
    mapa8mx-> key = "popenness8";
    mapa8mx -> required = NO;
    mapa8mx -> description = "8 directions positive openness map";

    mapaomx=G_define_standard_option(G_OPT_R_OUTPUT);
    mapaomx -> key = "popennesstot";
    mapaomx -> required = NO;
    mapaomx -> description = "All directions positive openness map";

    mapaupmx=G_define_standard_option(G_OPT_R_OUTPUT);
    mapaupmx -> key = "popennessup";
    mapaupmx -> required = NO;
    mapaupmx -> description = "Upslope directions positive openness map";

    mapadnmx=G_define_standard_option(G_OPT_R_OUTPUT);
    mapadnmx -> key = "popennessdn";
    mapadnmx -> required = NO;
    mapadnmx -> description = "Downslope directions positive openness map";

    mapabnmx=G_define_standard_option(G_OPT_R_OUTPUT);
    mapabnmx -> key = "popennessbn";
    mapabnmx -> required = NO;
    mapabnmx -> description = "Banks directions positive openness map";

    mapa8mn=G_define_standard_option(G_OPT_R_OUTPUT);
    mapa8mn-> key = "nopenness8";
    mapa8mn -> required = NO;
    mapa8mn -> description = "8 directions negative openness map";

    mapaomn=G_define_standard_option(G_OPT_R_OUTPUT);
    mapaomn -> key = "nopennesstot";
    mapaomn -> required = NO;
    mapaomn -> description = "All directions negative openness map";

    mapaupmn=G_define_standard_option(G_OPT_R_OUTPUT);
    mapaupmn -> key = "nopennessup";
    mapaupmn -> required = NO;
    mapaupmn -> description = "Upslope directions negative openness map";

    mapadnmn=G_define_standard_option(G_OPT_R_OUTPUT);
    mapadnmn -> key = "nopennessdn";
    mapadnmn -> required = NO;
    mapadnmn -> description = "Downslope directions negative openness map";

    mapabnmn=G_define_standard_option(G_OPT_R_OUTPUT);
    mapabnmn -> key = "nopennessbn";
    mapabnmn -> required = NO;
    mapabnmn -> description = "Banks directions negative openness map";

    ventana=G_define_option();
    ventana -> key = "size";
    ventana -> type = TYPE_INTEGER;
    ventana -> description = "window size";
    ventana -> required = YES;

    mode1=G_define_option();
    mode1 -> key = "mode";
    mode1 -> type = TYPE_INTEGER; 
    mode1 -> description = "1: square, 2: circle. If mode=2 size is the diameter"; 
    mode1 -> required = NO;
    mode1 -> answer = "1";

    if(G_parser(argc, argv)) exit(EXIT_FAILURE);

    name_b = mapab->answer;
    name_d = mapad->answer;
    name_8mx = mapa8mx->answer;
    name_omx = mapaomx->answer;
    name_upmx = mapaupmx->answer;
    name_dnmx = mapadnmx->answer;
    name_bnmx = mapabnmx->answer;
    name_8mn = mapa8mn->answer;
    name_omn = mapaomn->answer;
    name_upmn = mapaupmn->answer;
    name_dnmn = mapadnmn->answer;
    name_bnmn = mapabnmn->answer;

    if (name_omx == NULL && name_8mx == NULL && name_upmx == NULL && name_dnmx == NULL && name_bnmx == NULL && 
        name_omn == NULL && name_8mn == NULL && name_upmn == NULL && name_dnmn == NULL && name_bnmn == NULL){
	G_fatal_error(_("You must specify at least one of the parameters: "
			"<%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s>, <%s> or <%s>"),
                        mapa8mx->key, mapaomx->key, mapaupmx->key, mapadnmx->key, mapabnmx->key,
                        mapa8mn->key, mapaomn->key, mapaupmn->key, mapadnmn->key, mapabnmn->key);
    }
    if (name_b == NULL) G_fatal_error(_("You must specify an elevation map"));
    if ((name_upmx != NULL || name_dnmx != NULL || name_bnmx != NULL || name_upmn != NULL || name_dnmn != NULL || name_bnmn != NULL) && name_d==NULL)
        G_fatal_error(_("You must specify a flow direction map"));

    // Obtains basic info
    modo=atoi(mode1->answer);
    nn=atoi(ventana->answer);
    vv=(nn-1)/2;              
    tam=nn*nn;                 

    G_get_window(&region);

    // Reserves memory for values inside the window
    w = (FCELL *) G_malloc(tam* sizeof(FCELL));    
    ww = (FCELL *) G_malloc(tam* sizeof(FCELL));    
    wup = (FCELL *) G_malloc(tam* sizeof(FCELL));     
    wdown = (FCELL *) G_malloc(tam* sizeof(FCELL));      
    dxx =  (CELL *) G_malloc(tam* sizeof(CELL)); 
    dyy =  (CELL *) G_malloc(tam* sizeof(CELL)); 
    ll = (FCELL *) G_malloc(tam* sizeof(FCELL));   
    slo = (FCELL *) G_malloc(tam* sizeof(FCELL));  
    dist = (FCELL *) G_malloc(tam* sizeof(FCELL));   
    ang = (FCELL *) G_malloc(tam* sizeof(FCELL));   
    mxsl = (FCELL *) G_malloc(tam* sizeof(FCELL));  
    mnsl = (FCELL *) G_malloc(tam* sizeof(FCELL));  
    angr = (FCELL *) G_malloc(tam* sizeof(FCELL));   
    En8 = (CELL *) G_malloc(tam* sizeof(CELL));    
    borde = (CELL *) G_malloc(tam* sizeof(CELL));    
 
    // Open files and Reserve memory for input layers        
    fdb = Rast_open_old(name_b, "");
    tipob=Rast_get_map_type(fdb);
    
    if (name_d != NULL){
         fdd=Rast_open_old(name_d, "");		
         tipod=Rast_get_map_type(fdd);
         if (tipod!=0) G_fatal_error(_("%s=%s - must be CELL type"), mapad->key, name_d);
    }
    
    switch(tipob){
       case 0:{
          entrada_c=(CELL  *) G_malloc(region.rows*region.cols * sizeof(CELL));  
          datosb_c=(CELL  *) G_malloc(region.cols * sizeof(CELL));  break;
       }
       case 1:{
          entrada_f=(FCELL  *) G_malloc(region.rows*region.cols * sizeof(FCELL));
          datosb_f=(FCELL *) G_malloc(region.cols * sizeof(FCELL)); break;
       }
       case 2:{
          entrada_d=(DCELL  *) G_malloc(region.rows*region.cols * sizeof(DCELL));
          datosb_d=(DCELL *) G_malloc(region.cols * sizeof(DCELL)); break;
       }
    }
    if (name_d != NULL){
         datosd=(CELL *) G_malloc(region.cols * sizeof(CELL)); 
         entrada_dir=(CELL  *) G_malloc(region.rows*region.cols * sizeof(CELL));  
    }

    //Reserve memory for output layers
    if (name_omx != NULL) {
        fdomx=Rast_open_new(name_omx,1); 
        datosomx = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_8mx != NULL) {
        fd8mx=Rast_open_new(name_8mx,1); 
        datos8mx = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_upmx != NULL) {
        fdupmx=Rast_open_new(name_upmx,1); 
        datosupmx = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_dnmx != NULL) {
        fddnmx=Rast_open_new(name_dnmx,1); 
        datosdnmx = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_bnmx != NULL) {
        fdbnmx=Rast_open_new(name_bnmx,1); 
        datosbnmx = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_omn != NULL) {
        fdomn=Rast_open_new(name_omn,1); 
        datosomn  = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_8mn != NULL) {
        fd8mn=Rast_open_new(name_8mn,1); 
        datos8mn  = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_upmn != NULL) {
        fdupmn=Rast_open_new(name_upmn,1); 
        datosupmn  = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_dnmn != NULL) {
        fddnmn=Rast_open_new(name_dnmn,1); 
        datosdnmn  = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
    if (name_bnmn != NULL) {
        fdbnmn=Rast_open_new(name_bnmn,1); 
        datosbnmn  = (FCELL *) G_malloc(region.cols * sizeof(FCELL));
    }
 

    // Calculate weights matrix as a function of window size and mode (circle or square).
    // Border cells have 1 and the rest 0
    printf("Distances to center cell: \n");
    for (row=0; row<nn; row++){printf("   "); for (col=0; col<nn; col++){
        *(ll + row*nn +col) = sqrt((row-vv)*(row-vv)*region.ns_res*region.ns_res + (col-vv)*(col-vv)*region.ew_res*region.ew_res );
        if (modo==2){
           if (0.99 * sqrt((row-vv)*(row-vv) + (col-vv)*(col-vv)) <= vv && 0.99 * sqrt((row-vv)*(row-vv) + (col-vv)*(col-vv)) >= (vv-1)){
              *(ww + row*nn +col) = 1;
           } else {
              *(ww + row*nn +col) = 0;
           }
        } 
        if (modo==1){
            if ( (row==0) || (row==(nn-1)) || (col==0) || (col==(nn-1)) ) *(ww + row*nn +col) = 1; else *(ww + row*nn +col) = 0;
        }
        printf("%d ",(CELL) *(ll + row*nn +col));
    } printf("\n"); }

    printf("\nradius: %d\n",vv);
    for (row=-vv; row<=vv; row++){;for (col=-vv; col<=vv; col++){  
        if(row>0 && col==0) *(ang + (row+vv)*nn + col +vv) = 180;
        if(row==0 && col>0) *(ang + (row+vv)*nn + col +vv) = 90;
        if(row<0 && col==0) *(ang + (row+vv)*nn + col +vv) = 0;
        if(row==0 && col<0) *(ang + (row+vv)*nn + col +vv) = 270;
        if(row<0 && col>0) *(ang + (row+vv)*nn + col +vv) = 180*atan2(col,-row)/3.14159;
        if(row>0 && col>0) *(ang + (row+vv)*nn + col +vv) = 90 + 180*atan2(row,col)/3.14159;
        if(row>0 && col<0) *(ang + (row+vv)*nn + col +vv) = 180 + 180*atan2(-col,row)/3.14159;
        if(row<0 && col<0) *(ang + (row+vv)*nn + col +vv) = 270 + 180*atan2(-row,-col)/3.14159;
    } }

    printf("\nweights: \n");
    for (row=0; row<nn; row++){printf("   "); for (col=0; col<nn; col++){   
        printf("%d ",(CELL) *(ww + row*nn +col));
    } printf("\n"); }

    printf("\ncompass directions:\n");
    for (row=0; row<nn; row++){printf("   "); for (col=0; col<nn; col++){   
        printf("%3.3f ",  *(ang + row*nn +col));
    } printf("\n"); }
 
 
printf("Read input files: \n");
   // Read input files  
    for (row = 0; row < region.rows; row++) {
	    G_percent(row, region.rows - 1, 1);
        if(name_d!=NULL) Rast_get_row(fdd,datosd,row,tipod);
        switch(tipob){
           case 0: Rast_get_row(fdb,datosb_c,row,tipob);break;
           case 1: Rast_get_row(fdb,datosb_f,row,tipob);break;
           case 2: Rast_get_row(fdb,datosb_d,row,tipob);break;
        }
        for (col=0;col<region.cols;col++){ 
           if(name_d!=NULL) *(entrada_dir + row*region.cols + col)= abs(*(datosd + col));  
           switch(tipob){
              case 0: *(entrada_c + row*region.cols + col)= *(datosb_c + col);break;
              case 1: *(entrada_f + row*region.cols + col)= *(datosb_f + col);break;
              case 2: *(entrada_d + row*region.cols + col)= *(datosb_d + col);break;
           }
        }
    }

printf("Processing:\n");
    //Process
    for (row = 0; row < region.rows; row++) {
       G_percent(row, region.rows - 1, 1);
       for (col=0;col<region.cols;col++){  
          if (row>=vv && row<(region.rows-vv) && col>=vv && col<(region.cols-vv)){// && row==7 && col==7){

            // Get height of the central cell
            switch(tipob){
                case 0: z0= (FCELL) *(entrada_c + row*region.cols +col);break;
                case 1: z0= (FCELL) *(entrada_f + row*region.cols +col);break;
                case 2: z0= (FCELL) *(entrada_d + row*region.cols +col);break;
            } 

            // Compute slopes from the central cell for the whole window
            //printf("\nslopes:\n   ");
            for (dr=-vv;dr<=vv;dr++){ 
                //printf("   ");
                for (dc=-vv;dc<=vv;dc++){
                   switch(tipob){
                      case 0: zz= (FCELL) *(entrada_c + (row+dr)*region.cols + col+dc);break;
                      case 1: zz= (FCELL) *(entrada_f + (row+dr)*region.cols + col+dc);break;
                      case 2: zz= (FCELL) *(entrada_d + (row+dr)*region.cols + col+dc);break;
                   } 
                   *(slo + (dr+vv)*nn + (dc+vv)) = (zz - z0) / *(ll + (dr+vv)*nn + (dc+vv));
                   if (dr==0 && dc==0) *(slo + (dr+vv)*nn + (dc+vv)) = 0;             
                   //printf("   %d %d %f %f %f %f\n",dr,dc,z0,zz,*(ll + (dr+vv)*nn + (dc+vv)),*(slo + (dr+vv)*nn + (dc+vv)));
                  //printf ("%f ",*(slo + (dr+vv)*nn + (dc+vv)));
                } //printf("\n");
            } //printf("\n");

            // For each cell with w=1 (that is for each cell in the border)
            h=-1;
            for (dr=-vv;dr<=vv;dr++){
               for (dc=-vv;dc<=vv;dc++){
                  if (*(ww + (dr+vv)*nn + (dc+vv)) == 1){
                     h=h+1;
                     mxsl[h]=0;mnsl[h]=0;
                     if ( (abs(dr)==vv) && (abs(dc)==vv) || (abs(dr)==vv) && (dc==0 )|| (abs(dc)==vv) && (dr==0)) En8[h]=1; else En8[h]=0;
                     angr[h]=*(ang + (dr+vv)*nn + (dc+vv));
                     if (abs(dr)>=abs(dc)) mdr=TRUE; else mdr=FALSE;  
                     // Meridian and parallel lines of sight
                     if (dr==0 && dc>0) for (i=1;i<=dc;i++) {
                             if(mxsl[h]<*(slo + (dr+vv)*nn + (i+vv))) mxsl[h]=*(slo + (dr+vv)*nn + i + vv); 
                             if(mnsl[h]>*(slo + (dr+vv)*nn + (i+vv))) mnsl[h]=*(slo + (dr+vv)*nn + i + vv);
                        // printf("       dr=%d dc=%d   i=%d  slo=%f maxslo=%f  minslo=%f\n",dr,dc,i, j, rc, j*nn + rc + vv, *(slo + (dr+vv)*nn + i+vv),mxsl[h],mnsl[h]);
                     }
                     if (dr==0 && dc<0) for (i=-1;i>=dc;i--) {
                             if(mxsl[h]<*(slo + (dr+vv)*nn + (i+vv))) mxsl[h]=*(slo + (dr+vv)*nn + i + vv);
                             if(mnsl[h]>*(slo + (dr+vv)*nn + (i+vv))) mnsl[h]=*(slo + (dr+vv)*nn + i + vv);
                          //printf("       dr=%d dc=%d   i=%d  slo=%f maxslo=%f  minslo=%f\n",dr,dc,i, j, rc, j*nn + rc + vv, *(slo + (dr+vv)*nn + i+vv),mxsl[h],mnsl[h]);
                     }
                     if (dc==0 && dr>0) for (i=1;i<=dr;i++) {
                             if(mxsl[h]<*(slo + (i+vv)*nn + (dc+vv))) mxsl[h]=*(slo + (i+vv)*nn + dc + vv);
                             if(mnsl[h]>*(slo + (i+vv)*nn + (dc+vv))) mnsl[h]=*(slo + (i+vv)*nn + dc + vv);
                         // printf("       dr=%d dc=%d   i=%d  slo=%f maxslo=%f  minslo=%f\n",dr,dc,i, j, rc, j*nn + rc + vv, *(slo + (i+vv)*nn + dc+vv),mxsl[h],mnsl[h]);
                     }
                     if (dc==0 && dr<0) for (i=-1;i>=dr;i--) {
                             if(mxsl[h]<*(slo + (i+vv)*nn + (dc+vv))) mxsl[h]=*(slo + (i+vv)*nn + dc+vv);
                             if(mnsl[h]>*(slo + (i+vv)*nn + (dc+vv))) mnsl[h]=*(slo + (i+vv)*nn + dc+vv);
                          //printf("       dr=%d dc=%d   i=%d  slo=%f maxslo=%f  minslo=%f\n",dr,dc,i, j, rc, j*nn + rc + vv, *(slo + (i+vv)*nn + dc+vv),mxsl[h],mnsl[h]);
                     }
                     // The other lines of sight
                     if (dc>0 && dr<0 && mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=-vv;i<0;i++) {  
                                 j=round(i*dc/dr);  
                                 sl=*(slo + (i+vv)*nn + (j+vv));     
                                 if( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                 //printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);     
                           }
                     }
                     if (dc>0 && dr<0 && !mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=vv;i>0;i--) {  
                                 j=round(i*dr/dc);
                                 sl=*(slo + (j+vv)*nn + (i+vv));       
                                 if( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                 //printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);                    
                           }
                     }

                     if (dc>0 && dr>0 && mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=vv;i>0;i--) {  
                                 j=round(i*dc/dr); 
                                 sl=*(slo + (i+vv)*nn + (j+vv));    
                                 if( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                 //printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);                    
                           }
                     }
                     if (dc>0 && dr>0 && !mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=vv;i>0;i--) {  
                                 j=round(i*dr/dc); 
                                 sl=*(slo + (j+vv)*nn + (i+vv));    
                                 if( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                 //printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);                    
  }
                     }
                     if (dc<0 && dr>0 && !mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=-vv;i<0;i++) {  
                                 j=round(i*dr/dc);  
                                 sl=*(slo + (j+vv)*nn + (i+vv));  
                                 if ( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if ( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                 //printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);
                        }
                     }
                     if (dc<0 && dr>0 && mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=vv;i>0;i--) {  
                                 j=round(i*dc/dr);   
                                 sl=*(slo + (i+vv)*nn + (j+vv)); 
                                 if( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                 //printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);
                         }
                     }
                     if (dc<0 && dr<0 && mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=-vv;i<0;i++) {  
                                 j=round(i*dc/dr);
                                 sl=*(slo + (i+vv)*nn + (j+vv)); 
                                 if( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                 //printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);
                           }
                     }
                     if (dc<0 && dr<0 && !mdr) {
                           mxsl[h]=0;mnsl[h]=0;
                           for (i=-vv;i<0;i++) {  
                                 j=round(i*dr/dc);   
                                 sl=*(slo + (j+vv)*nn + (i+vv)); 
                                 if( (mxsl[h]<sl) && (sl>0) ) mxsl[h]=sl;
                                 if( (mnsl[h]>sl) && (sl<0) ) mnsl[h]=sl;
                                // printf("        dr=%d dc=%d   i=%d  j=%d    slo=%f maxslo=%f minslo=%f\n",dr,dc,i, j, sl,mxsl[h],mnsl[h]);
                          }
                     }

                     //printf("        h=%d %d %d    %f %f %f\n",h,dr,dc,mxsl[h],mnsl[h],(fabs(mnsl[h])));
                     mxsl[h]= 90 - 180*atan(mxsl[h])/PI;
                     mnsl[h]= 90 - 180*atan(fabs(mnsl[h]))/PI;

                     //printf("        h=%d %d %d    %f %f\n",h,dr,dc,mxsl[h],mnsl[h]);
                  }                       
            }}
//printf("gg\n");//a=40;

            if(name_d!=NULL & (name_omx != NULL || name_8mx != NULL || name_upmx != NULL || name_dnmx != NULL || name_bnmx != NULL || name_omn != NULL || name_8mn != NULL || name_upmn != NULL || name_dnmn != NULL || name_bnmn != NULL)){

                //Flow output angle
                 dr=row;dc=col;
                 while (dr<(row+vv)  && dr>(row-vv)  && dc<(col+vv)  && dc>(col-vv)){
                      asp2dxdy(*(entrada_dir + row*region.cols +col),&dx,&dy);
		      dr=dr+dy;dc=dc+dx;
                 }
                 angdown=dxdy2ang(dc-col,-(dr-row));
//printf("           %d %d    dc=%d dr=%d angdown=%f\n",col,row,dc-col,dr-row,angdown);   



                 //Flow input angle
                 pt1=0;pt2=-1;
                 *(dyy+pt1)=row;*(dxx+pt1)=col;*(slo+pt1)=999;*(dist+pt1)=0;
//printf(" pt1=%d pt2=%d \n",pt1,pt2);
                 while (pt1>pt2){
                     pt2=pt2+1;
//printf("pt1=%d  pt2=%d    dxx=%d dyy=%d\n",pt1,pt2,*(dxx+pt2), *(dyy+pt2));
                           // printf("%d ",pt2);
                     for (dx=-1;dx<=1;dx=dx+1){ for (dy=-1;dy<=1;dy=dy+1){
                         if ((dx!=0 || dy!=0) && (*(dyy+pt2)+dy)>=0 && (*(dxx+pt2)+dx)>=0 && (*(dyy+pt2)+dy)<region.rows && (*(dxx+pt2)+dx)<region.cols){
                             //printf("   %d %d   \n",dx,dy);
                            if (dx==0 | dy==0) resol=region.ns_res; else resol=sqrt(2*pow(region.ns_res,2));
                            //printf("   %d %d     %d %d %f ",*(dxx+pt2), *(dyy+pt2), dx,dy,resol);
                            //printf("   %d %d\n",dxdy2antiasp(dx,dy),*(entrada_dir + (*(dyy+pt2)+dy)*region.cols +*(dxx+pt2)+dx));
                            if(dxdy2antiasp(dx,dy)==*(entrada_dir + (*(dyy+pt2)+dy)*region.cols +*(dxx+pt2)+dx)){
                                if ((*(dyy+pt2)+dy)>=(row-vv) && (*(dyy+pt2)+dy)<=(row+vv) && (*(dxx+pt2)+dx)>=(col-vv) && (*(dxx+pt2)+dx)<=(col+vv)){
                                    pt1=pt1+1;
                                    *(dyy+pt1)=*(dyy+pt2)+dy ;*(dxx+pt1)=*(dxx+pt2)+dx;
                                    switch(tipob){
                                       case 0: MDE= (FCELL) *(entrada_c + *(dyy+pt1) * region.cols + *(dxx+pt1));break;
                                       case 1: MDE= (FCELL) *(entrada_f + *(dyy+pt1) * region.cols + *(dxx+pt1));break;
                                       case 2: MDE= (FCELL) *(entrada_d + *(dyy+pt1) * region.cols + *(dxx+pt1));break;
                                    } 
                                    *(dist+pt1)=*(dist+pt2)+resol;*(slo+pt1)=(MDE-z0) / *(dist+pt1);
                                    if (*(dyy+pt1)==row+vv || *(dyy+pt1)==row-vv || *(dxx+pt1)==col+vv || *(dxx+pt1)==col-vv) *(borde+pt1)=1; else *(borde+pt1)=0;
                                    //printf("            %d %d  %f %f \n",*(dxx+pt1), *(dyy+pt1), *(dist+pt1), *(slo+pt1));
                                } //else {printf("\n");}
                            } //else {printf("\n");}
                         }
                      }}
                 }
                 angup=9999;minslo=999;
//printf("gg   2 \n");
                 for (pt2=0;pt2<=pt1;pt2++){ 
                      //printf("X=%d Y=%d dx=%d dy=%d d=%f s=%f b=%d ang=%f\n",*(dxx+pt2), *(dyy+pt2),*(dxx+pt2)-col,-(*(dyy+pt2)-row), *(dist+pt2),*(slo+pt2),*(borde+pt2),dxdy2ang(*(dxx+pt2)-col,-(*(dyy+pt2)-row)));
                      if (*(borde+pt2)==1 && *(slo+pt2)<minslo){ minslo= *(slo+pt2); angup=dxdy2ang(*(dxx+pt2)-col,-(*(dyy+pt2)-row));}
                 }
                 if (angup==9999)  for (pt2=0;pt2<=pt1;pt2++){
                      if (*(slo+pt2)<minslo) {minslo=*(slo+pt2); angup=dxdy2ang(*(dxx+pt2)-col,-(*(dyy+pt2)-row));}
                 }
                 if (minslo==999) angup=180+angdown; if(angup>360) angup=angup-360;
                 //printf("r=%d c=%d z=%f minslo=%f  angup=%f angdown=%f\n",row,col,z0,minslo,angup,angdown);               
            }            

            //printf("angup=%f angdown=%f \n",angup,angdown);
            mediamx=0;mediamn=0;mediamx8=0;mediamn8=0;
            for (hh=0;hh<=h;hh++){
                //printf("               %d %d %f %f\n",hh,En8[hh], mxsl[hh],mnsl[hh]);
                mediamx=mediamx+mxsl[hh];mediamn=mediamn+mnsl[hh];
                if (En8[hh]==1) {
                    mediamx8=mediamx8+mxsl[hh];mediamn8=mediamn8+mnsl[hh];                    
                }
                //printf("%f+",mxsl[hh]);
            } 
            //printf("\n\n mediamx=%f mediamn=%f mediamx8=%f mediamn8=%f\n\n",mediamx,mediamn,mediamx8,mediamn8);
            //printf("\n\n mediamx=%f mediamn=%f mediamx8=%f mediamn8=%f\n\n",mediamx/(4*(nn-1)),mediamn/(4*(nn-1)),mediamx8/8,mediamn8/8);     
       
            ponderAngulos(angup, angdown, h, angr, wdown, wup, w);
            upmx=0;upmn=0;dnmx=0;dnmn=0;bnmn=0;bnmx=0;
            for (hh=0;hh<=h;hh++){
            // printf("      %d   a=%f mx=%f mn=%f    wu=%f wd=%f wb=%f\n",hh,*(angr+hh),*(mxsl+hh),*(mnsl+hh),*(wup+hh),*(wdown+hh),*(w+hh));
                   upmx = upmx + *(wup+hh) * mxsl[hh];
                   upmn = upmn + *(wup+hh) * mnsl[hh];
                   dnmx = dnmx + *(wdown+hh) * mxsl[hh];
                   dnmn = dnmn + *(wdown+hh) * mnsl[hh];
                   bnmx = bnmx + *(w+hh) * mxsl[hh];
                   bnmn = bnmn + *(w+hh) * mnsl[hh];
            }
          // printf("r=%d c=%d Smx=%f  Smn=%f  Smx8=%f Smn8=%f Smxup=%f  Smnup=%f Smxdn=%f Smndn=%f Smxbn=%f Smnbn=%f  \n\n",row,col,mediamx,mediamn,mediamx8,mediamn8,upmx,upmn,dnmx,dnmn,bnmx,bnmn);

            if (name_omx != NULL) *(datosomx + col) = mediamx / (4*(nn-1));
            if (name_8mx != NULL) *(datos8mx + col) = mediamx8 / 8;
            if (name_omn != NULL) *(datosomn + col) = mediamn / (4*(nn-1));
            if (name_8mn != NULL) *(datos8mn + col) = mediamn8 / 8;
            if (name_upmx != NULL) *(datosupmx + col) = upmx;
            if (name_dnmx != NULL) *(datosdnmx + col) = dnmx;
            if (name_bnmx != NULL) *(datosbnmx + col) = bnmx;
            if (name_upmn != NULL) *(datosupmn + col) = upmn;
            if (name_dnmn != NULL) *(datosdnmn + col) = dnmn;
            if (name_bnmn != NULL) *(datosbnmn + col) = bnmn;
          }else{ 
                if (name_omx != NULL) Rast_set_null_value (datosomx + col,1,1);
                if (name_omn != NULL) Rast_set_null_value (datosomn + col,1,1);
                if (name_8mx != NULL) Rast_set_null_value (datos8mx + col,1,1);
                if (name_8mn != NULL) Rast_set_null_value (datos8mn + col,1,1);
                if (name_upmx != NULL) Rast_set_null_value (datosupmx + col,1,1);
                if (name_dnmx != NULL) Rast_set_null_value (datosdnmx + col,1,1);
                if (name_bnmx != NULL) Rast_set_null_value (datosbnmx + col,1,1);
                if (name_upmn != NULL) Rast_set_null_value (datosupmn + col,1,1);
                if (name_dnmn != NULL) Rast_set_null_value (datosdnmn + col,1,1);
                if (name_bnmn != NULL) Rast_set_null_value (datosbnmn + col,1,1);
          }
       }
       if (name_omx != NULL) Rast_put_row(fdomx,datosomx,1);
       if (name_8mx != NULL) Rast_put_row(fd8mx,datos8mx,1);
       if (name_omn != NULL) Rast_put_row(fdomn,datosomn,1);
       if (name_8mn != NULL) Rast_put_row(fd8mn,datos8mn,1);
       if (name_upmx != NULL) Rast_put_row(fdupmx,datosupmx,1);
       if (name_dnmx != NULL) Rast_put_row(fddnmx,datosdnmx,1);
       if (name_bnmx != NULL) Rast_put_row(fdbnmx,datosbnmx,1);
       if (name_upmn != NULL) Rast_put_row(fdupmn,datosupmn,1);
       if (name_dnmn != NULL) Rast_put_row(fddnmn,datosdnmn,1);
       if (name_bnmn != NULL) Rast_put_row(fdbnmn,datosbnmn,1);
//printf("Cierro celda\n");
    }

    // CERRAR Y SALIR
printf("Closing\n");
    Rast_close(fdb);
    if (name_d != NULL){Rast_close(fdd); G_free(datosd);G_free(entrada_dir);}
    if (name_omx != NULL){Rast_close(fdomx);G_free(datosomx);}
    if (name_8mx != NULL){Rast_close(fd8mx);G_free(datos8mx);}
    if (name_upmx != NULL){Rast_close(fdupmx);G_free(datosupmx);}
    if (name_dnmx != NULL){Rast_close(fddnmx);G_free(datosdnmx);}
    if (name_bnmx != NULL){Rast_close(fdbnmx);G_free(datosbnmx);}
    if (name_omn != NULL){Rast_close(fdomn);G_free(datosomn);}
    if (name_8mn != NULL){Rast_close(fd8mn);G_free(datos8mn);}
    if (name_upmn != NULL){Rast_close(fdupmn);G_free(datosupmn);}
    if (name_dnmn != NULL){Rast_close(fddnmn);G_free(datosdnmn);}
    if (name_bnmn != NULL){Rast_close(fdbnmn);G_free(datosbnmn);}

    // Liberar memoria de los punteros de ventana
    switch(tipob){
           case 0:  G_free(datosb_c); G_free(entrada_c); break;
           case 1:  G_free(datosb_f); G_free(entrada_f); break;
           case 2:  G_free(datosb_d); G_free(entrada_d); break;
    }
    return (0);
}
