#include "athena_preprocessor.fh"
#ifdef INTERACTIVE_GRAPHICS
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>

void set_drawing_state(int lgraph, int lmovie);
#define NWINDOWS 1

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

#define DEBUG
#undef FILEOUT
#undef MARGINS
#define MYRGBMODE
#define DOUBLEBUFFER

#undef DUMP

#ifdef DUMP
#define MYRGBMODE 1		/* s. maher - for gltotiff */
#endif

#ifndef MYRGBMODE
#ifndef DOUBLEBUFFER
#else
#endif
#else /*MYRGBMODE*/
#endif /*MYRGBMODE*/

#	define NBAND	64
#	define IBASE	0
#	define BLACK	IBASE
#       define GRAY     1
#       define  GREEN   2
#       define  BLUE    3
#       define  YELLOW  4
#       define  CYAN    5
#       define  RED     6
#       define WHITE    7
#	define IBOT	(IBASE+8)
#	define ITOP	(IBASE+NBAND*4-1)


#undef WIDTH 
#undef HEIGHT
#define WIDTH  800
#define HEIGHT 800
#define NPIX (WIDTH*HEIGHT)


/*
 * Routine to convert from hue to red, green, blue.
 * Based upon Foley, Van Dam, Feiner & Hughes, p. 593 (fig. 13.34).
 * Units: h = degrees, everything else is 0..1.
 * Modified to remove first sextant and replace it with diminishing saturation,
 * and to add an extra sextant with diminishing value.
 */
/***********************************************************************/
/*                                                                     */
/*                         H u e T o R G B                             */
/*                                                                     */
/***********************************************************************/
void HueToRGB(float h, float *r, float *g, float *b) {
  int ih;    /* Integer part of hue */
  float fh;  /* Fractional part of hue */
  
  h = (360.-h) * 7. / 360.;
  ih = floor(h);
  fh = h - (float)ih;
  switch (ih % 7) {
  case 0: *r = 1.0;
    *g = 1.0-fh;
    *b = 1.0-fh;
    break;
  case 1: *r = 1.0;
    *g = fh;
    *b = 0.0;
    break;
  case 2: *r = 1.0-fh;
    *g = 1.0;
    *b = 0.0;
    break;
  case 3: *r = 0.0;
    *g = 1.0;
    *b = fh;
    break;
  case 4: *r = 0.0;
    *g = 1.0-fh;
    *b = 1.0;
    break;
  case 5: *r = fh;
    *g = 0.0;
    *b = 1.0;
    break;
  case 6: *r = 1.0-fh;
    *g = 0.0;
    *b = 1.0-fh;
    break;
  }
}


/***********************************************************************/
/*                                                                     */
/*                          q u k f i n i s h _                        */
/*                                                                     */
/***********************************************************************/
void qukfinish_(void)
{

}

#define NCOLORS 16000
static float my_color_table[NCOLORS][3];

void my_glutSetColor(int i, float r, float g, float b)
{
  if (i < 0 || i >= NCOLORS) return;
  my_color_table[i][0] = r;
  my_color_table[i][1] = g;
  my_color_table[i][2] = b;
}

void my_glutGetColor(int i, float *r, float *g, float *b)
{
  if (i < 0 || i >= NCOLORS) return;
  *r = my_color_table[i][0];
  *g = my_color_table[i][1];
  *b = my_color_table[i][2];
}

void my_glClearIndex(float xi)
{
  int i;
  float r, g, b;
  i = xi;
  my_glutGetColor(i, &r, &g, &b);
  glClearColor(r, g, b, 0.0);
}

void my_glIndexf(float xi)
{
  int i;
  float r, g, b;
  i = xi;
  my_glutGetColor(i, &r, &g, &b);
  glColor3f(r, g, b);
}
 
  

/***********************************************************************/
/*                                                                     */
/*                           q u k i n i t _                           */
/*                                                                     */
/***********************************************************************/
void qukinit_(void)
{
  int i;
  float j;

#ifdef DEBUG
  fprintf(stderr,"qukinit_\n");
#endif

#ifndef MYRGBMODE
  for (i = 0; i < NBAND; i++) {
    j = ((float) i)/((float) (NBAND-1));
    if (j > 1.0) j = 1.0;
    my_glutSetColor(IBASE+i,         0.0, j,     1.0);
    my_glutSetColor(IBASE+i+NBAND,   0.0, 1.0,   1.0-j);
    my_glutSetColor(IBASE+i+2*NBAND, j,   1.0,   0.0);
    my_glutSetColor(IBASE+i+3*NBAND, 1.0, 1.0-j, 0.0);
  }

  /* One color is taken out of the first band for black */
  my_glutSetColor(BLACK, 0.0, 0.0, 0.0);
  /*my_glutSetColor(1, 1.0, 0.0, 0.0); */
  my_glutSetColor(2, 0.0, 1.0, 0.0);
  my_glutSetColor(3, 0.0, 0.0, 1.0);
  my_glutSetColor(4, 1.0, 1.0, 0.0);
  my_glutSetColor(5, 0.0, 1.0, 1.0);
  my_glutSetColor(6, 1.0, 0.0, 0.0);
  my_glutSetColor(WHITE, 1.0, 1.0, 1.0);
  my_glutSetColor(GRAY, 0.5, 0.5, 0.5);

  my_glutSetColor(ITOP, 1.0, 0.0, 0.5);
  my_glClearIndex(7.0);
#else /*MYRGBMODE*/
  glClearColor(0.0, 0.0, 0.0, 0.0);
#endif /*MYRGBMODE*/

  glClear(GL_COLOR_BUFFER_BIT);

}


void impl_colors(int i,float *rv,float  *gv,float  *bv)
{
    float j;

    if(i>=0 && i<=7) {
        switch(i) {
           case 0:
               *rv=0.; *gv=0.; *bv=0.; return;
               break;
           case 1:
               *rv=0.5; *gv=0.5; *bv=0.5; return;
               break;
           case 2:
               *rv=0.; *gv=1.; *bv=0.; return;
               break;
           case 3:
               *rv=0.; *gv=0.; *bv=1.; return;
               break;
           case 4:
               *rv=1.; *gv=1.; *bv=0.; return;
               break;
           case 5:
               *rv=0.; *gv=1.; *bv=1.; return;
               break;
           case 6:
               *rv=1.; *gv=0.; *bv=1.; return;
               break;
           case 7:
               *rv=1.; *gv=1.; *bv=1.; return;
               break;
        }
    } else
    if(i>=IBOT && i<IBASE+NBAND) {
         j = (float) (i-IBASE)/((float) (NBAND-1));
         if (j > 1.0) j = 1.0;
         *rv=0.0; *gv=j; *bv=1.0; return;
    } else
    if(i>=IBASE+NBAND && i<IBASE+2*NBAND) {
         j = (float) (i-IBASE-NBAND)/((float) (NBAND-1));
         if (j > 1.0) j = 1.0;
         *rv=0.0; *gv=1.0; *bv=1.0-j; return;
    } else
    if(i>=IBASE+2*NBAND && i<IBASE+3*NBAND) {
         j = (float) (i-IBASE-2*NBAND)/((float) (NBAND-1));
         if (j > 1.0) j = 1.0;
         *rv=j; *gv=1.0; *bv=0.0; return;
    } else
    if(i>=IBASE+3*NBAND && i<IBASE+4*NBAND-1) {
         j = (float) (i-IBASE-3*NBAND)/((float) (NBAND-1));
         if (j > 1.0) j = 1.0;
         *rv=1.0; *gv=1.0-j; *bv=0.0; return;
    } else
    if(i==ITOP) {
         *rv=1.0; *gv=0.0; *bv=0.5; return;
    }
}


/***********************************************************************/
/*                                                                     */
/*                         q u k f r a m e 1 _                         */
/*                                                                     */
/***********************************************************************/
void qukframe1_(
		float *xmin,float *ymin,float *zmin,
		float *xmax,float *ymax,float *zmax,
		float *xtrans,float *ytrans,float *ztrans,
		float *xscale,float *yscale,float *zscale,
		float *anglex, float *angley, float *anglez, 
		float *xaxis ,float *yaxis, float *zaxis)
{
  /* glClear(GL_COLOR_BUFFER_BIT); */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(*xmin, *xmax, *ymin, *ymax, *zmin, *zmax);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);  /* GL_LESS GL_EQUAL GL_LEQUAL GL_GREATER GL_GEQUAL */

  /* Model coordinates is the real space coordinates for the object. */
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glRotatef(*anglex,1.0,0.0,0.0);
  glRotatef(*angley,0.0,1.0,0.0);
  glRotatef(*anglez,0.0,0.0,1.0);
  /* glRotatef(*angle,*xaxis,*yaxis,*zaxis); */
  glScalef(    *xscale,*yscale,*zscale);
  glTranslatef(*xtrans,*ytrans,*ztrans);


}


/***********************************************************************/
/*                                                                     */
/*                         q u k f r a m e 3 _                         */
/*                                                                     */
/***********************************************************************/
void qukframe3_(
		int *lnblocks,int nodetype[],
		float bbox[][6])
{
  int i,j,nelema,c;
  float x1,x2,x3,x4;
  float y1,y2,y3,y4;
  float rv, gv, bv;
  int mvar,mx,my;
  int nvar,nx,ny;

  float vector[2];

#ifdef DEBUG
  fprintf(stderr,"qukframe3_\n");
#endif
  nelema = *lnblocks;

  /* Draw bounding boxes for the blocks. */

  for(j=0;j<nelema;j++)
    { 
      if( nodetype[j] == 1 && nelema > 1)  /* leaf blocks only */
        {
	  register float f0,x0,y0, delx,dely;

	  x0   = bbox[j][0]; delx = (bbox[j][1] - bbox[j][0]);
	  y0   = bbox[j][2]; dely = (bbox[j][3] - bbox[j][2]);

          c = GRAY;
          impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);
          
	  glBegin(GL_LINE_LOOP);
	     glVertex2f(x0     ,y0);
	     glVertex2f(x0+delx,y0);
	     glVertex2f(x0+delx,y0+dely);
	     glVertex2f(x0     ,y0+dely);
	  glEnd();
	}
    }
}

void qukframe3_3d_(
           int *ijkplanep,int *iplanep,int *jplanep,int *kplanep,
           int *nv1,int *nv2,int *nv3,
           float *cutplanex,float *cutplaney,float *cutplanez,
           float *fac_n, float bbox[])
{
  /* Draw bounding boxes for the blocks. */
    float x1,x2,x3,x4;
    float y1,y2,y3,y4;
    float z1,z2,z3,z4;
    int c,ix,iy,iz;
    float rv, gv, bv, fact;
    float x0,y0,z0, delx,dely,delz;

    int nx, ny, nz;
    int ijkplane,iplane,jplane,kplane;

    nx = *nv1; ny = *nv2; nz = *nv3;
    ijkplane=*ijkplanep; iplane=*iplanep;
    jplane=*jplanep; kplane=*kplanep; fact= *fac_n;

    x0   = bbox[0]; delx = (bbox[1] - bbox[0])/nx;
    y0   = bbox[2]; dely = (bbox[3] - bbox[2])/ny;
    z0   = bbox[4]; delz = (bbox[5] - bbox[4])/nz;

    c = GRAY;
    impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);

    switch(ijkplane) {
   
        case 1:
             x1= *cutplanex * fact;x2=x1;x3=x1;x4=x1;

	          y1 = bbox[2]; y2 = bbox[3];
	          y3 = bbox[3]; y4 = bbox[2];
	      
	          z1 = bbox[4]; z2 = bbox[4];
	          z3 = bbox[5]; z4 = bbox[5];

	          glBegin(GL_LINE_LOOP);
	             glVertex3f(x1,y1,z1);
	             glVertex3f(x2,y2,z2);
	             glVertex3f(x3,y3,z3);
	             glVertex3f(x4,y4,z4);
	          glEnd();
             break;
   
        case 2:
             y1= *cutplaney * fact;y2=y1;y3=y1;y4=y1;

	          x1 = bbox[0]; x2 = bbox[1];
	          x3 = bbox[1]; x4 = bbox[0];
	    
	          z1 = bbox[4]; z2 = bbox[4];
	          z3 = bbox[5]; z4 = bbox[5];

	          glBegin(GL_LINE_LOOP);
	             glVertex3f(x1,y1,z1);
	             glVertex3f(x2,y2,z2);
	             glVertex3f(x3,y3,z3);
	             glVertex3f(x4,y4,z4);
	          glEnd();
             break;
   
        case 3:
             z1= *cutplanez * fact;z2=z1;z3=z1;z4=z1;

	          x1 = bbox[0]; x2 = bbox[1];
	          x3 = bbox[1]; x4 = bbox[0];
	      
	          y1 = bbox[2]; y2 = bbox[2];
	          y3 = bbox[3]; y4 = bbox[3];

	          glBegin(GL_LINE_LOOP);
	             glVertex3f(x1,y1,z1);
	             glVertex3f(x2,y2,z2);
	             glVertex3f(x3,y3,z3);
	             glVertex3f(x4,y4,z4);
	          glEnd();
             break;
    }
}

/***********************************************************************/
/*                                                                     */
/*                         q u k f r a m e 5 _                         */
/*                                                                     */
/***********************************************************************/
void qukframe5_(float bbox[6])
{
  float x0,y0, delx,dely;
  int c;
  float rv, gv, bv;

          x0   = bbox[0];
          delx = (bbox[1] - bbox[0]);
          y0   = bbox[2];
          dely = (bbox[3] - bbox[2]);

          c = BLACK;
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
          
          glBegin(GL_LINE_LOOP);
             glVertex2f(x0     ,y0);
             glVertex2f(x0+delx,y0);
             glVertex2f(x0+delx,y0+dely);
             glVertex2f(x0     ,y0+dely);
          glEnd();
}


/***********************************************************************/
/*                                                                     */
/*                         q u k f r a m e 6 _                         */
/*                                                                     */
/***********************************************************************/
void qukframe6_( int *nxb,int *nyb,int *nzb,
                 float bbox[6], int ibmask[])
{
  int ix,iy,iz;
  float x0,y0, delx,dely;
  float x,y;
  int n, nx, ny, nz;
  int c;
  float rv, gv, bv;

  nx=*nxb; ny=*nyb; nz=*nzb;
  x0   = bbox[0]; delx = (bbox[1] - bbox[0])/nx;
  y0   = bbox[2]; dely = (bbox[3] - bbox[2])/ny;

  for(iy=0;iy<ny;iy++) {
     for(ix=0;ix<nx;ix++) {
          if(ibmask[iy*nx+ix] == 0) { 

              c = BLACK;
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
          
              glBegin(GL_POLYGON);
                 x=x0+ix*delx; y=y0+iy*dely;
                 glVertex2f(x     ,y);
                 glVertex2f(x+delx,y);
                 glVertex2f(x+delx,y+dely);
                 glVertex2f(x     ,y+dely);
              glEnd();

          }
     }
  } 
}

void qukframe6_3d_( int *nxb,int *nyb,int *nzb,
                 float bbox[6], int ibmask[])
{
  int ix,iy,iz;
  float x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4;
  float x5,y5,z5, x6,y6,z6, x7,y7,z7, x8,y8,z8;
  float delx, dely, delz;
  int nx, ny, nz;
  int c;
  float rv, gv, bv;

  nx=*nxb; ny=*nyb; nz=*nzb;
  delx = (bbox[1] - bbox[0])/nx;
  dely = (bbox[3] - bbox[2])/ny;
  delz = (bbox[5] - bbox[4])/nz;
  for(iz=0;iz<nz;iz++) {
     for(iy=0;iy<ny;iy++) {
         for(ix=0;ix<nx;ix++) {
            if(ibmask[iz*ny*nx+iy*nx+ix] == 0) { 
 
                x1=bbox[0]+(ix+0)*delx;    x4=x1;x5=x1;x8=x1;
                x2=bbox[0]+(ix+1)*delx;    x3=x2;x6=x2;x7=x2;

                y1=bbox[2]+(iy+0)*dely;    y2=y1;y5=y1;y6=y1;
                y4=bbox[2]+(iy+1)*dely;    y3=y4;y7=y4;y8=y4;

                z1=bbox[4]+(iz+0)*delz;    z2=z1;z3=z1;z4=z1;
                z5=bbox[4]+(iz+1)*delz;    z6=z5;z7=z5;z8=z5;

          
                /* paint 6 faces */
                /* xfaces */
                rv=0.2;gv=0.2;bv=0.2;glColor3f(rv,gv,bv);
                glBegin(GL_POLYGON);
                   glVertex3f(x1 ,y1, z1); glVertex3f(x5 ,y5, z5);
                   glVertex3f(x8 ,y8, z8); glVertex3f(x4 ,y4, z4);
                glEnd();
                glBegin(GL_POLYGON);
                   glVertex3f(x2 ,y2, z2); glVertex3f(x3 ,y3, z3);
                   glVertex3f(x7 ,y7, z7); glVertex3f(x6 ,y6, z6);
                glEnd();

                /* yfaces */
                rv=0.4;gv=0.4;bv=0.4;glColor3f(rv,gv,bv);
                glBegin(GL_POLYGON);
                   glVertex3f(x1 ,y1, z1); glVertex3f(x2 ,y2, z2);
                   glVertex3f(x6 ,y6, z6); glVertex3f(x5 ,y5, z5);
                glEnd();
                glBegin(GL_POLYGON);
                   glVertex3f(x3 ,y3, z3); glVertex3f(x4 ,y4, z4);
                   glVertex3f(x8 ,y8, z8); glVertex3f(x7 ,y7, z7);
                glEnd();

                /* zfaces */
                rv=0.6;gv=0.6;bv=0.6;glColor3f(rv,gv,bv);
                glBegin(GL_POLYGON);
                   glVertex3f(x4 ,y4, z4); glVertex3f(x3 ,y3, z3);
                   glVertex3f(x2 ,y2, z2); glVertex3f(x1 ,y1, z1);
                glEnd();
                glBegin(GL_POLYGON);
                   glVertex3f(x5 ,y5, z5); glVertex3f(x6 ,y6, z6);
                   glVertex3f(x7 ,y7, z7); glVertex3f(x8 ,y8, z8);
                glEnd();
            }
         }
     }
  }
}

void qukframe7_( )
{
  int c;
  float rv, gv, bv;
/* plot cross hair */
          c = RED;
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
          
	      glBegin(GL_LINE_LOOP);
	        glVertex2f(-100,0.);
	        glVertex2f( 100,0.);
	      glEnd();
	      glBegin(GL_LINE_LOOP);
	        glVertex2f(0.     ,-100);
	        glVertex2f(0.     , 100);
	      glEnd();
}

/***********************************************************************/
/*                                                                     */
/*                 q u k _ f i e l d l i n e s _ 2 d _                 */
/*                                                                     */
/***********************************************************************/
void quk_fieldlines_2d_( int *npts, float xf[], float yf[], int *sign)
{
  int i,j,nelema;
  /*   Colorindex c; */
  int c;
  float rv, gv, bv;

  nelema = *npts;

/*  printf("nelema=%d\n",nelema); */

  /* Draw fieldlines. */
          c=0;
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
          
  glBegin(GL_LINES);
    for(j=0;j<nelema-1;j++)
      { 
	  register float x0,y0,x1,y1;

	  x0   = xf[j];
	  y0   = yf[j];
	  x1   = xf[j+1];
	  y1   = yf[j+1];

	  glVertex2f(x0     ,y0);
	  glVertex2f(x1     ,y1);
    }
  glEnd();

}

/***********************************************************************/
/*                                                                     */
/*                         q u k f r a m e 4 _                         */
/*                                                                     */
/***********************************************************************/
void qukframe4_(int *i2, int *lgraph, int *lmovie)
{
  int ik;
  char fname[100];
  static int icount = 0;

  /* Not the best solution, but, it'll do. */
  char *titles[NWINDOWS] = {"athena"};
  set_drawing_state(*lgraph, *lmovie);
#ifdef DEBUG
  fprintf(stderr,"qukframe4_.0\n");
  fprintf(stderr,"*i2=%d\n", *i2);
#endif
  ik=*i2 - 1;
#ifdef DEBUG
  fprintf(stderr,"qukframe4_.1\n");
#endif
  set_window(ik);
#ifdef DEBUG
  fprintf(stderr,"qukframe4_.2\n");
#endif
  /*   glutSwapBuffers(); */  /* Off screen hates this. */
  glFlush();
#ifdef DEBUG
  fprintf(stderr,"qukframe4_.3\n");
#endif
  sprintf(fname, "%s_%4.4d", titles[ik], icount);
  if (*lgraph) icount++;
#ifdef DEBUG
  fprintf(stderr,"qukframe4_.4\n");
#endif
  save_screen(fname);
#ifdef DEBUG
  fprintf(stderr,"qukframe4_.5\n");
#endif

}

/***********************************************************************/
/*                                                                     */
/*                         q u k f r a m e 2 _                         */
/*                                                                     */
/***********************************************************************/
void qukframe2_(
	     int *lnblocks,int nodetype[],
	     float bbox[][4],
	     int *Pmvar,int *Pmx,int *Pmy,
	     int *Pnvar,int *Pnx,int *Pny,
	     float f[],
	     float* fmin,float* fmax)
{
  int i,j,nelema;
  float x1,x2,x3,x4;
  float y1,y2,y3,y4;
  float f1,f2,f3,f4;
  float cc,bb;
  int c;
  float rv, gv, bv;
  int mvar,mx,my;
  int nvar,nx,ny;

#ifdef MARGINS
  float xmarg, ymarg;
#endif

  mvar = *Pmvar;
  mx = *Pmx;
  my = *Pmy;

  nvar = *Pnvar;
  nx = *Pnx;
  ny = *Pny;


#ifdef NEVER
   {
     GLint mm;
     glGetIntegerv(GL_MATRIX_MODE, &mm);
     glMatrixMode(GL_PROJECTION);
     glLoadIdentity();
     gluOrtho2D(xmin, xmax,  ymin, ymax);
     glMatrixMode(mm);
   }
#endif

#ifdef NEVER
#ifdef MARGINS
  xmarg = (*xmax - *xmin)/10.0;
  ymarg = (*ymax - *ymin)/10.0;
  glOrtho(*xmin-xmarg, *xmax+xmarg, *ymin-ymarg, *ymax+ymarg, -1.0, 1.0);
#else
  glOrtho(*xmin, *xmax, *ymin, *ymax, -1.0, 1.0);
#endif
#endif

  nelema = *lnblocks;

  bb = (ITOP-IBOT)/((*fmax - *fmin)+1.0e-30);
  cc = IBOT - bb*(*fmin);

#ifdef DEBUG
  printf("nelema=%d\n",nelema);
  printf("mvar,mx,my,nvar,nx,ny=%d %d %d %d %d %d \n",
	 mvar,mx,my,nvar,nx,ny);
  printf("fmin,fmax= %f %f \n",*fmin,*fmax);
#endif

  for(j=0;j<nelema;j++)
    { 
      if( nodetype[j] == 1)  /* leaf blocks only */
        {
	  register float f0,x0,y0, delx,dely;
	  register int ix,iy;
	  float fv,rv,gv,bv,hv;

#ifdef NEVER
	  x0 = coord[j][0]-0.5*size[j][0];
	  y0 = coord[j][1]-0.5*size[j][1];

	  delx = size[j][0]/nx;
	  dely = size[j][1]/ny;
#else
	  x0   = bbox[j][0];
	  delx = (bbox[j][1] - bbox[j][0])/nx;
	  y0   = bbox[j][2];
	  dely = (bbox[j][3] - bbox[j][2])/ny;
#endif


	  for( iy=0; iy<ny; iy++) 
	    for( ix=0; ix<nx; ix++) 
	    {

	      x1 = x0+delx*(ix+0);
	      x2 = x0+delx*(ix+1);
	      x3 = x0+delx*(ix+1);
	      x4 = x0+delx*(ix+0);
	      
	      y1 = y0+dely*(iy+0);
	      y2 = y0+dely*(iy+0);
	      y3 = y0+dely*(iy+1);
	      y4 = y0+dely*(iy+1);
	      
	      f0 = f[0+mvar*(ix+mx*(iy+my*(j)))];
	      
	      f1 = f0;
	      f2 = f0;
	      f3 = f0;
	      f4 = f0;

	      glBegin(GL_POLYGON);

	      c = max(IBOT,min(ITOP,(int)(cc+bb*f1)));
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
	      glVertex2f(x1,y1);
          
	  c = max(IBOT,min(ITOP,(int)(cc+bb*f2)));
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
	  glVertex2f(x2,y2);
	      
	  c = max(IBOT,min(ITOP,(int)(cc+bb*f3)));
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
	  glVertex2f(x3,y3);
          
	  c = max(IBOT,min(ITOP,(int)(cc+bb*f4)));
          impl_colors(c, &rv, &gv, &bv);
          glColor3f(rv,gv,bv);
	  glVertex2f(x4,y4);
          
	      glEnd();
	    }
	}
    }
  
  return;
}

void qukframe2_3d_(
	     int *ijkplanel,int *iplanel,int *jplanel,int *kplanel,
             float *cutplanex,float *cutplaney,float *cutplanez,
	     float bbox[6],
	     int *nv1,int *nv2,int *nv3,
	     float f[],
	     float* fmin,float* fmax)
{
    float x1,x2,x3,x4;
    float y1,y2,y3,y4;
    float z1,z2,z3,z4;
    float f0,f1,f2,f3,f4;
    float cc,bb;
    int c,ix,iy,iz, i0;
    float rv, gv, bv;
    float x0,y0,z0, delx,dely,delz;

    int nx, ny, nz;
    int ijkplane,iplane,jplane,kplane;

    nx = *nv1; ny = *nv2; nz = *nv3;
    ijkplane=*ijkplanel; iplane=*iplanel;
    jplane=*jplanel; kplane=*kplanel;
  
    bb = (ITOP-IBOT)/((*fmax - *fmin)+1.0e-30);
    cc = IBOT - bb*(*fmin);

    x0   = bbox[0]; delx = (bbox[1] - bbox[0])/nx;
    y0   = bbox[2]; dely = (bbox[3] - bbox[2])/ny;
    z0   = bbox[4]; delz = (bbox[5] - bbox[4])/nz;


    switch(ijkplane) {
   
        case 1:
             x1= *cutplanex; x2=x1;x3=x1;x4=x1;
	     for( iz=0; iz<nz; iz++) { 
	       for( iy=0; iy<ny; iy++) {

	          y1 = y0+dely*(iy+0); y2 = y0+dely*(iy+1);
	          y3 = y0+dely*(iy+1); y4 = y0+dely*(iy+0);
	      
	          z1 = z0+delz*(iz+0); z2 = z0+delz*(iz+0);
	          z3 = z0+delz*(iz+1); z4 = z0+delz*(iz+1);
	          i0 = iy+iz*ny; 
                  f0=f[i0];
	          c = max(IBOT,min(ITOP,(int)(cc+bb*f0)));
                  impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);

	          glBegin(GL_POLYGON);
	             glVertex3f(x1,y1,z1);
	             glVertex3f(x2,y2,z2);
	             glVertex3f(x3,y3,z3);
	             glVertex3f(x4,y4,z4);
	          glEnd();
               }
             }
             break;
   
        case 2:
             y1= *cutplaney;y2=y1;y3=y1;y4=y1;
	     for( iz=0; iz<nz; iz++) { 
	       for( ix=0; ix<nx; ix++) {

	          x1 = x0+delx*(ix+0); x2 = x0+delx*(ix+1);
	          x3 = x0+delx*(ix+1); x4 = x0+delx*(ix+0);
	      
	          z1 = z0+delz*(iz+0); z2 = z0+delz*(iz+0);
	          z3 = z0+delz*(iz+1); z4 = z0+delz*(iz+1);
	          i0 = ix+iz*nx; 
                  f0=f[i0];
	          c = max(IBOT,min(ITOP,(int)(cc+bb*f0)));
                  impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);

	          glBegin(GL_POLYGON);
	             glVertex3f(x1,y1,z1);
	             glVertex3f(x2,y2,z2);
	             glVertex3f(x3,y3,z3);
	             glVertex3f(x4,y4,z4);
	          glEnd();
               }
             }
             break;
   
        case 3:
             z1= *cutplanez;z2=z1;z3=z1;z4=z1;
	     for( iy=0; iy<ny; iy++) { 
	       for( ix=0; ix<nx; ix++) {

	          x1 = x0+delx*(ix+0); x2 = x0+delx*(ix+1);
	          x3 = x0+delx*(ix+1); x4 = x0+delx*(ix+0);
	      
	          y1 = y0+dely*(iy+0); y2 = y0+dely*(iy+0);
	          y3 = y0+dely*(iy+1); y4 = y0+dely*(iy+1);
	          i0 = ix+iy*nx; 
                  f0=f[i0];
	          c = max(IBOT,min(ITOP,(int)(cc+bb*f0)));
                  impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);

	          glBegin(GL_POLYGON);
	             glVertex3f(x1,y1,z1);
	             glVertex3f(x2,y2,z2);
	             glVertex3f(x3,y3,z3);
	             glVertex3f(x4,y4,z4);
	          glEnd();
               }
             }
             break;
    }
}

void qukframe8_(
		float *x1p,float *y1p,float *z1p,
		float *x2p,float *y2p,float *z2p)
{
    float x1,x2,y1,y2,z1,z2;
    float rv, gv, bv;
    int c;
     

    x1=*x1p; x2=*x2p; y1=*y1p; y2=*y2p; z1=*z1p; z2=*z2p;
          

    glBegin(GL_LINES);

    c = RED;
    impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);
       glVertex3f(x1,y1,z1);
       glVertex3f(x2,y1,z1);
    c = GREEN;
    impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);
       glVertex3f(x1,y1,z1);
       glVertex3f(x1,y2,z1);
    c = BLUE;
    impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);
       glVertex3f(x1,y1,z1);
       glVertex3f(x1,y1,z2);

    c = BLACK;
    impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);

       glVertex3f(x2,y1,z1);
       glVertex3f(x2,y2,z1);
       glVertex3f(x2,y1,z1);
       glVertex3f(x2,y1,z2);

    glEnd();

    glBegin(GL_LINE_LOOP);
    c = BLACK;
    impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);
       glVertex3f(x2,y2,z1);
       glVertex3f(x1,y2,z1);
       glVertex3f(x1,y2,z2);
       glVertex3f(x2,y2,z2);
    glEnd();

    glBegin(GL_LINE_LOOP);
    c = BLACK;
    impl_colors(c, &rv, &gv, &bv); glColor3f(rv,gv,bv);
       glVertex3f(x1,y1,z2);
       glVertex3f(x2,y1,z2);
       glVertex3f(x2,y2,z2);
       glVertex3f(x1,y2,z2);
    glEnd();
}
#endif
