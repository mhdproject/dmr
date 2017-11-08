#include "athena_preprocessor.fh"
#ifdef INTERACTIVE_GRAPHICS

/* The off screen graphic code. */
/*
 * Objective: 
 * Render the graphics to a specified window size, that is
 * directly drawn to the screen.  The images are from predetermined
 * files which are created by a parent process.
 */

/* By John Merritt - Emergent Information Technologies */
/* Software License: GPL */
/* Copyright (C) 2001 */

#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <GL/osmesa.h>
#include <signal.h>
#include <pthread.h>

void write_gif(char *fname);

#define NWINDOWS 1
int  my_onwindow[NWINDOWS];  /* Global onscreen for these routines. */

void *screen_buffer[NWINDOWS];  /* Buffer used for rendering off screen. */
int screen_width[NWINDOWS], screen_height[NWINDOWS];
static time_t last_t[NWINDOWS];

extern void *off_screen_buffer;

/* And from config.h in Mesa.../src/config.h */
#define MAX_WIDTH 2048
#define MAX_HEIGHT 2048

/*
  'stat' each file and only update the ones that are newer.
  Filenames are:
    last_athena.gif
*/

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void draw_image(int ifield, char *fname)
{

  struct stat sf;
  FILE *fp;

  if (stat(fname, &sf) != 0) {
    memset(screen_buffer[ifield], 0, 4*sizeof(char)*screen_width[ifield]*screen_height[ifield]);
  } else {
    if (sf.st_mtime > last_t[ifield]) {
      last_t[ifield] = sf.st_mtime;
      /*    printf("--- last_t = %d \n", last_t[ifield]);*/

      fp = fopen(fname, "r");
      if (fp != NULL) {
	fread(screen_buffer[ifield], 4*sizeof(char), screen_width[ifield]*screen_height[ifield], fp);
	fclose(fp);
      }
    }
  }
  glutSetWindow(my_onwindow[ifield]);

  glViewport(0, 0, screen_width[ifield], screen_height[ifield]);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D((double)0, (double)screen_width[ifield], (double)0, (double)screen_height[ifield]);
  glRasterPos2i(0, 0);
  glDrawPixels(screen_width[ifield], screen_height[ifield], GL_RGBA, GL_UNSIGNED_BYTE, screen_buffer[ifield]);
  glutSwapBuffers();

  /*
  off_screen_buffer = screen_buffer[ifield];
  write_gif(fname);
  */
}
void draw_athena(void)
{
#ifdef DEBUG_COLOR
  float r, g, b;
  r = (float)random()/65535.0/65535.0; /* 2^32 - 1 */
  g = (float)random()/65535.0/65535.0; /* 2^32 - 1 */
  b = (float)random()/65535.0/65535.0; /* 2^32 - 1 */
  fprintf(stderr, "draw_athena: r/g/b = %f/%f/%f\n", r,g,b);
  glClearColor(r, g, b, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glutSwapBuffers();
#endif
  draw_image(0, "last_athena.raw");
}

void visibility_athena(int state)
{
  /*  fprintf(stderr, "visibility_athena: %d\n", state); */
  if (state == GLUT_VISIBLE)
    draw_image(0, "last_athena.raw");
}

void reshape_athena(int w, int h)
{
  /*  fprintf(stderr, "reshape_athena\n"); */
  draw_image(0, "last_athena.raw");
}

void quit_graphics(int signum)
{
  fprintf(stderr, "quit_graphics\n");
  exit(0);
}

void new_image_handler(int signum)
{
  int i;
  for (i=0; i<NWINDOWS; i++) {
    glutSetWindow(my_onwindow[i]);
    glutPostRedisplay();
  }
  signal(SIGUSR1, new_image_handler);
}  

void toggle_view_of_windows(void)
{
  /* For now, jut iconify them. */
  /* FIXME: There is a BUG when the windows become visible again,
   * after iconfication, the display functions are not called!  Why?
   * It seems that the event loop is still intact because the
   * keyboard callbacks work.
   */
  int i;
  for (i=0; i<NWINDOWS; i++) {
    glutSetWindow(my_onwindow[i]);
    glutIconifyWindow();
  }
}

void anykey(unsigned char key, int x, int y)
{
  if (key == 'q' || key == 'Q') {
    system("rm *.raw 2> /dev/null");
    exit(0);
  }
  if (key == 'i' || key == 'i') {
    toggle_view_of_windows();
  }
  
}

void athena_graphics(int xdim, int ydim, 
		   int xpos[], int ypos[],
		   char *title[], int nfield) 
{
  int argc = 0 ;      /* Fake this for glutInit. */
  char **argv = NULL; /* Fake this for glutInit. */
  int i;

  /*  fprintf(stderr, "BEGIN -- athena_graphics\n"); */

  signal(SIGUSR1, new_image_handler);
  signal(SIGHUP, quit_graphics);

  argc = 1;
  argv = (char **) calloc(argc, sizeof(char *));
  argv[0] = "athena_graphics";
  argv[1] = NULL;
  
  for (i=0; i<nfield; i++) {
    screen_width[i] = xdim;
    screen_height[i] = ydim;
    screen_buffer[i] = calloc( xdim * ydim, 4*sizeof(char) );
    if (screen_buffer[i] == NULL) {
      perror("athena_graphics");
      return;
    }
  }

  glutInit(&argc, argv);
    
  for (i=0; i<nfield; i++) {
    glutInitWindowSize(xdim, ydim);
    glutInitWindowPosition(xpos[i], ypos[i]);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    my_onwindow[i] = glutCreateWindow(title[i]);
    glutKeyboardFunc(anykey);
    glutDisplayFunc(draw_athena);
    glutVisibilityFunc(visibility_athena);
    glutReshapeFunc(reshape_athena);
    /*    fprintf(stderr, "glClearColor\n"); */
    glClearColor(1.0, 1.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
    glClearColor(1.0, 1.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
  }

  glutMainLoop();
}
#endif
