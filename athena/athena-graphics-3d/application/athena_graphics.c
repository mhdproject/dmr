#include "athena_preprocessor.fh"
#include <stdio.h>
#include <stdlib.h>
#ifdef INTERACTIVE_GRAPHICS
#include <signal.h>

#define MAX_WIDTH 2048
#define MAX_HEIGHT 2048
void quit_graphics(int signum);

int main(int argc, char **argv)
{
  int i;

  int graphic_pid = 0;

#define NWINDOWS 1
#define XDIM 800
#define YDIM 800
  int xdim = XDIM; /* Should be argument. */
  int ydim = YDIM; /* Should be argument. */

  int xpos[NWINDOWS] = {0};
  int ypos[NWINDOWS] = {0};
  char *title[NWINDOWS] = {"Athena"};
  int nfield = 1;

  FILE *fp;
  char *home;
  char fname[200];

  fprintf(stderr, "BEGIN -- athena_graphics\n");
  fprintf(stderr, "   Use 'q' or 'Q', to quit.\n");
  fprintf(stderr, "   Use 'i' or 'I', to iconify.\n");

  system("rm -f *.raw 2> /dev/null");
  if (xdim > MAX_WIDTH) {
    fprintf(stderr, "Error xdim=%d exceeds MAX_WIDTH=%d\n", xdim, MAX_WIDTH);
    exit(-1);
  }

  if (ydim > MAX_HEIGHT) {
    fprintf(stderr, "Error ydim=%d exceeds MAX_HEIGHT=%d\n", ydim, MAX_HEIGHT);
    exit(-1);
  }

  signal(SIGHUP, quit_graphics);

  graphic_pid = getpid();
  home = getenv("HOME");
  sprintf(fname,"%s/%s",home,".athena_graphics.pid");
  fp = fopen(fname, "w");
  if (fp == NULL) {
    perror(fname);
    exit(-1);
  }
  fprintf(fp, "%d", graphic_pid);
  fclose(fp);

  setpgid(getpid(), getppid());
  athena_graphics(xdim, ydim, xpos, ypos, title, nfield); 

}

#else
int main(int ac, char *av)
{
  printf("Athena is built with '#undef INTERACTIVE_GRAPHICS'.\n");
  printf("Exiting....\n");
  exit(0);
}
#endif
