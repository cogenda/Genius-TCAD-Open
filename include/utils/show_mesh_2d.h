/********************************************************************************/
/*     888888    888888888   88     888  88888   888      888    88888888       */
/*   8       8   8           8 8     8     8      8        8    8               */
/*  8            8           8  8    8     8      8        8    8               */
/*  8            888888888   8   8   8     8      8        8     8888888        */
/*  8      8888  8           8    8  8     8      8        8            8       */
/*   8       8   8           8     8 8     8      8        8            8       */
/*     888888    888888888  888     88   88888     88888888     88888888        */
/*                                                                              */
/*       A Three-Dimensional General Purpose Semiconductor Simulator.           */
/*                                                                              */
/*                                                                              */
/*  Copyright (C) 2007-2008                                                     */
/*  Cogenda Pte Ltd                                                             */
/*                                                                              */
/*  Please contact Cogenda Pte Ltd for license information                      */
/*                                                                              */
/*  Author: Gong Ding   gdiso@ustc.edu                                          */
/*                                                                              */
/********************************************************************************/

#ifndef   __show_mesh_2d_h__
#define   __show_mesh_2d_h__

#include "config.h"

#ifdef HAVE_X11

#include <vector>
#include <map>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>

class SimulationSystem;

class ShowMesh2D
{
public:

  ShowMesh2D(const SimulationSystem & system, bool inv_y=true);

  int  show_mesh_init();

  void show_mesh(const char * TIFFFileName=NULL);

private:

  const SimulationSystem & _system;

  static char *prog_name;

  /**
   * X11 display
   */
  Display *display;

  /**
   * main window
   */
  Window  main_win;

  /**
   * draw window
   */
  Window  draw_win;

  /**
   * button windows
   */
  Window button[6];

  struct Button_Data
  {
    int          x0, y0;
    int          hght, wdth;
    int          border;
    const char   *caption;
    int          pressed;
  };

  static Button_Data butt_data[6];

  XTextProperty winname, iconame;

  Pixmap ico_pixm, back_pixm;

  int draw_wdth, draw_hght;

  int main_wdth, main_hght;

  int scr_num;

  GC gc_BoW, gc_WoB, gc_XOR, gc_THICK, gc_DASHED, gc_numb;

  /**
   * supported color, 16 VGA color + 20 rainbow color
   */
  GC gc_color[16+20];

  /**
   * rainbow colors
   */
  static const unsigned int rc[20][3];

  XFontStruct *text_font, *numb_font;

  Colormap colormap;

  /**
   * scaling and offset in screen
   */
  double  scl, X0, Y0;

  /**
   * mesh bound box
   */
  double  xmin, xmax, ymin, ymax;

  /**
   * map to screen
   */
  int Xmap(double x)
  { return int( ((x)-(xmin+xmax)/2)*scl+X0); }

  /**
   * map to screen
   */
  int Ymap(double y)
  { return int(-((y)-(ymin+ymax)/2)*scl+Y0); }


  void GetGCbyColorName(Window win, GC *gc, char *pszColorName );

  void GetGCbyRainbow(Window win, GC *gc, unsigned short R, unsigned short G,unsigned short B);

  void getGC(Window win);

  void load_fonts();

  void draw(Window win, GC gc, int win_x_dim, int win_y_dim);

  void place_text(Window win, GC gc, XFontStruct *text_font, int win_x_dim, int win_y_dim, char *string);

  void create_buttons(Window parent);

  void write_on_button(int b);

  void show_mesh_redwraw();

  int  show_mesh_save_screen (const char * tiff_file, Window window, int width, int height);

  void show_mesh_close();

  struct ScreenPoint
  {
    double          x, y;
    int             X, Y;
    double          doping;
    unsigned int    subdomain;
  };

  std::multimap<int, ScreenPoint> screen_points;

  ScreenPoint & screen_point(unsigned int region, int node)
  {
    typedef std::multimap<int, ScreenPoint>::iterator It;
    std::pair<It, It> pos = screen_points.equal_range(node);
    It it = pos.first;
    while(it!=pos.second)
    {
      if(it->second.subdomain == region)
        return it->second;
      ++it;
    }
    //prevent warning
    return screen_points.begin()->second;
  }

  struct ScreenFace
  {
    int A;
    int B;
    int C;
    int D;
    unsigned int subdomain;
  };

  std::vector<ScreenFace> screen_faces;

  void draw_junction(const ScreenPoint &A, const ScreenPoint &B, const ScreenPoint &C);

};

#endif // HAVE_X11

#endif // __show_mesh_2d_h__
