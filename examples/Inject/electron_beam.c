#include <math.h>


/**
 * calculate the current density of electron injected to semiconductor surface
 * given a surface point(x,y,z) with unit um, norm(nx, ny, nz) and time with unit s
 * return the electron current density in the unit of A/cm^2
 * @param x   the x location of surface point, unit um
 * @param y   the y location of surface point, unit um
 * @param z   the z location of surface point, unit um
 * @param nx   the norm vector
 * @param ny   the norm vector
 * @param nz   the norm vector
 * @param t   the current time, unit second
 * @return the electron current density. A per square cm
 */
double electron_inject_density(double x, double y, double z, double nx, double ny, double nz, double t)
{
  double center_x = 1.5;
  double center_y = 0.0;
  double center_z = 0.0;

  double d = sqrt( (x-center_x)*(x-center_x) + (y-center_y)*(y-center_y) + (z-center_z)*(z-center_z) );
  double d_sigma = 0.2;

  double t_center = 3e-9;
  double t_sigma  = 1e-9;
  return 1.0e3*exp(-d*d/d_sigma/d_sigma)*exp(-(t-t_center)*(t-t_center)/t_sigma/t_sigma);
}


/**
 * calculate the current density of hole injected to semiconductor surface
 * given a surface point(x,y,z) with unit um, norm(nx, ny, nz) and time with unit s
 * return the hole current density in the unit of A/cm^2
 * @param x   the x location of surface point, unit um
 * @param y   the y location of surface point, unit um
 * @param z   the z location of surface point, unit um
 * @param nx   the norm vector
 * @param ny   the norm vector
 * @param nz   the norm vector
 * @param t   the current time, unit second
 * @return the electron current density. A per square cm
 */
double hole_inject_density(double x, double y, double z, double nx, double ny, double nz, double t)
{
  return 0.0;
}
