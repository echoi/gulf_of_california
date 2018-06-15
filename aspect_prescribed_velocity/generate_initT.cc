#include <iostream>
#include <math.h>
#include <fstream>

int main()
{

    const double diffusivity = 1.0e-6;
    const double dip_FAR = 15.0 * M_PI / 180.0;
    const double tan_delta = tan(dip_FAR);
    const double cos_delta = cos(dip_FAR);
    const double tan_delta_sqrd = tan_delta * tan_delta;
    const double strike_Ridge = 7.0 * M_PI / 180.0;
    const double tan_theta = tan(strike_Ridge);
    const double cos_del_cos_theta = cos( dip_FAR * cos(strike_Ridge) );
    const double obliquity_FAR = 20.0 * M_PI / 180.0;
    const double cos_alpha_FAR = cos(obliquity_FAR);
    const double obliquity_PAC = 330.0 * M_PI / 180.0;
    const double cos_alpha_PAC = cos(obliquity_PAC);
    const double YEAR2SEC = 3.1536e7;
    const double v_half = 0.03 / YEAR2SEC; // m/yr to m/s.

    const double xt = 100.0e3;
    const double x0 = 0.0e3;
    const double surface_temperature = 293.0;
    const double mantle_temperature = 1593.0;
            
    const int nx = 701;
    const int nz = 351;
    const double xmin = 0.0e3;
    const double xmax = 700.0e3;
    const double zmin = 0.0e3;
    const double zmax = 350.0e3;
  
    std::ofstream initT_file;
    initT_file.open("goc_initT.txt");
    initT_file << "# POINTS: " << nx <<" "<< nz << std::endl;
    initT_file << "# Columns: x y temperature [K]" << std::endl;
    for( int j = 0; j < nz; ++j ) {
        for (int i = 0; i < nx; ++i) {
            int counter = i + j * nx;
            double xn = i*(xmax-xmin)/(nx-1);
            double zn = j*(zmax-zmin)/(nz-1);
            double yn = 0.0;
            double dn  = zmax - zn;
            double age = 0.0;
            double depth = 0.0;
            // If in or beneath unsubducted plates
            if( xn <= xt ) {
#if 0
                double xr = (x0 + tan_theta * (tan_theta*xn-yn))/(1.0+tan_theta*tan_theta);
                double yr = tan_theta * (xr-xn) + yn;
                double dist = sqrt( (xr-xn)*(xr-xn) + (yr-yn)*(yr-yn) );
                age = dist / (v_half * cos_alpha_FAR);
#endif
                age = xn / v_half ;
                depth = dn;
                // std::cerr << "Unsubducted plates: x="<< xn <<" age="<<age<<std::endl;
            }
            // If in the NA plate
            else if( xn > xt && dn < ( xn - xt ) * tan_delta ) {
                age = 100.0e6 * YEAR2SEC; // 100 My in seconds
                depth = dn;
                // std::cerr << "In NA plate: x="<< xn <<" age="<<age<<std::endl;
            }
            // If in or beneath the subducted plate
            else if( xn > xt && dn >= ( xn - xt ) * tan_delta ) {
#if 0
                double xs = (xn - tan_delta * zn + tan_delta_sqrd * xt) / (1.0 + tan_delta_sqrd );
                double ys = yn;
                double zs = -tan_delta * ( xs - xt );
                double xr = (x0 + tan_theta * (tan_theta*xs-ys))/(1.0+tan_theta*tan_theta);
                double yr = tan_theta * (xr-xs) + ys;
                double yt = tan_theta * (xt-xs) + ys;
                double d_tr = sqrt( (xr-xt)*(xr-xt) + (yr-yt)*(yr-yt) );
                double d_st = sqrt( (xt-xs)*(xt-xs) + (yt-ys)*(yt-ys) );
                age = d_tr / (v_half * cos_alpha_FAR) + d_st / ( v_half * cos_alpha_FAR * cos_del_cos_theta);
                depth = zmax - sqrt( (xs-xn)*(xs-xn) + (zs-zn)*(zs-zn) );
#endif
                age = xn / ( v_half * cos_delta );
                depth = ( dn - ( xn - xt ) * tan_delta ) * cos_delta;
                std::cerr << "In subducted Farallon plate: x="<< xn <<" age="<< age 
                          << " depth="<< depth <<" dn="<<dn<<" xn="<<xn<<" xt="<<xt 
                          << " tan_delta="<<tan_delta<<" cos_delta="<<cos_delta<<std::endl;
            }
           
            if(age==0.0) age = 0.1e6 * YEAR2SEC; 
            double w = 0.5 * depth / sqrt( diffusivity * age );
            double temperature = surface_temperature + ( mantle_temperature - surface_temperature ) * erf(w);
            initT_file << xn <<" "<< zn <<" "<< temperature << std::endl;
            // std::cerr << counter <<" "<< xn <<" "<< zn <<" "<< temperature << std::endl;
        }
    }
    initT_file.close();

    return 0;
}
