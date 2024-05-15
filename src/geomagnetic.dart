import 'dart:math';

const double twoPi = pi * 2;
const double degToRad = pi / 180.0;

class GeoMagResult{
    DateTime time;
    double alt;
    double glat;
    double glon;
    late double x;
    late double y;
    late double z;
    late double h;
    late double f;
    late double i;
    late double d;
    late double gv;
    bool in_blackout_zone = false;
    bool in_caution_zone = false;

    // The Magnetic Components values from ``GeoMag.calculate()``.

    // - **glat** *(float)* – Geodetic Latitude, -90.00 to +90.00 degrees (North positive, South negative)
    // - **glon** *(float)* – Geodetic Longitude, -180.00 to +180.00 degrees (East positive, West negative)
    // - **alt** *(float)* – Altitude, -1 to 850km referenced to the WGS 84 ellipsoid OR the Mean Sea Level (MSL)
    // - **time** *(float)* – Time (in decimal year), 2020.0 to 2025.0

    // - **f**, **.ti**, **.total_intensity** *(float)* – Total Intensity
    // - **h** *(float)* – Horizontal Intensity
    // - **x** *(float)* – North Component
    // - **y** *(float)* – East Component
    // - **z** *(float)* – Vertical Component
    // - **i**, **.dip**, **.inclination** *(float)* – Geomagnetic Inclination
    // - **d**, **.dec** *(float)* – Geomagnetic Declination (Magnetic Variation)
    // - **gv** *(float)* – Magnetic grid variation if the current geodetic position is in the arctic or antarctic

    GeoMagResult({
      required this.time, 
      required this.alt, 
      required this.glat, 
      required this.glon
    });

    double dec(){return this.d;}
    double dip(){return this.i;}
    double inclination(){return this.i;}
    double ti(){return this.f;}
    double total_intensity(){return this.f;}
  
    void calculate(bool raise_in_warning_zone){
        // Calculate extra result values.
        // COMPUTE X, Y, Z, AND H COMPONENTS OF THE MAGNETIC FIELD
        x = f * (cos(degToRad * d) * cos(degToRad * i));
        y = f * (cos(degToRad * i) * sin(degToRad * d));
        z = f * (sin(degToRad * i));
        h = f * (cos(degToRad * i));

        // Check if in Caution or Blackout Zones
        if (h < 2000){
            if (raise_in_warning_zone){
                // The Blackout Zones are defined as regions around the north and south magnetic poles where the horizontal intensity
                // of Earth's magnetic field (H) is less than 2000 nT. In these zones WMM declination values are inaccurate and
                // compasses are unreliable.
                Exception(
                    "The horizontal field strength at this location is $h. Compass readings have VERY LARGE "
                    "uncertainties in areas where H smaller than 2000 nT"
                );
            }
            in_blackout_zone = true;
        } else if (h < 6000){
            if (raise_in_warning_zone){
                // A Caution Zone is an areas around a Blackout Zone where caution must be exercised while using a compass. It is
                // defined where Earth's magnetic field (H) is between 2000 and 6000 nT. Compass accuracy may be degraded in this
                // region.
                Exception(
                    "The horizontal field strength at this location is $h. Compass readings have large "
                    "uncertainties in areas where H smaller than 6000 nT"
                );
            }
            in_caution_zone = true;
        }
    }
}
