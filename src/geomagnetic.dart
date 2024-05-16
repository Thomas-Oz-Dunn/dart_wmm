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

    double dec() => this.d;
    double dip() => this.i;
    double inclination() => this.i;
    double ti() => this.f;
    double total_intensity() => this.f;
  
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
        
    // Calculate the uncertainty values for this ``GeoMagResult``.
    // Uncertainty estimates provided by the **WMM2015** and **WMM2020** error model for the various field components.
    // H is expressed in nT in the formula providing the error in D.
    // These values can currently only be computed for ``GeoMagResult`` between 2015.0 and 2025.0 and using a value
    // outside this will raise an Exception
    GeoMagUncertaintyResult calculate_uncertainty() => GeoMagUncertaintyResult(result: this);

    }

}

    
class GeoMagUncertaintyResult{
    late double x;
    late double y;
    late double z;
    late double h;
    late double f;
    late double i;
    late double d;
    // The uncertainty values of a ``GeoMagResult``.

    // - **f** *(float)* – Uncertainty of the Total Intensity in nT
    // - **h** *(float)* – Uncertainty of the Horizontal Intensity in nT
    // - **x** *(float)* – Uncertainty of the North Component in nT
    // - **y** *(float)* – Uncertainty of the East Component in nT
    // - **z** *(float)* – Uncertainty of the Vertical Component in nT
    // - **i** *(float)* – Uncertainty of the Geomagnetic Inclination in degrees
    // - **d** *(float)* – Uncertainty of the Geomagnetic Declination (Magnetic Variation) in degrees

    GeoMagUncertaintyResult({
      required GeoMagResult result
    }){
      if (2020.0 <= result.time.year && result.time.year <= 2025.0){
          this._error_model_wmm_2020(result);
      } else if (2015.0 <= result.time.year && result.time.year < 2020.0){
          this._error_model_wmm_2015(result);
      } else{
        Exception("GeoMagResult outside of known uncertainty estimates.");
      }
    }

    void _error_model_wmm_2015(result){
        // Calculate uncertainty estimates for 2015.0 to 2020.0
        x = 138.0;
        y = 89.0;
        z = 165.0;
        h = 133.0;
        f = 152.0;
        i = 0.22;
        d = sqrt(pow(0.23,2) + pow(5430 / result.h, 2));
    }

    void _error_model_wmm_2020(result){
        // Calculate uncertainty estimates for 2020.0 to 2025.0.
        x = 131.0;
        y = 94.0;
        z = 157.0;
        h = 128.0;
        f = 148.0;
        i = 0.21;
        d = sqrt(pow(0.26, 2) + pow(5625 / result.h, 2));
    }

}
