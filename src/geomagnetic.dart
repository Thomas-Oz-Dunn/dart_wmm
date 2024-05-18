import 'dart:io';
import 'dart:math';
import 'dart:convert';

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

    GeoMagResult(
      this.time,
      this.alt,
      this.glat,
      this.glon
    );

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

class GeoMag{
    //     Dart port of Python port of the Legacy C code provided by NOAA for the World Magnetic Model (WMM).
    //     It defaults to using the WMM-2020 Coefficient file (WMM.COF) valid for 2020.0 - 2025.0.

    //     Included are the following coefficient files, if you have the need to calculate past values:
    //     .. table::
    //        :widths: auto
    //        ==============  ==========  ===============  ==========
    //        File            Model       Life Span        Creation
    //        ==============  ==========  ===============  ==========
    //        WMM.COF         WMM-2020    2020.0 - 2025.0  12/10/2019
    //        WMM_2015v2.COF  WMM-2015v2  2015.0 - 2020.0  09/18/2018
    //        WMM_2015.COF    WMM-2015    2015.0 - 2020.0  12/15/2014
    //        WMM_2010.COF    WMM-2010    2010.0 - 2015.0  11/20/2009
    //        ==============  ==========  ===============  ==========

  var _coefficients_data;
  var _coefficients_file;
  late int _maxord;
  late double _epoch;
  late String _model;
  late DateTime _release_date;
  late var _c;
  late var _cd;
  late var _p;
  late var _fn;
  late var _fm;
  late var _k;

  GeoMag({coefficients_file, coefficients_data}){

    if (coefficients_file != null && coefficients_data!= null){
        Exception("Both coefficients_file and coefficients_data supplied, supply none or only one.");
    }
    _coefficients_data = coefficients_data;
    _coefficients_file = coefficients_file;
    _maxord = 12;

    model(){
      // Return the life span for the selected coefficient file.
      if (_epoch == null){

      }
      return _model;
    }

    List<DateTime> lifespan(){
      // Return the model name for the selected coefficient file.
      if (_epoch == null){

      }
      return [_release_date, _release_date.add(Duration(days: 365 * 5))];
    }

    DateTime release_date(){
      // Return the release date for the selected coefficient file.
      if (_epoch == null){

      }
      return _release_date;
    }

    String _get_model_filename(){
      // Determine the model filename to load the coefficients from.
      if (_coefficients_file != null){
        if (_coefficients_file[0] == "/" || _coefficients_file[0] == "\\"){
          return _coefficients_file;
        }
      }

      File filepath = File('.');
      String sep = (filepath.path.contains("/"))? "/" : "\\";

      if (_coefficients_file != null){
          return filepath.parent.path + sep + _coefficients_file;
      }

      coefficients_file = "wmm${sep}WMM.COF";
      var wmm_filepath = filepath.parent.path + sep + coefficients_file;

      if (File(wmm_filepath).existsSync()){
          _coefficients_file = wmm_filepath;
          return wmm_filepath;
      }
      coefficients_file = "WMM.COF";
      var wmm_filepath2 = filepath.parent.path + sep + coefficients_file;

      if (File(wmm_filepath2).existsSync()){
          _coefficients_file = wmm_filepath2;
          return wmm_filepath2;
      } else{
        return wmm_filepath;
      }
    }

    ((double, String, DateTime), List<dynamic>) _read_coefficients_data_from_file(){
        // Read coefficients data from file to be processed by ``_load_coefficients``.
        var data = [];
        String model_filename = _get_model_filename(); 
        
        List<String> lines = new File(model_filename).readAsLinesSync();

        var header = LineSplitter().convert(lines[0]);
        if (header.length != 3){
            Exception("Corrupt header in model file");
        }
        double epoch = double.parse(header[0]);
        String model = header[1];
        DateTime release_date = DateTime.parse(header[2]);


        int i =0;
        for (var line in lines){
          if (i>0){
          
            var data_line = LineSplitter().convert(line);
            if (data_line.length != 6){
                Exception("Corrupt record in model file");
            }
            data += [
              int.parse(data_line[0]), 
              int.parse(data_line[1]), 
              double.parse(data_line[2]), 
              double.parse(data_line[3]), 
              double.parse(data_line[4]), 
              double.parse(data_line[5])
            ];
          }
          i += 1;
        }

        return ((epoch, model, release_date), data);

    }

    List<List<double>> _create_matrix(size, {double val=0}){
      return List<List<double>>.filled(size, List<double>.filled(size, val));
    }
    List<double> _create_list(size, {double val=0}){
      return List<double>.filled(size, val);
    }

    void _load_coefficients(){
        // Load the coefficients model to calculate the Magnetic Components from.

        if (_epoch != null){
            return;
        }

        List<List<double>> c = _create_matrix(13);
        List<List<double>> cd = _create_matrix(13);
        List<double> snorm = _create_list(169);
        List<double> fn = _create_list(13);
        List<double> fm = _create_list(13);
        List<List<double>> k = _create_matrix(13);

        Iterable coefficients;
        double epoch;
        String model;
        DateTime release_date;

        if (_coefficients_data != null){
            List vals = _coefficients_data;
            (epoch, model, release_date) = vals[0];
            coefficients = vals[1];
        }
        else{
            ((double, String, DateTime), List<dynamic>) vals = _read_coefficients_data_from_file();
            (epoch, model, release_date) = vals.$1;
            coefficients = vals.$2;
        }

        // READ WORLD MAGNETIC MODEL SPHERICAL HARMONIC COEFFICIENTS
        c[0][0] = 0.0;
        cd[0][0] = 0.0;
        for(var (n, m, gnm, hnm, dgnm, dhnm) in coefficients){
            if (m > _maxord){
                break;
            }
            else if (m > n || m < 0){
                Exception("Corrupt record in model file");
            }
            else if (m <= n){
                c[m][n] = gnm;
                cd[m][n] = dgnm;
                if (m != 0){
                    c[n][m - 1] = hnm;
                    cd[n][m - 1] = dhnm;
                  }
            }

        }

        // CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED
        snorm[0] = 1.0;
        fm[0] = 0.0;

        for (int n=1; n<_maxord+1; n+=1){
            snorm[n] = snorm[n - 1] * (2 * n - 1) / n;
            int j = 2;
            int m = 0;
            int D1 = 1;
            double D2 = (n - m + D1) / D1;
            while(D2>0){
                k[m][n] = ((n - 1) * (n - 1)) - (m * m) / ((2 * n - 1) * (2 * n - 3));
                if (m > 0){
                  double flnmj = sqrt((n - m + 1) * j / (n + m));
                  snorm[n + m * 13] = snorm[n + (m - 1) * 13] * flnmj;
                  j = 1;
                  c[n][m - 1] = snorm[n + m * 13] * c[n][m - 1];
                  cd[n][m - 1] = snorm[n + m * 13] * cd[n][m - 1];
                  c[m][n] = snorm[n + m * 13] * c[m][n];
                  cd[m][n] = snorm[n + m * 13] * cd[m][n];
                  D2 -= 1;
                  m += D1;
                }
            }
            fn[n] = n + 1;
            fm[n] = n * 1.0;
        }
        k[1][1] = 0.0;

        _epoch = epoch;
        _model = model;
        _release_date = release_date;
        _c = c;
        _cd = cd;
        _p = snorm;
        _fn = fn;
        _fm = fm;
        _k = k;
    }


    GeoMagResult calculate(
      double glat,
      double glon,
      double alt,
      DateTime time,
    {
        bool allow_date_outside_lifespan=false,
        bool raise_in_warning_zone=false,
    }){

//         Calculate the Magnetic Components from a latitude, longitude, altitude and date.

//         :param float glat: Geodetic Latitude, -90.00 to +90.00 degrees (North positive, South negative)
//         :param float glon: Geodetic Longitude, -180.00 to +180.00 degrees (East positive, West negative)
//         :param float alt: Altitude, -1 to 850km referenced to the WGS 84 ellipsoid OR the Mean Sea Level (MSL)
//         :param float time: Time (in decimal year), 2020.0 to 2025.0
//         :param bool allow_date_outside_lifespan: True, if you want an estimation outside the 5-year life span
//         :param bool raise_in_warning_zone: True if you want to raise a BlackoutZoneException or CautionZoneException
//             exception when the horizontal intensity is < 6000
//         :return type: GeoMagResult

//         Calculate the geomagnetic declination at the Space Needle in Seattle, WA:
//         >>> geo_mag = GeoMag()
//         >>> result = geo_mag.calculate(glat=47.6205, glon=-122.3493, alt=0, time=2023.75)
//         >>> print(result.d)
//         15.25942260585284

//         And calculate it for the same spot 10 years ago:
//         >>> geo_mag = GeoMag(coefficients_file='wmm/WMM_2010.COF')
//         >>> result = geo_mag.calculate(glat=47.6205, glon=-122.3493, alt=0, time=2013.75)
//         >>> print(result.d)
//         16.32554283003356

        List<List<double>> tc = _create_matrix(13);
        List<List<double>> dp = _create_matrix(13);
        List<double> sp = _create_list(169);
        List<double> cp = _create_list(13);
        List<double> pp = _create_list(13);

        // # INITIALIZE CONSTANTS
        sp[0] = 0.0;
        cp[0] = pp[0] = 1.0;
        dp[0][0] = 0.0;

        // Shape of earth
        double a = 6378.137;
        double b = 6356.7523142;
        double re = 6371.2;

        double a2 = a * a;
        double b2 = b * b;
        double c2 = a2 - b2;
        double a4 = a2 * a2;
        double b4 = b2 * b2;
        double c4 = a4 - b4;

        _load_coefficients();

        //  TODO #1: Legacy C code static vars for speed
        //   Decide to either:
        //    1. Pull out the tracking of previous values
        //         which are in the legacy c app for speeding it up for getting the Secular Change
        //    2. Remove them
        //  otime = oalt = olat = olon = -1000.0

        Duration dt = time.difference(DateTime(_epoch.toInt()));
        if (dt.inDays < 0.0 || dt.inDays > (5.0*365.25) && !allow_date_outside_lifespan){
            Exception("Time extends beyond model 5-year life span");
        }

        double rlon = degToRad * glon;
        double rlat = degToRad * glat;
        double srlon = sin(rlon);
        double srlat = sin(rlat);
        double crlon = cos(rlon);
        double crlat = cos(rlat);
        double srlat2 = srlat * srlat;
        double crlat2 = crlat * crlat;
        sp[1] = srlon;
        cp[1] = crlon;

        //  CONVERT FROM GEODETIC COORDINATES TO SPHERICAL COORDINATES
        double q = sqrt(a2 - c2 * srlat2);
        double q1 = alt * q;
        double q2 = ((q1 + a2) / (q1 + b2)) * ((q1 + a2) / (q1 + b2));
        double ct = srlat / sqrt(q2 * crlat2 + srlat2);
        double st = sqrt(1.0 - (ct * ct));
        double r2 = (alt * alt) + 2.0 * q1 + (a4 - c4 * srlat2) / (q * q);
        double r = sqrt(r2);
        double d = sqrt(a2 * crlat2 + b2 * srlat2);
        double ca = (alt + d) / r;
        double sa = c2 * crlat * srlat / (r * d);
        for (var m = 2; m<_maxord + 1; m+=1){
            sp[m] = sp[1] * cp[m - 1] + cp[1] * sp[m - 1];
            cp[m] = cp[1] * cp[m - 1] - sp[1] * sp[m - 1];
          
        }

        double aor = re / r;
        double ar = aor * aor;
        double br = 0.0; 
        double bt = 0.0;
        double bp = 0.0;
        double bpp = 0.0;

        for (int n=1; n<_maxord+1; n+=1){
          ar = ar * aor;
          var m = 0;
          var D3 = 1;
          var D4 = (n + m + D3) / D3;
          while (D4 > 0){
            //  COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
            //  AND DERIVATIVES VIA RECURSION RELATIONS
            if (n == m){
                _p[n + m * 13] = st * _p[n - 1 + (m - 1) * 13];
                dp[m][n] = st * dp[m - 1][n - 1] + ct * _p[n - 1 + (m - 1) * 13];
            } else if (n == 1 && m == 0){
                _p[n + m * 13] = ct * _p[n - 1 + m * 13];
                dp[m][n] = ct * dp[m][n - 1] - st * _p[n - 1 + m * 13];
            }
            else if (n > 1 && n != m){
                if (m > n - 2){
                    _p[n - 2 + m * 13] = 0.0;
                }
                if (m > n - 2){
                    dp[m][n - 2] = 0.0;
                }
                _p[n + m * 13] = ct * _p[n - 1 + m * 13] - _k[m][n] * _p[n - 2 + m * 13];
                dp[m][n] = ct * dp[m][n - 1] - st * _p[n - 1 + m * 13] - _k[m][n] * dp[m][n - 2];
            }

            // # TIME ADJUST THE GAUSS COEFFICIENTS
            tc[m][n] = _c[m][n] + dt * _cd[m][n];
            if (m != 0){
                tc[n][m - 1] = _c[n][m - 1] + dt * _cd[n][m - 1];
            }

            // ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
            var par = ar * _p[n + m * 13];
            var temp1;
            var temp2;
            if (m == 0){
                temp1 = tc[m][n] * cp[m];
                temp2 = tc[m][n] * sp[m];
            }
            else{

                temp1 = tc[m][n] * cp[m] + tc[n][m - 1] * sp[m];
                temp2 = tc[m][n] * sp[m] - tc[n][m - 1] * cp[m];
            }
            bt = bt - ar * temp1 * dp[m][n];
            bp += _fm[m] * temp2 * par;
            br += _fn[n] * temp1 * par;

            // SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
            if (st == 0.0 && m == 1){
                if (n == 1){
                    pp[n] = pp[n - 1];
                } else{
                    pp[n] = ct * pp[n - 1] - _k[m][n] * pp[n - 2];
                }
                bpp += _fm[m] * temp2 * ar * pp[n];
            }

            D4 -= 1;
            m += D3;

          }
        }

        if (st == 0.0){
            bp = bpp;
        } else{
            bp /= st;
        }

        //  ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
        //  GEODETIC COORDINATES
        var bx = -bt * ca - br * sa;
        var by = bp;
        var bz = bt * sa - br * ca;

        GeoMagResult result = GeoMagResult(time, alt, glat, glon);


        // COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
        // TOTAL INTENSITY (TI)
        var bh = sqrt((bx * bx) + (by * by));
        result.f = sqrt((bh * bh) + (bz * bz));
        result.d = atan2(by, bx) / degToRad;
        result.i = atan2(bz, bh) / degToRad;


        // COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
        // GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
        // (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)
        //  OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
        result.gv = -999.0;
        double gv_default = -999.0;

        if (glat.abs() >= 55.0){

            if (glat > 0.0 && glon >= 0.0){
                result.gv = result.d - glon;
            }
            if (glat > 0.0 && glon < 0.0){
                result.gv = result.d + glon.abs();
            }
            if (glat < 0.0 && glon >= 0.0){
                result.gv = result.d + glon;
            } 
            if (glat < 0.0 && glon < 0.0){
                result.gv = result.d - glon.abs();
            }
            if (result.gv > 180.0){
                result.gv -= 360.0;
            }
            if (result.gv < -180.0){
                result.gv += 360.0;
            }
        }

        if (result.gv == gv_default){
            // result.gv = null;
        }

        result.calculate(raise_in_warning_zone);
      return result;
    }

  }
}

