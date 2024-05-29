/*
Based on pygeomag

*/
import 'dart:io';
import 'dart:math';

const double twoPi = pi * 2;
const double degToRad = pi / 180.0;

/// The Magnetic Components values from ``GeoMag.calculate()``.
///
/// - **glat** *(float)* – Geodetic Latitude, -90.00 to +90.00 degrees (North positive, South negative)
/// - **glon** *(float)* – Geodetic Longitude, -180.00 to +180.00 degrees (East positive, West negative)
/// - **alt** *(float)* – Altitude, -1 to 850km referenced to the WGS 84 ellipsoid OR the Mean Sea Level (MSL)
/// - **time** *(float)* – Time (in decimal year), 2020.0 to 2025.0
/// - **f**, **.ti**, **.total_intensity** *(float)* – Total Intensity
/// - **h** *(float)* – Horizontal Intensity
/// - **north_comp** *(float)* – North Component
/// - **east_comp** *(float)* – East Component
/// - **up_comp** *(float)* – Vertical Component
/// - **i**, **.dip**, **.inclination** *(float)* – Geomagnetic Inclination
/// - **d**, **.dec** *(float)* – Geomagnetic Declination (Magnetic Variation)
/// - **gv** *(float)* – Magnetic grid variation if the current geodetic position is in the arctic or antarctic

class GeoMagResult {
  double dec_year;
  double alt;
  double glat;
  double glon;
  late double north_comp;
  late double east_comp;
  late double up_comp;
  late double h;
  late double f;
  late double i;
  late double d;
  late double gv;
  bool in_blackout_zone = false;
  bool in_caution_zone = false;

  GeoMagResult(this.dec_year, this.alt, this.glat, this.glon);

  double dec() => this.d;
  double dip() => this.i;
  double inclination() => this.i;
  double ti() => this.f;
  double total_intensity() => this.f;

  void calculate(bool raise_in_warning_zone) {
    // Calculate extra result values.
    // COMPUTE X, Y, Z, AND H COMPONENTS OF THE MAGNETIC FIELD
    double cosi = cos(degToRad * i);
    double cosd = cos(degToRad * d);
    north_comp = f * (cosd * cosi);
    east_comp = f * (cosi * sin(degToRad * d));
    up_comp = f * (sin(degToRad * i));
    // h should be e3 not e30
    h = f * cosi;

    // Check if in Caution or Blackout Zones
    if (h < 2000) {
      if (raise_in_warning_zone) {
        // The Blackout Zones are defined as regions around the north and south magnetic poles where the horizontal intensity
        // of Earth's magnetic field (H) is less than 2000 nT. In these zones WMM declination values are inaccurate and
        // compasses are unreliable.
        Exception(
            "The horizontal field strength at this location is $h. Compass readings have VERY LARGE "
            "uncertainties in areas where H smaller than 2000 nT");
      }
      in_blackout_zone = true;
    } else if (h < 6000) {
      if (raise_in_warning_zone) {
        // A Caution Zone is an areas around a Blackout Zone where caution must be exercised while using a compass. It is
        // defined where Earth's magnetic field (H) is between 2000 and 6000 nT. Compass accuracy may be degraded in this
        // region.
        Exception(
            "The horizontal field strength at this location is $h. Compass readings have large "
            "uncertainties in areas where H smaller than 6000 nT");
      }
      in_caution_zone = true;
    }
  }

  // Calculate the uncertainty values for this ``GeoMagResult``.
  // Uncertainty estimates provided by the **WMM2015** and **WMM2020** error model for the various field components.
  // H is expressed in nT in the formula providing the error in D.
  // These values can currently only be computed for ``GeoMagResult`` between 2015.0 and 2025.0 and using a value
  // outside this will raise an Exception
  GeoMagUncertaintyResult calculate_uncertainty() =>
      GeoMagUncertaintyResult(this);
}

/// The uncertainty values of a ``GeoMagResult``.
///
/// - **f** *(float)* – Uncertainty of the Total Intensity in nT
/// - **h** *(float)* – Uncertainty of the Horizontal Intensity in nT
/// - **north_comp_unc** *(float)* – Uncertainty of the North Component in nT
/// - **east_comp_unc** *(float)* – Uncertainty of the East Component in nT
/// - **up_comp_unc** *(float)* – Uncertainty of the Vertical Component in nT
/// - **i** *(float)* – Uncertainty of the Geomagnetic Inclination in degrees
/// - **d** *(float)* – Uncertainty of the Geomagnetic Declination (Magnetic Variation) in degrees
class GeoMagUncertaintyResult {
  late double north_comp_unc;
  late double east_comp_unc;
  late double up_comp_unc;
  late double h;
  late double f;
  late double i;
  late double d;
  GeoMagUncertaintyResult(GeoMagResult result) {
    if (2020.0 <= result.dec_year && result.dec_year <= 2025.0) {
      this._error_model_wmm_2020(result);
    } else if (2015.0 <= result.dec_year && result.dec_year < 2020.0) {
      this._error_model_wmm_2015(result);
    } else {
      Exception("GeoMagResult outside of known uncertainty estimates.");
    }
  }

  void _error_model_wmm_2015(GeoMagResult result) {
    // Calculate uncertainty estimates for 2015.0 to 2020.0
    north_comp_unc = 138.0;
    east_comp_unc = 89.0;
    up_comp_unc = 165.0;
    h = 133.0;
    f = 152.0;
    i = 0.22;
    print('result.h ${result.h}');
    d = sqrt(0.23 * 0.23 + 5430 / result.h * 5430 / result.h);
    print('d $d');
  }

  void _error_model_wmm_2020(GeoMagResult result) {
    // Calculate uncertainty estimates for 2020.0 to 2025.0.
    north_comp_unc = 131.0;
    east_comp_unc = 94.0;
    up_comp_unc = 157.0;
    h = 128.0;
    f = 148.0;
    i = 0.21;
    print('result.h ${result.h}');
    d = sqrt(0.26 * 0.26 + 5625 / result.h * 5625 / result.h);
    print('d $d');
  }
}

/// Dart port of Python port of the Legacy C code provided by NOAA for the World Magnetic Model (WMM).
///
/// It defaults to using the WMM-2020 Coefficient file (WMM.COF) valid for 2020.0 - 2025.0.
/// Included are the following coefficient files:
///
///     | File           | Model      | Life Span       | Creation |
///  -  | WMM_2020.csv   |  WMM-2020  | 2020.0 - 2025.0 | 12/10/2019 |
///  -  | WMM_2015.csv   |  WMM-2015  | 2015.0 - 2020.0 | 12/15/2014 |
class GeoMag {
  var _coefficients_data;
  var _coefficients_file;
  late int _maxord;
  double? _epoch;
  late String _model;
  late DateTime _release_date;
  late var _c;
  late var _cd;
  late var _p;
  late var _fn;
  late var _fm;
  late var _k;

  GeoMag({coefficients_file, coefficients_data}) {
    if (coefficients_file != null && coefficients_data != null) {
      Exception(
          "Both coefficients_file and coefficients_data supplied, supply none or only one.");
    }
    _coefficients_data = coefficients_data;
    _coefficients_file = coefficients_file;
    _maxord = 12;
  }

  model() {
    // Return the life span for the selected coefficient file.
    if (_epoch == null) {
      _load_coefficients();
    }
    return _model;
  }

  List<DateTime> lifespan() {
    // Return the model name for the selected coefficient file.
    if (_epoch == null) {
      _load_coefficients();
    }
    return [_release_date, _release_date.add(Duration(days: 365 * 5))];
  }

  DateTime release_date() {
    // Return the release date for the selected coefficient file.
    if (_epoch == null) {
      _load_coefficients();
    }
    return _release_date;
  }

  /// Determine the model filename to load the coefficients from.
  File _get_model_filename() {
    if (_coefficients_file != null) {
      if (_coefficients_file is File) {
        if (_coefficients_file.existsSync()) {
          return _coefficients_file;
        }
      } else {
        if (File(_coefficients_file).existsSync()) {
          return File(_coefficients_file);
        }
      }
    }

    String sep = (File('.').absolute.path.contains("/")) ? "/" : "\\";
    File filepath = File('./src');
    if (_coefficients_file != null) {
      return File(filepath.path + sep + _coefficients_file);
    }

    _coefficients_file = "wmm${sep}WMM.csv";
    File wmm_filepath = File(filepath.path + sep + _coefficients_file);

    if (wmm_filepath.existsSync()) {
      _coefficients_file = wmm_filepath;
      return wmm_filepath;
    }
    return File('src${sep}wmm${sep}WMM.csv');
  }

  /// Read coefficients data from file to be processed by ``_load_coefficients``.
  ((double, String, DateTime), List<dynamic>)
      _read_coefficients_data_from_file() {
    File model_filename = _get_model_filename();
    List<String> lines = model_filename.readAsLinesSync();

    // Header
    var header = lines[0].replaceAll(' ', '').split(',');
    if (header.length != 3) {
      Exception("Corrupt header in model file");
    }
    double epoch = double.parse(header[0]);
    String model = header[1];
    DateTime release_date = DateTime.parse(header[2]);

    // Coefficients
    var data = [];
    int i = 0;
    for (var line in lines) {
      if (i > 0) {
        var data_line = line.replaceAll(' ', '').split(',');
        if (data_line.length != 6) {
          Exception("Corrupt record in model file");
        }
        var line_list = [
          int.parse(data_line[0]),
          int.parse(data_line[1]),
          double.parse(data_line[2]),
          double.parse(data_line[3]),
          double.parse(data_line[4]),
          double.parse(data_line[5])
        ];
        data += [line_list];
      }
      i += 1;
    }

    return ((epoch, model, release_date), data);
  }

  List<List<double>> _create_matrix(size, {double val = 0}) {
    return List<List<double>>.filled(size, List<double>.filled(size, val));
  }

  List<double> _create_list(size, {double val = 0}) {
    return List<double>.filled(size, val);
  }

  /// Load the coefficients model to calculate the Magnetic Components from.
  void _load_coefficients() {
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

    if (_coefficients_data != null) {
      List vals = _coefficients_data;
      (epoch, model, release_date) = vals[0];
      coefficients = vals[1];
    } else {
      ((double, String, DateTime), List<dynamic>) vals =
          _read_coefficients_data_from_file();
      (epoch, model, release_date) = vals.$1;
      coefficients = vals.$2;
    }

    // READ WORLD MAGNETIC MODEL SPHERICAL HARMONIC COEFFICIENTS
    c[0][0] = 0.0;
    cd[0][0] = 0.0;
    for (var setsd in coefficients) {
      var n = setsd[0];
      var m = setsd[1];
      var gnm = setsd[2];
      var hnm = setsd[3];
      var dgnm = setsd[4];
      var dhnm = setsd[5];
      if (m > _maxord) {
        break;
      } else if (m > n || m < 0) {
        Exception("Corrupt record in model file");
      } else if (m <= n) {
        c[m][n] = gnm;
        cd[m][n] = dgnm;
        if (m != 0) {
          c[n][m - 1] = hnm;
          cd[n][m - 1] = dhnm;
        }
      }
    }

    // CONVERT SCHMIDT NORMALIZED GAUSS COEFFICIENTS TO UNNORMALIZED
    snorm[0] = 1.0;
    fm[0] = 0.0;

    for (int n = 1; n < _maxord + 1; n += 1) {
      snorm[n] = snorm[n - 1] * (2 * n - 1) / n;
      int j = 2;
      int m = 0;
      int D1 = 1;
      double D2 = (n - m + D1) / D1;
      while (D2 > 0.0) {
        k[m][n] = ((n - 1) * (n - 1)) - (m * m) / ((2 * n - 1) * (2 * n - 3));
        if (m > 0) {
          double flnmj = sqrt((n - m + 1) * j / (n + m));
          snorm[n + m * 13] = snorm[n + (m - 1) * 13] * flnmj;
          j = 1;
          c[n][m - 1] = snorm[n + m * 13] * c[n][m - 1];
          cd[n][m - 1] = snorm[n + m * 13] * cd[n][m - 1];
        }
        c[m][n] = snorm[n + m * 13] * c[m][n];
        cd[m][n] = snorm[n + m * 13] * cd[m][n];
        D2 -= 1;
        m += D1;
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

  /// Calculate the Magnetic Components from a latitude, longitude, altitude and date.
  ///
  /// - **float glat** Geodetic Latitude, -90.00 to +90.00 degrees (North positive, South negative)
  /// - **float glon** Geodetic Longitude, -180.00 to +180.00 degrees (East positive, West negative)
  /// - **float alt** Altitude, -1 to 850km referenced to the WGS 84 ellipsoid OR the Mean Sea Level (MSL)
  /// - **float time** Time (in decimal year), 2020.0 to 2025.0
  /// - **bool allow_date_outside_lifespan** True, if you want an estimation outside the 5-year life span
  /// - **bool raise_in_warning_zone** True if you want to raise a BlackoutZoneException or CautionZoneException
  ///   exception when the horizontal intensity is < 6000
  /// return type: GeoMagResult
  ///
  ///  Calculate the geomagnetic declination at the Space Needle in Seattle, WA:
  ///  ```dart
  ///   import 'dart_wmm/geomagnetic.dart'
  ///   geo_mag = GeoMag();
  ///   result = geo_mag.calculate(glat=47.6205, glon=-122.3493, alt=0, time=2023.75);
  ///   print(result.d);
  ///   15.25942260585284
  /// ```
  GeoMagResult calculate(
    double glat,
    double glon,
    double alt,
    double dec_year, {
    bool allow_date_outside_lifespan = false,
    bool raise_in_warning_zone = false,
  }) {
    List<List<double>> tc = _create_matrix(13);
    List<List<double>> dp = _create_matrix(13);
    List<double> sinp = _create_list(169);
    List<double> cosp = _create_list(13);
    List<double> pp = _create_list(13);

    // # INITIALIZE CONSTANTS
    sinp[0] = 0.0;
    cosp[0] = pp[0] = 1.0;
    dp[0][0] = 0.0;

    // Shape of earth
    double aE = 6378.137;
    double bE = 6356.7523142;
    double rE = 6371.2;
    double aE2 = aE * aE;
    double bE2 = bE * bE;
    double cE2 = aE2 - bE2;
    double aE4 = aE2 * aE2;
    double bE4 = bE2 * bE2;
    double cE4 = aE4 - bE4;

    _load_coefficients();
    double dt = dec_year - _epoch!;
    if ((dt < 0.0 || dt > 5.0) && !allow_date_outside_lifespan) {
      Exception("Time extends beyond model 5-year life span");
    }

    double rlon = degToRad * glon;
    double rlat = degToRad * glat;
    double sin_rlon = sin(rlon);
    double sin_rlat = sin(rlat);
    double cos_rlon = cos(rlon);
    double cos_rlat = cos(rlat);
    double sin_rlat2 = sin_rlat * sin_rlat;
    double cos_rlat2 = cos_rlat * cos_rlat;
    sinp[1] = sin_rlon;
    cosp[1] = cos_rlon;

    //  CONVERT FROM GEODETIC COORDINATES TO SPHERICAL COORDINATES
    double q = sqrt(aE2 - cE2 * sin_rlat2);
    double q1 = alt * q;
    double q2 = (q1 + aE2) / (q1 + bE2);
    double ct = sin_rlat / sqrt(q2 * q2 * cos_rlat2 + sin_rlat2);
    double st = sqrt(1.0 - (ct * ct));
    double r2 = (alt * alt) + 2.0 * q1 + (aE4 - cE4 * sin_rlat2) / (q * q);
    double r = sqrt(r2);  // km
    double d = sqrt(aE2 * cos_rlat2 + bE2 * sin_rlat2);
    double ca = (alt + d) / r;
    double sa = cE2 * cos_rlat * sin_rlat / (r * d);
    for (var m = 2; m < _maxord + 1; m += 1) {
      sinp[m] = sinp[1] * cosp[m - 1] + cosp[1] * sinp[m - 1];
      cosp[m] = cosp[1] * cosp[m - 1] - sinp[1] * sinp[m - 1];
    }

    double aor = rE / r;
    double ar = aor * aor;
    double b_rad = 0.0;
    double b_tan = 0.0;
    double b_pol = 0.0;
    double bpp = 0.0;

    for (int n = 1; n < _maxord + 1; n += 1) {
      ar = ar * aor;
      var m = 0;
      var D3 = 1;
      var D4 = (n + m + D3) / D3;
      while (D4 > 0) {
        //  COMPUTE UNNORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
        //  AND DERIVATIVES VIA RECURSION RELATIONS
        var i = n + m * 13;
        var j = n - 1 + (m - 1) * 13;
        if (n == m) {
          _p[i] = st * _p[j];
          dp[m][n] = st * dp[m - 1][n - 1] + ct * _p[j];
        } else if (n == 1 && m == 0) {
          _p[i] = ct * _p[i - 1];
          dp[m][n] = ct * dp[m][n - 1] - st * _p[i - 1];
        } else if (n > 1 && n != m) {
          if (m > n - 2) {
            _p[i - 2] = 0.0;
            dp[m][n - 2] = 0.0;
          }
          _p[i] =
              ct * _p[i - 1] - _k[m][n] * _p[i - 2];
          dp[m][n] = ct * dp[m][n - 1] -
              st * _p[i - 1] -
              _k[m][n] * dp[m][n - 2];
        }

        // # TIME ADJUST THE GAUSS COEFFICIENTS
        tc[m][n] = _c[m][n] + dt * _cd[m][n];
        if (m != 0) {
          tc[n][m - 1] = _c[n][m - 1] + dt * _cd[n][m - 1];
          print('tc[n][m-1] ${tc[n][m-1]}');
        }
        print('tc[m][n] ${tc[m][n]}');

        // ACCUMULATE TERMS OF THE SPHERICAL HARMONIC EXPANSIONS
        double par = ar * _p[n + m * 13];
        double temp1;
        double temp2;
        if (m == 0) {
          temp1 = tc[m][n] * cosp[m];
          temp2 = tc[m][n] * sinp[m];
        } else {
          var tc2 = tc[n][m - 1];
          temp1 = tc[m][n] * cosp[m] + tc2 * sinp[m];
          temp2 = tc[m][n] * sinp[m] - tc2 * cosp[m];
        }
        // investigate temp1
        // b_tan should be e3 not e30
        // print('ar ${ar} temp1 ${temp1} b_tan ${b_tan} dp[m][n] ${dp[m][n]}');
        b_tan = b_tan - ar * temp1 * dp[m][n];
        b_pol += _fm[m] * temp2 * par;
        b_rad += _fn[n] * temp1 * par;

        // SPECIAL CASE:  NORTH/SOUTH GEOGRAPHIC POLES
        if (st == 0.0 && m == 1) {
          if (n == 1) {
            pp[n] = pp[n - 1];
          } else {
            pp[n] = ct * pp[n - 1] - _k[m][n] * pp[n - 2];
          }
          bpp += _fm[m] * temp2 * ar * pp[n];
        }

        D4 -= 1;
        m += D3;
      }
    }

    if (st == 0.0) {
      b_pol = bpp;
    } else {
      b_pol /= st;
    }

    //  ROTATE MAGNETIC VECTOR COMPONENTS FROM SPHERICAL TO
    //  GEODETIC COORDINATES
    // bt should be e3 not e30
    print('bt ${b_tan}');
    print('br ${b_rad}');
    var bx = -b_tan * ca - b_rad * sa;
    print('bx ${bx}');
    var by = b_pol;
    var bz = b_tan * sa - b_rad * ca;

    GeoMagResult result = GeoMagResult(dec_year, alt, glat, glon);

    // COMPUTE DECLINATION (DEC), INCLINATION (DIP) AND
    // TOTAL INTENSITY (TI)
    var bh = sqrt((bx * bx) + (by * by));
    print('bh ${bh}');
    print('bz ${bz}');
    // f should be e3 not e30
    result.f = sqrt((bh * bh) + (bz * bz));
    print('result.f ${result.f}');
    result.d = atan2(by, bx) / degToRad;
    result.i = atan2(bz, bh) / degToRad;

    // COMPUTE MAGNETIC GRID VARIATION IF THE CURRENT
    // GEODETIC POSITION IS IN THE ARCTIC OR ANTARCTIC
    // (I.E. GLAT > +55 DEGREES OR GLAT < -55 DEGREES)
    //  OTHERWISE, SET MAGNETIC GRID VARIATION TO -999.0
    result.gv = -999.0;
    double gv_default = -999.0;

    if (glat.abs() >= 55.0) {
      if (glat > 0.0 && glon >= 0.0) {
        result.gv = result.d - glon;
      } 
      if (glat > 0.0 && glon < 0.0) {
        result.gv = result.d + glon.abs();
      } 
      if (glat < 0.0 && glon >= 0.0) {
        result.gv = result.d + glon;
      } 
      if (glat < 0.0 && glon < 0.0) {
        result.gv = result.d - glon.abs();
      }

      if (result.gv > 180.0) {
        result.gv -= 360.0;
      }
      if (result.gv < -180.0) {
        result.gv += 360.0;
      }
    }

    if (result.gv == gv_default) {
      // result.gv = null;
    }

    result.calculate(raise_in_warning_zone);
    return result;
  }
}
