import 'dart:io';
import 'package:test/test.dart';

import '../src/geomagnetic.dart';

bool test_calculate_uncertainty() {
  GeoMag geo_mag = GeoMag();
  GeoMagResult result = geo_mag.calculate(80, 0, 0, 2020);
  GeoMagUncertaintyResult uncertainty = result.calculate_uncertainty();
  print(uncertainty.d);
  return closeTo(uncertainty.d, 0.01).matches(0.89, {});
}

bool test_d_property() {
  var result = GeoMagResult(2020, 0, 0, 0);
  result.d = 5;
  return result.d == result.dec();
}

bool test_f_properties() {
  var result = GeoMagResult(2020, 0, 0, 0);
  result.f = 5;
  return (result.f == result.ti() && result.f == result.total_intensity());
}

bool test_i_property() {
  var result = GeoMagResult(2020, 0, 0, 0);
  result.i = 5;
  return (result.i == result.dip() && result.i == result.inclination());
}

bool test_static_values_2015() {
  var geo_mag = GeoMag(coefficients_file: "src/wmm/WMM_2015.csv");
  var result = geo_mag.calculate(80, 0, 0, 2015);
  var uncertainty = GeoMagUncertaintyResult(result);
  return (uncertainty.north_comp_unc == 138.0 &&
      uncertainty.east_comp_unc == 89.0 &&
      uncertainty.up_comp_unc == 165.0 &&
      uncertainty.h == 133.0 &&
      uncertainty.f == 152.0 &&
      uncertainty.i == 0.22);
}

bool test_static_values_2020() {
  var geo_mag = GeoMag();
  var result = geo_mag.calculate(80, 0, 0, 2020);
  var uncertainty = GeoMagUncertaintyResult(result);
  return (uncertainty.north_comp_unc == 131.0 &&
      uncertainty.east_comp_unc == 94.0 &&
      uncertainty.up_comp_unc == 157.0 &&
      uncertainty.h == 128.0 &&
      uncertainty.f == 148.0 &&
      uncertainty.i == 0.21);
}

bool test_uncertainty_degrees_2015() {
  var geo_mag = GeoMag(coefficients_file: "src/wmm/WMM_2015.csv");
  var result = geo_mag.calculate(80, 0, 0, 2015.0);
  var uncertainty = GeoMagUncertaintyResult(result);
  bool d1 = closeTo(uncertainty.d, 0.01).matches(0.85, {});

  var result2 = geo_mag.calculate(0, 120, 0, 2015.0);
  var uncertainty2 = GeoMagUncertaintyResult(result2);
  return d1 && closeTo(uncertainty2.d, 0.01).matches(0.27, {});
}

bool test_uncertainty_degrees_2022() {
  var geo_mag = GeoMag();
  var result = geo_mag.calculate(80, 0, 0, 2020.0);
  var uncertainty = GeoMagUncertaintyResult(result);
  bool d1 = closeTo(uncertainty.d, 0.01).matches(0.89, {});

  var result2 = geo_mag.calculate(0, 120, 0, 2020.0);
  var uncertainty2 = GeoMagUncertaintyResult(result2);
  return d1 && closeTo(uncertainty2.d, 0.01).matches(0.3, {});
}

List<double> get_test_values(String test_parameters) {
  List<String> splits = test_parameters.replaceAll(' ', '').split(',');
  return splits.map((e) => double.parse(e)).toList();
}

bool run_tests() {
  GeoMag geo_mag = GeoMag();
  List<bool> t = [];
  for (String test_filename in [
    'test\\WMM2015testvalues.txt',
    'test\\WMM2020testvalues.txt'
  ]) {
    File f = File(test_filename);
    List<String> lines = f.readAsLinesSync();
    int i = 0;
    List<bool> res = [];
    for (var line in lines) {
      if (line[0] != "#") {
        var vals = get_test_values(line);
        var time = vals[0];
        var alt = vals[1];
        var glat = vals[2];
        var glon = vals[3];
        var x = vals[4];
        var y = vals[5];
        var z = vals[6];
        var h = vals[7];
        var f = vals[8];
        var i = vals[9];
        var d = vals[10];
        var gv = vals[11];
        var result = geo_mag.calculate(glat, glon, alt, time);

        res += [
          closeTo(x, 0.1).matches(result.north_comp, {
                false: 'Row $i: X (nT) expected $x, result ${result.north_comp}'
              }) &&
              closeTo(y, 0.1).matches(result.east_comp, {
                false: 'Row $i: Y (nT) expected $y, result ${result.east_comp}'
              }) &&
              closeTo(z, 0.1).matches(result.up_comp,
                  {false: 'Row $i: Z (nT) expected $z, result ${result.up_comp}'}) &&
              closeTo(h, 0.1).matches(result.h,
                  {false: 'Row $i: H (nT) expected $h, result ${result.h}'}) &&
              closeTo(f, 0.1).matches(result.f,
                  {false: 'Row $i: F (nT) expected $f, result ${result.f}'}) &&
              closeTo(i, 0.01).matches(
                  result.i, {false: 'Row $i: I (Deg) expected $i, result ${result.i}'}) &&
              closeTo(d, 0.01).matches(result.d, {false: 'Row $i: D (Deg) expected $d, result ${result.d}'})
        ];
      }
      i += 1;
    }
    t += [!res.any((element) => (element == false))];
  }
  return !t.any((element) => element == false);
}

void main() {
  test('test_calculate_uncertainty', (() => expect(test_calculate_uncertainty(), true)));
  test('test_d_property', (() => expect(test_d_property(), true)));
  test('test_f_properties', (() => expect(test_f_properties(), true)));
  test('test_i_property', (() => expect(test_i_property(), true)));
  test('test_static_values_2015', (() => expect(test_static_values_2015(), true)));
  test('test_static_values_2020', (() => expect(test_static_values_2020(), true)));
  test('test_uncertainty_degrees_2015', (() => expect(test_uncertainty_degrees_2015(), true)));
  test('test_uncertainty_degrees_2022', (() => expect(test_uncertainty_degrees_2022(), true)));
  test('run_tests', (() => expect(run_tests(), true)));
}

// class TestGeoMag(TestCase):
//     def test_both_parameters_supplied_raises(self):
//         with self.assertRaisesRegex(
//             ValueError, "Both coefficients_file and coefficients_data supplied, supply none or only one."
//         ):
//             GeoMag(coefficients_file="wmm/WMM_2020.COF", coefficients_data=WMM_2020)

//     def test_calculate_declination_time_beyond_model_bypass(self):
//         geo_mag = GeoMag()
//         result = geo_mag.calculate(0, 80, 0, 2030, allow_date_outside_lifespan=True)
//         self.assertIsInstance(result, GeoMagResult)

//     def test_calculate_declination_time_beyond_model_raises(self):
//         geo_mag = GeoMag()
//         with self.assertRaisesRegex(ValueError, "Time extends beyond model 5-year life span"):
//             geo_mag.calculate(0, 80, 0, 2030)

//     def test_create_list(self):
//         self.assertEqual(GeoMag._create_list(2), [None, None])
//         self.assertEqual(GeoMag._create_list(3, 0), [0, 0, 0])
//         self.assertNotEqual(GeoMag._create_list(4, 0), [1, 1, 1, 1])

//     def test_create_matrix(self):
//         self.assertEqual(GeoMag._create_matrix(2, 2), [[None, None], [None, None]])
//         self.assertEqual(GeoMag._create_matrix(2, 3, 0), [[0, 0, 0], [0, 0, 0]])
//         self.assertNotEqual(GeoMag._create_matrix(2, 4), [[1, 1, 1, 1], [1, 1, 1, 1]])

//     def test_exception_blackout_zone_does_not_raise(self):
//         geo_mag = GeoMag()
//         result = geo_mag.calculate(90, 90, 0, 2020, raise_in_warning_zone=False)
//         self.assertEqual(result.in_blackout_zone, True)

//     def test_exception_blackout_zone_raises(self):
//         geo_mag = GeoMag()
//         with self.assertRaises(BlackoutZoneException):
//             geo_mag.calculate(90, 90, 0, 2020, raise_in_warning_zone=True)

//     def test_exception_caution_zone_does_not_raise(self):
//         geo_mag = GeoMag()
//         result = geo_mag.calculate(80, 80, 0, 2020, raise_in_warning_zone=False)
//         self.assertEqual(result.in_caution_zone, True)

//     def test_exception_caution_zone_raises(self):
//         geo_mag = GeoMag()
//         with self.assertRaises(CautionZoneException):
//             geo_mag.calculate(80, 80, 0, 2020, raise_in_warning_zone=True)

//     def test_get_model_filename_default(self):
//         geo_mag = GeoMag()
//         model_filename = geo_mag._get_model_filename()
//         self.assertEqual(model_filename[-20:], get_os_based_test_path("pygeomag/wmm/WMM.COF"))

//     def test_get_model_filename_default_not_in_wmm_path(self):
//         m = mock_open(read_data="")
//         m.side_effect = [OSError, DEFAULT]
//         with patch("pygeomag.geomag.open", m):
//             geo_mag = GeoMag()
//             model_filename = geo_mag._get_model_filename()
//             self.assertEqual(model_filename[-16:], get_os_based_test_path("pygeomag/WMM.COF"))
//         self.assertEqual(m.call_count, 2)

//     def test_get_model_filename_default_when_neither_file_exists(self):
//         m = mock_open()
//         m.side_effect = OSError
//         with patch("pygeomag.geomag.open", m):
//             geo_mag = GeoMag()
//             model_filename = geo_mag._get_model_filename()
//             self.assertEqual(model_filename[-20:], get_os_based_test_path("pygeomag/wmm/WMM.COF"))
//         self.assertEqual(m.call_count, 2)

//     def test_get_model_filename_different(self):
//         geo_mag = GeoMag(coefficients_file="wmm/WMM_NEW.COF")
//         model_filename = geo_mag._get_model_filename()
//         self.assertEqual(model_filename[-15:], "wmm/WMM_NEW.COF")

//     def test_get_model_filename_fullpath(self):
//         geo_mag = GeoMag(coefficients_file="/wmm/WMM_NEW.COF")
//         model_filename = geo_mag._get_model_filename()
//         self.assertEqual(model_filename, "/wmm/WMM_NEW.COF")

//     def test_load_coefficients_invalid_header(self):
//         geo_mag = GeoMag(coefficients_file="../test/test_files/INVALID_HEADER.COF")
//         with self.assertRaisesRegex(ValueError, "Invalid header in model file"):
//             geo_mag.calculate(0, 80, 0, 2030)

//     def test_load_coefficients_invalid_row(self):
//         geo_mag = GeoMag(coefficients_file="../test/test_files/INVALID_ROW.COF")
//         with self.assertRaisesRegex(ValueError, "Corrupt record in model file"):
//             geo_mag.calculate(0, 80, 0, 2030)

//     def test_load_coefficients_invalid_row_data(self):
//         geo_mag = GeoMag(coefficients_file="../test/test_files/INVALID_ROW_DATA.COF")
//         with self.assertRaisesRegex(ValueError, "Corrupt record in model file"):
//             geo_mag.calculate(0, 80, 0, 2030)

//     def test_load_coefficients_maxord(self):
//         maxord_11_value = -3.4655
//         maxord_12_value = -3.4599
//         self.assertNotEqual(maxord_11_value, maxord_12_value)
//         geo_mag = GeoMag()
//         geo_mag._maxord = 11
//         result = geo_mag.calculate(0, 80, 0, 2020)
//         self.assertAlmostEqual(result.d, maxord_11_value, 4)
//         geo_mag = GeoMag()
//         geo_mag._maxord = 12
//         result = geo_mag.calculate(0, 80, 0, 2020)
//         self.assertAlmostEqual(result.d, maxord_12_value, 4)

//     def test_load_coefficients_missing_file(self):
//         geo_mag = GeoMag(coefficients_file="missing.cof")
//         with self.assertRaisesRegex(FileNotFoundError, "No such file or directory"):
//             geo_mag.calculate(0, 80, 0, 2030)

//     def test_property_life_span(self):
//         geo_mag = GeoMag()
//         self.assertTupleEqual(geo_mag.life_span, (2020.0, 2025.0))

//     def test_property_model(self):
//         geo_mag = GeoMag()
//         self.assertEqual(geo_mag.model, "WMM-2020")

//     def test_property_release_date(self):
//         geo_mag = GeoMag()
//         self.assertEqual(geo_mag.release_date, "12/10/2019")
