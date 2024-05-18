import 'package:test/test.dart';

import '../src/geomagnetic.dart';


bool test_calculate_uncertainty(){
    GeoMag geo_mag = GeoMag();
    GeoMagResult result = geo_mag.calculate(80, 0, 0, 2020);
    GeoMagUncertaintyResult uncertainty = result.calculate_uncertainty();
    return closeTo(uncertainty.d, 0.01).matches(0.89, {});
}
bool test_d_property(){
    var result = GeoMagResult(DateTime(2020), 0, 0, 0);
    return result.d == result.dec;
}
bool test_f_properties(){
    var result = GeoMagResult(DateTime(2020), 0, 0, 0);
    return (result.f == result.ti && result.f == result.total_intensity);
}
bool test_i_property(){
    var result = GeoMagResult(DateTime(2020), 0, 0, 0);
    return (result.i == result.dip && result.i == result.inclination);
}

bool test_static_values_2015(self){

    var geo_mag = GeoMag(coefficients_file: "../src/wmm/WMM_2015.csv");
    var result = geo_mag.calculate(80, 0, 0, 2015);
    var uncertainty = GeoMagUncertaintyResult(result);
    return (
      uncertainty.north_comp_unc == 138.0 &&
      uncertainty.east_comp_unc == 89.0 &&
      uncertainty.up_comp_unc == 165.0 &&
      uncertainty.h == 133.0 &&
      uncertainty.f == 152.0 &&
      uncertainty.i == 0.22
    );
}

//     def test_static_values_2020(self):
//         geo_mag = GeoMag()
//         result = geo_mag.calculate(80, 0, 0, 2020)
//         uncertainty = GeoMagUncertaintyResult(result)
//         self.assertEqual(uncertainty.x, 131.0)
//         self.assertEqual(uncertainty.y, 94.0)
//         self.assertEqual(uncertainty.z, 157.0)
//         self.assertEqual(uncertainty.h, 128.0)
//         self.assertEqual(uncertainty.f, 148.0)
//         self.assertEqual(uncertainty.i, 0.21)

//     def test_uncertainty_degrees_2015(self):
//         geo_mag = GeoMag(coefficients_file="wmm/WMM_2015.COF")
//         result = geo_mag.calculate(80, 0, 0, 2015.0)
//         uncertainty = GeoMagUncertaintyResult(result)
//         self.assertAlmostEqual(uncertainty.d, 0.85, 2)

//         result = geo_mag.calculate(0, 120, 0, 2015.0)
//         uncertainty = GeoMagUncertaintyResult(result)
//         self.assertAlmostEqual(uncertainty.d, 0.27, 2)

//     def test_uncertainty_degrees_2022(self):
//         geo_mag = GeoMag()
//         result = geo_mag.calculate(80, 0, 0, 2020.0)
//         uncertainty = GeoMagUncertaintyResult(result)
//         self.assertAlmostEqual(uncertainty.d, 0.89, 2)

//         result = geo_mag.calculate(0, 120, 0, 2020.0)
//         uncertainty = GeoMagUncertaintyResult(result)
//         self.assertAlmostEqual(uncertainty.d, 0.30, 2)

//     def test_time_out_of_supported_range(self):
//         geo_mag = GeoMag()
//         result = geo_mag.calculate(80, 0, 0, 2026.0, allow_date_outside_lifespan=True)
//         with self.assertRaisesRegex(ValueError, "GeoMagResult outside of known uncertainty estimates."):
//             GeoMagUncertaintyResult(result)

void main() {
  
  test('test_calculate_uncertainty', test_calculate_uncertainty);
  test('test_d_property', test_d_property);
  test('test_f_properties', test_f_properties);
  test('test_i_property', test_i_property);



}

// class TestGeoMagCoefficients(TestCase):
//     def get_test_values(self, test_parameter, style):
//         if style == 0:
//             print(len(test_parameter.split()))
//             time, alt, glat, glon, d, i, h, x, y, z, f, gv, _, _, _, _, _, _, _ = (
//                 t(s)
//                 for t, s in zip(
//                     (
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                     ),
//                     test_parameter.split(),
//                 )
//             )

//         elif style == 1:
//             time, alt, glat, glon, x, y, z, h, f, i, d, gv, _, _, _, _, _, _, _ = (
//                 t(s)
//                 for t, s in zip(
//                     (
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                     ),
//                     test_parameter.split(),
//                 )
//             )

//         else:
//             time, alt, glat, glon, d, i, h, x, y, z, f, _, _, _, _, _, _, _ = (
//                 t(s)
//                 for t, s in zip(
//                     (
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                         float,
//                     ),
//                     test_parameter.split(),
//                 )
//             )
//             gv = None

//         return time, alt, glat, glon, x, y, z, h, f, i, d, gv

//     def run_tests(self, geo_mag, test_filename, style):
//         with open(get_test_filename(test_filename)) as test_values_file:
//             for row, test_parameter in enumerate(test_values_file):
//                 if test_parameter[0] == "#":
//                     continue

//                 time, alt, glat, glon, x, y, z, h, f, i, d, gv = self.get_test_values(test_parameter, style)

//                 result = geo_mag.calculate(glat, glon, alt, time)
//                 gv_test = -999 if result.gv is None else result.gv

//                 self.assertAlmostEqual(x, result.x, 1, f"Row {row}: X (nT) expected {x}, result {result.x}")
//                 self.assertAlmostEqual(y, result.y, 1, f"Row {row}: Y (nT) expected {y}, result {result.y}")
//                 self.assertAlmostEqual(z, result.z, 1, f"Row {row}: Z (nT) expected {z}, result {result.z}")
//                 self.assertAlmostEqual(h, result.h, 1, f"Row {row}: H (nT) expected {h}, result {result.h}")
//                 self.assertAlmostEqual(f, result.f, 1, f"Row {row}: F (nT) expected {f}, result {result.f}")
//                 self.assertAlmostEqual(i, result.i, 2, f"Row {row}: I (Deg) expected {i}, result {result.i}")
//                 self.assertAlmostEqual(d, result.d, 2, f"Row {row}: D (Deg) expected {d}, result {result.d}")
//                 if style != 2:
//                     self.assertAlmostEqual(gv, gv_test, 2, f"Row {row}: GV (Deg) expected {gv}, result {result.gv}")

//     def test_calculate_declination_from_2010_wmm_style_0_file(self):
//         self.run_tests(GeoMag(coefficients_file="wmm/WMM_2010.COF"), "test_values/WMM2010testvalues.txt", 0)

//     def test_calculate_declination_from_2015_wmm_style_1_file(self):
//         self.run_tests(GeoMag(coefficients_file="wmm/WMM_2015.COF"), "test_values/WMM2015testvalues.txt", 1)

//     def test_calculate_declination_from_2015_wmm_style_1_data(self):
//         self.run_tests(GeoMag(coefficients_data=WMM_2015), "test_values/WMM2015testvalues.txt", 1)

//     def test_calculate_declination_from_2015v2_wmm_style_1_file(self):
//         self.run_tests(GeoMag(coefficients_file="wmm/WMM_2015v2.COF"), "test_values/WMM2015v2testvalues.txt", 1)

//     def test_calculate_declination_from_2015v2_wmm_style_1_data(self):
//         self.run_tests(GeoMag(coefficients_data=WMM_2015v2), "test_values/WMM2015v2testvalues.txt", 1)

//     def test_calculate_declination_from_2020_wmm_style_1_file(self):
//         self.run_tests(GeoMag(), "test_values/WMM2020testvalues.txt", 1)

//     def test_calculate_declination_from_2020_wmm_style_1_data(self):
//         self.run_tests(GeoMag(coefficients_data=WMM_2020), "test_values/WMM2020testvalues.txt", 1)

//     def test_calculate_declination_from_2020_wmm_style_2(self):
//         self.run_tests(GeoMag(), "test_values/WMM2020_TEST_VALUES.txt", 2)


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