package com.ameliaWx.weatherUtils;

import static com.ameliaWx.weatherUtils.PrecipitationType.DRY_SNOW;
import static com.ameliaWx.weatherUtils.PrecipitationType.FREEZING_RAIN;
import static com.ameliaWx.weatherUtils.PrecipitationType.FREEZING_RAIN_SURFACE;
import static com.ameliaWx.weatherUtils.PrecipitationType.FRZR_ICEP_MIX;
import static com.ameliaWx.weatherUtils.PrecipitationType.FRZR_SNOW_MIX;
import static com.ameliaWx.weatherUtils.PrecipitationType.ICEP_SNOW_MIX;
import static com.ameliaWx.weatherUtils.PrecipitationType.ICE_PELLETS;
import static com.ameliaWx.weatherUtils.PrecipitationType.RAIN;
import static com.ameliaWx.weatherUtils.PrecipitationType.RAIN_ICEP_MIX;
import static com.ameliaWx.weatherUtils.PrecipitationType.RAIN_SNOW_MIX;
import static com.ameliaWx.weatherUtils.PrecipitationType.SNOW;
import static com.ameliaWx.weatherUtils.PrecipitationType.VERY_DRY_SNOW;
import static com.ameliaWx.weatherUtils.PrecipitationType.WET_SNOW;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class PtypeAlgorithms {
	private static final double G = 9.81; // m s^-2
	private static final double T0 = 273.15; // Kelvins

	public static void main(String[] args) {
//		double[] pressure = { 70000, 72500, 75000, 77500, 80000, 82500, 85000, 87500, 90000, 92500, 95000, 97500, 100000, 101225 };
//		double[] tmpIsoba = {   263,   270,   273,   275,   283,   280,   274,   272,   268,   263,   263,   266,    268,    270 };
//		double[] hgtIsoba = {  2929,  2656,  2390,  2132,  1881,  1638,  1401,  1172,   948,   731,    518,  268,     32,      2 };
		double[] pressure = { 70000, 72500, 75000, 77500, 80000, 82500, 85000, 87500, 90000, 92500, 95000, 97500,
				100000, 93700 };
		double[] tmpIsoba = { 253, 255, 257, 259, 261, 263, 265, 267, 269, 269, 269, 271, 271, 269 };
		double[] hgtIsoba = { 2929, 2656, 2390, 2132, 1881, 1638, 1401, 1172, 948, 731, 518, 268, 32, 2 };

		PrecipitationType ptype = bourgouinRevisedExtendedMethod(pressure, tmpIsoba, tmpIsoba, hgtIsoba, 101325, 0, 271,
				true);

		System.out.println(ptype);
	}

	public static PrecipitationType freezingPointMethod(double tmp2m) {
		if (tmp2m > 273.15) {
			return RAIN;
		} else {
			return SNOW;
		}
	}

	// make a bourgouin original method too just for shits

	public static PrecipitationType bourgouinRevisedMethod(double[] pressureLevels, double[] tmpIsobaric,
			double[] dptIsobaric, double[] hgtIsobaric, double presSurface, double hgtSurface) {
		return bourgouinRevisedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric, presSurface, hgtSurface,
				false);
	}

	public static PrecipitationType bourgouinRevisedMethod(float[] pressureLevels, float[] tmpIsobaric,
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface) {
		return bourgouinRevisedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric, presSurface, hgtSurface,
				false);
	}

	/**
	 * @param pressureLevels   All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric      All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric      All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric      All isobaric heights of the profile, in Meters
	 * @param dynamicInitLayer true = find highest level with a dewpoint depression
	 *                         of 3 K, false = always init at 500 mb
	 * @return The precipitation type diagnosed by the Bourgouin Revised 2021 Method
	 *         (https://doi.org/10.1175/WAF-D-20-0118.1)
	 */
	public static PrecipitationType bourgouinRevisedMethod(float[] pressureLevels, float[] tmpIsobaric,
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface, boolean dynamicInitLayer) {
//		System.out.println(tmpIsobaric[tmpIsobaric.length - 1]);
		if (tmpIsobaric[tmpIsobaric.length - 1] >= 283.15)
			return RAIN;

		List<RecordAtLevel> isobaricData = new ArrayList<>();
		for (int i = 0; i < pressureLevels.length; i++) {
//			System.out.println("pres: " + pressureLevels[i]);
//			System.out.println("tmp: " + tmpIsobaric[i]);
//			System.out.println("dpt: " + dptIsobaric[i]);
//			System.out.println("wbt: " + WeatherUtils.wetBulbTemperature(tmpIsobaric[i], dptIsobaric[i], pressureLevels[i]));
			isobaricData.add(new RecordAtLevel(pressureLevels[i], tmpIsobaric[i], dptIsobaric[i], hgtIsobaric[i]));
//			System.out.printf("pressure: %4d mb temperature: %6.2f K, wetbulb: %6.2f K, dewpoint: %6.2f K, hgtIsobaric: %6d m\n", (int) pressureLevels[i]/100, tmpIsobaric[i], isobaricData.get(i).wetbulb, dptIsobaric[i], (int) hgtIsobaric[i]);
		}

//		System.out.println("built isobaricData");

		Collections.sort(isobaricData); // should sort top of atmosphere to start of list, surface to bottom of list

		// effectively removes subsurface records
		for (int i = 0; i < pressureLevels.length; i++) {
			double heightAtLevel = isobaricData.get(i).height;

			if (heightAtLevel < hgtSurface) {
				isobaricData.get(i).height = hgtSurface;
				isobaricData.get(i).pressure = presSurface;
				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
			}
		}

//		System.out.println(Arrays.toString(isobaricData.toArray()));

		double initPressureLayer = 50000;

		if (dynamicInitLayer) {
			// looks for first level with a dewpoint depression less than or equal to 3
			// kelvins
			// assumes that's the precip initial layer
			for (int i = 0; i < pressureLevels.length; i++) {
				double dewpointDepression = isobaricData.get(i).temperature - isobaricData.get(i).dewpoint;
//				System.out.println(isobaricData.get(i).pressure/100 + "\t" + dewpointDepression);

				if (dewpointDepression <= 3) {
					initPressureLayer = isobaricData.get(i).pressure;
					break;
				}
			}
		}

//		System.out.println(initPressureLayer);

		int precipInitIndex = -1;
		for (int i = 0; i < pressureLevels.length; i++) {
			if (isobaricData.get(i).pressure >= initPressureLayer) {
				precipInitIndex = i;
				break;
			}
		}
//		System.out.println(precipInitIndex);

		// uses wet-bulb freezing point method if dynamic init layer is used and no
		// precip layers are detected
		if (precipInitIndex == -1)
			return freezingPointMethod(WeatherUtils.wetBulbTemperature(tmpIsobaric[tmpIsobaric.length - 1],
					dptIsobaric[dptIsobaric.length - 1], pressureLevels[pressureLevels.length - 1]));

		// checks if ice nucleation happens
		boolean someIceNucleation = false;
		boolean allIceNucleation = false;
		double probIce = 0; // fraction, NOT PERCENTAGE

		if (tmpIsobaric[precipInitIndex] < 258.15 || !dynamicInitLayer) {
			allIceNucleation = true;
			someIceNucleation = true;
			probIce = 1;
		} else if (tmpIsobaric[precipInitIndex] < 266.15) {
			someIceNucleation = true;

			double t = tmpIsobaric[precipInitIndex] - 273.15; // Celsius
			probIce = (-0.065 * t * t * t * t - 3.1544 * t * t * t - 56.414 * t * t - 449.6 * t - 1308) / 100.0;
		}

		// checks if all wet bulb is freezing, saves time on the vertical profile check
		boolean allFreezing = true;
//		System.out.println(precipInitIndex);
		for (int i = precipInitIndex; i < isobaricData.size(); i++) {
			if (isobaricData.get(i).wetbulb > 273.15) {
				allFreezing = false;
				break;
			}
		}

//		System.out.println(allFreezing);
//		System.out.println(allIceNucleation);
//		System.out.println(someIceNucleation);

		if (allFreezing) {
			if (allIceNucleation) {
				return SNOW;
			} else if (someIceNucleation) {
				return FRZR_SNOW_MIX;
			} else {
				return FREEZING_RAIN_SURFACE;
			}
		}

//		// checks if there's a sufficiently big surface melting layer, saves time on
//		// vertical profile check
//		if (isobaricData.get(isobaricData.size() - 1).wetbulb > 273.15) {
//			double meltingEnergy = 0.0;
//			for (int i = isobaricData.size() - 1; i >= 1; i--) {
//				double heightLower = isobaricData.get(i).height;
//				double heightUpper = isobaricData.get(i - 1).height;
//				double wetbulbLower = isobaricData.get(i).wetbulb;
//				double wetbulbUpper = isobaricData.get(i - 1).wetbulb;
//
//				if (wetbulbUpper > 273.15) {
//					// trapezoidal integration
//					double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//					double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);
//
//					double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
//					double trapezoidHeight = heightUpper - heightLower;
//
//					meltingEnergy += avgTrapezoidWidth * trapezoidHeight;
//				} else {
//					double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);
//
//					// trapezoidal integration (kinda)
//					double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//					double trapezoidWidthUpper = 0;
//
//					double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
//					double trapezoidHeight = freezingPointHeight - heightLower;
//
//					meltingEnergy += avgTrapezoidWidth * trapezoidHeight;
//
//					break;
//				}
//
//				if (meltingEnergy > 13.2) {
//					return RAIN;
//				}
//			}
//		}

		// uses wet-bulb
		double meltingEnergy = 0.0;
		double refreezingEnergy = 0.0;
		double remeltingEnergy = 0.0;

		for (int i = precipInitIndex; i < isobaricData.size() - 1; i++) {

//			System.out.printf("Pressure: %4d, Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, Remelting Energy: %4d J/kg\n", (int) isobaricData.get(i).pressure/100, (int) meltingEnergy, (int) refreezingEnergy, (int) remeltingEnergy);

			// do this tomorrow
			double heightLower = isobaricData.get(i + 1).height;
			double heightUpper = isobaricData.get(i).height;
			double wetbulbLower = isobaricData.get(i + 1).wetbulb;
			double wetbulbUpper = isobaricData.get(i).wetbulb;

			if (wetbulbUpper > 273.15 && wetbulbLower > 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
				double trapezoidHeight = heightUpper - heightLower;

				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;

				if (refreezingEnergy == 0.0) {
					meltingEnergy += Math.abs(trapezoidArea);
				} else {
					remeltingEnergy += Math.abs(trapezoidArea);
				}
			} else if (wetbulbUpper <= 273.15 && wetbulbLower <= 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
				double trapezoidHeight = heightUpper - heightLower;

				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;

				if (meltingEnergy != 0.0) {
					refreezingEnergy += Math.abs(trapezoidArea);
				}
			} else if (wetbulbUpper <= 273.15 && wetbulbLower > 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);

				double avgFreezingTrapezoidWidth = -(0 + trapezoidWidthUpper) / 2;
				double freezingTrapezoidHeight = heightUpper - freezingPointHeight;

				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;

				double avgMeltingTrapezoidWidth = (trapezoidWidthLower + 0) / 2;
				double meltingTrapezoidHeight = freezingPointHeight - heightLower;

				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;

				if (meltingEnergy == 0.0) {
					meltingEnergy += Math.abs(meltingTrapezoidArea);
				} else {
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
					remeltingEnergy += Math.abs(meltingTrapezoidArea);
				}
			} else if (wetbulbUpper > 273.15 && wetbulbLower <= 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);

				double avgFreezingTrapezoidWidth = -(trapezoidWidthLower + 0) / 2;
				double freezingTrapezoidHeight = freezingPointHeight - heightLower;

				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;

				double avgMeltingTrapezoidWidth = (0 + trapezoidWidthUpper) / 2;
				double meltingTrapezoidHeight = heightUpper - freezingPointHeight;

				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;

				if (refreezingEnergy == 0.0) {
					meltingEnergy += Math.abs(meltingTrapezoidArea);
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
				} else {
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
					remeltingEnergy += Math.abs(meltingTrapezoidArea);
				}
			}
		}

//		System.out.printf("Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, Remelting Energy: %4d J/kg\n", (int) meltingEnergy, (int) refreezingEnergy, (int) remeltingEnergy);

		double probSnowI = 15.4 * Math.exp(-0.29 * meltingEnergy); // fraction, NOT PERCENT
		if (probSnowI > 1)
			probSnowI = 1;
		if (probSnowI < 0)
			probSnowI = 0;
		double probSnow = probSnowI * probIce; // fraction, NOT PERCENT

		double probIcep = 0.0;
		if (meltingEnergy > 0) {
			double probIcepI = (2.3 * refreezingEnergy - 42 * Math.log(meltingEnergy + remeltingEnergy + 1) + 3)
					/ 100.0; // fraction, NOT PERCENT
			if (probIcepI > 1)
				probIcepI = 1;
			if (probIcepI < 0)
				probIcepI = 0;
			probIcep = probIcepI * probIce; // fraction, NOT PERCENT
		}

		double probRainI = (-2.1 * refreezingEnergy + 0.2 * (meltingEnergy + remeltingEnergy) + 458) / 100.0;
		if (probRainI > 1)
			probRainI = 1;
		if (probRainI < 0)
			probRainI = 0;
		if (meltingEnergy < 5) {
			probRainI = 0.2 * probRainI; // fraction, NOT PERCENT
		}

		double probRain = (1 - probIce) + probIce * probRainI; // fraction, NOT PERCENT
		if (probRain > 1)
			probRain = 1;

		int dominantType = 0; // 0 = snow, 1 = rain, 2 = ice pellets

		if (probSnow > probRain && probSnow > probIcep) {
			dominantType = 0;
		} else if (probRain > probSnow && probRain > probIcep) {
			dominantType = 1;
		} else if (probIcep > probRain && probIcep > probSnow) {
			dominantType = 2;
		}

		if (dominantType == 0) {
			if (probRain > probIcep && probRain > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_SNOW_MIX;
				} else {
					return FRZR_SNOW_MIX;
				}
			} else if (probIcep > probRain && probIcep > 0.5) {
				return ICEP_SNOW_MIX;
			} else {
				return SNOW;
			}
		} else if (dominantType == 1) {
			if (probSnow > probIcep && probSnow > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_SNOW_MIX;
				} else {
					return FRZR_SNOW_MIX;
				}
			} else if (probIcep > probSnow && probIcep > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_ICEP_MIX;
				} else {
					return FRZR_ICEP_MIX;
				}
			} else {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN;
				} else {
					return FREEZING_RAIN;
				}
			}
		} else if (dominantType == 2) {
			if (probSnow > probRain && probSnow > 0.5) {
				return ICEP_SNOW_MIX;
			} else if (probRain > probSnow && probRain > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_ICEP_MIX;
				} else {
					return FRZR_ICEP_MIX;
				}
			} else {
				return ICE_PELLETS;
			}
		} else {
			return RAIN;
		}
	}

	/**
	 * @param pressureLevels   All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric      All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric      All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric      All isobaric heights of the profile, in Meters
	 * @param dynamicInitLayer true = find highest level with a dewpoint depression
	 *                         of 3 K, false = always init at 500 mb
	 * @return The precipitation type diagnosed by the Bourgouin Revised 2021 Method
	 *         (https://doi.org/10.1175/WAF-D-20-0118.1)
	 */
	public static PrecipitationType bourgouinRevisedMethod(double[] pressureLevels, double[] tmpIsobaric,
			double[] dptIsobaric, double[] hgtIsobaric, double presSurface, double hgtSurface,
			boolean dynamicInitLayer) {
//		System.out.println(tmpIsobaric[tmpIsobaric.length - 1]);
		if (tmpIsobaric[tmpIsobaric.length - 1] >= 283.15)
			return RAIN;

		List<RecordAtLevel> isobaricData = new ArrayList<>();
		for (int i = 0; i < pressureLevels.length; i++) {
//			System.out.println("pres: " + pressureLevels[i]);
//			System.out.println("tmp: " + tmpIsobaric[i]);
//			System.out.println("dpt: " + dptIsobaric[i]);
//			System.out.println("wbt: " + WeatherUtils.wetBulbTemperature(tmpIsobaric[i], dptIsobaric[i], pressureLevels[i]));
			isobaricData.add(new RecordAtLevel(pressureLevels[i], tmpIsobaric[i], dptIsobaric[i], hgtIsobaric[i]));
//			System.out.printf("pressure: %4d mb temperature: %6.2f K, wetbulb: %6.2f K, dewpoint: %6.2f K, hgtIsobaric: %6d m\n", (int) pressureLevels[i]/100, tmpIsobaric[i], isobaricData.get(i).wetbulb, dptIsobaric[i], (int) hgtIsobaric[i]);
		}

//		System.out.println("built isobaricData");

		Collections.sort(isobaricData); // should sort top of atmosphere to start of list, surface to bottom of list

		// effectively removes subsurface records
		for (int i = 0; i < pressureLevels.length; i++) {
			double heightAtLevel = isobaricData.get(i).height;

			if (heightAtLevel < hgtSurface) {
				isobaricData.get(i).height = hgtSurface;
				isobaricData.get(i).pressure = presSurface;
				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
			}
		}

//		System.out.println(Arrays.toString(isobaricData.toArray()));

		double initPressureLayer = 50000;

		if (dynamicInitLayer) {
			// looks for first level with a dewpoint depression less than or equal to 3
			// kelvins
			// assumes that's the precip initial layer
			for (int i = 0; i < pressureLevels.length; i++) {
				double dewpointDepression = isobaricData.get(i).temperature - isobaricData.get(i).dewpoint;
//				System.out.println(isobaricData.get(i).pressure/100 + "\t" + dewpointDepression);

				if (dewpointDepression <= 3) {
					initPressureLayer = isobaricData.get(i).pressure;
					break;
				}
			}
		}

//		System.out.println(initPressureLayer);

		int precipInitIndex = -1;
		for (int i = 0; i < pressureLevels.length; i++) {
			if (isobaricData.get(i).pressure >= initPressureLayer) {
				precipInitIndex = i;
				break;
			}
		}
//		System.out.println(precipInitIndex);

		// uses wet-bulb freezing point method if dynamic init layer is used and no
		// precip layers are detected
		if (precipInitIndex == -1)
			return freezingPointMethod(WeatherUtils.wetBulbTemperature(tmpIsobaric[tmpIsobaric.length - 1],
					dptIsobaric[dptIsobaric.length - 1], pressureLevels[pressureLevels.length - 1]));

		// checks if ice nucleation happens
		boolean someIceNucleation = false;
		boolean allIceNucleation = false;
		double probIce = 0; // fraction, NOT PERCENTAGE

		if (tmpIsobaric[precipInitIndex] < 258.15 || !dynamicInitLayer) {
			allIceNucleation = true;
			someIceNucleation = true;
			probIce = 1;
		} else if (tmpIsobaric[precipInitIndex] < 266.15) {
			someIceNucleation = true;

			double t = tmpIsobaric[precipInitIndex] - 273.15; // Celsius
			probIce = (-0.065 * t * t * t * t - 3.1544 * t * t * t - 56.414 * t * t - 449.6 * t - 1308) / 100.0;
		}

		// checks if all wet bulb is freezing, saves time on the vertical profile check
		boolean allFreezing = true;
//		System.out.println(precipInitIndex);
		for (int i = precipInitIndex; i < isobaricData.size(); i++) {
			if (isobaricData.get(i).wetbulb > 273.15) {
				allFreezing = false;
				break;
			}
		}

//		System.out.println(allFreezing);
//		System.out.println(allIceNucleation);
//		System.out.println(someIceNucleation);

		if (allFreezing) {
			if (allIceNucleation) {
				return SNOW;
			} else if (someIceNucleation) {
				return FRZR_SNOW_MIX;
			} else {
				return FREEZING_RAIN_SURFACE;
			}
		}

//		// checks if there's a sufficiently big surface melting layer, saves time on
//		// vertical profile check
//		if (isobaricData.get(isobaricData.size() - 1).wetbulb > 273.15) {
//			double meltingEnergy = 0.0;
//			for (int i = isobaricData.size() - 1; i >= 1; i--) {
//				double heightLower = isobaricData.get(i).height;
//				double heightUpper = isobaricData.get(i - 1).height;
//				double wetbulbLower = isobaricData.get(i).wetbulb;
//				double wetbulbUpper = isobaricData.get(i - 1).wetbulb;
//
//				if (wetbulbUpper > 273.15) {
//					// trapezoidal integration
//					double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//					double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);
//
//					double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
//					double trapezoidHeight = heightUpper - heightLower;
//
//					meltingEnergy += avgTrapezoidWidth * trapezoidHeight;
//				} else {
//					double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);
//
//					// trapezoidal integration (kinda)
//					double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//					double trapezoidWidthUpper = 0;
//
//					double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
//					double trapezoidHeight = freezingPointHeight - heightLower;
//
//					meltingEnergy += avgTrapezoidWidth * trapezoidHeight;
//
//					break;
//				}
//
//				if (meltingEnergy > 13.2) {
//					return RAIN;
//				}
//			}
//		}

		// uses wet-bulb
		double meltingEnergy = 0.0;
		double refreezingEnergy = 0.0;
		double remeltingEnergy = 0.0;

		for (int i = precipInitIndex; i < isobaricData.size() - 1; i++) {

//			System.out.printf("Pressure: %4d, Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, Remelting Energy: %4d J/kg\n", (int) isobaricData.get(i).pressure/100, (int) meltingEnergy, (int) refreezingEnergy, (int) remeltingEnergy);

			// do this tomorrow
			double heightLower = isobaricData.get(i + 1).height;
			double heightUpper = isobaricData.get(i).height;
			double wetbulbLower = isobaricData.get(i + 1).wetbulb;
			double wetbulbUpper = isobaricData.get(i).wetbulb;

			if (wetbulbUpper > 273.15 && wetbulbLower > 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
				double trapezoidHeight = heightUpper - heightLower;

				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;

				if (refreezingEnergy == 0.0) {
					meltingEnergy += Math.abs(trapezoidArea);
				} else {
					remeltingEnergy += Math.abs(trapezoidArea);
				}
			} else if (wetbulbUpper <= 273.15 && wetbulbLower <= 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
				double trapezoidHeight = heightUpper - heightLower;

				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;

				if (meltingEnergy != 0.0) {
					refreezingEnergy += Math.abs(trapezoidArea);
				}
			} else if (wetbulbUpper <= 273.15 && wetbulbLower > 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);

				double avgFreezingTrapezoidWidth = -(0 + trapezoidWidthUpper) / 2;
				double freezingTrapezoidHeight = heightUpper - freezingPointHeight;

				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;

				double avgMeltingTrapezoidWidth = (trapezoidWidthLower + 0) / 2;
				double meltingTrapezoidHeight = freezingPointHeight - heightLower;

				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;

				if (meltingEnergy == 0.0) {
					meltingEnergy += Math.abs(meltingTrapezoidArea);
				} else {
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
					remeltingEnergy += Math.abs(meltingTrapezoidArea);
				}
			} else if (wetbulbUpper > 273.15 && wetbulbLower <= 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);

				double avgFreezingTrapezoidWidth = -(trapezoidWidthLower + 0) / 2;
				double freezingTrapezoidHeight = freezingPointHeight - heightLower;

				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;

				double avgMeltingTrapezoidWidth = (0 + trapezoidWidthUpper) / 2;
				double meltingTrapezoidHeight = heightUpper - freezingPointHeight;

				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;

				if (refreezingEnergy == 0.0) {
					meltingEnergy += Math.abs(meltingTrapezoidArea);
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
				} else {
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
					remeltingEnergy += Math.abs(meltingTrapezoidArea);
				}
			}
		}

//		System.out.printf("Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, Remelting Energy: %4d J/kg\n", (int) meltingEnergy, (int) refreezingEnergy, (int) remeltingEnergy);

		double probSnowI = 15.4 * Math.exp(-0.29 * meltingEnergy); // fraction, NOT PERCENT
		if (probSnowI > 1)
			probSnowI = 1;
		if (probSnowI < 0)
			probSnowI = 0;
		double probSnow = probSnowI * probIce; // fraction, NOT PERCENT

		double probIcep = 0.0;
		if (meltingEnergy > 0) {
			double probIcepI = (2.3 * refreezingEnergy - 42 * Math.log(meltingEnergy + remeltingEnergy + 1) + 3)
					/ 100.0; // fraction, NOT PERCENT
			if (probIcepI > 1)
				probIcepI = 1;
			if (probIcepI < 0)
				probIcepI = 0;
			probIcep = probIcepI * probIce; // fraction, NOT PERCENT
		}

		double probRainI = (-2.1 * refreezingEnergy + 0.2 * (meltingEnergy + remeltingEnergy) + 458) / 100.0;
		if (probRainI > 1)
			probRainI = 1;
		if (probRainI < 0)
			probRainI = 0;
		if (meltingEnergy < 5) {
			probRainI = 0.2 * probRainI; // fraction, NOT PERCENT
		}

		double probRain = (1 - probIce) + probIce * probRainI; // fraction, NOT PERCENT
		if (probRain > 1)
			probRain = 1;

		int dominantType = 0; // 0 = snow, 1 = rain, 2 = ice pellets

		if (probSnow > probRain && probSnow > probIcep) {
			dominantType = 0;
		} else if (probRain > probSnow && probRain > probIcep) {
			dominantType = 1;
		} else if (probIcep > probRain && probIcep > probSnow) {
			dominantType = 2;
		}

		if (dominantType == 0) {
			if (probRain > probIcep && probRain > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_SNOW_MIX;
				} else {
					return FRZR_SNOW_MIX;
				}
			} else if (probIcep > probRain && probIcep > 0.5) {
				return ICEP_SNOW_MIX;
			} else {
				return SNOW;
			}
		} else if (dominantType == 1) {
			if (probSnow > probIcep && probSnow > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_SNOW_MIX;
				} else {
					return FRZR_SNOW_MIX;
				}
			} else if (probIcep > probSnow && probIcep > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_ICEP_MIX;
				} else {
					return FRZR_ICEP_MIX;
				}
			} else {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN;
				} else {
					return FREEZING_RAIN;
				}
			}
		} else if (dominantType == 2) {
			if (probSnow > probRain && probSnow > 0.5) {
				return ICEP_SNOW_MIX;
			} else if (probRain > probSnow && probRain > 0.5) {
				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
					return RAIN_ICEP_MIX;
				} else {
					return FRZR_ICEP_MIX;
				}
			} else {
				return ICE_PELLETS;
			}
		} else {
			return RAIN;
		}
	}

	/**
	 * @param pressureLevels   All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric      All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric      All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric      All isobaric heights of the profile, in Meters
	 * @param dynamicInitLayer true = find highest level with a dewpoint depression
	 *                         of 3 K, false = always init at 500 mb
	 * @return The precipitation type diagnosed by the Bourgouin Revised 2021 Method
	 *         (https://doi.org/10.1175/WAF-D-20-0118.1), extended by the kuchera
	 *         method and a surface vs. elevated-only check on the freezing rain
	 */
	public static PrecipitationType bourgouinRevisedExtendedMethod(double[] pressureLevels, double[] tmpIsobaric,
			double[] dptIsobaric, double[] hgtIsobaric, double presSurface, double hgtSurface, double tmpSurface,
			boolean dynamicInitLayer) {
		PrecipitationType preExtension = bourgouinRevisedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric,
				presSurface, hgtSurface, true);

//		System.out.println(Arrays.toString(pressureLevels));
//		System.out.println(Arrays.toString(tmpIsobaric));
//		System.out.println(Arrays.toString(dptIsobaric));
//		System.out.println(Arrays.toString(hgtIsobaric));
//		System.out.println(presSurface);
//		System.out.println(hgtSurface);
//		System.out.println(tmpSurface);
//		System.out.println(dynamicInitLayer);

		if (preExtension == PrecipitationType.FREEZING_RAIN) {
			if (tmpSurface <= 273.15) {
				return PrecipitationType.FREEZING_RAIN_SURFACE;
			} else {
				return PrecipitationType.FREEZING_RAIN_ELEVATED;
			}
		}

		if (preExtension == PrecipitationType.SNOW) {
			double kucheraRatio = WeatherUtils.kucheraRatio(tmpIsobaric);

			if (kucheraRatio >= 25) {
				return VERY_DRY_SNOW;
			} else if (kucheraRatio >= 10) {
				return DRY_SNOW;
			} else {
				return WET_SNOW;
			}
		}

		return preExtension;
	}

	/**
	 * @param pressureLevels   All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric      All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric      All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric      All isobaric heights of the profile, in Meters
	 * @param dynamicInitLayer true = find highest level with a dewpoint depression
	 *                         of 3 K, false = always init at 500 mb
	 * @return The precipitation type diagnosed by the Bourgouin Revised 2021 Method
	 *         (https://doi.org/10.1175/WAF-D-20-0118.1), extended by the kuchera
	 *         method and a surface vs. elevated-only check on the freezing rain
	 */
	public static PrecipitationType bourgouinRevisedExtendedMethod(float[] pressureLevels, float[] tmpIsobaric,
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface, float tmpSurface,
			boolean dynamicInitLayer) {
		PrecipitationType preExtension = bourgouinRevisedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric,
				presSurface, hgtSurface, true);

//		System.out.println(Arrays.toString(pressureLevels));
//		System.out.println(Arrays.toString(tmpIsobaric));
//		System.out.println(Arrays.toString(dptIsobaric));
//		System.out.println(Arrays.toString(hgtIsobaric));
//		System.out.println(presSurface);
//		System.out.println(hgtSurface);
//		System.out.println(tmpSurface);
//		System.out.println(dynamicInitLayer);

		if (preExtension == PrecipitationType.FREEZING_RAIN) {
			if (tmpSurface <= 273.15) {
				return PrecipitationType.FREEZING_RAIN_SURFACE;
			} else {
				return PrecipitationType.FREEZING_RAIN_ELEVATED;
			}
		}

		if (preExtension == PrecipitationType.SNOW) {
			double kucheraRatio = WeatherUtils.kucheraRatio(tmpIsobaric);

			if (kucheraRatio >= 25) {
				return VERY_DRY_SNOW;
			} else if (kucheraRatio >= 10) {
				return DRY_SNOW;
			} else {
				return WET_SNOW;
			}
		}

		return preExtension;
	}

	private static double linScale(double preMin, double preMax, double postMin, double postMax, double value) {
		double slope = (postMax - postMin) / (preMax - preMin);

		return slope * (value - preMin) + postMin;
	}
}
