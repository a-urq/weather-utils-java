package com.ameliaWx.weatherUtils;

import static com.ameliaWx.weatherUtils.PrecipitationType.DRY_SNOW;
import static com.ameliaWx.weatherUtils.PrecipitationType.FREEZING_RAIN;
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

	private static final double SATURATION_THRESHOLD_RH = 0.8; // Fraction

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

	public static double meltingEnergy(double[] pressureLevels, double[] tmpIsobaric, double[] dptIsobaric,
			double[] hgtIsobaric, double presSurface, double hgtSurface) {
		List<RecordAtLevel> isobaricData = new ArrayList<>();
		for (int i = 0; i < dptIsobaric.length; i++) {
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
		for (int i = 0; i < dptIsobaric.length; i++) {
			double heightAtLevel = isobaricData.get(i).height;

			if (heightAtLevel < hgtSurface) {
				isobaricData.get(i).height = hgtSurface;
				isobaricData.get(i).pressure = presSurface;
				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
			}
		}

		// patch that removes 0 records in middle of profile
		for (int i = 0; i < isobaricData.size() - 1; i++) {
			double heightAtLevel = isobaricData.get(i).height - isobaricData.get(isobaricData.size() - 1).height;
			double heightAtLevelNext = isobaricData.get(i + 1).height - isobaricData.get(isobaricData.size() - 1).height;

			if (heightAtLevel < 0.1 && heightAtLevelNext > 1.9) {
				isobaricData.remove(i);
				i--; // ensures no entries are skipped after record removal
			}
		}

//		System.out.println(Arrays.toString(isobaricData.toArray()));

		double initPressureLayer = 50000;

		// looks for first level with RH >= 90
		for (int i = 0; i < dptIsobaric.length; i++) {
			double rh = WeatherUtils.relativeHumidity(isobaricData.get(i).temperature, isobaricData.get(i).dewpoint);

			if (rh >= SATURATION_THRESHOLD_RH) {
				initPressureLayer = isobaricData.get(i).pressure;
				break;
			}
		}

//		System.out.println(initPressureLayer);

		int precipInitIndex = -1;
		for (int i = 0; i < dptIsobaric.length; i++) {
			if (isobaricData.get(i).pressure >= initPressureLayer) {
				precipInitIndex = i;
				break;
			}
		}

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

		return meltingEnergy + remeltingEnergy;
	}

	public static double refreezingEnergy(double[] pressureLevels, double[] tmpIsobaric, double[] dptIsobaric,
			double[] hgtIsobaric, double presSurface, double hgtSurface) {
		List<RecordAtLevel> isobaricData = new ArrayList<>();
		for (int i = 0; i < dptIsobaric.length; i++) {
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
		for (int i = 0; i < dptIsobaric.length; i++) {
			double heightAtLevel = isobaricData.get(i).height;

			if (heightAtLevel < hgtSurface) {
				isobaricData.get(i).height = hgtSurface;
				isobaricData.get(i).pressure = presSurface;
				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
			}
		}

		// patch that removes 0 records in middle of profile
		for (int i = 0; i < isobaricData.size() - 1; i++) {
			double heightAtLevel = isobaricData.get(i).height - isobaricData.get(isobaricData.size() - 1).height;
			double heightAtLevelNext = isobaricData.get(i + 1).height - isobaricData.get(isobaricData.size() - 1).height;

			if (heightAtLevel < 0.1 && heightAtLevelNext > 1.9) {
				isobaricData.remove(i);
				i--; // ensures no entries are skipped after record removal
			}
		}

//		System.out.println(Arrays.toString(isobaricData.toArray()));

		double initPressureLayer = 50000;

		// looks for first level with RH >= 90
		for (int i = 0; i < dptIsobaric.length; i++) {
			double rh = WeatherUtils.relativeHumidity(isobaricData.get(i).temperature, isobaricData.get(i).dewpoint);

			if (rh >= SATURATION_THRESHOLD_RH) {
				initPressureLayer = isobaricData.get(i).pressure;
				break;
			}
		}

//		System.out.println(initPressureLayer);

		int precipInitIndex = -1;
		for (int i = 0; i < dptIsobaric.length; i++) {
			if (isobaricData.get(i).pressure >= initPressureLayer) {
				precipInitIndex = i;
				break;
			}
		}

		double meltingEnergy = 0.0;
		double refreezingEnergy = 0.0;
		double remeltingEnergy = 0.0;

		for (int i = precipInitIndex; i < isobaricData.size() - 1; i++) {
			// do this tomorrow
			double heightLower = isobaricData.get(i + 1).height;
			double heightUpper = isobaricData.get(i).height;
			double wetbulbLower = isobaricData.get(i + 1).wetbulb;
			double wetbulbUpper = isobaricData.get(i).wetbulb;

			if (wetbulbUpper > 273.15 && wetbulbLower > 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);

				double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
				double trapezoidHeight = heightUpper - heightLower;

				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;

				if (refreezingEnergy == 0.0) {
					meltingEnergy += Math.abs(trapezoidArea);
				} else {
					remeltingEnergy += Math.abs(trapezoidArea);
				}

				System.out.printf("%5.1f\t%5.1f\t%5d\t%5d\t", (wetbulbUpper - 273.15), (wetbulbLower - 273.15),
						(int) heightUpper, (int) heightLower);
				System.out.printf("%08.1fM/%8s ", trapezoidArea, "");
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

				System.out.printf("%5.1f\t%5.1f\t%5d\t%5d\t", (wetbulbUpper - 273.15), (wetbulbLower - 273.15),
						(int) heightUpper, (int) heightLower);
				System.out.printf("%08.1fF/%8s ", trapezoidArea, "");
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

				System.out.printf("%5.1f\t%5.1f\t%5d\t%5d\t", (wetbulbUpper - 273.15), (wetbulbLower - 273.15),
						(int) heightUpper, (int) heightLower);
				System.out.printf("%08.1fF/%08.1fM ", freezingTrapezoidArea, meltingTrapezoidArea);
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

				System.out.printf("%5.1f\t%5.1f\t%5d\t%5d\t", (wetbulbUpper - 273.15), (wetbulbLower - 273.15),
						(int) heightUpper, (int) heightLower);
				System.out.printf("%08.1fM/%08.1fF ", meltingTrapezoidArea, freezingTrapezoidArea);
			}

			System.out.printf("Pressure: %4d, Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, T_w: %4.1f C\n",
					(int) isobaricData.get(i).pressure / 100, (int) (meltingEnergy + remeltingEnergy),
					(int) refreezingEnergy, (wetbulbLower - 273.15));
		}

		return refreezingEnergy;
	}

	/**
	 * @param pressureLevels   All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric      All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric      All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric      All isobaric heights of the profile, in Meters
	 * @param dynamicInitLayer true = find highest level with a dewpoint depression
	 *                         of 3 K, false = always init at 500 mb
	 * @return The precipitation type diagnosed by the Bourgouin 2000 Method
	 *         (http://employees.oneonta.edu/blechmjb/JBpages/M361/BUFKITlabs/Bourgouin/Bourgouin2000.pdf)
	 */
	public static PrecipitationType bourgouinMethod(double[] pressureLevels, double[] tmpIsobaric, double[] dptIsobaric,
			double[] hgtIsobaric, double presSurface, double hgtSurface, boolean dynamicInitLayer) {
//		System.out.println("start bourgouin 2000");

		if (tmpIsobaric[tmpIsobaric.length - 1] >= 283.15)
			return RAIN;

		List<RecordAtLevel> isobaricData = new ArrayList<>();
		for (int i = 0; i < dptIsobaric.length; i++) {
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
		for (int i = 0; i < dptIsobaric.length; i++) {
			double heightAtLevel = isobaricData.get(i).height;

			if (heightAtLevel < hgtSurface) {
				isobaricData.get(i).height = hgtSurface;
				isobaricData.get(i).pressure = presSurface;
				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
			}
		}

		// patch that removes 0 records in middle of profile
		for (int i = 0; i < isobaricData.size() - 1; i++) {
			double heightAtLevel = isobaricData.get(i).height - isobaricData.get(isobaricData.size() - 1).height;
			double heightAtLevelNext = isobaricData.get(i + 1).height - isobaricData.get(isobaricData.size() - 1).height;

			if (heightAtLevel < 0.1 && heightAtLevelNext > 1.9) {
				isobaricData.remove(i);
				i--; // ensures no entries are skipped after record removal
			}
		}

//		System.out.println(Arrays.toString(isobaricData.toArray()));

		double initPressureLayer = 50000;

		if (dynamicInitLayer) {
			// looks for first level with a dewpoint depression less than or equal to 3
			// kelvins
			// assumes that's the precip initial layer
			for (int i = 0; i < isobaricData.size(); i++) {
//				System.out.println(isobaricData.get(i).pressure/100 + "\t" + dewpointDepression);

				double rh = WeatherUtils.relativeHumidity(isobaricData.get(i).temperature,
						isobaricData.get(i).dewpoint);

				if (rh >= SATURATION_THRESHOLD_RH) {
					initPressureLayer = isobaricData.get(i).pressure;
					break;
				}
			}
		}

//		System.out.println(initPressureLayer);

		int precipInitIndex = -1;
		for (int i = 0; i < dptIsobaric.length; i++) {
			if (isobaricData.get(i).pressure >= initPressureLayer) {
				precipInitIndex = i;
				break;
			}
		}
//		System.out.println(precipInitIndex);

		// uses freezing point method if dynamic init layer is used and no
		// precip layers are detected
		if (precipInitIndex == -1 || isobaricData.get(precipInitIndex).height - 3 < hgtSurface)
			return freezingPointMethod(tmpIsobaric[tmpIsobaric.length - 1]);

		// checks if all wet bulb is freezing, saves time on the vertical profile check
		boolean allFreezing = true;
		for (int i = precipInitIndex; i < isobaricData.size(); i++) {
			if (isobaricData.get(i).wetbulb > 273.15) {
				allFreezing = false;
				break;
			}
		}

		if (allFreezing) {
			return SNOW;
		}

		// checks if there's a sufficiently big surface melting layer, saves time on
		// vertical profile check
//		if (isobaricData.get(isobaricData.size() - 1).temperature > 273.15) {
//			double meltingEnergy = 0.0;
//			for (int i = isobaricData.size() - 1; i >= 1; i--) {
//				double heightLower = isobaricData.get(i).height;
//				double heightUpper = isobaricData.get(i - 1).height;
//				double temperatureLower = isobaricData.get(i).temperature;
//				double temperatureUpper = isobaricData.get(i - 1).temperature;
//
//				if (temperatureUpper > 273.15) {
//					// trapezoidal integration
//					double trapezoidWidthLower = G * ((temperatureLower - T0) / T0);
//					double trapezoidWidthUpper = G * ((temperatureUpper - T0) / T0);
//
//					double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
//					double trapezoidHeight = heightUpper - heightLower;
//
//					meltingEnergy += avgTrapezoidWidth * trapezoidHeight;
//				} else {
//					double freezingPointHeight = linScale(temperatureLower, temperatureUpper, heightLower, heightUpper, T0);
//
//					// trapezoidal integration (kinda)
//					double trapezoidWidthLower = G * ((temperatureLower - T0) / T0);
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
//				System.out.println("bourgouin ME:" + meltingEnergy);
//
//				if (meltingEnergy > 13.2) {
//					return RAIN;
//				} else if (meltingEnergy >= 5.6) {
//					return RAIN_SNOW_MIX;
//				}
//			}
//		}

		// uses wet-bulb
		double meltingEnergy = 0.0;
		double refreezingEnergy = 0.0;
		double remeltingEnergy = 0.0;

		for (int i = precipInitIndex; i < isobaricData.size() - 1; i++) {
			double heightLower = isobaricData.get(i + 1).height;
			double heightUpper = isobaricData.get(i).height;
			double temperatureLower = isobaricData.get(i + 1).temperature;
			double temperatureUpper = isobaricData.get(i).temperature;

			if (temperatureUpper > 273.15 && temperatureLower > 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((temperatureLower - T0) / T0);
				double trapezoidWidthUpper = G * ((temperatureUpper - T0) / T0);

				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
				double trapezoidHeight = heightUpper - heightLower;

				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;

//				System.out.println("TH          :  " + temperatureLower);
//				System.out.println("TH          :  " + temperatureUpper);
//				System.out.println("TH          :  " + heightLower);
//				System.out.println("TH          :  " + heightUpper);
//				System.out.println("TH          :  " + trapezoidWidthLower);
//				System.out.println("TH          :  " + trapezoidWidthUpper);
//				System.out.println("TH          :  " + avgTrapezoidWidth);
//				System.out.println("TH          :  " + trapezoidHeight);
//				System.out.println("PA trapezoid:  " + trapezoidArea);
//				System.out.println("bourgouin ME:  " + meltingEnergy);
//				System.out.println("bourgouin RFE: " + refreezingEnergy);
//				System.out.println("bourgouin RME: " + remeltingEnergy);

				if (refreezingEnergy == 0.0) {
					meltingEnergy += Math.abs(trapezoidArea);
				} else {
					remeltingEnergy += Math.abs(trapezoidArea);
				}
			} else if (temperatureUpper <= 273.15 && temperatureLower <= 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((temperatureLower - T0) / T0);
				double trapezoidWidthUpper = G * ((temperatureUpper - T0) / T0);

				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
				double trapezoidHeight = heightUpper - heightLower;

				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;

//				System.out.println("NA trapezoid:  " + trapezoidArea);
//				System.out.println("bourgouin ME:  " + meltingEnergy);
//				System.out.println("bourgouin RFE: " + refreezingEnergy);
//				System.out.println("bourgouin RME: " + remeltingEnergy);

				if (meltingEnergy != 0.0) {
					refreezingEnergy += Math.abs(trapezoidArea);
				}
			} else if (temperatureUpper <= 273.15 && temperatureLower > 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((temperatureLower - T0) / T0);
				double trapezoidWidthUpper = G * ((temperatureUpper - T0) / T0);

				double freezingPointHeight = linScale(temperatureLower, temperatureUpper, heightLower, heightUpper, T0);

				double avgFreezingTrapezoidWidth = -(0 + trapezoidWidthUpper) / 2;
				double freezingTrapezoidHeight = heightUpper - freezingPointHeight;

				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;

				double avgMeltingTrapezoidWidth = (trapezoidWidthLower + 0) / 2;
				double meltingTrapezoidHeight = freezingPointHeight - heightLower;

				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;

//				System.out.println("NA-> trapezoid:  " + freezingTrapezoidArea);
//				System.out.println("->PA trapezoid:  " + meltingTrapezoidArea);
//				System.out.println("bourgouin ME:    " + meltingEnergy);
//				System.out.println("bourgouin RFE:   " + refreezingEnergy);
//				System.out.println("bourgouin RME    " + remeltingEnergy);

				if (meltingEnergy == 0.0) {
					meltingEnergy += Math.abs(meltingTrapezoidArea);
				} else {
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
					remeltingEnergy += Math.abs(meltingTrapezoidArea);
				}
			} else if (temperatureUpper > 273.15 && temperatureLower <= 273.15) {
				// trapezoidal integration
				double trapezoidWidthLower = G * ((temperatureLower - T0) / T0);
				double trapezoidWidthUpper = G * ((temperatureUpper - T0) / T0);

				double freezingPointHeight = linScale(temperatureLower, temperatureUpper, heightLower, heightUpper, T0);

				double avgFreezingTrapezoidWidth = -(trapezoidWidthLower + 0) / 2;
				double freezingTrapezoidHeight = freezingPointHeight - heightLower;

				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;

				double avgMeltingTrapezoidWidth = (0 + trapezoidWidthUpper) / 2;
				double meltingTrapezoidHeight = heightUpper - freezingPointHeight;

				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;

//				System.out.println("PA-> trapezoid:  " + meltingTrapezoidArea);
//				System.out.println("->NA trapezoid:  " + freezingTrapezoidArea);
//				System.out.println("bourgouin ME:    " + meltingEnergy);
//				System.out.println("bourgouin RFE:   " + refreezingEnergy);
//				System.out.println("bourgouin RME    " + remeltingEnergy);

				if (refreezingEnergy == 0.0) {
					meltingEnergy += Math.abs(meltingTrapezoidArea);
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
				} else {
					refreezingEnergy += Math.abs(freezingTrapezoidArea);
					remeltingEnergy += Math.abs(meltingTrapezoidArea);
				}
			}
		}

//		System.out.println("bourgouin ME:" + meltingEnergy);
//		System.out.println("bourgouin RFE:" + refreezingEnergy);
//		System.out.println("bourgouin RME:" + remeltingEnergy);

		if (remeltingEnergy >= 5.6) {
			if (remeltingEnergy > 13.2) {
				return RAIN;
			} else {
				return RAIN_SNOW_MIX;
			}
		} else if (refreezingEnergy > 0) {
			if (refreezingEnergy > 66 + 0.66 * meltingEnergy) {
				return ICE_PELLETS;
			} else if (refreezingEnergy < 46 + 0.66 * meltingEnergy) {
				return FREEZING_RAIN;
			} else {
				return FRZR_ICEP_MIX;
			}
		} else {
			if (meltingEnergy > 13.2) {
				return RAIN;
			}
			if (meltingEnergy >= 5.6) {
				return RAIN_SNOW_MIX;
			} else {
				return SNOW;
			}
		}
	}

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

	public static PrecipitationType bourgouinRevisedMethod(float[] pressureLevels, float[] tmpIsobaric,
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface, boolean dynamicInitLayer) {
		return bourgouinRevisedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric, presSurface, hgtSurface, dynamicInitLayer, false);
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
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface, boolean dynamicInitLayer, boolean debug) {
		double[] pressureLevelsD = convFloatToDouble(pressureLevels);
		double[] tmpIsobaricD = convFloatToDouble(tmpIsobaric);
		double[] dptIsobaricD = convFloatToDouble(dptIsobaric);
		double[] hgtIsobaricD = convFloatToDouble(hgtIsobaric);

		if(debug) System.out.println("ptype debug on");
		return bourgouinRevisedMethod(pressureLevelsD, tmpIsobaricD, dptIsobaricD, hgtIsobaricD, presSurface,
				hgtSurface, dynamicInitLayer, debug);

////		System.out.println(tmpIsobaric[tmpIsobaric.length - 1]);
//		if (tmpIsobaric[tmpIsobaric.length - 1] >= 283.15)
//			return RAIN;
//
//		List<RecordAtLevel> isobaricData = new ArrayList<>();
//		for (int i = 0; i < dptIsobaric.length; i++) {
////			System.out.println("pres: " + pressureLevels[i]);
////			System.out.println("tmp: " + tmpIsobaric[i]);
////			System.out.println("dpt: " + dptIsobaric[i]);
////			System.out.println("wbt: " + WeatherUtils.wetBulbTemperature(tmpIsobaric[i], dptIsobaric[i], pressureLevels[i]));
//			isobaricData.add(new RecordAtLevel(pressureLevels[i], tmpIsobaric[i], dptIsobaric[i], hgtIsobaric[i]));
////			System.out.printf("pressure: %4d mb temperature: %6.2f K, wetbulb: %6.2f K, dewpoint: %6.2f K, hgtIsobaric: %6d m\n", (int) pressureLevels[i]/100, tmpIsobaric[i], isobaricData.get(i).wetbulb, dptIsobaric[i], (int) hgtIsobaric[i]);
//		}
//
////		System.out.println("built isobaricData");
//
//		Collections.sort(isobaricData); // should sort top of atmosphere to start of list, surface to bottom of list
//
//		// effectively removes subsurface records
//		for (int i = 0; i < dptIsobaric.length; i++) {
//			double heightAtLevel = isobaricData.get(i).height;
//
//			if (heightAtLevel < hgtSurface) {
//				isobaricData.get(i).height = hgtSurface;
//				isobaricData.get(i).pressure = presSurface;
//				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
//				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
//				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
//			}
//		}
//
////		System.out.println(Arrays.toString(isobaricData.toArray()));
//
//		double initPressureLayer = 50000;
//
//		if (dynamicInitLayer) {
//			// looks for first level with RH >= 90
//			for (int i = 0; i < dptIsobaric.length; i++) {
//				double rh = WeatherUtils.relativeHumidity(isobaricData.get(i).temperature, isobaricData.get(i).dewpoint);
//				
//				if (rh >= SATURATION_THRESHOLD_RH) {
//					initPressureLayer = isobaricData.get(i).pressure;
//					break;
//				}
//			}
//		}
//
////		System.out.println(initPressureLayer);
//
//		int precipInitIndex = -1;
//		for (int i = 0; i < dptIsobaric.length; i++) {
//			if (isobaricData.get(i).pressure >= initPressureLayer) {
//				precipInitIndex = i;
//				break;
//			}
//		}
////		System.out.println(precipInitIndex);
//
//		// uses wet-bulb freezing point method if dynamic init layer is used and no
//		// precip layers are detected
//		if (precipInitIndex == -1 || isobaricData.get(precipInitIndex).height - 3 < hgtSurface)
//			return freezingPointMethod(WeatherUtils.wetBulbTemperature(tmpIsobaric[tmpIsobaric.length - 1],
//					dptIsobaric[dptIsobaric.length - 1], pressureLevels[pressureLevels.length - 1]));
//
//		// checks if ice nucleation happens
//		boolean someIceNucleation = false;
//		boolean allIceNucleation = false;
//		double probIce = 0; // fraction, NOT PERCENTAGE
//
//		if (tmpIsobaric[precipInitIndex] < 258.15 || !dynamicInitLayer) {
//			allIceNucleation = true;
//			someIceNucleation = true;
//			probIce = 1;
//		} else if (tmpIsobaric[precipInitIndex] < 266.15) {
//			someIceNucleation = true;
//
//			double t = tmpIsobaric[precipInitIndex] - 273.15; // Celsius
//			probIce = (-0.065 * t * t * t * t - 3.1544 * t * t * t - 56.414 * t * t - 449.6 * t - 1308) / 100.0;
//		}
//
//		// checks if all wet bulb is freezing, saves time on the vertical profile check
//		boolean allFreezing = true;
////		System.out.println(precipInitIndex);
//		for (int i = precipInitIndex; i < isobaricData.size(); i++) {
//			if (isobaricData.get(i).wetbulb > 273.15) {
//				allFreezing = false;
//				break;
//			}
//		}
//
////		System.out.println(allFreezing);
////		System.out.println(allIceNucleation);
////		System.out.println(someIceNucleation);
//
//		if (allFreezing) {
//			if (allIceNucleation) {
//				return SNOW;
//			} else if (someIceNucleation) {
//				return FRZR_SNOW_MIX;
//			} else {
//				return FREEZING_RAIN;
//			}
//		}
//
////		// checks if there's a sufficiently big surface melting layer, saves time on
////		// vertical profile check
////		if (isobaricData.get(isobaricData.size() - 1).wetbulb > 273.15) {
////			double meltingEnergy = 0.0;
////			for (int i = isobaricData.size() - 1; i >= 1; i--) {
////				double heightLower = isobaricData.get(i).height;
////				double heightUpper = isobaricData.get(i - 1).height;
////				double wetbulbLower = isobaricData.get(i).wetbulb;
////				double wetbulbUpper = isobaricData.get(i - 1).wetbulb;
////
////				if (wetbulbUpper > 273.15) {
////					// trapezoidal integration
////					double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
////					double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);
////
////					double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
////					double trapezoidHeight = heightUpper - heightLower;
////
////					meltingEnergy += avgTrapezoidWidth * trapezoidHeight;
////				} else {
////					double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);
////
////					// trapezoidal integration (kinda)
////					double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
////					double trapezoidWidthUpper = 0;
////
////					double avgTrapezoidWidth = (trapezoidWidthLower + trapezoidWidthUpper) / 2;
////					double trapezoidHeight = freezingPointHeight - heightLower;
////
////					meltingEnergy += avgTrapezoidWidth * trapezoidHeight;
////
////					break;
////				}
////
////				if (meltingEnergy > 13.2) {
////					return RAIN;
////				}
////			}
////		}
//
//		// uses wet-bulb
//		double meltingEnergy = 0.0;
//		double refreezingEnergy = 0.0;
//		double remeltingEnergy = 0.0;
//
//		for (int i = precipInitIndex; i < isobaricData.size() - 1; i++) {
//
////			System.out.printf("Pressure: %4d, Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, Remelting Energy: %4d J/kg\n", (int) isobaricData.get(i).pressure/100, (int) meltingEnergy, (int) refreezingEnergy, (int) remeltingEnergy);
//
//			// do this tomorrow
//			double heightLower = isobaricData.get(i + 1).height;
//			double heightUpper = isobaricData.get(i).height;
//			double wetbulbLower = isobaricData.get(i + 1).wetbulb;
//			double wetbulbUpper = isobaricData.get(i).wetbulb;
//
//			if (wetbulbUpper > 273.15 && wetbulbLower > 273.15) {
//				// trapezoidal integration
//				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);
//
//				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
//				double trapezoidHeight = heightUpper - heightLower;
//
//				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;
//
//				if (refreezingEnergy == 0.0) {
//					meltingEnergy += Math.abs(trapezoidArea);
//				} else {
//					remeltingEnergy += Math.abs(trapezoidArea);
//				}
//			} else if (wetbulbUpper <= 273.15 && wetbulbLower <= 273.15) {
//				// trapezoidal integration
//				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);
//
//				double avgTrapezoidWidth = -(trapezoidWidthLower + trapezoidWidthUpper) / 2;
//				double trapezoidHeight = heightUpper - heightLower;
//
//				double trapezoidArea = avgTrapezoidWidth * trapezoidHeight;
//
//				if (meltingEnergy != 0.0) {
//					refreezingEnergy += Math.abs(trapezoidArea);
//				}
//			} else if (wetbulbUpper <= 273.15 && wetbulbLower > 273.15) {
//				// trapezoidal integration
//				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);
//
//				double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);
//
//				double avgFreezingTrapezoidWidth = -(0 + trapezoidWidthUpper) / 2;
//				double freezingTrapezoidHeight = heightUpper - freezingPointHeight;
//
//				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;
//
//				double avgMeltingTrapezoidWidth = (trapezoidWidthLower + 0) / 2;
//				double meltingTrapezoidHeight = freezingPointHeight - heightLower;
//
//				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;
//
//				if (meltingEnergy == 0.0) {
//					meltingEnergy += Math.abs(meltingTrapezoidArea);
//				} else {
//					refreezingEnergy += Math.abs(freezingTrapezoidArea);
//					remeltingEnergy += Math.abs(meltingTrapezoidArea);
//				}
//			} else if (wetbulbUpper > 273.15 && wetbulbLower <= 273.15) {
//				// trapezoidal integration
//				double trapezoidWidthLower = G * ((wetbulbLower - T0) / T0);
//				double trapezoidWidthUpper = G * ((wetbulbUpper - T0) / T0);
//
//				double freezingPointHeight = linScale(wetbulbLower, wetbulbUpper, heightLower, heightUpper, T0);
//
//				double avgFreezingTrapezoidWidth = -(trapezoidWidthLower + 0) / 2;
//				double freezingTrapezoidHeight = freezingPointHeight - heightLower;
//
//				double freezingTrapezoidArea = avgFreezingTrapezoidWidth * freezingTrapezoidHeight;
//
//				double avgMeltingTrapezoidWidth = (0 + trapezoidWidthUpper) / 2;
//				double meltingTrapezoidHeight = heightUpper - freezingPointHeight;
//
//				double meltingTrapezoidArea = avgMeltingTrapezoidWidth * meltingTrapezoidHeight;
//
//				if (refreezingEnergy == 0.0) {
//					meltingEnergy += Math.abs(meltingTrapezoidArea);
//					refreezingEnergy += Math.abs(freezingTrapezoidArea);
//				} else {
//					refreezingEnergy += Math.abs(freezingTrapezoidArea);
//					remeltingEnergy += Math.abs(meltingTrapezoidArea);
//				}
//			}
//		}
//
////		System.out.printf("Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, Remelting Energy: %4d J/kg\n", (int) meltingEnergy, (int) refreezingEnergy, (int) remeltingEnergy);
//
//		double probSnowI = 15.4 * Math.exp(-0.29 * meltingEnergy); // fraction, NOT PERCENT
//		if (probSnowI > 1)
//			probSnowI = 1;
//		if (probSnowI < 0)
//			probSnowI = 0;
//		double probSnow = probSnowI * probIce; // fraction, NOT PERCENT
//
//		double probIcep = 0.0;
//		if (meltingEnergy > 0) {
//			double probIcepI = (2.3 * refreezingEnergy - 42 * Math.log(meltingEnergy + remeltingEnergy + 1) + 3)
//					/ 100.0; // fraction, NOT PERCENT
//			if (probIcepI > 1)
//				probIcepI = 1;
//			if (probIcepI < 0)
//				probIcepI = 0;
//			probIcep = probIcepI * probIce; // fraction, NOT PERCENT
//		}
//
//		double probRainI = (-2.1 * refreezingEnergy + 0.2 * (meltingEnergy + remeltingEnergy) + 458) / 100.0;
//		if (probRainI > 1)
//			probRainI = 1;
//		if (probRainI < 0)
//			probRainI = 0;
//		if (meltingEnergy < 5) {
//			probRainI = 0.2 * probRainI; // fraction, NOT PERCENT
//		}
//
//		double probRain = (1 - probIce) + probIce * probRainI; // fraction, NOT PERCENT
//		if (probRain > 1)
//			probRain = 1;
//
//		int dominantType = 0; // 0 = snow, 1 = rain, 2 = ice pellets
//
//		if (probSnow > probRain && probSnow > probIcep) {
//			dominantType = 0;
//		} else if (probRain > probSnow && probRain > probIcep) {
//			dominantType = 1;
//		} else if (probIcep > probRain && probIcep > probSnow) {
//			dominantType = 2;
//		}
//
//		if (dominantType == 0) {
//			if (probRain > probIcep && probRain > 0.5) {
//				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
//					return RAIN_SNOW_MIX;
//				} else {
//					return FRZR_SNOW_MIX;
//				}
//			} else if (probIcep > probRain && probIcep > 0.5) {
//				return ICEP_SNOW_MIX;
//			} else {
//				return SNOW;
//			}
//		} else if (dominantType == 1) {
//			if (probSnow > probIcep && probSnow > 0.5) {
//				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
//					return RAIN_SNOW_MIX;
//				} else {
//					return FRZR_SNOW_MIX;
//				}
//			} else if (probIcep > probSnow && probIcep > 0.5) {
//				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
//					return RAIN_ICEP_MIX;
//				} else {
//					return FRZR_ICEP_MIX;
//				}
//			} else {
//				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
//					return RAIN;
//				} else {
//					return FREEZING_RAIN;
//				}
//			}
//		} else if (dominantType == 2) {
//			if (probSnow > probRain && probSnow > 0.5) {
//				return ICEP_SNOW_MIX;
//			} else if (probRain > probSnow && probRain > 0.5) {
//				if (tmpIsobaric[tmpIsobaric.length - 1] > 273.15) {
//					return RAIN_ICEP_MIX;
//				} else {
//					return FRZR_ICEP_MIX;
//				}
//			} else {
//				return ICE_PELLETS;
//			}
//		} else {
//			return RAIN;
//		}
	}
	
	public static PrecipitationType bourgouinRevisedMethod(double[] pressureLevels, double[] tmpIsobaric,
			double[] dptIsobaric, double[] hgtIsobaric, double presSurface, double hgtSurface,
			boolean dynamicInitLayer) {
		return bourgouinRevisedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric, presSurface, hgtSurface, dynamicInitLayer, false);
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
			boolean dynamicInitLayer, boolean debug) {
		if(debug) System.out.println("ptype debug on");
		
//		System.out.println(tmpIsobaric[tmpIsobaric.length - 1]);
		if (tmpIsobaric[tmpIsobaric.length - 1] >= 283.15)
			return RAIN;

		List<RecordAtLevel> isobaricData = new ArrayList<>();
		for (int i = 0; i < dptIsobaric.length; i++) {
//			System.out.println("pres: " + pressureLevels[i]);
//			System.out.println("tmp: " + tmpIsobaric[i]);
//			System.out.println("dpt: " + dptIsobaric[i]);
//			System.out.println("wbt: " + WeatherUtils.wetBulbTemperature(tmpIsobaric[i], dptIsobaric[i], pressureLevels[i]));
			isobaricData.add(new RecordAtLevel(pressureLevels[i], tmpIsobaric[i], dptIsobaric[i], hgtIsobaric[i]));
//			System.out.printf("pressure: %4d mb temperature: %6.2f K, wetbulb: %6.2f K, dewpoint: %6.2f K, hgtIsobaric: %6d m\n", (int) pressureLevels[i]/100, tmpIsobaric[i], isobaricData.get(i).wetbulb, dptIsobaric[i], (int) hgtIsobaric[i]);
		}
		
		if(debug) {
			for (int i = 0; i < isobaricData.size(); i++) {
				System.out.printf("%7.1f\t%7.1f\t%7.2f\t%7.2f\t%7.2f\n", isobaricData.get(i).height, 
						isobaricData.get(i).pressure, isobaricData.get(i).temperature, 
						isobaricData.get(i).wetbulb, isobaricData.get(i).dewpoint);
			}
			System.out.println();
		}

//		System.out.println("built isobaricData");

		Collections.sort(isobaricData); // should sort top of atmosphere to start of list, surface to bottom of list

		// effectively removes subsurface records
		for (int i = 0; i < isobaricData.size(); i++) {
			double heightAtLevel = isobaricData.get(i).height;

			if (heightAtLevel < hgtSurface) {
				isobaricData.get(i).height = hgtSurface;
				isobaricData.get(i).pressure = presSurface;
				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
			}
		}

		// patch that removes 0 records in middle of profile
		for (int i = 0; i < isobaricData.size() - 1; i++) {
			double heightAtLevel = isobaricData.get(i).height - isobaricData.get(isobaricData.size() - 1).height;
			double heightAtLevelNext = isobaricData.get(i + 1).height - isobaricData.get(isobaricData.size() - 1).height;

			if (heightAtLevel < 0.1 && heightAtLevelNext > 1.9) {
				isobaricData.remove(i);
				i--; // ensures no entries are skipped after record removal
			}
		}
		
		if(debug) {
			for (int i = 0; i < isobaricData.size(); i++) {
				System.out.printf("%7.1f\t%7.1f\t%7.1f\t%7.1f\t%7.1f\n", isobaricData.get(i).height, 
						isobaricData.get(i).pressure, isobaricData.get(i).temperature, 
						isobaricData.get(i).wetbulb, isobaricData.get(i).dewpoint);
			}
		}

//		System.out.println(Arrays.toString(isobaricData.toArray()));

		double initPressureLayer = 35000;

		if (dynamicInitLayer) {
			// looks for first level with RH >= 90
			for (int i = 0; i < isobaricData.size(); i++) {
				double rh = WeatherUtils.relativeHumidity(isobaricData.get(i).temperature,
						isobaricData.get(i).dewpoint);
				
				if(debug) System.out.println("checking for init layer @ <" + isobaricData.get(i).pressure + ">: " + rh);

				if (rh >= SATURATION_THRESHOLD_RH) {
					initPressureLayer = isobaricData.get(i).pressure;
					break;
				}
			}
		}
		
//		initPressureLayer = 35000; // for testing purposes

//		System.out.println(initPressureLayer);

		int precipInitIndex = -1;
		for (int i = 0; i < isobaricData.size(); i++) {
			if (isobaricData.get(i).pressure >= initPressureLayer) {
				precipInitIndex = i;
				break;
			}
		}
//		System.out.println(precipInitIndex);

		// uses wet-bulb freezing point method if dynamic init layer is used and no
		// precip layers are detected
		if (precipInitIndex == -1 || isobaricData.get(precipInitIndex).height - 3 < hgtSurface)
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
				return FREEZING_RAIN;
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

		if(debug) System.out.printf("Melting Energy: %4d J/kg, Refreezing Energy: %4d J/kg, Remelting Energy: %4d J/kg\n", (int) meltingEnergy, (int) refreezingEnergy, (int) remeltingEnergy);
		
		double probSnowI = 15.4 * Math.exp(-0.29 * (meltingEnergy + remeltingEnergy)); // fraction, NOT PERCENT
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
			probRainI = 0.2 * probRainI * (meltingEnergy + remeltingEnergy); // fraction, NOT PERCENT
		}

		double probRain = (1 - probIce) + probIce * probRainI; // fraction, NOT PERCENT
		if (probRain > 1)
			probRain = 1;
		
		if(debug) {
			System.out.println("prob rain: " + probRain);
			System.out.println("prob icep: " + probIcep);
			System.out.println("prob snow: " + probSnow);
			System.out.println("prob ice: " + probIce);
		}

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
	 * @param presSurface      Surface pressure, in Pascals
	 * @param hgtSurface       Surface height, in Meters
	 * @param tmpSurface       Surface temperature, in Kelvins
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

	public static PrecipitationType bourgouinRevisedExtendedMethod(float[] pressureLevels, float[] tmpIsobaric,
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface, float tmpSurface,
			boolean dynamicInitLayer) {
		return bourgouinRevisedExtendedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric, presSurface, hgtSurface,
				tmpSurface, dynamicInitLayer, false);
	}

	/**
	 * @param pressureLevels   All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric      All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric      All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric      All isobaric heights of the profile, in Meters
	 * @param dynamicInitLayer true = find highest level with a dewpoint depression
	 *                         of 3 K, false = always init at 500 mb
	 * @return The precipitation type diagnosed by the Bourgouin Revised 2021 Method
	 *         (https://doi.org/10.1175/WAF-D-20-0118.1), extended by the Kuchera
	 *         Method and a surface vs. elevated-only check on the freezing rain
	 */
	public static PrecipitationType bourgouinRevisedExtendedMethod(float[] pressureLevels, float[] tmpIsobaric,
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface, float tmpSurface,
			boolean dynamicInitLayer, boolean debug) {
		if(debug) System.out.println("ptype debug on");
		if(debug) System.out.println(tmpIsobaric[tmpIsobaric.length - 1] + "\t" + tmpSurface);
		PrecipitationType preExtension = bourgouinRevisedMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric,
				presSurface, hgtSurface, true, debug);

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

		if(debug) System.out.println("tmpSurface: " + tmpSurface);
		if (preExtension == PrecipitationType.SNOW) {
			for(int i = 0; i < tmpIsobaric.length - 1; i++) {
				if(pressureLevels[i] > presSurface) {
					tmpIsobaric[i] = tmpIsobaric[tmpIsobaric.length - 1];
				}
			}
			
			double kucheraRatio = WeatherUtils.kucheraRatio(tmpIsobaric, dptIsobaric.length);
			
			if(debug) WeatherUtils.kucheraRatio(tmpIsobaric, true);
			if(debug) System.out.println("kuchera ratio: " + kucheraRatio);

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
	 * Described in
	 * https://journals.ametsoc.org/view/journals/apme/61/9/JAMC-D-21-0202.1.xml?tafloatb_body=pdf
	 * at Section 2.b.1
	 * 
	 * @param pressure   All pressure levels of the profile, in Pascals
	 * @param height     All isobaric heights of the profile, in Meters
	 * @param tmp2m      2 m Temperature, in Kelvins
	 * @param tmpSurface Surface Temperature, in Kelvins
	 */
	public static PrecipitationType cantinBachandMethod(double[] pressure, double[] height, double tmp2m,
			double tmpSurface) {
		double surfacePressure = pressure[pressure.length - 1];

		double hgt1000mb = logInterp(pressure, height, 100000);
		double hgt850mb = logInterp(pressure, height, 85000);
		double hgt700mb = logInterp(pressure, height, 70000);

		if (surfacePressure < 100000) {
			hgt1000mb = WeatherUtils.heightAtPressure(surfacePressure, 100000, tmp2m);
		}

		if (surfacePressure < 85000) {
			hgt850mb = WeatherUtils.heightAtPressure(surfacePressure, 85000, tmp2m);
		}

		if (surfacePressure < 70000) {
			hgt700mb = WeatherUtils.heightAtPressure(surfacePressure, 70000, tmp2m);
		}

		double thickness1000_850 = hgt850mb - hgt1000mb;
		double thickness850_700 = hgt700mb - hgt850mb;

//		for(int i = 0; i < height.length; i++) {
//			System.out.printf("%6.1f\t%6.0f\n", pressure[i]/100.0, height[i]);
//		}

//		System.out.println(surfacePressure);
//		System.out.println(hgt1000mb);
//		System.out.println(hgt850mb);
//		System.out.println(hgt700mb);
//		System.out.println();
//		
//		System.out.println(thickness1000_850);
//		System.out.println(thickness850_700);

		if (thickness1000_850 >= 1310) {
			return PrecipitationType.RAIN;
		} else if (thickness1000_850 >= 1290) {
			if (thickness850_700 >= 1540) {
				if (tmpSurface <= 273.15) {
					return PrecipitationType.FREEZING_RAIN;
				} else {
					return PrecipitationType.RAIN;
				}
			} else {
				return PrecipitationType.SNOW; // confirm in paper
			}
		} else {
			if (thickness850_700 >= 1540) {
				return PrecipitationType.ICE_PELLETS;
			} else {
				return PrecipitationType.SNOW;
			}
		}
	}

	public static PrecipitationType freezingPointMethod(double tmp2m) {
		if (tmp2m > 273.15) {
			return RAIN;
		} else {
			return SNOW;
		}
	}
	
	public static PrecipitationType ramerMethod(float[] pressureLevels, float[] tmpIsobaric, float[] dptIsobaric,
			float[] hgtIsobaric, float presSurface, float hgtSurface) {
		double[] pressureLevelsD = convFloatToDouble(pressureLevels);
		double[] tmpIsobaricD = convFloatToDouble(tmpIsobaric);
		double[] dptIsobaricD = convFloatToDouble(dptIsobaric);
		double[] hgtIsobaricD = convFloatToDouble(hgtIsobaric);
		
		return ramerMethod(pressureLevelsD, tmpIsobaricD, dptIsobaricD, hgtIsobaricD, (double) presSurface, (double) hgtSurface);
	}

	/**
	 * @param pressureLevels All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric    All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric    All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric    All isobaric heights of the profile, in Meters
	 * @return The precipitation type diagnosed by the Ramer 1993 Method
	 *         (https://ams.confex.com/ams/pdfpapers/30176.pdf)
	 */
	public static PrecipitationType ramerMethod(double[] pressureLevels, double[] tmpIsobaric, double[] dptIsobaric,
			double[] hgtIsobaric, double presSurface, double hgtSurface) {
//		System.out.println("start ramer");

		if (tmpIsobaric[tmpIsobaric.length - 1] >= 283.15) {
//			System.out.println("sfc >10");
			return RAIN;
		}

		List<RecordAtLevel> isobaricData = new ArrayList<>();
		for (int i = 0; i < dptIsobaric.length; i++) {
			isobaricData.add(new RecordAtLevel(pressureLevels[i], tmpIsobaric[i], dptIsobaric[i], hgtIsobaric[i]));
		}

		Collections.sort(isobaricData);

		// effectively removes subsurface records
		for (int i = 0; i < dptIsobaric.length; i++) {
			double heightAtLevel = isobaricData.get(i).height;

			if (heightAtLevel < hgtSurface) {
				isobaricData.get(i).height = hgtSurface;
				isobaricData.get(i).pressure = presSurface;
				isobaricData.get(i).temperature = isobaricData.get(i - 1).temperature;
				isobaricData.get(i).wetbulb = isobaricData.get(i - 1).wetbulb;
				isobaricData.get(i).dewpoint = isobaricData.get(i - 1).dewpoint;
			}
		}

		// patch that removes 0 records in middle of profile
		for (int i = 0; i < isobaricData.size() - 1; i++) {
			double heightAtLevel = isobaricData.get(i).height - isobaricData.get(isobaricData.size() - 1).height;
			double heightAtLevelNext = isobaricData.get(i + 1).height - isobaricData.get(isobaricData.size() - 1).height;

			if (heightAtLevel < 0.1 && heightAtLevelNext > 1.9) {
				isobaricData.remove(i);
				i--; // ensures no entries are skipped after record removal
			}
		}

		double initPressureLayer = 50000;

		for (int i = 0; i < isobaricData.size(); i++) {
			double rh = WeatherUtils.relativeHumidity(isobaricData.get(i).temperature, isobaricData.get(i).dewpoint);

			if (rh >= SATURATION_THRESHOLD_RH) {
				initPressureLayer = isobaricData.get(i).pressure;
				break;
			}
		}

		if (isobaricData.get(isobaricData.size() - 1).wetbulb > 273.15 + 2) {
//			System.out.println("sfc >2");
			return RAIN;
		}

		int precipInitIndex = -1;
		for (int i = 0; i < dptIsobaric.length; i++) {
			if (isobaricData.get(i).pressure >= initPressureLayer) {
				precipInitIndex = i;
				break;
			}
		}

		boolean allFreezing = true;
		for (int i = precipInitIndex; i < isobaricData.size() - 1; i++) {
			if (isobaricData.get(i).wetbulb > 273.15 - 6.6) {
				allFreezing = false;
				break;
			}
		}
		if (allFreezing && isobaricData.get(isobaricData.size() - 1).wetbulb <= 273.15 + 2) {
//			System.out.println("all freezing");
			return SNOW;
		}

		double hydrometeorIceFraction = 0;
		if (isobaricData.get(precipInitIndex).getWetbulb() < 273.15 - 6.6) {
			hydrometeorIceFraction = 1;
		}

		double minimumIceFraction = hydrometeorIceFraction;

		for (int i = precipInitIndex + 1; i < isobaricData.size(); i++) {
			double wetbulb = isobaricData.get(i).getWetbulb();

			double temperature = isobaricData.get(i).getTemperature();
			double dewpoint = isobaricData.get(i).getDewpoint();
			double rh = WeatherUtils.relativeHumidity(temperature, dewpoint);

			double pressure1 = isobaricData.get(i - 1).getPressure();
			double pressure0 = isobaricData.get(i).getPressure();

			double dLnP = Math.log(pressure0) - Math.log(pressure1);
			double e = 0.45 * rh;

			double dI = (273.15 - wetbulb) / e * dLnP;

			hydrometeorIceFraction += dI;

			if (hydrometeorIceFraction < 0) {
				hydrometeorIceFraction = 0;
			}

			if (hydrometeorIceFraction > 1) {
				hydrometeorIceFraction = 1;
			}

			minimumIceFraction = Double.min(minimumIceFraction, hydrometeorIceFraction);

//			System.out.printf("%4.2f\t%4.2f\n", hydrometeorIceFraction, minimumIceFraction);
		}

		if (hydrometeorIceFraction > 0.85) {
			double refreezingFraction = hydrometeorIceFraction - minimumIceFraction;

			if (refreezingFraction > 0.85) {
				return ICE_PELLETS;
			} else if (refreezingFraction >= 0.15) {
				return ICEP_SNOW_MIX;
			} else {
				return SNOW;
			}
		} else if (hydrometeorIceFraction > 0.04) {
			double refreezingFraction = hydrometeorIceFraction - minimumIceFraction;
			double wbt2m = isobaricData.get(isobaricData.size() - 1).getWetbulb();

			if (refreezingFraction > minimumIceFraction) {
				if (wbt2m < 273.15) {
					return FRZR_ICEP_MIX;
				} else {
					return RAIN_ICEP_MIX;
				}
			} else {
				if (wbt2m < 273.15) {
					return FRZR_SNOW_MIX;
				} else {
					return RAIN_SNOW_MIX;
				}
			}
		} else {
			double wbt2m = isobaricData.get(isobaricData.size() - 1).getWetbulb();

			if (wbt2m < 273.15) {
				return FREEZING_RAIN;
			} else {
				return RAIN;
			}
		}
	}

	/**
	 * @param pressureLevels   All pressure levels of the profile, in Pascals
	 * @param tmpIsobaric      All isobaric temperatures of the profile, in Kelvins
	 * @param dptIsobaric      All isobaric dewpoints of the profile, in Kelvins
	 * @param hgtIsobaric      All isobaric heights of the profile, in Meters
	 * @param presSurface      Surface pressure, in Pascals
	 * @param hgtSurface       Surface height, in Meters
	 * @param tmpSurface       Surface temperature, in Kelvins
	 * @param dynamicInitLayer true = find highest level with a dewpoint depression
	 *                         of 3 K, false = always init at 500 mb
	 * @return The precipitation type diagnosed by the Bourgouin Revised 2021 Method
	 *         (https://doi.org/10.1175/WAF-D-20-0118.1), extended by the Kuchera
	 *         Method and a surface vs. elevated-only check on the freezing rain
	 */
	public static PrecipitationType ramerExtendedMethod(float[] pressureLevels, float[] tmpIsobaric,
			float[] dptIsobaric, float[] hgtIsobaric, float presSurface, float hgtSurface, float tmpSurface) {
		PrecipitationType preExtension = ramerMethod(pressureLevels, tmpIsobaric, dptIsobaric, hgtIsobaric,
				presSurface, hgtSurface);

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
			double kucheraRatio = WeatherUtils.kucheraRatio(tmpIsobaric, dptIsobaric.length);

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

	private static double[] convFloatToDouble(float[] arr) {
		double[] ret = new double[arr.length];

		for (int i = 0; i < ret.length; i++) {
			ret[i] = (double) arr[i];
		}

		return ret;
	}

	// inputArr assumed to already be sorted and increasing
	private static double logInterp(double[] inputArr, double[] outputArr, double input) {
		if (input < inputArr[0]) {
			return outputArr[0];
		} else if (input >= inputArr[inputArr.length - 1]) {
			return outputArr[outputArr.length - 1];
		} else {
			for (int i = 0; i < inputArr.length - 1; i++) {
				if (i + 1 == outputArr.length) {
					return outputArr[outputArr.length - 1];
				}

				double input1 = inputArr[i];
				double input2 = inputArr[i + 1];

				if (input == input1) {
					return outputArr[i];
				} else if (input < input2) {
					double logInput1 = Math.log(input1);
					double logInput2 = Math.log(input2);
					double logInput = Math.log(input);

					double output1 = outputArr[i];
					double output2 = outputArr[i + 1];

					double weight1 = (logInput2 - logInput) / (logInput2 - logInput1);
					double weight2 = (logInput - logInput1) / (logInput2 - logInput1);

					return output1 * weight1 + output2 * weight2;
				} else {
					continue;
				}
			}

			return -1024.0;
		}
	}

	private static double linScale(double preMin, double preMax, double postMin, double postMax, double value) {
		double slope = (postMax - postMin) / (preMax - preMin);

		return slope * (value - preMin) + postMin;
	}
}
