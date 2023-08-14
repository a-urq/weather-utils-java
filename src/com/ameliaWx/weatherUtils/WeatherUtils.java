package com.ameliaWx.weatherUtils;

import java.util.ArrayList;

/**
 * @author Amelia Urquhart
 * @apiNote All methods use SI units, except for angles which use degrees. A
 *          separate class, 'UnitConversions' has been provided (WIP).
 */
public class WeatherUtils {
	/**
	 * Section 1 --- Equivalent to MetBasics.cpp in C++ version of the library
	 */

	/** Units: J K^-1 mol^-1 */
	public static final double molarGasConstant = 8.314;
	/** Units: kg mol^-1 */
	public static final double avgMolarMass = 0.029;

	/** Units: J kg^-1 K^-1 */
	public static final double dryAirGasConstant = 287;
	/** Units: J kg^-1 K^-1 */
	public static final double waterVaporGasConstant = 461.5;

	/** Units: J kg^-1 */
	public static final double latentHeatOfVaporization = 2500000;
	/** Units: J kg^-1 */
	public static final double latentHeatOfMelting = 334000;
	/** Units: J kg^-1 */
	public static final double latentHeatOfSublimation = 2830000;

	/** Units: J kg^-1 K^-1 */
	public static final double specificHeatCapacityDryAir = 1007.5; // constant pressure
	/** Units: J kg^-1 K^-1 */
	public static final double specificHeatCapacityWaterVapor = 1864; // constant pressure

	/** Units: m s^-2 */
	public static final double gravAccel = 9.81;

	/**
	 * Computes the partial density of water vapor, also called absolute humidity.
	 * 
	 * @param vaporPressure Units: Pascals
	 * @param temperature   Units: Kelvins
	 * @return <b>absoluteHumidity</b> Units: kg m^-3
	 */
	public static double absoluteHumidity(double vaporPressure, double temperature) {
		double waterVaporDensity = vaporPressure / (waterVaporGasConstant * temperature);

		return waterVaporDensity;
	}

	/**
	 * Computes the partial density of dry air.
	 * 
	 * @param dryAirPressure Units: Pascals
	 * @param temperature    Units: Kelvins
	 * @return <b>dryAirDensity</b> Units: kg m^-3
	 */
	public static double dryAirDensity(double dryAirPressure, double temperature) {
		double dryAirDensity = dryAirPressure / (dryAirGasConstant * temperature);

		return dryAirDensity;
	}

	/**
	 * Computes the dewpoint using temperature and relative humidity.
	 * 
	 * @param temperature      Units: Kelvins
	 * @param relativeHumidity Units: Fraction (not Percent!)
	 * @return <b>dewpoint</b> Units: Kelvins
	 */
	public static double dewpoint(double temperature, double relativeHumidity) {
		double e0 = 611; // Pascals
		double t0 = 273.15; // Kelvins

		double es = vaporPressure(temperature);

		double dewpointReciprocal = 1 / t0
				- waterVaporGasConstant / latentHeatOfVaporization * Math.log(es * relativeHumidity / e0);

		return 1 / dewpointReciprocal;
	}

	/**
	 * Computes the equivalent potential temperature. Formula comes from David
	 * Bolton (1980): "The Computation of Equivalent Potential Temperature" <a href=
	 * 'https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2'></a>
	 * 
	 * @param temperature Units: Kelvins
	 * @param dewpoint    Units: Kelvins
	 * @param pressure    Units: Pascals
	 * @return <b>equivalentPotentialTemperature</b> Units: Kelvins
	 */
	public static double equivalentPotentialTemperature(double temperature, double dewpoint, double pressure) {
		double mixingRatio = mixingRatio(pressure, dewpoint); // units of g g^-1, NOT g kg^-1
		double relativeHumidity = relativeHumidity(temperature, dewpoint);

		double T_L = 1 / ((1 / (temperature - 55)) - Math.log(relativeHumidity) / 2840) + 55;

		double thetaE = temperature * Math.pow(100000 / pressure, 0.2854 * (1 - 0.28 * mixingRatio))
				* Math.exp((3.376 / T_L - 0.00254) * (1000.0 * mixingRatio) * (1 + 0.81 * mixingRatio));

		return thetaE;
	}

	/**
	 * Computes height above sea level at a given pressure. Note that this assumes a
	 * constant scale height of 8500 m, regardless of the air temperatures.
	 * 
	 * @param seaLevelPres  Units: Pascals
	 * @param pressureAtAlt Units: Pascals
	 * @return <b>heightAtPressure</b> Units: Meter
	 */
	public static double heightAtPressure(double seaLevelPres, double pressureAtAlt) {
		double scaleHeight = 8500; // Meters

		return scaleHeight * Math.log(seaLevelPres / pressureAtAlt);
	}

	/**
	 * Computes height above sea level at a given pressure. Takes a temperature to estimate scale height.
	 * 
	 * @param seaLevelPres  Units: Pascals
	 * @param pressureAtAlt Units: Pascals
	 * @param temperature   Units: Kelvins
	 * @return <b>heightAtPressure</b> Units: Meter
	 */
	public static double heightAtPressure(double seaLevelPres, double pressureAtAlt, double temperature) {
		double scaleHeight = (molarGasConstant * temperature) / (avgMolarMass * gravAccel); // Meters

		return scaleHeight * Math.log(seaLevelPres / pressureAtAlt);
	}

	/**
	 * Computes the mixing ratio using total pressure, vapor pressure, and
	 * temperature.
	 * 
	 * @param pressure Units: Pascals
	 * @param dewpoint Units: Kelvins
	 * @return <b>mixingRatio</b> Units: Fraction
	 */
	public static double mixingRatio(double pressure, double dewpoint) {
		double vaporPressure = WeatherUtils.vaporPressure(dewpoint);
		double mixingRatio = 0.62197 * (vaporPressure) / (pressure - vaporPressure);

		return mixingRatio;
	}
	
	public static double moistStaticEnergy(double temperature, double dewpoint, double height, double pressure) {
		double vaporPressure = WeatherUtils.vaporPressure(dewpoint);
		double specificHumidity = WeatherUtils.specificHumidity(pressure, vaporPressure, temperature);
		
		double moistStaticEnergy = specificHeatCapacityDryAir * temperature + latentHeatOfVaporization * specificHumidity + gravAccel * height;
		
		return moistStaticEnergy;
	}

	private static final double REFERENCE_PRESSURE = 100000; // Pascals

	public static double potentialTemperature(double temperature, double pressure) {
		double potentialTemperature = temperature * Math.pow(REFERENCE_PRESSURE / pressure, 0.286);

		return potentialTemperature;
	}

	/**
	 * Computes pressure at a given height above sea level. Note that this assumes a
	 * constant scale height of 7290 m, regardless of the air temperatures.
	 * 
	 * @param seaLevelPres        Units: Pascals
	 * @param heightAboveSeaLevel Units: Meters
	 * @return <b>pressureAtHeight</b> Units: Pascals
	 */
	public static double pressureAtHeight(double seaLevelPres, double heightAboveSeaLevel) {
		double scaleHeight = 8500; // Meters

		return seaLevelPres * Math.exp(-heightAboveSeaLevel / scaleHeight);
	}

	/**
	 * Computes pressure at a given height above sea level. Uses temperature to
	 * figure out scale height.
	 * 
	 * @param seaLevelPres        Units: Pascals
	 * @param heightAboveSeaLevel Units: Meters
	 * @param temperature         Units: Kelvins
	 * @return <b>pressureAtHeight</b> Units: Pascals
	 */
	public static double pressureAtHeight(double seaLevelPres, double heightAboveSeaLevel, double temperature) {
		double scaleHeight = (molarGasConstant * temperature) / (avgMolarMass * gravAccel); // Meters

		return seaLevelPres * Math.exp(-heightAboveSeaLevel / scaleHeight);
	}

	/**
	 * Computes precipitable water given a vertical profile of height, temperature,
	 * and dewpoint.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure/decreasing
	 * height.
	 * 
	 * @param height      Units: Meters
	 * @param temperature Units: Kelvins
	 * @param dewpoint    Units: Kelvins
	 * @return <b>precipitableWater</b> Units: kg m^-2 (equivalent to mm of water)
	 */
	public static double precipitableWater(double[] height, double[] temperature, double[] dewpoint) {
		double pwat = 0.0; // kg m^-2

		for (int i = 0; i < dewpoint.length - 1; i++) {
			double height1 = height[i];
			double temperature1 = temperature[i];
			double dewpoint1 = dewpoint[i];

			double vaporPressure1 = WeatherUtils.vaporPressure(dewpoint1);
			double absoluteHumidity1 = WeatherUtils.absoluteHumidity(vaporPressure1, temperature1);

			double height2 = height[i + 1];
			double temperature2 = temperature[i + 1];
			double dewpoint2 = dewpoint[i + 1];

			double vaporPressure2 = WeatherUtils.vaporPressure(dewpoint2);
			double absoluteHumidity2 = WeatherUtils.absoluteHumidity(vaporPressure2, temperature2);

			double absoluteHumidity = (absoluteHumidity1 + absoluteHumidity2) / 2.0; // kg m^-3
			double dz = height1 - height2; // m

			pwat += absoluteHumidity * dz;
		}

		return pwat;
	}

	/**
	 * Computes relative humidity using temperature and dewpoint.
	 * 
	 * @param temperature Units: Kelvins
	 * @param dewpoint    Units: Kelvins
	 * @return <b>relativeHumidity</b> Units: Fraction (not Percent!)
	 */
	public static double relativeHumidity(double temperature, double dewpoint) {
		double vaporPres = vaporPressure(dewpoint);
		double satVaporPres = vaporPressure(temperature);

		return vaporPres / satVaporPres;
	}

	/**
	 * Computes the specific heat capacity (constant pressure) of a parcel of air
	 * with the given pressure, temperature, and dewpoint.
	 * 
	 * @param pressure    Units: Pascals
	 * @param temperature Units: Kelvins
	 * @param dewpoint    Units: Kelvins
	 * @return <b>specificHeatCapacity</b> Units: J kg^-1 K^-1
	 */
	public static double specificHeatCapacityAtDewpoint(double pressure, double temperature, double dewpoint) {
		double vaporPressure = vaporPressure(dewpoint);
		double specificHumidity = specificHumidity(pressure, vaporPressure, temperature);

		return specificHeatCapacityAtSpecificHumidity(specificHumidity);
	}

	/**
	 * Computes the specific heat capacity (constant pressure) of a parcel of air
	 * with the given specific humidity.
	 * 
	 * @param specificHumidity Units: Fraction
	 * @return <b>specificHeatCapacity</b> Units: J kg^-1 K^-1
	 */
	public static double specificHeatCapacityAtSpecificHumidity(double specificHumidity) {
		double specificHeatCapacity = specificHumidity * specificHeatCapacityWaterVapor
				+ (1 - specificHumidity) * specificHeatCapacityDryAir;

		return specificHeatCapacity;
	}

	/**
	 * Computes the specific humidity using total pressure, vapor pressure, and
	 * temperature.
	 * 
	 * @param pressure      Units: Pascals
	 * @param vaporPressure Units: Pascals
	 * @param temperature   Units: Kelvins
	 * @return <b>specificHumidity</b> Units: Fraction
	 */
	public static double specificHumidity(double pressure, double vaporPressure, double temperature) {
		double waterVaporDensity = absoluteHumidity(vaporPressure, temperature); // kg m^-3
		double airDensity = dryAirDensity(pressure - vaporPressure, temperature); // kg m^-3

		return waterVaporDensity / (waterVaporDensity + airDensity);
	}

	/**
	 * Computes the equivalent potential temperature. Formula comes from David
	 * Bolton (1980): "The Computation of Equivalent Potential Temperature" <a href=
	 * 'https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2'></a>
	 * 
	 * @param temperature Units: Kelvins
	 * @param dewpoint    Units: Kelvins
	 * @param pressure    Units: Pascals
	 * @return <b>thetaE</b> Units: Kelvins
	 */
	public static double thetaE(double temperature, double dewpoint, double pressure) {
		return equivalentPotentialTemperature(temperature, dewpoint, pressure);
	}

	/**
	 * Computes the partial pressure of water vapor in air with the dewpoint given.
	 * 
	 * @param dewpoint Units: Kelvins
	 * @return <b>vaporPressure</b> Units: Pascals
	 */
	public static double vaporPressure(double dewpoint) {
		double e0 = 611; // Pascals
		double t0 = 273.15; // Kelvins

		return e0 * Math.exp(latentHeatOfVaporization / waterVaporGasConstant * (1 / t0 - 1 / dewpoint));
	}

	// make virt temp
	public static double virtualTemperature(double temperature, double dewpoint, double pressure) {
		double vaporPressure = vaporPressure(dewpoint);

		double eOverP = vaporPressure / pressure;
		double virtualTemperature = temperature / (1 - eOverP * (1 - dryAirGasConstant / waterVaporGasConstant));

		return virtualTemperature;
	}

	/**
	 * Computes the wet-bulb temperature from temperature, dewpoint, and station
	 * pressure. Known to have gone into infinite loops in the past. If this
	 * happens, please tell me and also tell me the method inputs used.
	 * 
	 * @param tmp  Units: Kelvins
	 * @param dpt  Units: Kelvins
	 * @param pres Units: Pascals
	 * @return <b>wetBulbTemperature</b> Units: Kelvins
	 */
	static final double WBT_TOLERANCE = 0.004; // this value is arbitrary; it's worked well for me
	static final int ITER_LIMIT = 100; // sets an upper bound on amount of iterations, prevents infinite loops

	public static double wetBulbTemperature(double tmp, double dpt, double pres) {
		// tbh i have no clue where this implementation came from or how it works, but
		// it's been faithfully serving me for years.

		double tWB = tmp - 273.15;

		double vaporPres = vaporPressure(dpt);
		double mixingRatio1 = 0.62197 * (vaporPres / (pres - vaporPres));

		for (int i = 0; i < ITER_LIMIT; i++) {
			double satVapPres = vaporPressure(tWB + 273.15);
			double wSWB = 0.62197 * (satVapPres / (pres - satVapPres));
			double mixingRatio2 = ((2501 - (4.186 - 1.805) * tWB) * wSWB - 1.006 * ((tmp - 273.15) - tWB))
					/ (2501 + 1.805 * (tmp - 273.15) - 4.186 * tWB);

			if (Math.abs(Math.log(mixingRatio1) - Math.log(mixingRatio2)) < WBT_TOLERANCE) { // logarithm used to ensure
																								// equal accuracy at low
																								// and high tmps and
																								// dpts
				break;
			}

//			System.out.println(mixingRatio1);
//			System.out.println(mixingRatio2);
//			System.out.println();

			if (mixingRatio1 > mixingRatio2) {
				tWB += 0.01;
//				System.out.println("inc wbt to " + (tWB + 273.15));
			} else {
//				System.out.println(tWB + 273.15);
				tWB -= 500 * Math.abs(mixingRatio1 - mixingRatio2);
//				System.out.println("dec wbt to " + (tWB + 273.15) + " (" + Math.abs(mixingRatio1 - mixingRatio2) + ")");
			}
		}

		return tWB + 273.15;
	}

	/**
	 * Section 2 --- Equivalent to PhysMet.cpp in C++ version of the library
	 */

	/**
	 * Section 3 --- Equivalent to Kinematics.cpp in C++ version of the library
	 */

	/**
	 * Units: s^-1
	 */
	public static final double coriolisConstant = 0.00007292;

	/**
	 * Computes the coriolis parameter.
	 * 
	 * @param latitude Units: Degrees
	 * @return <b>coriolisParameter</b> Units: s^-1
	 */
	public static double coriolisParameter(double latitude) {
		double f = 2 * coriolisConstant * Math.sin(Math.toRadians(latitude));

		return f;
	}

	/**
	 * Computes the coriolis force on air with a given velocity. Expected order of
	 * components: zonal (+east), meridional (+north), altitudinal (+up).
	 * 
	 * @param velocity Units: m s^-1
	 * @param latitude Units: Degrees
	 * @return <b>coriolisForce</b> Units: m s^-2
	 */
	public static double[] coriolisForce(double[] velocity, double latitude) {
		int dimensions = velocity.length;

		if (dimensions <= 1 || dimensions >= 4) {
			return new double[dimensions]; // returns array of zeroes
		}

		double u = velocity[0];
		double v = velocity[1];
		double f = coriolisParameter(latitude);

		if (dimensions == 2) {
			double[] force = { f * v, -f * u };

			return force;
		} else { // will only execute if dimensions == 3
			double[] force = { f * v, -f * u, 0 };

			return force;
		}
	}

	/**
	 * Section 4 - Convective Meteorology Equivalent to ConvMet.cpp in C++ version
	 * of the library
	 */

	/** Units: K m^-1 */
	public static final double dryAdiabaticLapseRate = 0.0098;
	/** Units: K m^-1 */
	public static final double dewpointLapseRate = 0.0018;

	/**
	 * Computes the bulk wind shear
	 * 
	 * @param pressure   Units: array of Pascals
	 * @param uWind      Units: array of m s^-1
	 * @param vWind      Units: array of m s^-1
	 * @param lowerLimit Units: Meters
	 * @param upperLimit Units: Meters
	 * @return <b>bulkShearMagnitude</b> Units: m s^-1
	 */
	public static double bulkShearMagnitude(double[] height, double[] uWind, double[] vWind, double lowerLimit,
			double upperLimit) {
		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < vWind.length; i++) {
			heightRev[vWind.length - 1 - i] = height[i];
			uWindRev[vWind.length - 1 - i] = uWind[i];
			vWindRev[vWind.length - 1 - i] = vWind[i];
		}

		double uWindLower = linearInterp(heightRev, uWindRev, lowerLimit);
		double vWindLower = linearInterp(heightRev, vWindRev, lowerLimit);
		double uWindUpper = linearInterp(heightRev, uWindRev, upperLimit);
		double vWindUpper = linearInterp(heightRev, vWindRev, upperLimit);

//		System.out.println(Arrays.toString(heightRev));
//		System.out.println(Arrays.toString(uWindRev));
//		System.out.printf("%4.1f\t%4.1f\n", lowerLimit, upperLimit);
//		System.out.printf("%4.1f\t%4.1f\t%4.1f\t%4.1f\n\n", uWindLower, vWindLower, uWindUpper, vWindUpper);

		double bulkShearMagnitude = Math.hypot(uWindUpper - uWindLower, vWindUpper - vWindLower);

		return bulkShearMagnitude;
	}

	/**
	 * Computes the convective available potential energy given a parcel path.
	 * Assumes that all arrays are sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @param parcelPath     Type: ArrayList of RecordAtLevel, should be type
	 *                       returned by computeParcelPath()
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeCape(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			ArrayList<RecordAtLevel> parcelPath) {
		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		return computeCape(envPressure, envTemperature, envDewpoint, parcelPressure, parcelHeight, parcelTemperature,
				parcelDewpoint);
	}

	/**
	 * Computes the convective available potential energy given a parcel path.
	 * Assumes that all arrays are sorted in order of increasing pressure.
	 * 
	 * @param envPressure       Units: array of Pascals
	 * @param envTemperature    Units: array of Kelvins
	 * @param envDewpoint       Units: array of Kelvins
	 * @param parcelPressure    Units: array of Pascals
	 * @param parcelTemperature Units: array of Kelvins
	 * @param parcelDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeCape(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			double[] parcelPressure, double[] parcelHeight, double[] parcelTemperature, double[] parcelDewpoint) {
		double potentialEnergy = 0.0; // J kg^-1

		double equilibriumLevel = equilibriumLevel(envPressure, envTemperature, envDewpoint, parcelPressure,
				parcelHeight, parcelTemperature, parcelDewpoint);
		double levelOfFreeConvection = levelOfFreeConvection(envPressure, envTemperature, envDewpoint, parcelPressure,
				parcelHeight, parcelTemperature, parcelDewpoint);

		for (int i = 0; i < parcelPressure.length - 1; i++) {
			double pPres = parcelPressure[i];
			double height = parcelHeight[i];

//			System.out.printf("%6.1f\t", height);
//			System.out.printf("%6.1f\t", equilibriumLevel);

			if (height < equilibriumLevel && height > levelOfFreeConvection && pPres >= envPressure[0]) {
				double pTemp = parcelTemperature[i];
				double pDwpt = parcelDewpoint[i];

				double ePres = parcelPressure[i];
				double eTemp = logInterp(envPressure, envTemperature, ePres);
				double eDwpt = logInterp(envPressure, envDewpoint, ePres);

				double eVirtTemp = virtualTemperature(eTemp, eDwpt, ePres);
				double pVirtTemp = virtualTemperature(pTemp, pDwpt, pPres);

				double energyAdded = gravAccel * (pVirtTemp - eVirtTemp) / eVirtTemp;
//				System.out.println(dz);
//				System.out.printf("%6.1f\t", ePres/100.0);
//				System.out.print(energyAdded + " energyAdded\t");
//				System.out.print((energyAdded > 0) + " energyAdded\t");

				if (energyAdded > 0 && pDwpt >= pTemp) {
					double dz = parcelHeight[i] - parcelHeight[i + 1];

					potentialEnergy += energyAdded * dz;

//					System.out.printf("%7.1f\t%5.1f\t%5.1f\n", ePres/100.0, eVirtTemp, pVirtTemp);
//					System.out.println(energyAdded * dz + " J kg^-1 added");
//					System.out.println(potentialEnergy + " J kg^-1 total");
				} else {
					return potentialEnergy;
				}
			}
		}

		return potentialEnergy;
	}

	/**
	 * Computes SBCAPE given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeSbcape(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.SURFACE_BASED, false);

		return computeCape(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes ML50CAPE given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMl50cape(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.MIXED_LAYER_50MB, false);

		return computeCape(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes ML100CAPE given an environmental sounding Assumes that all arrays
	 * are sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMl100cape(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.MIXED_LAYER_100MB, false);

		return computeCape(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes MLCAPE given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure     Units: array of Pascals
	 * @param envTemperature  Units: array of Kelvins
	 * @param envDewpoint     Units: array of Kelvins
	 * @param mixedLayerDepth Units: Pascals
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMlcape(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			double mixedLayerDepth) {
		ArrayList<RecordAtLevel> parcelPath = computeMixedLayerParcelPath(envPressure, envTemperature, envDewpoint,
				mixedLayerDepth, false);

		return computeCape(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes MUCAPE given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMucape(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.MOST_UNSTABLE, false);

		return computeCape(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes the convective inhibition given a parcel path. Assumes that all
	 * arrays are sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @param parcelPath     Type: ArrayList of RecordAtLevel, should be type
	 *                       returned by computeParcelPath()
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeCinh(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			ArrayList<RecordAtLevel> parcelPath) {
		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		return computeCinh(envPressure, envTemperature, envDewpoint, parcelPressure, parcelHeight, parcelTemperature,
				parcelDewpoint);
	}

	/**
	 * Computes the convective available potential energy given a parcel path.
	 * Assumes that all arrays are sorted in order of increasing pressure.
	 * 
	 * @param envPressure       Units: array of Pascals
	 * @param envTemperature    Units: array of Kelvins
	 * @param envDewpoint       Units: array of Kelvins
	 * @param parcelPressure    Units: array of Pascals
	 * @param parcelTemperature Units: array of Kelvins
	 * @param parcelDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeCinh(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			double[] parcelPressure, double[] parcelHeight, double[] parcelTemperature, double[] parcelDewpoint) {
		double potentialEnergy = 0.0; // J kg^-1

		double levelOfFreeConvection = levelOfFreeConvection(envPressure, envTemperature, envDewpoint, parcelPressure,
				parcelHeight, parcelTemperature, parcelDewpoint);

		for (int i = 0; i < parcelPressure.length - 1; i++) {
			double pPres = parcelPressure[i];
			double height = parcelHeight[i];

//			System.out.printf("%6.1f\t", height);
//			System.out.printf("%6.1f\t", equilibriumLevel);

			if (height < levelOfFreeConvection && pPres >= envPressure[0]) {
				double pTemp = parcelTemperature[i];
				double pDwpt = parcelDewpoint[i];

				double ePres = parcelPressure[i];
				double eTemp = logInterp(envPressure, envTemperature, ePres);
				double eDwpt = logInterp(envPressure, envDewpoint, ePres);

				double eVirtTemp = virtualTemperature(eTemp, eDwpt, ePres);
				double pVirtTemp = virtualTemperature(pTemp, pDwpt, pPres);

				double energyAdded = gravAccel * (pVirtTemp - eVirtTemp) / eVirtTemp;
//				System.out.println(dz);
//				System.out.printf("%6.1f\t", ePres/100.0);
//				System.out.print(energyAdded + " energyAdded\t");
//				System.out.print((energyAdded > 0) + " energyAdded\t");

				if (energyAdded < 0) {
					double dz = parcelHeight[i] - parcelHeight[i + 1];

					potentialEnergy += energyAdded * dz;

//					System.out.printf("%7.1f\t%5.1f\t%5.1f\n", ePres/100.0, eVirtTemp, pVirtTemp);
//					System.out.println(energyAdded * dz + " J kg^-1 added");
//					System.out.println(potentialEnergy + " J kg^-1 total");
				} else {
					return potentialEnergy;
				}
			}
		}

		return potentialEnergy;
	}

	/**
	 * Computes SBCINH given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeSbcinh(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.SURFACE_BASED, false);

		return computeCinh(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes ML50CINH given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMl50cinh(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.MIXED_LAYER_50MB, false);

		return computeCinh(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes ML100CINH given an environmental sounding Assumes that all arrays
	 * are sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMl100cinh(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.MIXED_LAYER_100MB, false);

		return computeCinh(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes MLCINH given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure     Units: array of Pascals
	 * @param envTemperature  Units: array of Kelvins
	 * @param envDewpoint     Units: array of Kelvins
	 * @param mixedLayerDepth Units: Pascals
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMlcinh(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			double mixedLayerDepth) {
		ArrayList<RecordAtLevel> parcelPath = computeMixedLayerParcelPath(envPressure, envTemperature, envDewpoint,
				mixedLayerDepth, false);

		return computeCinh(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes MUCINH given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeMucinh(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.MOST_UNSTABLE, false);

		return computeCinh(envPressure, envTemperature, envDewpoint, parcelPath);
	}

	/**
	 * Computes MUCINH given an environmental sounding Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeDcape(double[] envPressure, double[] envTemperature, double[] envDewpoint) {
		ArrayList<RecordAtLevel> parcelPath = computeDcapeParcelPath(envPressure, envTemperature, envDewpoint);

		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		double potentialEnergy = 0.0; // J kg^-1

		for (int i = 0; i < parcelPressure.length - 1; i++) {
			double pTemp = parcelTemperature[i];

			double ePres = parcelPressure[i];
			double eTemp = logInterp(envPressure, envTemperature, ePres);

			double energyAdded = gravAccel * (eTemp - pTemp) / eTemp;

			if (energyAdded > 0) {
				double dz = parcelHeight[i + 1] - parcelHeight[i];

				potentialEnergy += energyAdded * dz;
			} else {
				return potentialEnergy;
			}
		}

		return potentialEnergy;
	}

	/**
	 * Computes the convective parcel path.
	 * 
	 * If entrainment is turned on, assumes a constant entrainment rate of 0.00001 m^-1.
	 * Not accurate to the study, but a fairly good approximation.
	 * 
	 * Also entrains air based on temperature and mixing ratio rather than moist static energy as in the study.
	 * Unsure if equivalent to the study.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @return <b>parcelPath</b> RecordAtLevel {pressure: Pascals, temperature:
	 *         Kelvins, wetbulb: Kelvins, dewpoint: Kelvins, height: Meters}
	 */
	public static ArrayList<RecordAtLevel> computeParcelPath(double[] pressure, double[] temperature, double[] dewpoint,
			ParcelPath pathType, boolean doEntrainment) {
		ArrayList<RecordAtLevel> parcelPath = new ArrayList<>();
		
		double entrainmentRate = 0.0;
		if(doEntrainment) {
			entrainmentRate = 0.00006; // m^-1
		}

		final double ITER_HEIGHT_CHANGE = 20; // change per iteration

		double height = 0;

		double parcelPressure = -1024.0;
		double parcelTemperature = -1024.0;
		double parcelDewpoint = -1024.0;

		switch (pathType) {
		case SURFACE_BASED:
			parcelPressure = pressure[dewpoint.length - 1];
			parcelTemperature = temperature[dewpoint.length - 1];
			parcelDewpoint = dewpoint[dewpoint.length - 1];

			break;
		case MIXED_LAYER_50MB:
			return computeMixedLayerParcelPath(pressure, temperature, dewpoint, 5000, doEntrainment);
		case MIXED_LAYER_100MB:
			return computeMixedLayerParcelPath(pressure, temperature, dewpoint, 10000, doEntrainment);
		case MOST_UNSTABLE:
			double surfPresMinus300Mb = pressure[pressure.length - 1] - 30000.0;

			int maxThetaEIndex = -1;
			double maxThetaEValue = 0.0;

			for (int i = 0; i < Double.min(Double.min(temperature.length, dewpoint.length), pressure.length); i++) {
				if (pressure[i] > surfPresMinus300Mb) {
					double thetaE = WeatherUtils.thetaE(temperature[i], dewpoint[i], pressure[i]);

					if (thetaE > maxThetaEValue) {
						maxThetaEIndex = i;
						maxThetaEValue = thetaE;
					}
				}
			}

			parcelPressure = pressure[maxThetaEIndex];
			parcelTemperature = temperature[maxThetaEIndex];
			parcelDewpoint = dewpoint[maxThetaEIndex];
		}

		{
			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, height);
			parcelPath.add(record);
		}

		while (parcelPressure > 10000) {
			parcelPressure = WeatherUtils.pressureAtHeight(parcelPressure, ITER_HEIGHT_CHANGE, parcelTemperature);
			height += ITER_HEIGHT_CHANGE;

			double ePres = parcelPressure;
			double eTemp = logInterp(pressure, temperature, ePres);
			double eDwpt = logInterp(pressure, dewpoint, ePres);

			if (parcelDewpoint >= parcelTemperature) {
				parcelTemperature -= moistAdiabaticLapseRate(parcelTemperature, parcelPressure)
						* ITER_HEIGHT_CHANGE;

				double envWetBulb = wetBulbTemperature(eTemp, eDwpt, ePres);
				parcelTemperature += -entrainmentRate * (parcelTemperature - envWetBulb) * ITER_HEIGHT_CHANGE;
				
				parcelDewpoint = parcelTemperature;
			} else {
				parcelTemperature -= dryAdiabaticLapseRate * ITER_HEIGHT_CHANGE;
				parcelDewpoint -= dewpointLapseRate * ITER_HEIGHT_CHANGE;
				
				double parcelMixingRatio = mixingRatio(parcelPressure, parcelDewpoint);
				double envMixingRatio = mixingRatio(ePres, eDwpt);
				
				parcelTemperature += -entrainmentRate * (parcelTemperature - eTemp) * ITER_HEIGHT_CHANGE;
				parcelMixingRatio += -entrainmentRate * (parcelMixingRatio - envMixingRatio) * ITER_HEIGHT_CHANGE;

				double parcelVaporPressure = (parcelPressure * parcelMixingRatio) / (0.62197 + parcelMixingRatio);
				parcelDewpoint = 1 / (1 / 273.15
						- (waterVaporGasConstant / latentHeatOfVaporization * Math.log(parcelVaporPressure / 611)));
			}

			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, height);
			parcelPath.add(record);
		}

		return parcelPath;
	}

	/**
	 * Computes the downdraft parcel path
	 * 
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @return <b>parcelPath</b> RecordAtLevel {pressure: Pascals, temperature:
	 *         Kelvins, wetbulb: Kelvins, dewpoint: Kelvins, height: Meters}
	 */
	public static ArrayList<RecordAtLevel> computeDcapeParcelPath(double[] pressure, double[] temperature,
			double[] dewpoint) {
		ArrayList<RecordAtLevel> parcelPath = new ArrayList<>();

		final double ITER_HEIGHT_CHANGE = -20; // change per iteration

		double height = 0;

		double parcelPressure = -1024.0;
		double parcelTemperature = -1024.0;

		int minThetaEIndex = -1;
		double minThetaEValue = 1024.0;

		for (int i = 0; i < Double.min(Double.min(temperature.length, dewpoint.length), pressure.length); i++) {
			if (pressure[i] > 50000) {
				double thetaE = WeatherUtils.thetaE(temperature[i], dewpoint[i], pressure[i]);

				if (thetaE < minThetaEValue) {
					minThetaEIndex = i;
					minThetaEValue = thetaE;
				}
			}
		}

		parcelPressure = pressure[minThetaEIndex];
		parcelTemperature = WeatherUtils.wetBulbTemperature(temperature[minThetaEIndex], dewpoint[minThetaEIndex],
				pressure[minThetaEIndex]);

		{
			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelTemperature, height);
			parcelPath.add(record);
		}

		while (parcelPressure < pressure[pressure.length - 1]) {
			parcelPressure = WeatherUtils.pressureAtHeight(parcelPressure, ITER_HEIGHT_CHANGE, parcelTemperature);

			parcelTemperature -= WeatherUtils.moistAdiabaticLapseRate(parcelTemperature, parcelPressure)
					* ITER_HEIGHT_CHANGE;

			height += ITER_HEIGHT_CHANGE;

			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelTemperature, height);
			parcelPath.add(record);
		}

		return parcelPath;
	}

	/**
	 * Computes the mixed-layer convective parcel path. Entrainment cape not yet
	 * implemented because I either haven't found or haven't read the study that
	 * derived it.
	 * 
	 * 
	 * @param pressure        Units: array of Pascals
	 * @param temperature     Units: array of Kelvins
	 * @param dewpoint        Units: array of Kelvins
	 * @param mixedLayerDepth Units: Pascals
	 * @return <b>parcelPath</b> RecordAtLevel {pressure: Pascals, temperature:
	 *         Kelvins, wetbulb: Kelvins, dewpoint: Kelvins, height: Meters}
	 */
	public static ArrayList<RecordAtLevel> computeMixedLayerParcelPath(double[] pressure, double[] temperature,
			double[] dewpoint, double mixedLayerDepth, boolean doEntrainment) {
		ArrayList<RecordAtLevel> parcelPath = new ArrayList<>();

		final double ITER_HEIGHT_CHANGE = 20; // change per iteration
		
		double entrainmentRate = 0.0;
		if(doEntrainment) {
			entrainmentRate = 0.00006; // m^-1
		}

		double height = 0;

		double parcelPressure = -1024.0;
		double parcelTemperature = -1024.0;
		double parcelDewpoint = -1024.0;

		double[] potentialTemperature = new double[dewpoint.length];
		double[] mixingRatio = new double[dewpoint.length];

		for (int i = 0; i < dewpoint.length; i++) {
			potentialTemperature[i] = WeatherUtils.potentialTemperature(temperature[i], pressure[i]);
			mixingRatio[i] = WeatherUtils.mixingRatio(pressure[i], dewpoint[i]);
		}

		double topOfMixedLayer = pressure[pressure.length - 2] - mixedLayerDepth;

		// computes the starting mixed-layer parcel
		double averageThetaSum = 0.0;
		double averageMixingRatioSum = 0.0;
		double weightSum = 0.0;

		for (int i = dewpoint.length - 1; i > 0; i--) {
			double pressure1 = pressure[i];
			double theta1 = potentialTemperature[i];
			double mixingRatio1 = mixingRatio[i];

			double pressure2 = pressure[i - 1];
			double theta2 = potentialTemperature[i - 1];
			double mixingRatio2 = mixingRatio[i - 1];

			if (pressure2 > topOfMixedLayer) {
				double weight = (pressure1 - pressure2) / 100.0;

				averageThetaSum += weight * (theta1 + theta2) / 2.0;
				averageMixingRatioSum += weight * (mixingRatio1 + mixingRatio2) / 2.0;
				weightSum += weight;
			} else {
				double weight = (pressure1 - topOfMixedLayer) / 100.0;

				double thetaInterp = linearInterp(pressure, potentialTemperature, topOfMixedLayer);
				double mixingRatioInterp = linearInterp(pressure, mixingRatio, topOfMixedLayer);

				double thetaAddend = weight * (theta1 + thetaInterp) / 2.0;

				averageThetaSum += thetaAddend;
				averageMixingRatioSum += weight * (mixingRatio1 + mixingRatioInterp) / 2.0;
				weightSum += weight;
				break;
			}
		}

		double averageTheta = averageThetaSum / weightSum;
		double averageMixingRatio = averageMixingRatioSum / weightSum;

//		System.out.println("theta:\t" + (averageTheta) + " K");
//		System.out.println("w:\t" + averageMixingRatio + " g g^-1");

		parcelPressure = pressure[pressure.length - 1];
		parcelTemperature = averageTheta / Math.pow(100000 / parcelPressure, 0.286);
		double averageVaporPressure = (parcelPressure * averageMixingRatio) / (0.62197 + averageMixingRatio);
		parcelDewpoint = 1 / (1 / 273.15
				- (waterVaporGasConstant / latentHeatOfVaporization * Math.log(averageVaporPressure / 611)));

//		System.out.println(parcelPressure / 100.0 + " mb");
//		System.out.println(parcelTemperature - 273.15 + " C");
//		System.out.println(parcelDewpoint - 273.15 + " C");

		{
			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, height);
			parcelPath.add(record);
		}

		while (parcelPressure > 10000) {
			parcelPressure = WeatherUtils.pressureAtHeight(parcelPressure, ITER_HEIGHT_CHANGE, parcelTemperature);
			height += ITER_HEIGHT_CHANGE;

			double ePres = parcelPressure;
			double eTemp = logInterp(pressure, temperature, ePres);
			double eDwpt = logInterp(pressure, dewpoint, ePres);

			if (parcelDewpoint >= parcelTemperature) {
				parcelTemperature -= moistAdiabaticLapseRate(parcelTemperature, parcelPressure)
						* ITER_HEIGHT_CHANGE;

				double envWetBulb = wetBulbTemperature(eTemp, eDwpt, ePres);
				parcelTemperature += -entrainmentRate * (parcelTemperature - envWetBulb) * ITER_HEIGHT_CHANGE;
				
				parcelDewpoint = parcelTemperature;
			} else {
				parcelTemperature -= dryAdiabaticLapseRate * ITER_HEIGHT_CHANGE;
				parcelDewpoint -= dewpointLapseRate * ITER_HEIGHT_CHANGE;
				
				double parcelMixingRatio = mixingRatio(parcelPressure, parcelDewpoint);
				double envMixingRatio = mixingRatio(ePres, eDwpt);
				
				parcelTemperature += -entrainmentRate * (parcelTemperature - eTemp) * ITER_HEIGHT_CHANGE;
				parcelMixingRatio += -entrainmentRate * (parcelMixingRatio - envMixingRatio) * ITER_HEIGHT_CHANGE;

				double parcelVaporPressure = (parcelPressure * parcelMixingRatio) / (0.62197 + parcelMixingRatio);
				parcelDewpoint = 1 / (1 / 273.15
						- (waterVaporGasConstant / latentHeatOfVaporization * Math.log(parcelVaporPressure / 611)));
			}

			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, height);
			parcelPath.add(record);
		}

		return parcelPath;
	}

	/**
	 * Computes the convective available potential energy in the Assumes that all
	 * arrays are sorted in order of increasing pressure.
	 * 
	 * @param envPressure       Units: array of Pascals
	 * @param envTemperature    Units: array of Kelvins
	 * @param envDewpoint       Units: array of Kelvins
	 * @param parcelPressure    Units: array of Pascals
	 * @param parcelTemperature Units: array of Kelvins
	 * @param parcelDewpoint    Units: array of Kelvins
	 * @return <b>potentialEnergy</b> Units: J kg^-1
	 */
	public static double computeThreeCape(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			ParcelPath pathType) {
		ArrayList<RecordAtLevel> parcelPath = computeParcelPath(envPressure, envTemperature, envDewpoint, pathType,
				false);

		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		double potentialEnergy = 0.0; // J kg^-1

//		double surfacePressure = envPressure[envPressure.length - 1]; // Pascals

		double equilibriumLevel = equilibriumLevel(envPressure, envTemperature, envDewpoint, parcelPressure,
				parcelHeight, parcelTemperature, parcelDewpoint);
		double levelOfFreeConvection = levelOfFreeConvection(envPressure, envTemperature, envDewpoint, parcelPressure,
				parcelHeight, parcelTemperature, parcelDewpoint);

		for (int i = 0; i < parcelPressure.length - 1; i++) {
			double pPres = parcelPressure[i];
			double height = parcelHeight[i];

//			System.out.printf("%6.1f\t", height);
//			System.out.printf("%6.1f\t", equilibriumLevel);

			if (height < equilibriumLevel && height > levelOfFreeConvection && pPres >= envPressure[0]) {
				double pTemp = parcelTemperature[i];
				double pDwpt = parcelDewpoint[i];

				double ePres = parcelPressure[i];
				double eTemp = logInterp(envPressure, envTemperature, ePres);
				double eDwpt = logInterp(envPressure, envDewpoint, ePres);

				double eVirtTemp = virtualTemperature(eTemp, eDwpt, ePres);
				double pVirtTemp = virtualTemperature(pTemp, pDwpt, pPres);

				double energyAdded = gravAccel * (pVirtTemp - eVirtTemp) / eVirtTemp;
//				System.out.println(dz);
//				System.out.printf("%6.1f\t", ePres/100.0);
//				System.out.print(energyAdded + " energyAdded\t");
//				System.out.print((energyAdded > 0) + " energyAdded\t");

				if (energyAdded > 0 && pDwpt >= pTemp) {
					double dz = parcelHeight[i] - parcelHeight[i + 1];

					if (parcelHeight[i] <= 3000) {
						potentialEnergy += energyAdded * dz;

//						System.out.printf("%7.1f\t%5.1f\t%5.1f\n", ePres / 100.0, eVirtTemp, pVirtTemp);
//						System.out.println(height);
//						System.out.println(energyAdded * dz + " J kg^-1 added");
//						System.out.println(potentialEnergy + " J kg^-1 total");
//						System.out.println();
					}
				} else {
					return potentialEnergy;
				}
			}
		}

		return potentialEnergy;
	}

	/**
	 * Computes the level of free convection given a parcel path. Assumes that all
	 * arrays are sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @param parcelPath     Type: ArrayList of RecordAtLevel, should be type
	 *                       returned by computeParcelPath()
	 * @return <b>levelOfFreeConvection</b> Units: Meters
	 */
	public static double convectiveCondensationLevel(double[] envPressure, double[] envTemperature,
			double[] envDewpoint, ArrayList<RecordAtLevel> parcelPath) {
		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		return convectiveCondensationLevel(envPressure, envTemperature, envDewpoint, parcelPressure, parcelHeight,
				parcelTemperature, parcelDewpoint);
	}

	/**
	 * Computes the level of free convection given the environment and a parcel
	 * path. Assumes that all arrays are sorted by increasing pressure.
	 * 
	 * @param envPressure       Units: Pascals
	 * @param envTemperature    Units: Kelvins
	 * @param envDewpoint       Units: Kelvins
	 * @param parcelPressure    Units: Pascals
	 * @param parcelTemperature Units: Kelvins
	 * @param parcelDewpoint    Units: Kelvins
	 * @return <b>levelOfFreeConvection</b> Units: Meters
	 */
	public static double convectiveCondensationLevel(double[] envPressure, double[] envTemperature,
			double[] envDewpoint, double[] parcelPressure, double[] parcelHeight, double[] parcelTemperature,
			double[] parcelDewpoint) {
		double ccl = -1024.0; // meters

		double surfaceDewpoint = parcelDewpoint[parcelDewpoint.length - 1];

		for (int i = parcelPressure.length - 1; i >= 0; i--) {

			double ePres = parcelPressure[i];
			double eTemp = logInterp(envPressure, envTemperature, ePres);

			double dewpointAtHeight = surfaceDewpoint - dewpointLapseRate * parcelHeight[i];

			if (dewpointAtHeight >= eTemp) {
				return parcelHeight[i];
			}
		}

		return ccl;
	}

	/**
	 * Computes the Corfidi Downshear vector, which estimates the motion of a
	 * bow echo MCS propelled by a cold pool.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param height      Units: array of Meters
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @return <b>corfidiUpshear</b> Units: array of m s^-1
	 */
	public static double[] corfidiDownshear(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind) {
		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		double inflowLayerBase = effectiveInflowLayer(pressure, height, temperature, dewpoint)[0];
		ArrayList<RecordAtLevel> muParcel = computeParcelPath(pressure, temperature, dewpoint, ParcelPath.MOST_UNSTABLE,
				false);
		double muEl = equilibriumLevel(pressure, temperature, dewpoint, muParcel);

		double[] meanWindStormMotion = stormRelativeMeanWind(pressure, height, uWind, vWind, new double[] { 0, 0 },
				inflowLayerBase, muEl);

		double[] coldPoolRelativeFlow = corfidiUpshear(pressure, height, temperature, dewpoint, uWind, vWind);

		return new double[] { meanWindStormMotion[0] + coldPoolRelativeFlow[0], meanWindStormMotion[1] + coldPoolRelativeFlow[1] };
	}

	/**
	 * Computes the Corfidi Upshear vector, which estimates the motion of a
	 * backbuilding MCS.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param height      Units: array of Meters
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @return <b>corfidiUpshear</b> Units: array of m s^-1
	 */
	public static double[] corfidiUpshear(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind) {
		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		double inflowLayerBase = effectiveInflowLayer(pressure, height, temperature, dewpoint)[0];
		ArrayList<RecordAtLevel> muParcel = computeParcelPath(pressure, temperature, dewpoint, ParcelPath.MOST_UNSTABLE,
				false);
		double muEl = equilibriumLevel(pressure, temperature, dewpoint, muParcel);

		double[] meanWindStormMotion = stormRelativeMeanWind(pressure, height, uWind, vWind, new double[] { 0, 0 },
				inflowLayerBase, muEl);

		double[] lowLevelJet = stormRelativeMeanWind(pressure, height, uWind, vWind, new double[] { 0, 0 }, 0, 2000);

		return new double[] { meanWindStormMotion[0] - lowLevelJet[0], meanWindStormMotion[1] - lowLevelJet[1] };
	}

	/**
	 * Computes the deviant tornado motion vector according to <a href =
	 * 'https://cameronnixonphotography.wordpress.com/research/anticipating-deviant-tornado-motion/'>Dr.
	 * Cameron Nixon's method</a>.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure Units: array of Pascals
	 * @param uWind    Units: array of m s^-1
	 * @param vWind    Units: array of m s^-1
	 * @return <b>devTorMotion</b> Units: array of m s^-1, first entry: uWind,
	 *         second entry: vWind
	 */
	public static double[] deviantTornadoMotion(double[] pressure, double[] height, double[] uWind, double[] vWind) {
		double[] stormMotion = stormMotionBunkersIDRightMoving(pressure, height, uWind, vWind);

		return deviantTornadoMotion(pressure, height, uWind, vWind, stormMotion);
	}

	/**
	 * Computes the deviant tornado motion vector according to <a href =
	 * 'https://cameronnixonphotography.wordpress.com/research/anticipating-deviant-tornado-motion/'>Dr.
	 * Cameron Nixon's method</a>.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure Units: array of Pascals
	 * @param uWind    Units: array of m s^-1
	 * @param vWind    Units: array of m s^-1
	 * @return <b>devTorMotion</b> Units: array of m s^-1, first entry: uWind,
	 *         second entry: vWind
	 */
	public static double[] deviantTornadoMotion(double[] pressure, double[] height, double[] uWind, double[] vWind,
			double[] stormMotion) {
		double[] devTorMotion = new double[2];

		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		double uWindSum = 0.0;
		double vWindSum = 0.0;
		double weightSum = 0.0;

		for (int i = uWind.length - 1; i >= 1; i--) {
			double uWind1 = uWind[i]; // m s^-1
			double vWind1 = vWind[i]; // m s^-1
			double height1 = height[i] - height[height.length - 1]; // meters AGL
//			System.out.println("0-500 MW summing at: " + (height[i] - height[height.length - 1]) + " m");

			double uWind2 = uWind[i - 1]; // m s^-1
			double vWind2 = vWind[i - 1]; // m s^-1
			double height2 = height[i - 1] - height[height.length - 1]; // meters AGL

			if (height2 < 500) {
				height2 = 500;
				uWind2 = linearInterp(heightRev, uWindRev, 500 + height[height.length - 1]);
				vWind2 = linearInterp(heightRev, vWindRev, 500 + height[height.length - 1]);
			}

			double weight = (height2 - height1) / 100.0;

			uWindSum += weight * (uWind1 + uWind2) / 2.0;
			vWindSum += weight * (vWind1 + vWind2) / 2.0;
			weightSum += weight;

			if (height[i - 1] > 500) {
				break;
			}
		}

//		System.out.println(uWindSum);
//		System.out.println(vWindSum);

		double meanUWind0_500m = uWindSum / weightSum;
		double meanVWind0_500m = vWindSum / weightSum;

//		System.out.println(meanUWind0_500m);
//		System.out.println(meanVWind0_500m);

		devTorMotion[0] = (stormMotion[0] + meanUWind0_500m) / 2.0;
		devTorMotion[1] = (stormMotion[1] + meanVWind0_500m) / 2.0;

		return devTorMotion;
	}

	/**
	 * Computes the effective bulk wind difference given an environmental sounding.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_eshr.html
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @param uWind          Units: array of m s^-1
	 * @param vWind          Units: array of m s^-1
	 * @return <b>effectiveBulkWindDifference</b> Units: array of m s^-1, first
	 *         entry: u component, second entry: v component
	 */
	public static double[] effectiveBulkWindDifference(double[] envPressure, double[] envHeight,
			double[] envTemperature, double[] envDewpoint, double[] uWind, double[] vWind) {
		double[] inflowLayer = effectiveInflowLayer(envPressure, envHeight, envTemperature, envDewpoint);

		if (inflowLayer[0] == -1024)
			return new double[] { 0.0, 0.0 };

		ArrayList<RecordAtLevel> muParcel = computeParcelPath(envPressure, envTemperature, envDewpoint,
				ParcelPath.MOST_UNSTABLE, false);
		double muEl = equilibriumLevel(envPressure, envTemperature, envDewpoint, muParcel);

		double lowerPressure = WeatherUtils.pressureAtHeight(envPressure[envPressure.length - 1], inflowLayer[0]);
		double upperPressure = WeatherUtils.pressureAtHeight(envPressure[envPressure.length - 1], muEl / 2.0);

		double uWindLower = logInterp(envPressure, uWind, lowerPressure);
		double vWindLower = logInterp(envPressure, vWind, lowerPressure);

		double uWindUpper = logInterp(envPressure, uWind, upperPressure);
		double vWindUpper = logInterp(envPressure, vWind, upperPressure);

		return new double[] { uWindUpper - uWindLower, vWindUpper - vWindLower };
	}

	/**
	 * Computes the effective inflow layer given an environmental sounding.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/soundings/help/effective.html
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @return <b>inflowLayer</b> Units: array of Meters, first entry: lower limit,
	 *         second entry: upper limit
	 */
	public static double[] effectiveInflowLayer(double[] envPressure, double[] envHeight, double[] envTemperature,
			double[] envDewpoint) {
		double[] inflowLayer = new double[2];

		double mucape = computeMucape(envPressure, envTemperature, envDewpoint); // J kg^-1
		double mucinh = computeMucinh(envPressure, envTemperature, envDewpoint); // J kg^-1

//		double surfacePressure = envPressure[envPressure.length - 1]; // Pascals

//		System.out.println("WeatherUtils::effectiveInflowLayer() - MUCAPE " + (int) (mucape) + " J kg^-1");
//		System.out.println("WeatherUtils::effectiveInflowLayer() - MUCINH " + (int) (mucinh) + " J kg^-1");

		// checking to make sure there is enough CAPE to create an inflow layer
		if (mucape <= 100 || mucinh <= -250) {
			inflowLayer[0] = -1024.0;
			inflowLayer[1] = -1024.0;

			return inflowLayer;
		} else {
			boolean foundBottomOfLayer = false;

			for (int i = envDewpoint.length - 1; i >= 1; i--) {
				double[] truncEnvPressure = truncateArray(envPressure, i + 1);
				double[] truncEnvHeight = truncateArray(envHeight, i + 1);
				double[] truncEnvTemperature = truncateArray(envTemperature, i + 1);
				double[] truncEnvDewpoint = truncateArray(envDewpoint, i + 1);

//				System.out.println("WeatherUtils::effectiveInflowLayer() - doing CAPE at " + truncEnvPressure[truncEnvPressure.length - 1]/100.0 + " mb");

				// a little messy but should get the job done right with low effort
				double capeAtLevel = computeSbcape(truncEnvPressure, truncEnvTemperature, truncEnvDewpoint);
				double cinhAtLevel = computeSbcinh(truncEnvPressure, truncEnvTemperature, truncEnvDewpoint);

//				System.out.printf("%4d\t%4d\t", (int) capeAtLevel, (int) cinhAtLevel);
//				System.out.println(foundBottomOfLayer);

				if (capeAtLevel > 100 && cinhAtLevel > -250) {
					if (!foundBottomOfLayer) {
						double hgtAtLayer = truncEnvHeight[truncEnvHeight.length - 1] - envHeight[envHeight.length - 1];

						inflowLayer[0] = hgtAtLayer;

						foundBottomOfLayer = true;
					}
				} else {
					if (!foundBottomOfLayer) {
						continue;
					} else {
						double hgtAtLayer = truncEnvHeight[truncEnvHeight.length - 1] - envHeight[envHeight.length - 1];

						inflowLayer[1] = hgtAtLayer;

						return inflowLayer;
					}
				}
			}

			// statement should never be reached, is present so code will compile
			return inflowLayer;
		}
	}

	/**
	 * Computes the equilibrium given a parcel path. Assumes that all arrays are
	 * sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @param parcelPath     Type: ArrayList of RecordAtLevel, should be type
	 *                       returned by computeParcelPath()
	 * @return <b>liftedCondensationLevel</b> Units: Meters
	 */
	public static double equilibriumLevel(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			ArrayList<RecordAtLevel> parcelPath) {
		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		return equilibriumLevel(envPressure, envTemperature, envDewpoint, parcelPressure, parcelHeight,
				parcelTemperature, parcelDewpoint);
	}

	/**
	 * Computes the equilibrium level given the environment and a parcel path.
	 * Assumes that all arrays are sorted by increasing pressure.
	 * 
	 * @param envPressure       Units: Pascals
	 * @param envTemperature    Units: Kelvins
	 * @param envDewpoint       Units: Kelvins
	 * @param parcelPressure    Units: Pascals
	 * @param parcelTemperature Units: Kelvins
	 * @param parcelDewpoint    Units: Kelvins
	 * @return <b>liftedCondensationLevel</b> Units: Meters
	 */
	public static double equilibriumLevel(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			double[] parcelPressure, double[] parcelHeight, double[] parcelTemperature, double[] parcelDewpoint) {
		double el = -1024.0; // meters

		for (int i = 0; i < parcelPressure.length; i++) {
			double pPres = parcelPressure[i];
			double pTemp = parcelTemperature[i];
			double pDwpt = parcelDewpoint[i];

			double ePres = parcelPressure[i];
			double eTemp = logInterp(envPressure, envTemperature, ePres);
			double eDwpt = logInterp(envPressure, envDewpoint, ePres);

			double eVirtTemp = virtualTemperature(eTemp, eDwpt, ePres);
			double pVirtTemp = virtualTemperature(pTemp, pDwpt, pPres);

//			System.out.printf("%7.1f\t%5.1f\t%5.1f\n", ePres/100.0, eVirtTemp, pVirtTemp);

			if (pVirtTemp >= eVirtTemp) {
//				System.out.println(ePres + " Pa");
				return parcelHeight[i];
			}
		}

		return el;
	}

	/**
	 * Computes the level of free convection given a parcel path. Assumes that all
	 * arrays are sorted in order of increasing pressure.
	 * 
	 * @param envPressure    Units: array of Pascals
	 * @param envTemperature Units: array of Kelvins
	 * @param envDewpoint    Units: array of Kelvins
	 * @param parcelPath     Type: ArrayList of RecordAtLevel, should be type
	 *                       returned by computeParcelPath()
	 * @return <b>levelOfFreeConvection</b> Units: Meters
	 */
	public static double levelOfFreeConvection(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			ArrayList<RecordAtLevel> parcelPath) {
		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		return levelOfFreeConvection(envPressure, envTemperature, envDewpoint, parcelPressure, parcelHeight,
				parcelTemperature, parcelDewpoint);
	}

	/**
	 * Computes the level of free convection given the environment and a parcel
	 * path. Assumes that all arrays are sorted by increasing pressure.
	 * 
	 * @param envPressure       Units: Pascals
	 * @param envTemperature    Units: Kelvins
	 * @param envDewpoint       Units: Kelvins
	 * @param parcelPressure    Units: Pascals
	 * @param parcelTemperature Units: Kelvins
	 * @param parcelDewpoint    Units: Kelvins
	 * @return <b>levelOfFreeConvection</b> Units: Meters
	 */
	public static double levelOfFreeConvection(double[] envPressure, double[] envTemperature, double[] envDewpoint,
			double[] parcelPressure, double[] parcelHeight, double[] parcelTemperature, double[] parcelDewpoint) {
		double lfc = 0.0; // meters
		double el = equilibriumLevel(envPressure, envTemperature, envDewpoint, parcelPressure, parcelHeight,
				parcelTemperature, parcelDewpoint);

		for (int i = 0; i < parcelPressure.length; i++) {
			double height = parcelHeight[i];

			if (height < el) {
				double pPres = parcelPressure[i];
				double pTemp = parcelTemperature[i];
				double pDwpt = parcelDewpoint[i];

				double ePres = parcelPressure[i];
				double eTemp = logInterp(envPressure, envTemperature, ePres);
				double eDwpt = logInterp(envPressure, envDewpoint, ePres);

				double eVirtTemp = virtualTemperature(eTemp, eDwpt, ePres);
				double pVirtTemp = virtualTemperature(pTemp, pDwpt, pPres);

//				System.out.printf("%7.1f\t%5.1f\t%5.1f\n", ePres/100.0, eVirtTemp, pVirtTemp);

				if (pVirtTemp <= eVirtTemp) {
//					System.out.println(ePres + " Pa");
//					System.out.println(parcelHeight[i] + " m");
					return parcelHeight[i];
				}
			}
		}

		return lfc;
	}

	/**
	 * Computes the level of free convection given a parcel path. Assumes that all
	 * arrays are sorted in order of increasing pressure.
	 * 
	 * @param surfacePressure Units: Pascals
	 * @param parcelPath      Type: ArrayList of RecordAtLevel, should be type
	 *                        returned by computeParcelPath()
	 * @return <b>liftedCondensationLevel</b> Units: Meters
	 */
	public static double liftedCondensationLevel(double surfacePressure, ArrayList<RecordAtLevel> parcelPath) {
		double[] parcelPressure = new double[parcelPath.size()];
		double[] parcelHeight = new double[parcelPath.size()];
		double[] parcelTemperature = new double[parcelPath.size()];
		double[] parcelDewpoint = new double[parcelPath.size()];

		for (int i = 0; i < parcelPath.size(); i++) {
			parcelPressure[i] = parcelPath.get(parcelPath.size() - 1 - i).pressure;
			parcelHeight[i] = parcelPath.get(parcelPath.size() - 1 - i).height;
			parcelTemperature[i] = parcelPath.get(parcelPath.size() - 1 - i).temperature;
			parcelDewpoint[i] = parcelPath.get(parcelPath.size() - 1 - i).dewpoint;
		}

		return liftedCondensationLevel(surfacePressure, parcelPressure, parcelHeight, parcelTemperature,
				parcelDewpoint);
	}

	/**
	 * Computes the height above the surface at which the air in a parcel condenses.
	 * Assumes that all arrays are sorted by increasing pressure.
	 * 
	 * @param surfacePressure   Units: Pascals
	 * @param parcelPressure    Units: Pascals
	 * @param parcelTemperature Units: Kelvins
	 * @param parcelDewpoint    Units: Kelvins
	 * @return <b>liftedCondensationLevel</b> Units: Meters
	 */
	public static double liftedCondensationLevel(double surfacePressure, double[] parcelPressure, double[] parcelHeight,
			double[] parcelTemperature, double[] parcelDewpoint) {
		double lcl = -1024.0; // meters

		for (int i = parcelPressure.length - 1; i >= 0; i--) {
//			System.out.printf("%6.1f\t%6.1f\t%5.1f\t%5.1f\n", parcelHeight[i], parcelPressure[i], parcelTemperature[i], parcelDewpoint[i]);
			if (parcelDewpoint[i] >= parcelTemperature[i]) {
//				System.out.println("LCL:\t" + parcelHeight[i]);
				return parcelHeight[i];
			}
		}

		return lcl;
	}

	/**
	 * Moist Adiabatic Lapse Rate intended for use in computing parcel paths. Found
	 * it in my old code, actually no clue how it works and only a faint idea where
	 * I got it.
	 * 
	 * @param temperature Units: Kelvins
	 * @param pressure    Units: Pascals
	 * @return <b>moistAdiabaticLapseRate</b> Units: K m^-1
	 */
	public static double moistAdiabaticLapseRate(double temperature, double pressure) {
		double mixingRatio = mixingRatio(pressure, temperature);

		double lapseRate = 9.81 * (1 + ((2501000 * mixingRatio) / (287 * temperature)))
				/ (1003.5 + 10000 * ((250100 * 2501 * mixingRatio) / (461.5 * temperature * temperature)));

		return lapseRate;
	}

	/**
	 * Computes the storm motion according to <a href=
	 * 'https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2'>Bunkers et.
	 * al. 2000</a>. Mean wind component.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure Units: array of Pascals
	 * @param uWind    Units: array of m s^-1
	 * @param vWind    Units: array of m s^-1
	 * @return <b>stormMotion</b> Units: array of m s^-1, first entry: uWind, second
	 *         entry: vWind
	 */
	public static double[] stormMotionBunkersIDMeanWindComponent(double[] pressure, double[] height, double[] uWind,
			double[] vWind) {
		double[] stormMotion = new double[2];

		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		double uWindSum = 0.0;
		double vWindSum = 0.0;
		double weightSum = 0.0;

		for (int i = uWind.length - 1; i >= 1; i--) {
			double uWind1 = uWind[i]; // m s^-1
			double vWind1 = vWind[i]; // m s^-1
			double height1 = height[i]; // meters

			double uWind2 = uWind[i - 1]; // m s^-1
			double vWind2 = vWind[i - 1]; // m s^-1
			double height2 = height[i - 1]; // meters

			if (height2 < 6000) {
				height2 = 6000;
				uWind2 = linearInterp(heightRev, uWindRev, 6000);
				vWind2 = linearInterp(heightRev, vWindRev, 6000);
			}

			double weight = (height2 - height1) / 100.0;

			uWindSum += weight * (uWind1 + uWind2) / 2.0;
			vWindSum += weight * (vWind1 + vWind2) / 2.0;
			weightSum += weight;

//			System.out.println("uWindSum:\t" + uWindSum);
//			System.out.println("vWindSum:\t" + vWindSum);

			if (height[i - 1] > 6000) {
				break;
			}
		}

		stormMotion[0] = uWindSum / weightSum;
		stormMotion[1] = vWindSum / weightSum;

//		System.out.println("storm motion mw: " + Arrays.toString(stormMotion));

		return stormMotion;
	}

	/**
	 * Computes the storm motion according to <a href=
	 * 'https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2'>Bunkers et.
	 * al. 2000</a>. Deviation component. Added to right-movers, subtracted from
	 * left-movers.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure Units: array of Pascals
	 * @param uWind    Units: array of m s^-1
	 * @param vWind    Units: array of m s^-1
	 * @return <b>stormMotion</b> Units: array of m s^-1, first entry: uWind, second
	 *         entry: vWind
	 */
	public static double[] stormMotionBunkersIDDeviationComponent(double[] pressure, double[] height, double[] uWind,
			double[] vWind) {
		double[] stormMotionDev = new double[2];

		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		double uWind0km = uWind[uWind.length - 1];
		double vWind0km = vWind[vWind.length - 1];

		double uWind0_5km = linearInterp(heightRev, uWindRev, 500 + heightRev[0]);
		double vWind0_5km = linearInterp(heightRev, vWindRev, 500 + heightRev[0]);

		double uWind5_5km = linearInterp(heightRev, uWindRev, 5500 + heightRev[0]);
		double vWind5_5km = linearInterp(heightRev, vWindRev, 5500 + heightRev[0]);

		double uWind6km = linearInterp(heightRev, uWindRev, 6000 + heightRev[0]);
		double vWind6km = linearInterp(heightRev, vWindRev, 6000 + heightRev[0]);

		double uWindHead = (uWind5_5km + uWind6km) / 2.0;
		double vWindHead = (vWind5_5km + vWind6km) / 2.0;

		double uWindTail = (uWind0km + uWind0_5km) / 2.0;
		double vWindTail = (vWind0km + vWind0_5km) / 2.0;

		double uBulkShear = uWindHead - uWindTail;
		double vBulkShear = vWindHead - vWindTail;

//		System.out.println("wind 0km: " + uWind0km + "\t" + vWind0km);
//		System.out.println("wind 0.5km: " + uWind0_5km + "\t" + vWind0_5km);
//		System.out.println("wind 5.5km: " + uWind5_5km + "\t" + vWind5_5km);
//		System.out.println("wind 6km: " + uWind6km + "\t" + vWind6km);
//		System.out.println("head: " + uWindHead + "\t" + vWindHead);
//		System.out.println("tail: " + uWindTail + "\t" + vWindTail);
//		System.out.println("bulk shear: " + uBulkShear + "\t" + vBulkShear);
//		System.out.println();

		double uBulkShearRot = vBulkShear;
		double vBulkShearRot = -uBulkShear;

		double bulkShearMag = Math.hypot(uBulkShearRot, vBulkShearRot);

		uBulkShearRot /= bulkShearMag;
		vBulkShearRot /= bulkShearMag;

		final double DEVIATION_MAGNITUDE = 7.5; // m s^-1

		stormMotionDev[0] = DEVIATION_MAGNITUDE * uBulkShearRot;
		stormMotionDev[1] = DEVIATION_MAGNITUDE * vBulkShearRot;

//		System.out.println("storm motion dev: " + Arrays.toString(stormMotionDev));

		return stormMotionDev;
	}

	/**
	 * Computes the storm motion according to <a href=
	 * 'https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2'>Bunkers et.
	 * al. 2000</a>. For left-movers.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure Units: array of Pascals
	 * @param uWind    Units: array of m s^-1
	 * @param vWind    Units: array of m s^-1
	 * @return <b>stormMotion</b> Units: array of m s^-1, first entry: uWind, second
	 *         entry: vWind
	 */
	public static double[] stormMotionBunkersIDLeftMoving(double[] pressure, double[] height, double[] uWind,
			double[] vWind) {
		double[] stormMotionMw = stormMotionBunkersIDMeanWindComponent(pressure, height, uWind, vWind);
		double[] stormMotionDev = stormMotionBunkersIDDeviationComponent(pressure, height, uWind, vWind);

		double[] stormMotion = { stormMotionMw[0] - stormMotionDev[0], stormMotionMw[1] - stormMotionDev[1] };

//		System.out.println("storm motion left: " + Arrays.toString(stormMotion));

		return stormMotion;
	}

	/**
	 * Computes the storm motion according to <a href=
	 * 'https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2'>Bunkers et.
	 * al. 2000</a>. For right-movers.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure Units: array of Pascals
	 * @param uWind    Units: array of m s^-1
	 * @param vWind    Units: array of m s^-1
	 * @return <b>stormMotion</b> Units: array of m s^-1, first entry: uWind, second
	 *         entry: vWind
	 */
	public static double[] stormMotionBunkersIDRightMoving(double[] pressure, double[] height, double[] uWind,
			double[] vWind) {
		double[] stormMotionMw = stormMotionBunkersIDMeanWindComponent(pressure, height, uWind, vWind);
		double[] stormMotionDev = stormMotionBunkersIDDeviationComponent(pressure, height, uWind, vWind);

		double[] stormMotion = { stormMotionMw[0] + stormMotionDev[0], stormMotionMw[1] + stormMotionDev[1] };

//		System.out.println("storm motion right: " + Arrays.toString(stormMotion));

		return stormMotion;
	}

	/**
	 * I UNDERSTOOD THE CONCEPT INCORRECTLY WHEN I WROTE THIS
	 * 
	 * Computes the storm relative helicity using the method of taking the area
	 * between the hodograph, the storm motion vector, the lower limit SR wind
	 * vector, and the upper limit SR wind vector.
	 * 
	 * Method detailed on
	 * https://cameronnixonphotography.wordpress.com/research/the-storm-relative-hodograph/
	 * CTRL+F -> "Storm-Relative Helicity"
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @param stormMotion Units: array of m s^-1 {first: uWind, second: vWind}
	 * @param lowerLimit  Units: Meters
	 * @param upperLimit  Units: Meters
	 * @return <b>stormRelativeHelicity</b> Units: m^2 s^-2
	 */
	public static double stormRelativeHelicityLegacy(double[] pressure, double[] uWind, double[] vWind,
			double[] stormMotion, double lowerLimit, double upperLimit) {
		double stormRelativeHelicity = 0.0; // m^2 s^-2

		double lowerPressure = WeatherUtils.pressureAtHeight(pressure[pressure.length - 1], lowerLimit);
		double upperPressure = WeatherUtils.pressureAtHeight(pressure[pressure.length - 1], upperLimit);

		for (int i = 0; i < uWind.length - 1; i++) {
			double pressure1 = pressure[i];
			double uWind1 = uWind[i];
			double vWind1 = vWind[i];

			double pressure2 = pressure[i + 1];
			double uWind2 = uWind[i + 1];
			double vWind2 = vWind[i + 1];

			double height1 = heightAtPressure(pressure[pressure.length - 1], pressure1);
			double height2 = heightAtPressure(pressure[pressure.length - 1], pressure2);

			if (height1 <= lowerLimit || height2 >= upperLimit) {
				continue;
			} else {
				if (upperLimit < height1 && upperLimit > height2) {
					uWind1 = logInterp(pressure, uWind, upperPressure);
					vWind1 = logInterp(pressure, vWind, upperPressure);
				}

				if (lowerLimit < height1 && lowerLimit > height2) {
					uWind2 = logInterp(pressure, uWind, lowerPressure);
					vWind2 = logInterp(pressure, vWind, lowerPressure);
				}

				double[] x = { stormMotion[0], uWind1, uWind2 };
				double[] y = { stormMotion[1], vWind1, vWind2 };

//				System.out.printf("storm motion: %4.1f\t%4.1f\n", x[0], y[0]);
//				System.out.printf("upper wind:   %4.1f\t%4.1f\n", x[1], y[1]);
//				System.out.printf("lower wind:   %4.1f\t%4.1f\n", x[2], y[2]);
//
//				System.out.println("lowerLimit: " + lowerLimit);
//				System.out.println("upperLimit: " + upperLimit);
//				System.out.println("height1: " + height1);
//				System.out.println("height2: " + height2);
//				System.out.println("triangle: " + triangleArea(x, y));
//				System.out.println("srh-b: \t" + stormRelativeHelicity);
				stormRelativeHelicity += 1.835 * triangleArea(x, y);
//				System.out.println("srh-a: \t" + stormRelativeHelicity);
//				System.out.println();
			}
		}

		return stormRelativeHelicity;
	}

	/**
	 * Computes the storm relative helicity given a storm motion vector.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @param stormMotion Units: array of m s^-1 {first: uWind, second: vWind}
	 * @param lowerLimit  Units: Meters
	 * @param upperLimit  Units: Meters
	 * @return <b>stormRelativeHelicity</b> Units: m^2 s^-2
	 */
	public static double stormRelativeHelicity(double[] pressure, double height[], double[] uWind, double[] vWind,
			double[] stormMotion, double lowerLimit, double upperLimit) {
		double stormRelativeHelicity = 0.0; // m^2 s^-2

		final double ITER_HEIGHT_CHANGE = 20; // change per iteration

		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];

		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		double zSurface = heightRev[0];

		for (double z = lowerLimit; z < upperLimit; z += ITER_HEIGHT_CHANGE) {
			double uWindAtLevel = linearInterp(heightRev, uWindRev, z + zSurface);
			double vWindAtLevel = linearInterp(heightRev, vWindRev, z + zSurface);

			double uWindAtLevelAbove = linearInterp(heightRev, uWindRev, z + zSurface + ITER_HEIGHT_CHANGE);
			double vWindAtLevelAbove = linearInterp(heightRev, vWindRev, z + zSurface + ITER_HEIGHT_CHANGE);

			double uWindSR = uWindAtLevel - stormMotion[0];
			double vWindSR = vWindAtLevel - stormMotion[1];

			double horizontalVorticityU = -(vWindAtLevelAbove - vWindAtLevel) / ITER_HEIGHT_CHANGE;
			double horizontalVorticityV = (uWindAtLevelAbove - uWindAtLevel) / ITER_HEIGHT_CHANGE;

			double integrand = uWindSR * horizontalVorticityU + vWindSR * horizontalVorticityV;

			stormRelativeHelicity += integrand * ITER_HEIGHT_CHANGE;

//			System.out.println(z);
//			System.out.printf("%4.1f\t%4.1f\n", uWindSR, vWindSR);
//			System.out.printf("%5.3f\t%5.3f\n", horizontalVorticityU, horizontalVorticityV);
//			System.out.println(integrand + " m s^-2");
//			System.out.println(stormRelativeHelicity + " m^2 s^-2");
//			System.out.println();
		}

		return stormRelativeHelicity;
	}

	/**
	 * Computes the storm relative mean wind in a given layer
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @param stormMotion Units: array of m s^-1 {first: uWind, second: vWind}
	 * @param lowerLimit  Units: Meters
	 * @param upperLimit  Units: Meters
	 * @return <b>stormRelativeMeanWind</b> Units: array of m s^-1, first entry:
	 *         uWind, second entry: vWind
	 */
	public static double[] stormRelativeMeanWind(double[] pressure, double[] height, double[] uWind, double[] vWind,
			double[] stormMotion, double lowerLimit, double upperLimit) {
		double meanWindUSum = 0.0; // m s^-1
		double meanWindVSum = 0.0; // m s^-1
		double weightSum = 0.0;

		double lowerPressure = WeatherUtils.pressureAtHeight(pressure[pressure.length - 1], lowerLimit);
		double upperPressure = WeatherUtils.pressureAtHeight(pressure[pressure.length - 1], upperLimit);

		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		for (int i = 0; i < uWind.length - 1; i++) {
			double uWind1 = uWind[i];
			double vWind1 = vWind[i];

			double uWind2 = uWind[i + 1];
			double vWind2 = vWind[i + 1];

			double height1 = height[i];
			double height2 = height[i + 1];

			if (height1 <= lowerLimit || height2 >= upperLimit) {
				continue;
			} else {
				if (upperLimit < height1 && upperLimit > height2) {
					uWind1 = logInterp(pressure, uWind, upperPressure);
					vWind1 = logInterp(pressure, vWind, upperPressure);
					height1 = upperLimit;
				}

				if (lowerLimit < height1 && lowerLimit > height2) {
					uWind2 = logInterp(pressure, uWind, lowerPressure);
					vWind2 = logInterp(pressure, vWind, lowerPressure);
					height2 = lowerLimit;
				}

				double weight = height1 - height2;

				meanWindUSum += weight * (uWind1 + uWind2) / 2.0;
				meanWindVSum += weight * (vWind1 + vWind2) / 2.0;
				weightSum += weight;
			}
		}

		double meanWindU = meanWindUSum / weightSum;
		double meanWindV = meanWindVSum / weightSum;

		return new double[] { meanWindU - stormMotion[0], meanWindV - stormMotion[1] };
	}

	/**
	 * Computes the streamwiseness of the vorticity available in the inflow of a
	 * convective storm.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @param stormMotion Units: array of m s^-1 {first: uWind, second: vWind}
	 * @param lowerLimit  Units: Meters
	 * @param upperLimit  Units: Meters
	 * @return <b>streamwisenessOfVorticity</b> Units: Fraction
	 */
	public static double streamwisenessOfVorticity(double[] pressure, double[] height, double[] uWind, double[] vWind,
			double[] stormMotion, double lowerLimit, double upperLimit) {
		double streamwiseVorticity = streamwiseVorticity(pressure, height, uWind, vWind, stormMotion, lowerLimit,
				upperLimit);
		double totalVorticity = totalHorizontalVorticity(pressure, height, uWind, vWind, lowerLimit, upperLimit);

		return streamwiseVorticity / totalVorticity;
	}

	/**
	 * Computes the streamwise vorticity available in the inflow of a convective
	 * storm.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @param stormMotion Units: array of m s^-1 {first: uWind, second: vWind}
	 * @param lowerLimit  Units: Meters
	 * @param upperLimit  Units: Meters
	 * @return <b>streamwiseVorticity</b> Units: s^-1
	 */
	public static double streamwiseVorticity(double[] pressure, double[] height, double[] uWind, double[] vWind,
			double[] stormMotion, double lowerLimit, double upperLimit) {
		double streamwiseVorticityWeightedSum = 0.0; // s^-1
		double weightSum = 0.0;

		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		for (int i = 0; i < uWind.length - 1; i++) {
			double uWind1 = uWind[i];
			double vWind1 = vWind[i];

			double uWind2 = uWind[i + 1];
			double vWind2 = vWind[i + 1];

			double height1 = height[i];
			double height2 = height[i + 1];

			if (height1 <= lowerLimit || height2 >= upperLimit) {
				continue;
			} else {
				if (upperLimit < height1 && upperLimit > height2) {
					uWind1 = linearInterp(heightRev, uWindRev, upperLimit);
					vWind1 = linearInterp(heightRev, vWindRev, upperLimit);
					height1 = upperLimit;
				}

				if (lowerLimit < height1 && lowerLimit > height2) {
					uWind1 = linearInterp(heightRev, uWindRev, lowerLimit);
					vWind1 = linearInterp(heightRev, vWindRev, lowerLimit);
					height2 = lowerLimit;
				}

				double weight = height1 - height2;

				if (weight != 0) {
					double horizontalVorticityU = -(vWind1 - vWind2) / (height1 - height2);
					double horizontalVorticityV = (uWind1 - uWind2) / (height1 - height2);

					double meanWindU = (uWind1 + uWind2) / 2.0;
					double meanWindV = (vWind1 + vWind2) / 2.0;

					double stormInflowU = meanWindU - stormMotion[0];
					double stormInflowV = meanWindV - stormMotion[1];

					double stormInflowNormU = stormInflowU / Math.hypot(stormInflowU, stormInflowV);
					double stormInflowNormV = stormInflowV / Math.hypot(stormInflowU, stormInflowV);

					double streamwiseVorticity = horizontalVorticityU * stormInflowNormU
							+ horizontalVorticityV * stormInflowNormV;

//					System.out.printf("%6.3f\t%6.3f\n", horizontalVorticityU, horizontalVorticityV);
//					System.out.printf("%6.3f\t%6.3f\n", stormInflowNormU, stormInflowNormV);
//					System.out.println("streamwiseVorticityWeightedSum-b: " + streamwiseVorticityWeightedSum);

					streamwiseVorticityWeightedSum += weight * streamwiseVorticity;
					weightSum += weight;

//					System.out.println("streamwiseVorticityWeightedSum-a: " + streamwiseVorticityWeightedSum);
				}
			}
		}

		return streamwiseVorticityWeightedSum / weightSum;
	}

	/**
	 * Computes the total wind shear induced horizontal vorticity available in the
	 * inflow of a convective storm.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @param stormMotion Units: array of m s^-1 {first: uWind, second: vWind}
	 * @param lowerLimit  Units: Meters
	 * @param upperLimit  Units: Meters
	 * @return <b>streamwiseVorticity</b> Units: s^-1
	 */
	public static double totalHorizontalVorticity(double[] pressure, double[] height, double[] uWind, double[] vWind,
			double lowerLimit, double upperLimit) {
		double horizontalVorticityWeightedSum = 0.0; // s^-1
		double weightSum = 0.0;

		double[] heightRev = new double[vWind.length];
		double[] uWindRev = new double[vWind.length];
		double[] vWindRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			uWindRev[heightRev.length - 1 - i] = uWind[i];
			vWindRev[heightRev.length - 1 - i] = vWind[i];
		}

		for (int i = 0; i < uWind.length - 1; i++) {
			double uWind1 = uWind[i];
			double vWind1 = vWind[i];

			double uWind2 = uWind[i + 1];
			double vWind2 = vWind[i + 1];

			double height1 = height[i];
			double height2 = height[i + 1];

			if (height1 <= lowerLimit || height2 >= upperLimit) {
				continue;
			} else {
				if (upperLimit < height1 && upperLimit > height2) {
					uWind1 = linearInterp(heightRev, uWindRev, upperLimit);
					vWind1 = linearInterp(heightRev, vWindRev, upperLimit);
					height1 = upperLimit;
				}

				if (lowerLimit < height1 && lowerLimit > height2) {
					uWind1 = linearInterp(heightRev, uWindRev, lowerLimit);
					vWind1 = linearInterp(heightRev, vWindRev, lowerLimit);
					height2 = lowerLimit;
				}

				double weight = height1 - height2;

				if (weight != 0) {
					double horizontalVorticityU = -(vWind1 - vWind2) / (height1 - height2);
					double horizontalVorticityV = (uWind1 - uWind2) / (height1 - height2);

					double horizontalVorticity = Math.hypot(horizontalVorticityU, horizontalVorticityV);

					horizontalVorticityWeightedSum += weight * horizontalVorticity;
					weightSum += weight;
				}
			}
		}

		return horizontalVorticityWeightedSum / weightSum;
	}

	/**
	 * Section 5 - Composite Indices Equivalent to CompInd.cpp in C++ version of the
	 * library
	 */

	/**
	 * An experimental parameter currently in development by Amelia Urquhart. It
	 * calculates the angle between the ground-relative storm motion vector and the
	 * Deviant Tornado Motion vector (developed by Dr. Cameron Nixon). It then
	 * assigns a value to that angle ranging between -4 and +4. It is primarily
	 * intended for use by storm spotters and weather photographers to keep them
	 * aware of the possibility of direction changes at the end of a tornado's life,
	 * but may also be useful for some operational forecasting purposes. <br>
	 * <br>
	 * 
	 * If the significant tornado parameter is less than 0.5, DEVTOR is set to
	 * 0.<br>
	 * <br>
	 * 
	 * Interpretation:<br>
	 * <br>
	 * 
	 * -4 - Reversal <br>
	 * -3 - Extreme Right Deviation (135 degrees) <br>
	 * -2 - Major Right Deviation (90 degrees) <br>
	 * -1 - Significant Right Deviation (45 degrees) <br>
	 * +0 - No Significant Deviation Expected <br>
	 * +1 - Significant Left Deviation (45 degrees) <br>
	 * +2 - Major Left Deviation (90 degrees) <br>
	 * +3 - Extreme Left Deviation (135 degrees) <br>
	 * +4 - Reversal
	 * 
	 * @param sigtor               Units: dimensionless
	 * @param stormMotion          Units: array of m s^-1
	 * @param deviantTornadoMotion Units: array of m s^-1
	 * @return <b>deviantTornadoParameter</b> Units: dimensionless
	 */
	public static double deviantTornadoParameter(double sigtor, double[] stormMotion, double[] deviantTornadoMotion) {
		if (Math.abs(sigtor) < 0.5)
			return 0.0;

		double stormMotionMagnitude = Math.hypot(stormMotion[0], stormMotion[1]);
		double deviantTornadoMotionMagnitude = Math.hypot(deviantTornadoMotion[0], deviantTornadoMotion[1]);

		double[] stormMotionNorm = new double[2];
		double[] deviantTornadoMotionNorm = new double[2];

		stormMotionNorm[0] = stormMotion[0] / stormMotionMagnitude;
		stormMotionNorm[1] = stormMotion[1] / stormMotionMagnitude;

		deviantTornadoMotionNorm[0] = deviantTornadoMotion[0] / deviantTornadoMotionMagnitude;
		deviantTornadoMotionNorm[1] = deviantTornadoMotion[1] / deviantTornadoMotionMagnitude;

		double dotProduct = stormMotionNorm[0] * deviantTornadoMotionNorm[0]
				+ stormMotionNorm[1] * deviantTornadoMotionNorm[1];
		double crossProduct = stormMotionNorm[0] * deviantTornadoMotionNorm[1]
				- stormMotionNorm[1] * deviantTornadoMotionNorm[0];

//		System.out.println("devtor:");
//		System.out.println(dotProduct);
//		System.out.println(crossProduct);

		double angleSegment = Math.abs(Math.toDegrees(Math.acos(dotProduct)));

//		System.out.println(angleSegment);

		double angle = 0;

		if (dotProduct >= 0 && crossProduct >= 0) {
			angle = angleSegment;
		} else if (dotProduct < 0 && crossProduct >= 0) {
			angle = angleSegment;
		} else if (dotProduct >= 0 && crossProduct < 0) {
			angle = -angleSegment;
		} else if (dotProduct < 0 && crossProduct < 0) {
			angle = -angleSegment;
		}

//		System.out.println(angle);

		double deviantTornadoParameter = angle / 45.0;
		return deviantTornadoParameter;
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_stpc.html
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @return <b>significantTornadoParameter</b> Units: dimensionless
	 */

	public static double significantHailParameter(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint, double[] uWind, double[] vWind) {
		ArrayList<RecordAtLevel> muParcel = computeParcelPath(pressure, temperature, dewpoint, ParcelPath.MOST_UNSTABLE,
				false);

		double[] heightRev = new double[vWind.length];
		double[] temperatureRev = new double[vWind.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i];
			temperatureRev[heightRev.length - 1 - i] = temperature[i];
		}

		double mucape = computeCape(pressure, temperature, dewpoint, muParcel);

		RecordAtLevel muParcelOrigin = muParcel.get(muParcel.size() - 1);
		double muParcelMixingRatio = mixingRatio(muParcelOrigin.pressure, muParcelOrigin.dewpoint);

		double height700mb = logInterp(pressure, height, 70000);
		double temperature700mb = linearInterp(heightRev, temperatureRev, height700mb);

		double height500mb = logInterp(pressure, height, 50000);
		double temperature500mb = linearInterp(heightRev, temperatureRev, height500mb);

		double lapseRate700_500 = (temperature700mb - temperature500mb) / (height500mb - height700mb);

		double shear0_6km = WeatherUtils.bulkShearMagnitude(pressure, uWind, vWind, 0, 6000);

		return significantHailParameter(mucape, muParcelMixingRatio, lapseRate700_500, temperature500mb, shear0_6km);
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_sigh.html
	 * 
	 * @param mucape              Units: J kg^-1
	 * @param muParcelMixingRatio Units: g g^-1
	 * @param lapseRate700_500    Units: K m^-1
	 * @param tmp500mb            Units: Kelvins
	 * @param shear0_6km          Units: m s^-1
	 * @return <b>significantHailParameter</b> Units: dimensionless
	 */

	public static double significantHailParameter(double mucape, double muParcelMixingRatio, double lapseRate700_500,
			double tmp500mb, double shear0_6km) {
		muParcelMixingRatio *= 1000;
		lapseRate700_500 *= 1000;

		if (shear0_6km < 7)
			shear0_6km = 7;
		if (shear0_6km > 27)
			shear0_6km = 27;

		if (muParcelMixingRatio < 11)
			muParcelMixingRatio = 11;
		if (muParcelMixingRatio > 13.6)
			muParcelMixingRatio = 13.6;

		double tmp500mbC = tmp500mb - 273.15;
		if (tmp500mbC > -5.5)
			muParcelMixingRatio = -5.5;

		double significantHailParameter = mucape * muParcelMixingRatio * lapseRate700_500 * -tmp500mbC * shear0_6km
				/ 42000000.0;

		if (mucape < 1300)
			significantHailParameter *= (mucape / 1300);
		if (lapseRate700_500 < 5.8)
			significantHailParameter *= (lapseRate700_500 / 5.8);

		return significantHailParameter;
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_stpc.html
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @return <b>significantTornadoParameter</b> Units: dimensionless
	 */

	public static double significantTornadoParameterEffective(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint, double[] uWind, double[] vWind) {
		double[] stormMotion = stormMotionBunkersIDRightMoving(pressure, height, uWind, vWind);

		return significantTornadoParameterEffective(pressure, height, temperature, dewpoint, uWind, vWind, stormMotion);
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_stpc.html
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @return <b>significantTornadoParameter</b> Units: dimensionless
	 */

	public static double significantTornadoParameterEffective(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint, double[] uWind, double[] vWind, double[] stormMotion) {
		double[] inflowLayer = effectiveInflowLayer(pressure, height, temperature, dewpoint);

		if (inflowLayer[0] == -1024.0 || inflowLayer[0] > 2.1)
			return 0.0;

		double mlcape = computeMl100cape(pressure, temperature, dewpoint);

		ArrayList<RecordAtLevel> mlParcel = WeatherUtils.computeParcelPath(pressure, temperature, dewpoint,
				ParcelPath.MIXED_LAYER_100MB, false);
		double mlLcl = WeatherUtils.liftedCondensationLevel(pressure[pressure.length - 1], mlParcel);

		double effectiveSrh = stormRelativeHelicity(pressure, height, uWind, vWind, stormMotion, inflowLayer[0],
				inflowLayer[1]);

		double[] effectiveBwdVector = effectiveBulkWindDifference(pressure, height, temperature, dewpoint, uWind,
				vWind);
		double effectiveBwd = Math.hypot(effectiveBwdVector[0], effectiveBwdVector[1]);

		double mlcinh = computeMl100cinh(pressure, temperature, dewpoint);

		return significantTornadoParameterEffective(mlcape, inflowLayer[0], mlLcl, effectiveSrh, effectiveBwd, mlcinh);
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_stpc.html
	 * 
	 * @param mlcape       Units: J kg^-1
	 * @param inflowBase   Units: meters
	 * @param mlLcl        Units: meters
	 * @param effectiveSrh Units: m^2 s^-2
	 * @param effectiveBwd Units: m s^-1
	 * @param mlcinh       Units: J kg^-1
	 * @return <b>significantTornadoParameter</b> Units: dimensionless
	 */

	public static double significantTornadoParameterEffective(double mlcape, double inflowBase, double mlLcl,
			double effectiveSrh, double effectiveBwd, double mlcinh) {
		if (inflowBase == -1024.0 || inflowBase > 2.1)
			return 0.0;

		double mlcapeTerm = mlcape / 1500.0;
		double mlLclTerm = (2000.0 - mlLcl) / 1000.0;
		double esrhTerm = effectiveSrh / 150.0;
		double ebwdTerm = effectiveBwd / 20.0;
		double mlcinhTerm = (200.0 + mlcinh) / 150.0;

		if (mlLcl > 2000) {
			mlLclTerm = 0.0;
		} else if (mlLcl < 1000) {
			mlLclTerm = 1.0;
		}

		if (effectiveBwd > 30) {
			ebwdTerm = 1.5;
		} else if (effectiveBwd < 12.5) {
			ebwdTerm = 0.0;
		}

		if (mlcinh > -50) {
			mlcinhTerm = 1.0;
		} else if (mlcinh < -200) {
			mlcinhTerm = 0.0;
		}

		double significantTornadoParameter = mlcapeTerm * mlLclTerm * esrhTerm * ebwdTerm * mlcinhTerm;
		return significantTornadoParameter;
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_stpc.html
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @return <b>significantTornadoParameter</b> Units: dimensionless
	 */

	public static double significantTornadoParameterFixed(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint, double[] uWind, double[] vWind) {
		double[] stormMotion = stormMotionBunkersIDRightMoving(pressure, height, uWind, vWind);

		return significantTornadoParameterFixed(pressure, height, temperature, dewpoint, uWind, vWind, stormMotion);
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_stpc.html
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @param uWind       Units: array of m s^-1
	 * @param vWind       Units: array of m s^-1
	 * @return <b>significantTornadoParameter</b> Units: dimensionless
	 */

	public static double significantTornadoParameterFixed(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint, double[] uWind, double[] vWind, double[] stormMotion) {
		double[] inflowLayer = effectiveInflowLayer(pressure, height, temperature, dewpoint);

		if (inflowLayer[0] == -1024.0 || inflowLayer[0] > 2.1)
			return 0.0;

		double sbcape = computeSbcape(pressure, temperature, dewpoint);

		ArrayList<RecordAtLevel> sbParcel = WeatherUtils.computeParcelPath(pressure, temperature, dewpoint,
				ParcelPath.SURFACE_BASED, false);
		double sbLcl = WeatherUtils.liftedCondensationLevel(pressure[pressure.length - 1], sbParcel);

		double srh1km = stormRelativeHelicity(pressure, height, uWind, vWind, stormMotion, inflowLayer[0],
				inflowLayer[1]);

		double sixBwd = bulkShearMagnitude(height, uWind, vWind, 0, 6000);

		double sbcinh = computeSbcinh(pressure, temperature, dewpoint);

		return significantTornadoParameterFixed(sbcape, inflowLayer[0], sbLcl, srh1km, sixBwd, sbcinh);
	}

	/**
	 * Computes the significant tornado parameter using the SPC's formula.
	 * 
	 * Reference: https://www.spc.noaa.gov/exper/mesoanalysis/help/help_stpc.html
	 * 
	 * @param sbcape     Units: J kg^-1
	 * @param inflowBase Units: meters
	 * @param sbLcl      Units: meters
	 * @param srh1km     Units: m^2 s^-2
	 * @param sixBwd     Units: m s^-1
	 * @param sbcinh     Units: J kg^-1
	 * @return <b>significantTornadoParameter</b> Units: dimensionless
	 */

	public static double significantTornadoParameterFixed(double sbcape, double inflowBase, double sbLcl, double srh1km,
			double sixBwd, double sbcinh) {
		if (inflowBase == -1024.0 || inflowBase > 2.1)
			return 0.0;

		double mlcapeTerm = sbcape / 1500.0;
		double mlLclTerm = (2000.0 - sbLcl) / 1000.0;
		double esrhTerm = srh1km / 150.0;
		double ebwdTerm = sixBwd / 20.0;
		double mlcinhTerm = (200.0 + sbcinh) / 150.0;

		if (sbLcl > 2000) {
			mlLclTerm = 0.0;
		} else if (sbLcl < 1000) {
			mlLclTerm = 1.0;
		}

		if (sixBwd > 30) {
			ebwdTerm = 1.5;
		} else if (sixBwd < 12.5) {
			ebwdTerm = 0.0;
		}

		if (sbcinh > -50) {
			mlcinhTerm = 1.0;
		} else if (sbcinh < -200) {
			mlcinhTerm = 0.0;
		}

		double significantTornadoParameter = mlcapeTerm * mlLclTerm * esrhTerm * ebwdTerm * mlcinhTerm;
		return significantTornadoParameter;
	}

	/**
	 * Computes the supercell composite using the SPC's formula.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param mucape       Units: J kg^-1
	 * @param effectiveSrh Units: m^2 s^-2
	 * @param effectiveBwd Units: m s^-1
	 * @param mucinh       Units: J kg^-1
	 * @return <b>supercellComposite</b> Units: dimensionless
	 */

	public static double supercellComposite(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind) {
		double[] stormMotion = stormMotionBunkersIDRightMoving(pressure, height, uWind, vWind);

		return supercellComposite(pressure, height, temperature, dewpoint, uWind, vWind, stormMotion);
	}

	/**
	 * Computes the supercell composite using the SPC's formula.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param mucape       Units: J kg^-1
	 * @param effectiveSrh Units: m^2 s^-2
	 * @param effectiveBwd Units: m s^-1
	 * @param mucinh       Units: J kg^-1
	 * @return <b>supercellComposite</b> Units: dimensionless
	 */

	public static double supercellComposite(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind, double[] stormMotion) {
		double mucape = computeMucape(pressure, temperature, dewpoint);

		double[] inflowLayer = effectiveInflowLayer(pressure, height, temperature, dewpoint);

		double effectiveSrh = stormRelativeHelicity(pressure, height, uWind, vWind, stormMotion, inflowLayer[0],
				inflowLayer[1]);

		double[] effectiveBwdVector = effectiveBulkWindDifference(pressure, height, temperature, dewpoint, uWind,
				vWind);
		double effectiveBwd = Math.hypot(effectiveBwdVector[0], effectiveBwdVector[1]);

		double mucinh = computeMucinh(pressure, temperature, dewpoint);

		return supercellComposite(mucape, effectiveSrh, effectiveBwd, mucinh);
	}

	/**
	 * Computes the supercell composite using the SPC's formula.
	 * 
	 * @param mucape       Units: J kg^-1
	 * @param effectiveSrh Units: m^2 s^-2
	 * @param effectiveBwd Units: m s^-1
	 * @param mucinh       Units: J kg^-1
	 * @return <b>supercellComposite</b> Units: dimensionless
	 */

	public static double supercellComposite(double mucape, double effectiveSrh, double effectiveBwd, double mucinh) {
		double mucapeTerm = mucape / 1000.0;
		double esrhTerm = effectiveSrh / 50.0;
		double ebwdTerm = 0.0;

		if (effectiveBwd > 20) {
			ebwdTerm = 1.0;
		} else if (effectiveBwd > 10) {
			ebwdTerm = effectiveBwd / 20.0;
		}

		double mucinhTerm = 1.0;

		if (mucinh < -40) {
			mucinhTerm = -40 / mucinh;
		}

		double supercellComposite = mucapeTerm * esrhTerm * ebwdTerm * mucinhTerm;
		return supercellComposite;
	}

	/**
	 * Section 5 - Winter Weather Meteorology Equivalent to WinterWxMet.cpp in C++ version
	 * of the library
	 */

	/**
	 * Computes the dendritic growth zone given an environmental sounding.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * Layer between -12C and -17C
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param height      Units: array of Meters
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @return <b>dgzLayer</b> Units: array of Meters, first entry: lower limit,
	 *         second entry: upper limit
	 */
	public static double[] dendriticGrowthZoneLayer(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint) {
		double[] dgzLayer = new double[2];

		for (int i = dewpoint.length - 1; i >= 1; i--) {
			double height1 = height[i];
			double temperature1 = temperature[i];
			
			double height2 = height[i - 1];
			double temperature2 = temperature[i - 1];
			
			if(temperature1 >= 273.15 - 12 && temperature2 < 273.15 - 12) {
				if(temperature1 - temperature2 != 0) {
					double weight2 = (temperature1 - (273.15 - 12))/(temperature1 - temperature2);
					double weight1 = 1 - weight2;
					
					dgzLayer[0] = weight1 * height1 + weight2 * height2 - height[height.length - 1];
				}
			}
			
			if(temperature1 >= 273.15 - 17 && temperature2 < 273.15 - 17) {
				if(temperature1 - temperature2 != 0) {
					double weight2 = (temperature1 - (273.15 - 17))/(temperature1 - temperature2);
					double weight1 = 1 - weight2;
					
					dgzLayer[1] = weight1 * height1 + weight2 * height2 - height[height.length - 1];
				}
			}
		}

		return dgzLayer;
	}


	/**
	 * Computes the freezing level given an environmental sounding.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param pressure    Units: array of Pascals
	 * @param height      Units: array of Meters
	 * @param temperature Units: array of Kelvins
	 * @param dewpoint    Units: array of Kelvins
	 * @return <b>freezingLevel</b> Units: Meters
	 */
	public static double freezingLevel(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint) {
		double freezingLevel = -1024.0;

		for (int i = dewpoint.length - 1; i >= 1; i--) {
			double height1 = height[i];
			double temperature1 = temperature[i];
			
			double height2 = height[i - 1];
			double temperature2 = temperature[i - 1];
			
			if(temperature1 >= 273.15 && temperature2 < 273.15) {
				if(temperature1 - temperature2 != 0) {
					double weight2 = (temperature1 - (273.15))/(temperature1 - temperature2);
					double weight1 = 1 - weight2;
					
					freezingLevel = weight1 * height1 + weight2 * height2;
				}
			}
		}

		return freezingLevel - height[height.length - 1];
	}

	/**
	 * Returns the snow-to-liquid ratio according to the Kuchera method
	 * 
	 * @param temperatures Array of temperatures, Units: Kelvins
	 * @return <b>kucheraRatio:<b/> Units: dimensionless
	 */
	public static double kucheraRatio(double[] temperatures) {
		if (temperatures.length == 0)
			return 10; // assumes 10:1 ratio if null input

		double maxColTemp = temperatures[0];
		for (int i = 1; i < temperatures.length; i++) {
			maxColTemp = Double.max(maxColTemp, temperatures[i]);
		}

		double kucheraRatio = (maxColTemp < 271.16 ? 12 + (271.16 - maxColTemp) : 12 + 2 * (271.16 - maxColTemp));
		return kucheraRatio;
	}

	/**
	 * Returns the snow-to-liquid ratio according to the Kuchera method
	 * 
	 * @param temperatures Array of temperatures, Units: Kelvins
	 * @return <b>kucheraRatio:<b/> Units: dimensionless
	 */
	public static float kucheraRatio(float[] temperatures) {
		if (temperatures.length == 0)
			return 10; // assumes 10:1 ratio if null input

		float maxColTemp = temperatures[0];
		for (int i = 1; i < temperatures.length; i++) {
			maxColTemp = Float.max(maxColTemp, temperatures[i]);
		}

		float kucheraRatio = (maxColTemp < 271.16f ? 12 + (271.16f - maxColTemp) : 12 + 2 * (271.16f - maxColTemp));
		return kucheraRatio;
	}

	/**
	 * Miscellaneous
	 */

	/**
	 * Computes the average parameter over a layer using height above ground level.
	 * 
	 * Assumes all arrays are sorted in order of increasing pressure.
	 * 
	 * @param height     Units: array of Meters
	 * @param parameter  Units: user's choice
	 * @param lowerLimit Units: Meters
	 * @param upperLimit Units: Meters
	 * @return <b>averageParameterOverLayer</b> Units: user's choice
	 */
	public static double averageParameterOverLayer(double[] height, double[] parameter, double lowerLimit,
			double upperLimit) {
		double parameterWeightedSum = 0.0;
		double weightSum = 0.0;

		double surfaceHeight = height[height.length - 1];

		double[] heightRev = new double[parameter.length];
		double[] parameterRev = new double[parameter.length];
		for (int i = 0; i < heightRev.length; i++) {
			heightRev[heightRev.length - 1 - i] = height[i] - surfaceHeight;
			parameterRev[heightRev.length - 1 - i] = parameter[i];
		}

		for (int i = 0; i < parameter.length - 1; i++) {
			double height1 = height[i];
			double parameter1 = parameter[i];

			double height2 = height[i + 1];
			double parameter2 = parameter[i + 1];

			if (height1 <= lowerLimit || height2 >= upperLimit) {
				continue;
			} else {
				if (upperLimit < height1 && upperLimit > height2) {
					parameter1 = linearInterp(heightRev, parameterRev, upperLimit);
					height1 = upperLimit;
				}

				if (lowerLimit < height1 && lowerLimit > height2) {
					parameter2 = linearInterp(heightRev, parameterRev, lowerLimit);
					height2 = lowerLimit;
				}

				double weight = (height1 - height2) / 1000.0;

				if (weight != 0) {
					double prm = (parameter1 + parameter2) / 2.0;
					
//					System.out.println(parameter1);
//					System.out.println(parameter2);
//					System.out.println(prm);
//					System.out.println();

					parameterWeightedSum += weight * prm;
					weightSum += weight;
				}
			}
		}

//		System.out.println("SUM: " + (parameterWeightedSum / weightSum));
		return parameterWeightedSum / weightSum;
	}

	public static double heatIndex(double tmp, double dpt) {
		double rh = WeatherUtils.relativeHumidity(tmp, dpt);

		double t = 1.8 * (tmp - 273.15) + 32;
		double c1 = -42.379;
		double c2 = 2.04901523;
		double c3 = 10.14333127;
		double c4 = -0.22475541;
		double c5 = -0.00683783;
		double c6 = -0.05481717;
		double c7 = 0.00122874;
		double c8 = 0.00085282;
		double c9 = -0.00000199;
		double ret = c1 + c2 * t + c3 * rh + c4 * t * rh + c5 * t * t + c6 * rh * rh + c7 * t * t * rh
				+ c8 * t * rh * rh + c9 * t * t * rh * rh;

		if (tmp > 300 && rh > 0.4) {
			return (ret - 32) * 5.0 / 9.0 + 273.15;
		} else {
			return tmp;
		}
	}

	public static double windChill(double tmp, double wndSpd) {
		double tmpF = 1.8 * (tmp - 273.15) + 32;
		double wMph = 2.24 * wndSpd;

//        System.out.println(tmp);
//        System.out.println(tmpF);
		double ret = 0;
		if (tmpF > 50 || wMph < 1.5) {
			ret = tmpF;
		} else {
			double v1 = ((0.6215 * tmpF));
//            System.out.println(v1);
			double v2 = ((35.75 * Math.pow(wMph, 0.16)));
//            System.out.println(v2);
			double v3 = ((0.4275 * tmpF * Math.pow(wMph, 0.16)));
//            System.out.println(v3);
			ret = 35.74 + v1 - v2 + v3;
//            System.out.println(ret);
		}

		return (ret - 32) * 5.0 / 9.0 + 273.15;
	}

	// private backend methods

	// inputArr assumed to already be sorted and increasing
	private static double linearInterp(double[] inputArr, double[] outputArr, double input) {
		if (input < inputArr[0]) {
			return outputArr[0];
		} else if (input >= inputArr[inputArr.length - 1]) {
			return outputArr[outputArr.length - 1];
		} else {
			for (int i = 0; i < inputArr.length - 1; i++) {
				double input1 = inputArr[i];
				double input2 = inputArr[i + 1];

				if (input == input1) {
					return outputArr[i];
				} else if (input < input2) {
					double output1 = outputArr[i];
					double output2 = outputArr[i + 1];

					double weight1 = (input2 - input) / (input2 - input1);
					double weight2 = (input - input1) / (input2 - input1);

					return output1 * weight1 + output2 * weight2;
				} else {
					continue;
				}
			}

			return -1024.0;
		}
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

	// for use in eff inflow layer calculations
	private static double[] truncateArray(double[] arr, int length) {
		double[] ret = new double[Integer.min(arr.length, length)];

		for (int i = 0; i < ret.length; i++) {
			ret[i] = arr[i];
		}

		return ret;
	}

	// for use in helicity calculations. first point is storm motion
	private static double triangleArea(double[] x, double[] y) {
		double stormMotionU = x[0];
		double stormMotionV = y[0];

		double windUpperU = x[1];
		double windUpperV = y[1];

		double windLowerU = x[2];
		double windLowerV = y[2];

		double basisVectorU = windUpperU - stormMotionU;
		double basisVectorV = windUpperV - stormMotionV;

		double rightAngleU = -basisVectorV / Math.hypot(basisVectorU, basisVectorV);
		double rightAngleV = basisVectorU / Math.hypot(basisVectorU, basisVectorV);

		double diffVectorU = windLowerU - windUpperU;
		double diffVectorV = windLowerV - windUpperV;

		double dotProduct = diffVectorU * rightAngleU + diffVectorV * rightAngleV;

		double area = 0.5 * Math.abs(stormMotionU * (windUpperV - windLowerV) + windUpperU * (windLowerV - stormMotionV)
				+ windLowerU * (stormMotionV - windUpperV));

		return Math.signum(dotProduct) * area;
	}
}
