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
	 * Computes the mixing ratio using total pressure, vapor pressure, and
	 * temperature.
	 * 
	 * @param pressure      Units: Pascals
	 * @param dewpoint      Units: Kelvins
	 * @return <b>mixingRatio</b> Units: Fraction
	 */
	public static double mixingRatio(double pressure, double dewpoint) {
		double vaporPressure = WeatherUtils.vaporPressure(dewpoint);
		double mixingRatio = 0.62197 * (vaporPressure)/(pressure - vaporPressure);

		return mixingRatio;
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
		double scaleHeight = (molarGasConstant * temperature) / (avgMolarMass * 9.81); // Meters

		return seaLevelPres * Math.exp(-heightAboveSeaLevel / scaleHeight);
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

	/** Units: K m^-1 */
	public static final double dryAdiabaticLapseRate = 0.0098;
	/** Units: K m^-1 */
	public static final double dewpointLapseRate = 0.0018;

	/**
	 * Computes the convective parcel path. Entrainment cape not yet implemented
	 * because I either haven't found or haven't read the study that derived it.
	 * 
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

		final double ITER_HEIGHT_CHANGE = 20; // change per iteration

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
			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, -1024.0);
			parcelPath.add(record);
		}

		while (parcelPressure > 10000) {
			parcelPressure = WeatherUtils.pressureAtHeight(parcelPressure, ITER_HEIGHT_CHANGE);

			if (parcelDewpoint >= parcelTemperature) {
				parcelTemperature -= WeatherUtils.moistAdiabaticLapseRate(parcelTemperature, parcelPressure)
						* ITER_HEIGHT_CHANGE;
				parcelDewpoint = parcelTemperature;
			} else {
				parcelTemperature -= dryAdiabaticLapseRate * ITER_HEIGHT_CHANGE;
				parcelDewpoint -= dewpointLapseRate * ITER_HEIGHT_CHANGE;
			}

			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, -1024.0);
			parcelPath.add(record);
		}

		return parcelPath;
	}

	/**
	 * Computes the mixed-layer convective parcel path. Entrainment cape not yet implemented
	 * because I either haven't found or haven't read the study that derived it.
	 * 
	 * 
	 * @param pressure        Units: array of Pascals
	 * @param temperature     Units: array of Kelvins
	 * @param dewpoint        Units: array of Kelvins
	 * @param mixedLayerDepth Units: Pascals
	 * @return <b>parcelPath</b> RecordAtLevel {pressure: Pascals, temperature:
	 *         Kelvins, wetbulb: Kelvins, dewpoint: Kelvins, height: Meters}
	 */
	public static ArrayList<RecordAtLevel> computeMixedLayerParcelPath(double[] pressure, double[] temperature, double[] dewpoint,
			double mixedLayerDepth, boolean doEntrainment) {
		ArrayList<RecordAtLevel> parcelPath = new ArrayList<>();

		final double ITER_HEIGHT_CHANGE = 20; // change per iteration

		double parcelPressure = -1024.0;
		double parcelTemperature = -1024.0;
		double parcelDewpoint = -1024.0;
		
		double[] potentialTemperature = new double[dewpoint.length];
		double[] mixingRatio = new double[dewpoint.length];
		
		for(int i = 0; i < dewpoint.length; i++) {
			potentialTemperature[i] = WeatherUtils.potentialTemperature(temperature[i], pressure[i]);
			mixingRatio[i] = WeatherUtils.mixingRatio(pressure[i], dewpoint[i]);
		}

		double topOfMixedLayer = pressure[pressure.length - 2] - mixedLayerDepth;
		
		// computes the starting mixed-layer parcel
		double averageThetaSum = 0.0;
		double averageMixingRatioSum = 0.0;
		double weightSum = 0.0;
		
		for(int i = dewpoint.length - 1; i > 0; i--) {
			double pressure1 = pressure[i];
			double theta1 = potentialTemperature[i];
			double mixingRatio1 = mixingRatio[i];
			
			double pressure2 = pressure[i - 1];
			double theta2 = potentialTemperature[i - 1];
			double mixingRatio2 = mixingRatio[i - 1];
			
			if(pressure2 > topOfMixedLayer) {
				double weight = (pressure1 - pressure2)/100.0;
				
				averageThetaSum += weight * (theta1 + theta2) / 2.0;
				averageMixingRatioSum += weight * (mixingRatio1 + mixingRatio2) / 2.0;
				weightSum += weight;
			} else {
				double weight = (pressure1 - topOfMixedLayer)/100.0;
				
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
		
		System.out.println("theta:\t" + (averageTheta) + " K");
		System.out.println("w:\t" + averageMixingRatio + " g g^-1");
		
		parcelPressure = pressure[pressure.length - 1];
		parcelTemperature = averageTheta / Math.pow(100000 / parcelPressure, 0.286);
		double averageVaporPressure = (parcelPressure * averageMixingRatio)  / (0.62197 + averageMixingRatio);
		parcelDewpoint = 1 / (1/273.15 - (waterVaporGasConstant/latentHeatOfVaporization * Math.log(averageVaporPressure/611)));
		
		System.out.println(parcelPressure/100.0 + " mb");
		System.out.println(parcelTemperature - 273.15 + " C");
		System.out.println(parcelDewpoint - 273.15 + " C");
		
		{
			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, -1024.0);
			parcelPath.add(record);
		}

		while (parcelPressure > 10000) {
			parcelPressure = WeatherUtils.pressureAtHeight(parcelPressure, ITER_HEIGHT_CHANGE);

			if (parcelDewpoint >= parcelTemperature) {
				parcelTemperature -= WeatherUtils.moistAdiabaticLapseRate(parcelTemperature, parcelPressure)
						* ITER_HEIGHT_CHANGE;
				parcelDewpoint = parcelTemperature;
			} else {
				parcelTemperature -= dryAdiabaticLapseRate * ITER_HEIGHT_CHANGE;
				parcelDewpoint -= dewpointLapseRate * ITER_HEIGHT_CHANGE;
			}

			RecordAtLevel record = new RecordAtLevel(parcelPressure, parcelTemperature, parcelDewpoint, -1024.0);
			parcelPath.add(record);
		}

		return parcelPath;
	}

	/**
	 * Miscellaneous
	 */

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
	
	// private backend methods
	
	// inputArr assumed to already be sorted and increasing
	private static double linearInterp(double[] inputArr, double[] outputArr, double input) {
		if(input < inputArr[0]) {
			return outputArr[0];
		} else if (input >= inputArr[inputArr.length - 1]) {
			return outputArr[outputArr.length - 1];
		} else {
			for(int i = 0; i < inputArr.length - 1; i++) {
				double input1 = inputArr[i];
				double input2 = inputArr[i + 1];
				
				if (input == input1) {
					return outputArr[i];
				} else if (input < input2) {
					double output1 = outputArr[i];
					double output2 = outputArr[i + 1];
					
					double weight1 = (input2 - input)/(input2 - input1);
					double weight2 = (input - input1)/(input2 - input1);
					
					return output1 * weight1 + output2 * weight2;
				} else {
					continue;
				}
			}
			
			return -1024.0;
		}
	}
}
