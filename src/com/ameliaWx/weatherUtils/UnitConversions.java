package com.ameliaWx.weatherUtils;

/**
 * @author Zander Urquhart
 * @apiNote Constructed for use with WeatherUtils.java.
 */
public class UnitConversions {
	/**
	 * Section 1 --- Equivalent to UnitConversions.cpp in C++ version of the library
	 */

	public static double fahrenheitToCelsius(double temperatureF) {
		double temperatureC = (temperatureF - 32) / 1.8;

		return temperatureC;
	}

	public static double celsiusToKelvin(double temperatureC) {
		double temperatureK = temperatureC + 273.15;

		return temperatureK;
	}

	public static double fahrenheitToKelvin(double temperatureF) {
		double temperatureK = celsiusToKelvin(fahrenheitToCelsius(temperatureF));

		return temperatureK;
	}

	public static double kelvinToCelsius(double temperatureK) {
		double temperatureC = temperatureK - 273.15;

		return temperatureC;
	}

	public static double celsiusToFahrenheit(double temperatureC) {
		double temperatureF = 1.8 * temperatureC + 32;

		return temperatureF;
	}

	public static double kelvinToFahrenheit(double temperatureK) {
		double temperatureF = celsiusToFahrenheit(kelvinToCelsius(temperatureK));

		return temperatureF;
	}

	public static double hectopascalsToPascals(double pressureHPa) {
		double pressurePa = pressureHPa / 100.0;

		return pressurePa;
	}

	public static double pascalsToHectopascals(double pressurePa) {
		double pressureHPa = pressurePa * 100.0;

		return pressureHPa;
	}

}
