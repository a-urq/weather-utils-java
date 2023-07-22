package com.ameliaWx.weatherUtils;

public class RecordAtLevel implements Comparable<RecordAtLevel> {
	double pressure; // Pascals
	double temperature; // Kelvins
	double wetbulb; // Kelvins
	double dewpoint; // Kelvins
	double height; // Meters

	public RecordAtLevel(double pres, double tmp, double dpt, double hgt) {
		pressure = pres;
		temperature = tmp;
		dewpoint = dpt;
		wetbulb = WeatherUtils.wetBulbTemperature(temperature, dewpoint, pressure);
		height = hgt;
	}

	@Override
	public int compareTo(RecordAtLevel r) {
		return (int) Math.signum(this.pressure - r.pressure);
	}
	
	@Override
	public String toString() {
		return "[p: " + pressure + ",\tt:" + temperature + ",\tw:" + wetbulb + ",\td:" + dewpoint + ",\th:" + height + "]";
	}

	public double getPressure() {
		return pressure;
	}

	public double getTemperature() {
		return temperature;
	}

	public double getWetbulb() {
		return wetbulb;
	}

	public double getDewpoint() {
		return dewpoint;
	}

	public double getHeight() {
		return height;
	}
}
