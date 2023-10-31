package com.ameliaWx.weatherUtils;

public class RadarQC {
	// TODO:
	// IMPLEMENT REGION BASED ALGORITHM FROM PYART (https://github.com/ARM-DOE/pyart/blob/main/pyart/correct/region_dealias.py)

	// algorithm is of my own invention
	// some points inspired by UNRAVEL (https://doi.org/10.1175/JTECH-D-19-0020.1)
	// but hopefully faster
	private double[][][] dealiasVelocity(double[][][] data) {
		double[][][] dealiasedVelocity = new double[data.length][data[0].length][data[0][0].length];

		for (int i = 0; i < data.length; i++) {
			dealiasedVelocity[i] = dealiasVelocityAtTilt(data[i]);
		}

		return dealiasedVelocity;
	}

	// zero-isodop method similar to storm relative motion calc
	// call multiple times if possible to get rid of higher order aliasing such as
	// when the nyquist velocity is extremely low (looking at you KMPX)
	private static final int INNER_MID_GATE = 4000;
	private static final int MID_OUTER_GATE = 8000;

	private double[][] dealiasVelocityAtTilt(double[][] data) {
		double[][] dealiasedVelocity = new double[data.length][data[0].length];

		double nyquistVelocity = -1.0;

		for (int a = 0; a < data.length; a++) {
			for (int r = 0; r < data[a].length; r++) {
				double record = data[a][r];

				dealiasedVelocity[a][r] = record;

				if (Math.abs(record) > nyquistVelocity && record != -1024 && record != -2048) {
					nyquistVelocity = Math.abs(record);
				}
			}
		}

		nyquistVelocity += 0.5;

		// maybe make a max too
		final double NEEDS_DEALIASING_FRACTION = 0.8;
		double[][] needsDealiasingMinRange = new double[data.length][3];
		for (int a = 0; a < data.length; a++) {
			for (int r = 0; r < 3; r++) {
				needsDealiasingMinRange[a][r] = 2000;
			}
		}

		double[][] velSum = new double[data.length][3];
		int[][] velCount = new int[data.length][3];
		for (int a = 0; a < data.length; a++) {
			for (int r = 0; r < data[a].length; r++) {
				int rId = 0;
				if (r > INNER_MID_GATE)
					rId = 1;
				if (r > MID_OUTER_GATE)
					rId = 2;

				double record = data[a][r];

				if (record != -1024 && record != -2048) {
					velSum[a][rId] += record;
					velCount[a][rId]++;

					if (r >= 3) {
						double record2 = data[a][r - 2];
						double record1 = data[a][r - 1];
						double record0 = data[a][r];

						if (record0 != -1024 && record0 != -2048 && record1 != -1024 && record1 != -2048
								&& record2 != -1024 && record2 != -2048) {
							double addend0 = (Math.abs(record0) > NEEDS_DEALIASING_FRACTION * nyquistVelocity) ? 0.25
									: 0;
							double addend1 = (Math.abs(record1) > NEEDS_DEALIASING_FRACTION * nyquistVelocity) ? 0.25
									: 0;
							double addend2 = (Math.abs(record2) > NEEDS_DEALIASING_FRACTION * nyquistVelocity) ? 0.25
									: 0;

							double sum = addend0 + addend1 + addend2;

							if (sum >= 0.5 && r < needsDealiasingMinRange[a][rId]) {
								needsDealiasingMinRange[a][rId] = r - 2;
							}
						}
					}
				}
			}
		}

		for (int a = 0; a < velSum.length; a++) {
			for (int r = 0; r < 3; r++) {
				if (velCount[a][r] != 0) {
					velSum[a][r] /= velCount[a][r];
				} else {
					velSum[a][r] = -1024.0;
				}
			}
		}

		double[][] smoothedSums = new double[velSum.length][3];

		final int KERNEL_SIZE = 180;
		for (int a = KERNEL_SIZE / 2; a < 720 + KERNEL_SIZE / 2; a++) {
			for (int r = 0; r < 3; r++) {
				double sum = 0.0;

				for (int da = -KERNEL_SIZE / 2; da <= KERNEL_SIZE / 2; da++) {
					int _a = (a + da + 720) % 720;

					sum += velSum[_a][r];
				}

				smoothedSums[a % 720][r] = sum / 180.0;
			}
		}

		velSum = smoothedSums;

		for (int a = 0; a < data.length; a++) {
			for (int r = 0; r < data[a].length; r++) {
				if (data[a][r] != -1024 && data[a][r] != -2048) {

					int rId = 0;
					if (r > INNER_MID_GATE)
						rId = 1;
					if (r > MID_OUTER_GATE)
						rId = 2;

					if (data[a][r] > 0 && velSum[a][rId] < 0 && r >= needsDealiasingMinRange[a][rId]) {
						dealiasedVelocity[a][r] = data[a][r] - 2 * nyquistVelocity;
					} else if (data[a][r] < 0 && velSum[a][rId] > 0 && r >= needsDealiasingMinRange[a][rId]) {
						dealiasedVelocity[a][r] = data[a][r] + 2 * nyquistVelocity;
					} else {
						dealiasedVelocity[a][r] = data[a][r];
					}
				}
			}
		}

		// box filter to remove straggler pixels
		final int BOX_SIZE = 5;
		for (int a = 0; a < data.length; a++) {
			for (int r = BOX_SIZE + 300; r < data[a].length - BOX_SIZE; r++) {
				int _a = (a + 720) % 720;

				double record = dealiasedVelocity[_a][r];

				if (record != -1024 && record != -2048) {
					double sum = 0;
					int count = 0;

					for (int da = -BOX_SIZE; da <= BOX_SIZE; da++) {
						for (int dr = -BOX_SIZE; dr <= BOX_SIZE; dr++) {
							int _da = (a + da + 720) % 720;

							double record0 = dealiasedVelocity[_da][r + dr];

							if (record0 != -1024 && record0 != -2048) {
								sum += record0;
								count++;
							}
						}
					}

					if (count != 0) {
						double boxAvg = sum / count;

						if (boxAvg - record > nyquistVelocity) {
							dealiasedVelocity[_a][r] += 2 * nyquistVelocity;
						} else if (boxAvg - record < -nyquistVelocity) {
							dealiasedVelocity[_a][r] -= 2 * nyquistVelocity;
						}
					}
				}
			}
		}

		return dealiasedVelocity;
	}
}
