package score.math;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import score.distribution.StructuredNetworkEvent;

public class Euler2ndOrder implements Euler2ndOrderBase {

    double epsilon;
    double max_step;

    double[] migration_rates; // flattened square matrix of migration rates
    int n; // dimension of migration rate matrix and indicators matrix
    int[] indicators;
    double[] coalescent_rates;
    double[] reassortment_rates;
    double probs;
    int lineages;
    int types;
    int dimension;
    double[] sumTypes;
    double[] tCR;
    double[] sumDotTypes;
    List<Integer> n_segs;
	int subIntervalID;
	double durationCopy;

    int iterations;

    public Euler2ndOrder() {
    };

    @Override
    public void init(double[] migration_rates, double[] coalescent_rates, double[] reassortment_rates, int lineages,
			List<Integer> n_segs) {
	this.migration_rates = migration_rates;
	n = (int) (Math.sqrt(migration_rates.length) + 0.5);
	this.coalescent_rates = coalescent_rates;
	this.reassortment_rates = reassortment_rates;
	this.lineages = lineages;
	this.dimension = this.lineages * this.types;
	sumTypes = new double[types];
	tCR = new double[types];
	sumDotTypes = new double[types];
	this.n_segs = n_segs;

	iterations = 0;

    }

    double[] linProbs_tmpdt;
    double[] linProbs_tmpddt;
    double[] linProbs_tmpdddt;

    @Override
	public void setup(int maxSize, int types, double epsilon, double max_step) {
	linProbs_tmpdt = new double[maxSize];
	linProbs_tmpddt = new double[maxSize];
	linProbs_tmpdddt = new double[maxSize];

	this.max_step = max_step;
	this.epsilon = epsilon;
	this.types = types;
    }

    public double[][] coalescentRates;
    double[][] migrationRates;
    double[][] reassortmentRates;
    int[][] indicators_;
    double[] nextRateShift;

    @Override
    public void setUpDynamics(double[][] coalescentRates, double[][] migrationRates, double[][] reassortmentRates,
	    int[][] indicators, double[] nextRateShift) {
	this.coalescentRates = coalescentRates;
	this.migrationRates = migrationRates;
	this.reassortmentRates = reassortmentRates;
	this.indicators_ = indicators;
	this.nextRateShift = nextRateShift;
    }

    @Override
    public void initAndcalculateValues(int ratesInterval, int lineages, double duration, double[] p, int length,
			List<Integer> n_segs, StructuredNetworkEvent startEvent) {
	double nextRateShiftTime = ratesInterval == nextRateShift.length ? Double.POSITIVE_INFINITY
		: nextRateShift[ratesInterval];
	if (ratesInterval >= nextRateShift.length) {
	    ratesInterval = nextRateShift.length - 1;
	}
	migration_rates = migrationRates[ratesInterval];
	coalescent_rates = coalescentRates[ratesInterval];
	reassortment_rates = reassortmentRates[ratesInterval];
	indicators = indicators_[ratesInterval];

	iterations = 0;

	durationCopy = duration;
	n = (int) (Math.sqrt(migration_rates.length) + 0.5);
	this.lineages = lineages;
	this.dimension = this.lineages * this.types;
	this.n_segs = n_segs;

	sumTypes = new double[types];
	tCR = new double[types];
	sumDotTypes = new double[types];

		if (startEvent != null)
			subIntervalID = startEvent.numRecords;

		calculateValues(duration, p, startEvent, length);
    }


    @Override
	public void calculateValues(double duration, double[] p, StructuredNetworkEvent startEvent, int length) {
	double[] pDot = linProbs_tmpdt;
	double[] pDotDot = linProbs_tmpddt;
	double[] pDotDotDot = linProbs_tmpdddt;

		calculateValues(duration, p, pDot, pDotDot, pDotDotDot, startEvent, length);
    }


    public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot,
			StructuredNetworkEvent startEvent, int length) {


	clearArray(pDotDot, length);
	clearArray(pDotDotDot, length);

	while (duration > 0) {
	    iterations++;
	    clearArray(pDot, length);
	    computeDerivatives(p, pDot, pDotDot, pDotDotDot, length);
	    computeSecondDerivate(p, pDot, pDotDot, length);
	    approximateThirdDerivate(p, pDot, pDotDot, pDotDotDot, length);

			if (startEvent != null
					&& (duration < (durationCopy * subIntervalID) / startEvent.numRecords || iterations == 1)) {
				int pos = startEvent.numRecords - subIntervalID;
				startEvent.intermediateTimeStored[pos] -= duration;
				startEvent.p_stored[pos] = Arrays.copyOfRange(p, 0, p.length);
				startEvent.pDot_stored[pos] = Arrays.copyOfRange(pDot, 0, p.length);
				subIntervalID -= 1;
			}

	    duration = updateP(duration, p, pDot, pDotDot, pDotDotDot, length - 1);

			if (iterations > 100000) {
				System.err.println("too many iterations, return negative infinity");
				p[length - 1] = Double.NEGATIVE_INFINITY;
				break;
	    }
	}

		if (startEvent != null) {
			if (duration < (durationCopy * subIntervalID) / startEvent.numRecords) {
				int pos = startEvent.numRecords - subIntervalID;
				startEvent.intermediateTimeStored[pos] -= duration;
				startEvent.p_stored[pos] = Arrays.copyOfRange(p, 0, p.length);
				startEvent.pDot_stored[pos] = Arrays.copyOfRange(pDot, 0, p.length);
				subIntervalID -= 1;
			}

			startEvent.numRecords -= subIntervalID;
			startEvent.intermediateTimeStored = Arrays.copyOfRange(startEvent.intermediateTimeStored, 0,
					startEvent.numRecords);
			startEvent.p_stored = Arrays.copyOfRange(startEvent.p_stored, 0, startEvent.numRecords);
		}
    }

    void clearArray(double[] v, int n) {
	for (int i = 0; i < n; i++) {
	    v[i] = 0.0;
	}
    }

    double updateP(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
	final double max_dotdotdot = maxAbs(pDotDotDot, length);


	double timeStep = FastMath.min(FastMath.cbrt(epsilon * 6 / max_dotdotdot), FastMath.min(duration, max_step));
	double timeStepSquare = timeStep * timeStep * 0.5;

	for (int i = 0; i < length; i++) {
	    double new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
	    double diff = FastMath.abs(new_val - p[i]);
	    int its = 0;
	    while (new_val > 1 || new_val < 0 || diff > 0.2) {
		timeStep *= 0.9;
		timeStepSquare = timeStep * timeStep * 0.5;
		new_val = p[i] + pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
		diff = FastMath.abs(new_val - p[i]);
		its++;

				if (its > 100000) {
		    p[length - 1] = Double.NEGATIVE_INFINITY;
		    break;
		}
	    }
	}

	if (p[length - 1] == Double.NEGATIVE_INFINITY)
	    return 0.0;

	updateP2(timeStep, timeStepSquare, p, length + 1, pDot, pDotDot);

	// normalize to ensure stability
	for (int i = 0; i < lineages; i++) {
	    normalise(i, p);
	}

	duration -= timeStep;
	return duration;
    }

    double maxAbs(double[] pDotDotDot, int length) {
	double max_dotdotdot = 0.0;
	for (int i = 0; i < length; i++) {
	    max_dotdotdot = FastMath.max(max_dotdotdot, FastMath.abs(pDotDotDot[i]));
	}
	return max_dotdotdot;
    }

    void normalise(final int i, final double[] p) {
	final int k = types * i;
	double linSum = 0;

	int u = k;
	int q;
	double x;

	for (q = 0; q < types; q++) {
	    x = p[u++];
	    linSum += x;
	    if (x < 0.0) {
		p[p.length - 1] = Double.NEGATIVE_INFINITY;
		return;
	    } // XXX

	}
	u = k;
	for (q = 0; q < types; q++) {
	    p[u++] /= linSum;
	}
    }

    void bailout(double[] p) {
	System.err.println(Arrays.toString(p));
	System.exit(0);
    }

    void updateP2(final double timeStep, final double timeStepSquare, final double[] p, final int length,
	    final double[] pDot, final double[] pDotDot) {
	for (int i = 0; i < length; i++) {
	    p[i] += pDot[i] * timeStep + pDotDot[i] * timeStepSquare;
	}
    }

    public void computeDerivatives(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {

	double migrates;
	// Compute the sum of line state probabilities for each state
	clearArray(sumTypes, types);
	calcSumStates(sumTypes, p);

	// Calculate the change in the lineage state probabilities for every lineage in
	// every state
	int currlin = 0, j, k;
	for (int i = 0; i < lineages; i++) {

	    double sumCoal = 0;
	    k = currlin;
	    for (j = 0; j < types; j++) {
		tCR[j] = coalescent_rates[j] * (sumTypes[j] - p[k]);
		sumCoal += p[k] * tCR[j];
		k++;
	    }
	    pDot[length - 1] -= sumCoal;

	    k = currlin;
	    for (j = 0; j < types; j++) {
		// Calculate the Derivate of p:
		double coal = sumCoal - tCR[j];
		pDotDot[k] = coal;
		pDotDotDot[k] = coal;
		pDot[k] += p[k] * coal;
		k++;
	    } // j
	    currlin += types;
	}

	// Calculate the probability of a lineage changing states

	int u = 0, v;
	double pj;
	for (int i = 0; i < lineages; i++) {
	    // Calculate the probability of a lineage changing states
	    for (j = 0; j < types; j++) {
		v = u;
		pj = p[u];
		for (k = j + 1; k < types; k++) {
		    v++;
		    // the probability of lineage i being in state j is p[i*nr_states +j]
		    migrates = p[v] * migration_rates[k * n + j] - pj * migration_rates[j * n + k];
		    pDot[u] += migrates;
		    pDot[v] -= migrates;
		} // j XXX
		u++;
	    } // j
	} // lineages

	pDot[length - 1] /= 2;

	currlin = 0;
	double[] reassort = new double[types];

	for (int i = 0; i < lineages; i++) {

	    double sumReassort = 0;
	    k = currlin;

	    for (j = 0; j < types; j++) {
		reassort[j] = reassortment_rates[j] * (1 - Math.pow(0.5, n_segs.get(i) - 1));
		sumReassort += p[k] * reassort[j];
		k++;
	    }

	    pDot[length - 1] -= sumReassort;

	    k = currlin;
	    for (j = 0; j < types; j++) {
		double r = sumReassort - reassort[j];
		pDot[k] += p[k] * r;
		pDotDot[k] += r;
		pDotDotDot[k] += r;
		k++;
	    }
	    currlin += types;
	}

    }

    private void calcSumStates(final double[] sumStates, final double[] p) {
	int u = 0, j;
	for (int i = 0; i < lineages; i++) {
	    for (j = 0; j < types; j++) {
		sumStates[j] += p[u++];
	    }
	}
    }

    public void computeSecondDerivate(double[] p, double[] pDot, double[] pDotDot, int length) {
	clearArray(sumDotTypes, types);
	calcSumStates(sumDotTypes, pDot);

	// Calculate the change in the lineage state probabilities for every lineage in
	// every state
	int currlin = 0, j;
	for (int i = 0; i < lineages; i++) {
	    double pCoalRate = 0.0;
	    int k = currlin;
	    for (j = 0; j < types; j++) {
		pCoalRate += coalescent_rates[j] * (pDot[k] * (sumTypes[j] - 2 * p[k]) + p[k] * (sumDotTypes[j]));
		k++;
	    }

	    k = currlin;
	    for (j = 0; j < types; j++) {
		pDotDot[k] = pDotDot[k] * pDot[k]
			+ p[k] * (pCoalRate - coalescent_rates[j] * (sumDotTypes[j] - pDot[k]));
		k++;
	    } // j

	    pDotDot[length - 1] -= pCoalRate;

	    currlin += types;
	} // lineages

	double migrates;

	// Calculate the probability of a lineage changing states

	int u = 0;
	for (int i = 0; i < lineages; i++) {
	    // Calculate the probability of a lineage changing states
	    for (j = 0; j < types; j++) {
		double pj = pDot[u];
		int v = u + 1;
		for (int k = j + 1; k < types; k++) {

		    // the probability of lineage i being in state j is p[i*nr_states +j]
		    migrates = pDot[v] * migration_rates[k * n + j] - pj * migration_rates[j * n + k];
		    pDotDot[u] += migrates;
		    pDotDot[v] -= migrates;
		    v++;
		} // j XXX
		u++;
	    } // j
	} // lineages

	pDotDot[length - 1] /= 2;

	currlin = 0;

	double[] reassort = new double[types];
	for (int i = 0; i < lineages; i++) {
	    double sumReassort = 0;
	    int k = currlin;
	    for (j = 0; j < types; j++) {
		reassort[j] = reassortment_rates[j] * (1 - Math.pow(0.5, n_segs.get(i) - 1));
		sumReassort += pDot[k] * reassort[j];
		k++;
	    }

	    pDotDot[length - 1] -= sumReassort;

	    k = currlin;
	    for (j = 0; j < types; j++) {
		pDotDot[k] += p[k] * sumReassort;
		k++;
	    }
	    currlin += types;
	}

    }

    public void approximateThirdDerivate(double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
	double migrates;

	// Calculate the change in the lineage state probabilities for every lineage in
	// every state
	for (int u = 0; u < length - 1; u++) {
	    pDotDotDot[u] *= pDotDot[u];
	}

	// Calculate the probability of a lineage changing states

	int k;
	for (int j = 0; j < types; j++) {
	    for (k = 0; k < types; k++) {
		double mrate = migration_rates[j * n + k];
		int u = j;
		int v = k;
		for (int i = 0; i < lineages; i++) {
		    migrates = pDotDot[u] * mrate;
		    pDotDotDot[v] += migrates;
		    pDotDotDot[u] -= migrates;
		    u += types;
		    v += types;
		} // XXX

	    }
	}

	int currlin = 0, j;
	double[] reassort = new double[types];
	for (int i = 0; i < lineages; i++) {
	    double sumReassort_1 = 0;
	    double sumReassort_2 = 0;
	    k = currlin;
	    for (j = 0; j < types; j++) {
		reassort[j] = reassortment_rates[j] * (1 - Math.pow(0.5, n_segs.get(i) - 1));
		sumReassort_1 += pDot[k] * reassort[j];
		sumReassort_2 += pDotDot[k] * reassort[j];
		k++;
	    }

	    k = currlin;
	    for (j = 0; j < types; j++) {
		pDotDotDot[k] += (2 * pDot[k] * sumReassort_1) + (p[k] * sumReassort_2);
		k++;
	    }
	    currlin += types;
	}
    }
}
