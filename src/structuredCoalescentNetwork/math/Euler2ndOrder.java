package structuredCoalescentNetwork.math;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.util.FastMath;

import beast.mascot.distribution.Mascot;

public class Euler2ndOrder implements Euler2ndOrderBase {

    double epsilon;
    double max_step;

    double[] migration_rates; // flattened square matrix of migration rates
    int n; // dimension of migration rate matrix and indicators matrix
    int n2 = 2;
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
    public void setup(int maxSize, int states, double epsilon, double max_step) {
	linProbs_tmpdt = new double[maxSize];
	linProbs_tmpddt = new double[maxSize];
	linProbs_tmpdddt = new double[maxSize];

	this.max_step = max_step;
	this.epsilon = epsilon;
	this.types = states;
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

    public Euler2ndOrder(double[] migration_rates, double[] coalescent_rates, double[] reassortment_rates, int lineages,
	    int states, double epsilon, double max_step, List<Integer> n_segs) {
	this.max_step = max_step;
	this.epsilon = epsilon;
	this.migration_rates = migration_rates;
	n = (int) (Math.sqrt(migration_rates.length) + 0.5);
	this.coalescent_rates = coalescent_rates;
	this.reassortment_rates = reassortment_rates;
	this.lineages = lineages;
	this.types = states;
	this.dimension = this.lineages * this.types;
	sumTypes = new double[states];
	tCR = new double[states];
	sumDotTypes = new double[states];
	this.n_segs = n_segs;

	iterations = 0;
    }

    @Override
    public void initAndcalculateValues(int ratesInterval, int lineages, double duration, double[] p, int length,
	    List<Integer> n_segs) {
	double nextRateShiftTime = ratesInterval == nextRateShift.length ? Double.POSITIVE_INFINITY
		: nextRateShift[ratesInterval];
	if (ratesInterval >= nextRateShift.length) {
	    ratesInterval = nextRateShift.length - 1;
	}
	migration_rates = migrationRates[ratesInterval];
	coalescent_rates = coalescentRates[ratesInterval];
	indicators = indicators_[ratesInterval];

	iterations = 0;
	n = (int) (Math.sqrt(migration_rates.length) + 0.5);
	this.lineages = lineages;
	this.dimension = this.lineages * this.types;
	this.n_segs = n_segs;

	sumTypes = new double[types];
	tCR = new double[types];
	sumDotTypes = new double[types];

	calculateValues(duration, p, length);
    }

    @Override
    public void calculateValues(double duration, double[] p, int length) {
	double[] pDot = linProbs_tmpdt;
	double[] pDotDot = linProbs_tmpddt;
	double[] pDotDotDot = linProbs_tmpdddt;
	calculateValues(duration, p, pDot, pDotDot, pDotDotDot, length);
    }

    public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot) {
	calculateValues(duration, p, pDot, pDotDot, pDotDotDot, pDot.length);
    }

    public void calculateValues(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot,
	    int length) {

	if (Mascot.debug && false) {
	    System.err.println(duration);
	    System.err.println("caol " + Arrays.toString(coalescent_rates));
	    System.err.println("imgr " + Arrays.toString(migration_rates));
	    System.err.println("p " + Arrays.toString(p));
	}

	clearArray(pDotDot, length);
	clearArray(pDotDotDot, length);

	while (duration > 0) {
	    iterations++;
	    // pDot = new double[length];
	    clearArray(pDot, length);
	    computeDerivatives(p, pDot, pDotDot, pDotDotDot, length);
	    computeSecondDerivate(p, pDot, pDotDot, length);
	    approximateThirdDerivate(pDotDot, pDotDotDot, length);
	    duration = updateP(duration, p, pDot, pDotDot, pDotDotDot, length - 1);

	    if (iterations > 10000) {
		System.err.println("too many iterations, return negative infinity");
		p[length - 1] = Double.NEGATIVE_INFINITY;
		break;
	    }
	}

    }

    void clearArray(double[] v, int n) {
	for (int i = 0; i < n; i++) {
	    v[i] = 0.0;
	}
    }

    double updateP(double duration, double[] p, double[] pDot, double[] pDotDot, double[] pDotDotDot, int length) {
	final double max_dotdotdot = maxAbs(pDotDotDot, length);

	// double timeStep = FastMath.min(FastMath.pow(epsilon*6/max_dotdotdot, C),
	// FastMath.min(duration, max_step));

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

		if (its > 10000) {
//					System.err.println("cannot find proper time step, skip these parameter values");
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
	double sumReassort = 0;
	double[] reassort = new double[types];
	for (int i = 0; i < lineages; i++) {
	    k = currlin;
	    for (j = 0; j < types; j++) {
		reassort[j] = reassortment_rates[j] * (1 - Math.pow(0.5, n_segs.get(i) - 1));
		sumReassort += p[k] * reassort[j];
		k++;
	    }

	    k = currlin;
	    for (j = 0; j < types; j++) {
		double r = sumReassort - reassort[j];
		pDot[k] += p[k] * r;
		pDotDot[k] = r;
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

    }

    public void approximateThirdDerivate(double[] pDotDot, double[] pDotDotDot, int length) {
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
    }
}
