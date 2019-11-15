package structuredCoalescentNetwork.mapping;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import beast.core.Description;


@Description("")
public class ode implements FirstOrderDifferentialEquations {

	double[] migration_rates;
	double[] coalescent_rates;
	double[] reassortment_rates;
	int n; // dimension of migration rate matrix and indicators matrix
	double probs;
    int lineages;
	int length;
    int types;
    int dimension;
	double[] tCR;
	double[] sumTypes;
    Integer[][] connectivity;
    Integer[][] sums;
    List<List<Integer>> lineage_type;
    List<Integer> n_segs;
    
    boolean belowzero = false;

    // constructor
    public ode(double[] migration_rates, double[] coalescent_rates, double[] reassortment_rates, int lineages,
			int types, List<List<Integer>> lineage_type,
			List<Integer> n_segs) {
    	this.migration_rates = migration_rates;
        this.coalescent_rates = coalescent_rates;
        this.reassortment_rates = reassortment_rates;
        this.lineages = lineages;
        this.types = types;
        belowzero = false;
		this.dimension = lineages * types;
        this.n_segs = n_segs;
        this.lineage_type = lineage_type;
		this.length = lineages * types;
		n = (int) (Math.sqrt(migration_rates.length) + 0.5);
		sumTypes = new double[types];
		tCR = new double[types];

    }

    @Override
	public int getDimension() {
        return dimension;
    }

    
    @Override
	public void computeDerivatives(double t, double[] p, double[] pDot) {

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
//			pDot[length - 1] -= sumCoal;

			k = currlin;
			for (j = 0; j < types; j++) {
				// Calculate the Derivate of p:
				double coal = sumCoal - tCR[j];
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

//		pDot[length - 1] /= 2;

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

//			pDot[length - 1] -= sumReassort;

			k = currlin;
			for (j = 0; j < types; j++) {
				double r = sumReassort - reassort[j];
				pDot[k] += p[k] * r;
				k++;
			}
			currlin += types;
		}

	}

	void clearArray(double[] v, int n) {
		for (int i = 0; i < n; i++) {
			v[i] = 0.0;
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
        
    public static void main(String[] args) throws Exception{
        // 2d test
		double[] migration_rates = { 0.0, 0.5, 0.5, 0.0 };
    	double[] coalescent_rates = {1.0, 1.0};
		double[] reassortment_rates = { 0.5, 0.5 };
		int lineages = 1;
        int types = 2;
        /*
         * 0	1	1	0
         *	1	0	0	1
         *	1	0	0	1
         *	0	1	1	0
         *
         */
        Integer[][] con = {{null,0,0,null},{1,null,null,0},{1,null,null,0},{null,1,1,null}};
        Integer[][] sums = {{2,0},{1,1},{1,1},{0,2}};

		final double BACKWARD_INTEGRATION_MIN_STEP = 1e-100;
		final double BACKWARD_INTEGRATION_MAX_STEP = 0.1;
		final double BACKWARD_INTEGRATION_ABS_TOLERANCE = 1e-50;
		final double BACKWARD_INTEGRATION_REL_TOLERANCE = 1e-5;
		final int RATE_CHANGE_CHECKS_PER_EDGE = 100;
		final double RATE_CHANGE_CHECK_CONVERGENCE = 1e-5;
		final int RATE_CHANGE_MAX_ITERATIONS = 1000;

		FirstOrderIntegrator odeIntegrator = new DormandPrince54Integrator(0.1 * BACKWARD_INTEGRATION_MIN_STEP,
				0.1 * BACKWARD_INTEGRATION_MAX_STEP, BACKWARD_INTEGRATION_ABS_TOLERANCE,
				BACKWARD_INTEGRATION_REL_TOLERANCE);

        List<List<Integer>> lineage_type = new ArrayList<>();
        lineage_type.add(new ArrayList<>(Arrays.asList(0, 0)));
        lineage_type.add(new ArrayList<>(Arrays.asList(1, 0)));
        lineage_type.add(new ArrayList<>(Arrays.asList(0, 1)));
        lineage_type.add(new ArrayList<>(Arrays.asList(1, 1)));
        
		List<Integer> n_segs = new ArrayList<>(Arrays.asList(2));

//		FirstOrderIntegrator integrator = new ClassicalRungeKuttaIntegrator(0.00001);
		FirstOrderDifferentialEquations ode = new ode(migration_rates, coalescent_rates, reassortment_rates, lineages,
				types, lineage_type, n_segs);
		double[] y0 = new double[] { 1, 0 };
		double[] y = new double[2];
//		integrator.integrate(ode, 0, y0, 0.1, y);
		odeIntegrator.integrate(ode, 0, y0, 0.1, y);

		System.out.println("Solution: " + y[0] + " " + y[1]);
    }

}