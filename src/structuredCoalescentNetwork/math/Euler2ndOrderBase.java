package structuredCoalescentNetwork.math;

import java.util.List;

import structuredCoalescentNetwork.distribution.StructuredNetworkEvent;

public interface Euler2ndOrderBase {

    public void setup(int maxSize, int types, double epsilon, double max_step);

    public void init(double[] migration_rates, double[] coalescent_rates, double[] reassortment_rates, int lineages,
			List<Integer> n_segs, int numRecords);

//    public void calculateValues(double duration, double[] p, int length);

	public void calculateValues(double duration, double[] p, StructuredNetworkEvent startEvent, int length);

//    public void initAndcalculateValues(int ratesInterval, int lineages, double duration, double[] p, int length,
//	    List<Integer> n_segs);

	public void initAndcalculateValues(int ratesInterval, int lineages, double duration, double[] p, int length,
			List<Integer> n_segs, StructuredNetworkEvent startEvent);

    public void setUpDynamics(double[][] coalescentRates, double[][] migrationRates, double[][] reassortment_rates,
	    int[][] indicators, double[] nextRateShift);

}
