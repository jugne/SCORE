package score.distribution;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;

import beast.core.Input;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import score.dynamics.ConstantReassortment;
import score.math.Euler2ndOrder;
import score.math.Euler2ndOrderBase;

public class SCORE extends StructuredNetworkDistribution {

    public Input<ConstantReassortment> dynamicsInput = new Input<>("dynamics", "Input of rates",
	    Input.Validate.REQUIRED);
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration", 0.000001);
    public Input<Double> maxStepInput = new Input<>("maxStep", "max step for the RK4 integration",
			0.1);

	public Input<Integer> nRecordsInput = new Input<>("nRecords",
			"maximum number of records to keep per interval for stochastic mapping", 100);

    public int samples;
    public int nrSamples;
    public DoubleMatrix[] nodeStateProbabilities;
    public List<NetworkNode> nodes = new ArrayList<>();

    private int nrLineages;

    private double[] coalescentRates;
    private double[] reassortmentRates;

    // Set up for lineage state probabilities
    List<NetworkEdge> activeLineages;
    private double[] linProbs;
    private double[] linProbsNew;
    private int linProbsLength;
	public int types;

    // check if this is the first calculation
    private int first = 0;

    Euler2ndOrderBase euler;
	public Network network;
	public ConstantReassortment dynamics;
    StructuredNetworkIntervals networkIntervals;
	public List<StructuredNetworkEvent> networkEventList;
    double[] linProbs_tmp;
    int[] parents;

    @Override
    public void initAndValidate() {
	dynamics = dynamicsInput.get();
	networkIntervals = networkIntervalsInput.get();
	network = networkIn.get();
	if (network == null) {
	    network = networkIntervals.networkInput.get();
	}
	networkEventList = networkIntervals.getNetworkEventList();
	nodeStateProbabilities = new DoubleMatrix[network.getInternalNodes().size()];
	nrSamples = network.getLeafNodes().size();
		types = dynamics.getNrTypes();

	int intCount = networkEventList.size();

	parents = new int[intCount];

	activeLineages = new ArrayList<>();

	int MAX_SIZE = intCount * types;
	linProbs_tmp = new double[MAX_SIZE];
	linProbs = new double[MAX_SIZE];
	linProbsNew = new double[MAX_SIZE];

	euler = new Euler2ndOrder();
	euler.setup(MAX_SIZE, types, epsilonInput.get(), maxStepInput.get());
    }

    @Override
	public double calculateLogP() {
	networkIntervals = networkIntervalsInput.get();
	networkEventList = networkIntervals.getNetworkEventList();
	nodeStateProbabilities = new DoubleMatrix[networkIntervals.networkInput.get().getInternalNodes().size()];
	nrSamples = networkIntervals.networkInput.get().getLeafNodes().size();
	nodes = new ArrayList<>(networkIntervals.networkInput.get().getInternalNodes());

	int intCount = networkEventList.size();

	parents = new int[intCount];

	int MAX_SIZE = intCount * types;
	linProbs_tmp = new double[MAX_SIZE];
	linProbs = new double[MAX_SIZE];
	linProbsNew = new double[MAX_SIZE];

	euler.setup(MAX_SIZE, types, epsilonInput.get(), maxStepInput.get());

	// Set up for lineage state probabilities
	activeLineages.clear();

	// Compute likelihood at each integration time and tree event starting at final
	// sampling time and moving backwards
	logP = 0;
	nrLineages = 0;

	linProbsLength = 0;
	int networkInterval = 0, ratesInterval = 0;
	double nextEventTime = 0.0;
	double prevEventTime = 0.0;

	// Time to the next rate shift or event on the tree
	StructuredNetworkEvent nextNetworkEvent = networkEventList.get(networkInterval);
	double nextNetworkEventTime = nextNetworkEvent.time;
	double nextRateShift = dynamics.getInterval(ratesInterval);

	if (first == 0 || !dynamics.areDynamicsKnown()) {
	    setUpDynamics();
	}

	coalescentRates = dynamics.getCoalescentRate(ratesInterval);
	reassortmentRates = dynamics.getReassortmentRate(ratesInterval);
	nrLineages = activeLineages.size();
	linProbsLength = nrLineages * types;


	
	// Calculate the likelihood
	do {
	    nextEventTime = Math.min(nextNetworkEventTime, nextRateShift);
	    if (nextEventTime > 0) { // if true, calculate the interval contribution
				logP += doEuler(prevEventTime, nextEventTime, ratesInterval);
	    }

	    if (nextNetworkEventTime <= nextRateShift) {
		switch (nextNetworkEvent.type) {
		case COALESCENCE:
		    nrLineages--;
		    logP += coalesce(nextNetworkEvent);
		    break;

		case SAMPLE:
		    nrLineages++;
		    sample(nextNetworkEvent);
		    break;

		case REASSORTMENT:
		    logP += reassortment(nextNetworkEvent);
		    nrLineages++;
		    break;
		}

		networkInterval++;
		nextRateShift -= nextNetworkEventTime;
		try {
		    nextNetworkEvent = networkEventList.get(networkInterval);
		    nextNetworkEventTime = nextNetworkEvent.time;

		} catch (Exception e) {
		    break;
		}
	    } else {
		ratesInterval++;
		coalescentRates = dynamics.getCoalescentRate(ratesInterval);

		nextNetworkEventTime -= nextRateShift;
		nextRateShift = dynamics.getInterval(ratesInterval);
	    }
	    prevEventTime = nextEventTime;
	    if (logP == Double.NEGATIVE_INFINITY) {
		return logP;
	    }
	} while (nextNetworkEventTime <= Double.POSITIVE_INFINITY);

	first++;
//	setNodeTypes();
	return logP;
    }

    private void sample(StructuredNetworkEvent event) {

	List<NetworkEdge> incomingLines = event.lineagesAdded;
	int sampleState = 0;
	int newLength = linProbsLength + 1 * types;
	int currPosition = linProbsLength;

	/*
	 * If there is no trait given as Input, the model will simply assume that the
	 * last value of the taxon name, the last value after a _, is an integer that
	 * gives the type of that taxon
	 */
	if (dynamics.typeTraitInput.get() != null) {
	    for (NetworkEdge l : incomingLines) {
		activeLineages.add(l);
		sampleState = dynamics.getValue(l.childNode.getTaxonLabel());
		
	    }

			if (sampleState >= dynamics.getNrTypes()) {
		System.err.println("sample discovered with higher state than dimension");
	    }

	    for (int i = 0; i < types; i++) {
		if (i == sampleState) {
		    linProbs[currPosition] = 1.0;
		    currPosition++;
		} else {
		    linProbs[currPosition] = 0.0;
		    currPosition++;
		}
	    }

	} else {
	    /*
	     * If there is no trait given as Input, the model will simply assume that the
	     * last value of the taxon name, the last value after a _, is an integer that
	     * gives the type of that taxon
	     */
	    for (NetworkEdge l : incomingLines) {
		activeLineages.add(l);
		String sampleID = l.childNode.getTaxonLabel();
		String[] splits = sampleID.split("_");
		sampleState = Integer.parseInt(splits[splits.length - 1]); // samples types (or priors) should
									   // eventually be specified in the XML
	    }
	    for (int i = 0; i < types; i++) {
		if (i == sampleState) {
		    linProbs[currPosition] = 1.0;
		    currPosition++;
		} else {
		    linProbs[currPosition] = 0.0;
		    currPosition++;
		}
	    }
	}
	linProbsLength = newLength;
    }

    private double coalesce(StructuredNetworkEvent event) {
	List<NetworkEdge> coalLines = event.lineagesRemoved;
	if (coalLines.size() > 2) {
	    System.err.println("Unsupported coalescent at non-binary node");
	    System.exit(0);
	}
	if (coalLines.size() < 2) {
	    System.out.println();
	    System.out.println("WARNING: Less than two lineages found at coalescent event!");
	    System.out.println();
	    return Double.NaN;
	}

	// get the indices of the two daughter lineages
	final int daughterIndex1 = activeLineages.indexOf(coalLines.get(0));
	final int daughterIndex2 = activeLineages.indexOf(coalLines.get(1));
	if (daughterIndex1 == -1 || daughterIndex2 == -1) {
	    System.out.println("daughter lineages at coalescent event not found");
	    return Double.NaN;
	}

	DoubleMatrix lambda = DoubleMatrix.zeros(types);

	/*
	 * Calculate the overall probability for two strains to coalesce independent of
	 * the state at which this coalescent event is supposed to happen
	 */
	for (int k = 0; k < types; k++) {
	    Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1 * types + k]
		    * linProbs[daughterIndex2 * types + k];
	    if (!Double.isNaN(pairCoalRate)) {
		lambda.put(k, pairCoalRate);
	    } else {
		return Double.NEGATIVE_INFINITY;
	    }
	}

	// add the new parent lineage as an active lineage
	activeLineages.add(event.lineagesAdded.get(0));

	// get the node state probabilities
	DoubleMatrix pVec = new DoubleMatrix();
	pVec.copy(lambda);
	pVec = pVec.div(pVec.sum());

	nodeStateProbabilities[nodes.indexOf(coalLines.get(0).parentNode)] = pVec;

	int linCount = 0;
	// add all lineages execpt the daughter lineage to the new p array
	for (int i = 0; i <= nrLineages; i++) {
	    if (i != daughterIndex1 && i != daughterIndex2) {
		for (int j = 0; j < types; j++) {
		    linProbsNew[linCount * types + j] = linProbs[i * types + j];
		}
		linCount++;
	    }
	}
	// add the parent lineage
	for (int j = 0; j < types; j++) {
	    linProbsNew[linCount * types + j] = pVec.get(j);
	}
	// set p to pnew
	linProbs = linProbsNew;
	linProbsNew = linProbs;
	linProbsLength = linProbsLength - types;

	// check which index is large such that the removing starts
	// with the one with the larger value
	if (daughterIndex1 > daughterIndex2) {
	    activeLineages.remove(daughterIndex1);
	    activeLineages.remove(daughterIndex2);
	} else {
	    activeLineages.remove(daughterIndex2);
	    activeLineages.remove(daughterIndex1);
	}

	if (lambda.min() < 0.0) {
	    System.err.println("Coalescent probability is: " + lambda.min());
	    return Double.NEGATIVE_INFINITY;
	}

	if (lambda.sum() == 0)
	    return Double.NEGATIVE_INFINITY;
	else
	    return Math.log(lambda.sum());
    }

    private double reassortment(StructuredNetworkEvent event) {
	List<NetworkEdge> reassortLines = event.lineagesAdded;
	if (reassortLines.size() > 2) {
	    System.out.println();
	    System.err.println("WARNING: More than two parent lineages at reassortment event!");
	    System.out.println();
	    return Double.NaN;
	}
	if (reassortLines.size() < 2) {
	    System.out.println();
	    System.err.println("WARNING: Less than two parent lineages at reassortment event!");
	    System.out.println();
	    return Double.NaN;
	}

	if (event.lineagesRemoved.size() > 1) {
	    System.out.println("More than one daughter lineage at reassortment event");
	    return Double.NaN;
	}

	// get the indices of the daughter lineage
	final int daughterIndex = activeLineages.indexOf(event.lineagesRemoved.get(0));
	if (daughterIndex == -1) {
	    System.out.println("Daughter lineage at reassortment event not found");
	    return Double.NaN;
	}

	DoubleMatrix lambda = DoubleMatrix.zeros(types);

	for (int k = 0; k < types; k++) {
	    Double typeProb = reassortmentRates[k] * linProbs[daughterIndex * types + k]
		    * Math.pow(0.5, event.segsSortedLeft) * Math.pow(0.5, (event.segsToSort - event.segsSortedLeft))
		    * 2.0;

	    if (!Double.isNaN(typeProb)) {
		lambda.put(k, typeProb);
	    } else {
		return Double.NEGATIVE_INFINITY;
	    }
	}

	// remove daughter lineage from active lineages
	activeLineages.remove(daughterIndex);

	// add two new parent lineages as an active lineages
	activeLineages.add(event.lineagesAdded.get(0));
	activeLineages.add(event.lineagesAdded.get(1));

	// get the node state probabilities
	DoubleMatrix pVec = new DoubleMatrix();
	pVec.copy(lambda);
	pVec = pVec.div(pVec.sum());

	nodeStateProbabilities[nodes.indexOf(reassortLines.get(0).childNode)] = pVec;

	int linCount = 0;
	// add all lineages execpt the daughter lineage to the new p array
	for (int i = 0; i < nrLineages; i++) {
	    if (i != daughterIndex) {
		for (int j = 0; j < types; j++) {
		    try {
			linProbsNew[linCount * types + j] = linProbs[i * types + j];
		    } catch (Exception e) {
			System.out.println(network.getExtendedNewick());
		    }
		}
		linCount++;
	    }
	}
	// add the parent lineage
	for (int l = 0; l < event.lineagesAdded.size(); l++) {
	    for (int j = 0; j < types; j++) {
		linProbsNew[linCount * types + j] = pVec.get(j);
	    }
	    linCount++;
	}

	// set p to pnew
	linProbs = linProbsNew;
	linProbsNew = linProbs;
	linProbsLength = linProbsLength + types;

	if (lambda.min() < 0.0) {
	    System.err.println("Reassortment probability is: " + lambda.min());
	    return Double.NEGATIVE_INFINITY;
	}

	if (lambda.sum() == 0)
	    return Double.NEGATIVE_INFINITY;
	else
	    return Math.log(lambda.sum());
    }

    private void setUpDynamics() {
	int n = dynamics.getEpochCount();
	double[][] coalescentRates = new double[n][];
	double[][] migrationRates = new double[n][];
	double[][] reassortmentRates = new double[n][];
	int[][] indicators = new int[n][];
	double[] nextRateShift = dynamics.getIntervals();
	for (int i = 0; i < n; i++) {
	    coalescentRates[i] = dynamics.getCoalescentRate(i);
	    migrationRates[i] = dynamics.getBackwardsMigration(i);
	    reassortmentRates[i] = dynamics.getReassortmentRate(i);
	    indicators[i] = dynamics.getIndicators(i);
	}
	dynamics.setDynamicsKnown();
	euler.setUpDynamics(coalescentRates, migrationRates, reassortmentRates, indicators, nextRateShift);
    }

	private double doEuler(double start, double end, int ratesInterval) {
		double duration = end - start;

	if (linProbs_tmp.length != linProbsLength + 1) {
	    linProbs_tmp = new double[linProbsLength + 1];
	}
	System.arraycopy(linProbs, 0, linProbs_tmp, 0, linProbsLength);
	linProbs_tmp[linProbsLength] = 0;

	linProbs[linProbsLength - 1] = 0;

	List<Integer> n_segs = new ArrayList<>();
	for (NetworkEdge l : activeLineages) {
	    n_segs.add(l.hasSegments.cardinality());
	}

		euler.initAndcalculateValues(ratesInterval, nrLineages, duration, linProbs_tmp, linProbsLength + 1,
				n_segs, null);

	System.arraycopy(linProbs_tmp, 0, linProbs, 0, linProbsLength);

	return linProbs_tmp[linProbsLength];
    }
    
    public DoubleMatrix getRootTypes(){
    	return nodeStateProbabilities[nodeStateProbabilities.length-1];
    }
    
    public void setNodeTypes() {    	
    	for (NetworkNode n : nodes) {
    		DoubleMatrix m = nodeStateProbabilities[nodes.indexOf(n)];
    		int i = m.argmax();
    		n.setTypeLabel(dynamics.getStringStateValue(i));
    	}
    	
    }

}
