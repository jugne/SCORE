package structuredCoalescentNetwork.mapping;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.jblas.DoubleMatrix;

import beast.core.Input;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import structuredCoalescentNetwork.distribution.StructuredNetworkEvent;
import structuredCoalescentNetwork.distribution.StructuredNetworkIntervals;
import structuredCoalescentNetwork.dynamics.ConstantReassortment;
import structuredCoalescentNetwork.math.Euler2ndOrder;
import structuredCoalescentNetwork.math.Euler2ndOrderBase;


/**
 * @author Ugne Stolz. Created on 15 Nov 2019. Mostly Uses ideas of stochastic
 *         mapping from Tim Vaughan's derivations for structured coalescent
 */
public class MappedNetwork extends Network {

	public Input<Boolean> mapOnInitInput = new Input<>("mapOnInit",
			"If true, mapping will be performed when object is " + "first initialize.", true);

	public Input<Network> netwokInput = new Input<>("untypedNetwork", "Network on which to apply mapping.",
			Input.Validate.REQUIRED);

	public Input<ConstantReassortment> dynamicsInput = new Input<>("dynamics", "Input of rates",
			Input.Validate.REQUIRED);

	public Input<Integer> nRecordsInput = new Input<>("nRecords",
			"maximum number of records to keep per interval for stochastic mapping", 1000);

	public Input<Boolean> rejectionInput = new Input<>("rejection",
			"If true, mapper will reject simulation if parent lineage types of reassortment event are different. "
					+ "Can significantly increase the runtime",
			true);

	public Input<Boolean> remapOnLogInput = new Input<>("remapOnLog",
			"If true, mapping will be regenerated when this object " +
					"is logged.",
			true);

	public Input<Boolean> inheritReaTypeInput = new Input<>("inheritReaType",
			"Does child lineage of reassortment event inherits the type from at least one parent",
			true);

	// DEBUG purposes
	boolean rejection;
	boolean inheritReaType;


	StructuredNetworkIntervals intervals = new StructuredNetworkIntervals();
	public ConstantReassortment dynamics;

	private Network untypedNetwork;
	private List<StructuredNetworkEvent> eventList;
	private DoubleMatrix[] nodeStateProbabilities;
	public int nrSamples;
	public List<NetworkNode> nodes = new ArrayList<>();
	private int types;
	private int nrLineages;
	private double[] linProbs;
	private double[] linProbsNew;
	double[] linProbs_tmp;
	private int linProbsLength;
	Euler2ndOrderBase euler;

	/**
	 * Maximum number of steps in each waiting time calculation in forward
	 * simulation.
	 */
	private final int FORWARD_INTEGRATION_STEPS = 100;

	private final double STEP_SIZE_BACKWARD_INTEGRATION = 0.000001;

	private final double MAX_STEP_FOR_BACKWARD_INTEGRATION = 0.001;

	List<NetworkEdge> activeLineages;

	private double[] coalescentRates;
	private double[] reassortmentRates;
	int[] parents;


	@Override
	public void initAndValidate() {

		rejection = rejectionInput.get();
		inheritReaType = inheritReaTypeInput.get();

		dynamics = dynamicsInput.get();
		types = dynamics.getNrTypes();

		activeLineages = new ArrayList<>();

		if (mapOnInitInput.get())
			doStochasticMapping();
	}

	public void doStochasticMapping() {
		untypedNetwork = (Network) netwokInput.get().copy();
//		System.out.println("untyped");
//		System.out.println(untypedNetwork.getExtendedNewick());
		this.setRootEdge(untypedNetwork.getRootEdge());

		intervals.initAndValidate(untypedNetwork);
		eventList = intervals.getNetworkEventList(untypedNetwork);

		backwardIntegration();
		Object[] net = forwardSimulateNetwork();

		while (!(Boolean) net[0]) {
			untypedNetwork = (Network) netwokInput.get().copy();

			intervals.initAndValidate(untypedNetwork);
			eventList = intervals.getNetworkEventList(untypedNetwork);
			backwardIntegration();
			net = forwardSimulateNetwork();
		}

		this.setRootEdge(((NetworkNode) net[1]).getParentEdges().get(0));
	}

	private void backwardIntegration() {
		nodeStateProbabilities = new DoubleMatrix[untypedNetwork.getInternalNodes().size()];
		nrSamples = untypedNetwork.getLeafNodes().size();
		nodes = new ArrayList<>(untypedNetwork.getInternalNodes());

		int nIntervals = eventList.size();

		parents = new int[nIntervals];

		int MAX_SIZE = nIntervals * types;
		linProbs_tmp = new double[MAX_SIZE];
		linProbs = new double[MAX_SIZE];
		linProbsNew = new double[MAX_SIZE];

		euler = new Euler2ndOrder();
		euler.setup(MAX_SIZE, types, STEP_SIZE_BACKWARD_INTEGRATION, MAX_STEP_FOR_BACKWARD_INTEGRATION);

		activeLineages.clear();
		nrLineages = 0;
		linProbsLength = 0;
		int networkInterval = 0, ratesInterval = 0;
		double nextEventTime = 0.0;
		double prevEventTime = 0.0;

		StructuredNetworkEvent nextNetworkEvent = eventList.get(networkInterval);
		StructuredNetworkEvent startEvent = new StructuredNetworkEvent();
		double nextRateShift = dynamics.getInterval(ratesInterval);
		double nextNetworkEventTime = nextNetworkEvent.node.getHeight();

		setUpDynamics();

		// Get rates
		coalescentRates = dynamics.getCoalescentRate(ratesInterval);
		reassortmentRates = dynamics.getReassortmentRate(ratesInterval);

		nrLineages = activeLineages.size();
		linProbsLength = nrLineages * types;

		// Calculate the likelihood
		do {
			nextEventTime = Math.min(nextNetworkEventTime, nextRateShift);
			if (nextEventTime > 0) { // if true, calculate the interval contribution
				startEvent = eventList.get(networkInterval - 1);
				startEvent.activeLineages = new ArrayList<NetworkEdge>();
				startEvent.activeLineages.addAll(activeLineages);
				doEuler(prevEventTime, nextEventTime, ratesInterval, startEvent);
			}

			if (nextNetworkEventTime <= nextRateShift) {
				switch (nextNetworkEvent.type) {
				case COALESCENCE:
					nrLineages--;
					coalesce(nextNetworkEvent);
					break;

				case SAMPLE:
					nrLineages++;
					sample(nextNetworkEvent);
					break;

				case REASSORTMENT:
					reassortment(nextNetworkEvent);
					nrLineages++;
					break;
				}

				networkInterval++;
				nextRateShift -= nextNetworkEventTime;
				try {
					nextNetworkEvent = eventList.get(networkInterval);
					nextNetworkEventTime = nextNetworkEvent.node.getHeight();
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
		} while (nextNetworkEventTime <= Double.POSITIVE_INFINITY);
	}



	private Object[] forwardSimulateNetwork() {
		int nIntervals = eventList.size() - 1;
		StructuredNetworkEvent currentEvent = eventList.get(nIntervals);
		StructuredNetworkEvent nextEvent = eventList.get(nIntervals - 1);
		int startType = Randomizer.randomChoicePDF(nodeStateProbabilities[nodeStateProbabilities.length - 1].toArray());

		HashMap<NetworkEdge, Integer> lineageType = new HashMap<NetworkEdge, Integer>();
		for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
			lineageType.put(nextEvent.activeLineages.get(j), startType);
		}

		NetworkEdge rootEdge = currentEvent.node.getParentEdges().get(0);
		NetworkNode root = rootEdge.childNode;
		root.setTypeIndex(startType);
		root.setTypeLabel(dynamics.getStringStateValue(startType));
		root.setMetaData(currentEvent.node.getMetaData());

		NetworkNode currentNode = root;
		double currentTime = currentEvent.node.getHeight();
		double endTime = nextEvent.node.getHeight();

		double[] rates = new double[types];
		double[] ratesNext = new double[types];
		double totalRate;
		double totalRateNext;
		int its = 0;

		while (nIntervals > 0) {
			if (its > 1000)
				return new Object[] { false, null };
			while (true) {
				double dt = (endTime - currentTime) / FORWARD_INTEGRATION_STEPS;
				// If the length of interval between events is 0, jump to next event
				if (dt == 0.0)
					break;
				NetworkEdge minEdge = new NetworkEdge();
				double minTime = Double.NEGATIVE_INFINITY;
				double[] minRates = new double[types];

				for (NetworkEdge e : nextEvent.activeLineages)
				{
					double K = Math.log(Randomizer.nextDouble());
					double I = 0.0;
					double t = currentTime;
					double currentTime_l = currentTime;
					int idx = nextEvent.activeLineages.indexOf(e);

					totalRate = getTotalForwardsRate(lineageType.get(e), currentTime, idx,
							rates, nextEvent);

					int integrationStep;
					for (integrationStep = 0; integrationStep < FORWARD_INTEGRATION_STEPS; integrationStep++) {
						double tnext = currentTime_l
								- (currentTime_l - endTime) * (integrationStep + 1) / FORWARD_INTEGRATION_STEPS;
						totalRateNext = getTotalForwardsRate(lineageType.get(e), tnext, idx,
								ratesNext, nextEvent);

						I += 0.5 * (totalRateNext + totalRate) * dt;

						if (I <= K) {
							currentTime_l = t + 0.5 * dt;
							if (currentTime_l > minTime) {
								minTime = currentTime_l;
								minEdge = e;
								minRates = ratesNext.clone();
							}
							break;
						}
						totalRate = totalRateNext;
						t = tnext;
					}
				}

				if (minTime > endTime) {


					NetworkEdge lineage = minEdge;
					NetworkNode parent = lineage.parentNode;
					NetworkEdge newEdge = new NetworkEdge();
					NetworkNode newNode = new NetworkNode();
					newEdge.hasSegments = lineage.hasSegments;

					int oldType = lineageType.get(lineage);
					
					int newType = 0;
					if (minRates[0] == 0.0 && minRates[1] == 0.0) {
						break;
					}
					newType = Randomizer.randomChoicePDF(minRates);

					if (oldType == newType) {// happens only when rate of change is very close to zero
						its++;
						Network test = new Network();
						test.setRootEdge(root.getParentEdges().get(0));
						System.out.println(test.getExtendedNewick());
						System.out.println("Type matched fired");
						System.exit(0);
						continue;
					}

					newNode.setTypeIndex(oldType);
					newNode.setTypeLabel(dynamics.getStringStateValue(oldType));

					// track type change on a lineage
					lineageType.put(lineage, newType);

					parent.removeChildEdge(lineage);
					newNode.addChildEdge(lineage);
					newNode.addParentEdge(newEdge);
					parent.addChildEdge(newEdge);

					currentNode = newNode;
					newNode.setHeight(minTime);

					currentTime = minTime;

				}
				else
					break;
			}

			switch (nextEvent.type) {
			case SAMPLE:
				NetworkEdge sampleLineage = nextEvent.node.getParentEdges().get(0);
				NetworkNode sampleNode = new NetworkNode();
				sampleNode.setMetaData(nextEvent.node.getMetaData());

				sampleNode.setHeight(endTime);
				sampleNode.setTaxonLabel(nextEvent.node.getTaxonLabel());
				sampleNode.setTaxonIndex(nextEvent.node.getTaxonIndex());
				String label = nextEvent.node.getTaxonLabel();
				int labelIdx = dynamics.getValue(label);
				int idx = lineageType.get(sampleLineage);
				if (idx != labelIdx) {
					return new Object[] { false, null };
				}

				sampleNode.setTypeIndex(idx);
				sampleNode.setTypeLabel(dynamics.getStringStateValue(idx));

				nextEvent.node.removeParentEdge(sampleLineage);
				sampleNode.addParentEdge(sampleLineage);

				lineageType.remove(sampleLineage);

				break;
			case COALESCENCE:
				NetworkEdge coalLineage = nextEvent.node.getParentEdges().get(0);
				int coalType = lineageType.get(coalLineage);
				NetworkNode coalNode = new NetworkNode();
				coalNode.setMetaData(nextEvent.node.getMetaData());

				coalNode.setHeight(endTime);
				coalNode.setTypeIndex(coalType);
				coalNode.setTypeLabel(dynamics.getStringStateValue(coalType));

				nextEvent.node.removeParentEdge(coalLineage);
				NetworkEdge child1 = nextEvent.node.getChildEdges().get(0);
				NetworkEdge child2 = nextEvent.node.getChildEdges().get(1);

				nextEvent.node.removeChildEdge(child1);
				nextEvent.node.removeChildEdge(child2);

				coalNode.addParentEdge(coalLineage);
				coalNode.addChildEdge(child1);
				coalNode.addChildEdge(child2);

				lineageType.remove(coalLineage);
				lineageType.put(child1, coalType);
				lineageType.put(child2, coalType);

				break;
			case REASSORTMENT:
				NetworkEdge parent1 = nextEvent.node.getParentEdges().get(0);
				NetworkEdge parent2 = nextEvent.node.getParentEdges().get(1);
				int type1 = lineageType.get(parent1);
				int type2 = lineageType.get(parent2);
				int reassType = Integer.MAX_VALUE;
				
				if (type1==type2)
					reassType = type1;
				else {
					if (rejection)
					{ // reject simulation if reassortment parents have different types
						return new Object[] { false, null };
					}
					else {
							double[] r1 = new double[types];
							double rates1 = getTotalForwardsRate(lineageType.get(parent1),
									parent1.childNode.getHeight(), nextEvent.activeLineages.indexOf(parent1),
									r1, nextEvent);
							double[] r2 = new double[types];
							double rates2 = getTotalForwardsRate(lineageType.get(parent2),
									parent1.childNode.getHeight(), nextEvent.activeLineages.indexOf(parent2),
									r2, nextEvent);

						if (inheritReaType) { // Reassortment child lineage must inherit the type of at least one parent
							while (reassType != type1 && reassType != type2) {
								reassType = rates1 > rates2 ? Randomizer.randomChoicePDF(r1)
										: Randomizer.randomChoicePDF(r2);
							}

						} else // Reassortment child lineage may or may not inherit the type of at least one
								// parent
							reassType = rates1 > rates2 ? Randomizer.randomChoicePDF(r1)
									: Randomizer.randomChoicePDF(r2);
					}
				}

				// Need to chose the type for child from cumulative rates of the parents
				// Or migrate one child to the others type based on rates
				NetworkNode reassNode = new NetworkNode();

				reassNode.setHeight(endTime);
				reassNode.setTypeIndex(reassType);
				reassNode.setTypeLabel(dynamics.getStringStateValue(reassType));
				
				nextEvent.node.removeParentEdge(parent1);
				nextEvent.node.removeParentEdge(parent2);
				NetworkEdge child = nextEvent.node.getChildEdges().get(0);
				nextEvent.node.removeChildEdge(child);

				if (reassType != type2) {
						NetworkNode fakeMigNode2 = new NetworkNode();
						NetworkEdge fakeMigEdge2 = new NetworkEdge();

						fakeMigNode2.setHeight(endTime);
						fakeMigNode2.setTypeIndex(type2);
						fakeMigNode2.setTypeLabel(dynamics.getStringStateValue(type2));

						fakeMigNode2.addParentEdge(parent2);
						fakeMigEdge2.hasSegments = parent2.hasSegments;
						fakeMigNode2.addChildEdge(fakeMigEdge2);

						reassNode.addParentEdge(fakeMigEdge2);
					}
					if (reassType != type1) {
						NetworkNode fakeMigNode1 = new NetworkNode();
						NetworkEdge fakeMigEdge1 = new NetworkEdge();

					fakeMigNode1.setHeight(endTime);
					fakeMigNode1.setTypeIndex(type1);
					fakeMigNode1.setTypeLabel(dynamics.getStringStateValue(type1));

						fakeMigNode1.addParentEdge(parent1);
						fakeMigEdge1.hasSegments = parent1.hasSegments;
						fakeMigNode1.addChildEdge(fakeMigEdge1);

						reassNode.addParentEdge(fakeMigEdge1);

					}
				if (reassType == type1)
					reassNode.addParentEdge(parent1);
				if (reassType == type2)
					reassNode.addParentEdge(parent2);

				reassNode.addChildEdge(child);

				lineageType.remove(parent1);
				lineageType.remove(parent2);
				lineageType.put(child, reassType);

				break;
			default:
				throw new IllegalArgumentException("Switch fell through in forward simulation.");
			}

			nIntervals -= 1;
			its = 0;
			if (nIntervals > 0) {
				currentEvent = eventList.get(nIntervals);
				nextEvent = eventList.get(nIntervals - 1);
				currentTime = currentEvent.node.getHeight();
				endTime = nextEvent.node.getHeight();
			}
		}

		return new Object[] { true, root };
	}


	private double getTotalForwardsRate(int fromType, double t, int lineageIdx,
			double[] rates, StructuredNetworkEvent nextEvent) {
		double totalRate = 0.0;
		getForwardsRates(fromType, t, lineageIdx, rates, nextEvent);
		for (int type = 0; type < types; type++)
			totalRate += rates[type];

		return totalRate;
	}

	private double[] getForwardsRates(int fromType, double time, int lineageIdx,
			double[] result, StructuredNetworkEvent nextEvent) {

		// Assert if x2 is not 0
		int x1 = -1;
		double time1 = 0;
		double time2 = 0;
		boolean interpolate = true;
		int x2 = bigger(nextEvent.intermediateTimeStored, time);
		if (x2 > nextEvent.intermediateTimeStored.length - 1) {
			x2 = nextEvent.intermediateTimeStored.length - 1;
			interpolate = false;
		}
		else if (x2 == 0) {
			interpolate = false;
		}

		if (interpolate && Math.abs(nextEvent.intermediateTimeStored[x2] - time) > 1e-16) {
			x1 += x2;

			time1 = nextEvent.intermediateTimeStored[x1];
			time2 = nextEvent.intermediateTimeStored[x2];
		}
		else
			interpolate = false;


		int ratesInterval = getIntervalIndex(time);
		double[] migMatrix = dynamics.getBackwardsMigration(ratesInterval);
		int n = (int) (Math.sqrt(migMatrix.length) + 0.5);

		for (int type = 0; type < types; type++) {
			if (type == fromType) {
				result[type] = 0.0;
				continue;
			}

//			UnivariateInterpolator interpolator2;
			double pTo;
//			double[] points = new double[nextEvent.p_stored.length];
			if (interpolate) {
//				interpolator2 = new SplineInterpolator();
//				for (int ss = 0; ss < nextEvent.p_stored.length; ss++) {
//					interpolator.addSamplePoint(nextEvent.intermediateTimeStored[ss],
//							nextEvent.p_stored[ss]);
//					points[ss] = nextEvent.p_stored[ss][lineageIdx * types + type];
//
////							nextEvent.pDot_stored[ss]);
//
//				}
//				UnivariateFunction function = interpolator2.interpolate(nextEvent.intermediateTimeStored, points);
//				System.out.println(function.value(time));

				pTo = (nextEvent.p_stored[x1][lineageIdx * types + type] * (time - time1)
						+ nextEvent.p_stored[x2][lineageIdx * types + type] * (time2 - time))
						/ (time2 - time1);
//				System.out.println(pTo);
			}
			else {
				if (x2 == -1) {
					System.out.println("Something wrong in mapping: getForwardsRates");
					System.exit(1);
				}
				pTo = nextEvent.p_stored[x2][lineageIdx * types + type];
			}

			result[type] = migMatrix[type * n + fromType] * pTo; // p[lineageIdx * score.types + type];
		}

		double pFrom;
		if (interpolate)
			pFrom = (nextEvent.p_stored[x1][lineageIdx * types + fromType] * (time - time1)
					+ nextEvent.p_stored[x2][lineageIdx * types + fromType] * (time2 - time))
					/ (time2 - time1);
		else
			pFrom = nextEvent.p_stored[x2][lineageIdx * types + fromType];

		if (pFrom <= 0.0) {
			// The source type prob approaches zero as the integration closes
			// in on a node with a defined type. This causes the transition
			// rate to this type to become infinite. What follows is a hack
			// to ensure that this important situation is handled properly.

			int maxRateIdx = -1;
			double maxRate = Double.NEGATIVE_INFINITY;
			for (int type = 0; type < types; type++) {
				if (result[type] > maxRate) {
					maxRateIdx = type;
					maxRate = result[type];
				}
			}

			for (int type = 0; type < types; type++)
				result[type] = type == maxRateIdx ? 1.0 : 0.0;

		} else {
			// Apply source type prob as rate denominator:

			for (int type = 0; type < types; type++)
				result[type] /= pFrom;

		}

		return result;
	}


	/**
	 * Finds the index of the rates shift interval t lies in. Note that if t lies on
	 * a boundary between intervals, the interval returned will be the _earlier_ of
	 * these two intervals.
	 *
	 * @param t time for which to identify interval
	 * @return index identifying interval.
	 */
	public int getIntervalIndex(double t) {

		int index = Arrays.binarySearch(dynamics.getIntervals(), t);

		if (index < 0)
			index = -index - 1;

		// return at most the index of the last interval (m-1)
		return Math.max(0, Math.min(index, dynamics.getIntervals().length - 1));
	}

	// XXX Backwards stuff
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

			if (sampleState >= types) {
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
					linProbsNew[linCount * types + j] = linProbs[i * types + j];
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

	private double doEuler(double start, double end, int ratesInterval, StructuredNetworkEvent startEvent) {
		double duration = end - start;

		startEvent.numRecords = nRecordsInput.get();
		startEvent.p_stored = new double[startEvent.numRecords][linProbsLength];
		startEvent.pDot_stored = new double[startEvent.numRecords][linProbsLength];
		startEvent.intermediateTimeStored = new double[startEvent.numRecords];
		Arrays.fill(startEvent.intermediateTimeStored, end);

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

		euler.initAndcalculateValues(ratesInterval, nrLineages, duration, linProbs_tmp, linProbsLength + 1, n_segs,
				startEvent);

		System.arraycopy(linProbs_tmp, 0, linProbs, 0, linProbsLength);

		return linProbs_tmp[linProbsLength];
	}

	// XXX logging

	private long lastRemapSample = -1;

	/**
	 * Remap the tree. Intended to be called by loggers requiring a mapped tree.
	 * Supplying the sample number allows the result of the remapping to be cached
	 * and used for other loggers.
	 *
	 * @param sample sample number at log
	 */
	public void remapForLog(long sample) {
		if (!remapOnLogInput.get() || sample == lastRemapSample)
			return;

		doStochasticMapping();
		lastRemapSample = sample;
	}

	@Override
	public void init(PrintStream out) {
		untypedNetwork.init(out);
	}

	@Override
	public void log(long sample, PrintStream out) {
		remapForLog(sample);

		out.print("tree STATE_" + sample + " = ");
		// Don't sort, this can confuse CalculationNodes relying on the tree
		// tree.getRoot().sort();
		final String newick = this.getExtendedNewick();
		out.print(newick);
	}

	@Override
	public void close(PrintStream out) {
		untypedNetwork.close(out);
	}

	// XXX utils

	public static int bigger(double[] arr, double target) {
		int idx = Arrays.binarySearch(arr, target);
		if (idx < 0) {
			// target not found, so return index of insertion point
			return -idx - 1;
		}
		// target found, so return that
		return idx;
	}
}
