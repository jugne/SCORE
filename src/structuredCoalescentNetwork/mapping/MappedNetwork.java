package structuredCoalescentNetwork.mapping;


import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import beast.core.Input;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import structuredCoalescentNetwork.distribution.SCORE;
import structuredCoalescentNetwork.distribution.StructuredNetworkEvent;


/**
 * @author Ugne Stolz. Created on 15 Nov 2019. Mostly Uses ideas of stochastic
 *         mapping implementation from Tim Vaughan's BDMM Prime implementation
 */
public class MappedNetwork extends Network {

	public Input<SCORE> scoreDistribInput = new Input<>("scoreDistribution",
			"If provided, extract the parameterization from here.", Input.Validate.REQUIRED);

	public Input<Boolean> mapOnInitInput = new Input<>("mapOnInit",
			"If true, mapping will be performed when object is " + "first initialize.", true);

	private SCORE score;

	private Network network;
	private List<StructuredNetworkEvent> eventList;

	/**
	 * Maximum number of steps in each waiting time calculation in forward
	 * simulation.
	 */
	private final int FORWARD_INTEGRATION_STEPS = 100;

	@Override
	public void initAndValidate() {
		score = scoreDistribInput.get();
		network = score.network;


		eventList = score.networkEventList;

//		if (mapOnInitInput.get())
//			doStochasticMapping();
	}

	public void doStochasticMapping() {

//		if (mapOnInitInput.get())
//			int rootNodeType = Randomizer.randomChoicePDF(score.getRootTypes().toArray());

		network = score.network;
		eventList = score.networkEventList;


//		this.storedNet = (Network) network.copy();

		NetworkNode net = forwardSimulateNetwork();
		Network netw = new Network();
		netw.setRootEdge(net.getParentEdges().get(0));
//		System.out.println(netw.getExtendedNewick());


	}

	private NetworkNode forwardSimulateNetwork() {
		int nIntervals = eventList.size() - 1;
		StructuredNetworkEvent currentEvent = eventList.get(nIntervals);
		StructuredNetworkEvent nextEvent = eventList.get(nIntervals - 1);
		int startType = Randomizer.randomChoicePDF(score.getRootTypes().toArray());

		HashMap<NetworkEdge, Integer> lineageType = new HashMap<NetworkEdge, Integer>();
		for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
			lineageType.put(nextEvent.activeLineages.get(j), startType);
		}

		NetworkEdge rootEdge = currentEvent.node.getParentEdges().get(0);
		NetworkNode root = rootEdge.childNode;
		root.setTypeIndex(startType);
		root.setTypeLabel(score.dynamics.getStringStateValue(startType));
		root.setMetaData(currentEvent.node.getMetaData());
		
		NetworkNode currentNode = root;
		double currentTime = currentEvent.time;
		double endTime = nextEvent.time;
		
		double[] rates = new double[score.types];
		double[] ratesNext = new double[score.types];
		double totalRate;
		double totalRateNext;
		
		while (nIntervals > 0) {
			// construct interpolator for each interval
			// here we are going forwards in time, therefore our lasts stored value is now
			// our first
			// this doesn't matter for interpolator construction,
			// but is written down in such way for clarity
//			HermiteInterpolator interpolator = new HermiteInterpolator();
//			for (int i = nextEvent.intermediateTimeStored.length - 1; i >= 0; i--) {
////				double[] p_sqrt = root(Arrays.copyOfRange(nextEvent.p_stored[i], 0, nextEvent.p_stored[i].length - 1));
////				double[] p_sqrt_Dot = pDot(p_sqrt,
////						Arrays.copyOfRange(nextEvent.pDot_stored[i], 0, nextEvent.pDot_stored[i].length - 1));
//				interpolator.addSamplePoint(nextEvent.intermediateTimeStored[i],
//
////						
//						Arrays.copyOfRange(root(nextEvent.p_stored[i]), 0, nextEvent.p_stored[i].length - 1));
////						Arrays.copyOfRange(nextEvent.pDot_stored[i], 0, nextEvent.pDot_stored[i].length - 1),
////						Arrays.copyOfRange(nextEvent.pDotDot_stored[i], 0, nextEvent.pDotDot_stored[i].length - 1));
//																												// nextEvent.pDotDot_stored[i]);
//			}
			
//			LinearInterpolator linearInterpolator = new LinearInterpolator();
//			HashMap<Integer, HashMap<Integer, PolynomialSplineFunction>> lineagePolynomials = new HashMap<Integer, HashMap<Integer,PolynomialSplineFunction>>();
//			HashMap<Integer, PolynomialSplineFunction> typePolynomials = new HashMap<Integer, PolynomialSplineFunction>();
//			for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
//				for (int i=0; i<score.types; i++) {
//					typePolynomials.put(i, linearInterpolator.interpolate(nextEvent.intermediateTimeStored, arg1))
//					
//					
//				}
//				
//			}
			


			while (true) {

				double dt = (endTime - currentTime) / FORWARD_INTEGRATION_STEPS;
				int minLineageIdx = 0;
				double minTime = Double.NEGATIVE_INFINITY;

				for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
					double K = Math.log(Randomizer.nextDouble());
					double I = 0.0;
					double t = currentTime;
					double currentTime_l = currentTime;

					totalRate = getTotalForwardsRate(lineageType.get(nextEvent.activeLineages.get(j)), currentTime, j,
							rates, nextEvent);

					int integrationStep;
					for (integrationStep = 0; integrationStep < FORWARD_INTEGRATION_STEPS; integrationStep++) {
						// TODO debug if time behaves correctly
						double tnext = currentTime_l
								- (currentTime_l - endTime) * (integrationStep + 1) / FORWARD_INTEGRATION_STEPS;
						totalRateNext = getTotalForwardsRate(lineageType.get(nextEvent.activeLineages.get(j)), tnext, j,
								ratesNext, nextEvent);
//						System.out.println(totalRateNext);

						I += 0.5 * (totalRateNext + totalRate) * dt;

						if (I <= K) {
							currentTime_l = t + 0.5 * dt;
							if (currentTime_l > minTime) {
								minTime = currentTime_l;
								minLineageIdx = j;
							}
							break;
						}
						totalRate = totalRateNext;
						double[] tmp = rates;
						rates = ratesNext;
						ratesNext = tmp;

						t = tnext;
					}


				}

				if (minTime > endTime) {
//					currentNode.setHeight(currentEvent.time - currentTime);
					Network test = new Network();
					test.setRootEdge(root.getParentEdges().get(0));
					System.out.println(network.getExtendedNewick());
					System.out.println("test: " + test.getExtendedNewick());

					NetworkEdge lineage = nextEvent.activeLineages.get(minLineageIdx);
					NetworkNode parent = lineage.parentNode;
					NetworkEdge newEdge = new NetworkEdge();
					NetworkNode newNode = new NetworkNode();
					newEdge.hasSegments = lineage.hasSegments;

					int oldType = lineageType.get(lineage);
					int newType = Randomizer.randomChoicePDF(ratesNext);
					newNode.setTypeIndex(oldType);
					newNode.setTypeLabel(score.dynamics.getStringStateValue(oldType));
					// track type change on a lineage
//					nextEvent.node.setTypeIndex(newType);
//					nextEvent.node.setTypeLabel(score.dynamics.getStringStateValue(newType));
//					lineage.childNode.setTypeIndex(newType);
//					lineage.childNode.setTypeLabel(score.dynamics.getStringStateValue(newType));

					lineageType.put(lineage, newType);

					parent.addChildEdge(newEdge);
					parent.removeChildEdge(lineage);

					newNode.addChildEdge(lineage);
					newNode.addParentEdge(newEdge);

					currentNode = newNode;
//					newNode.setHeight(currentEvent.time - currentTime);
					newNode.setHeight(minTime);

//					Network test = new Network();
					test.setRootEdge(root.getParentEdges().get(0));

					currentTime = minTime;

//					System.out.println(network.getExtendedNewick());
					System.out.println(test.getExtendedNewick());
				}
				else
					break;
			}
			// how to deal with events???
//			currentNode.setHeight(root.getHeight() - endTime);
			switch (nextEvent.type) {
			case SAMPLE:
				NetworkEdge sampleLineage = nextEvent.node.getParentEdges().get(0);

				NetworkNode sampleNode = new NetworkNode();
				// TODO debug this!
				sampleNode.setMetaData(nextEvent.node.getMetaData());

				sampleNode.setHeight(endTime);
				sampleNode.setTaxonLabel(nextEvent.node.getTaxonLabel());
				sampleNode.setTaxonIndex(nextEvent.node.getTaxonIndex());
				String label = nextEvent.node.getTaxonLabel();
				sampleNode.setTypeIndex(score.dynamics.getValue(label));
				sampleNode.setTypeLabel(score.dynamics.getStringStateValue(sampleNode.getTypeIndex()));

				nextEvent.node.removeParentEdge(sampleLineage);
				sampleNode.addParentEdge(sampleLineage);

				lineageType.remove(sampleLineage);
				break;
			case COALESCENCE:
				Network test = new Network();
				test.setRootEdge(root.getParentEdges().get(0));
				System.out.println(network.getExtendedNewick());
				System.out.println(test.getExtendedNewick());

				NetworkEdge coalLineage = nextEvent.node.getParentEdges().get(0);
				int coalType = lineageType.get(coalLineage);
				NetworkNode coalNode = new NetworkNode();
				// TODO debug this!

				coalNode.setMetaData(nextEvent.node.getMetaData());

				coalNode.setHeight(endTime);
				coalNode.setTypeIndex(coalType);
				coalNode.setTypeLabel(score.dynamics.getStringStateValue(coalType));

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

				test = new Network();
				test.setRootEdge(root.getParentEdges().get(0));

//				System.out.println(network.getExtendedNewick());
				System.out.println(test.getExtendedNewick());
				break;
			case REASSORTMENT:
				NetworkEdge parent1 = nextEvent.node.getParentEdges().get(0);
				NetworkEdge parent2 = nextEvent.node.getParentEdges().get(1);
				int type1 = lineageType.get(parent1);
				int type2 = lineageType.get(parent2);
				int reassType = 0;
				
				if (type1==type2)
					reassType = type1;
				else {
					test = new Network();
					test.setRootEdge(root.getParentEdges().get(0));
					System.out.println(test.getExtendedNewick());
					System.out.println("What now?");
					System.exit(0);
//					double[] r1 = new double[score.types];
//					double[] r2 = new double[score.types];
//					getForwardsRates(type1, endTime, nextEvent.activeLineages.indexOf(parents.get(0)), interpolator,
//							r1);
				}
				NetworkNode reassNode = new NetworkNode();
				reassNode.setHeight(endTime);
				reassNode.setTypeIndex(reassType);
				reassNode.setTypeLabel(score.dynamics.getStringStateValue(reassType));
				
				nextEvent.node.removeParentEdge(parent1);
				nextEvent.node.removeParentEdge(parent2);
				NetworkEdge child = nextEvent.node.getChildEdges().get(0);
				nextEvent.node.removeChildEdge(child);

				reassNode.addParentEdge(parent1);
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
			Network test = new Network();
			test.setRootEdge(root.getParentEdges().get(0));

//			System.out.println(network.getExtendedNewick());
			System.out.println(test.getExtendedNewick());
			if (nIntervals > 0) {
				currentEvent = eventList.get(nIntervals);
				nextEvent = eventList.get(nIntervals - 1);
				currentTime = currentEvent.time;
				endTime = nextEvent.time;

//				for (int j = 0; j < nextEvent.lineagesRemoved.size(); j++) {
//					lineageType.put(nextEvent.lineagesRemoved.get(j),
//							nextEvent.activeLineages.get(j).parentNode.getTypeIndex());
//				}
			}


		}
//		score.restore();
		return root;

	}


	private double getTotalForwardsRate(int fromType, double t, int lineageIdx,
			double[] rates, StructuredNetworkEvent nextEvent) {
		double totalRate = 0.0;
		getForwardsRates(fromType, t, lineageIdx, rates, nextEvent);
		for (int type = 0; type < score.types; type++)
			totalRate += rates[type];

		return totalRate;
	}

	private double[] getForwardsRates(int fromType, double time, int lineageIdx,
			double[] result, StructuredNetworkEvent nextEvent) {

		// Assert if x2 is not 0
		int x1 = -1;
		double time1 = 0;
		double time2 = 0;
		boolean interpolate = false;
		int x2 = bigger(nextEvent.intermediateTimeStored, time);
		if (x2 == 0)
			System.out.println(x2);
		if (Math.abs(nextEvent.intermediateTimeStored[x2] - time) > 1e-10) {
			interpolate = true;
			x1 += x2;

			time1 = nextEvent.intermediateTimeStored[x1];
			time2 = nextEvent.intermediateTimeStored[x2];
		}



//		double[] p = interpolator.value(time);
		int interval = getIntervalIndex(time);
		double[] migMatrix = score.dynamics.getBackwardsMigration(interval);
		int n = (int) (Math.sqrt(migMatrix.length) + 0.5);

		for (int type = 0; type < score.types; type++) {
			if (type == fromType) {
				result[type] = 0.0;
				continue;
			}

			double pTo;
			if (interpolate)
				pTo = (nextEvent.p_stored[x1][lineageIdx * score.types + type] * (time - time1)
					+ nextEvent.p_stored[x2][lineageIdx * score.types + type] * (time2 - time)) / (time2 - time1);
			else
				pTo = nextEvent.p_stored[x2][lineageIdx * score.types + type];

			result[type] = migMatrix[type * n + fromType] * pTo; // p[lineageIdx * score.types + type];
		}

		double pFrom;
		if (interpolate)
			pFrom = (nextEvent.p_stored[x1][lineageIdx * score.types + fromType] * (time - time1)
				+ nextEvent.p_stored[x2][lineageIdx * score.types + fromType] * (time2 - time)) / (time2 - time1);
		else
			pFrom = nextEvent.p_stored[x2][lineageIdx * score.types + fromType];

		if (pFrom <= 0.0) {
			// The source type prob approaches zero as the integration closes
			// in on a node with a defined type. This causes the transition
			// rate to this type to become infinite. What follows is a hack
			// to ensure that this important situation is handled properly.

			int maxRateIdx = -1;
			double maxRate = Double.NEGATIVE_INFINITY;
			for (int type = 0; type < score.types; type++) {
				if (result[type] > maxRate) {
					maxRateIdx = type;
					maxRate = result[type];
				}
			}

			for (int type = 0; type < score.types; type++)
				result[type] = type == maxRateIdx ? 1.0 : 0.0;

		} else {
			// Apply source type prob as rate denominator:

			for (int type = 0; type < score.types; type++)
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

		int index = Arrays.binarySearch(score.dynamics.getIntervals(), t);

		if (index < 0)
			index = -index - 1;

		// return at most the index of the last interval (m-1)
		return Math.max(0, Math.min(index, score.dynamics.getIntervals().length - 1));
	}

	/**
	 * Return time of node, i.e. N - node_age.
	 *
	 * @param node node whose time to query.
	 * @return time of node.
	 */
	public double getNodeTime(NetworkNode node) {
		return network.getRootEdge().childNode.getHeight() - node.getHeight();
	}

	private double[] root(double[] array) {

		double[] squared = new double[array.length];

		for (int i = 0; i < array.length; i++)
			squared[i] = Math.sqrt(array[i]);

		return squared;
	}

	private static double[] pDot(final double[] p_sqrt, final double[] pDot) {
		double[] p_sqrt_dot = new double[p_sqrt.length];
		for (int i = 0; i < p_sqrt.length; i++)
			p_sqrt_dot[i] = pDot[i] / (2 * p_sqrt[i]);

		return p_sqrt_dot;
	}

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
