package structuredCoalescentNetwork.mapping;


import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.analysis.interpolation.HermiteInterpolator;

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
	private double[] rootNodeType;
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

		if (mapOnInitInput.get())
			doStochasticMapping();
	}

	private void doStochasticMapping() {

		int rootNodeType = Randomizer.randomChoicePDF(score.getRootTypes().toArray());

		forwardSimulateNetwork();

	}

	private void forwardSimulateNetwork() {
		int nEvents = eventList.size();
		int nIntervals = eventList.size() - 1;
		StructuredNetworkEvent currentEvent = eventList.get(nIntervals);
		StructuredNetworkEvent nextEvent = eventList.get(nIntervals - 1);
		int startType = currentEvent.node.getTypeIndex();

		HashMap<NetworkEdge, Integer> lineageType = new HashMap<NetworkEdge, Integer>();
		for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
			lineageType.put(nextEvent.activeLineages.get(j), startType);
		}

		NetworkNode root = currentEvent.node;
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
			HermiteInterpolator interpolator = new HermiteInterpolator();
			for (int i = 1; i <= nextEvent.intermediateTimeStored.length; i++) {
				interpolator.addSamplePoint(nextEvent.intermediateTimeStored[-i], nextEvent.p_stored[-i],
						nextEvent.pDot_stored[-i], nextEvent.pDotDot_stored[-i]);
			}

			while (true) {

				double dt = (endTime - currentTime) / FORWARD_INTEGRATION_STEPS;
				int minLineageIdx = 0;
				double minTime = Double.POSITIVE_INFINITY;

				for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
					double K = -Math.log(Randomizer.nextDouble());
					double I = 0.0;
					double t = currentTime;

					totalRate = getTotalForwardsRate(lineageType.get(nextEvent.activeLineages.get(j)), currentTime, j,
							interpolator, rates);

					int integrationStep;
					for (integrationStep = 0; integrationStep < FORWARD_INTEGRATION_STEPS; integrationStep++) {
						// TODO debug if time behaves correctly
						double tnext = currentTime
								- (currentTime - endTime) * (integrationStep + 1) / FORWARD_INTEGRATION_STEPS;
						totalRateNext = getTotalForwardsRate(lineageType.get(nextEvent.activeLineages.get(j)), tnext, j,
								interpolator, ratesNext);

						I += 0.5 * (totalRateNext + totalRate) * dt;

						if (I >= K) {
							currentTime = t + 0.5 * dt;
							if (currentTime < minTime) {
								minTime = currentTime;
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
//					// If we here because we ran through max number of steps
//					// continue with the next lineage
//					if (integrationStep == FORWARD_INTEGRATION_STEPS)
//						continue;

				}

				if (minTime < endTime) {
//					currentNode.setHeight(currentEvent.time - currentTime);

					NetworkEdge lineage = nextEvent.activeLineages.get(minLineageIdx);
					NetworkEdge newEdge = new NetworkEdge();
					NetworkNode newNode = new NetworkNode();
					int newType = Randomizer.randomChoicePDF(ratesNext);
					newNode.setTypeIndex(newType);
					newNode.setTypeLabel(score.dynamics.getStringStateValue(newType));
					// track type change on a lineage
					lineageType.put(lineage, newType);

					currentNode.addChildEdge(newEdge);
					currentNode.removeChildEdge(lineage);

					newNode.addChildEdge(lineage);
					newNode.addParentEdge(newEdge);

					currentNode = newNode;
					newNode.setHeight(currentEvent.time - currentTime);
				}
				else
					break;
			}
			// how to deal with events???
//			currentNode.setHeight(root.getHeight() - endTime);
			switch (nextEvent.type) {
			case SAMPLE:
				NetworkEdge sampleLineage = nextEvent.node.getParentEdges().get(0);

				// for debug only
				int sampleType = lineageType.get(sampleLineage);

				NetworkNode sampleNode = new NetworkNode();
				// TODO debug this!
				sampleNode.setMetaData(nextEvent.node.getMetaData());

				sampleNode.setHeight(endTime);
				sampleNode.setTaxonLabel(nextEvent.node.getTaxonLabel());
				sampleNode.setTaxonIndex(nextEvent.node.getTaxonIndex());

				nextEvent.node.removeParentEdge(sampleLineage);
				sampleNode.addParentEdge(sampleLineage);
			case COALESCENCE:
				NetworkEdge coalLineage = nextEvent.node.getParentEdges().get(0);
				int coalType = lineageType.get(coalLineage);
				NetworkNode coalNode = new NetworkNode();
				// TODO debug this!

				coalNode.setMetaData(nextEvent.node.getMetaData());

				coalNode.setHeight(endTime);
				coalNode.setTaxonIndex(coalType);
				coalNode.setTaxonLabel(score.dynamics.getStringStateValue(coalType));

				nextEvent.node.removeParentEdge(coalLineage);
				List<NetworkEdge> children = nextEvent.node.getChildEdges();
				nextEvent.node.removeChildEdge(children.get(0));
				nextEvent.node.removeChildEdge(children.get(1));

				coalNode.addParentEdge(coalLineage);
				coalNode.addChildEdge(children.get(0));
				coalNode.addChildEdge(children.get(1));
			case REASSORTMENT:
				List<NetworkEdge> parents = nextEvent.node.getParentEdges();
				int type1 = lineageType.get(parents.get(0));
				int type2 = lineageType.get(parents.get(1));
				int reassType = 0;
				
				if (type1==type2)
					reassType = type1;
				else {
					System.out.println("What now?");
					System.exit(0);
//					double[] r1 = new double[score.types];
//					double[] r2 = new double[score.types];
//					getForwardsRates(type1, endTime, nextEvent.activeLineages.indexOf(parents.get(0)), interpolator,
//							r1);
				}
				NetworkNode reassNode = new NetworkNode();
				reassNode.setHeight(endTime);
				reassNode.setTaxonIndex(reassType);
				reassNode.setTaxonLabel(score.dynamics.getStringStateValue(reassType));
				
				nextEvent.node.removeParentEdge(parents.get(0));
				nextEvent.node.removeParentEdge(parents.get(1));
				NetworkEdge child = nextEvent.node.getChildEdges().get(0);
				nextEvent.node.removeChildEdge(child);

				reassNode.addParentEdge(parents.get(0));
				reassNode.addParentEdge(parents.get(1));
				reassNode.addChildEdge(child);

			}

			nIntervals -= 1;
			currentEvent = eventList.get(nIntervals);
			nextEvent = eventList.get(nIntervals - 1);

			for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
				lineageType.put(nextEvent.activeLineages.get(j),
						nextEvent.activeLineages.get(j).parentNode.getTypeIndex());
			}

		}

	}

	private NetworkNode forwardSimulateSubtree(NetworkNode subNetRoot, double startTime, int startType) {
		int networkInterval = eventList.size();
		StructuredNetworkEvent currentEvent = eventList.get(networkInterval-1);
		StructuredNetworkEvent nextEvent = eventList.get(networkInterval-2);
		
		NetworkNode root = currentEvent.node;
		root.setTypeIndex(startType);
		root.setTypeLabel(score.dynamics.getStringStateValue(startType));

		NetworkNode currentNode = root;
//		currentNode.setTypeIndex(startType);
//		currentNode.setTypeLabel(score.dynamics.getStringStateValue(startType));
		double currentTime = currentEvent.time;

		double endTime = getNodeTime(nextEvent.node);
		double[] rates = new double[score.types];
		double[] ratesNext = new double[score.types];
		double totalRate;
		double totalRateNext;

		// Start going by the interval
		while (networkInterval >= 0) {
			// Make an interpolator for lineage type probabilities in current interval
			HermiteInterpolator interpolator = new HermiteInterpolator();
			for (int i = 1; i <= nextEvent.intermediateTimeStored.length; i++) {
				interpolator.addSamplePoint(nextEvent.intermediateTimeStored[-i], nextEvent.p_stored[-i],
						nextEvent.pDot_stored[-i], nextEvent.pDotDot_stored[-i]);
			}

			double dt = (endTime - currentTime) / FORWARD_INTEGRATION_STEPS;
			int minLineageIdx = 0;
			double minTime = Double.POSITIVE_INFINITY;

			while (true) {

			for (int j = 0; j < nextEvent.activeLineages.size(); j++) {
				double K = -Math.log(Randomizer.nextDouble());
				double I = 0.0;
				double t = currentTime;

				
					totalRate = getTotalForwardsRate(currentNode.getTypeIndex(), currentTime, j, interpolator,
						rates);
				int integrationStep;
				for (integrationStep = 0; integrationStep < FORWARD_INTEGRATION_STEPS; integrationStep++) {
					double tnext = currentTime
							+ (endTime - currentTime) * (integrationStep + 1) / FORWARD_INTEGRATION_STEPS;
						totalRateNext = getTotalForwardsRate(currentNode.getTypeIndex(), tnext, j, interpolator,
							ratesNext);

					I += 0.5 * (totalRateNext + totalRate) * dt;

					if (I >= K) {
						currentTime = t + 0.5*dt;
						if (currentTime < minTime) {
							minTime = currentTime;
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
					if (integrationStep == FORWARD_INTEGRATION_STEPS)
						continue;

				}

			if (minTime < endTime) {
				currentNode.setHeight(currentEvent.time - currentTime);

				NetworkEdge lineage = nextEvent.activeLineages.get(minLineageIdx);
				NetworkEdge newEdge = new NetworkEdge();
				NetworkNode newNode = new NetworkNode();
				int newType = Randomizer.randomChoicePDF(ratesNext);
				newNode.setTypeIndex(newType);
				newNode.setTypeLabel(score.dynamics.getStringStateValue(newType));

				currentNode.addChildEdge(newEdge);
				currentNode.removeChildEdge(lineage);

				newNode.addChildEdge(lineage);
				newNode.addParentEdge(newEdge);

				currentNode = newNode;
			}
				else {
					break;
				}
		}
			currentNode.setHeight(root.getHeight() - endTime);
			switch (nextEvent.type) {
			case SAMPLE:
				currentNode.setTaxonLabel(nextEvent.node.getTaxonLabel());
				currentNode.setTaxonIndex(nextEvent.node.getTaxonIndex());
				currentNode.setMetaData(nextEvent.node.getMetaData());
			case COALESCENCE:

			}
		}
		return null;
	}


	private double getTotalForwardsRate(int fromType, double t, int lineageIdx,
			HermiteInterpolator interpolator, double[] rates) {
		double totalRate = 0.0;
		getForwardsRates(fromType, t, lineageIdx, interpolator, rates);
		for (int type = 0; type < score.types; type++)
			totalRate += rates[type];

		return totalRate;
	}

	private double[] getForwardsRates(int fromType, double time, int lineageIdx,
			HermiteInterpolator interpolator, double[] result) {
		double[] p = interpolator.value(time);
		int interval = getIntervalIndex(time);
		double[] migMatrix = score.dynamics.getBackwardsMigration(interval);
		int n = (int) (Math.sqrt(migMatrix.length) + 0.5);

		for (int type = 0; type < score.types; type++) {
			if (type == fromType) {
				result[type] = 0.0;
				continue;
			}

			result[type] = migMatrix[type * n + fromType] * p[lineageIdx * score.types + type];
		}

		if (p[lineageIdx * score.types + fromType] <= 0.0) {
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
				result[type] /= p[lineageIdx * score.types + fromType];

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
}
