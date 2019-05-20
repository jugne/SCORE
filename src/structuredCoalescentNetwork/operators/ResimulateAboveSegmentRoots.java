package structuredCoalescentNetwork.operators;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import coalre.distribution.NetworkEvent;
import coalre.distribution.NetworkIntervals;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.operators.NetworkOperator;

public class ResimulateAboveSegmentRoots extends NetworkOperator {

    public Input<RealParameter> reassortmentRateInput = new Input<>("reassortmentRate",
	    "Rate of reassortment (per lineage per unit time)", Validate.REQUIRED);

    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);

    public Input<NetworkIntervals> networkIntervalsInput = new Input<>("networkIntervals",
	    "Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);

    private int nSegments;

    private PopulationFunction populationFunction;
    private RealParameter reassortmentRate;
    private RealParameter Ne;
    private NetworkIntervals intervals;
    double meanReassotmentRate;
    double meanNe;
    double logP;

    @Override
    public void initAndValidate() {
	nSegments = segmentTreesInput.get().size();

	reassortmentRate = reassortmentRateInput.get();
	Ne = NeInput.get();
	intervals = networkIntervalsInput.get();
	super.initAndValidate();
    }

    @Override
    public double networkProposal() {
	network.startEditing(this);
	return resimulate();

    }

    double unstructuredSubNetworkProb(double startTime) {

	logP = 0.0;
	meanReassotmentRate = calculate_average_of(reassortmentRate.getValues());
	meanNe = calculate_average_of(Ne.getValues());

	List<NetworkEvent> subNetEventList = intervals.getNetworkEventList();

	subNetEventList = subNetEventList.stream().filter(e -> e.time > startTime).collect(Collectors.toList());
	NetworkEvent prevEvent = null;

	for (NetworkEvent event : subNetEventList) {
	    if (prevEvent != null)
		logP += intervalContribution(prevEvent, event);

	    switch (event.type) {
	    case COALESCENCE:
		logP += coalesce(event);
		break;

	    case SAMPLE:
		break;

	    case REASSORTMENT:
		logP += reassortment(event);
		break;
	    }

	    if (logP == Double.NEGATIVE_INFINITY)
		break;

	    prevEvent = event;
	}

	return logP;
    }

    private double reassortment(NetworkEvent event) {

	return Math.log(meanReassotmentRate) + event.segsSortedLeft * Math.log(intervals.getBinomialProb())
		+ (event.segsToSort - event.segsSortedLeft) * Math.log(1 - intervals.getBinomialProb()) + Math.log(2.0);
    }

    private double coalesce(NetworkEvent event) {

	return Math.log(1.0 / meanNe);
    }

    private double intervalContribution(NetworkEvent prevEvent, NetworkEvent nextEvent) {

	double result = 0.0;

	result += -reassortmentRate.getArrayValue() * prevEvent.totalReassortmentObsProb
		* (nextEvent.time - prevEvent.time);

	result += -0.5 * prevEvent.lineages * (prevEvent.lineages - 1) * (1.0 / meanNe)
		* (nextEvent.time - prevEvent.time);

	return result;
    }

    private double calculate_average_of(Double[] array) {

	double total = 0;
	for (double element : array) {
	    total += element;
	}

	return total / array.length;
    }

    double resimulate() {
	network.startEditing(this);

	// get the place where to cut
	double maxHeight = getMaxSegmentMRCA();

	// get all network edges
	List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

	// keep only those that coexist at the time of maxHeight
	List<NetworkEdge> startingEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
		.filter(e -> e.parentNode.getHeight() > maxHeight).filter(e -> e.childNode.getHeight() <= maxHeight)
		.collect(Collectors.toList());

//        System.out.println("max " + maxHeight);
//        for (int i = 0; i < startingEdges.size(); i++)
//        	System.out.println(startingEdges.get(i).childNode.getHeight());

	if (startingEdges.size() == 0)
	    return Double.NEGATIVE_INFINITY;

//        System.out.println(network.getExtendedNewick());
	// simulate the rest of the network starting from mxHeight
	double currentTime = maxHeight;
	double timeUntilNextSample = Double.POSITIVE_INFINITY;
	do {

	    // get the current propensities
	    int k = startingEdges.size();

	    double currentTransformedTime = populationFunction.getIntensity(currentTime);
	    double transformedTimeToNextCoal = k >= 2 ? Randomizer.nextExponential(0.5 * k * (k - 1))
		    : Double.POSITIVE_INFINITY;
	    double timeToNextCoal = populationFunction
		    .getInverseIntensity(transformedTimeToNextCoal + currentTransformedTime) - currentTime;

	    double timeToNextReass = k >= 1 ? Randomizer.nextExponential(k * reassortmentRate.getValue())
		    : Double.POSITIVE_INFINITY;

	    // next event time
	    double timeUntilNextEvent = Math.min(timeToNextCoal, timeToNextReass);
	    if (timeUntilNextEvent < timeUntilNextSample) {
		currentTime += timeUntilNextEvent;
		if (timeUntilNextEvent == timeToNextCoal)
		    coalesce(currentTime, startingEdges);
		else
		    reassort(currentTime, startingEdges);
	    }

	} while (startingEdges.size() > 1);

	network.setRootEdge(startingEdges.get(0));
//        System.out.println(network.getExtendedNewick());

	return Double.POSITIVE_INFINITY;

    }

    double getMaxSegmentMRCA() {
	double maxHeight = 0.0;
	for (int i = 0; i < segmentTreesInput.get().size(); i++) {
	    double height = segmentTreesInput.get().get(i).getRoot().getHeight();
	    if (height > maxHeight)
		maxHeight = height;
	}

	return maxHeight;
    }

    private void coalesce(double coalescentTime, List<NetworkEdge> extantLineages) {
	// Sample the pair of lineages that are coalescing:
	NetworkEdge lineage1 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
	NetworkEdge lineage2;
	do {
	    lineage2 = extantLineages.get(Randomizer.nextInt(extantLineages.size()));
	} while (lineage1 == lineage2);

	// Create coalescent node
	NetworkNode coalescentNode = new NetworkNode();
	coalescentNode.setHeight(coalescentTime).addChildEdge(lineage1).addChildEdge(lineage2);
	lineage1.parentNode = coalescentNode;
	lineage2.parentNode = coalescentNode;

	// Merge segment flags:
	BitSet hasSegments = new BitSet();
	hasSegments.or(lineage1.hasSegments);
	hasSegments.or(lineage2.hasSegments);

	// Create new lineage
	NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
	coalescentNode.addParentEdge(lineage);

	extantLineages.remove(lineage1);
	extantLineages.remove(lineage2);
	extantLineages.add(lineage);
    }

    private void sample(List<NetworkNode> remainingSampleNodes, List<NetworkEdge> extantLineages) {
	// sample the network node
	NetworkNode n = remainingSampleNodes.get(0);

	// Create corresponding lineage
	BitSet hasSegs = new BitSet();
	hasSegs.set(0, nSegments);
	NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
	extantLineages.add(lineage);
	n.addParentEdge(lineage);

	remainingSampleNodes.remove(0);
    }

    private void reassort(double reassortmentTime, List<NetworkEdge> extantLineages) {
	NetworkEdge lineage = extantLineages.get(Randomizer.nextInt(extantLineages.size()));

	BitSet hasSegs_left = new BitSet();
	BitSet hasSegs_right = new BitSet();

	for (int segIdx = lineage.hasSegments.nextSetBit(0); segIdx != -1; segIdx = lineage.hasSegments
		.nextSetBit(segIdx + 1)) {
	    if (Randomizer.nextBoolean()) {
		hasSegs_left.set(segIdx);
	    } else {
		hasSegs_right.set(segIdx);
	    }
	}

	// Stop here if reassortment event is unobservable
	if (hasSegs_left.cardinality() == 0 || hasSegs_right.cardinality() == 0)
	    return;

	// Create reassortment node
	NetworkNode node = new NetworkNode();
	node.setHeight(reassortmentTime).addChildEdge(lineage);

	// Create reassortment lineages
	NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
	NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
	node.addParentEdge(leftLineage);
	node.addParentEdge(rightLineage);

	extantLineages.remove(lineage);
	extantLineages.add(leftLineage);
	extantLineages.add(rightLineage);
    }

}
