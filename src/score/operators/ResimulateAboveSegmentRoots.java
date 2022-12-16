package score.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.stream.Collectors;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.coalescent.PopulationFunction;
import beast.base.util.Randomizer;
import coalre.distribution.NetworkEvent;
import coalre.distribution.NetworkIntervals;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.operators.NetworkOperator;
import score.distribution.StructuredNetworkEvent;
import score.distribution.StructuredNetworkIntervals;

public class ResimulateAboveSegmentRoots extends NetworkOperator {

    public Input<RealParameter> reassortmentRateInput = new Input<>("reassortmentRates",
	    "Rate of reassortment (per lineage per unit time)", Validate.REQUIRED);

    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes", Validate.REQUIRED);
    
    public Input<StructuredNetworkIntervals> structuredNetworkIntervalsInput = new Input<>("networkIntervals",
            "Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);

    private int nSegments;

    private PopulationFunction populationFunction;
    private RealParameter reassortmentRate;
    private RealParameter Ne;
    private StructuredNetworkIntervals structuredIntervals;
    double meanReassotmentRate;
    double meanNe;
    int i=0;

    @Override
    public void initAndValidate() {
	nSegments = segmentTreesInput.get().size();

	reassortmentRate = reassortmentRateInput.get();
	Ne = NeInput.get();
	structuredIntervals = structuredNetworkIntervalsInput.get();
	super.initAndValidate();
    }

    @Override
    public double networkProposal() {
//    System.out.println("before: ");
//    System.out.println(network.getExtendedNewick());
	network.startEditing(this);
	// get the place where to cut
	double networkRootHeight = network.getRootEdge().childNode.getHeight();
	double maxHeight = getMaxSegmentMRCA();
//	if (networkRootHeight == maxHeight)
//		return Double.NEGATIVE_INFINITY;
	double currentSubNetProb = unstructuredSubNetworkProb(maxHeight);
	
//	i += 1;
//	System.out.println(i);
//	if (i==50) {
//		System.out.println("stop");
//	}
	double r = resimulate(maxHeight);
	if (r == Double.NEGATIVE_INFINITY) {
//		i += 1;
//		System.out.println(i);
		return r;
	}
	
	double newSubNetProb = unstructuredSubNetworkProb(maxHeight);
	return  currentSubNetProb-newSubNetProb;
    }

    double unstructuredSubNetworkProb(double startTime) {

	double prob = 0.0;
	meanReassotmentRate = calculate_average_of(reassortmentRate.getValues());
	meanNe = calculate_average_of(Ne.getValues());
	
//	System.out.println("Ne_full: "+Arrays.deepToString(Ne.getValues()));
//	System.out.println("Ne: "+meanNe);
//	System.out.println("Rea: "+meanReassotmentRate);
	
	List<StructuredNetworkEvent> subNetEventList = structuredIntervals.getNetworkEventList();
//	subNetEventList = subNetEventList.stream()
//			.filter(e -> e.time > startTime)
//			.collect(Collectors.toList());
	
	StructuredNetworkEvent prevEvent = null;

	for (StructuredNetworkEvent event : subNetEventList) {
	    if (prevEvent != null)
		prob += intervalContribution(prevEvent, event);

	    switch (event.type) {
	    case COALESCENCE:
		prob += coalesce(event);
		break;

	    case SAMPLE:
		break;

	    case REASSORTMENT:
		prob += reassortment(event);
		break;
	    }

	    if (prob == Double.NEGATIVE_INFINITY)
		break;

	    prevEvent = event;
	}

	return prob;
    }

    private double reassortment(StructuredNetworkEvent event) {

	return Math.log(meanReassotmentRate) + event.segsSortedLeft * Math.log(structuredIntervals.getBinomialProb())
		+ (event.segsToSort - event.segsSortedLeft) * Math.log(1 - structuredIntervals.getBinomialProb()) + Math.log(2.0);
    }

    private double coalesce(StructuredNetworkEvent event) {

	return Math.log(1.0 / meanNe);
    }

    private double intervalContribution(StructuredNetworkEvent prevEvent, StructuredNetworkEvent nextEvent) {

	double result = 0.0;

	result += -meanReassotmentRate * prevEvent.totalReassortmentObsProb
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

    double resimulate(double startTime) {
//	network.startEditing(this);

	// get the place where to cut
	double maxHeight = startTime;

	// get all network edges
	List<NetworkEdge> networkEdges = new ArrayList<>(network.getEdges());

	// keep only those that coexist at the time of maxHeight
	List<NetworkEdge> startingEdges = networkEdges.stream().filter(e -> !e.isRootEdge())
		.filter(e -> e.parentNode.getHeight() > maxHeight).filter(e -> e.childNode.getHeight() <= maxHeight)
		.collect(Collectors.toList());


	if (startingEdges.size() == 0)
	    return Double.NEGATIVE_INFINITY;

//        System.out.println(network.getExtendedNewick());
	// simulate the rest of the network starting from mxHeight
	double currentTime = maxHeight;
	double timeUntilNextSample = Double.POSITIVE_INFINITY;
	do {

	    // get the current propensities
	    int k = startingEdges.size();
	    
	    double timeToNextCoal = k >= 2 ? Randomizer
			    .nextExponential(0.5 * k * (k - 1) * (1/meanNe))
			    : Double.POSITIVE_INFINITY;;

	    double timeToNextReass = k >= 1 ? Randomizer.nextExponential(k * meanReassotmentRate)
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
//	System.out.println(network.getExtendedNewick());

	return Double.POSITIVE_INFINITY;

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

    
    double getMaxSegmentMRCA() {
	double maxHeight = 0.0;
	for (int i = 0; i < segmentTreesInput.get().size(); i++) {
	    double height = segmentTreesInput.get().get(i).getRoot().getHeight();
	    if (height > maxHeight)
		maxHeight = height;
	}

	return maxHeight;
    }
}
