package score.simulator;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;

import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;


public class SimulateNetworkGivenRates extends Network {

    // array of reassortment rates for each state
	public Input<Function> reassortmentRatesInput = new Input<>("reassortmentRate",
	    "Rate of reassortment for each state (per lineage per unit time)", Validate.REQUIRED);

    // array of migration rates for each state
	public Input<Function> migrationRatesInput = new Input<>("migrationRate",
	    "Rate of migration for each state (per lineage per unit time)", Validate.REQUIRED);

//    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
//	    "Population model to use.", Validate.REQUIRED);

	public Input<Function> coalescentRatesInput = new Input<>("coalescentRate",
	    "Rate of migration for each state (per lineage per unit time)", Validate.REQUIRED);

    public Input<Integer> dimensionInput = new Input<>("dimension",
	    "the number of different states." + " if -1, it will use the number of different types ", -1);

    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set. ", Validate.REQUIRED);


    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree", "One or more segment trees to initialize.",
	    new ArrayList<>());

    public Input<Integer> nSegmentsInput = new Input<>("nSegments",
	    "Number of segments. Used if no segment trees are supplied.");

    public Input<TraitSet> traitSetInput = new Input<>("traitSet", "Trait set used to assign leaf ages.");

    public Input<Boolean> enableSegTreeUpdateInput = new Input<>("enableSegmentTreeUpdate",
	    "If false, segment tree objects won't be updated to agree with simulated " + "network. (Default true.)",
	    true);

    final public Input<Boolean> ignoreMigrationNodes = new Input<>("ignoreMigrationNodes",
	    "Do not include migration nodes in network output", false);

    private RealParameter reassortmentRates;
    private RealParameter migrationRates;

    private RealParameter coalescentRates;
    
	private List<NetworkNode> sampleNodes = new ArrayList<>();

    private final HashMap<String, Integer> typeNameToIndex = new HashMap<>();
	public final HashMap<Integer, String> typeIndexToName = new HashMap<>();

    private ArrayList<String> uniqueTypes;

    private enum MigrationType {
	symmetric, asymmetric
    }

    private MigrationType migrationType;

    private int nSegments;

    @Override
    public void initAndValidate() {

	if (nSegmentsInput.get() != null)
	    nSegments = nSegmentsInput.get();
	else
	    nSegments = segmentTreesInput.get().size();

//	populationFunction = populationFunctionInput.get();
		reassortmentRates = new RealParameter(doubleToDouble(reassortmentRatesInput.get().getDoubleValues()));
		reassortmentRates.setDimension(dimensionInput.get());
		migrationRates = new RealParameter(doubleToDouble(migrationRatesInput.get().getDoubleValues()));
		migrationRates.setDimension(dimensionInput.get());
		coalescentRates = new RealParameter(doubleToDouble(coalescentRatesInput.get().getDoubleValues()));
		coalescentRates.setDimension(dimensionInput.get());

	if (nSegments == 0) {
	    throw new IllegalArgumentException("Need at least one segment!");
	}

	// Set up sample nodes:



	final TraitSet leafAgeTraitSet = traitSetInput.get();
	final TraitSet typeTraitSet = typeTraitInput.get();
	TaxonSet taxonSet;
	if (leafAgeTraitSet != null)
	    taxonSet = leafAgeTraitSet.taxaInput.get();
	else
	    taxonSet = typeTraitSet.taxaInput.get();

	if (taxonSet == null)
	    throw new IllegalArgumentException("Must define either a " + "trait set, type set or a taxon set.");

	final SortedSet<String> typeNameSet = new TreeSet<>();
	taxonSet.asStringList().forEach(n -> typeNameSet.add(typeTraitSet.getStringValue(n)));
	uniqueTypes = new ArrayList<>(typeNameSet);

	for (int i = 0; i < uniqueTypes.size(); i++) {
	    typeNameToIndex.put(uniqueTypes.get(i), i);
	    typeIndexToName.put(i, uniqueTypes.get(i));
	}
	
	if (coalescentRates.getDimension() != uniqueTypes.size())
		coalescentRates.setDimension(uniqueTypes.size());
	
	if (reassortmentRates.getDimension() != uniqueTypes.size())
		reassortmentRates.setDimension(uniqueTypes.size());

		final int migDim = dimensionInput.get() != -1 ? dimensionInput.get() * (dimensionInput.get() - 1)
		: uniqueTypes.size() * (uniqueTypes.size() - 1);


	if (migDim == migrationRates.getDimension()) {
	    migrationType = MigrationType.asymmetric;
	} else if (migDim / 2 == migrationRates.getDimension()) {
	    migrationType = MigrationType.symmetric;
	} else {
	    migrationType = MigrationType.asymmetric;
	    System.err.println("Wrong number of migration elements, assume asymmetric migration:");
	    System.err.println("the dimension of " + migrationRates.getID() + " is set to " + migDim);
	    migrationRates.setDimension(migDim);
	}

	for (int taxonIndex = 0; taxonIndex < taxonSet.getTaxonCount(); taxonIndex++) {
	    final String taxonName = taxonSet.getTaxonId(taxonIndex);

			final NetworkNode sampleNode = new NetworkNode();
	    sampleNode.setTaxonLabel(taxonName);
	    sampleNode.setTaxonIndex(taxonIndex);

	    if (leafAgeTraitSet != null)
		sampleNode.setHeight(leafAgeTraitSet.getValue(taxonName));
	    else
		sampleNode.setHeight(0.0);

	    final String typeName = typeTraitSet.getStringValue(taxonName);
	    sampleNode.setTypeLabel(typeName);
	    sampleNode.setTypeIndex(typeNameToIndex.get(typeTraitSet.getStringValue(taxonName)));

	    sampleNodes.add(sampleNode);
	}

//	// Perform network simulation:
//	simulateNetwork(sampleNodes);
//
//	// Update segment trees:
//	if (enableSegTreeUpdateInput.get()) {
//	    for (int segIdx = 0; segIdx < nSegments; segIdx++) {
//		Tree segmentTree = segmentTreesInput.get().get(segIdx);
//		updateSegmentTree(segmentTree, segIdx);
//		segmentTree.setEverythingDirty(false);
//	    }
//	}
//
//	if (ignoreMigrationNodes.get())
//	    removeMigrationNodes();
    }

    /**
     * Simulate network under structured coalescent with reassortment model.
     * 
     * @param sampleNodes network nodes corresponding to samples.
     */
    public void simulateNetwork(List<NetworkNode> sampleNodes) {

		reassortmentRates = new RealParameter(doubleToDouble(reassortmentRatesInput.get().getDoubleValues()));
		reassortmentRates.setDimension(dimensionInput.get());
		migrationRates = new RealParameter(doubleToDouble(migrationRatesInput.get().getDoubleValues()));
		migrationRates.setDimension(dimensionInput.get());
		coalescentRates = new RealParameter(doubleToDouble(coalescentRatesInput.get().getDoubleValues()));
		coalescentRates.setDimension(dimensionInput.get());

	final List<NetworkNode> remainingSampleNodes = new ArrayList<>(sampleNodes);

	// #####################################
	// extant lineages have to be sorted by state id
	// #####################################

	final HashMap<Integer, List<NetworkEdge>> extantLineages = new HashMap<Integer, List<NetworkEdge>>(
		uniqueTypes.size());
	for (int i = 0; i < uniqueTypes.size(); i++) {
	    extantLineages.put(i, new ArrayList<>());
	}

	remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

	double currentTime = 0;
	double timeUntilNextSample;
	List<List<NetworkEdge>> remaining;
	do {
	    // get the timing of the next sampling event
	    if (!remainingSampleNodes.isEmpty()) {
		timeUntilNextSample = remainingSampleNodes.get(0).getHeight() - currentTime;
	    } else {
		timeUntilNextSample = Double.POSITIVE_INFINITY;
	    }

	    // TODO make work for different pop models
	    // assume fixed population for now, so transformation like this not needed:
//	     double currentTransformedTime = populationFunction.getIntensity(currentTime);
//	     double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
//	     double timeToNextCoal = populationFunction.getInverseIntensity(
//	     transformedTimeToNextCoal + currentTransformedTime) - currentTime;

	    double minCoal = Double.POSITIVE_INFINITY;
	    double minReassort = Double.POSITIVE_INFINITY;
	    double minMigration = Double.POSITIVE_INFINITY;

	    int typeIndexCoal = -1, typeIndexReassortment = -1, typeIndexMigrationFrom = -1, typeIndexMigrationTo = -1;

	    int c = 0;
	    for (int i = 0; i < uniqueTypes.size(); i++) {
		// how many lineages are in this state
		final int k_ = extantLineages.get(i).size();

		if (k_ >= 2) {
//		    double currentTransformedTime = populationFunction.getIntensity(currentTime);
//		    double transformedTimeToNextCoal = Randomizer.nextExponential(0.5 * k_ * (k_ - 1));
//		    double timeToNextCoal = populationFunction
//			    .getInverseIntensity(transformedTimeToNextCoal + currentTransformedTime) - currentTime;
		    final double timeToNextCoal = Randomizer
			    .nextExponential(0.5 * k_ * (k_ - 1) * coalescentRates.getArrayValue(i));
		    if (timeToNextCoal < minCoal) {
			minCoal = timeToNextCoal;
			typeIndexCoal = i;
		    }
		}

		if (k_ >= 1) {
		    final double timeToNextReass = Randomizer.nextExponential(k_ * reassortmentRates.getArrayValue(i));
		    if (timeToNextReass < minReassort) {
			minReassort = timeToNextReass;
			typeIndexReassortment = i;
		    }

		    for (int j = 0; j < uniqueTypes.size(); j++) {
			if (i != j) {
			    final double timeToNextMigration = Randomizer.nextExponential(k_ * migrationRates.getArrayValue(c));
			    c++;
			    if (migrationType == MigrationType.symmetric)
				c %= migrationRates.getDimension();
			    if (timeToNextMigration < minMigration) {
				minMigration = timeToNextMigration;
				typeIndexMigrationFrom = i;
				typeIndexMigrationTo = j;
			    }
			}
		    }
		}
	    }

	    // next event time
	    double timeUntilNextEvent = Math.min(minCoal, minReassort);
	    timeUntilNextEvent = Math.min(timeUntilNextEvent, minMigration);
	    if (timeUntilNextEvent < timeUntilNextSample) {
		currentTime += timeUntilNextEvent;
		if (timeUntilNextEvent == minCoal)
		    coalesce(currentTime, extantLineages, typeIndexCoal);
		else if (timeUntilNextEvent == minReassort)
		    reassort(currentTime, extantLineages, typeIndexReassortment);
		else
		    migrate(currentTime, extantLineages, typeIndexMigrationFrom, typeIndexMigrationTo);
	    } else {
		currentTime += timeUntilNextSample;
		sample(remainingSampleNodes, extantLineages);
	    }

	    remaining = extantLineages.values().stream().filter(l -> l.size() >= 1).collect(Collectors.toList());

	} while ((remaining.size() > 1 || remaining.get(0).size() > 1) || !remainingSampleNodes.isEmpty());

	final List<List<NetworkEdge>> root = extantLineages.values().stream().filter(l -> l.size() == 1)
		.collect(Collectors.toList());
	if (root.size() > 1)
	    System.err.println("More than one root edge");
	setRootEdge(root.get(0).get(0));
    }

    private void sample(List<NetworkNode> remainingSampleNodes, HashMap<Integer, List<NetworkEdge>> extantLineages) {
	// sample the network node
	final NetworkNode n = remainingSampleNodes.get(0);

	// Create corresponding lineage
	final BitSet hasSegs = new BitSet();
	hasSegs.set(0, nSegments);
	final NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
	final int id = n.getTypeIndex();
	extantLineages.get(id).add(lineage);
	n.addParentEdge(lineage);

	remainingSampleNodes.remove(0);
    }

    private void coalesce(double coalescentTime, HashMap<Integer, List<NetworkEdge>> extantLineages, int stateIdCoal) {
	// Sample the pair of lineages that are coalescing:
	final NetworkEdge lineage1 = extantLineages.get(stateIdCoal)
		.get(Randomizer.nextInt(extantLineages.get(stateIdCoal).size()));
	NetworkEdge lineage2;
	do {
	    lineage2 = extantLineages.get(stateIdCoal).get(Randomizer.nextInt(extantLineages.get(stateIdCoal).size()));
	} while (lineage1 == lineage2);

	// Create coalescent node
	final NetworkNode coalescentNode = new NetworkNode();
	coalescentNode.setHeight(coalescentTime).addChildEdge(lineage1).addChildEdge(lineage2);
	coalescentNode.setTypeIndex(stateIdCoal);
	coalescentNode.setTypeLabel(uniqueTypes.get(stateIdCoal));
	lineage1.parentNode = coalescentNode;
	lineage2.parentNode = coalescentNode;

	// Merge segment flags:
	final BitSet hasSegments = new BitSet();
	hasSegments.or(lineage1.hasSegments);
	hasSegments.or(lineage2.hasSegments);

	// Create new lineage
	final NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
	coalescentNode.addParentEdge(lineage);

	extantLineages.get(stateIdCoal).remove(lineage1);
	extantLineages.get(stateIdCoal).remove(lineage2);
	extantLineages.get(stateIdCoal).add(lineage);
    }

    private void reassort(double reassortmentTime, HashMap<Integer, List<NetworkEdge>> extantLineages,
	    int stateIdReassortment) {
	final NetworkEdge lineage = extantLineages.get(stateIdReassortment)
		.get(Randomizer.nextInt(extantLineages.get(stateIdReassortment).size()));

	final BitSet hasSegs_left = new BitSet();
	final BitSet hasSegs_right = new BitSet();

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
	final NetworkNode node = new NetworkNode();
	node.setHeight(reassortmentTime).addChildEdge(lineage);
	node.setTypeIndex(lineage.childNode.getTypeIndex());
	node.setTypeLabel(lineage.childNode.getTypeLabel());

	// Create reassortment lineages
	final NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
	final NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
	node.addParentEdge(leftLineage);
	node.addParentEdge(rightLineage);

	extantLineages.get(stateIdReassortment).remove(lineage);
	extantLineages.get(stateIdReassortment).add(leftLineage);
	extantLineages.get(stateIdReassortment).add(rightLineage);
    }

    private void migrate(double migrationTime, HashMap<Integer, List<NetworkEdge>> extantLineages,
	    int stateIdMigrationFrom, int stateIdMigrationTo) {
	// Sample the lineage for migration:
	final NetworkEdge lineage = extantLineages.get(stateIdMigrationFrom)
		.get(Randomizer.nextInt(extantLineages.get(stateIdMigrationFrom).size()));

	final NetworkNode migrationPoint = new NetworkNode();
	final NetworkEdge newParentEdge = new NetworkEdge();
	// NetworkEdge newParentEdge = lineage.getCopy();
	newParentEdge.hasSegments = lineage.hasSegments;

	migrationPoint.setHeight(migrationTime);
	migrationPoint.addParentEdge(newParentEdge);

	migrationPoint.addChildEdge(lineage);

	migrationPoint.setTypeIndex(stateIdMigrationTo);
	migrationPoint.setTypeLabel(uniqueTypes.get(stateIdMigrationTo));

	extantLineages.get(stateIdMigrationFrom).remove(lineage);
	extantLineages.get(stateIdMigrationTo).add(newParentEdge);
    }

    private void removeMigrationNodes() {

	List<NetworkNode> migrationNodes = this.getNodes().stream().filter(n -> n.getParentCount() == 1)
		.filter(n -> n.getChildCount() == 1).sorted(Comparator.comparing(NetworkNode::getHeight))
		.collect(Collectors.toList());

	Collections.reverse(migrationNodes);

	for (NetworkNode m : migrationNodes) {
	    NetworkEdge parentEdge = m.getParentEdges().get(0);
	    NetworkEdge childEdge = m.getChildEdges().get(0);

	    NetworkNode newChildNode = null;
//			while (newChildNode == null) {
	    if (childEdge.childNode.getChildCount() > 1 || childEdge.childNode.getParentCount() > 1
		    || childEdge.childNode.isLeaf()) {
		newChildNode = childEdge.childNode;
		m.removeChildEdge(childEdge);
		m.removeParentEdge(parentEdge);
		newChildNode.addParentEdge(parentEdge);
		newChildNode.removeParentEdge(childEdge);

	    } else {

		childEdge = childEdge.childNode.getChildEdges().get(0);
		m.removeChildEdge(childEdge.parentNode.getParentEdges().get(0));
		m.removeParentEdge(parentEdge);
		childEdge.parentNode.removeParentEdge(childEdge.parentNode.getParentEdges().get(0));
		childEdge.parentNode.addParentEdge(parentEdge);
		m = childEdge.parentNode;
		parentEdge = m.getParentEdges().get(0);
	    }
//			}

	}

    }

	private Double[] doubleToDouble(double[] arr) {
		final Double[] convertArray = new Double[arr.length];

		for (int index = 0; index < arr.length; index++)
			convertArray[index] = Double.valueOf(arr[index]);

		return convertArray;
	}
	
	@Override
	public void log (long sample, PrintStream out) {
		// Perform network simulation:
		simulateNetwork(sampleNodes);

		// Update segment trees:
		if (enableSegTreeUpdateInput.get()) {
		    for (int segIdx = 0; segIdx < nSegments; segIdx++) {
			Tree segmentTree = segmentTreesInput.get().get(segIdx);
			updateSegmentTree(segmentTree, segIdx);
			segmentTree.setEverythingDirty(false);
		    }
		}

		if (ignoreMigrationNodes.get())
		    removeMigrationNodes();
		
		out.print("tree STATE_" + sample + " = ");
		// Don't sort, this can confuse CalculationNodes relying on the tree
		// tree.getRoot().sort();
		final String newick = this.getExtendedNewick();
		out.print(newick);
	}

}
