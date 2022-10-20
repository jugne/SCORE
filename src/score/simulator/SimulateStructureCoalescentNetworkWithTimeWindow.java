package score.simulator;

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

import java.util.*;
import java.util.stream.Collectors;
import org.javatuples.Pair;


public class SimulateStructureCoalescentNetworkWithTimeWindow extends Network {

    // array of reassortment rates for each state
    public Input<RealParameter> reassortmentRatesInput = new Input<>("reassortmentRate",
	    "Rate of reassortment for each state (per lineage per unit time)", Validate.REQUIRED);

	public Input<RealParameter> scalerInput = new Input<>("scaler", "scalar for reassortment rate", Validate.REQUIRED);

	public Input<Double> timeWindowInput = new Input<>("timeWindow", "length of time window after migration event (backwards in time)", Validate.REQUIRED);

    // array of migration rates for each state
    public Input<RealParameter> migrationRatesInput = new Input<>("migrationRate",
	    "Rate of migration for each state (per lineage per unit time)", Validate.REQUIRED);

//    public Input<PopulationFunction> populationFunctionInput = new Input<>("populationModel",
//	    "Population model to use.", Validate.REQUIRED);

    public Input<RealParameter> coalescentRatesInput = new Input<>("coalescentRate",
			"Rate of coalescence for each state (per lineage per unit time)", Validate.REQUIRED);

    public Input<Integer> dimensionInput = new Input<>("dimension",
	    "the number of different states." + " if -1, it will use the number of different types ", -1);

    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set. ", Validate.REQUIRED);

    public Input<String> typesInput = new Input<>("types",
	    "input of the different types, can be helpful for multilocus data");

    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree", "One or more segment trees to initialize.",
	    new ArrayList<>());

    public Input<Integer> nSegmentsInput = new Input<>("nSegments",
	    "Number of segments. Used if no segment trees are supplied.");

    public Input<TraitSet> traitSetInput = new Input<>("traitSet", "Trait set used to assign leaf ages.");

    public Input<Boolean> enableSegTreeUpdateInput = new Input<>("enableSegmentTreeUpdate",
	    "If false, segment tree objects won't be updated to agree with simulated " + "network. (Default true.)",
	    true);

    public Input<String> fileNameInput = new Input<>("fileName", "Name of file to write simulated network to.");

    public Input<Double> ploidyInput = new Input<>("ploidy",
	    "Ploidy (copy number) for this gene," + "typically a whole number or half (default is 1).", 1.0);

    final public Input<Boolean> ignoreMigrationNodes = new Input<>("ignoreMigrationNodes",
	    "Do not include migration nodes in network output", false);

    private RealParameter reassortmentRates;

	private RealParameter scaler;
	private double timeWindow;

    private RealParameter migrationRates;
//    private PopulationFunction populationFunction;
    private RealParameter coalescentRates;
    private final HashMap<String, Integer> typeNameToIndex = new HashMap<>();
	public final HashMap<Integer, String> typeIndexToName = new HashMap<>();

    private ArrayList<String> uniqueTypes;

    private enum MigrationType {
	symmetric, asymmetric
    }

    private MigrationType migrationType;

    private int nSegments;

	double minChangeTime;

    @Override
    public void initAndValidate() {

		minChangeTime = Double.POSITIVE_INFINITY;

	if (nSegmentsInput.get() != null)
	    nSegments = nSegmentsInput.get();
	else
	    nSegments = segmentTreesInput.get().size();

//	populationFunction = populationFunctionInput.get();
	reassortmentRates = reassortmentRatesInput.get();
	scaler = scalerInput.get();
	timeWindow = timeWindowInput.get();
	migrationRates = migrationRatesInput.get();
	coalescentRates = coalescentRatesInput.get();

	if (nSegments == 0) {
	    throw new IllegalArgumentException("Need at least one segment!");
	}

	// Set up sample nodes:

	final List<NetworkNode> sampleNodes = new ArrayList<>();

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

    }

    /**
     * Simulate network under structured coalescent with reassortment model.
     * 
     * @param sampleNodes network nodes corresponding to samples.
     */
    public void simulateNetwork(List<NetworkNode> sampleNodes) {

	final List<NetworkNode> remainingSampleNodes = new ArrayList<>(sampleNodes);

	// #####################################
	// extant lineages have to be sorted by state id
	// #####################################

	final HashMap<Integer, List<Pair>> extantLineages = new HashMap<>(
			uniqueTypes.size() * 2);

	for (int i = 0; i < uniqueTypes.size() * 2; i++) {
	    extantLineages.put(i, new ArrayList<>());
	}

	remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

	double currentTime = 0;
	double timeUntilNextSample;
	List<List<Pair>> remaining;
	do {

		clearTimeWindow(extantLineages, uniqueTypes.size(), currentTime);

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

	    int c = 0; int d = 0;
	    for (int i = 0; i < uniqueTypes.size() * 2; i++) {
			boolean psudoTypeFlag;
			if (i % 2 == 0) {
				psudoTypeFlag = false;
			} else {
				psudoTypeFlag = true;
			}

			// how many lineages are in this state
			final int k_ = extantLineages.get(i).size();
			int k_pseudoType = 0;
			if (!psudoTypeFlag) {
				k_pseudoType = extantLineages.get(i + 1).size(); // Get the # of lineages in type i'
			}


			if ((k_ + k_pseudoType) >= 2 && !psudoTypeFlag) { // !psudoTypeFlag is added here to prevent sampling coalescent twice (one for type i and one for type i')
				// coalescent event for event i as a whole (including type i and type i')
				final double timeToNextCoal = Randomizer
					.nextExponential(0.5 * (k_ + k_pseudoType) * ((k_ + k_pseudoType) - 1) * coalescentRates.getArrayValue(i/2));
				if (timeToNextCoal < minCoal) {
				minCoal = timeToNextCoal;
				typeIndexCoal = i;
				}
			}

			if ((k_ + k_pseudoType) >= 1 && !psudoTypeFlag) {
				double timeToNextReass;
				double r = k_ * reassortmentRates.getArrayValue(i / 2);
				double r_pseudo = k_pseudoType * reassortmentRates.getArrayValue(i / 2) * scaler.getArrayValue(i / 2);
				timeToNextReass = Randomizer.nextExponential(r + r_pseudo); // type i, use normal reassortment rate


				if (timeToNextReass < minReassort) {
					minReassort = timeToNextReass;

					double u = Randomizer.nextDouble() * (r + r_pseudo);
					typeIndexReassortment = u < r ? i : (i + 1);
				}


				for (int j = 1; j < uniqueTypes.size()*2; j+=2) { // one can only migrate to pseudotype
					if (i != (j - 1)) {
						final double timeToNextMigration = Randomizer.nextExponential((k_ + k_pseudoType) * migrationRates.getArrayValue(c));

						if (timeToNextMigration < minMigration) {
							minMigration = timeToNextMigration;
							double u = Randomizer.nextDouble() * ((k_ + k_pseudoType) * migrationRates.getArrayValue(c));
							typeIndexMigrationFrom = u < k_*migrationRates.getArrayValue(c) ? i : (i+1);
							typeIndexMigrationTo = j;
						}
						c++;
						if (migrationType == MigrationType.symmetric)
							c %= migrationRates.getDimension();
					}
				}
			}
		}



	    // next event time
	    double timeUntilNextEvent = Math.min(minCoal, minReassort);
	    timeUntilNextEvent = Math.min(timeUntilNextEvent, minMigration);
		timeUntilNextEvent = Math.min(timeUntilNextEvent, minChangeTime);


		if (timeUntilNextEvent < timeUntilNextSample) {
			currentTime += timeUntilNextEvent;
			if (timeUntilNextEvent == minCoal)
				coalesce(currentTime, extantLineages, typeIndexCoal);
			else if (timeUntilNextEvent == minReassort)
				reassort(currentTime, extantLineages, typeIndexReassortment);
			else if (timeUntilNextEvent == minMigration)
				migrate(currentTime, extantLineages, typeIndexMigrationFrom, typeIndexMigrationTo, timeWindow);
	    } else {
			currentTime += timeUntilNextSample;
			sample(remainingSampleNodes, extantLineages);
	    }

	    remaining = extantLineages.values().stream().filter(l -> l.size() >= 1).collect(Collectors.toList());

	} while ((remaining.size() > 1 || remaining.get(0).size() > 1) || !remainingSampleNodes.isEmpty());

	final List<List<Pair>> root = extantLineages.values().stream().filter(l -> l.size() == 1)
		.collect(Collectors.toList());
	if (root.size() > 1)
	    System.err.println("More than one root edge");
	setRootEdge((NetworkEdge) root.get(0).get(0).getValue0());
    }

    private void sample(List<NetworkNode> remainingSampleNodes, HashMap<Integer, List<Pair>> extantLineages) {
		// sample the network node
		final NetworkNode n = remainingSampleNodes.get(0);

		// Create corresponding lineage
		final BitSet hasSegs = new BitSet();
		hasSegs.set(0, nSegments);
		final NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
		final int id = n.getTypeIndex();
		final Pair<NetworkEdge, Double> lineage_pair = Pair.with(lineage, -1.0);
		extantLineages.get(id * 2).add(lineage_pair);
		n.addParentEdge(lineage);

		remainingSampleNodes.remove(0);
    }

    private void coalesce(double coalescentTime, HashMap<Integer, List<Pair>> extantLineages, int stateIdCoal) {
		// Sample the pair of lineages that are coalescing:

		// Uniformly sample from type stateIdCoal and stateIdCoal'
		final int randomNumber = Randomizer.nextInt(extantLineages.get(stateIdCoal).size() + extantLineages.get(stateIdCoal + 1).size());

		boolean lineage1PseudoTypeFlag = false;

		Pair<NetworkEdge, Double> lineage1;
		if (randomNumber < extantLineages.get(stateIdCoal).size()) {
			lineage1 = extantLineages.get(stateIdCoal).get(randomNumber);
		} else {
			lineage1 = extantLineages.get(stateIdCoal + 1).get(randomNumber - extantLineages.get(stateIdCoal).size());
		}

		if (lineage1.getValue1() > -1) {
			lineage1PseudoTypeFlag = true;
		}

		int randomNumber2;
		boolean linege2PseudoTypeFlag = false;
		Pair<NetworkEdge, Double> lineage2;
		do {
			randomNumber2 = Randomizer.nextInt(extantLineages.get(stateIdCoal).size() + extantLineages.get(stateIdCoal + 1).size());
			if(randomNumber2 < extantLineages.get(stateIdCoal).size()) {
				lineage2 = extantLineages.get(stateIdCoal).get(randomNumber2);
			} else {
				lineage2 = extantLineages.get(stateIdCoal + 1).get(randomNumber2 - extantLineages.get(stateIdCoal).size());
			}

		} while (lineage1.equals(lineage2));

		if (lineage2.getValue1() > -1) {
			linege2PseudoTypeFlag = true;
		}

		// Create coalescent node
		final NetworkNode coalescentNode = new NetworkNode();
		coalescentNode.setHeight(coalescentTime).addChildEdge(lineage1.getValue0()).addChildEdge(lineage2.getValue0());
		coalescentNode.setTypeIndex(stateIdCoal / 2);
		coalescentNode.setTypeLabel(uniqueTypes.get(stateIdCoal / 2));
		lineage1.getValue0().parentNode = coalescentNode;
		lineage2.getValue0().parentNode = coalescentNode;

		// Merge segment flags:
		final BitSet hasSegments = new BitSet();
		hasSegments.or(lineage1.getValue0().hasSegments);
		hasSegments.or(lineage2.getValue0().hasSegments);

		// Create new lineage
		final NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
		coalescentNode.addParentEdge(lineage);
		double newLineageTimeWindow = Math.max(lineage1.getValue1(), lineage2.getValue1());
		boolean newLineagePsudoTypeFlag = false;
		if (newLineageTimeWindow > -1) {
			newLineagePsudoTypeFlag = true;
		}
		final Pair<NetworkEdge, Double> newLineage = Pair.with(lineage, newLineageTimeWindow);

		if (lineage1PseudoTypeFlag) {
			extantLineages.get(stateIdCoal + 1).remove(lineage1);
		} else {
			extantLineages.get(stateIdCoal).remove(lineage1);
		}

		if (linege2PseudoTypeFlag) {
			extantLineages.get(stateIdCoal + 1).remove(lineage2);
		} else {
			extantLineages.get(stateIdCoal).remove(lineage2);
		}

		if (newLineagePsudoTypeFlag) {
			extantLineages.get(stateIdCoal + 1).add(newLineage);
		} else {
			extantLineages.get(stateIdCoal).add(newLineage);
		}
    }

    private void reassort(double reassortmentTime, HashMap<Integer, List<Pair>> extantLineages,
	    int stateIdReassortment) {
		final Pair<NetworkEdge, Double> lineage = extantLineages.get(stateIdReassortment).get(Randomizer.nextInt(extantLineages.get(stateIdReassortment).size()));

		final BitSet hasSegs_left = new BitSet();
		final BitSet hasSegs_right = new BitSet();

		for (int segIdx = lineage.getValue0().hasSegments.nextSetBit(0); segIdx != -1; segIdx = lineage.getValue0().hasSegments
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
		node.setHeight(reassortmentTime).addChildEdge(lineage.getValue0());
		node.setTypeIndex(lineage.getValue0().childNode.getTypeIndex());
		node.setTypeLabel(lineage.getValue0().childNode.getTypeLabel());

		// Create reassortment lineages
		final NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
		final NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
		node.addParentEdge(leftLineage);
		node.addParentEdge(rightLineage);
		final Pair<NetworkEdge, Double> leftLineagePair = Pair.with(leftLineage, lineage.getValue1());
		final Pair<NetworkEdge, Double> rightLineagePair = Pair.with(rightLineage, lineage.getValue1());

		extantLineages.get(stateIdReassortment).remove(lineage);
		extantLineages.get(stateIdReassortment).add(leftLineagePair);
		extantLineages.get(stateIdReassortment).add(rightLineagePair);
    }

    private void migrate(double migrationTime, HashMap<Integer, List<Pair>> extantLineages,
	    int stateIdMigrationFrom, int stateIdMigrationTo, double timeWindow) {
		// Sample the lineage for migration:
		final Pair<NetworkEdge, Double> lineage = extantLineages.get(stateIdMigrationFrom).get(Randomizer.nextInt(extantLineages.get(stateIdMigrationFrom).size()));

		final NetworkNode migrationPoint = new NetworkNode();
		final NetworkEdge newParentEdge = new NetworkEdge();
		// NetworkEdge newParentEdge = lineage.getCopy();
		newParentEdge.hasSegments = lineage.getValue0().hasSegments;

		migrationPoint.setHeight(migrationTime);
		migrationPoint.addParentEdge(newParentEdge);

		migrationPoint.addChildEdge(lineage.getValue0());

		int typeIndex = (stateIdMigrationTo - 1) / 2;
		migrationPoint.setTypeIndex(typeIndex);
		migrationPoint.setTypeLabel(uniqueTypes.get(typeIndex));

		Pair<NetworkEdge, Double> newParentEdgePair = Pair.with(newParentEdge, migrationTime + timeWindow);
/*		if (minChangeTime>(migrationTime + timeWindow))
			minChangeTime = migrationTime + timeWindow;*/
		extantLineages.get(stateIdMigrationFrom).remove(lineage);
		extantLineages.get(stateIdMigrationTo).add(newParentEdgePair);
    }

	// clear up those lineages whose time window has passed
	private void clearTimeWindow(HashMap<Integer, List<Pair>> extantLineages, double typeNumber, double currentTime) {
		Double minWindow = Double.POSITIVE_INFINITY;
		for(int i = 1; i < typeNumber * 2; i += 2){
			int size_i = extantLineages.get(i).size();
			List<Pair> toBeDeleted = new ArrayList<Pair>();

			for(int j = 0; j < size_i; j++) {
				Pair<NetworkEdge, Double> lineage = extantLineages.get(i).get(j);
				final double timeWindow = lineage.getValue1();
				if (currentTime >= timeWindow && timeWindow != -1){
					Pair<NetworkEdge, Double> lineageNew = Pair.with(lineage.getValue0(), -1.0);
					// extantLineages.get(i).remove(lineage);
					toBeDeleted.add(lineage);
					extantLineages.get(i - 1).add(lineageNew);
				} else if (timeWindow != -1 && (timeWindow-currentTime)<minWindow){
					minWindow = timeWindow-currentTime;
				}
			}

			extantLineages.get(i).removeAll(toBeDeleted);
		}
		minChangeTime = minWindow;

	}

}