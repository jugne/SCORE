package structured.simulator;

import java.util.*;
import java.util.stream.Collectors;

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

//TODO check length for rate of reassort and migrattion rates to be the same as number of states
//TODO figure out number of unique states like in Dynamics.java


public class SimulateStructureCoalescentNetwork extends Network {
	
	// array of reassortment rates for each state
    public Input<RealParameter> reassortmentRatesInput = new Input<>("reassortmentRate",
            "Rate of reassortment for each state (per lineage per unit time)", Validate.REQUIRED);
    
	// array of migration rates for each state
    public Input<RealParameter> migrationRatesInput = new Input<>("migrationRate",
            "Rate of migration for each state (per lineage per unit time)", Validate.REQUIRED);
    
    
    public Input<Integer> dimensionInput = new Input<>("dimension", "the number of different states." + 
    		" if -1, it will use the number of different types ", -1);
    
    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set. ", Validate.REQUIRED);
    
    public Input<String> typesInput = new Input<>(
    		"types", "input of the different types, can be helpful for multilocus data");
    
    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "One or more segment trees to initialize.", new ArrayList<>());

    public Input<Integer> nSegmentsInput = new Input<>("nSegments",
            "Number of segments. Used if no segment trees are supplied.");

    public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set used to assign leaf ages.");

    public Input<Boolean> enableSegTreeUpdateInput = new Input<>("enableSegmentTreeUpdate",
            "If false, segment tree objects won't be updated to agree with simulated " +
                    "network. (Default true.)", true);

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file to write simulated network to.");
    
    public Input<RealParameter> NeInput = new Input<>("Ne", "input of effective population sizes");
    
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene,"
    		+ "typically a whole number or half (default is 1).", 1.0);
	
    
    private RealParameter reassortmentRates;
    private RealParameter migrationRates;
    private RealParameter Ne;
    private HashMap<String, Integer> typeNameToIndex = new HashMap<>();
    private HashMap<Integer, String> typeIndexToName = new HashMap<>();

    private ArrayList<String> uniqueTypes;
    
	private enum MigrationType {
	    symmetric, asymmetric 
	}

    private MigrationType migrationType;
    
    

    private int nSegments;
    
    public void initAndValidate() {
    	
    	Ne = NeInput.get();
        if (nSegmentsInput.get() != null)
            nSegments = nSegmentsInput.get();
        else
            nSegments = segmentTreesInput.get().size();

//        populationFunction = populationFunctionInput.get();
        reassortmentRates = reassortmentRatesInput.get();
        migrationRates = migrationRatesInput.get();

        if (nSegments==0) {
            throw new IllegalArgumentException("Need at least one segment!");
        }
        
        int migDim = Ne.getDimension()*(Ne.getDimension()-1);
        
        if (migDim == migrationRates.getDimension()){
			migrationType = MigrationType.asymmetric;
		}else if ((int) migDim/2 == migrationRates.getDimension()){
			migrationType = MigrationType.symmetric;
		}else{
			migrationType = MigrationType.asymmetric;
			System.err.println("Wrong number of migration elements, assume asymmetric migration:");
    		System.err.println("the dimension of " + migrationRates.getID() + " is set to " + migDim);
    		migrationRates.setDimension(migDim);       		
		}
//        System.out.println(migrationType);
        // Set up sample nodes:
        
        List<NetworkNode> sampleNodes = new ArrayList<>();
        
        TraitSet leafAgeTraitSet = traitSetInput.get();
        TraitSet typeTraitSet = typeTraitInput.get();
        TaxonSet taxonSet;
        if (leafAgeTraitSet != null)
            taxonSet = leafAgeTraitSet.taxaInput.get();
        else
        	taxonSet = typeTraitSet.taxaInput.get();

        if (taxonSet == null)
                throw new IllegalArgumentException("Must define either a " +
                        "trait set, type set or a taxon set.");

        uniqueTypes = new ArrayList<>(new TreeSet<>(taxonSet.asStringList()));

		for (int i = 0; i < uniqueTypes.size(); i++) {
			typeNameToIndex.put(uniqueTypes.get(i), i);
			typeIndexToName.put(i, uniqueTypes.get(i));
		}
			
        for (int taxonIndex=0; taxonIndex<taxonSet.getTaxonCount(); taxonIndex++) {
            String taxonName = taxonSet.getTaxonId(taxonIndex);

            NetworkNode sampleNode = new NetworkNode();
            sampleNode.setTaxonLabel(taxonName);
            sampleNode.setTaxonIndex(taxonIndex);

            if (leafAgeTraitSet != null)
                sampleNode.setHeight(leafAgeTraitSet.getValue(taxonName));
            else
                sampleNode.setHeight(0.0);
            
            String typeName = typeTraitSet.getStringValue(taxonName);
            sampleNode.setTypeLabel(typeName);
            sampleNode.setTypeIndex(typeNameToIndex.get(typeTraitSet.getStringValue(taxonName)));

            sampleNodes.add(sampleNode);
        }
        
        // Perform network simulation:
        simulateNetwork(sampleNodes);
    }
    
    
    
    /**
     * Simulate network under structured coalescent with reassortment model.
     * @param sampleNodes network nodes corresponding to samples.
     */
    public void simulateNetwork(List<NetworkNode> sampleNodes) {

        List<NetworkNode> remainingSampleNodes = new ArrayList<>(sampleNodes);
        
//        #####################################
//        extant lineages have to be sorted by state id
//        #####################################
        
        HashMap<Integer, List<NetworkEdge>> extantLineages = new HashMap<Integer, List<NetworkEdge>>(uniqueTypes.size());
        for (int i = 0; i < uniqueTypes.size(); i++ ) {
        	extantLineages.put(i, new ArrayList<>());        	
        }

        remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

        double currentTime = 0;
        double timeUntilNextSample;
        boolean allEmpty;
        List<List <NetworkEdge>> remaining;
        do {
            // get the timing of the next sampling event
            if (!remainingSampleNodes.isEmpty()) {
                timeUntilNextSample = remainingSampleNodes.get(0).getHeight() - currentTime;
            } else {
                timeUntilNextSample = Double.POSITIVE_INFINITY;
            }
            
            
//          TODO make work for different pop models
//            assume fixed population for now, so transformation like this not needed:
//            double currentTransformedTime = populationFunction.getIntensity(currentTime);
//            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
//            double timeToNextCoal = populationFunction.getInverseIntensity(
//                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;
            
            double minCoal = Double.POSITIVE_INFINITY;
            double minReassort = Double.POSITIVE_INFINITY;
            double minMigration = Double.POSITIVE_INFINITY;

            int typeIndexCoal = -1, typeIndexReassortment = -1, typeIndexMigrationFrom = -1, typeIndexMigrationTo = -1;

            int c = 0;
            for (int i = 0; i < uniqueTypes.size(); i++) {
            	//how many lineages are in this state
            	int k_ = extantLineages.get(i).size();
            	
            	if (k_ >= 2) {
            		double timeToNextCoal = Randomizer.nextExponential(0.5*k_*(k_-1));
            		if (timeToNextCoal < minCoal) {
            			minCoal = timeToNextCoal;
            			typeIndexCoal = i;
            		}
            	}
            	
            	if (k_ >= 1) {
            		double timeToNextReass = Randomizer.nextExponential(k_*reassortmentRates.getArrayValue(i));
                	if (timeToNextReass < minReassort) {
                		minReassort = timeToNextReass;
                		typeIndexReassortment = i;
                	}


//                	if (migrationType == MigrationType.asymmetric) {
                      	for (int j = 0; j< uniqueTypes.size(); j++) {
                    		if (i!=j) {
                    			double timeToNextMigration = Randomizer.nextExponential(k_*migrationRates.getArrayValue(c));
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
//                	}else {
//                      	for (int j=0; j<uniqueStates.size(); j++) {
//                        		if (i!=j) {
//                        			double timeToNextMigration = Randomizer.nextExponential(k_*migrationRates.getArrayValue(c));
//                              		c++;
//                              		c %= migrationRates.getDimension();
////                        					.getArrayValue(i*uniqueStates.size()+j));
//                                	if (timeToNextMigration < minMigration) {
//                                		minMigration = timeToNextMigration;
//                                		stateIdMigrationFrom = i;
//                                		stateIdMigrationTo = j;
//                                	}
//                        		}
//                      		}
//                      	}
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
                	migrate(currentTime, extantLineages,
                            typeIndexMigrationFrom, typeIndexMigrationTo);
            } else {
                currentTime += timeUntilNextSample;
                sample(remainingSampleNodes, extantLineages);
            }
        
//            allEmpty = extantLineages.values().stream().allMatch(l -> l == null || l.size() <= 1);
            remaining = extantLineages.values().stream().filter(l -> l.size() >= 1).collect(Collectors.toList());
            

        }
        while ((remaining.size() > 1 || remaining.get(0).size() >1) || !remainingSampleNodes.isEmpty());
        

        
        final List<List <NetworkEdge>> root = extantLineages.values().stream().filter(l -> l.size()==1).collect(Collectors.toList());
        if (root.size() > 1)
        	System.err.println("More than one root edge");
        setRootEdge(root.get(0).get(0));
    }
    
    
    private void sample(List<NetworkNode> remainingSampleNodes, HashMap<Integer, List<NetworkEdge>> extantLineages) {
        // sample the network node
        NetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BitSet hasSegs = new BitSet();
        hasSegs.set(0, nSegments);
        NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
        int id = n.getTypeIndex();
        extantLineages.get(id).add(lineage);
        n.addParentEdge(lineage);

        remainingSampleNodes.remove(0);
    }
    
    private void coalesce(double coalescentTime, HashMap<Integer, List<NetworkEdge>> extantLineages, int stateIdCoal) {
        // Sample the pair of lineages that are coalescing:
        NetworkEdge lineage1 = extantLineages.get(stateIdCoal).
        		get(Randomizer.nextInt(extantLineages.get(stateIdCoal).size()));
        NetworkEdge lineage2;
        do {
            lineage2 = extantLineages.get(stateIdCoal).
            		get(Randomizer.nextInt(extantLineages.get(stateIdCoal).size()));
        } while (lineage1 == lineage2);

        // Create coalescent node
        NetworkNode coalescentNode = new NetworkNode();
        coalescentNode.setHeight(coalescentTime)
                .addChildEdge(lineage1)
                .addChildEdge(lineage2);
        coalescentNode.setTypeIndex(stateIdCoal);
        coalescentNode.setTypeLabel(uniqueTypes.get(stateIdCoal));
        lineage1.parentNode = coalescentNode;
        lineage2.parentNode = coalescentNode;

        // Merge segment flags:
        BitSet hasSegments = new BitSet();
        hasSegments.or(lineage1.hasSegments);
        hasSegments.or(lineage2.hasSegments);

        // Create new lineage
        NetworkEdge lineage = new NetworkEdge(null, coalescentNode, hasSegments);
        coalescentNode.addParentEdge(lineage);

        extantLineages.get(stateIdCoal).remove(lineage1);
        extantLineages.get(stateIdCoal).remove(lineage2);
        extantLineages.get(stateIdCoal).add(lineage);
    }
    
    private void reassort(double reassortmentTime, HashMap<Integer, List<NetworkEdge>> extantLineages, int stateIdReassortment) {
        NetworkEdge lineage = extantLineages.get(stateIdReassortment).
        		get(Randomizer.nextInt(extantLineages.get(stateIdReassortment).size()));

        BitSet hasSegs_left = new BitSet();
        BitSet hasSegs_right = new BitSet();

        for (int segIdx = lineage.hasSegments.nextSetBit(0);
             segIdx != -1; segIdx = lineage.hasSegments.nextSetBit(segIdx+1)) {
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
        node.setTypeIndex(lineage.childNode.getTypeIndex());
        node.setTypeLabel(lineage.childNode.getTypeLabel());

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.get(stateIdReassortment).remove(lineage);
        extantLineages.get(stateIdReassortment).add(leftLineage);
        extantLineages.get(stateIdReassortment).add(rightLineage);
    }
    
    private void migrate(double migrationTime, HashMap<Integer, List<NetworkEdge>> extantLineages, int stateIdMigrationFrom, int stateIdMigrationTo) {
        // Sample the lineage for migration:
        NetworkEdge lineage = extantLineages.get(stateIdMigrationFrom).
        		get(Randomizer.nextInt(extantLineages.get(stateIdMigrationFrom).size()));

        NetworkNode migrationPoint = new NetworkNode();
//        NetworkNode newParentNode = lineage.parentNode;
        NetworkEdge newParentEdge = lineage.getCopy();
        
        migrationPoint.setHeight(migrationTime);
        migrationPoint.addParentEdge(newParentEdge);
//        newParentNode.addChildEdge(newParentEdge);
//        
//        newParentNode.removeChildEdge(lineage);
        migrationPoint.addChildEdge(lineage);
        
        migrationPoint.setTypeIndex(stateIdMigrationTo);
        migrationPoint.setTypeLabel(uniqueTypes.get(stateIdMigrationTo));
        
        extantLineages.get(stateIdMigrationFrom).remove(lineage);
        extantLineages.get(stateIdMigrationTo).add(newParentEdge);
    }

}
