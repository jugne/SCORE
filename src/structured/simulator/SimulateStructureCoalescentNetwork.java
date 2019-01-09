package structured.simulator;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.TraitSet;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.PopulationFunction;
import beast.util.Randomizer;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

//TODO check length for rate of reassort and migrattion rates to be the same as number of states
//TODO figure out number of unique states like in Dynamics.java


public class SimulateStructureCoalescentNetwork extends Network{
	
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
    		"types", "input of the different types, can be helpful for multilocus data", Validate.OPTIONAL);
    
    public Input<List<Tree>> segmentTreesInput = new Input<>("segmentTree",
            "One or more segment trees to initialize.", new ArrayList<>());

    public Input<Integer> nSegmentsInput = new Input<>("nSegments",
            "Number of segments. Used if no segment trees are supplied.");

    public Input<TraitSet> traitSetInput = new Input<>("traitSet",
            "Trait set used to assign leaf ages.");

    public Input<TaxonSet> taxonSetInput = new Input<>("taxonSet",
            "Taxon set used to define leaves");

    public Input<Boolean> enableSegTreeUpdateInput = new Input<>("enableSegmentTreeUpdate",
            "If false, segment tree objects won't be updated to agree with simulated " +
                    "network. (Default true.)", true);

    public Input<String> fileNameInput = new Input<>("fileName",
            "Name of file to write simulated network to.");
    
    public Input<RealParameter> NeImput = new Input<>("Ne", "input of effective population sizes"); 
    
    public Input<Double> ploidyInput = new Input<>("ploidy", "Ploidy (copy number) for this gene,"
    		+ "typically a whole number or half (default is 1).", 1.0);
	
    
//    private PopulationFunction populationFunction;
    private RealParameter reassortmentRates;
    private RealParameter migrationRates;
    private RealParameter Ne;
    private HashMap<String, Integer> traitToType = new HashMap<>(); 
    private HashMap<Integer, String> reverseTraitToType = new HashMap<>();
    private ArrayList<String> uniqueStates = new ArrayList<>();
    
	private enum MigrationType {
	    symmetric, asymmetric 
	}
	
	MigrationType migrationType;
    
    

    private int nSegments;
    
    public void initAndValidate() {
    	
    	Ne = NeImput.get();    	
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
        
        // Set up sample nodes:
        
        List<NetworkNode> sampleNodes = new ArrayList<>();
        
        TraitSet leafAgeSet = traitSetInput.get();
        TraitSet stateSet = typeTraitInput.get();
        TaxonSet taxonSet;
        if (leafAgeSet != null)
            taxonSet = leafAgeSet.taxaInput.get();
        else if (stateSet != null)
        	taxonSet = stateSet.taxaInput.get();
        else
            taxonSet = taxonSetInput.get();

        if (taxonSet == null)
                throw new IllegalArgumentException("Must define either a " +
                        "trait set, type set or a taxon set.");
    	
		List<String> taxa = stateSet.taxaInput.get().asStringList();
		for (int i = 0; i < taxa.size(); i++)
			uniqueStates.add(stateSet.getStringValue(taxa.get(i)));
		
		Collections.sort(uniqueStates);
		for (int i = uniqueStates.size()-2; i > -1; i--)
			if(uniqueStates.get(i+1).equals(uniqueStates.get(i)))
				uniqueStates.remove(i+1);
		  
		for (int i = 0; i < uniqueStates.size(); i++)
			traitToType.put(uniqueStates.get(i), i);
		for (int i = 0; i < uniqueStates.size(); i++) {
			reverseTraitToType.put(i, uniqueStates.get(i)); 
		}
			
        for (int taxonIndex=0; taxonIndex<taxonSet.getTaxonCount(); taxonIndex++) {
            String taxonName = taxonSet.getTaxonId(taxonIndex);

            NetworkNode sampleNode = new NetworkNode();
            sampleNode.setTaxonLabel(taxonName);
            sampleNode.setTaxonIndex(taxonIndex);

            if (leafAgeSet != null)
                sampleNode.setHeight(leafAgeSet.getValue(taxonName));
            else
                sampleNode.setHeight(0.0);
            
            String stateName = stateSet.getStringValue(taxonName);
            sampleNode.setStateLabel(stateName);
            sampleNode.setStateIndex(traitToType.get(stateSet.getStringValue(taxonName)));

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
        
        List<List<NetworkEdge>> extantLineages = new ArrayList<List<NetworkEdge>>(uniqueStates.size());

        remainingSampleNodes.sort(Comparator.comparingDouble(NetworkNode::getHeight));

        double currentTime = 0;
        double timeUntilNextSample;
        boolean allEmpty;
        do {
            // get the timing of the next sampling event
            if (!remainingSampleNodes.isEmpty()) {
                timeUntilNextSample = remainingSampleNodes.get(0).getHeight() - currentTime;
            } else {
                timeUntilNextSample = Double.POSITIVE_INFINITY;
            }
            
            // get the current propensities
//            int k = extantLineages.size();
            
//          TODO make work for different pop models
//            assume fixed population for now, so transformation like this not needed:
//            double currentTransformedTime = populationFunction.getIntensity(currentTime);
//            double transformedTimeToNextCoal = k>=2 ? Randomizer.nextExponential(0.5*k*(k-1)) : Double.POSITIVE_INFINITY;
//            double timeToNextCoal = populationFunction.getInverseIntensity(
//                    transformedTimeToNextCoal + currentTransformedTime) - currentTime;
            
//            int[] k_ = new int[uniqueStates.size()];
            double minCoal = Double.POSITIVE_INFINITY;
            double minReassort = Double.POSITIVE_INFINITY;
            double minMigration = Double.POSITIVE_INFINITY;
            Integer stateIdCoal = null;
            Integer stateIdReassortment = null;
            Integer stateIdMigrationFrom = null;
            Integer stateIdMigrationTo = null;
            for (int i = 0; i < uniqueStates.size(); i++) {
            	//how many lineages are in this state
            	int k_ = extantLineages.get(i).size();
            	
            	if (k_ >= 2) {
            		double timeToNextCoal = Randomizer.nextExponential(0.5*k_*(k_-1));
            		if (timeToNextCoal < minCoal) {
            			minCoal = timeToNextCoal;
            			stateIdCoal = i;
            		}
            	}
            	
            	if (k_ >= 1) {
            		double timeToNextReass = Randomizer.nextExponential(k_*reassortmentRates.getArrayValue(i));
                	if (timeToNextReass < minReassort) {
                		minReassort = timeToNextReass;
                		stateIdReassortment = i;
                	}


                	if (migrationType == MigrationType.asymmetric) {
                      	for (int j=0; j<uniqueStates.size(); j++) {
                    		if (i!=j) {
                    			double timeToNextMigration = Randomizer.nextExponential(k_*migrationRates
                    					.getArrayValue(i*uniqueStates.size()+j));
                            	if (timeToNextMigration < minMigration) {
                            		minMigration = timeToNextMigration;
                            		stateIdMigrationFrom = i;
                            		stateIdMigrationTo = j;
                            	}
                    		}
                      	}
                	}else {
                      	for (int j=i+1; j<uniqueStates.size(); j++) {
                    		if (i!=j) {
                    			double timeToNextMigration = Randomizer.nextExponential(k_*migrationRates
                    					.getArrayValue(i*uniqueStates.size()+j));
                            	if (timeToNextMigration < minMigration) {
                            		minMigration = timeToNextMigration;
                            		stateIdMigrationFrom = i;
                            		stateIdMigrationTo = j;
                            	}
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
                    coalesce(currentTime, extantLineages, stateIdCoal);
                else if (timeUntilNextEvent == minReassort)
                    reassort(currentTime, extantLineages, stateIdReassortment);
                else
                	migrate(currentTime, extantLineages, stateIdMigrationFrom, stateIdMigrationTo);
            } else {
                currentTime += timeUntilNextSample;
                sample(remainingSampleNodes, extantLineages);
            }
        
            allEmpty = extantLineages.stream().allMatch(l -> l == null || l.size() <= 1);
        }
        while (!allEmpty || !remainingSampleNodes.isEmpty());

        final List<List <NetworkEdge>> root = extantLineages.stream().filter(l -> l.size()==1).collect(Collectors.toList());
        if (root.size() > 1 || root.get(0).size() > 1)
        	System.err.println("More than one root edge");
        setRootEdge(root.get(0).get(0));
    }
    
    
    private void sample(List<NetworkNode> remainingSampleNodes, List<List<NetworkEdge>> extantLineages) {
        // sample the network node
        NetworkNode n = remainingSampleNodes.get(0);

        // Create corresponding lineage
        BitSet hasSegs = new BitSet();
        hasSegs.set(0, nSegments);
        NetworkEdge lineage = new NetworkEdge(null, n, hasSegs);
        int id = n.getStateIndex();
        extantLineages.get(id).add(lineage);
        n.addParentEdge(lineage);

        remainingSampleNodes.remove(0);
    }
    
    private void coalesce(double coalescentTime, List<List<NetworkEdge>> extantLineages, int stateIdCoal) {
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
    
    private void reassort(double reassortmentTime, List<List<NetworkEdge>> extantLineages, int stateIdReassortment) {
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
        node.setStateIndex(lineage.childNode.getStateIndex());
        node.setStateLabel(lineage.childNode.getStateLabel());

        // Create reassortment lineages
        NetworkEdge leftLineage = new NetworkEdge(null, node, hasSegs_left);
        NetworkEdge rightLineage = new NetworkEdge(null, node, hasSegs_right);
        node.addParentEdge(leftLineage);
        node.addParentEdge(rightLineage);

        extantLineages.get(stateIdReassortment).remove(lineage);
        extantLineages.get(stateIdReassortment).add(leftLineage);
        extantLineages.get(stateIdReassortment).add(rightLineage);
    }
    
    private void migrate(double migrationTime, List<List<NetworkEdge>> extantLineages, int stateIdMigrationFrom, int stateIdMigrationTo) {
        // Sample the lineage for migration:
        NetworkEdge lineage = extantLineages.get(stateIdMigrationFrom).
        		get(Randomizer.nextInt(extantLineages.get(stateIdMigrationFrom).size()));

        NetworkNode migrationPoint = new NetworkNode();
        NetworkNode newParentNode = lineage.parentNode;
        NetworkEdge newParentEdge = lineage.getCopy();
        
        migrationPoint.setHeight(migrationTime);
        migrationPoint.addParentEdge(newParentEdge);
        newParentNode.addChildEdge(newParentEdge);
        
        newParentNode.removeChildEdge(lineage);
        migrationPoint.addChildEdge(lineage);
        
        migrationPoint.setStateIndex(stateIdMigrationTo);
        migrationPoint.setStateLabel(uniqueStates.get(stateIdMigrationTo));
        
        extantLineages.get(stateIdMigrationFrom).remove(lineage);
        extantLineages.get(stateIdMigrationTo).add(newParentEdge);
    }

}
