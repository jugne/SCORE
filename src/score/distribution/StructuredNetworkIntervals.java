package score.distribution;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.util.Precision;

import beast.base.inference.CalculationNode;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import coalre.network.Network;

/**
 * Extended for structured case by Ugne Jankauskaite
 * @author Nicola Felix Mueller
 */
public class StructuredNetworkIntervals extends CalculationNode {
    public Input<Network> networkInput = new Input<>("network",
            "network for which to calculate the intervals", Validate.REQUIRED);

    public Input<Function> binomialProbInput = new Input<>("binomialProb",
            "Probability of a given segment choosing a particular parent.");

    private Network network;

    private List<StructuredNetworkEvent> networkEventList, storedNetworkEventList;

    public boolean eventListDirty = true;

    @Override
    public void initAndValidate() {
        network = networkInput.get();

        storedNetworkEventList = new ArrayList<>();
    }

	public void initAndValidate(Network network) {
		this.network = network;

		storedNetworkEventList = new ArrayList<>();
	}

    public List<StructuredNetworkEvent> getNetworkEventList() {
		network = networkInput.get();
		update();

		return networkEventList;
	}

	public List<StructuredNetworkEvent> getNetworkEventList(Network network) {
		this.network = network;
        update();

        return networkEventList;
    }

    public double getBinomialProb() {
        return binomialProbInput.get() != null
                ? binomialProbInput.get().getArrayValue()
                : 0.5;
    }

    void update() {
//        if (!eventListDirty)
//            return;

//		System.out.println(network.getExtendedNewick());
        networkEventList = network.getNodes().stream().map(n -> {
            StructuredNetworkEvent event = new StructuredNetworkEvent();
			event.time = Precision.round(n.getHeight(), 10);
//			event.time = n.getHeight();
            event.node = n;
            switch(n.getChildCount()) {
                case 0:
                    event.type = StructuredNetworkEvent.NetworkEventType.SAMPLE;
                    break;

                case 1:
                    event.type = StructuredNetworkEvent.NetworkEventType.REASSORTMENT;
                    break;

                case 2:
                    event.type = StructuredNetworkEvent.NetworkEventType.COALESCENCE;
                    break;

                default:
                    throw new RuntimeException("Network node has illegal number of children.");
            }
            return event;
        }).sorted(Comparator.comparingDouble(e -> e.time)).collect(Collectors.toList());

        int lineages = 0;
        double totalReassortmentObsProb = 0;

        for (StructuredNetworkEvent event : networkEventList) {
            switch(event.type) {
                case SAMPLE:
                    lineages += 1;
                    event.lineagesAdded.add(event.node.getParentEdges().get(0));
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    break;

                case REASSORTMENT:
                    lineages += 1;
                    
                    event.lineagesAdded.add(event.node.getParentEdges().get(0));
                    event.lineagesAdded.add(event.node.getParentEdges().get(1));
                    event.lineagesRemoved.add(event.node.getChildEdges().get(0));
                    
                    totalReassortmentObsProb -= event.node.getChildEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(1).getReassortmentObsProb(getBinomialProb());

                    event.segsToSort = event.node.getChildEdges().get(0).hasSegments.cardinality();
                    event.segsSortedLeft = event.node.getParentEdges().get(0).hasSegments.cardinality();
                    break;

                case COALESCENCE:
                    lineages -= 1;
                    
                    event.lineagesAdded.add(event.node.getParentEdges().get(0));
                    event.lineagesRemoved.add(event.node.getChildEdges().get(0));
                    event.lineagesRemoved.add(event.node.getChildEdges().get(1));
                    
                    totalReassortmentObsProb -= event.node.getChildEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb -= event.node.getChildEdges().get(1).getReassortmentObsProb(getBinomialProb());
                    totalReassortmentObsProb += event.node.getParentEdges().get(0).getReassortmentObsProb(getBinomialProb());
                    break;
            }

            event.lineages = lineages;
            event.totalReassortmentObsProb = totalReassortmentObsProb;
        }

//        eventListDirty = false;
    }

    @Override
    protected boolean requiresRecalculation() {
        eventListDirty = true;

        return true;
    }

    @Override
    protected void restore() {
        List<StructuredNetworkEvent> tmp = networkEventList;
        networkEventList = storedNetworkEventList;
        storedNetworkEventList = tmp;
        
        eventListDirty = true;
        
        super.restore();
    }

    @Override
    protected void store() {
        storedNetworkEventList.clear();
        storedNetworkEventList.addAll(networkEventList);

        super.store();
    }
}