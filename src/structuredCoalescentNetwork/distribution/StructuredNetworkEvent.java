package structuredCoalescentNetwork.distribution;

import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

public class StructuredNetworkEvent {
    public enum NetworkEventType { SAMPLE, COALESCENCE, REASSORTMENT }

    public NetworkEventType type;
    public double time;

    /**
     * Number of segments on a reassorting lineage.
     */
    public int segsToSort;

    /**
     * Number of segments sent to the first parent.
     */
    public int segsSortedLeft;

    public int lineages;
    public List<NetworkEdge> lineagesAdded = new ArrayList<NetworkEdge>();
    public List<NetworkEdge> lineagesRemoved = new ArrayList<NetworkEdge>();;
    
    
    public double totalReassortmentObsProb;

    /**
     * Only used when setting up event list.
     * May not point to a compatible node at other times.
     */
    public NetworkNode node;
}
