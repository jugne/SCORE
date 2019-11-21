package structuredCoalescentNetwork.distribution;

import java.util.ArrayList;
import java.util.List;

import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;

public class StructuredNetworkEvent {
    public enum NetworkEventType {
	SAMPLE, COALESCENCE, REASSORTMENT
    }

    public NetworkEventType type;
    public double time;

	// Needed for stochastic mapping
	public double[][] p_stored;
	public double[][] pDot_stored;
	public double[][] pDotDot_stored;
	public double[] intermediateTimeStored;
	public int numRecords;
	public ArrayList<NetworkEdge> activeLineages;

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
     * Only used when setting up event list. May not point to a compatible node at
     * other times.
     */
    public NetworkNode node;
}
