package score.logger;

import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import score.dynamics.ConstantReassortment;
import score.mapping.MappedNetwork;
import score.simulator.SimulateStructureCoalescentNetwork;

/**
 * Logger for generating statistics from type mapped trees.
 */
public class TypedNetworkStatsLogger extends BEASTObject implements Loggable {

	public Input<MappedNetwork> typedNetworkInput = new Input<>("network",
			"Network with type changes mapped.");

	public Input<SimulateStructureCoalescentNetwork> simulatedTypedNetworkInput = new Input<>("simulatedNetwork",
			"Simulated network with type changes mapped.",
			Input.Validate.XOR, typedNetworkInput);


    int[][] countMatrix;

	MappedNetwork network;
	SimulateStructureCoalescentNetwork simNetwork;
	ConstantReassortment dynamics;
    int nTypes;
    String typeLabel;
	boolean simulation;

    @Override
    public void initAndValidate() {
		if (simulatedTypedNetworkInput.get() != null) {
			simulation = true;
			simNetwork = simulatedTypedNetworkInput.get();
			nTypes = simNetwork.typeIndexToName.keySet().size();
		} else {
			simulation = false;
			network = typedNetworkInput.get();
			dynamics = network.dynamics;
			nTypes = dynamics.getNrTypes();
		}
        countMatrix = new int[nTypes][nTypes];
    }


	private void countChanges() {
		for (int i = 0; i < nTypes; i++) {
			for (int j = 0; j < nTypes; j++) {
				countMatrix[i][j] = 0;
			}
		}
		List<NetworkNode> migrationNodes;

		if (simulation)
			migrationNodes = simNetwork.getNodes().stream().filter(n -> n.getParentCount() == 1)
				.filter(n -> n.getChildCount() == 1).sorted(Comparator.comparing(NetworkNode::getHeight))
				.collect(Collectors.toList());
		else
			migrationNodes = network.getNodes().stream().filter(n -> n.getParentCount() == 1)
					.filter(n -> n.getChildCount() == 1).sorted(Comparator.comparing(NetworkNode::getHeight))
					.collect(Collectors.toList());

		Collections.reverse(migrationNodes);

		for (NetworkNode m : migrationNodes) {
			NetworkEdge childEdge = m.getChildEdges().get(0);

			int fromType = m.getTypeIndex();
			int toType = childEdge.childNode.getTypeIndex();
			if (fromType == toType)
			{
				System.out.println("log");
				System.out.println(network.getExtendedNewick());
				throw new IllegalArgumentException("Not valid migration node");
			}

			countMatrix[fromType][toType] += 1;
		}
    }


    @Override
    public void init(PrintStream out) {

		String prefix;
		if (simulation)
			prefix = simNetwork.getID() != null
					? simNetwork.getID() + "."
					: "";
		else
			prefix = network.getID() != null
					? network.getID() + "."
					: "";

        for (int type=0; type<nTypes; type++) {
            for (int typeP=0; typeP<nTypes; typeP++) {
                if (type == typeP)
                    continue;

				if (simulation)
					out.print(prefix + "count_" + simNetwork.typeIndexToName.get(type)
							+ "_to_" + simNetwork.typeIndexToName.get(typeP) + "\t");

				else
					out.print(prefix + "count_" + dynamics.getStringStateValue(type)
						+ "_to_" + dynamics.getStringStateValue(typeP) + "\t");
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {

		if (!simulation) {
			network = (MappedNetwork) network.getCurrent();
			network.remapForLog(sample);
		}
		else
			simNetwork = (SimulateStructureCoalescentNetwork) simNetwork.getCurrent();

		countChanges();

        for (int type=0; type<nTypes; type++) {
            for (int typeP = 0; typeP < nTypes; typeP++) {
                if (type == typeP)
                    continue;

                out.print(countMatrix[type][typeP] + "\t");
            }
        }
    }

    @Override
    public void close(PrintStream out) { }
}
