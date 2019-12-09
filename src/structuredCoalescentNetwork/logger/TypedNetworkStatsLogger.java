package structuredCoalescentNetwork.logger;

import java.io.PrintStream;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import structuredCoalescentNetwork.dynamics.ConstantReassortment;
import structuredCoalescentNetwork.mapping.MappedNetwork;

/**
 * Logger for generating statistics from type mapped trees.
 */
public class TypedNetworkStatsLogger extends BEASTObject implements Loggable {

	public Input<MappedNetwork> typedNetworkInput = new Input<>("typedNetwork",
            "Tree with type changes mapped.",
            Input.Validate.REQUIRED);

	public Input<ConstantReassortment> dynamicsInput = new Input<>("dynamics", "Input of rates",
			Input.Validate.REQUIRED);


    int[][] countMatrix;

	MappedNetwork network;
	ConstantReassortment dynamics;
    int nTypes;
    String typeLabel;
    boolean includeRootEdgeChanges;

    @Override
    public void initAndValidate() {
		network = typedNetworkInput.get();
		dynamics = dynamicsInput.get();
		nTypes = dynamics.getNrTypes();

        countMatrix = new int[nTypes][nTypes];
    }


	private void countChanges() {
		for (int i = 0; i < nTypes; i++) {
			for (int j = 0; j < nTypes; j++) {
				countMatrix[i][j] = 0;
			}
		}

		List<NetworkNode> migrationNodes = network.getNodes().stream().filter(n -> n.getParentCount() == 1)
				.filter(n -> n.getChildCount() == 1).sorted(Comparator.comparing(NetworkNode::getHeight))
				.collect(Collectors.toList());

		Collections.reverse(migrationNodes);

		for (NetworkNode m : migrationNodes) {
			NetworkEdge parentEdge = m.getParentEdges().get(0);

			int fromType = parentEdge.parentNode.getTypeIndex();
			int toType = m.getTypeIndex();
			if (fromType == toType)
				throw new IllegalArgumentException("Not valid migration node");
			countMatrix[fromType][toType] += 1;
		}
    }


    @Override
    public void init(PrintStream out) {

		String prefix = network.getID() != null
				? network.getID() + "."
                : "";

        for (int type=0; type<nTypes; type++) {
            for (int typeP=0; typeP<nTypes; typeP++) {
                if (type == typeP)
                    continue;

				out.print(prefix + "count_" + dynamics.getStringStateValue(type)
						+ "_to_" + dynamics.getStringStateValue(typeP) + "\t");
            }
        }
    }

    @Override
    public void log(long sample, PrintStream out) {

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
