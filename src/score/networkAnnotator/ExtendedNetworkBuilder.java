package score.networkAnnotator;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.antlr.v4.runtime.CharStream;
import org.antlr.v4.runtime.CharStreams;
import org.antlr.v4.runtime.CommonTokenStream;
import org.antlr.v4.runtime.tree.ParseTree;

import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;
import coalre.network.parser.NetworkBaseVisitor;
import coalre.network.parser.NetworkLexer;
import coalre.network.parser.NetworkParser;

/**
 * Code copied from CoalRe package by N.F. Muller et al Adjustment for
 * structured network by Ugne Stolz <ugne.stolz@protonmail.com>
 * 
 * @date 12 May 2020
 */
public class ExtendedNetworkBuilder extends Network {

	@Override
	public void fromExtendedNewick(String newickStr) {

		CharStream inputStream = CharStreams.fromString(newickStr);
		NetworkLexer lexer = new NetworkLexer(inputStream);
		CommonTokenStream tokenStream = new CommonTokenStream(lexer);
		NetworkParser parser = new NetworkParser(tokenStream);
		ParseTree tree = parser.network();

		StructuredNetworkBuilderVisitor builder = new StructuredNetworkBuilderVisitor();
		rootEdge = builder.visit(tree);

		List<NetworkNode> leafNodes = new ArrayList<>(getLeafNodes());
	}

	class StructuredNetworkBuilderVisitor extends NetworkBaseVisitor<NetworkEdge> {
		Map<Integer, NetworkNode> seenHybrids;
		Map<NetworkEdge, Double> edgeLengths;



		double getMaxRootToLeafTime(NetworkNode node, Set<NetworkNode> seenNodes) {

			if (seenNodes.contains(node))
				return 0.0;

			seenNodes.add(node);

			double maxTime = 0.0;
			for (NetworkEdge childEdge : node.getChildEdges()) {
				NetworkNode childNode = childEdge.childNode;
				childNode.setHeight(node.getHeight() - edgeLengths.get(childEdge));

				double thisTime = edgeLengths.get(childEdge) +
						getMaxRootToLeafTime(childNode, seenNodes);
				if (thisTime > maxTime)
					maxTime = thisTime;
			}

			return maxTime;
		}

		void shiftNodeHeights(double maxRTLT, NetworkNode node, Set<NetworkNode> seenNodes) {
			if (seenNodes.contains(node))
				return;

			seenNodes.add(node);

			node.setHeight(node.getHeight() + maxRTLT);

			for (NetworkEdge childEdge : node.getChildEdges())
				shiftNodeHeights(maxRTLT, childEdge.childNode, seenNodes);
		}

		@Override
		public NetworkEdge visitNetwork(NetworkParser.NetworkContext ctx) {
			seenHybrids = new HashMap<>();
			edgeLengths = new HashMap<>();

			NetworkEdge rootEdge = visit(ctx.node());

			Set<NetworkNode> seenNodes = new HashSet<>();
			NetworkNode rootNode = rootEdge.childNode;
			rootNode.setHeight(0.0);
			double maxRTLT = getMaxRootToLeafTime(rootNode, seenNodes);

			seenNodes.clear();
			shiftNodeHeights(maxRTLT, rootEdge.childNode, seenNodes);

			return rootEdge;
		}

		private String removeQuotes(String str) {

			String[] quoteChars = { "\"", "'" };

			for (String quoteChar : quoteChars) {
				if (str.startsWith(quoteChar) && str.endsWith(quoteChar) && str.length() >= 2)
					str = str.substring(1, str.length() - 1);
			}

			return str;
		}

		@Override
		public NetworkEdge visitNode(NetworkParser.NodeContext ctx) {

			visit(ctx.post());

			NetworkNode node;

			if (ctx.post().hybrid() != null) {
				int hybridID = Integer.valueOf(ctx.post().hybrid().id.getText());

				if (seenHybrids.containsKey(hybridID)) {
					node = seenHybrids.get(hybridID);
				} else {
					node = new NetworkNode();
					seenHybrids.put(hybridID, node);
				}
			} else {
				node = new NetworkNode();
			}

			if (ctx.post().label() != null)
				node.setTaxonLabel(removeQuotes(ctx.post().label().getText()));

			for (NetworkParser.NodeContext childNodeCtx : ctx.node()) {
				NetworkEdge childEdge = visit(childNodeCtx);
				childEdge.parentNode = node;
				node.addChildEdge(childEdge);
			}

			boolean segmentsProcessed = false;
			BitSet hasSegments = new BitSet();
			if (ctx.post().meta() != null
					&& ctx.post().meta().attrib() != null) {

				for (NetworkParser.AttribContext attribCtx : ctx.post().meta().attrib()) {
					if (removeQuotes(attribCtx.attribKey.getText()).equals("state"))
						node.setTypeLabel(attribCtx.attribValue().getText());

					if (!removeQuotes(attribCtx.attribKey.getText()).equals("segments"))
						continue;

					if (attribCtx.attribValue().vector() == null)
						continue;

					for (NetworkParser.AttribValueContext attribValueCtx : attribCtx.attribValue().vector()
							.attribValue())
						hasSegments.set(Integer.valueOf(attribValueCtx.getText()));

					segmentsProcessed = true;
				}

			}

			if (!segmentsProcessed) {
				throw new RuntimeException("Segment attribute missing/malformed " +
						"for edge in input network string.");
			}

			NetworkEdge edge = new NetworkEdge(null, node, hasSegments);
			node.addParentEdge(edge);

			if (ctx.post().length == null) {
				throw new RuntimeException("Edge missing length in input " +
						"network string.");
			}

			edgeLengths.put(edge, Double.valueOf(ctx.post().length.getText()));

			return edge;
		}
	}

}
