/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package score.networkAnnotator;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.swing.BoxLayout;
import javax.swing.GroupLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.WindowConstants;
import javax.swing.border.EtchedBorder;

import beast.core.util.Log;
import coalre.network.Network;
import coalre.network.NetworkEdge;
import coalre.network.NetworkNode;


/**
 * Number of reassortment events preceding migration events on fit and unfit
 * network edges.
 * 
 * 
 * @author Ugne Stolz <ugne.stolz@protonmail.com>
 * @date 20 Apr 2020
 */
public class ReassortmentAndMigration extends SCoReAnnotator {

	List<NetworkNode> allFitNodes;

	List<NetworkNode> allUnfitNodes;
	List<Double> leaveDistance;
    
	List<NetworkNode> reaAfterList;
	List<NetworkNode> reaBeforeList;
	List<NetworkNode> migrationForwardList;

	List<NetworkNode> migrationBackwardList;
	List<NetworkNode> migrationUnasignedList;

	List<Double> reaToTipLengths;

	List<Double> migToTipLengths;
	Double lengthForward;
	Double lengthBackward;
	private boolean firstNet;
	Set<String> typeKeys;

    private static class NetworkAnnotatorOptions {
        File inFile;
		File outFile = new File("reassortment_and_migration.txt");
        double burninPercentage = 10.0;

		double maxMigrationDistance;

		double minTipDistance;
        int[] removeSegments = new int[0];


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
					"Max distance between reassortment and migration event " + maxMigrationDistance + "%\n" +
					"minimal distance to a tip to be considered trunk\n" +
					"(ignored in MostRecentSample Case): " + minTipDistance;
        }
    }

	public ReassortmentAndMigration(NetworkAnnotatorOptions options) throws Exception {

        // Display options:
        System.out.println(options + "\n");

        // Initialise reader
		StructuredReassortmentLogReader logReader = new StructuredReassortmentLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

	      System.out.println("\nWriting output to " + options.outFile.getName()
	      + "...");
        
		// compute number of migration within maxMigrationDistance to reassotment
        try (PrintStream ps = new PrintStream(options.outFile)) {

			ps.print("n_reassortment" + "\t" +
					"n_migration" + "\t" +
					"n_migration_fit" + "\t" +
					"n_migration_unfit" + "\t" +
					"n_reassortment_total_fit" + "\t" +
					"n_reassortment_total_unfit" + "\t" +
					"n_reassortment_window" + "\t" +
					"n_reassortment_window_fit" + "\t" +
					"n_reassortment_window_unfit" + "\t" +
					"n_reassortment_off_window" + "\t" +
					"n_reassortment_off_window_fit" + "\t" +
					"n_reassortment_off_window_unfit" + "\t" +
					"network_length" + "\t" + 
					"network_length_fit" + "\t" +
					"network_length_unfit" + "\t" +
					"network_length_window" + "\t" +
					"network_length_window_fit" + "\t" +
					"network_length_window_unfit" + "\t" +
					"network_length_off_window_fit" + "\t" +
					"network_length_off_window_unfit");

			// + "\t" + "mig_tip_lengths" + "\t"
//					+ "mig_tip_heights");
//			ps.print("\n");
			firstNet = true;
			for (Network network : logReader) {
	        	pruneNetwork(network, options.removeSegments);
				computeTrunkReassortmentLeaveDist(network, options.minTipDistance);
				computeReassortmentAndMigration(network, ps, options.maxMigrationDistance);
	        	
	        	ps.print("\n");
				firstNet = false;
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }

	/**
	 * gets how many reticulation events happen on the trunk vs. not on the trunk
	 * The trunk is define as any edge on the network that has descendents that are
	 * more than minTipDistance away from that node
	 * 
	 * @param network
	 * @param ps
	 * @param minTipDistance
	 */
	private void computeTrunkReassortmentLeaveDist(Network network, double minTipDistance) {

		// get the length of the network
//		List<NetworkEdge> allEdges = network.getEdges().stream()
//				.filter(e -> !e.isRootEdge())
//				.collect(Collectors.toList());

//		double fullLength = 0.0;
//		for (NetworkEdge edge : allEdges)
//			fullLength += edge.getLength();

		// compute which nodes are on the trunk of the network defined as any connection
		// between
		// the root and the most recent sampled individual. To do so, get all zero
		// height edges
		List<NetworkEdge> zeroHeightEdges = network.getEdges().stream()
				.filter(e -> e.isLeafEdge())
				.collect(Collectors.toList());

		allFitNodes = new ArrayList<>();
		allUnfitNodes = new ArrayList<>();
		leaveDistance = new ArrayList<>();

		if (zeroHeightEdges.size() == 0)
			throw new IllegalArgumentException("no leaf node with 0 height found");

		for (NetworkEdge zeroEdge : zeroHeightEdges) {
			double dist = 0.0;
			getAllAncestralEdgesLeaveDist(zeroEdge.childNode, dist, minTipDistance);
		}

		// drop every node in allTrunkNodes that is less far away from a leave than some
		// threshold
		for (int i = allFitNodes.size() - 1; i >= 0; i--) {
			if (leaveDistance.get(i) < minTipDistance) {
				allUnfitNodes.add(allFitNodes.get(i));
				allFitNodes.remove(i);
			}
		}

		if (network.getNodes().size() != (allUnfitNodes.size() + allFitNodes.size()))
			System.out.println("False");

		// calculate how many reassortment events are on the trunk and how many aren't
//		int onTrunk = allFitNodes.stream()
//				.filter(e -> e.isReassortment())
//				.collect(Collectors.toList()).size();
//
//		int offTrunk = allEdges.stream()
//				.filter(e -> !e.isRootEdge())
//				.filter(e -> e.parentNode.isReassortment())
//				.collect(Collectors.toList()).size() - onTrunk;
//
//		// calculate the length of the trunk
//		double trunkLength = 0.0;
//		for (NetworkNode node : allFitNodes) {
//			for (NetworkEdge edge : node.getParentEdges())
//				if (!edge.isRootEdge())
//					trunkLength += edge.getLength();
//		}

//		ps.print(onTrunk + "\t" + offTrunk + "\t" + trunkLength + "\t" + (fullLength - trunkLength));

	}
        

	private void computeReassortmentAndMigration(Network network, PrintStream ps, double maxMigrationDistance)
			throws Exception {
		// get all reassortment nodes of the network

		List<NetworkNode> reaNodes = network.getNodes().stream()
				.filter(n -> n.isReassortment())
				.collect(Collectors.toList());
		
		List<NetworkNode> migNodes = network.getNodes().stream()
				.filter(n -> isMigrationNode(n))
				.sorted(Comparator.comparingDouble(n -> n.getHeight()))
				.collect(Collectors.toList());
		


		// clear the lists and sums
		migrationForwardList = new ArrayList<>();
		migrationBackwardList = new ArrayList<>();
		migrationUnasignedList = new ArrayList<>();
		reaToTipLengths = new ArrayList<>();
		migToTipLengths = new ArrayList<>();
		reaAfterList = new ArrayList<NetworkNode>();
		reaBeforeList = new ArrayList<NetworkNode>();
		lengthBackward = 0.0;
		lengthForward = 0.0;
		List<Double> migToTipHeight = new ArrayList<Double>();
		List<NetworkNode> reaCounts = new ArrayList<NetworkNode>();
		Double length = 0.0;

		// Map NetworkEdge -> length. Here length is taken up to the point
		// where maxMigrationDistance (measured in network height, not length)
		// to the migration node is reached
		// (which can also be its full length).
		// It is used to calculate the total length of the migration node parent
		// network, within the maxMigrationDistance window
		HashMap<NetworkEdge, Double> edgeLengthMap = new HashMap<NetworkEdge, Double>();


		for (NetworkNode node : migNodes) {
			reaCounts.addAll(reaNodesBefore(node, 0.0, maxMigrationDistance, length, edgeLengthMap));
			migToTipLengths.add(lengthToTip(node));
			migToTipHeight.add(node.getHeight() - closesDescendantTip(node).getHeight());
		}
		
		// Colect reassortment nodes within the windows above migration nodes.
		// Create a set to avoid double counting.
		Set<NetworkNode> reaNodesInWindowSet = new HashSet<NetworkNode>();
		for (NetworkEdge edge : edgeLengthMap.keySet()) {
			if (edge.childNode.isReassortment())
				reaNodesInWindowSet.add(edge.childNode);
		}
		
		// Total edge length within the windows above migration nodes
		double totalLengthWindow = 0.0;
		for (double l : edgeLengthMap.values())
			totalLengthWindow += l;

		// get the length of the network
		List<NetworkEdge> allEdges = network.getEdges().stream()
				.filter(e -> !e.isRootEdge())
				.collect(Collectors.toList());

		HashMap<String, Double> typesToLength = new HashMap<>();
		double fullLength = 0.0;
		for (NetworkEdge edge : allEdges) {
			fullLength += edge.getLength();
			typesToLength.merge(edge.childNode.getTypeLabel(), edge.getLength(), Double::sum);
		}

		// Fit and unfit edges separation

		// On fit edges
		List<NetworkNode> reaNodesFit = new ArrayList<NetworkNode>();
		List<NetworkNode> reaNodesFitInWindow = new ArrayList<NetworkNode>(); // on window
		List<NetworkNode> reaNodesFitOffWindow = new ArrayList<NetworkNode>(); // off window
		List<NetworkNode> migNodesFit = new ArrayList<NetworkNode>();
		// Off fit edges
		List<NetworkNode> reaNodesUnfit = new ArrayList<NetworkNode>();
		List<NetworkNode> reaNodesUnfitInWindow = new ArrayList<NetworkNode>();
		List<NetworkNode> reaNodesUnfitOffWindow = new ArrayList<NetworkNode>();
		List<NetworkNode> migNodesUnfit = new ArrayList<NetworkNode>();
		
		// for loops should speed up the multiple streams on the same list below
		for (NetworkNode n : allFitNodes) {
			if (n.isReassortment()) {
				reaNodesFit.add(n);
			} else if(isMigrationNode(n)) {
				migNodesFit.add(n);
			}
		}
		
		for (NetworkNode n : reaNodesFit) {
			if (reaNodesInWindowSet.contains(n)) {
				reaNodesFitInWindow.add(n);
			} else {
				reaNodesFitOffWindow.add(n);
			}
		}

		for (NetworkNode n : allUnfitNodes) {
			if (n.isReassortment()) {
				reaNodesUnfit.add(n);
			} else if (isMigrationNode(n)) {
				migNodesUnfit.add(n);
			}
		}

		for (NetworkNode n : reaNodesUnfit) {
			if (reaNodesInWindowSet.contains(n)) {
				reaNodesUnfitInWindow.add(n);
			} else {
				reaNodesUnfitOffWindow.add(n);
			}
		}


//		reaNodesFit = allFitNodes.stream()
//				.filter(e -> e.isReassortment())
//				.collect(Collectors.toList());
		
		
//		reaNodesFitInWindow = reaNodesFit.stream()
//				.filter(n -> reaNodesInWindowSet.contains(n))
//				.collect(Collectors.toList());
//		reaNodesFitOffWindow = reaNodesFit.stream()
//				.filter(n -> !reaNodesInWindowSet.contains(n))
//				.collect(Collectors.toList());

//		migNodesFit = allFitNodes.stream()
//				.filter(e -> isMigrationNode(e))
//				.collect(Collectors.toList());

//		reaNodesUnfit = allUnfitNodes.stream()
//				.filter(e -> e.isReassortment())
//				.collect(Collectors.toList());
//		reaNodesUnfitInWindow = reaNodesUnfit.stream()
//				.filter(n -> reaNodesInWindowSet.contains(n))
//				.collect(Collectors.toList());
//		reaNodesUnfitOffWindow = reaNodesUnfit.stream()
//				.filter(n -> !reaNodesInWindowSet.contains(n))
//				.collect(Collectors.toList());
		
//		migNodesUnfit = allUnfitNodes.stream()
//				.filter(e -> isMigrationNode(e))
//				.collect(Collectors.toList());

		// calculate the length of the trunk
		double fitOnWindowLength = 0.0;
		double unfitOnWindowLength = 0.0;
		for (NetworkEdge edge : edgeLengthMap.keySet()) {
			if (allFitNodes.contains(edge.childNode))
				fitOnWindowLength += edgeLengthMap.get(edge);
			else
				unfitOnWindowLength += edgeLengthMap.get(edge);
		}

		

		double fitOffWindowLength = 0.0;
		double unfitOffWindowLength = 0.0;
		double totalFitLength = 0.0;
		double totalUnfitLength = 0.0;
		for (NetworkNode node : allFitNodes) {
			for (NetworkEdge edge : node.getParentEdges())
				if (!edge.isRootEdge())
					totalFitLength += edge.getLength();
		}
		
		fitOffWindowLength = totalFitLength - fitOnWindowLength;

		totalUnfitLength = fullLength - totalFitLength;
		unfitOffWindowLength = totalUnfitLength - unfitOnWindowLength;
		
		if (firstNet) {
			ps.print("\t");
			for (String key : typesToLength.keySet()) {
				ps.print("network_length_" + key + "\t");
			}
			typeKeys = typesToLength.keySet();
		}

		ps.print("\n");
		String s = reaNodes.size() + "\t" +
				migNodes.size() + "\t" +
				migNodesFit.size() + "\t" +
				migNodesUnfit.size() + "\t" +
				reaNodesFit.size() + "\t" +
				reaNodesUnfit.size() + "\t" +
				reaNodesInWindowSet.size() + "\t" +
				reaNodesFitInWindow.size() + "\t" +
				reaNodesUnfitInWindow.size() + "\t" +
				(reaNodes.size() - reaNodesInWindowSet.size()) + "\t" +
				reaNodesFitOffWindow.size() + "\t" +
				reaNodesUnfitOffWindow.size() + "\t" +
				fullLength + "\t" +
				totalFitLength + "\t" +
				totalUnfitLength + "\t" +
				totalLengthWindow + "\t" +
				fitOnWindowLength + "\t" +
				unfitOnWindowLength + "\t" +
				fitOffWindowLength + "\t" +
				unfitOffWindowLength + "\t";

		int l = 1;
		for (String key : typeKeys) {
			s += typesToLength.get(key);
			if (l < typeKeys.size())
				s += "\t";
			l++;
		}

		ps.print(s);

		// + "\t" + Arrays.toString(migToTipLengths.toArray()) + "\t"
//				+ Arrays.toString(migToTipHeight.toArray()));
	}

	private NetworkNode findDownStreamReas(NetworkNode node) {
		NetworkNode left = null;
		NetworkNode right = null;
		
		if (!node.isLeaf()) {
			left = node.getChildEdges().get(0).childNode;
			if (!left.isReassortment()) {
				left = findDownStreamReas(left);
			}

			// check if right sibling exists
			if (node.getChildCount() > 1) {
				right = node.getChildEdges().get(1).childNode;
				if (!right.isReassortment()) {
					right = findDownStreamReas(right);
				}
			}
			if (right == null) {
				if (left == null)
					return null;
				else
					return left;
			} else if (left == null)
				return right;
			else if (left.getHeight() > right.getHeight())
				return left;
			else
				return right;
		}
		return null;
	}

	private NetworkNode findUpStreamReas(NetworkNode node) {
		NetworkNode left = null;
		NetworkNode right = null;

		if (!node.getParentEdges().get(0).isRootEdge()) {
			left = node.getParentEdges().get(0).parentNode;
			if (!left.isReassortment()) {
				left = findUpStreamReas(left);
			}

			// check if right sibling exists
			if (node.getParentCount() > 1) {
				right = node.getParentEdges().get(0).parentNode;
				if (!right.isReassortment()) {
					right = findUpStreamReas(right);
				}
			}
			if (right == null) {
				if (left == null)
					return null;
				else
					return left;
			} else if (left == null)
				return right;
			else if (left.getHeight() < right.getHeight())
				return left;
			else
				return right;
		}
		return null;
	}

	private double findMedian(List<Double> list) {
		int size = list.size();
		list = list.stream().sorted().collect(Collectors.toList());

		if ((size - 1) % 2 == 0)
			return (list.get((size - 1) / 2) + list.get((size - 1) / 2 - 1)) / 2.0;

		return list.get((size - 1) / 2);

	}

	private double findMin(List<Double> list) {
		list = list.stream().sorted().collect(Collectors.toList());

		return list.get(0);

	}

	private List<NetworkNode> getMigrationNodesForward(NetworkNode node, double offset, double threshold) {
		List<NetworkNode> migrationForward = new ArrayList<>();
		if (!node.isLeaf()) {
		NetworkNode left = node.getChildEdges().get(0).childNode;
			if (node.getHeight() - left.getHeight() < threshold - offset) {
			if (isMigrationNode(left))
				migrationForward.add(left);
				migrationForward.addAll(
						getMigrationNodesForward(left, node.getHeight() - left.getHeight() + offset, threshold));
			}
		
		
		// check if right sibling exists
		if (node.getChildCount() > 1) {
			NetworkNode right = node.getChildEdges().get(1).childNode;
				if (node.getHeight() - right.getHeight() < threshold - offset) {
				if (isMigrationNode(right))
					migrationForward.add(right);
					migrationForward.addAll(
							getMigrationNodesForward(right, node.getHeight() - right.getHeight() + offset, threshold));
				}

		}
		
			// no double counting
			LinkedHashSet<NetworkNode> hashSet = new LinkedHashSet<>(migrationForward);
			migrationForward = new ArrayList<>(hashSet);
			return migrationForward.stream()
					.sorted(Comparator.comparingDouble(n -> (node.getHeight() - n.getHeight() + offset))).distinct()
					.collect(Collectors.toList());
		} else {
			return migrationForward;
		}
	}

	private double lengthToTip(NetworkNode node) {
		double leftLength = 0.0;
		double rightLength = 0.0;

		if (!node.isLeaf()) {
			NetworkNode left = node.getChildEdges().get(0).childNode;
			leftLength += node.getHeight() - left.getHeight();
			if (!left.isLeaf())
				leftLength += lengthToTip(left);

			// check if right sibling exists
			if (node.getChildCount() > 1) {
				NetworkNode right = node.getChildEdges().get(1).childNode;
				rightLength += node.getHeight() - right.getHeight();
				if (!right.isLeaf())
					leftLength += lengthToTip(right);
			}
		}
		return leftLength + rightLength;
	}

	private NetworkNode closesDescendantTip(NetworkNode node) {

		if (!node.isLeaf()) {
			NetworkNode left = node.getChildEdges().get(0).childNode;
			if (!left.isLeaf())
				left = closesDescendantTip(left);

			// check if right sibling exists
			if (node.getChildCount() > 1) {
				NetworkNode right = node.getChildEdges().get(1).childNode;
				if (!right.isLeaf())
					right = closesDescendantTip(right);

				if (left.isLeaf() && right.isLeaf())
					return left.getHeight() > right.getHeight() ? left : right;
			}
			return left.isLeaf() ? left : null;
		}
		return null;
	}

	private double lengthAfterRea(NetworkNode reassortment, double offset, double threshold) {
		double leftLength = 0.0;
		double rightLength = 0.0;

		if (!reassortment.isLeaf()) {
			NetworkNode left = reassortment.getChildEdges().get(0).childNode;
			if (reassortment.getHeight() - left.getHeight() < threshold - offset) {
				leftLength += reassortment.getHeight() - left.getHeight();
				if (!left.isReassortment())
					leftLength += lengthAfterRea(left, reassortment.getHeight() - left.getHeight() + offset, threshold);
			} else
				leftLength += threshold - offset;

			// check if right sibling exists
			if (reassortment.getChildCount() > 1) {
				NetworkNode right = reassortment.getChildEdges().get(1).childNode;
				if (reassortment.getHeight() - right.getHeight() < threshold - offset) {
					rightLength += reassortment.getHeight() - right.getHeight();
					if (!right.isReassortment())
						rightLength += lengthAfterRea(right, reassortment.getHeight() - right.getHeight() + offset,
							threshold);
				} else
					rightLength += threshold - offset;
			}
		}
		return leftLength + rightLength;
	}

	private double lengthBeforeRea(NetworkNode reassortment, double offset, double threshold) {
		double leftLength = 0.0;
		double rightLength = 0.0;

		if (!reassortment.getParentEdges().get(0).isRootEdge()) {
			NetworkNode left = reassortment.getParentEdges().get(0).parentNode;
			if (left.getHeight() - reassortment.getHeight() < threshold - offset) {
				leftLength += left.getHeight() - reassortment.getHeight();
				if (!left.isReassortment())
					leftLength += lengthBeforeRea(left, left.getHeight() - reassortment.getHeight() + offset,
							threshold);
			} else
				leftLength += threshold - offset;

			// check if right sibling exists
			if (reassortment.getParentCount() > 1) {
				NetworkNode right = reassortment.getParentEdges().get(1).parentNode;
				if (right.getHeight() - reassortment.getHeight() < threshold - offset) {
					rightLength += right.getHeight() - reassortment.getHeight();
					if (!right.isReassortment())
						rightLength += lengthBeforeRea(right, right.getHeight() - reassortment.getHeight() + offset,
							threshold);
				} else
					rightLength += threshold - offset;
			}
		}
		return leftLength + rightLength;
	}

	private NetworkNode closestReaBeforeMig(NetworkNode node, double offset, double threshold) {
		List<NetworkNode> reaBackward = new ArrayList<>();
		if (!node.getParentEdges().get(0).isRootEdge()) {
			NetworkNode left = node.getParentEdges().get(0).parentNode;
			if (left.getHeight() - node.getHeight() < threshold - offset) {
				if (left.isReassortment())
					reaBackward.add(left);
				else
					reaBackward.add(closestReaBeforeMig(left, left.getHeight() - node.getHeight() + offset, threshold));
			}

			// check if right parent exists
			if (node.getParentCount() > 1) {
				NetworkNode right = node.getParentEdges().get(1).parentNode;
				if (right.getHeight() - node.getHeight() < threshold - offset) {
					if (right.isReassortment())
						reaBackward.add(right);
					else
						reaBackward.add(
								closestReaBeforeMig(right, right.getHeight() - node.getHeight() + offset, threshold));
				}

			}

			// no double counting
			if (!reaBackward.isEmpty())
				return reaBackward.stream()
						.sorted(Comparator.comparingDouble(n -> (n.getHeight() - node.getHeight() + offset))).distinct()
						.collect(Collectors.toList()).get(0);
			else
				return null;

		} else {
			return null;
		}
	}

	private List<NetworkNode> reaNodesBefore(NetworkNode node, double offset, Double threshold, double length,
			HashMap<NetworkEdge, Double> edgeLengthMap)
			throws Exception {
		List<NetworkNode> reaBeforeList = new ArrayList<>();
		if (!node.getParentEdges().get(0).isRootEdge()) {
			NetworkEdge leftParentdEdge = node.getParentEdges().get(0);
			NetworkNode left = node.getParentEdges().get(0).parentNode;
			double leftEdgeLength = leftParentdEdge.getLength();

			if (leftEdgeLength < threshold - offset) {
				length += leftEdgeLength;
				
				if (!edgeLengthMap.containsKey(leftParentdEdge)
						|| (edgeLengthMap.containsKey(leftParentdEdge) && edgeLengthMap.get(leftParentdEdge) < leftEdgeLength))
					edgeLengthMap.put(leftParentdEdge, leftEdgeLength);

				if (isMigrationNode(left)) {
					// reset window from new migration node
//					reaBeforeList.addAll(reaNodesBefore(left, 0.0, threshold, length));
				}
				if (left.isReassortment()) {
					// add the reassortment node
					reaBeforeList.add(left);
					// move further up
					reaBeforeList.addAll(
							reaNodesBefore(left, leftEdgeLength + offset, threshold, length,
									edgeLengthMap));
				} else if (left.isCoalescence()) {
					// move further up
					reaBeforeList.addAll(
							reaNodesBefore(left, leftEdgeLength + offset, threshold, length,
									edgeLengthMap));
				}
			} else {
				if ((threshold - offset) > 0) {
					length += threshold - offset;
					if (!edgeLengthMap.containsKey(leftParentdEdge)
							|| (edgeLengthMap.containsKey(leftParentdEdge)
									&& edgeLengthMap.get(leftParentdEdge) < threshold - offset))
						edgeLengthMap.put(leftParentdEdge, threshold - offset); // add only the edge length up until the
																			// window ends
				}

//					throw new Exception("offset greater than the threshold!!");
//				length += threshold - offset;
			}

				if (node.getParentCount() > 1) {
				NetworkEdge rightParentEdge = node.getParentEdges().get(1);
					NetworkNode right = node.getParentEdges().get(1).parentNode;
				double rightEdgeLength = rightParentEdge.getLength();

				if (rightEdgeLength < threshold - offset) {
					length += rightEdgeLength;

					if (!edgeLengthMap.containsKey(rightParentEdge)
							|| (edgeLengthMap.containsKey(rightParentEdge)
									&& edgeLengthMap.get(rightParentEdge) < rightEdgeLength))
						edgeLengthMap.put(rightParentEdge, rightEdgeLength);

					if (isMigrationNode(right)) {
							// reset window from new migration node
//						reaBeforeList.addAll(reaNodesBefore(right, 0.0, threshold, length));
					}
					if (right.isReassortment()) {
						// add the reassortment node
						reaBeforeList.add(right);
						// move further up
						reaBeforeList.addAll(
								reaNodesBefore(right, rightEdgeLength + offset, threshold,
										length, edgeLengthMap));
					} else if (right.isCoalescence()) {
						// move further up
						reaBeforeList.addAll(
								reaNodesBefore(right, rightEdgeLength + offset, threshold,
										length, edgeLengthMap));
						}
				} else
				if ((threshold - offset) > 0) {
					length += threshold - offset;

					if (!edgeLengthMap.containsKey(rightParentEdge)
							|| (edgeLengthMap.containsKey(rightParentEdge)
									&& edgeLengthMap.get(rightParentEdge) < threshold - offset))
						edgeLengthMap.put(rightParentEdge, threshold - offset);
				}

//					throw new Exception("offset greater than the threshold!!");

					}
				}
		return reaBeforeList;
	}


	private NetworkNode closestReaAfterMig(NetworkNode node, double offset, double threshold) {
		List<NetworkNode> reaForward = new ArrayList<>();
		if (!node.isLeaf()) {
			NetworkNode left = node.getChildEdges().get(0).childNode;
			if (node.getHeight() - left.getHeight() < threshold - offset) {
				if (left.isReassortment())
					reaForward.add(left);
				else {
					NetworkNode temp = closestReaAfterMig(left, node.getHeight() - left.getHeight() + offset,
							threshold);
					if (temp != null)
						reaForward.add(temp);
				}
//					reaForward.add(
//							closestReaAfterMig(left, node.getHeight() - left.getHeight() + offset, threshold));
			}

			// check if right sibling exists
			if (node.getChildCount() > 1) {
				NetworkNode right = node.getChildEdges().get(1).childNode;
				if (node.getHeight() - right.getHeight() < threshold - offset) {
					if (right.isReassortment())
						reaForward.add(right);
					else {
						NetworkNode temp = closestReaAfterMig(right, node.getHeight() - right.getHeight() + offset,
								threshold);
						if (temp != null)
							reaForward.add(temp);
					}
//						reaForward.add(
//								closestReaAfterMig(right, node.getHeight() - right.getHeight() + offset, threshold));
				}

			}
			// no double counting
			if (!reaForward.isEmpty())
				return reaForward.stream()
						.sorted(Comparator.comparingDouble(n -> (node.getHeight() - n.getHeight() + offset))).distinct()
						.collect(Collectors.toList()).get(0);
			else
				return null;
		} else {
			return null;
		}
	}

	private List<NetworkNode> getMigrationNodesBackward(NetworkNode node, double offset, double threshold) {
		List<NetworkNode> migrationBackwards = new ArrayList<>();
		if (!node.getParentEdges().get(0).isRootEdge()) {
			NetworkNode left = node.getParentEdges().get(0).parentNode;
			if (left.getHeight() - node.getHeight() < threshold - offset) {
				if (isMigrationNode(left))
					migrationBackwards.add(left);
			migrationBackwards
					.addAll(getMigrationNodesBackward(left, offset + left.getHeight() - node.getHeight(), threshold));
			}

			// check if right parent exists
			if (node.getParentCount() > 1) {
				NetworkNode right = node.getParentEdges().get(1).parentNode;
				if (right.getHeight() - node.getHeight() < threshold - offset) {
					if (isMigrationNode(right))
						migrationBackwards.add(right);
				migrationBackwards
						.addAll(getMigrationNodesBackward(right, offset + right.getHeight() - node.getHeight(),
								threshold));
				}

			}

			// no double counting
			LinkedHashSet<NetworkNode> hashSet = new LinkedHashSet<>(migrationBackwards);
			migrationBackwards = new ArrayList<>(hashSet);
			return migrationBackwards.stream()
					.sorted(Comparator.comparingDouble(n -> (n.getHeight() - node.getHeight() + offset))).distinct()
					.collect(Collectors.toList());
		} else {
			return migrationBackwards;
		}

			
	}

	private boolean isMigrationNode(NetworkNode node) {
		if (!node.isLeaf() && node.getChildCount() == 1 && node.getParentCount() == 1)
			return true;
		else
			return false;
	}

	/**
	 * Copied from CoalRe TrunkReassortment.java by N.F. Muller and T. Vaughan
	 * 
	 * @param node
	 * @param dist
	 * @param threshold
	 */
	private void getAllAncestralEdgesLeaveDist(NetworkNode node, double dist, double threshold) {
		int index = allFitNodes.indexOf(node);
		if (index == -1) {
			allFitNodes.add(node);
			leaveDistance.add(dist);
		} else {
			if (leaveDistance.get(index) > dist)
				return;
			else if (leaveDistance.get(index) > threshold)
				return;
			else
				leaveDistance.set(index, dist);
		}

		for (NetworkEdge parentEdge : node.getParentEdges()) {
			if (parentEdge.isRootEdge()) {
				return;
			} else {
				getAllAncestralEdgesLeaveDist(parentEdge.parentNode, dist + parentEdge.getLength(), threshold);
			}
		}
	}

        
    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(NetworkAnnotatorOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Reassortment Event Trunk Mapper");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
		JLabel maxMigDistLabel = new JLabel("Max reassortment-to-migration distance:");
		JLabel minTipDistanceLabel = new JLabel("Min distance to tip for fit edge:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

		JTextField maxMigrationDistance = new JTextField(20);
		maxMigrationDistance.setEditable(true);
		JTextField minTipDistance = new JTextField(20);
		minTipDistance.setEnabled(true);

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);


        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel)
						.addComponent(burninLabel)
						.addComponent(maxMigDistLabel)
						.addComponent(minTipDistanceLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
						.addComponent(maxMigrationDistance)
						.addComponent(minTipDistance))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton))
                );

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
						.addComponent(maxMigDistLabel)
						.addComponent(maxMigrationDistance))
                .addGroup(layout.createParallelGroup()
						.addComponent(minTipDistanceLabel)
						.addComponent(minTipDistance))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
			options.maxMigrationDistance = Double.parseDouble(maxMigrationDistance.getText());
			options.minTipDistance = Double.parseDouble(minTipDistance.getText());
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select Reassortment Network log file to summarize");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
                outFileChooser.setCurrentDirectory(options.inFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });

        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Reassortment Event Locator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
			"Reassortment and migration analyzer - analyses reassortment and migration relation.\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-trunkDefinition {MostRecentSample, TipDistance} Choose trunk definition method.\n"
                    + "                         (default MostRecentSample)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
					+ "-maxReaMigDistance     	maximum distance between reassortment node and migration node\n"
					+ "                         In other words, maximum window width to consider.\n"
                    + "-minTipDistance     		minimum distance between internal network node\n"
                    + "                         and tip node such that the internal node is considered trunk.\n"
                    + "                         If not  specified, the trunk is any node between samples\n"
                    + "                         height=0 and the root.\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'reassortment_distances.txt'.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve TrunkReassortment options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, NetworkAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;


			case "-maxReaMigDistance":
                    if (args.length<=i+1) {
					printUsageAndError("-maxReaMigDistance must be followed by a number.");
                    }

                    try {
					options.maxMigrationDistance =
                                Double.parseDouble(args[i + 1]);
                    } catch (NumberFormatException e) {
					printUsageAndError("maxReaMigDistance must be a positive number. ");
                     }

                    i += 1;
                    break;
                    
			case "-minTipDistance":
				if (args.length <= i + 1) {
					printUsageAndError("-minTipDistance must be followed by a number.");
				}

				try {
					options.minTipDistance = Double.parseDouble(args[i + 1]);
				} catch (NumberFormatException e) {
					printUsageAndError("minTipDistance must be a positive number. ");
				}

				i += 1;
				break;

			case "-removeSegments":
                    if (args.length<=i+1) {
                        printUsageAndError("-removeSegments must be followed by at least one number.");
                    }

                    try {
                    	String[] argarray = args[i + 1].split(",");
                    	options.removeSegments = new int[argarray.length];
                    	for (int j = 0; j < argarray.length; j++)
                    		options.removeSegments[j] = Integer.parseInt(argarray[j]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("removeSegments must be an array of integers separated by commas if more than one");
                     }

                    i += 1;
                    break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	NetworkAnnotatorOptions options = new NetworkAnnotatorOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }


        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
			new ReassortmentAndMigration(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
}