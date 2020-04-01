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

package structuredCoalescentNetwork.networkAnnotator;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.Arrays;
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
import coalre.networkannotator.ReassortmentLogReader;

/**
 * A rewrite of TreeAnnotator that outputs how often reassortment events happen on trunk branches vs. other branches 
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class ReassortmentAndMigration extends ReassortmentAnnotator {

    
	List<NetworkNode> reaAfterList;
	List<NetworkNode> reaBeforeList;
	List<NetworkNode> migrationForwardList;

	List<NetworkNode> migrationBackwardList;
	List<NetworkNode> migrationUnasignedList;

	List<Double> reaToTipLengths;

	List<Double> migToTipLengths;
	Double lengthForward;
	Double lengthBackward;

    private static class NetworkAnnotatorOptions {
        File inFile;
		File outFile = new File("reassortment_and_migration.txt");
        double burninPercentage = 10.0;

		double minMigrationDistance;
        int[] removeSegments = new int[0];


        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "%\n" +
					"Minimal distance before/after reassortment event " + minMigrationDistance;
        }
    }

	public ReassortmentAndMigration(NetworkAnnotatorOptions options) throws Exception {

        // Display options:
        System.out.println(options + "\n");

        // Initialise reader
        ReassortmentLogReader logReader = new ReassortmentLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");

	      System.out.println("\nWriting output to " + options.outFile.getName()
	      + "...");
        
		// compute number of migration within minMigrationDistance to reassotment
        try (PrintStream ps = new PrintStream(options.outFile)) {
			ps.print("n_reassortment \t" + "n_migration \t" + "n_after \t"
					+ "n_before \t" + "length_after \t" + "sum_rate_after \t" + "length_before \t"
					+ "sum_rate_before \t" + "network_length \t"
					+ "min_rea_tip_length \t" + "median_rea_tip_length \t" + "min_mig_tip_length \t"
					+ "median_mig_tip_length\t" + "rea_count" + "length" + "mig_tip_height");
			ps.print("\n");
			for (Network network : logReader) {
	        	pruneNetwork(network, options.removeSegments);
				computeReassortmentAndMigration(network, ps, options.minMigrationDistance);
	        	
	        	ps.print("\n");
	        }
	        ps.close();
        }
        System.out.println("\nDone!");
    }
        

	private void computeReassortmentAndMigration(Network network, PrintStream ps, double minMigrationDistance)
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
		double sumRateForward = 0.0;
		double sumRateBackward = 0.0;
		List<Double> migToTipHeight = new ArrayList<Double>();
		List<NetworkNode> reaCounts = new ArrayList<NetworkNode>();
		Double length = 0.0;
		HashMap<NetworkEdge, Double> lineageMap = new HashMap<NetworkEdge, Double>();
		

//		for (NetworkNode node : reaNodes) {
//			migrationForwardList.addAll(getMigrationNodesForward(node, 0.0, minMigrationDistance));
//			migrationBackwardList.addAll(getMigrationNodesBackward(node, 0.0, minMigrationDistance));
//		}
//		migrationForwardList = migrationForwardList.stream().distinct().collect(Collectors.toList());
//		migrationBackwardList = migrationBackwardList.stream().distinct().collect(Collectors.toList());
//		
//		ps.print(reaNodes.size() + "\t" + migNodes.size() + "\t" + migrationForwardList.size() + "\t"
//				+ migrationBackwardList.size());

		for (NetworkNode node : migNodes) {
			reaCounts.addAll(reaNodesBefore(node, 0.0, minMigrationDistance, length, lineageMap));


			
			NetworkNode reaBefore = closestReaBeforeMig(node, 0.0, minMigrationDistance);
			NetworkNode reaAfter = closestReaAfterMig(node, 0.0, minMigrationDistance);
			double heigthDiffBefore = Double.POSITIVE_INFINITY;
			double heigthDiffAfter = Double.POSITIVE_INFINITY;

			migToTipLengths.add(lengthToTip(node));
			migToTipHeight.add(node.getHeight() - closesDescendantTip(node).getHeight());
			double middlePoint = Double.POSITIVE_INFINITY;

			if (reaBefore != null)
				heigthDiffBefore = reaBefore.getHeight() - node.getHeight();
			if (reaAfter != null)
				heigthDiffAfter = node.getHeight() - reaAfter.getHeight();
			if (reaBefore != null && reaAfter != null)
				middlePoint = (reaBefore.getHeight() - reaAfter.getHeight()) / 2.0;
			
			if (heigthDiffBefore == Double.POSITIVE_INFINITY && heigthDiffAfter == Double.POSITIVE_INFINITY) {
				continue;
			}
			else if (heigthDiffBefore < heigthDiffAfter) {
				migrationForwardList.add(node);

//				double middlePoint = Double.POSITIVE_INFINITY;
//				NetworkNode nextRea = findDownStreamReas(reaBefore);
//				if (nextRea != null) {
//					middlePoint = (reaBefore.getHeight() - nextRea.getHeight()) / 2.0;
//				}
//				middlePoint = (reaBefore.getHeight() - reaAfter.getHeight()) /2.0;

				double minDist = middlePoint < minMigrationDistance ? middlePoint : minMigrationDistance;
				double l1 = lengthAfterRea(reaBefore, 0.0, minDist);
				if (!reaBeforeList.contains(reaBefore)) {
					reaBeforeList.add(reaBefore);
					lengthForward += l1;
				}
				sumRateForward += 1 / l1;

			} else {
				migrationBackwardList.add(node);

//				double middlePoint = Double.POSITIVE_INFINITY;
//				NetworkNode nextRea = findUpStreamReas(reaAfter);
//				if (nextRea != null) {
//					middlePoint = (nextRea.getHeight() - reaAfter.getHeight()) / 2.0;
//				}

				double minDist = middlePoint < minMigrationDistance ? middlePoint : minMigrationDistance;

				double l2 = lengthBeforeRea(reaAfter, 0.0, minDist);
				if (!reaAfterList.contains(reaAfter)) {
					reaAfterList.add(reaAfter);
					lengthBackward += l2;
				}
				sumRateBackward += 1 / l2;
			}
		}
		
		for (NetworkNode node : reaNodes) {
			reaToTipLengths.add(lengthToTip(node));
		}

		// get the length of the network
		List<NetworkEdge> allEdges = network.getEdges().stream()
				.filter(e -> !e.isRootEdge())
				.collect(Collectors.toList());

		double fullLength = 0.0;
		for (NetworkEdge edge : allEdges)
			fullLength += edge.getLength();

		Set<NetworkNode> set2 = new HashSet<NetworkNode>(reaCounts);
		if (set2.size() < reaCounts.size())
			throw new Exception("Duplicates detected!");

		ps.print(reaNodes.size() + "\t" + migNodes.size() + "\t" + migrationForwardList.size() + "\t"
				+ migrationBackwardList.size() + "\t" + lengthForward + "\t" + "\t"
				+ sumRateForward + "\t" + lengthBackward + "\t" + "\t" + sumRateBackward + "\t" +
				+ fullLength + "\t" + findMin(reaToTipLengths) + "\t" + findMedian(reaToTipLengths) + "\t"
				+ findMin(migToTipLengths) + "\t" + findMedian(migToTipLengths) + "\t"
				+ reaCounts.size() + "\t"
				+ length + "\t"
				+ Arrays.toString(migToTipHeight.toArray()));
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

	private List<NetworkNode> reaNodesBefore(NetworkNode node, double offset, double threshold, double length,
			HashMap<NetworkEdge, Double> edgeList)
			throws Exception {
		List<NetworkNode> reaBeforeList = new ArrayList<>();
		if (!node.getParentEdges().get(0).isRootEdge()) {
			NetworkEdge leftChildEdge = node.getParentEdges().get(0);
			NetworkNode left = node.getParentEdges().get(0).parentNode;

			if (left.getHeight() - node.getHeight() < threshold - offset) {
				length += left.getHeight() - node.getHeight();

				if (isMigrationNode(left)) {
					// reset window from new migration node
//					reaBeforeList.addAll(reaNodesBefore(left, 0.0, threshold, length));
				}
				if (left.isReassortment()) {
					// add the reassortment node
					reaBeforeList.add(left);
					// move further up
					reaBeforeList.addAll(
							reaNodesBefore(left, left.getHeight() - node.getHeight() + offset, threshold, length,
									edgeList));
				} else if (left.isCoalescence()) {
					// move further up
					reaBeforeList.addAll(
							reaNodesBefore(left, left.getHeight() - node.getHeight() + offset, threshold, length,
									edgeList));
				}
			} else {
				if ((threshold - offset) > 0)
					length += threshold - offset;
//					throw new Exception("offset greater than the threshold!!");
//				length += threshold - offset;
			}

				if (node.getParentCount() > 1) {
					NetworkNode right = node.getParentEdges().get(1).parentNode;
					if (right.getHeight() - node.getHeight() < threshold - offset) {
					length += right.getHeight() - node.getHeight();
						if (isMigrationNode(right)) {
							// reset window from new migration node
//						reaBeforeList.addAll(reaNodesBefore(right, 0.0, threshold, length));
					}
					if (right.isReassortment()) {
						// add the reassortment node
						reaBeforeList.add(right);
						// move further up
						reaBeforeList.addAll(
								reaNodesBefore(right, right.getHeight() - node.getHeight() + offset, threshold,
										length, edgeList));
					} else if (right.isCoalescence()) {
						// move further up
						reaBeforeList.addAll(
								reaNodesBefore(right, right.getHeight() - node.getHeight() + offset, threshold,
										length, edgeList));
						}
				} else
				if ((threshold - offset) > 0)
					length += threshold - offset;
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
        JLabel trunkDefinitionLabel = new JLabel("Trunk definition:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");

		JTextField minMigrationDistance = new JTextField(20);
		minMigrationDistance.setEditable(true);
//        minTipDistance.setEnabled(false);        

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
                        .addComponent(trunkDefinitionLabel))
//                        .addComponent(thresholdLabel)
//                        .addComponent(geneFlowCheckBox))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(burninSlider)
//                        .addComponent(thresholdSlider)
						.addComponent(minMigrationDistance))
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
						.addComponent(trunkDefinitionLabel))
                .addGroup(layout.createParallelGroup()
						.addComponent(minMigrationDistance))
                );

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
			options.minMigrationDistance = Double.parseDouble(minMigrationDistance.getText());
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
            "TrunkReassortment - counts how many reassortment events happened on trunk and non-trunk nodes.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-trunkDefinition {MostRecentSample, TipDistance} Choose trunk definition method.\n"
                    + "                         (default MostRecentSample)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
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


			case "-minMigrationDistance":
                    if (args.length<=i+1) {
					printUsageAndError("-minMigrationDistance must be followed by a number.");
                    }

                    try {
					options.minMigrationDistance =
                                Double.parseDouble(args[i + 1]);
                    } catch (NumberFormatException e) {
					printUsageAndError("minMigrationDistance must be a positive number. ");
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