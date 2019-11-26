package structuredCoalescentNetwork.mapping;

import java.io.PrintStream;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;

@Description("logs the state of the root, i.e. the probability the root being in "+ 
			"any of m states based on the MASCOT density class")
public class MappedLogger extends CalculationNode implements Loggable {
	public Input<MappedNetwork> mappedInput = new Input<MappedNetwork>(
			"mappedNetwork",
			"");

	
	private int states;
	MappedNetwork mappedNet;
	@Override
	public void init(PrintStream out) {
		mappedNet = mappedInput.get();
//		DoubleMatrix RootStates = scoreInput.get().getRootTypes();
//		states = RootStates.length;
//		for (int i = 0 ; i < states; i++){
//			out.print("RootProbability." + reassortmentDynamics.get().getStringStateValue(i) + "\t");
//		}

	}

	@Override
	public void log(int sample, PrintStream out) {
		mappedNet.doStochasticMapping();
//		DoubleMatrix RootStates = scoreInput.get().getRootTypes();
//		states = RootStates.length;
//		for (int i = 0 ; i < states; i++){
//			out.print(RootStates.get(i) + "\t");
//		}
	}

	@Override
	public void close(PrintStream out) {	
	}

	@Override
	public void initAndValidate() {	
	}	
	

	
}
