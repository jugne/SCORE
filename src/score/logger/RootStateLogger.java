package score.logger;

import java.io.PrintStream;

import org.jblas.DoubleMatrix;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import mascot.dynamics.Dynamics;
import score.distribution.SCORE;

@Description("logs the state of the root, i.e. the probability the root being in "+ 
			"any of m states based on the MASCOT density class")
public class RootStateLogger extends CalculationNode implements Loggable {
	public Input<SCORE> scoreInput = new Input<SCORE>(
			"score",
			"");
	public Input<Dynamics> reassortmentDynamics = new Input<Dynamics>("reassortmentDynamics","");

	
	private int states;
	@Override
	public void init(PrintStream out) {
		DoubleMatrix RootStates = scoreInput.get().getRootTypes();
		states = RootStates.length;
		for (int i = 0 ; i < states; i++){
			out.print("RootProbability." + reassortmentDynamics.get().getStringStateValue(i) + "\t");
		}

	}


	@Override
	public void log(long sample, PrintStream out) {
		DoubleMatrix RootStates = scoreInput.get().getRootTypes();
		states = RootStates.length;
		for (int i = 0 ; i < states; i++){
			out.print(RootStates.get(i) + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {	
	}

	@Override
	public void initAndValidate() {	
	}	
	

	
}
