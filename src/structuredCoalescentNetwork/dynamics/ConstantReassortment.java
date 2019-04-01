package structuredCoalescentNetwork.dynamics;

import java.io.PrintStream;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;
import beast.mascot.dynamics.Dynamics;

public class ConstantReassortment extends Dynamics implements Loggable {

    public Input<RealParameter> reassortmentRates = new Input<>("reassortmentRates", "input of reassortment rates",
	    Validate.REQUIRED);

    enum DynamicsTypes {
	constant, bssvs, glm
    };

    public Input<DynamicsTypes> implementationInput = new Input<>("implementation",
	    "implementation, one of " + DynamicsTypes.values().toString(), DynamicsTypes.constant);

    @Override
    public void initAndValidate() {
	super.initAndValidate();

    }

    public double[] getReassortmentRate(int i) {
	double[] coal = new double[NeInput.get().getDimension()];
	int c = 0;
	for (int j = 0; j < NeInput.get().getDimension(); j++) {
	    reassort[c] = 1 / (ploidyInput.get() * NeInput.get().getArrayValue(j));
	    c++;
	}
	return reassort;

    }

    @Override
    public void close(PrintStream out) {
	// TODO Auto-generated method stub

    }

    @Override
    public void recalculate() {
	// TODO Auto-generated method stub

    }

    @Override
    public double getInterval(int i) {
	// TODO Auto-generated method stub
	return 0;
    }

    @Override
    public double[] getIntervals() {
	// TODO Auto-generated method stub
	return null;
    }

    @Override
    public boolean intervalIsDirty(int i) {
	// TODO Auto-generated method stub
	return false;
    }

    @Override
    public double[] getCoalescentRate(int i) {
	// TODO Auto-generated method stub
	return null;
    }

    @Override
    public double[] getBackwardsMigration(int i) {
	// TODO Auto-generated method stub
	return null;
    }

    @Override
    public int getEpochCount() {
	// TODO Auto-generated method stub
	return 0;
    }

}
