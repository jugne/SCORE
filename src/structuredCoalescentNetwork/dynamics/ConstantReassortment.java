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

    public Input<Dynamics> structuredCoalescentDynamics = new Input<>("structuredCoalescentDynamics",
	    "input of dynamics to use for structured coalescent. See MASCOT dynamics package.", Validate.REQUIRED);

    Dynamics scDynamics;

//    enum DynamicsTypes {
//	constant, bssvs, glm
//    };
//
//    public Input<DynamicsTypes> implementationInput = new Input<>("implementation",
//	    "implementation, one of " + DynamicsTypes.values().toString(), DynamicsTypes.constant);

    @Override
    public void initAndValidate() {
	super.initAndValidate();

	if (dimensionInput.get() < 1)
	    dimensionInput.set(getNrTypes());
	scDynamics = structuredCoalescentDynamics.get();
    }

    public double[] getReassortmentRate(int i) {
	if (reassortmentRates.get().getDimension() != dimensionInput.get()) {
	    System.err.println("Wrong dimension of reassortment rates input. "
		    + "Reassortment rates for all types set to the value: " + reassortmentRates.get().getArrayValue(0));
	    reassortmentRates.get().setDimension(dimensionInput.get());
	}
	double[] reassort = new double[dimensionInput.get()];
	for (int k = 0; k < dimensionInput.get(); k++) {
	    reassort[k] = reassortmentRates.get().getArrayValue(k);
	}

	return reassort;
    }

    @Override
    public void close(PrintStream out) {
	scDynamics.recalculate();
    }

    @Override
    public void recalculate() {
	scDynamics.recalculate();
    }

    @Override
    public double getInterval(int i) {
	return scDynamics.getInterval(i);
    }

    @Override
    public double[] getIntervals() {
	return scDynamics.getIntervals();
    }

    @Override
    public boolean intervalIsDirty(int i) {
	return scDynamics.intervalIsDirty(i);
    }

    @Override
    public double[] getCoalescentRate(int i) {
	return scDynamics.getCoalescentRate(i);
    }

    @Override
    public double[] getBackwardsMigration(int i) {
	return scDynamics.getBackwardsMigration(i);
    }

    @Override
    public int getEpochCount() {
	return scDynamics.getEpochCount();
    }

    @Override
    public void init(PrintStream out) {
	scDynamics.init(out);

    }

}
