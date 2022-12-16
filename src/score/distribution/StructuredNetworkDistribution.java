package score.distribution;

import java.util.List;
import java.util.Random;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;

public class StructuredNetworkDistribution extends Distribution {
//	public Input<Network> networkIn = new Input<>("network", "nework over which to calculate a prior or likelihood");
    public Input<StructuredNetworkIntervals> networkIntervalsInput = new Input<>("networkIntervals",
            "Structured Intervals for a phylogenetic beast tree", Validate.REQUIRED);

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
    
    @Override
    protected boolean requiresRecalculation() {
        final StructuredNetworkIntervals ti = networkIntervalsInput.get();
        if (ti != null) {
            assert ti.isDirtyCalculation();
            return true;
        }
		return ti.networkInput.get().somethingIsDirty();
//        return networkIn.get().somethingIsDirty();
//    	return true; problem not here
    }

}
