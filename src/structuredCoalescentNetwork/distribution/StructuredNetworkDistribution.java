package structuredCoalescentNetwork.distribution;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import coalre.network.Network;
import beast.core.State;

import java.util.List;
import java.util.Random;

public class StructuredNetworkDistribution extends Distribution {
	public Input<Network> networkIn = new Input<>("network", "nework over which to calculate a prior or likelihood");
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
        return networkIn.get().somethingIsDirty();
//    	return true; problem not here
    }

}
