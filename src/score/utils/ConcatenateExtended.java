package score.utils;

import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import feast.function.Concatenate;

@Description("A Function whose elements are the result of concatenating the elements of the input Functions.")
public class ConcatenateExtended extends Concatenate {

	public Input<List<RealParameter>> functionsInput = new Input<>("arg",
            "One or more functions to concatenate.",
            new ArrayList<>());

	List<RealParameter> functions;
    int totalDim;

    @Override
    public void initAndValidate() {
        functions = functionsInput.get();

        if (functions.isEmpty())
            throw new IllegalArgumentException("Concatenate expects at least one input function.");

        totalDim = 0;
        for (Function func : functionsInput.get())
            totalDim += func.getDimension();
    }

    @Override
    public int getDimension() {
        return totalDim;
    }

    @Override
    public double getArrayValue(int dim) {
        int paramIdx=0;
        while (dim >= functions.get(paramIdx).getDimension()) {
            dim -= functions.get(paramIdx).getDimension();
            paramIdx += 1;
        }

        return functions.get(paramIdx).getArrayValue(dim);
    }
}
