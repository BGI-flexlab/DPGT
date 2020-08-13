package org.bgi.flexlab.gaea.tools.jointcalling.util;

import java.util.Collections;

import org.apache.commons.math3.special.Gamma;

public class Dirichlet {
    final double[] alpha;

    public Dirichlet(final double... alpha) {
    	JointCallingUtils.nonNull(alpha);
    	JointCallingUtils.validateArg(alpha.length >= 1, "Dirichlet parameters must have at least one element");
    	JointCallingUtils.validateArg(GvcfMathUtils.allMatch(alpha, x -> x >= 0), "Dirichlet parameters may not be negative");
    	JointCallingUtils.validateArg(GvcfMathUtils.allMatch(alpha, Double::isFinite), "Dirichlet parameters must be finite");
        this.alpha = alpha.clone();
    }

    /**
     * Create a symmetric distribution Dir(a/K, a/K, a/K . . .) where K is the number of states and
     * a is the concentration.
     */
    public static Dirichlet symmetricDirichlet(final int numStates, final double concentration) {
    	JointCallingUtils.validateArg(numStates > 0, "Must have at leat one state");
    	JointCallingUtils.validateArg(concentration > 0, "concentration must be positive");
        return new Dirichlet(Collections.nCopies(numStates, concentration/numStates).stream().mapToDouble(x->x).toArray());
    }

    // in variational Bayes one often needs the effective point estimate of a multinomial distribution with a
    // Dirichlet prior.  This value is not the mode or mean of the Dirichlet but rather the exp of the expected log weights.
    // note that these effective weights do not add up to 1.  This is fine because in any probabilistic model scaling all weights
    // amounts to an arbitrary normalization constant, but it's important to keep in mind because some classes may expect
    // normalized weights.  In that case the calling code must normalize the weights.
    public double[] effectiveMultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(GvcfMathUtils.sum(alpha));
        return GvcfMathUtils.applyToArray(alpha, a -> Math.exp(Gamma.digamma(a) - digammaOfSum));
    }

    public double[] effectiveLog10MultinomialWeights() {
        final double digammaOfSum = Gamma.digamma(GvcfMathUtils.sum(alpha));
        return GvcfMathUtils.applyToArray(alpha, a -> (Gamma.digamma(a) - digammaOfSum) * GvcfMathUtils.LOG10_OF_E);
    }

    public double[] meanWeights() {
        final double sum = GvcfMathUtils.sum(alpha);
        return GvcfMathUtils.applyToArray(alpha, x -> x / sum);
    }

    public double[] log10MeanWeights() {
        final double sum = GvcfMathUtils.sum(alpha);
        return GvcfMathUtils.applyToArray(alpha, x -> Math.log10(x / sum));
    }
    
    public int size() { return alpha.length; }
}
