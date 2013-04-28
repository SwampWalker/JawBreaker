package ca.tonita.math.numerical;

/**
 *
 * @author atonita
 */
public class ODEIndexer1D {

    private int[] rank;
    private int[] variables;
    private int parameters;
    private int domains;

    /**
     * Creates an indexer. When all the variables are collected into a single
     * vector, this class is capable of indexing that variable.
     *
     * @param rank The rank of the basis in each domain.
     * @param variables The number of variables per domain.
     * @param parameters The number of parameters of the problem.
     * @param domains The number of domains.
     */
    public ODEIndexer1D(int[] rank, int[] variables, int parameters, int domains) {
        this.rank = rank;
        this.variables = variables;
        this.parameters = parameters;
        this.domains = domains;
    }

    /**
     * Returns the index into the data vector for the specific variable.
     * Specifically, this returns V_iVariable(x_iAbscissa) \in \Omega_iDomain.
     *
     * @param iDomain The domain of the variable.
     * @param iVariable The index of the variable.
     * @param iAbscissa The index of the coordinate location of the variable.
     * @return The index of the variable in the data vector.
     */
    public int index(int iDomain, int iVariable, int iAbscissa) {
        int index = 0;
        for (int i = 0; i < iDomain; i++) {
            index += variables[i] * rank[i];
        }
        index += iVariable * rank[iDomain] + iAbscissa;
        return index;
    }

    /**
     * Returns the index into the data vector of the parameter/constraint.
     * @param iParameter The parameter index.
     * @return the corresponding index in the data vector.
     */
    public int index(int iParameter) {
        int index = 0;
        for (int i = 0; i < domains; i++) {
            index += variables[i] * rank[i];
        }
        return index + iParameter;
    }

    public int[] indices(int iVector) {
        int iDomain = 0;
        int iVariable = 0;
        int iAbscissa = 0;
        int index = 0;
        while (iDomain < domains && iVector > index + variables[iDomain] * rank[iDomain]) {
            iDomain++;
            index += variables[iDomain] * rank[iDomain];
        }
        if (iDomain == domains) {
            // iVector corresponds to a parameter.
            return new int[]{iVector - index};
        }
        throw new UnsupportedOperationException("Not finished yet.");
    }

    public int getNVariables(int iDomain) {
        return variables[iDomain];
    }
}
