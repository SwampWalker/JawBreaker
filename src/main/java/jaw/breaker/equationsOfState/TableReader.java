/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jaw.breaker.equationsOfState;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import atonita.unitconversion.dimensionalanalysis.Dimension;
import atonita.unitconversion.dimensionalanalysis.SIConstants;
import atonita.unitconversion.dimensionalanalysis.UnitSystem;
import java.io.*;
import java.util.ArrayList;

/**
 *
 * @author tonita
 */
public class TableReader {

    /**
     * Reads an equation of state from file, returns a <code>TabulateHermite</code> object.
     * @param eosFile The file to open.
     * @return the <code>TabulatedHermite</code> from the file
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public static TabulatedHermite readTableFromFile(File eosFile) throws FileNotFoundException, IOException {
        // Read the file.
        FileInputStream fis = new FileInputStream(eosFile);
        BufferedReader br = new BufferedReader(new InputStreamReader(fis));
        String line = br.readLine();
        ArrayList<String[]> tokenSet = new ArrayList<String[]>();
        int nLineOne = -1;
        while (line != null) {
            if (!line.startsWith("#")) {
                String[] tokens = line.split("\\s+");
                if (nLineOne == -1) {
                    nLineOne = tokens.length;
                }
                if (tokens.length == nLineOne) {
                    tokenSet.add(tokens);
                }
            }
            line = br.readLine();
        }
        fis.close();

        double[] logn = new double[tokenSet.size()];
        double[] logp = new double[tokenSet.size()];
        double[] energyPerParticle = new double[tokenSet.size()];
        
        // Unit conversion
        double nConversion = Math.log10(UnitSystem.convert(1, Dimension.NUMBERDENSITY, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO));
        double pConversion = Math.log10(UnitSystem.convert(1, Dimension.PRESSURE, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO));
        double eConversion = UnitSystem.convert(1, Dimension.ENERGY, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO);
        double particleMass = UnitSystem.convert(931.494, Dimension.MASS, CommonUnits.MEV, CommonUnits.GEOMETRICASTRO); // value from shen guide
        
        for (int i = 0; i < tokenSet.size(); i++) {
            String[] tokens = tokenSet.get(i);
            logn[i] = Double.valueOf(tokens[1]) + nConversion;
            logp[i] = Double.valueOf(tokens[13]) + pConversion;
            energyPerParticle[i] = Double.valueOf(tokens[4])*eConversion;
        }
        
        TabulatedHermite eos = new TabulatedHermite(logn, logp, energyPerParticle, particleMass);
        return eos;
    }
}
