package jaw.breaker.shenzerotemperature.drivers;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import jaw.breaker.equationsOfState.ElectronPositronPlasmaEOS;

/**
 * Reads the equation of state file, and adds leptons to it.
 *
 * @author atonita
 *
 */
public class AddLeptons implements InputOutput {

    private ElectronPositronPlasmaEOS eos;

    public AddLeptons() {
        eos = new ElectronPositronPlasmaEOS();
        eos.setUnits(CommonUnits.MEV);
    }

    @Override
    public int readAndWrite(File input, File output) {
        try {
            process(input, output);
        } catch (IOException e) {
            // This should never happen.
            e.printStackTrace();
        }
        return SUCCESS;
    }

    private void process(File input, File output) throws IOException {
        // A lot of work just to readLine() >:( That's boiler plate for you.
        FileInputStream rawInput = new FileInputStream(input);
        DataInputStream inData = new DataInputStream(rawInput);
        BufferedReader in = new BufferedReader(new InputStreamReader(inData));

        FileWriter fstream = new FileWriter(output);
        BufferedWriter out = new BufferedWriter(fstream);

        // First 4 lines are comments;
        for (int i = 0; i < 4; i++) {
            out.write(in.readLine() + "\n");
        }

        // 66 Yp blocks
        for (int j = 0; j < 66; j++) {
            // with 110 density lines
            for (int i = 0; i < 110; i++) {
                try {
                    String lineout = processLine(in.readLine());
                    out.write(lineout);
                } catch (ArrayIndexOutOfBoundsException e) {
                    System.err.println("i = " + i + " j = " + j);
                    throw e;
                }
            }
            // block seperator.
            out.write(in.readLine() + "\n");
        }

        // Just to be on the safe side.
        rawInput.close();
        inData.close();
        in.close();

        out.flush();
        fstream.flush();
        out.close();
        fstream.close();
    }

    private String processLine(String readLine) {
        String[] tokens = readLine.split("\\s+");
        // Token 0 is some leading spaces.
        // rho, n_B and Y_p don't change
        String outline = tokens[1] + "  " + tokens[2] + "  " + tokens[3] + "  ";

        // Get the density.
        double baryonDensity = Double.valueOf(tokens[2]);    // fm^-3
        double protonFraction = Double.valueOf(tokens[3]);
        double freeEnergy = Double.valueOf(tokens[4]);
        double internalEnergy = Double.valueOf(tokens[5]);
        double pressure = Double.valueOf(tokens[14]);        // MeV fm^-3

        double chargeDensity = baryonDensity * protonFraction; // fm^-3
        double muElectron = eos.chemicalPotential(chargeDensity);
        double deltaEnergy = eos.internalEnergyDensity(chargeDensity) / baryonDensity;
        freeEnergy += deltaEnergy; // Free energy is internal energy at T=0.
        internalEnergy += deltaEnergy;

        outline += freeEnergy + "  " + internalEnergy + "  ";
        for (int i = 6; i < 14; i++) {
            outline += tokens[i] + "  ";
        }

        pressure += eos.pressure(chargeDensity);
        outline += pressure + "  " + tokens[15] + "  " + tokens[16] + "  ";
        if (tokens.length < 18) { // No hyperons.
            outline += 0.0 + "  " + 0.0 + "  " + muElectron + "\n";
        } else {
            outline += tokens[17] + "  " + tokens[18] + "  " + muElectron + "\n";
        }
        return outline;
    }
}
