package ca.tonita.jawbreaker.shenzerotemperature.drivers;

import atonita.unitconversion.dimensionalanalysis.CommonUnits;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import ca.tonita.jawbreaker.equationsOfState.StrohmayerShearModulus;

public class AddShear implements InputOutput {

    StrohmayerShearModulus shear;

    public AddShear() {
        shear = new StrohmayerShearModulus();
        shear.setUnits(CommonUnits.MEV);
    }

    /**
     * Entry point of the routine.
     *
     * @param input The file to read input from, must have had leptons added.
     * @param output The file to store output in.
     * @return Any return value other than SUCCESS is clearly not one.
     */
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

        // First line is comment;
        for (int i = 0; i < 1; i++) {
            out.write(in.readLine() + "\n");
        }

        // 1 Yp blocks
        // with 110 density lines
        for (int i = 0; i < 110; i++) {
            try {
                String lineout = processLine(in.readLine());
                out.write(lineout);
            } catch (ArrayIndexOutOfBoundsException e) {
                System.err.println("i = " + i);
                throw e;
            }
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
        String outline = "";
        for (int i = 0; i < tokens.length; i++) {
            outline += tokens[i] + "  ";
        }

        // Now compute the shear modulus.
        double density = Math.pow(10., Double.valueOf(tokens[InputOutput.numberdensity]));
        double Z = Double.valueOf(tokens[InputOutput.chargeNumber]);
        double A = Double.valueOf(tokens[InputOutput.massNumber]);
        double muOverN = this.shear.shearModulusOverNumberDensity(density, Z, A);

        return outline + muOverN + "\n";
    }
}
