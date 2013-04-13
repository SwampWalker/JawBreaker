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
import ca.tonita.jawbreaker.equationsOfState.ElectronPositronPlasmaEOS;

/**
 * Reads the equation of state file, and adds leptons to it.
 * @author atonita
 *
 */
public class Logarithms implements InputOutput {
	
	private ElectronPositronPlasmaEOS eos;
	
	public Logarithms() {
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
		for (int i = 0; i < 4; i++) out.write(in.readLine() + "\n");
		
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
		// rho is already logarithm
		String outline = tokens[0] + "  ";
		
		// Get the density.
		double baryonDensity = Double.valueOf(tokens[1]);    // fm^-3
		outline += Math.log10(baryonDensity) + "  " + tokens[2] + "  ";
		
		for (int i = 3; i < 13; i++) outline += tokens[i] + "  ";
		
		double pressure = Math.abs(Double.valueOf(tokens[13])); // The yp=0 branch has some inconsistencies.
		outline += Math.log10(pressure) + "  ";
		
		for (int i = 14; i < tokens.length; i++) outline += tokens[i] + "  ";
		return outline + "\n";
	}
}
