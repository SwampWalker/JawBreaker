package ca.tonita.jawbreaker.shenzerotemperature.drivers;

import java.io.File;

/**
 * Just outputs the files to screen.
 * @author atonita
 *
 */
public class Debug implements InputOutput {

	@Override
	public int readAndWrite(File input, File output) {
		System.out.println("Opening " + input.getName() + " and writing to " + output.getName());
		return 0;
	}

}
