package jaw.breaker.shenzerotemperature;

import java.io.File;
import javax.swing.UIManager;
import jaw.breaker.shenzerotemperature.drivers.CppClassMaker;
import jaw.breaker.shenzerotemperature.drivers.InputOutput;
import jaw.breaker.shenzerotemperature.gui.InputOutputView;

public class CppClassWriter {
	
	public static void usage() {
		System.out.println("Usage: CppClassMaker [table.t00 output]");
	}
	
	static String btnLabel = "Save data into class file.";
	static String helpText = "Pushing that button will open two dialog boxes successively. The first dialog box will ask for the zero temperature equation of state file, at beta equilibrium, with leptons and shear. The second dialog box will prompt for the output filename.";


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		InputOutput payload = new CppClassMaker();
		String file1 = "C:\\projects\\elasticity\\data\\eos\\eos3.shear.betaEquilibrium.2.t00";
		String file2 = "C:\\projects\\elasticity\\data\\eos\\Shen3BetaEquilibrium.cpp";
		
		String[] files = {file1, file2};
		
		if (files.length == 2) {
			// Then don't make the gui
			File input = new File(files[0]);
			File output = new File(files[1]);
			payload.readAndWrite(input, output);
		} else if (args.length == 0) {
			// Make the gui.
			try {
				UIManager.setLookAndFeel(
						UIManager.getSystemLookAndFeelClassName());
				InputOutputView gui = new InputOutputView(payload, btnLabel, helpText);
				gui.setVisible(true);
			} catch (Exception e) {
				e.printStackTrace();
			}
		} else {
			usage();
		}
	}

}
