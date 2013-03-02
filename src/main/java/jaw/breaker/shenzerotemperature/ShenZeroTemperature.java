package jaw.breaker.shenzerotemperature;

import java.io.File;
import javax.swing.UIManager;
import jaw.breaker.shenzerotemperature.drivers.AddLeptons;
import jaw.breaker.shenzerotemperature.drivers.InputOutput;
import jaw.breaker.shenzerotemperature.gui.InputOutputView;

public class ShenZeroTemperature {
	
	public static void usage() {
		System.out.println("Usage: ShenZeroTemperature [table.t00 output]");
		System.out.println("\tOpens Shen equation of state table table.t00, adds leptons and saves this to output.");
		System.out.println("\tCalling with no command line arguments will open the gui.");
	}
	

	static String btnLabel = "Adds leptons to zero-temperature equation of state table";
	static String helpText = "Pushing that button will open two dialog boxes successively."
                +" The first dialog box will ask for the zero temperature equation of state file, "
                +"which is the combination of the T=0 and y_p=0 files from Shen's website using "
                +"mergeT0YP0.py. The second dialog box will prompt for the output filename. Once "
                +"both files are selected, the table will be opened and the equation of state "
                +" contributions due to leptons will be added to the equation of state. The result"
                +"will be written to the output file. Take the logarithm next.";

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		InputOutput payload = new AddLeptons();
		
		if (args.length == 2) {
			// Then don't make the gui
			File input = new File(args[0]);
			File output = new File(args[1]);
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
