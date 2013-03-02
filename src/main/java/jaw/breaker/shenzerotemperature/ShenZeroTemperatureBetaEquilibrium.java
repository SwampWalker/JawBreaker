package jaw.breaker.shenzerotemperature;

import java.io.File;
import javax.swing.UIManager;
import jaw.breaker.shenzerotemperature.drivers.BetaEquilibrium;
import jaw.breaker.shenzerotemperature.drivers.InputOutput;
import jaw.breaker.shenzerotemperature.gui.InputOutputView;

public class ShenZeroTemperatureBetaEquilibrium {

    public static void usage() {
        System.out.println("Usage: ShenZeroTemperatureBetaEquilibrium [table.t00 output]");
        System.out.println("\tOpens Shen equation of state table logtable.t00, finds beta equilibrium and saves this to output.");
        System.out.println("\tCalling with no comman line arguments will open the gui.");
    }
    static String btnLabel = "Find beta equilibrium in zero-temperature equation of state";
    static String helpText = "Pushing that button will open two dialog boxes successively. "
            + "The first dialog box will ask for the logarithmic zero temperature equation of state file"
            + " which has had the leptons added to it, the beta equilibrium of this equation of state will"
            + " then be computed and saved into the second file.";

    /**
     * @param args
     */
    public static void main(String[] args) {
        InputOutput payload = new BetaEquilibrium();
        String file1 = "D:\\elasticity\\data\\eos\\eos3.log.leptons.yp0.t00";
        String file2 = "D:\\elasticity\\data\\eos\\eos3.betaEquilibrium.2.t00";

        String[] files = {/*file1, file2*/};

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
