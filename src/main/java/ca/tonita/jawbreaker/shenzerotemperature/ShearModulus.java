package ca.tonita.jawbreaker.shenzerotemperature;

import java.io.File;
import javax.swing.UIManager;
import ca.tonita.jawbreaker.shenzerotemperature.drivers.AddShear;
import ca.tonita.jawbreaker.shenzerotemperature.drivers.InputOutput;
import ca.tonita.jawbreaker.shenzerotemperature.gui.InputOutputView;

public class ShearModulus {

    public static void usage() {
        System.out.println("Usage: ShearModulus [table.t00 output]");
        System.out.println("\tOpens Shen beta equilibrium equation of state table table.t00, adds shear modulus.");
        System.out.println("\tCalling with no comman line arguments will open the gui.");
    }
    static String btnLabel = "Add shear to beta equilibrium equation of state table";
    static String helpText = "Pushing that button will open two dialog boxes successively. "
            + "The first dialog box will ask for the beta equilibrium equation of state file,"
            + " as output by the beta-equilibrium code. The output is saved in the second file chosen.";

    /**
     * @param args
     */
    public static void main(String[] args) {
        InputOutput payload = new AddShear();
        String file1 = "C:\\projects\\elasticity\\data\\eos\\eos3.betaEquilibrium.2.t00";
        String file2 = "C:\\projects\\elasticity\\data\\eos\\eos3.shear.betaEquilibrium.2.t00";

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
